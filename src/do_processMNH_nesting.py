from __future__ import print_function
from builtins import str
import sys
import numpy as np
import pdb 
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.interpolate
import os
import subprocess
from scipy import interpolate
import imp
import asciitable 
import glob 
import math 
import socket 
import argparse
import f90nml 
import shutil 

path_processMNH = os.environ['PATH_SRC_PYTHON_LOCAL']+'/Process_MNH/src/'
sys.path.append(path_processMNH)
sys.path.append(path_processMNH+'writePlt/')
#sys.path.append(path_processMNH+'../../Compare_FDS_MNH/src/')

#homebrewed
import MNH_tools
#import writePlt

#####################################################
def str_to_boolean(input):
    if  input in ['true', 'TRUE' , 'True' , '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh']:
        return True
    else:
        return False


#######################################
if __name__ == '__main__':
#######################################

    #reload(writePlt)
    reload(MNH_tools)

    parser = argparse.ArgumentParser(description='this is the driver of the processMNH algo.')
    parser.add_argument('-i','--input', help='Input Dir Simulation',required=True)
    parser.add_argument('-d','--domain', help='Number of the domain simulation behind the EXSEG, 1 for Father, 2 for Son1, 3 for Son2 ...',required=True)    
    parser.add_argument('-t','--time_ref', help='Reference time',required=True)
    parser.add_argument('-f','--format', help='format Output, availabale format are tecplot and netcdf. example: -f plt,nc. Default -f nc',required=False)
    parser.add_argument('-v','--verbose', help='verbose, default -v 0, minimum print',required=False)
    parser.add_argument('-r','--remove', help='remove previous output file, default -r False, minimum print',required=False)
    parser.add_argument('-dp','--dx_along_plume', help='step in meter used in the creation of the plume hull, default -dp 200',required=False)
    args = parser.parse_args()
    
 
    #in param
    flag_plot_lambda2 = True
    dxyz_factor = 1 

    #input directory 
    if args.input.isdigit():
        if args.input == '1':
            dir_data = path_processMNH + '../data/' # test directory
        else:
            print('this number not defined')
            sys.exit()
    else:
        dir_data = args.input
    if dir_data[-1] != '/': 
        dir_data = dir_data + '/'

    #reference time & reference domain

    time_ref = float(args.time_ref)
    domain   = str(args.domain)

    #format output 
    flag_plt = False
    flag_nc  = True
    if args.format is not None:
        if 'plt' in args.format:
            flag_plt = True
        else:
            flag_plt = False
        if 'nc' in args.format:
            flag_nc = True
        else:
            flag_nc = False
   
    dx_along_plume = args.dx_along_plume
    if dx_along_plume is None: 
        dx_along_plume = 200

    list_outputFormat = []
    if flag_plt: 
        list_outputFormat.append('plt')
    if flag_nc:
        list_outputFormat.append('nc')
    print('output format are: ', ' ,'.join(list_outputFormat))

    #verbose 
    verbose = 0
    if args.verbose is not None:
        verbose = np.int(args.verbose)
    
    #bookkeeping 
    flag_cleanStart = False
    if args.remove is not None:
        flag_cleanStart = str_to_boolean(args.remove)

    if not(os.path.isfile(dir_data+'/EXSEG1.nam')):
        print('missing EXSEG1.nam file in ', dir_data)
        print('stop here')
        sys.exit()

    #required input to compute the hull. use same as in compare_xxx code
    #inputConfig = imp.load_source('inputConfig' ,\
    #                               path_processMNH+'../../Compare_FDS_MNH/src/config.py')
    #dx_along_plume = inputConfig.inputdata['section_step']

    dir_postproc = dir_data + 'Postproc/'
    MNH_tools.ensure_dir(dir_postproc)
    outputDir    = dir_postproc + 'outputFile/'
    outputDir_nc = dir_postproc + 'outputFile/Netcdf/'
    if (flag_cleanStart) & (os.path.isdir(outputDir)):
        shutil.rmtree(outputDir)
    MNH_tools.ensure_dir(outputDir)
    MNH_tools.ensure_dir(outputDir_nc)
 

    if flag_plt:
        #init tecplot env
        writePlt.global_var.timestr_id = 1


    #load mnh namlist
    mnh_namelist  = f90nml.read(dir_data+'/EXSEG1.nam')

    #read namelist info
    expr_name = mnh_namelist['NAM_CONF']['CEXP'][:5]
    #set up passive tracer name 
    try:
        prep_namelist = f90nml.read(dir_data+'../0'+domain+'_pgd'+domain+'/PRE_PGD1.nam')
        nsv_user = prep_namelist['NAM_CONFn']['NSV_USER']
    except: 
        nsv_user = 0
    try:
        nsv_ff = mnh_namelist['nam_forefire']['nffscalars']
    except KeyError: 
        nsv_ff = 0
    ffsv_name = ''; ff_passiveTracer_info = None
    if nsv_ff > 0:
        ffsv_name = mnh_namelist['nam_forefire']['ffsvnames']
        item = ('mm',0)
        ff_passiveTracer_info = np.array([item]*nsv_ff,dtype=np.dtype([('name','S100'),('id_mnh',int)]))
        ff_passiveTracer_info = ff_passiveTracer_info.view(np.recarray)
        ff_passiveTracer_info.name = ffsv_name
        ff_passiveTracer_info.id_mnh   = nsv_user + np.arange(nsv_ff) + 1

    #get mnh output file
    mesonh_filenames  = sorted(glob.glob(dir_data+expr_name+'.'+domain+'*.nc4'))
    if len(mesonh_filenames) == 0: 
        print('could not find mnh files')
        print('in ', dir_data+expr_name+'.'+domain+'*.nc4')
        sys.exit()


    #loop over the files
    print('**')
    print('warning: first time iteration as shown below should be reference state')
    print('**')
    print('loop over {:d} MNH files in {:s}'.format(len(mesonh_filenames[1:]), dir_data))
    indent = '    '
    if flag_plt:
        time_scene_time_strand_tecplot = np.zeros([2,len(mesonh_filenames[1:])])
   
    i_mnh_file = 0
    for mesonh_filename in mesonh_filenames[1:]:
       
        #load mesonh data
        #-----------------
        time_scene, FireScene, FireScene2D, z_mesh = MNH_tools.load_centered_MesoNHField(mesonh_filename,     \
                                                                                 ff_sv=ff_passiveTracer_info, \
                                                                                 dxyz_factor= dxyz_factor,    \
                                                                                 verbose=verbose,indent=indent)
        if time_scene < time_ref:
            continue
        print(indent,os.path.basename(mesonh_filename), end=' ') 

        #set y axis origin at the center
        extent_y = (FireScene.yc+.5*FireScene.dy).max() -  (FireScene.yc-.5*FireScene.dy).min()
        FireScene.yc   = FireScene.yc   - .5* extent_y 
        FireScene2D.yc = FireScene2D.yc - .5* extent_y 
        
        print(time_scene.data)      
        
        if flag_plt:
            time_scene_time_strand_tecplot[0,i_mnh_file] = time_scene.data
            time_scene_time_strand_tecplot[1,i_mnh_file] = writePlt.global_var.timestr_id


        if i_mnh_file == 0: 
            
            # save vertical mesh
            #assume homogeneous mesh here
            f = open(outputDir+'vertical_mesh_edges.txt','w')
            temp_ = []
            for z_ in z_mesh[0,0,:]: # assume homogeneous mesh 
                temp_.append("{:18.6f}\n".format(z_))
            f.writelines(temp_)
            f.close()
            
            # and set reference scene
            FireScene_ref = np.copy(FireScene)
            FireScene_ref = FireScene_ref.view(np.recarray)

        if flag_plt:
            if verbose > 0: print(indent, 'write plt')
            # write tecplot file
            #-----------------
            ngx, ngy, ngz = FireScene.shape
            xg = FireScene.xc
            yg = FireScene.yc
            zg = FireScene.zc
            
            #list of variable
            vara = list(FireScene.dtype.names)
            vel_var = ['u','v','w']
            var_to_keep = []
            for ivar, var in enumerate(vara):
                if var not in ['xc','yc','zc']: # remove dim
                    var_to_keep.append(var)
            var_to_keep.append('plume_mask')

            for ivar,var in enumerate(var_to_keep):
                if var in vel_var:
                    i_vel_var = vel_var.index(var)
                    var_to_keep.insert(i_vel_var, var_to_keep.pop(var_to_keep.index(var)))
            
            str_var_to_keep = ' '.join(var_to_keep)
            nvar = len(var_to_keep)
            
            QQ = np.zeros([ngx, ngy, ngz,nvar])
            for ivar,var in enumerate(var_to_keep):
                if var in FireScene.dtype.names:
                    QQ[:,:,:,ivar] = FireScene[var]
                
                elif var == 'plume_mask':
                   if i_mnh_file == 0: # no plume in ref state
                        QQ[:,:,:,ivar] = np.zeros([ngx, ngy, ngz])
                   else:
                        QQ[:,:,:,ivar] = MNH_tools.plume_hull(dx_along_plume,FireScene, FireScene_ref, \
                                                                'SFObs', 0.02 , ref_val=1.e-3,flag='use xc,yc,zc keys')
                else:
                    print('issue in var selection, var =', var)
                    pdb.set_trace()
            
            writePlt.teclib.save_flow(expr_name,len(expr_name),flag_plot_lambda2,        \
                                  len(outputDir),outputDir,             \
                                  time_scene,                           \
                                  ngx, ngy, ngz,nvar,xg,yg,zg,          \
                                  QQ,                                   \
                                  len(str_var_to_keep),str_var_to_keep  )
            

            var_to_keep = ['heatFlux']
            QQ2d = np.zeros([ngx, ngy, nvar])
            for ivar,var in enumerate(var_to_keep):
                QQ2d[:,:,ivar] = FireScene2D[var]
            
            str_var_to_keep = ' '.join(var_to_keep)
            nvar2d = len(var_to_keep)

            writePlt.teclib.save_surface(expr_name,len(expr_name),        \
                                  len(outputDir),outputDir,             \
                                  time_scene,                           \
                                  ngx, ngy,nvar2d,FireScene2D.xc,FireScene2D.yc,FireScene2D.orography,          \
                                  QQ2d,                                   \
                                  len(str_var_to_keep),str_var_to_keep  )
            
            #increment time strandID
            writePlt.global_var.timestr_id = i_mnh_file + 2
       
            
        if flag_nc:
            if verbose > 0: print(indent, 'write netcdf')
            # write netcdf file
            #-----------------
            if i_mnh_file == 0: 
                flag_write = 'init'
            else:
                flag_write = 'append'
            MNH_tools.dump_netcdf(outputDir_nc,expr_name+domain,time_scene,FireScene,FireScene2D,\
                                  ff_sv=ff_passiveTracer_info,\
                                  flag_write=flag_write)
   
        #end loop
        i_mnh_file+=1

    if flag_plt:
        np.save(outputDir+'time_plt_strand',time_scene_time_strand_tecplot)
     
