from __future__ import print_function
from builtins import zip
import sys
import numpy as np
import pdb 
#import matplotlib.pyplot as plt
#from matplotlib import cm
import scipy.interpolate
import os
import subprocess
from scipy import interpolate
import imp
import glob 
import math 
import socket 
import argparse
import f90nml 
import shutil 
import importlib 
import datetime
import xarray as xr 


if 'PATH_SRC_PYTHON_LOCAL' not in os.environ: 
    print('script need to be call from where is is the python script as PATH_SRC_PYTHON_LOCAL is not set in env variable')
    os.environ['PATH_SRC_PYTHON_LOCAL'] = '../../'

path_processMNH = os.environ['PATH_SRC_PYTHON_LOCAL']+'/process_MNH/src/'
sys.path.append(path_processMNH)
sys.path.append(path_processMNH+'writePlt/')
#sys.path.append(path_processMNH+'../../Compare_FDS_MNH/src/')

#homebrewed
import MNH_tools
import writePlt

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
    importlib.reload(MNH_tools)

    parser = argparse.ArgumentParser(description='this is the driver of the processMNH algo.')
    parser.add_argument('-i','--input', help='Input Dir Simulation',required=True)
    parser.add_argument('-var','--variable', nargs="*", help='list of variable to save. need to give pair of mnh_name ouput_name, ex: -var HBLTOP,hbltop,lambda2,lambda2  lambda2 is set 2 times as it is computed here, hbltop need to be in dirDiag, see below',required=False)
    parser.add_argument('-tbeg','--time_beg', help='MNH time where concatenation start, all file before are skipped',required=False)
    parser.add_argument('-t0','--time_ref', help='Reference time to substract to MNH time',required=False)
    parser.add_argument('-tend','--time_end', help='end time, time elapsed after start time',required=False)
    parser.add_argument('-f','--format', help='format Output, availabale format are tecplot and netcdf. example: -f plt,nc. Default -f nc',required=False)
    parser.add_argument('-v','--verbose', help='verbose, default -v 0, minimum print',required=False)
    parser.add_argument('-r','--remove', help='remove previous output file, default -r False, minimum print',required=False)
    parser.add_argument('-dp','--dx_along_plume', help='step in meter used in the creation of the plume hull, default -dp 200',required=False)
    parser.add_argument('-m','--model', help='model number: 1(default), 2, 3, ..,',required=False)
    parser.add_argument('-zlim','--zlimit', help='max altitude to keep in the nc',required=False)
    parser.add_argument('-dirDiag','--dirDiag', help='directory where to look for the associated diag file, default is ../08_diag/, ',required=False)
    parser.add_argument('-out','--outFile', help='set to True to use high freq files',required=False)
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

    if args.model is None:
        modelN = 1
    else:
        modelN = int(args.model)
   
    if args.dirDiag is None:
        dirDiag =dir_data.replace('07_mnh', '08_diag')
    else: 
        dirDiag = args.dirDiag

    if args.zlimit is None:
        zlimit = None
    else:
        zlimit = float(args.zlimit)

    extraPattern = '*'
    if args.outFile is not None: 
        if str_to_boolean(args.outFile): 
            extraPattern='*OUT*'


    #reference time
    time_beg = float(args.time_beg) if (args.time_beg is not None) else 0
    time_ref = float(args.time_ref) if (args.time_ref is not None) else 0.
    time_end  = float(args.time_end) if (args.time_end is not None) else 1.e12

    #extar var
    if args.variable is not None:
        #listExtravariable = list(zip(args.variable[::2],args.variable[1::2]))
        variable_ = args.variable[0].split(',')
        listExtravariable = list(zip(variable_[::2],variable_[1::2]))
    else:
        listExtravariable = []

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
    else:
        dx_along_plume = int(dx_along_plume)
    
    list_outputFormat = []
    if flag_plt: 
        list_outputFormat.append('plt')
    if flag_nc:
        list_outputFormat.append('nc')
    print('output format are: ', ' ,'.join(list_outputFormat))

    #verbose 
    verbose = 0
    if args.verbose is not None:
        verbose = int(args.verbose)
    
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
    expr_name = mnh_namelist['NAM_CONF']['CEXP']
    expr_name_save = expr_name
    if extraPattern=='*OUT*': 
        expr_name_save += '_out'
    

    #Deal with extra var
    extraVar_info = None
    nextra = len(listExtravariable)

    #set up passive tracer name 
    try:
        prep_namelist = f90nml.read(dir_data+'../01_prep_ideal_case/PRE_IDEA1.nam')
        nsv_user = prep_namelist['NAM_CONFn']['NSV_USER']
    except: 
        nsv_user = 0
    ffsv_name = ''
    try:
        nsv_ff = mnh_namelist['nam_forefire']['nffscalars']-1
        ffsv_name = mnh_namelist['nam_forefire']['ffsvnames']
        ffsv_name.remove('BRatio')
    except KeyError: 
        nsv_ff = 0
   
    if nsv_ff + nextra > 0:
    
        item = ('mm','mm',int,'mm')
        extraVar_info = np.array([item]*(nsv_ff+nextra),dtype=np.dtype([('name','U100'),('name_mnh','U100'),('name_dtype',np.dtype),('dim','U100')]))
        extraVar_info = extraVar_info.view(np.recarray)
        
        if nsv_ff > 0:
            for i_, ffsv_name_ in enumerate(ffsv_name):
                extraVar_info.name[:nsv_ff][i_]       = ffsv_name_
                extraVar_info.name_mnh[:nsv_ff][i_]   = 'SVT{:03d}'.format(nsv_user + i_ + 2)
                extraVar_info.name_dtype[:nsv_ff][i_]      = float
                extraVar_info.dim[:nsv_ff][i_]      = '3D'

        if nextra > 0:
            for i_, var_ in enumerate(listExtravariable):
                extraVar_info.name[nsv_ff:][i_]       = var_[1]
                extraVar_info.name_mnh[nsv_ff:][i_]   = var_[0]
                extraVar_info.name_dtype[nsv_ff:][i_] = float
                extraVar_info.dim[nsv_ff:][i_] = 'na'


    #get mnh output file
    mesonh_filenames  = sorted(glob.glob(dir_data+expr_name+'*.{:1d}.{:s}.nc'.format(modelN,extraPattern)))
    if extraPattern == '*':
        mesonh_filenames_copy = mesonh_filenames[:]
        for file_ in mesonh_filenames_copy: 
            if '.OUT.' in file_:
                mesonh_filenames.remove(file_)

    if len(mesonh_filenames) == 0: 
        print('could not find mnh files')
        print('in ', dir_data)
        sys.exit()

    #remove init simulation files
    mesonh_filenames_ = []
    mesonh_filenames_time = []
    for mesonh_filename in mesonh_filenames:
        if '.spa' in mesonh_filename: continue 
        if int(os.path.basename(mesonh_filename).split('.')[-2]) == 0: continue
        mesonh_filenames_.append(mesonh_filename)
        ds = xr.open_dataset(mesonh_filename)
        mesonh_filenames_time.append( ds["time"].values[0]  )
    
    mesonh_filenames = np.array(mesonh_filenames_)[ np.argsort(mesonh_filenames_time)]


    #loop over the files
    if extraVar_info is not None:
        if 'plume_mask' in extraVar_info.name: 
            print('**')
            print('warning: first time iteration as shown below should be reference state')
            print('**')
            print('dx_along_plume = ', dx_along_plume)
            print('')
    print('loop over {:d} MNH files in {:s}'.format(len(mesonh_filenames), dir_data))
    indent = '    '
    if flag_plt:
        time_scene_time_strand_tecplot = np.zeros([2,len(mesonh_filenames)])
   
    i_mnh_file = 0
    time_scene_arr = []
    for ifile,  mesonh_filename in enumerate(mesonh_filenames):
      
        #load mesonh data
        #-----------------
        time_scene, datetime_MNHref, FireScene, FireScene2D, FireSceneMeshExtra = MNH_tools.load_centered_MesoNHField(mesonh_filename,      \
                                                                                                               ff_sv=extraVar_info, \
                                                                                                               dxyz_factor= dxyz_factor,    \
                                                                                                               verbose=verbose,indent=indent,dirDiag=dirDiag)
        
        if ifile == 0: 
            if extraVar_info is not None:
                for ii_ in range(len(extraVar_info)):
                    if extraVar_info.dim[ii_] == 'na': 
                        if extraVar_info.name[ii_] in FireScene.dtype.names: 
                           extraVar_info.dim[ii_] = '3D'
                        elif extraVar_info.name[ii_] in FireScene2D.dtype.names:
                           extraVar_info.dim[ii_] = '2D'
                        else:
                            pdb.set_trace()
        
        if zlimit is not None: 
            #do a clip on altitude
            klim = np.where(FireScene.zc<=zlimit)[2].max()
            FireScene = FireScene[:,:,:klim+1]
            FireSceneMeshExtra['zhat'] = FireSceneMeshExtra['zhat'][:klim+3]

        print(indent,os.path.basename(mesonh_filename), end=' ') 
        if 'SFObs' in FireScene.dtype.names:  
            print (' | max SVT001 {:.3e}'.format(FireScene.SFObs.max()), end=' ' )
        
        if time_scene < time_beg:
            print ('skip', end='\n')
            continue
        time_scene = float(time_scene) - time_ref
        print(' | t={:6.1f} '.format(time_scene), end= ' ')      
        datetime_ref =datetime_MNHref + datetime.timedelta(0,time_ref) 

        #print(indent,os.path.basename(mesonh_filename), end=' ') 
        if time_scene > time_end: 
            print ('')
            print (' ** stop here time_scene > time_end=',time_end )
            print ('time_scene = ', time_scene)
            print ('with time_ref = ', time_ref)
            break        

        if i_mnh_file > 0: 
            if time_scene == time_scene_arr[-1]: 
                print(' restart file, skip')
                continue

        time_scene_arr.append(time_scene)
        #set y axis origin at the center
        #extent_y = (FireScene.yc+.5*FireScene.dy).max() -  (FireScene.yc-.5*FireScene.dy).min()
        #FireScene.yc   = FireScene.yc   - .5* extent_y 
        #FireScene2D.yc = FireScene2D.yc - .5* extent_y 
        
        #print(time_scene, end=' ')      
        print(' | dimension: {:d}x{:d}x{:d}'.format(*list(FireScene.shape)) , end='\n')      
        
        if flag_plt:
            time_scene_time_strand_tecplot[0,i_mnh_file] = time_scene
            time_scene_time_strand_tecplot[1,i_mnh_file] = writePlt.global_var.timestr_id


        if i_mnh_file == 0: 
            
            # save vertical mesh
            #assume homogeneous mesh here over horizontal
            f = open(outputDir+'vertical_mesh_edge_model{:1d}s.txt'.format(modelN),'w')
            temp_ = []
            for z_ in FireSceneMeshExtra['zhat']: # assume homogeneous mesh 
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
            
            #if nextra >0:
            #    if 'SFObs' in extraVar_info.name : var_to_keep.append('plume_mask')

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
                   if (i_mnh_file == 0) | (nsv_ff == 0): # no plume in ref state
                        QQ[:,:,:,ivar] = np.zeros([ngx, ngy, ngz])
                   else:
                        QQ[:,:,:,ivar] = MNH_tools.plume_hull(dx_along_plume,FireScene, FireScene_ref, \
                                                                'SFObs', 0.02 , ref_val=1.e-3,flag='use xc,yc,zc keys')
                else:
                    print('issue in var selection, var =', var)
                    pdb.set_trace()
            
            writePlt.teclib.save_flow(expr_name_save,len(expr_name_save),flag_plot_lambda2,        \
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

            writePlt.teclib.save_surface(expr_name_save,len(expr_name_save),        \
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
            MNH_tools.dump_netcdf(outputDir_nc,expr_name_save,modelN,time_scene,datetime_ref,FireScene,FireScene2D,FireSceneMeshExtra,\
                                  ff_sv=extraVar_info,\
                                  flag_write=flag_write)
   
        #end loop
        i_mnh_file+=1

    if flag_plt:
        np.save(outputDir+'time_plt_strand',time_scene_time_strand_tecplot)
     
