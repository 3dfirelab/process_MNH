from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
import sys
import numpy as np
#import vtk
import pdb 
#import matplotlib.pyplot as plt
#from matplotlib import cm
from netCDF4 import Dataset
import scipy.interpolate
import os
import subprocess
from scipy import interpolate
import imp
#from vtk import *
#from vtk.util import numpy_support as VN
import glob 
import math 
import socket 
import itertools
from functools import reduce
import datetime
from scipy import ndimage

path_processMNH = os.environ['PATH_SRC_PYTHON_LOCAL']+'/process_MNH/src/'

sys.path.append(path_processMNH)

#homebrewed
sys.path.append(path_processMNH+'computeLambda2/')
import postproc
#sys.path.append(path_processMNH+'3dInterp/')
#import interp
    
#constant 
c_p = 1005.   # J  / kg / K
R   =  287.05 #  J / kg / K
kappa = old_div(R, c_p)
pressure_ref = 1.e5 # Pa


#####################################################
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


#####################################################################
def extrapolate_Nan(arr):

    # Mask of NaNs
    mask = np.isnan(arr)

    # Get indices of nearest non-NaN values
    idx = ndimage.distance_transform_edt(mask,
                                         return_distances=False,
                                         return_indices=True)

    # Replace NaNs with nearest neighbor values
    arr_filled = arr[tuple(idx)]

    return arr_filled


#####################################################################
def read_mesoNH(filename,extraVar=None, dirDiag='../08_diag/'):

    #from MNH: ch_ini_orilam.f90
    XBOLTZ = 1.380658E-23
    XAVOGADRO = 6.0221367E+23
    #XG = 9.80665
    #XP00 = 1.E5
    XMD = 28.9644E-3
    XRD = old_div(XAVOGADRO * XBOLTZ, XMD)
    XCPD = 7.* XRD /2.


    #open file
    ncfile = Dataset(filename,'r')
    if '.OUT.' in filename : 
        tmp = filename.replace('.OUT.','.').split('.')
        ncfile_grid_filename = '.'.join(tmp[:-2])+'.001.nc'
        ncfile_grid = Dataset(ncfile_grid_filename,'r')
        nctime = ncfile.variables['time']
    else:
        ncfile_grid_filename = filename
        ncfile_grid = ncfile 
        nctime = ncfile.variables['DTCUR']
    
    datetime_MNHref_str = nctime.units.split('since')[1].strip().split(' ')[0]
    datetime_MNHref = datetime.datetime.strptime(datetime_MNHref_str, '%Y-%m-%d')

    try:
        nclat  = ncfile_grid.variables['LAT']
        nclon  = ncfile_grid.variables['LON']
    except: 
        nclat = None
        nclon = None 

    ncx  = ncfile_grid.variables['XHAT']
    ncy  = ncfile_grid.variables['YHAT']
    ncz  = ncfile_grid.variables['ZHAT']
    ncrhodrefz = ncfile_grid.variables['RHOREFZ']
    
    nczs  = ncfile_grid.variables['ZS']

    ncTheta      = ncfile.variables['THT']
    ncTke      = ncfile.variables['TKET']
    ncpressure   = ncfile.variables['PABST']
    ncu   = ncfile.variables['UT']
    ncv   = ncfile.variables['VT']
    ncw   = ncfile.variables['WT']
    ncrvap   = ncfile.variables['RVT']
    try:
        ncrcloud = ncfile.variables['RCT']
        ncrrain = ncfile.variables['RRT']
        ncrice = ncfile.variables['RIT']
    except: 
        ncrcloud = np.zeros_like(ncrvap) 
        ncrrain  = np.zeros_like(ncrvap) 
        ncrice   = np.zeros_like(ncrvap) 
    try:
        ncthw        =  ncfile.variables['THW_FLX']
    except: 
        ncthw = np.zeros_like(ncTheta) 


    #0D var
    time = np.round(nctime[:],2)
    
    #1D var
    x_hat = ncx[:]
    y_hat = ncy[:]
    z_hat = ncz[:]
    
    rhodrefz = ncrhodrefz[:]
    
    #2D var
    zs = nczs[:,:].T

    #3D var
    #
    #if ncfile.variables['MASDEV'][0].data==53:
    #    #dynamics
    #    theta_raw    = ncTheta[:,:,:].T
    #    tke_raw      = ncTke[:,:,:].T
    #    thw_raw      = ncthw[:,:,:].T
    #    p_raw        = ncpressure[:,:,:].T
    #    u_raw        = ncu[:,:,:].T
    #    v_raw        = ncv[:,:,:].T
    #    w_raw        = ncw[:,:,:].T
    #    #microphysics
    #    rvap_raw      = ncrvap[:,:,:].T
    #    rcloud_raw    = ncrcloud[:,:,:].T
    #    rrain_raw     = ncrrain[:,:,:].T
    #    rice_raw      = ncrice[:,:,:].T
    # 
    #else: 
    #dynamics
    theta_raw    = ncTheta[0,:,:,:].T
    tke_raw      = ncTke[0,:,:,:].T
    thw_raw      = ncthw[0,:,:,:].T
    p_raw        = ncpressure[0,:,:,:].T
    u_raw        = ncu[0,:,:,:].T
    v_raw        = ncv[0,:,:,:].T
    w_raw        = ncw[0,:,:,:].T
    #microphysics
    rvap_raw      = ncrvap[0,:,:,:].T
    rcloud_raw    = ncrcloud[0,:,:,:].T
    rrain_raw     = ncrrain[0,:,:,:].T
    rice_raw      = ncrice[0,:,:,:].T


    nx= ncx.shape[0]-2
    ny= ncy.shape[0]-2
    nz= ncz.shape[0]-2
   
    #edge
    xe = x_hat[1:] 
    ye = y_hat[1:] 
    z_top = z_hat[-1]
    zb = np.zeros([nx+2,ny+2,nz+2])
    for k in range(nz+2):
        zb[:,:,k] = (old_div(z_hat[k], z_top)) * (z_top - zs[:,:]) + zs[:,:]

    #centered
    xc = np.zeros(nx+1); xc[:-1] = .5*(x_hat[1:-1]+x_hat[2:]) ; xc[-1] = xe[-1] + .5*(xe[-1]-xe[-2])
    yc = np.zeros(ny+1); yc[:-1] = .5*(y_hat[1:-1]+y_hat[2:]) ; yc[-1] = ye[-1] + .5*(ye[-1]-ye[-2])
    zc = np.zeros([nx+1,ny+1,nz+1])
    for k in range(nz):
        zc[:,:,k] = zb[1:,1:,k+1] + .5*(zb[1:,1:,k+2]-zb[1:,1:,k+1])
    zc[:,:,-1] = zb[1:,1:,-1] + .5*(zb[1:,1:,-1]-zb[1:,1:,-2])

    #grid the data 
    dtype_here = [('xe',float),('ye',float),('zb',float), ('xc',float),('yc',float),('zc',float),('u',float),('v',float),('w',float),('theta',float),('temp',float),('rhod',float),('tke',float),('rvap',float),('rcloud',float),('rrain',float),('rice',float)]
    if extraVar is not None: 
        for name_,dtype_ in zip(extraVar.name,extraVar.name_dtype): # add ff passive tracer
            dtype_here.append((name_,dtype_))
    #print('extar var:')
    #print(extraVar.name) 
    #print('')
    #out = np.zeros([nx+1,ny+1,nz+1],dtype=np.dtype(dtype_here))   
    #out = out.view(np.recarray)
    try:
        out = np.empty((nx+1, ny+1, nz+1), dtype=np.dtype(dtype_here))
        for name in out.dtype.names:       
            out[name].fill(np.nan)
        out = out.view(np.recarray)
            
    except: 
        pdb.set_trace()

    out.ye, out.xe, out.zb = np.meshgrid(ye,xe,z_hat[1:]) # use fake vertical grid here 
    out.zb = zb[1:,1:,1:]                                 # just to fill the horizontal grid
    out.yc, out.xc, out.zc = np.meshgrid(yc,xc,z_hat[1:]) # use fake vertical grid here 
    out.zc = zc                                           # just to fill the horizontal grid
    
    if theta_raw.shape[0] == ncx.shape[0]:
        out.theta = theta_raw[1:,1:,1:]
        out.tke = tke_raw[1:,1:,1:]
        out.temp = theta_raw[1:,1:,1:] * (old_div(p_raw[1:,1:,1:],pressure_ref))**(kappa)  
        out.u = u_raw[1:,1:,1:]
        out.v = v_raw[1:,1:,1:]
        out.w = w_raw[1:,1:,1:]
        out.rvap    = rvap_raw[1:,1:,1:]
        out.rcloud = rcloud_raw[1:,1:,1:]
        out.rrain   = rrain_raw[1:,1:,1:]
        out.rice    = rice_raw[1:,1:,1:]
 
        rhodref = np.zeros_like(thw_raw)
        for i,j in list(itertools.product(list(range(rhodref.shape[0])),list(range(rhodref.shape[1])))):
            rhodref[i,j,:] = rhodrefz[:]
    
        out.rhod = old_div(p_raw[1:,1:,1:], (R * out.temp))
    else: 
        out.theta[:-1,:-1,:-1] = theta_raw[0:,0:,0:]
        out.tke[:-1,:-1,:-1] = tke_raw[0:,0:,0:]
        out.temp[:-1,:-1,:-1] = theta_raw[0:,0:,0:] * (old_div(p_raw[0:,0:,0:],pressure_ref))**(kappa)  
        out.u[:,:-1,:-1] = u_raw[0:,0:,0:]
        out.v[:-1,:,:-1] = v_raw[0:,0:,0:]
        out.w[:-1,:-1,:] = w_raw[0:,0:,0:]
        out.rvap[:-1,:-1,:-1]    = rvap_raw[0:,0:,0:]
        out.rcloud[:-1,:-1,:-1] = rcloud_raw[0:,0:,0:]
        out.rrain[:-1,:-1,:-1]   = rrain_raw[0:,0:,0:]
        out.rice[:-1,:-1,:-1]    = rice_raw[0:,0:,0:]
    
        rhodref = np.zeros_like(thw_raw)
        for i,j in list(itertools.product(list(range(rhodref.shape[0])),list(range(rhodref.shape[1])))):
            rhodref[i,j,:] = rhodrefz[1:-1]

        
        out.rhod[:-1,:-1,:-1] = old_div(p_raw[:,:,:], (R * out.temp[1:,1:,1:]))


    extraVar2D = [] #var to find in diag
    if extraVar is not None: 
        #deal with passive tracer
        for name_,name_mnh_ in zip(extraVar.name,extraVar.name_mnh):
            if name_mnh_ == 'lambda2':
                
                postproc.lambda2.compute(out.shape[1], out.shape[2], out.shape[0], out.xc, out.yc, out.zc, 
                                         extrapolate_Nan(out.u), 
                                         extrapolate_Nan(out.v), 
                                         extrapolate_Nan(out.w)) 
                out[name_]    = postproc.lambda2.values
                postproc.lambda2.clear()
            elif name_mnh_ in ['HBLTOP', 'UM10', 'VM10', 'FF10MAX', 'T2M_ISBA']:
                extraVar2D.append([name_,name_mnh_])   
                continue 
            else:
                try:
                    ncVar        = ncfile.variables[name_mnh_]
                    #if ncfile.variables['MASDEV'][0].data==53:
                    #    var_raw      = ncVar[:,:,:].T
                    #else:
                    if len(ncVar.shape)==4:
                        var_raw      = ncVar[0,:,:,:].T
                        out[name_]    = var_raw[1:,1:,1:]
                    elif len(ncVar.shape)==3:
                        var_raw      = ncVar[0,:,:].T
                        out[name_][:,:,0]    = var_raw[1:,1:]
                        out[name_][:,:,1:]    = -999 

                    #print(name_mnh_)
                    #if name_mnh_ == 'SVT002': 
                    #    print('mmmm', out[name_].max(), end=' ')
                except:
                    if (name_mnh_ == 'SVT002') | (name_mnh_ == 'SFObs')  | (name_mnh_ == 'SVFF001'):
                        print('no '+name_mnh_)
                        pass
                    else:
                        pdb.set_trace()

    #ground = np.zeros([nx+1,ny+1],dtype=np.dtype([('xc',float),('yc',float),('heatFlux',float),('orography',float)]))   
    ground = np.empty(
                (nx+1, ny+1),
                dtype=[('xc', float), ('yc', float), ('heatFlux', float), ('orography', float)]
            )
    for name in ground.dtype.names:
        ground[name].fill(np.nan)
    ground = ground.view(np.recarray)

    ground.yc, ground.xc = np.meshgrid(yc,xc)
    ground.orography = zs[1:,1:]
    
    if theta_raw.shape[0] == ncx.shape[0]:
        ground.heatFlux = XCPD * rhodref[1:,1:,1] * thw_raw[1:,1:,1] * 1.e-3 # kW/m2
    else:
        ground.heatFlux[:-1,:-1] = XCPD * rhodref[0:,0:,1] * thw_raw[0:,0:,1] * 1.e-3 # kW/m2


    if len(extraVar2D)>0: #load variable from 08_diag
     
        # Define the new field
        new_field_name = [xx[0] for xx in extraVar2D]
        new_field_dtype = [float for ii in range( len(new_field_name))]
        new_dtype = ground.dtype.descr + \
                    [(new_field_name_, new_field_dtype_) \
                                       for new_field_name_, new_field_dtype_ in zip(new_field_name,new_field_dtype)]

        # Create the new recarray
        new_ground = np.recarray(ground.shape, dtype=new_dtype)

        # Copy existing data
        for field in ground.dtype.names:
            new_ground[field] = ground[field]

        # Initialize the new field with HBLTOP from 08_diag
        try: 
            ncfileDIAG = Dataset(dirDiag+os.path.basename(filename).replace('.nc','.d.nc'),'r')
        except: 
            print('no DIAG file found in ', dirDiag)
            print(dirDiag+os.path.basename(filename).replace('.nc','.d.nc'))
            print('stop here')
            sys.exit()

        for field, fiedl_mnh in extraVar2D:
            try:
                if len(ncfileDIAG[fiedl_mnh].shape)==3:
                    new_ground[field] = (ncfileDIAG[fiedl_mnh][0,:,:].T)[1:,1:]
                elif len(ncfileDIAG[fiedl_mnh].shape)==2:
                    new_ground[field] = (ncfileDIAG[fiedl_mnh][:,:].T)[1:,1:]

            except: 
                pdb.set_trace()

        ground = new_ground
        
    #mm  = (1./chf_to_HRR_ratio) *\
    #                       c_p/R * p_raw[1:,1:,1] / theta_raw[1:,1:,1] * thw_raw[1:,1:,1] * 1.e-3 # kW/m2

    #plt.imshow((ground.heatFlux).T,origin='lower',interpolation='nearest'); plt.show()
    
    if nclat != None:
        FireSceneMeshExtra = {'zhat': z_hat, 'lat':np.array(nclat[:,:]), 'lon': np.array(nclon[:,:])}
    else: 
        FireSceneMeshExtra = {'zhat': z_hat,}

    # close file
    ncfile.close()
    if ncfile_grid_filename != filename:
        ncfile_grid.close()
    
    return  out, time, datetime_MNHref, ground, FireSceneMeshExtra


######################################################################
def interp_mesonh(MesoNH, ground, ff_sv=None, dxyz_factor=1, verbose=0, indent='    '):

    '''
    there isnt any interpolation at the end here, appart from the C-grid wind data set to the center.
    the grid is the one from the input MNH file
    '''

    yc, xc, zc = MesoNH.yc[:-1,:-1,:-1], MesoNH.xc[:-1,:-1,:-1], MesoNH.zc[:-1,:-1,:-1]

    #output 3D flow
    nx = xc.shape[0]; ny = xc.shape[1]; nz = xc.shape[2]
    dtype_here = [('xc',float),   ('yc',float),    ('zc',float),  \
                  ('dx',float),   ('dy',float),    ('dz',float),  \
                  ('theta',float),('temp',float), ('rhod',float),('tke',float),  \
                  ('u',float),    ('v',float),     ('w',float), \
                  ('rvap',float),('rcloud',float),('rrain',float),('rice',float) ,\
                  #diag variable below
                  ('rh',float)  ]  

    if ff_sv is not None:
        for name,dtype in zip(ff_sv.name,ff_sv.name_dtype): # add ff passive tracer
            
            if name in ground.dtype.names : 
                continue
            
            try: 
                dtype_here.append((name.decode(),dtype))
            except: 
                dtype_here.append((name,dtype))

    #out = np.zeros([nx,ny,nz],dtype=np.dtype(dtype_here))   
    try: 
        out = np.empty((nx, ny, nz), dtype=np.dtype(dtype_here))
        for name in out.dtype.names:
            out[name].fill(np.nan)
        out = out.view(np.recarray)
    
    except: 
        pdb.set_trace()
    out = out.view(np.recarray)
    out.yc, out.xc, out.zc = yc, xc, zc
    if verbose > 0: print(2*indent, '3D output field is ', xc.shape)

    out.dx = MesoNH.xe[1:,:-1,:-1] -  MesoNH.xe[:-1,:-1,:-1]
    out.dy = MesoNH.ye[:-1,1:,:-1] -  MesoNH.ye[:-1,:-1,:-1]
    out.dz = MesoNH.zb[:-1,:-1,1:] -  MesoNH.zb[:-1,:-1,:-1]

    #output 2D
    var2D = ['pressure','theta','rvap']
    for name_ in ground.dtype.names:
        var2D.append(name_)

    out2D = np.zeros([nx,ny],dtype=np.dtype([ (xx,float) for xx in var2D ] ))  
    #('xc',float),('yc',float),\
    #('orography',float),('heatFlux',float), ('pressure',float), ('theta',float), ('rvap', float)]))   
    out2D = out2D.view(np.recarray)
    out2D.yc, out2D.xc = yc[:,:,0], xc[:,:,0] 

    ##
    ##3D field
    ##
    coord_pts = np.dstack((MesoNH.xc[:-1,:-1,:-1].flatten(),MesoNH[:-1,:-1,:-1].yc.flatten(),MesoNH[:-1,:-1,:-1].zc.flatten()))
   
    #for u
    if verbose > 0: print(2*indent,'u ...')
    u_c = .5*(MesoNH.u[:-1,:-1,:-1] + MesoNH.u[1:,:-1,:-1])
    #out.u = interpolate.griddata(coord_pts.reshape(-1,3) , u_c.flatten() , (xc, yc, zc),fill_value=0, method='nearest' )
    out.u = u_c

    #for v
    if verbose > 0: print(2*indent,'v ...')
    v_c = .5*(MesoNH.v[:-1,:-1,:-1] + MesoNH.v[:-1,1:,:-1])
    #out.v = interpolate.griddata(coord_pts.reshape(-1,3) , v_c.flatten() , (xc, yc, zc),fill_value=0, method='nearest' )
    out.v = v_c

    #for w
    if verbose > 0: print(2*indent,'w ...')
    w_c = .5*(MesoNH.w[:-1,:-1,:-1] + MesoNH.w[:-1,:-1,1:])
    #out.w = interpolate.griddata(coord_pts.reshape(-1,3) , w_c.flatten() , (xc, yc, zc),fill_value=0, method='nearest' )
    out.w = w_c
    
    #for mass point data
    if ff_sv is None: 
        tmp_list_var_ =  ['theta','temp','rhod','tke','rvap','rcloud','rrain','rice']
    else: 
        var3D = list(ff_sv.name)
        if 'hbltop' in var3D:
            var3D.remove('hbltop')
        if 'um10' in var3D:
            var3D.remove('um10')
        if 'vm10' in var3D:
            var3D.remove('vm10')
        if 'gust10' in var3D:
            var3D.remove('gust10')
        if 't2m' in var3D:
            var3D.remove('t2m')
        
        try: 
            tmp_list_var_ =  ['theta','temp','rhod','tke','rvap','rcloud','rrain','rice'] + [xx.decode() for xx in var3D]
        except: 
            tmp_list_var_ =  ['theta','temp','rhod','tke','rvap','rcloud','rrain','rice'] + [xx for xx in var3D]

    for key in tmp_list_var_:
        if verbose > 0: print(2*indent,key,' ...')
        #coord_pts = np.dstack((MesoNH.xc.flatten(),MesoNH.yc.flatten(),MesoNH.zc.flatten()))
        #out[key] = interpolate.griddata(coord_pts.reshape(-1,3) , MesoNH[key].flatten() , (xc, yc, zc),fill_value=0, method='nearest' ) # hear nearest as centered
        try: 
            out[key] = MesoNH[key][:-1,:-1,:-1]
        except: 
            pdb.set_trace()
    #3D diagnotic variables
    rvap_     =  out['rvap']
    temp_     =  out['temp']
    theta_    =  out['theta']
    pressure_ = pressure_ref*(old_div(temp_,theta_))**(1./kappa)
    out.rh = relativeHumidity(rvap_,pressure_,temp_)

    ##
    ##2D field
    ##
    xin = ground.xc.flatten()
    yin = ground.yc.flatten()
    coord_pts = np.dstack((xin,yin))
 
    #run  interpolation
    #out2D.orography = interpolate.griddata(coord_pts.reshape(-1,2) , ground.orography.flatten(), \
    #                                       (out2D.xc, out2D.yc),fill_value=0, method='nearest' )
    #out2D.heatFlux  = interpolate.griddata(coord_pts.reshape(-1,2) , ground.heatFlux.flatten(),  \
    #                                       (out2D.xc, out2D.yc),fill_value=0, method='nearest' )
    for field in ground.dtype.names: 
        if (field == 'xc') | (field == 'yc'): continue
        out2D[field] = ground[field][:-1,:-1]
        #out2D[field]  = interpolate.griddata(coord_pts.reshape(-1,2) , ground[field].flatten(),  \
        #                                       (out2D.xc, out2D.yc),fill_value=0, method='nearest' )
    
    #out2D.orography = ground.orography[:-1,:-1]
    #out2D.heatFlux = ground.heatFlux[:-1,:-1]
    
    rvap_     =  MesoNH['rvap'][:,:,0]
    temp_     =  MesoNH['temp'][:,:,0]
    theta_    =  MesoNH['theta'][:,:,0]
    pressure_ = pressure_ref*(old_div(temp_,theta_))**(1./kappa)
    
    #out2D.pressure = interpolate.griddata(coord_pts.reshape(-1,2) , pressure_.flatten(), (out2D.xc, out2D.yc),fill_value=0, method='nearest' ) 
    out2D.pressure = pressure_[:-1,:-1]
    out2D.theta = MesoNH['theta'][:-1,:-1,0]
    out2D.rvap = MesoNH['rvap'][:-1,:-1,0]
   
    

    if dxyz_factor == 1:
        return out, out2D 

    else:

        dxyz_factor_x = list(factors(out.shape[0]))[::-1] [ np.abs( np.array(list(factors(out.shape[0])))[::-1] - dxyz_factor ).argmin() ]
        dxyz_factor_y = list(factors(out.shape[1]))[::-1] [ np.abs( np.array(list(factors(out.shape[1])))[::-1] - dxyz_factor ).argmin() ]
        dxyz_factor_z = list(factors(out.shape[2]))[::-1] [ np.abs( np.array(list(factors(out.shape[2])))[::-1] - dxyz_factor ).argmin() ]
        
        nxr,nyr,nzr = old_div(nx,dxyz_factor_x), old_div(ny,dxyz_factor_y), old_div(nz,dxyz_factor_z)

        #new reduced 3D fields
        out_r = np.zeros([nxr,nyr,nzr],dtype= out.dtype ) 
        out_r = out_r.view(np.recarray)
        for key in out.dtype.names:
            out_r[key] = shrink_average_3d(out[key], *out_r.shape)

        #new reduced 2D fields
        out2D_r = np.zeros([nxr,nyr],dtype=  out2D.dtype )
        out2D_r = out2D_r.view(np.recarray)
        for key in out2D.dtype.names:
            out2D_r[key] = shrink_average_2d(out2D[key], *out2D_r.shape)


        return out_r, out2D_r


#############################################################
def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))


#############################################################
def shrink_average_3d(data, rows, cols, levels):
    return (old_div((old_div((old_div(data.reshape(rows,  old_div(data.shape[0],rows),    \
                          cols,   old_div(data.shape[1],cols),    \
                          levels, old_div(data.shape[2],levels),  ).sum(axis=1),(old_div(data.shape[0],rows)))  \
                                                         ).sum(axis=2),(old_div(data.shape[1],cols)))  \
                                                         ).sum(axis=3),(old_div(data.shape[2],levels)))\
            )

#############################################################
def shrink_average_2d(data, rows, cols):
    return  (old_div((old_div(data.reshape(rows,  old_div(data.shape[0],rows),    \
                          cols,   old_div(data.shape[1],cols),    ).sum(axis=1),(old_div(data.shape[0],rows)))  \
                                                         ).sum(axis=2),(old_div(data.shape[1],cols)))  \
            )

#############################################################
def dump_netcdf(outputDir,filename_in,modelN,time,datetimeRef,Firescene,Firescene2D,FiresceneMeshExtra,ff_sv=None,flag_write='append'):

    def for_display_in_netcdf(arr):
        if len(arr.shape) == 3:
            return np.swapaxes(arr,0,2)
        if len(arr.shape) == 2:
            return np.swapaxes(arr,0,1)


    filename_out = outputDir+filename_in+"_model{:1d}.nc".format(modelN)

    if flag_write == 'init' :
        if os.path.isfile(filename_out):
            os.remove(filename_out)
        ncfile = Dataset(filename_out,'w')
        
        ncfile.description = 'MNH Output for simulation ' + filename_in + ' and model{:1d}'.format(modelN)
       
        # Global attributes
        setattr(ncfile, 'created', 'R. Paugam') 
        setattr(ncfile, 'company', 'DoWeNeedOne')
        setattr(ncfile, 'title', 'MNH Output')

        # dimensions
        ncfile.createDimension('x',Firescene.xc.shape[0])
        ncfile.createDimension('y',Firescene.yc.shape[1])
        ncfile.createDimension('z',Firescene.zc.shape[2])
        ncfile.createDimension('time',None)

        ncx = ncfile.createVariable('x', 'f4', ('x',))
        setattr(ncx, 'long_name', 'x')
        setattr(ncx, 'standard_name', 'x')
        setattr(ncx, 'units','m')

        ncy = ncfile.createVariable('y', 'f4', ('y',))
        setattr(ncy, 'long_name', 'y')
        setattr(ncy, 'standard_name', 'y')
        setattr(ncy, 'units','m')
        
        ncz = ncfile.createVariable('z', 'f4', ('z',))
        setattr(ncz, 'long_name', 'z')
        setattr(ncz, 'standard_name', 'z')
        setattr(ncz, 'units','m')
        
        ncTime = ncfile.createVariable('time', 'f8', ('time',))
        setattr(ncTime, 'long_name', 'time')
        setattr(ncTime, 'standard_name', 'time')
        setattr(ncTime, 'units','seconds since fire ignition: {:s}'.format(datetimeRef.strftime('%Y-%m-%d %H:%M:%S.%f') ))
        
        #mesh
        if 'lat' in FiresceneMeshExtra.keys():
            nclat = ncfile.createVariable('lat', 'f8', ('y','x'))
            setattr(nclat, 'long_name', 'latitude center cell')
            setattr(nclat, 'standard_name', 'lat')
            setattr(nclat, 'units','degree')
            
            nclon = ncfile.createVariable('lon', 'f8', ('y','x',))
            setattr(nclon, 'long_name', 'longitude center cell')
            setattr(nclon, 'standard_name', 'lon')
            setattr(nclon, 'units','degree')
        
        ncdx = ncfile.createVariable('dx', 'f8', ('z', 'y', 'x',))
        setattr(ncdx, 'long_name', 'dx')
        setattr(ncdx, 'standard_name', 'dx')
        setattr(ncdx, 'units','m')

        ncdy = ncfile.createVariable('dy', 'f8', ('z', 'y', 'x',))
        setattr(ncdy, 'long_name', 'dy')
        setattr(ncdy, 'standard_name', 'dy')
        setattr(ncdy, 'units','m')
        
        ncdz = ncfile.createVariable('dz', 'f8', ('z', 'y', 'x',))
        setattr(ncdz, 'long_name', 'dz')
        setattr(ncdz, 'standard_name', 'dz')
        setattr(ncdz, 'units','m')
       
        # 2D variables
        nctopo    = ncfile.createVariable('topography',    'f8', ('y', 'x',))
        setattr(nctopo, 'long_name', 'topography, altitude at z_hat=0') 
        setattr(nctopo, 'standard_name', 'topo') 
        setattr(nctopo, 'units', 'm') 
        
        ncHeatFlux    = ncfile.createVariable('heatFlux',    'f8', ('time','y', 'x',))
        setattr(ncHeatFlux, 'long_name', 'Sensible Heat Flux') 
        setattr(ncHeatFlux, 'standard_name', 'heatFlux') 
        setattr(ncHeatFlux, 'units', 'kW/m2') 
        
        ncp2D    = ncfile.createVariable('surfPressure',    'f8', ('time','y', 'x',))
        setattr(ncp2D, 'long_name', 'pressure, altitude at z_hat=0') 
        setattr(ncp2D, 'standard_name', 'surfPressure') 
        setattr(ncp2D, 'units', 'Pa') 
        
        ncth2D    = ncfile.createVariable('surfTheta',    'f8', ('time','y', 'x',))
        setattr(ncth2D, 'long_name', 'theta, altitude at z_hat=0') 
        setattr(ncth2D, 'standard_name', 'surfTheta') 
        setattr(ncth2D, 'units', 'K') 
        
        ncrv2D    = ncfile.createVariable('surfRvap',    'f8', ('time','y', 'x',))
        setattr(ncrv2D, 'long_name', 'rvap, altitude at z_hat=0') 
        setattr(ncrv2D, 'standard_name', 'surfRvap') 
        setattr(ncrv2D, 'units', 'kg/kg') 
        
        # 3D variables
        ncalt    = ncfile.createVariable('altitude',    'f4', ('z','y', 'x',))
        setattr(ncalt, 'long_name', 'altitude level at mesh center') 
        setattr(ncalt, 'standard_name', 'alt') 
        setattr(ncalt, 'units', 'm') 

        ncu    = ncfile.createVariable('u',    'f8', ('time','z', 'y', 'x',))
        setattr(ncu, 'long_name', 'wind field u') 
        setattr(ncu, 'standard_name', 'u') 
        setattr(ncu, 'units', 'm/s') 
        
        ncv    = ncfile.createVariable('v',    'f8', ('time','z', 'y', 'x',))
        setattr(ncv, 'long_name', 'wind field v') 
        setattr(ncv, 'standard_name', 'v') 
        setattr(ncv, 'units', 'm/s') 
        
        ncw    = ncfile.createVariable('w',    'f8', ('time','z', 'y', 'x',))
        setattr(ncw, 'long_name', 'wind field w') 
        setattr(ncw, 'standard_name', 'w') 
        setattr(ncw, 'units', 'm/s') 

        ncTemp    = ncfile.createVariable('Temp',    'f8', ('time','z', 'y', 'x',))
        setattr(ncTemp, 'long_name', 'Kinetic Temperature') 
        setattr(ncTemp, 'standard_name', 'Temp') 
        setattr(ncTemp, 'units', 'K') 
        
        ncTheta    = ncfile.createVariable('Theta',    'f8', ('time','z', 'y', 'x',))
        setattr(ncTheta, 'long_name', 'Potential Temperature') 
        setattr(ncTheta, 'standard_name', 'Theta') 
        setattr(ncTheta, 'units', 'K') 
        
        ncRho    = ncfile.createVariable('Rho',    'f8', ('time','z', 'y', 'x',))
        setattr(ncRho, 'long_name', 'Dry air density') 
        setattr(ncRho, 'standard_name', 'Rho') 
        setattr(ncRho, 'units', 'kg/m3') 
        
        ncTke    = ncfile.createVariable('Tke',    'f8', ('time','z', 'y', 'x',))
        setattr(ncTke, 'long_name', 'Turbulent Kinetic Enery') 
        setattr(ncTke, 'standard_name', 'tke') 
        setattr(ncTke, 'units', 'm2/s2') 
        
        ncRvap    = ncfile.createVariable('Rvap',    'f8', ('time','z', 'y', 'x',))
        setattr(ncRvap, 'long_name', 'Vapor Mixing ratio') 
        setattr(ncRvap, 'standard_name', 'rvap') 
        setattr(ncRvap, 'units', 'kg/kg') 
        
        ncRcloud    = ncfile.createVariable('Rcloud',    'f8', ('time','z', 'y', 'x',))
        setattr(ncRcloud, 'long_name', 'Cloud Mixing ratio') 
        setattr(ncRcloud, 'standard_name', 'rcloud') 
        setattr(ncRcloud, 'units', 'kg/kg')               
        
        ncRrain    = ncfile.createVariable('Rrain',    'f8', ('time','z', 'y', 'x',))
        setattr(ncRrain, 'long_name', 'Rain Mixing ratio') 
        setattr(ncRrain, 'standard_name', 'rrain') 
        setattr(ncRrain, 'units', 'kg/kg')               
        
        ncRice    = ncfile.createVariable('Rice',    'f8', ('time','z', 'y', 'x',))
        setattr(ncRice, 'long_name', 'Ice Mixing ratio') 
        setattr(ncRice, 'standard_name', 'rice') 
        setattr(ncRice, 'units', 'kg/kg')               
        
        ncRH    = ncfile.createVariable('RH',    'f8', ('time','z', 'y', 'x',))
        setattr(ncRH, 'long_name', 'Relative Humidity') 
        setattr(ncRH, 'standard_name', 'RH') 
        setattr(ncRH, 'units', '%') 
        
        ncTracer = []
        if ff_sv is not None: 
            for name,dimType in zip(ff_sv.name,ff_sv.dim):
                if dimType == '3D': 
                    ncTracer.append(ncfile.createVariable(name,    'f8', ('time','z', 'y', 'x',)) )
                    setattr(ncTracer[-1], 'long_name', name ) 
                    setattr(ncTracer[-1], 'standard_name', name) 
                    setattr(ncTracer[-1], 'units', '-') 
                elif dimType == '2D': 
                    ncTracer.append(ncfile.createVariable(name,    'f8', ('time','y', 'x',)) )
                    setattr(ncTracer[-1], 'long_name', name ) 
                    setattr(ncTracer[-1], 'standard_name', name) 
                    setattr(ncTracer[-1], 'units', '-') 
                else: 
                    pdb.set_trace()

        zhat = FiresceneMeshExtra['zhat']

        ncx[:]  = np.array(Firescene.xc,dtype=np.float32)[:,0,0] #for_display_in_netcdf(np.array(Firescene.xc,dtype=np.float32))
        ncy[:]  = np.array(Firescene.yc,dtype=np.float32)[0,:,0] #for_display_in_netcdf(np.array(Firescene.yc,dtype=np.float32))
        ncz[:]  = np.array(.5*(zhat[1:-1]+zhat[2:]),dtype=np.float32) #for_display_in_netcdf(np.array(Firescene.zc,dtype=np.float32))       
        nctopo[:,:] = for_display_in_netcdf(Firescene2D.orography)
        ncalt[:,:,:] = for_display_in_netcdf(Firescene.zc)
      
        if 'lat' in FiresceneMeshExtra.keys():
            nclat[:,:] = 0.5*(FiresceneMeshExtra['lat'][1:-1,1:-1] + FiresceneMeshExtra['lat'][2:,2:])
            nclon[:,:] = 0.5*(FiresceneMeshExtra['lon'][1:-1,1:-1] + FiresceneMeshExtra['lon'][2:,2:])

        ncdx[:,:,:]  = for_display_in_netcdf(Firescene.dx)
        ncdy[:,:,:]  = for_display_in_netcdf(Firescene.dy)
        ncdz[:,:,:]  = for_display_in_netcdf(Firescene.dz)       

        i_time = 0

    elif flag_write == 'append':
        ncfile = Dataset(filename_out,'a')
        #1D
        ncTime = ncfile.variables['time']
        #2D
        ncHeatFlux    = ncfile.variables['heatFlux']
        ncp2D         = ncfile.variables['surfPressure']
        ncth2D       = ncfile.variables['surfTheta']
        ncrv2D       = ncfile.variables['surfRvap']
        #3D
        ncu           = ncfile.variables['u']
        ncv           = ncfile.variables['v']
        ncw           = ncfile.variables['w']
        ncTheta       = ncfile.variables['Theta']
        ncTke         = ncfile.variables['Tke']
        ncRvap        = ncfile.variables['Rvap']
        ncRcloud      = ncfile.variables['Rcloud']
        ncRrain       = ncfile.variables['Rrain']
        ncRice        = ncfile.variables['Rice']
        ncRH          = ncfile.variables['RH']
        ncRho         = ncfile.variables['Rho']
        ncTemp        = ncfile.variables['Temp']
        
        #2D or 3D
        ncTracer      = []
        if ff_sv is not None: 
            for name,dimType in zip( ff_sv.name, ff_sv.dim):
                try: 
                    ncTracer.append(ncfile.variables[name.decode()])
                except: 
                    ncTracer.append(ncfile.variables[name])
        
        
        i_time = ncTime.shape[0]

    else:
        print('issue with flag_write_netcdf')
        pdb.set_trace()
    
    #1D
    ncTime[i_time]          = time
    #2D
    ncHeatFlux[i_time,:,:]         = for_display_in_netcdf(Firescene2D.heatFlux)
    ncp2D[i_time,:,:]              = for_display_in_netcdf(Firescene2D.pressure)
    ncth2D[i_time,:,:]              = for_display_in_netcdf(Firescene2D.theta)
    ncrv2D[i_time,:,:]              = for_display_in_netcdf(Firescene2D.rvap)
    #3D
    ncu[i_time,:,:,:]       = for_display_in_netcdf(Firescene.u[:,:,:])
    ncv[i_time,:,:,:]       = for_display_in_netcdf(Firescene.v[:,:,:])
    ncw[i_time,:,:,:]       = for_display_in_netcdf(Firescene.w[:,:,:])
    ncTemp[i_time,:,:,:]    = for_display_in_netcdf(Firescene.temp[:,:,:])
    ncTheta[i_time,:,:,:]   = for_display_in_netcdf(Firescene.theta[:,:,:])
    ncTke[i_time,:,:,:]     = for_display_in_netcdf(Firescene.tke[:,:,:])
    ncRvap[i_time,:,:,:]    = for_display_in_netcdf(Firescene.rvap[:,:,:])
    ncRcloud[i_time,:,:,:]    = for_display_in_netcdf(Firescene.rcloud[:,:,:])
    ncRrain[i_time,:,:,:]    = for_display_in_netcdf(Firescene.rrain[:,:,:])
    ncRice[i_time,:,:,:]    = for_display_in_netcdf(Firescene.rice[:,:,:])
    ncRH[i_time,:,:,:]      = for_display_in_netcdf(Firescene.rh[:,:,:])
    ncRho[i_time,:,:,:]     = for_display_in_netcdf(Firescene.rhod[:,:,:])
    if ff_sv is not None: 
        for i_tracer, (name, dimType) in enumerate(zip(ff_sv.name, ff_sv.dim)):
            if dimType == '3D': 
                try: 
                    ncTracer[i_tracer][i_time,:,:,:]  = for_display_in_netcdf(Firescene[name.decode()][:,:,:])
                except: 
                    ncTracer[i_tracer][i_time,:,:,:]  = for_display_in_netcdf(Firescene[name][:,:,:])

            elif dimType == '2D': 
                try: 
                    ncTracer[i_tracer][i_time,:,:]  = for_display_in_netcdf(Firescene2D[name.decode()][:,:])
                except: 
                    ncTracer[i_tracer][i_time,:,:]  = for_display_in_netcdf(Firescene2D[name][:,:])

            else:
                pdb.set_trace()

    ncfile.close()
   
    #and copy a version without hdf for paraview
    #filename_out_2 = outputDir+filename_in+"_nohdf.nc"
    #subprocess.call(["nccopy", "-k", "1", filename_out, filename_out_2 ])
    
    return 0

#######################################
def load_centered_MesoNHField(mesonh_filename,        \
                              ff_sv=None,             \
                              dxyz_factor=1,          \
                              verbose=0, indent='    ',\
                              dirDiag='../08_diag/'):

    #load raw MesoNH Data from NetCDF
    if verbose > 0 : print(indent, 'load raw data')
    MesoNH_raw, time_frame, datetime_MNHref, ground, FireSceneMeshExtra = read_mesoNH(mesonh_filename,extraVar=ff_sv,dirDiag=dirDiag)

    #interpolated on the mesh centered
    if verbose > 0 : print(indent, 'interpolate data')
    FireScene, FireScene2D = interp_mesonh(MesoNH_raw, ground, ff_sv=ff_sv, dxyz_factor=dxyz_factor, \
                                           verbose=verbose, indent=indent)

    return time_frame, datetime_MNHref, FireScene, FireScene2D, FireSceneMeshExtra


#############################
def relativeHumidity(y_v,pressure,temperature):
    '''
    in:
    yv: vapor mixing ratio kg/kg
    pressure Pa
    temperature K 
    
    out:relative humidity in %
    '''
    Wvap = 18.01
    Wair = 28.85
    # pressure at saturation form sonntag http://cires.colorado.edu/~voemel/vp.html
    p_vsat =   np.exp(old_div(-6096.9385,temperature) + 16.635794  - 2.711193e-2 * temperature  + 1.673952e-5 * temperature**2 + 2.433502 * np.log(temperature) ) * 100 
    x_vsat = old_div(p_vsat,pressure)
    y_vsat = x_vsat*Wvap / (x_vsat*Wvap+ (1-x_vsat)*Wair) # kg/kg
    return old_div(y_v,y_vsat) * 100.


################################################
def plume_hull(dx_plume, data,data_ref,key,threshold_percent,ref_val=1.e-3,flag='use x,y,z keys'):
    
    if flag=='use x,y,z keys':
        x_key = 'x'
        y_key = 'y'
        z_key = 'z'
    elif flag=='use xc,yc,zc keys':
        x_key = 'xc'
        y_key = 'yc'
        z_key = 'zc'
    else :
        print('bad key in plume_hull, flag =', flag)
        pdb.set_trace()

    hull = np.zeros(data.shape)
    pts_grid = np.array([data[x_key].flatten(), data[y_key].flatten(), data[z_key].flatten()]).transpose()

    #ref_val_arr = np.where( np.abs(data_ref[key]) > ref_val, np.abs(data_ref[key]) , np.zeros_like(data_ref[key])+ref_val)
    test = np.where( data_ref[key]>=0, np.abs(data[key]-data_ref[key]),  0)#/ ref_val_arr
    test = old_div(test,test.max())
    idx=np.where(np.round(test,2) >= threshold_percent)
        
    #y_center = data.shape[1]/2
    #plt.figure()
    #plt.imshow(np.ma.masked_where(test[:,y_center,:]<threshold_percent,test[:,y_center,:]).T,origin='lower',interpolation='nearest')
    
    #plt.figure()
    #ax = plt.subplot(111)
    #ax.imshow(data[key][:,y_center,:].T,origin='lower')
    ##ax = plt.subplot(122)
    ##ax.imshow(condition[:,y_center,:].T,origin='lower')
    #plt.show()
    #pdb.set_trace()
 
    if len(idx[0])==0: 
        return np.zeros(data.shape)

    #dx_plume = 200
    range_x_plume_min = data[x_key][idx].min()
    range_x_plume_max = data[x_key][idx].max()*1.01
    range_x_plume     = np.arange(range_x_plume_min,range_x_plume_max+dx_plume,dx_plume)
    
    for i_range, (xb,xe) in enumerate(zip(range_x_plume[:-1],range_x_plume[1:])):
       
        idx_sub=np.where( (test >= threshold_percent) & (data[x_key]>=xb) & (data[x_key]<xe))
        pts_plume = np.array([data[x_key][idx_sub], data[y_key][idx_sub], data[z_key][idx_sub]]).transpose()
       
        if len(idx_sub[0]) == 0: 
            continue

        ii = 0
        while (np.unique(pts_plume[:,0]).shape[0] == 1) | (len(idx_sub[0]) < 5):
            xb = range_x_plume[i_range-ii] # assume this does not happen at the start of the loop
            idx_sub=np.where( (test >= threshold_percent) & (data[x_key]>=xb) & (data[x_key]<xe))
            pts_plume = np.array([data[x_key][idx_sub], data[y_key][idx_sub], data[z_key][idx_sub]]).transpose()
            if ii == i_range: 
                return np.zeros(data.shape) 
                #pdb.set_trace()
            ii+=1

        try: 
            hull_sub_flat = in_hull(pts_grid,pts_plume )
            idx_ = np.where(hull_sub_flat.reshape(data.shape) == 1 )
            hull[idx_] = 1
        except :
            print('qhull: no hull between {:.1f} - {:.1f}'.format(xb,xe))

    #idx=np.where((hull==1)&(test >= threshold_percent))
    #hull[idx]=1

    return hull
   


#######################################
if __name__ == '__main__':
#######################################
    mesonh_filename  = '../data/FI20m.1.SEG01.005.nc4'
    time, FireScene, FireScene2D, FireSceneMeshExtra = load_centered_MesoNHField(mesonh_filename, dxyz_factor=2)

    '''
    fig = plt.figure()
    ax = plt.subplot(121)
    im = ax.imshow(FireScene.w[:,old_div(FireScene.shape[1],2),:].T,origin='lower',interpolation='nearest')
    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel('w (m/s)')

    ax = plt.subplot(122)
    im = ax.imshow(FireScene2D.heatFlux.T,origin='lower',interpolation='nearest')
    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel('heatFlux kW/m2')

    plt.show()
    '''
