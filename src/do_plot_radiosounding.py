import numpy as np
import sys
import os
import argparse
import matplotlib as mpl
from matplotlib import cm 
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
from matplotlib.patches import Polygon
import shutil
import pdb

#####################################################
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
    
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='to plot nc data from real case simulation')
    parser.add_argument('-i','--input', help='path to input nc file',required=True)
    parser.add_argument('-o','--output', help='name of the radiosounding',required=True)
    parser.add_argument('-t','--time', help='time in seconds in the reference of the input nc',required=True)
    parser.add_argument('-ll','--latlon', help='latlon of point to show on map',required=True)
    args = parser.parse_args()

    '''
    for shabeni1
    run do_plot_nc.py -i /mnt/dataEstrella/MNH_Real/shabeni/SHABE_model2.nc -o shabeniM2_rsou -t 11:00:00 -ll [-25.11596,31.23530]
    '''
    name_in = args.input
    name_out = args.output
    timersou = np.array(np.array(args.time.split(':'),dtype=float)*np.array([3600,60,1])).sum()
    latlonpt = [ float(xx) for xx in args.latlon[1:-1].split(',') ] 
    
    nc = Dataset(name_in,'r')
    
    dir_in = os.path.dirname(name_in)
    dir_out = dir_in+'/rsou_png/'
    ensure_dir(dir_out)

    topo = nc.variables['topography'][:,:].T
    u = np.transpose(nc.variables['u'][:,:,:,:],[0,1,3,2])
    v = np.transpose(nc.variables['v'][:,:,:,:],[0,1,3,2])
    theta = np.transpose(nc.variables['Theta'][:,:,:,:],[0,1,3,2])
    rh = np.transpose(nc.variables['RH'][:,:,:,:],[0,1,3,2])

    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    timenc = nc.variables['time'][:]
    itrsou = np.abs(timenc.data-timersou).argmin()
    print(itrsou)
    yy,xx = np.meshgrid(x,y)
    extent=(xx.min(),xx.max(),yy.min(),yy.max())

    lats = nc.variables['lat'][:,:].T
    lons = nc.variables['lon'][:,:].T
    latpt_, lonpt_ = latlonpt
    idxpt = np.unravel_index( ( (lats-latpt_)**2 + (lons-lonpt_)**2 ).argmin(), lats.shape)
   

    refdatetime = datetime.datetime(2000,1,1,0,0,0)

    rsou = np.array([0]*nc['z'][:].size, dtype=np.dtype([('z_ground',float),('theta',float),('rh',float),('u',float),('v',float)])) 
    rsou = rsou.view(np.recarray)

    rsou.z_ground = nc['z'][:]*1.e-3 
    rsou.u = u[itrsou,:,idxpt[0],idxpt[1]]
    rsou.v = v[itrsou,:,idxpt[0],idxpt[1]]
    rsou.rh = rh[itrsou,:,idxpt[0],idxpt[1]]
    rsou.theta = theta[itrsou,:,idxpt[0],idxpt[1]]

    mpl.rcdefaults()
    mpl.rcParams['figure.subplot.left'] = .14
    mpl.rcParams['figure.subplot.right'] = .95
    mpl.rcParams['figure.subplot.top'] = .8
    mpl.rcParams['figure.subplot.bottom'] = .1
    mpl.rcParams['figure.subplot.hspace'] = 0.02
    mpl.rcParams['figure.subplot.wspace'] = 0.02

    fig = plt.figure(figsize=(4,8))
    ax1 = plt.subplot(111)
    ax1.plot(rsou.theta, rsou.z_ground, c='k')
    ax1.set_xlabel('potential temperature (K)')
    ax1.set_ylabel('altitude from ground level (km)')
    
    ax2=ax1.twiny()
    ax2.spines['top'].set_position(("axes", 1))
    ax2.spines['top'].set_visible(True)
    ax2.spines['top'].set_color('r')
    ax2.xaxis.set_label_position('top')
    ax2.xaxis.set_ticks_position('top')
    ax2.yaxis.label.set_color('r')
    ax2.tick_params(axis='x', colors='r',)
    ax2.set_xlabel('relative humidity (%)')
    ax2.xaxis.label.set_color('r')
    
    ax2.plot(rsou.rh, rsou.z_ground, c='r', ls=':')
    
    ax3=ax1.twiny()
    ax3.spines['top'].set_position(("axes", 1.12))
    ax3.spines['top'].set_visible(True)
    ax3.spines['top'].set_color('b')
    ax3.xaxis.set_label_position('top')
    ax3.xaxis.set_ticks_position('top')
    ax3.yaxis.label.set_color('b')
    ax3.tick_params(axis='x', colors='b',)
    ax3.set_xlabel('wind speed (m/s)')
    ax3.xaxis.label.set_color('b')

    
    ax3.plot(np.sqrt(rsou.u**2+rsou.v**2), rsou.z_ground, c='b', ls='--')
    
    ax1.set_ylim(0,6)
    ax1.set_title('shabeni1 22-08-2014 '+args.time, pad=13)

    ax1.set_zorder(1)
    ax2.set_zorder(2)
    ax3.set_zorder(3)

    fig.savefig(dir_out+name_out+'.png')
    plt.close(fig)
    sys.exit()








    for it in range(timenc.shape[0]):

        # determine range to print based on min, max lat and lon of the data
        margin = 0 # buffer to add to the range
        lat_min = np.min(lat) - margin
        lat_max = np.max(lat) + margin
        lon_min = np.min(lon) - margin
        lon_max = np.max(lon) + margin
        extent=(lon_min, lon_max, lat_min, lat_max)

        fig =plt.figure(figsize=(12., ratioCorrection*12.*topo.shape[1]/topo.shape[0]))

        #first plot
        ax = plt.subplot(121)
        m = Basemap(llcrnrlon=lon_min,
                    llcrnrlat=lat_min,
                    urcrnrlon=lon_max,
                    urcrnrlat=lat_max,
                    lat_0=(lat_max - lat_min)/2,
                    lon_0=(lon_max-lon_min)/2,
                    projection='merc',
                    resolution = 'h',
                    area_thresh=10000.,
                    )
        try: 
            m.drawcoastlines()
        except: 
            pass
        m.drawcountries()
        
        if it == 0:
            # convert lat and lon to map projection coordinates
            lons, lats = m(lon.T, lat.T)
            lonpt_, latpt_ = m(*latlonpt[::-1])

            if flag_plotpt: 
                idxpt = np.unravel_index( ( (lats-latpt_)**2 + (lons-lonpt_)**2 ).argmin(), lats.shape)

            scale, width = get_qv_prop()


        im=m.imshow(topo.T, origin='lower',extent=extent, cmap=cm.terrain)
        qv=m.quiver(lons[::skipW,::skipW].flatten(),lats[::skipW,::skipW].flatten(),u[it,0,::skipW,::skipW].flatten(),v[it,0,::skipW,::skipW].flatten(), scale=scale, width=width)
        ax.quiverkey(qv, 0.9, 1.02, 5, '5 m/s', coordinates='axes')
        
        if flag_plotpt: 
            ax.scatter(lons[idxpt],lats[idxpt],c='r',s=5)
            qv2=m.quiver([lons[idxpt]],[lats[idxpt]],[u[it,0,idxpt[0],idxpt[1]]],[v[it,0,idxpt[0],idxpt[1]]],color='r',scale=scale, width=width)
            print(u[it,0,idxpt[0],idxpt[1]],v[it,0,idxpt[0],idxpt[1]]) 

        ax.set_title('wind + topography')

        divider = make_axes_locatable(ax)
        cbaxes = divider.append_axes("bottom", size="5%", pad=0.05)
        cbar = fig.colorbar(im,orientation='horizontal',cax = cbaxes)
        cbar.set_label('topography (m)')
        
        if flag_plotbox:
            bbx1,bby1 = m(bb[0,0],bb[0,1])
            bbx2,bby2 = m(bb[1,0],bb[1,1])
            bbx3,bby3 = m(bb[2,0],bb[2,1])
            bbx4,bby4 = m(bb[3,0],bb[3,1])
            poly = Polygon([(bbx1,bby1),(bbx2,bby2),(bbx3,bby3),(bbx4,bby4)],facecolor='none',edgecolor='k',linewidth=1, linestyle=':')
            ax.add_patch(poly)
            
        #seconf plot
        ax = plt.subplot(122)
        m = Basemap(llcrnrlon=lon_min,
                    llcrnrlat=lat_min,
                    urcrnrlon=lon_max,
                    urcrnrlat=lat_max,
                    lat_0=(lat_max - lat_min)/2,
                    lon_0=(lon_max-lon_min)/2,
                    projection='merc',
                    resolution = 'h',
                    area_thresh=10000.,
                    )
        try: 
            m.drawcoastlines()
        except: 
            pass
        m.drawcountries()

        im=m.imshow(surfTheta[it,:,:].T, origin='lower',extent=extent)
        ax.set_title('surface theta')

        divider = make_axes_locatable(ax)
        cbaxes = divider.append_axes("bottom", size="5%", pad=0.05)
        cbar = fig.colorbar(im,orientation='horizontal',cax = cbaxes)
        cbar.set_label('theta (K)')

        datetime_ = refdatetime + datetime.timedelta(0,timenc[it])
        fig.suptitle('time = {:s}'.format(datetime_.strftime('%H:%M')) )

        print( 't={:s}'.format(datetime_.strftime('%H%M')) )
        plt.savefig('{:s}/MNHReal_{:s}_{:s}.png'.format(dir_out, name_out, datetime_.strftime('%H%M')))    
        plt.close(fig)


