import numpy as np
import sys
import os
import argparse
from matplotlib import cm 
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
from matplotlib.patches import Polygon
#homebrewed
import MNH_tools
import shutil
import pdb
from matplotlib.ticker import FormatStrFormatter

def get_qv_prop():
    fig =plt.figure(figsize=(12,6))
    ax = plt.subplot(121)
    m = Basemap(llcrnrlon=lon_min,
                llcrnrlat=lat_min,
                urcrnrlon=lon_max,
                urcrnrlat=lat_max,
                lat_0=(lat_max - lat_min)/2,
                lon_0=(lon_max-lon_min)/2,
                projection='merc',
                resolution = 'l',
                area_thresh=10000.,
                )
    try: 
        m.drawcoastlines()
    except: 
        pass
    m.drawcountries()
   
    im=m.pcolor(lons,lats, topo, cmap=cm.terrain)
    #im=m.imshow(topo.T, origin='lower',extent=extent, cmap=cm.terrain)
    qv=m.quiver(lons[::skipW,::skipW].flatten(),lats[::skipW,::skipW].flatten(),u[it,0,::skipW,::skipW].flatten(),v[it,0,::skipW,::skipW].flatten(),)# scale=1.e-5)
   
    fig.savefig('mm.png')
    scale = qv.scale
    width = qv.width
    plt.close(fig)
    os.remove('mm.png')

    return scale, width


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='to plot nc data from real case simulation')
    parser.add_argument('-i','--input', help='path to input nc file',required=True)
    parser.add_argument('-o','--output', help='base name of the pngs',required=True)
    parser.add_argument('-ll','--latlon', help='latlon of point to show on map',required=False)
    parser.add_argument('-bb','--domainbox', help='relative path to domain box to overlay',required=False)
    parser.add_argument('-zoom','--zoom', help='True to plot 0.05 degre around latlon',required=False)
    args = parser.parse_args()

    '''
    for shabeni1
    run do_plot_nc.py -i /mnt/dataEstrella/MNH_Real/shabeni/SHABE_model2.nc -o shabeniM2 -ll [-25.11596,31.23530]
    
    for skukuza4
    run do_plot_nc.py -i /mnt/dataEstrella/MNH_Real/skukuza/SKUKU_model1.nc -o skukuza4M1 -ll [-25.09969,31.44104] -bb SKUKU_model2.nc
t=0000
    '''
    name_in = args.input
    name_out = args.output
  
    if 'model1' in name_in: 
        ratioCorrection = .6
    if 'model2' in name_in: 
        ratioCorrection = .5
    if 'model3' in name_in: 
        ratioCorrection = .5

    flag_plotpt = False
    if args.latlon is not None:
        latlonpt = [ float(xx) for xx in args.latlon.split(',') ] 
        flag_plotpt = True
   
    print('latlon = ', latlonpt)

    flag_zoom = False
    if args.zoom is not None:
        flag_zoom =  args.zoom


    nc = Dataset(name_in,'r')
    
    dir_in = os.path.dirname(name_in)
    MNH_tools.ensure_dir(dir_in+'/png/')
    dir_out = dir_in+'/png/{:s}/'.format(name_out)
    MNH_tools.ensure_dir(dir_out)

    
    topo = nc.variables['topography'][:,:].T
    u =  np.transpose(nc.variables['u' ][:,:,:,:],[0,1,3,2])
    v =  np.transpose(nc.variables['v' ][:,:,:,:],[0,1,3,2])
    RH = np.transpose(nc.variables['RH'][:,:,:,:],[0,1,3,2])
    surfTheta = np.transpose(nc.variables['surfTheta'][:,:,:],[0,2,1])


    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    timenc = nc.variables['time'][:]

    yy,xx = np.meshgrid(x,y)
    extent=(xx.min(),xx.max(),yy.min(),yy.max())

    lat = nc.variables['lat'][:,:]
    lon = nc.variables['lon'][:,:]
    

    #refdatetime = datetime.datetime(2000,1,1,0,0,0)
    refdatetime = datetime.datetime.strptime(nc.variables['time'].units.split('ion:')[1].strip(), '%Y-%m-%d %H:%M:%S.%f') 
    skipW = 10

    flag_plotbox = False
    if args.domainbox is not None:
        flag_plotbox = True
        bb_list = [] 
        for name_box_ in args.domainbox.split(','):
            name_box = os.path.dirname(name_in) +'/'+ name_box_
            ncbox = Dataset(name_box,'r')
            bb = np.zeros([4,2])
            bb[0,:] = ncbox['lon'][:,:].T[0,0], ncbox['lat'][:,:].T[0,0]
            bb[1,:] = ncbox['lon'][:,:].T[-1,0], ncbox['lat'][:,:].T[-1,0]
            bb[2,:] = ncbox['lon'][:,:].T[-1,-1], ncbox['lat'][:,:].T[-1,-1]
            bb[3,:] = ncbox['lon'][:,:].T[0,-1], ncbox['lat'][:,:].T[0,-1]
            bb_list.append(bb)
            

    if flag_zoom:
        widthPlot =[ 0.015, 0.025]
        ilat1 = np.abs(lat[:,0]-(latlonpt[0]-widthPlot[0])).argmin()
        ilat2 = np.abs(lat[:,0]-(latlonpt[0]+widthPlot[0])).argmin()
        jlon1 = np.abs(lon[0,:]-(latlonpt[1]-widthPlot[1])).argmin()
        jlon2 = np.abs(lon[0,:]-(latlonpt[1]+widthPlot[1])).argmin()

    
        lat = lat[ilat1:ilat2,jlon1:jlon2]
        lon = lon[ilat1:ilat2,jlon1:jlon2]
        topo = topo[ilat1:ilat2,jlon1:jlon2]
        surfTheta = surfTheta[:,ilat1:ilat2,jlon1:jlon2]
        u         = u[:,:,ilat1:ilat2,jlon1:jlon2] 
        v         = v[:,:,ilat1:ilat2,jlon1:jlon2] 
        RH         = RH[:,:,ilat1:ilat2,jlon1:jlon2] 

        skipW = 3
    

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
                    resolution = 'f',
                    area_thresh=10000.,
                    )
        try: 
            m.drawcoastlines()
        except: 
            pass
        m.drawcountries()
        
        if it == 0:
            # convert lat and lon to map projection coordinates
            lons, lats = m(lon, lat)
            lonpt_, latpt_ = m(*latlonpt[::-1])

            if flag_plotpt: 
                idxpt = np.unravel_index( ( (lats-latpt_)**2 + (lons-lonpt_)**2 ).argmin(), lats.shape)

            scale, width = get_qv_prop()

        #scale = 1000
        #width = 0.01

        im=m.pcolor(lons,lats, topo.T, cmap=cm.terrain)
        qv=m.quiver(lons[::skipW,::skipW].flatten(),lats[::skipW,::skipW].flatten(),u[it,0,::skipW,::skipW].flatten(),v[it,0,::skipW,::skipW].flatten(), scale=2*scale, width=width)
        ax.quiverkey(qv, 0.9, 1.02, 10, '10 m/s', coordinates='axes')

        if flag_plotpt: 
            ax.scatter(lons[idxpt],lats[idxpt],c='r',s=5)
            qv2=m.quiver([lons[idxpt]],[lats[idxpt]],[u[it,0,idxpt[0],idxpt[1]]],[v[it,0,idxpt[0],idxpt[1]]],color='r',scale=2*scale, width=width)
            #print(u[it,0,idxpt[0],idxpt[1]],v[it,0,idxpt[0],idxpt[1]]) 

        ax.set_title('wind + topography')

        ax.set_xticks(lons[0,::2*skipW])
        ax.set_xlabel('x (m)')
        ax.set_yticks(lats[::2*skipW,0])
        ax.set_ylabel('y (m)')
        
        divider = make_axes_locatable(ax)
        cbaxes = divider.append_axes("bottom", size="5%", pad=0.5)
        cbar = fig.colorbar(im,orientation='horizontal',cax = cbaxes)
        cbar.set_label('topography (m)')
        
        if flag_plotbox:
            for bb in bb_list:
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
                    resolution = 'f',
                    area_thresh=10000.,
                    )
        try: 
            m.drawcoastlines()
        except: 
            pass
        m.drawcountries()

        #im=m.imshow(surfTheta[it,:,:].T, origin='lower',extent=extent)
        #im=m.pcolor(lons,lats, surfTheta[it,:,:])
        #ax.set_title('surface theta')
        im=m.pcolor(lons,lats, RH[it,0,:,:].T, vmin=20,vmax=60)
        ax.set_title('Relative Humidity')
        
        ax.set_xticks(lons[0,::2*skipW])
        ax.set_xticklabels(['{:.2f}'.format(xx) for xx in  lon[0,::2*skipW]])
        ax.set_xlabel('longitude')

        ax.set_yticks(lats[::2*skipW,0])
        ax.set_yticklabels(['{:.2f}'.format(xx) for xx in  lat[::2*skipW,0]])
        ax.set_ylabel('latitude')

        divider = make_axes_locatable(ax)
        cbaxes = divider.append_axes("bottom", size="5%", pad=0.5)
        cbar = fig.colorbar(im,orientation='horizontal',cax = cbaxes)
        #cbar.set_label('theta (K)')
        cbar.set_label('RH (%)')
        
        if flag_plotpt: 
            ax.scatter(lons[idxpt],lats[idxpt],c='r',s=5)

        datetime_ = refdatetime + datetime.timedelta(0,timenc[it])
        fig.suptitle('time = {:s}'.format(datetime_.strftime('%Y-%m-%d %H:%M')) )

        print( 't={:s}'.format(datetime_.strftime('%H%M')) )
        plt.savefig('{:s}/MNHReal_{:s}_{:s}.png'.format(dir_out, name_out, datetime_.strftime('%Y-%m-%d-%H%M')))    
        plt.close(fig)


