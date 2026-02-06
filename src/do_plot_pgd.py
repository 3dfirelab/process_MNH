import numpy as np
import sys
import os
import glob
import argparse
from matplotlib import cm 
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
from matplotlib.patches import Polygon
from shapely.geometry import Polygon as PolygonShapely

import geojson 
import geopandas as gpd

#homebrewed
import MNH_tools
import shutil
import pdb
from matplotlib.ticker import FormatStrFormatter

###################################
if __name__ == '__main__':
###################################
    
    parser = argparse.ArgumentParser(description='to plot nc data from real case simulation')
    parser.add_argument('-i','--input', help='path to input 01_prep_pgd dir',required=True)
    parser.add_argument('-ll','--latlon', help='latlon of point to show on map',required=False)
    args = parser.parse_args()
    flag_plotpt = False
    if args.latlon is not None:
        latlonpt = [ float(xx) for xx in args.latlon[1:-1].split(',') ] 
        flag_plotpt = True

    dir_in = args.input

    ncfiles = glob.glob(dir_in+'/PGD*.nested.nc')
    reso = []
    for ncfile in ncfiles:
        reso.append(float(ncfile.split('_D')[1].split('mA')[0]))
        
    bb = np.zeros([len(ncfiles),4,2])

    for  inc, ncfile in enumerate(ncfiles):
        ncbox = Dataset(ncfile,'r')
        bb[inc,0,:] = ncbox['longitude'][:,:].T[0,0], ncbox['latitude'][:,:].T[0,0]
        bb[inc,1,:] = ncbox['longitude'][:,:].T[-1,0], ncbox['latitude'][:,:].T[-1,0]
        bb[inc,2,:] = ncbox['longitude'][:,:].T[-1,-1], ncbox['latitude'][:,:].T[-1,-1]
        bb[inc,3,:] = ncbox['longitude'][:,:].T[0,-1], ncbox['latitude'][:,:].T[0,-1]

    ii = np.array(reso).argmax()
    ncbox = Dataset(ncfiles[ii],'r')
    lat = ncbox['latitude'][:,:].T
    lon = ncbox['longitude'][:,:].T
    topo = ncbox.variables['ZS'][:,:].T

    # determine range to print based on min, max lat and lon of the data
    margin = 0 # buffer to add to the range
    lat_min = np.min(lat) - margin
    lat_max = np.max(lat) + margin
    lon_min = np.min(lon) - margin
    lon_max = np.max(lon) + margin
    extent=(lon_min, lon_max, lat_min, lat_max)

    ratioCorrection = .6
    fig =plt.figure(figsize=(8., 8))

    #first plot
    ax = plt.subplot(111)
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
    
    lons, lats = m(lon, lat)
    
    im = m.pcolor(lons,lats, topo, cmap=cm.terrain)
    
    # Example of your data
    polygons = []
    names = ["pgd1", "pgd2", "pgd3"]  # Names for the polygons

    for inc in range(len(ncfiles)):
        llbbx1,llbby1 = bb[inc,0,0],bb[inc,0,1]
        llbbx2,llbby2 = bb[inc,1,0],bb[inc,1,1]
        llbbx3,llbby3 = bb[inc,2,0],bb[inc,2,1]
        llbbx4,llbby4 = bb[inc,3,0],bb[inc,3,1]
            
        # Create polygon geometry
        coords = [(llbbx1, llbby1), (llbbx2, llbby2), (llbbx3, llbby3), (llbbx4, llbby4), (llbbx1, llbby1)]  # Closed loop
        polygons.append(PolygonShapely(coords))
        
        if inc == 0 : continue
        bbx1,bby1 = m(bb[inc,0,0],bb[inc,0,1])
        bbx2,bby2 = m(bb[inc,1,0],bb[inc,1,1])
        bbx3,bby3 = m(bb[inc,2,0],bb[inc,2,1])
        bbx4,bby4 = m(bb[inc,3,0],bb[inc,3,1])
        poly = Polygon([(bbx1,bby1),(bbx2,bby2),(bbx3,bby3),(bbx4,bby4)],
                        facecolor='none',edgecolor='k',linewidth=1, linestyle=':')
        ax.add_patch(poly)
        
    
    # Create GeoDataFrame with CRS set to EPSG:4326
    gdf = gpd.GeoDataFrame({"name": names[:len(ncfiles)], "geometry": polygons})
    gdf = gdf.set_geometry("geometry")
    gdf = gdf.set_crs("EPSG:4326")

    # Export to GeoJSON
    gdf.to_file("pdfFootPrints.geojson", driver="GeoJSON")

    if flag_plotpt: 
        lonpt_, latpt_ = m(*latlonpt[::-1])
        idxpt = np.unravel_index( ( (lats-latpt_)**2 + (lons-lonpt_)**2 ).argmin(), lats.shape)
        ax.scatter(lons[idxpt],lats[idxpt],c='r',s=5)

    divider = make_axes_locatable(ax)
    cbaxes = divider.append_axes("bottom", size="5%", pad=0.5)
    cbar = fig.colorbar(im,orientation='horizontal',cax = cbaxes)
    cbar.set_label('topography (m)')

    plt.savefig(dir_in+'pgds.png')
    plt.close(fig)
