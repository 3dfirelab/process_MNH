import numpy as np 
from netCDF4 import Dataset
import argparse
import os, sys
import pandas as pd
from datetime import datetime, timedelta

#homebrewed
#import MNH_tools

def calculate_wind_direction(u, v):
    # Calculate the wind direction in degrees
    wind_direction = (270 - (np.degrees(np.arctan2(v, u)))) % 360
    return wind_direction


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='to extart point info from MNH netcdf postroc file') 
    parser.add_argument('-i','--input', help='path to input nc file',required=True)
    parser.add_argument('-ll','--latlon', help='latlon of point to show on map',required=False)
    parser.add_argument('-o','--dirout', help='path dir out',required=False)
    args = parser.parse_args()

    name_in = args.input

    if args.latlon is not None:
        latlonpt = [ float(xx) for xx in args.latlon.split(',') ]

    name_out = 'weather_MNH_{:04.2f}-{:04.2f}'.format(latlonpt[0],latlonpt[1]).replace('.','_')+'.wxs'

    nc = Dataset(name_in,'r')

    dir_in = os.path.dirname(name_in)
    if args.dirout is None:
        dir_out = dir_in+'/FasrsiteWxsFile/'
    else: 
        dir_out = args.dirout
    #MNH_tools.ensure_dir(dir_out)
    os.makedirs(dir_out, exist_ok=True)
   
    #load data
    topo = nc.variables['topography'][:,:].T
    u = np.transpose(nc.variables['u'][:,:,:,:],[0,1,3,2])
    v = np.transpose(nc.variables['v'][:,:,:,:],[0,1,3,2])
    
    temp = np.transpose(nc.variables['Temp'][:,:,:],[0,1,3,2])
    rh = np.transpose(nc.variables['RH'][:,:,:,:],[0,1,3,2])
    pcp = np.transpose(nc.variables['acprecip'][:,:,:,:],[0,1,3,2])
    cldF = np.transpose(nc.variables['cloudfrac'][:,:,:,:],[0,1,3,2])

    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    timenc = nc.variables['time'][:]
    str_time_ref = nc.variables['time'].units.split('ion:')[1].strip()
    time_ref = datetime.strptime(str_time_ref, "%Y-%m-%d %H:%M:%S.%f")
    time =  []
    for time_sec in timenc:
        time.append(time_ref +  timedelta(seconds=time_sec))
    time= np.array(time)

    yy,xx = np.meshgrid(x,y)
    extent=(xx.min(),xx.max(),yy.min(),yy.max())

    lat = nc.variables['lat'][:,:]
    lon = nc.variables['lon'][:,:]

    ilat = np.abs(lat[:,0]-latlonpt[0]).argmin()
    jlon = np.abs(lon[0,:]-latlonpt[1]).argmin()

    #print(jlon,ilat)

    # Example time series data
    time_series =  pd.to_datetime(time) #pd.date_range(start=time[0].strftime("%Y-%m-%d %H:%M"), end=time[-1].strftime("%Y-%m-%d %H:%M"), freq='h')
    temperature = temp[:,0,ilat,jlon] -273.15   # Example data, replace with your actual time series data
    humidity    = rh[:,0,ilat,jlon] 
    hourly_precipitation = pcp[:,0,ilat,jlon]
    wind_speed = np.sqrt(u[:,0,ilat,jlon]**2 + v[:,0,ilat,jlon]**2)
    wind_direction = calculate_wind_direction(u,v)[:,0,ilat,jlon]
    cloud_coverage = cldF[:,0,ilat,jlon]

# Create a DataFrame
    data = {
        'Year': [dt.year for dt in time_series],
        'Mth': [dt.month for dt in time_series],
        'Day': [dt.day for dt in time_series],
        'Time': [dt.strftime('%H%M') for dt in time_series],
        'Temp': [f'{temp:.0f}' for temp in temperature],
        'RH': [f'{rh:.0f}' for rh in humidity],
        'HrlyPcp': [f'{pcp:.2f}' for pcp in hourly_precipitation],
        'WindSpd': [f'{ws:.0f}' for ws in wind_speed],
        'WindDir': [f'{wd:.0f}' for wd in wind_direction],
        'CloudCov': [f'{cc:.0f}' for cc in cloud_coverage]
    }
    df = pd.DataFrame(data)
    
    '''
    df_wxs = df.copy()
    df_wxs['Temp'] = df_wxs['Temp'].round(0).astype(int)
    df_wxs['RH'] = df_wxs['RH'].round(0).astype(int)
    df_wxs['WindSpd'] = df_wxs['WindSpd'].round(0).astype(int)
    df_wxs['WindDir'] = df_wxs['WindDir'].round(0).astype(int)
    df_wxs['CloudCov'] = df_wxs['CloudCov'].round(0).astype(int)
    '''
    # Write the DataFrame to a CSV file
    with open(dir_out+'/'+name_out, 'w') as f:
        f.write('RAWS_UNITS: Metric\n')
        f.write('RAWS_ELEVATION: {:3.1f}\n'.format(topo[ilat,jlon]))
        f.write('RAWS: {:d}\n'.format(len(df)))
        df.to_csv(f, index=False, sep='\t')

    #write a dataframe with timestamp as well
    # Combine 'Year', 'Mth', 'Day', and 'Time' into a single 'datetime' column
    df['Timestamp'] = pd.to_datetime(df['Year'].astype(str) + '-' +
                                     df['Mth'].astype(str) + '-' +
                                     df['Day'].astype(str) + ' ' +
                                     df['Time'].str.zfill(4).str[:2] + ':' +
                                     df['Time'].str.zfill(4).str[2:])

    # Set 'Timestamp' as the index
    df.set_index('Timestamp', inplace=True)

    # Drop the old Year, Mth, Day, and Time columns (optional)
    df.drop(['Year', 'Mth', 'Day', 'Time'], axis=1, inplace=True)
    
    #convert to numeric
    df = df.apply(pd.to_numeric, errors='coerce')

    df.to_json(dir_out+'/'+name_out.replace('wxs','json'), orient='records', date_format='iso')
