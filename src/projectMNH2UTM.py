import xarray as xr
import numpy as np
import pyproj
from scipy.interpolate import griddata
import argparse
import sys
import multiprocessing
import os 
################################################
def cpu_count():
    try:
        return int(os.environ['ntask'])
    except:
        print('env variable ntask is not defined')
        sys.exit()
        #return multiprocessing.cpu_count()


#####################################################
def griddata_(data2d):
    return griddata((x.flatten(), y.flatten()), data2d.flatten(), (x_mesh, y_mesh), method="linear")

def griddata_star(arg):
    return griddata_(*arg)  

#####################################################
if __name__ == "__main__":
#####################################################

    parser = argparse.ArgumentParser(description='convert 2d latlon netcdf to 1d xy')
    parser.add_argument('-i','--input', help='input nc file',required=True)
    parser.add_argument('-utmproj','--utmprojection', help='projection of the output file',required=False)
    parser.add_argument('-p','--parallel', help='run with mulriprocessing',required=False, default=False)
    args = parser.parse_args()
    
    flag_parallel = args.parallel

# Load the NetCDF file
    ds = xr.open_dataset(args.input, decode_times=False)

# Extract 2D latitude, longitude, and a data variable (e.g., temperature)
    lat = ds["lat"].values  # Replace with actual variable name
    lon = ds["lon"].values  # Replace with actual variable name
    
    variables_4d = [var for var in ds.data_vars if ds[var].ndim == 4]
    variables_3d = [var for var in ds.data_vars if ds[var].ndim == 3]
    variables_2d = [var for var in ds.data_vars if ds[var].ndim == 2]
    
    nx = ds[variables_2d[0]].shape[1]
    ny = ds[variables_2d[0]].shape[0]
    nz = ds[variables_3d[0]].shape[0]
    nt = ds[variables_4d[0]].shape[0]

#create list of variables

# Define the projections
    proj_wgs84 = pyproj.CRS("EPSG:4326")  # WGS84
    if args.utmprojection is None: 
        utm_zone = int((lon.mean() + 180) / 6) + 1  # Calculate UTM zone
        proj_utm = pyproj.CRS(f"EPSG:326{utm_zone}" if lat.mean() >= 0 else f"EPSG:327{utm_zone}")  # Northern or Southern Hemisphere
    else:  
        proj_utm = pyproj.CRS("EPSG:{:s}".format(args.utmprojection)) 

# Create a transformer
    transformer = pyproj.Transformer.from_crs(proj_wgs84, proj_utm, always_xy=True)


# Convert lat/lon to UTM
    x, y = transformer.transform(lon, lat)

# Define a regular UTM grid (change resolution as needed)
    x_grid = np.linspace(x.min(), x.max(), nx)
    y_grid = np.linspace(y.min(), y.max(), ny)
    x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)

# Initialize the dataset with coordinates (empty for now)
    ds_utm = xr.Dataset(
        coords={
            "time": ds.time,
            "z": ds.z,
            "x": x_grid,
            "y": y_grid,
        }
    )
    ds_utm = ds_utm.rio.write_crs(args.utmprojection)

# Interpolate data to the regular UTM grid
    for varlist,nbdim in zip([variables_4d,variables_3d,variables_2d],[4,3,2]):
        print(varlist)
        for ivar, var in enumerate(varlist):
            print('{:d} {:d}/{:d} {:s}'.format(nbdim,ivar,len(varlist),var) )
            data = ds[var]
            
            if nbdim == 4: 
                args_here = []
                for it in range(nt):
                  for k in range(nz):
                      #print(it,k)
                      data2d = data.values[it,k] 
                      args_here.append([data2d])
                
                if flag_parallel: 
             
                    # set up a pool to run the parallel processing
                    cpus = cpu_count()
                    pool = multiprocessing.Pool(processes=cpus)

                    # then the map method of pool actually does the parallelisation  
                    data_utm = pool.map(griddata_star, args_here)
                    pool.close()
                    pool.join()
                    
                else:
                    data_utm = [] 
                    for arg in args_here:
                         data_utm.append(griddata_(arg[0]))

                ds_utm[var] = (["time", "z", "y", "x"], np.array(data_utm).reshape(nt,nz,ny,nx) )
            
            if nbdim == 3: 
                
                args_here = []
                if 'time' in list(ds[var].sizes): 
                    nn = nt
                    dimname = 'time'
                if 'z' in list(ds[var].sizes): 
                    nn = nz
                    dimname = 'z'

                for k in range(nn):
                    data2d = data.values[k] 
                    args_here.append([data2d])

                if flag_parallel: 
             
                    # set up a pool to run the parallel processing
                    cpus = cpu_count()
                    pool = multiprocessing.Pool(processes=cpus)

                    # then the map method of pool actually does the parallelisation  
                    data_utm = pool.map(griddata_star, args_here)
                    pool.close()
                    pool.join()
                    
                else:
                    data_utm = [] 
                    for arg in args_here:
                         data_utm.append(griddata_(arg[0]))

                ds_utm[var] = ([dimname, "y", "x"], np.array(data_utm).reshape(nn,ny,nx) )
            
            if nbdim == 2: 
                data2d = data.values 
                data_utm =  griddata((x.flatten(), y.flatten()), data2d.flatten(), (x_mesh, y_mesh), method="linear")
            
                ds_utm[var] = ([ "y", "x"], data_utm )

# Save to a new NetCDF file
    ds_utm.to_netcdf(args.input.replace('.nc','_utm.nc'))
    
    crs_wkt = ds_utm.rio.crs.to_wkt()
    # Save the CRS to a .prj file
    prj_file_path = args.input.replace('.nc','_utm.prj')
    with open(prj_file_path, 'w') as prj_file:
        prj_file.write(crs_wkt)
    
    print("Projection complete! Saved as", args.input.replace('.nc','_utm.nc'))
