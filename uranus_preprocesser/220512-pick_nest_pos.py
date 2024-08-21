import numpy as np
import xarray as xr

def get_closest_idx(lat2d, lon2d, lat0, lon0):
    dis_lat2d=lat2d-lat0
    dis_lon2d=lon2d-lon0
    dis=abs(dis_lat2d)+abs(dis_lon2d)
    idx=np.argwhere(dis==dis.min()).tolist()[0] # x, y position
    return idx 

outerfile_dir='/home/lzhenn/Njord_dev/Projects/Njord_t1t2/roms_swan_grid/roms_d02_lp0d1.3m.nc'
innerfile_dir='/home/lzhenn/cooperate/data/dtm_bathymetry_100m3.nc'
#outer domain
ds_out=xr.open_dataset(outerfile_dir)
# inner domain
ds_in=xr.open_dataset(innerfile_dir)

lat2d_in=ds_in['lat_rho'].values
lon2d_in=ds_in['lon_rho'].values

sw_point=lat2d_in[0,0], lon2d_in[0,0]
ne_point=lat2d_in[-1,-1], lon2d_in[-1,-1]


lat2d_out=ds_out['lat_rho'].values
lon2d_out=ds_out['lon_rho'].values

sw_id=get_closest_idx(lat2d_out, lon2d_out, sw_point[0], sw_point[1])
print(sw_id)
print(get_closest_idx(lat2d_out, lon2d_out, ne_point[0], ne_point[1]))

