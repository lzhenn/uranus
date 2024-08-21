import numpy as np
import xarray as xr
from PIL import Image, ImageDraw 

ds_in=xr.open_dataset('/home/lzhenn/array74/workspace/njord_pipeline/njord-implement/domaindb/njord_t1t2/roms_d02_test.nc')

h=ds_in['h']
h[1,:]=(h[0,:]+h[2,:])/2
h[-2,:]=(h[-1,:]+h[-3,:])/2
h[:,1]=(h[:,0]+h[:,2])/2
h[:,-2]=(h[:,-1]+h[:,-3])/2

ds_in.to_netcdf('/home/lzhenn/array74/workspace/njord_pipeline/njord-implement/domaindb/njord_t1t2/roms_d02_omp.nc')


