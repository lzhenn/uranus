import numpy as np
import xarray as xr
from PIL import Image, ImageDraw 

ds_in=xr.open_dataset('/home/lzhenn/Njord_dev/njord_his_d02.nc')
ds_out=xr.open_dataset('/home/lzhenn/array74/workspace/njord_pipeline/njord-implement/domaindb/njord_t1t2/roms_d02_omp.nc')

ds_out['h'].values=ds_in['h'].values

ds_out.to_netcdf('/home/lzhenn/array74/workspace/njord_pipeline/njord-implement/domaindb/njord_t1t2/roms_d02_test.nc')


