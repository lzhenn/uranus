import sys, datetime
import numpy as np
import xarray as xr
from PIL import Image, ImageDraw 

ROMS_FILE='/home/lzhenn/array74/workspace/njord_pipeline/njord-implement/domaindb/njord_t1t2/roms_d02_omp.nc'
BMP_FILE='/home/lzhenn/array74/workspace/njord_pipeline/preprocess/fig/bitmap_roms_d02.new.bmp'
OUT_DIR='/home/lzhenn/array74/workspace/njord_pipeline/njord-implement/domaindb/njord_t1t2/'
ds=xr.open_dataset(ROMS_FILE)

im_w, im_h=ds['mask_rho'].shape

imask=ds['mask_rho']
im = Image.open(BMP_FILE)
values=list(im.getdata())
for x in range(im_w):
    for y in range(im_h):
        imask.values[x,y]=values[y*im_w+x]

imask.values=imask.values/255
ds.to_netcdf(OUT_DIR+'roms_d02_test.nc')
