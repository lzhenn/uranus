import numpy as np
import xarray as xr
from PIL import Image, ImageDraw 

ds=xr.open_dataset('/home/lzhenn/array74/workspace/njord_pipeline/njord-implement/domaindb/njord_t1t2/roms_d01_omp.nc')

imask=ds['mask_rho']
im_w,im_h=imask.shape

imask=imask.values.reshape((im_w,im_h))

image = Image.new('1', (im_w, im_h), 0)
draw = ImageDraw.Draw(image)
for x in range(im_w):
    for y in range(im_h):
        draw.point((x, y), fill=int(imask[x,y]))
image.save('../fig/bitmap_roms_d01.bmp', 'bmp')
