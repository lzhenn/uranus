import numpy as np
import xarray as xr
from PIL import Image, ImageDraw 

ncfile='/home/lzhenn/array74/workspace/uranus/uranus/domaindb/njord_noluzon/geo_em.d01.nc'
bmpfile='geo_em.d01.bmp'
ds=xr.open_dataset(ncfile)

imask=ds['LANDMASK'].values[0,::-1,:]
print(imask)
im_w, im_h = imask.shape
# Convert the mask to a PIL image
im = Image.fromarray(imask.astype(np.uint8) * 255, mode="L")
# Save the image as a BMP file
im.save(bmpfile)
