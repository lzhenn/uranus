import xarray as xr
import numpy as np
from PIL import Image
work_dir='/home/lzhenn/array74/workspace/uranus/uranus/domaindb/njord_noluzon'
ncfile=f'{work_dir}/geo_em.d01.bck.nc'
bmpfile=f'{work_dir}/geo_em.d01.noluzon.bmp'
outfile=f'{work_dir}/geo_em.d01.nc'


ds=xr.open_dataset(ncfile)
im = Image.open(bmpfile)
bmp_array = np.array(im)
bmp_array = bmp_array[::-1,:] 
print(bmp_array)
ds['LANDMASK'].values[0,:,:]=bmp_array
ds.to_netcdf(outfile)
