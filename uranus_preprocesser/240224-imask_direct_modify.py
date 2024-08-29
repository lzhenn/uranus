import numpy as np
import xarray as xr
from PIL import Image, ImageDraw 

ncfile='/home/lzhenn/array74/Njord_Calypso/domaindb/poseidon_1500m_L12/roms_d03_omp.nc.bck'
bmpfile='../fig/roms_d03_omp.bmp'
ds=xr.open_dataset(ncfile)

