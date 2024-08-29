import xarray as xr
import numpy as np
from PIL import Image
work_dir='/home/lzhenn/array74/workspace/uranus/uranus/domaindb/njord_noluzon/'
infile=f'{work_dir}/geo_em.d01.nc'
reffile=f'{work_dir}/geo_em.d01.bck.nc'
outfile=f'{work_dir}/geo_em.d01.new.nc'

var2d={'CLAYFRAC':0,'HGT_M':0,'LU_INDEX':17,'OA1':0,'OA2':0,'OA3':0,'OL1':0,'OL2':0,'OL3':0,'OL4':0,
       'SANDFRAC':0,'SCB_DOM':14,'SCT_DOM':14,'SNOALB':0,'SOILTEMP':0,'VAR':0,'VAR_SSO':0}
var3d={'ALBEDO12M':8,'GREENFRAC':0,'LAI12M':0} # all 12-dim
var3dsp={'LANDUSEF':16,'SOILCBOT':13,'SOILCTOP':13} #[num,:,:]=1.0, others=0.0


ds=xr.open_dataset(infile)
imask=ds['LANDMASK'].values

ds_ref=xr.open_dataset(reffile)
maskchg=np.logical_xor(ds_ref['LANDMASK'].values[0,:,:],imask[0,:,:])

for var in var2d:
    print('var2d=',var)
    ds[var].values=xr.where(maskchg==1,var2d[var],ds[var].values)
for var in var3d:
    print('var3d=',var)
    for i in range(0,12):
        ds[var].values[0,i,:,:]=xr.where(maskchg==1,var3d[var],ds[var].values[0,i,:,:]) 
for var in var3dsp:
    print('var3dsp=',var)
    val=var3dsp[var]
    _,len_var,_,_=ds[var].shape
    for i in range(0,len_var):
        ds[var].values[0,i,:,:]=xr.where(maskchg==1,0,ds[var].values[0,i,:,:])
    ds[var].values[0,val,:,:]=xr.where(maskchg==1,1.0,ds[var].values[0,val,:,:])
ds.to_netcdf(outfile)