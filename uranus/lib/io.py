#!/usr/bin/env python3
"""specific module for IO"""
# ---imports---
import os, glob, shutil
import numpy as np
import pandas as pd
import struct
from . import utils, const

# ---Module regime consts and variables---
print_prefix='lib.io>>'

def copy_files(src_path, dest_path):
    utils.write_log(
        f'{print_prefix}Copy files from {src_path} to {dest_path}...')
    for file in glob.glob(rf'{src_path}'):
        shutil.copy(file, dest_path)
def move_files(src_path, dest_path):
    utils.write_log(
        f'{print_prefix}Move files from {src_path} to {dest_path}...')
    for file in glob.glob(rf'{src_path}'):
        shutil.move(file, dest_path)
        
def symlink_files(src_path, dest_path):
    utils.write_log(
        f'{print_prefix}Symlink files from {src_path} to {dest_path}...')
    os.system(f'ln -sf {src_path} {dest_path}')        
# ---Classes and Functions---
def del_files(tgt_path, fnpatterns):
    file_list = os.listdir(tgt_path)
    utils.write_log(f'{print_prefix}Clean workspace for {tgt_path}...')
    for filename in file_list:
        if filename.startswith(tuple(fnpatterns)):
            os.remove(os.path.join(tgt_path, filename))

def get_wrf_fn(tgt_time, wrf_domain):
    '''
    return aimed wrf file name given tgt_time and wrf_domain
    '''
    return 'wrfout_'+wrf_domain+'_'+tgt_time.strftime('%Y-%m-%d_%H:00:00')
        
def sel_frm(da_src, tf, itf):
    try:
        da=da_src.sel(time=tf, method='nearest')
    except TypeError:
        da=da_src.isel(time=itf)
    return da

def gen_patternfn_lst(drv_root, drv_dic, inittime, 
        endtime, kw='atm',special=''):
    
    fn_lst=[]
    # do not include the last time frame
    tfs=pd.date_range(start=inittime, end=endtime, 
                      inclusive='left', freq=drv_dic[f'{kw}_file_nfrq'])
    pattern=drv_dic[f'{kw}_naming_pattern']
    if '$I' in pattern:
        init_fmt=pattern.split('$I')[1]
        init_str = inittime.strftime(init_fmt)
        pattern=pattern.replace(f'$I{init_fmt}$I', init_str)
    if '$S' in pattern:
        pattern=pattern.replace(f'$S', special)
        
    frm_fmt=pattern.split('$F')[1]
    for tf in tfs: 
        frm_str = tf.strftime(frm_fmt)
        # Replace the placeholders with the formatted strings
        fn=pattern.replace(f'$F{frm_fmt}$F', frm_str)
        fn_full=os.path.join(drv_root, fn)
        if not(os.path.exists(fn_full)):
            utils.throw_error(f'file not exist: {fn_full}')
        fn_lst.append(fn_full)
    return fn_lst, tfs
def gen_roms_forc_template(ts, lat2d, lon2d):
    '''generate the forcing template file for ROMS
        ts: time series, in days relative to basetime
        lat2d, lon2d: 2D arrays
    '''
    import xarray as xr
    
    nt, (nx,ny)=len(ts), lat2d.shape
    # Create the dimensions
    xr_dim = xr.IndexVariable('xr', np.arange(nx))
    er_dim = xr.IndexVariable('er', np.arange(ny))
    time_dim = xr.IndexVariable('time', ts)
    coords={'xr': xr_dim, 'er': er_dim, 'time': time_dim}


    # Create the variables
    ds_lon = xr.Variable(('er', 'xr'), np.zeros((ny, nx)), attrs={
        'long_name': 'longitude', 'units': 'degrees_east', 'field': 'xp, scalar, series'})
    ds_lat = xr.Variable(('er', 'xr'), np.zeros((ny, nx)), attrs={
        'long_name': 'latitude', 'units': 'degrees_north', 'field': 'yp, scalar, series'})
    ocean_time = xr.Variable('time', ts, attrs={
        'long_name': 'atmospheric forcing time', 'units': 'days', 'field': 'time, scalar, series'})

    var_xarr={}
    for var in const.FORC_TIMEVAR_LIST:
        var_xarr[var+'_time']=xr.Variable(var+'_time', ts,
            attrs={'long_name':var+'_time', 
                'units':'days since 1858-11-17 00:00:00 UTC',
                'field':'wind_time, scalar, series'})
        
    for key, attr in const.FORC_VAR_DIC.items():
        timedim=attr['time']
        var_xarr[key]=xr.Variable((timedim, 'er', 'xr'), np.zeros((nt, ny, nx)), attrs=attr)
    # Assign Values
    ds_lon.values=lon2d.values
    ds_lat.values=lat2d.values
    var_xarr['ocean_time']=ocean_time
    var_xarr['lat'],var_xarr['lon']=ds_lat,ds_lon

    # Create the dataset
    ds = xr.Dataset(var_xarr, coords=coords)

    # Add global attributes
    ds.attrs['type'] = 'bulk fluxes forcing file'
    ds.attrs['gridid'] = 'combined grid'

    return ds

def gen_wrf_mid_template(drv_key):
    # init the dictionary
    drv_dic=const.DRV_DIC[drv_key]
    lats,lons=drv_dic['lats'],drv_dic['lons']
    NLAT,NLON=len(lats),len(lons)
    slab_dict={
        'IFV':5, 'HDATE':'0000-00-00_00:00:00:0000',
        'XFCST':0.0, 'MAP_SOURCE':drv_key.ljust(32),
        'FIELD':'', 'UNIT':'', 'DESC':'', 
        'XLVL':0.0, 'NX':NLON, 'NY':NLAT,
        'IPROJ':0,'STARTLOC':'SWCORNER',
        'STARTLAT':lats[0], 'STARTLON':lons[0],
        'DELTLAT':lats[1]-lats[0], 'DELTLON':lons[1]-lons[0], 'EARTH_RAD':6371.229,
        'IS_WIND_EARTH_REL': 0, 
        'SLAB':np.array(np.zeros((NLAT,NLON)), dtype=np.float32),
        'key_lst':['IFV', 'HDATE', 'XFCST', 'MAP_SOURCE', 'FIELD', 'UNIT', 
        'DESC', 'XLVL', 'NX', 'NY', 'IPROJ', 'STARTLOC', 
        'STARTLAT', 'STARTLON', 'DELTLAT', 'DELTLON', 
        'EARTH_RAD', 'IS_WIND_EARTH_REL', 'SLAB']
    }
    return slab_dict

def write_record(out_file, slab_dic):
    '''
    Write a record to a WRF intermediate file
    '''
    slab_dic['FIELD']=slab_dic['FIELD'].ljust(9)
    slab_dic['UNIT']=slab_dic['UNIT'].ljust(25)
    slab_dic['DESC']=slab_dic['DESC'].ljust(46)
    
    # IFV header
    out_file.write_record(struct.pack('>I',slab_dic['IFV']))
    
    # HDATE header
    pack=struct.pack('>24sf32s9s25s46sfIII', 
        slab_dic['HDATE'].encode(), slab_dic['XFCST'],
        slab_dic['MAP_SOURCE'].encode(), slab_dic['FIELD'].encode(),
        slab_dic['UNIT'].encode(), slab_dic['DESC'].encode(),
        slab_dic['XLVL'], slab_dic['NX'], slab_dic['NY'],
        slab_dic['IPROJ'])
    out_file.write_record(pack)

    # STARTLOC header
    pack=struct.pack('>8sfffff',
        slab_dic['STARTLOC'].encode(), slab_dic['STARTLAT'],
        slab_dic['STARTLON'], slab_dic['DELTLAT'], slab_dic['DELTLON'],
        slab_dic['EARTH_RAD'])
    out_file.write_record(pack)

    # IS_WIND_EARTH_REL header
    pack=struct.pack('>I', slab_dic['IS_WIND_EARTH_REL'])
    out_file.write_record(pack)

    # Let's play with the SLAB
    out_file.write_record(
        slab_dic['SLAB'].astype('>f'))

# ---Unit test---
if __name__ == '__main__':
    pass
