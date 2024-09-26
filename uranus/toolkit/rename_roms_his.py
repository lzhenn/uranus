#/usr/bin/env python3
"""
Rename ROMS history files 
and move to archive
"""
import os, subprocess
import xarray as xr

domid='d03'
src_dir='/home/lzhenn/array130/poseidon/2018091200_noluzon/roms_his'
arch_dir='/home/lzhenn/array130/poseidon/2018091200_noluzon'

fn_stream=subprocess.check_output(
    f'ls {arch_dir}/roms_his_{domid}*', shell=True).decode('utf-8')
archfn_list=fn_stream.split()
last_fn=archfn_list[-1]
ds_roms=xr.open_dataset(last_fn)
arch_tf=ds_roms['ocean_time'][0].dt
ymdhms_last=arch_tf.strftime('%Y%m%d%H%M%S').values
last_fn=os.path.basename(last_fn)
last_num=int(last_fn.split('_')[-1].split('.')[0])
print(f'archive end: {ymdhms_last} in {last_fn}')

fn_stream=subprocess.check_output(
    f'ls {src_dir}/roms_his_{domid}*', shell=True).decode('utf-8')
srcfn_list=fn_stream.split()

for fn in srcfn_list:
    ds_roms=xr.open_dataset(fn)
    arch_tf=ds_roms['ocean_time'][0].dt
    ymdhms=arch_tf.strftime('%Y%m%d%H%M%S').values
    base=os.path.basename(fn)
    if ymdhms>ymdhms_last:
        last_num+=1
        print(f'mv {fn} to {arch_dir}/roms_his_{domid}_{last_num:05}.nc')
        subprocess.run(f'cp {fn} {arch_dir}/roms_his_{domid}_{last_num:05}.nc', shell=True)
    
