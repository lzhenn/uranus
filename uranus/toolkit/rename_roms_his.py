#/usr/bin/env python3
"""
Rename ROMS history files 
and move to archive
"""
import subprocess
import xarray as xr

domid='d03'
src_dir='/home/lzhenn/array130/poseidon/2018091200_noluzon_fc1.2/roms_his'
arch_dir='/home/lzhenn/array130/poseidon/2018091200_noluzon_fc1.2'

fn_stream=subprocess.check_output(
    f'ls {arch_dir}/roms_his_{domid}*', shell=True).decode('utf-8')
archfn_list=fn_stream.split()

fn_stream=subprocess.check_output(
    f'ls {src_dir}/roms_his_{domid}*', shell=True).decode('utf-8')
srcfn_list=fn_stream.split()

ds_roms=xr.load_dataset(srcfn_list[0])
src_stf=ds_roms['ocean_time'][0].dt
ymdhms0=int(src_stf.strftime('%Y%m%d%H%M%S').values)

print(f'src start: {ymdhms0}')
for fn in archfn_list[::-1]:
    ds_roms=xr.open_dataset(fn)
    arch_tf=ds_roms['ocean_time'][0].dt
    ymdhms=int(arch_tf.strftime('%Y%m%d%H%M%S').values)
    print(f'archive end: {ymdhms}')
    if ymdhms>ymdhms0:
        print(f'archive {fn}')
    else:
        print(f'not {fn}')
    