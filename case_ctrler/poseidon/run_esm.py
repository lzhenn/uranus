import numpy as np
import xarray as xr
import os, subprocess 
import uranus

nesm=4
istrt=3
inipath='/home/lzhenn/Poseidon/Projects/poseidon_LTtmr_L12/'
pert_flag=True
pert_var='temp'
domid='d03'

int_file=os.path.join(inipath,f'roms_ini_{domid}.nc')
for i in range(istrt, nesm):
    if pert_flag:
        ds=xr.load_dataset(int_file)
        pert = np.random.normal(loc=0.0, scale=0.01, size=ds[pert_var].shape)
        ds[pert_var]=ds[pert_var]+pert
        outfn=os.path.join(inipath,f'roms_ini_{domid}_{i:02d}.nc')
        print('perturbed init: '+outfn)
        ds.to_netcdf(outfn)
for i in range(istrt, nesm):
    srcfn=os.path.join(inipath,f'roms_ini_{domid}_{i:02d}.nc')
    aimfn=os.path.join(inipath,f'roms_ini_{domid}.nc')
    cmd=f'mv {aimfn} roms_ini_{domid}.nc.bck; cp {srcfn} {aimfn}'
    print('mv perturbed init: '+cmd)
    rcode=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    agent=uranus.Uranus(cfgfn=f'config.poseidon.{domid}.ini')
    agent.arch_root=os.path.join(agent.arch_root, f'esm_{i:02d}')
    agent.cfg['ROMS']['gen_roms_icbc']='False'
    agent.waterfall()
