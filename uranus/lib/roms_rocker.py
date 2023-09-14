#/usr/bin/env python3
"""Build ROMS workflow"""
import os 
import pandas as pd
import xarray as xr
import numpy as np
from . import utils, io, const, mathlib

# ---Module regime consts and variables---
print_prefix='lib.ROMSRocker>>'


# ---Classes and Functions---
class ROMSRocker:

    '''
    Construct icbcmaker to generate initial and boundary file for ROMS
    '''
    
    def __init__(self, uranus):
        """ construct maker obj """
        self.uranus=uranus
        self.cfg=uranus.cfg
        romscfg=self.cfg['ROMS']
        self.run_maker=romscfg.getboolean('gen_roms_icbc')
        self.strt_time,self.end_time=uranus.sim_strt_time,uranus.sim_end_time
        self.drv_type=romscfg['drv_type']
        self.drv_dic=const.DRV_DIC[self.drv_type+'_roms']
        self.drv_root=utils.parse_fmt_timepath(self.strt_time, romscfg['drv_root'])
        self.proj_root=uranus.proj_root
        self.domdb_root=uranus.domdb_root
        self.dt=int(romscfg['dt'])
        self.nrst=int(romscfg['nrst'])
        self.nhis=int(romscfg['nhis'])
        if not(os.path.exists(self.proj_root)):
            utils.write_log(print_prefix+'mkdir '+self.proj_root)
            os.mkdir(self.proj_root)
        
        # build meta
        self.frm_time_series=pd.date_range(
            start=self.strt_time, end=self.end_time, freq=self.drv_dic['ocn_nfrq'])
        self.ocnfn_lst, self.file_time_series=io.gen_patternfn_lst(
            self.drv_root, self.drv_dic, self.strt_time, self.end_time, kw='ocn')

    def prepare_rock(self):
        nml_temp=self.uranus.nml_temp
        domfn=os.path.join(
                self.domdb_root, nml_temp, '*omp.nc')
        io.symlink_files(domfn, self.proj_root)
            
        roms_in=os.path.join(self.proj_root, 'roms_d01.in')
        run_seconds=int((self.end_time-self.strt_time).total_seconds())
        ntimes=run_seconds//self.dt
        nrst=self.nrst*const.DAY_IN_SEC//self.dt
        nhis=self.nhis*const.HR_IN_SEC//self.dt 
        ndefhis=const.DAY_IN_SEC//self.dt
        
        grdname=os.path.join('Projects', nml_temp, 'roms_d01_omp.nc')
        ininame=os.path.join('Projects', nml_temp, 'coawst_ini.nc')
        
        tf_str=self.file_time_series[0].strftime('%Y%m%d')
        clmname=os.path.join('Projects', nml_temp, f'coawst_clm_{tf_str}.nc')
        bryname=os.path.join('Projects', nml_temp, f'coawst_bdy_{tf_str}.nc')
        
        sed_dic={
            'NtileI':self.uranus.ntasks_iocn, 'NtileJ':self.uranus.ntasks_jocn,
            'NTIMES':ntimes, 'DT':f"{self.dt:.1f}d0", 'NRST':nrst, 'NHIS':nhis, 
            'NDEFHIS':ndefhis, 'NDIA':nhis, 'NAVG':nhis, 
            'GRDNAME':grdname, 'ININAME':ininame, 'CLMNAME':clmname, 'BRYNAME':bryname
        }
        for key, itm in sed_dic.items():
            utils.sedline(key, f'{key} == {itm}', roms_in)
  
    def load_domain(self):
        """ load domain file """
        utils.write_log(print_prefix+'Load domain file...')
        domfn=os.path.join(
            self.domdb_root, self.uranus.nml_temp, 'roms_d01_omp.nc')
        self.ds_static=xr.load_dataset(domfn)
        ds_static=self.ds_static
        self.mask=ds_static['mask_rho'].values
        self.h=ds_static['h'].values
        self.lat1d,self.lon1d=ds_static['lat_rho'][:,0].values,ds_static['lon_rho'][0,:].values
        self.lat_u, self.lon_u=ds_static['lat_u'][:,0].values, ds_static['lon_u'][0,:].values
        self.lat_v, self.lon_v=ds_static['lat_v'][:,0].values, ds_static['lon_v'][0,:].values
       
        # load sample file
        sampfn=os.path.join(
            self.uranus.domdb_root, self.uranus.nml_temp, 'roms_d01_inismp.nc')
        self.ds_smp=xr.load_dataset(sampfn)
        self.hc=self.ds_smp['hc']
        
        # ONLY FOR REGULAR LAT LON!!! and Mercator projection
        self.ds_smp=self.ds_smp.rename_dims({
            'erho':'lat','xrho':'lon'})
        self.ds_smp=self.ds_smp.assign_coords({
            'lon': self.lon1d, 
            'lat': self.lat1d}) 
        # generate data template
        self.roms_3dtemplate=self.ds_smp['temp']
  
    def build_icbc(self):
        """ build icbcs for ROMS """
        if not(self.run_maker):
            return
        utils.write_log(
            print_prefix+'build icbcs from %s to %s...'%(
                self.strt_time.strftime('%Y%m%d'),self.end_time.strftime('%Y%m%d')))
        
        # load domain file
        self.load_domain()
        
                
        for it, tf in enumerate(self.frm_time_series):
            utils.write_log(
                print_prefix+'build icbcs@%s...'% tf.strftime('%Y%m%d%H'))

            if tf in self.file_time_series:
                # load raw file  
                self.load_raw(it) 
            
            # first deal with zeta
            roms_var='zeta'
            if self.drv_type=='cfs':
                utils.write_log(
                    print_prefix+'cfs data set zeta=0...')
                self.ds_smp['zeta'].values[:]=0.0
            else:
                self.inter2d(roms_var)
            # calculate vertical coordinate before 3d interpolation
            zeta = self.ds_smp['zeta'] 
            h=zeta # obtain dimensional information of zeta  
            h= self.ds_static.h.values
            z_rho = mathlib.sigma2depth(zeta, h, self.ds_smp)
            self.dz=z_rho[1:,:,:]-z_rho[:-1,:,:]
            self.dz=np.concatenate((self.dz,self.dz[-1:,:,:]),axis=0)
            self.dp_idx=mathlib.assign_depth_idx(
                z_rho, self.mask)

            for roms_var in ['temp','salt','u','v']:
            #for roms_var in ['temp']:
                self.inter3d(roms_var)
            if tf==self.strt_time:
                # pkg time
                time_offset=self.strt_time - const.BASE_TIME
                self.ds_smp['ocean_time'].values[:]=int(time_offset.total_seconds())*const.S2NS
                self.ds_smp=self.ds_smp.assign_coords({'ocean_time':self.ds_smp['ocean_time']})

                # output
                inifn=os.path.join(self.proj_root, 'coawst_ini.nc')
                self.ds_smp.to_netcdf(inifn)
                utils.write_log(print_prefix+'build initial conditions done!')    
            self.build_clm(tf)
            self.build_bdy(tf)
   
    def load_raw(self, it):
        """ load raw GCM files """
        fn = self.ocnfn_lst[it]
        self.ds_raw=xr.load_dataset(fn)
        ds_raw=self.ds_raw
        self.varname_remap()
        # deal with previous hycom
        if self.drv_type=='hycom':
            if self.varmap['lat_rho']=='Y':
                ds_raw=ds_raw.rename_dims({
                    'Depth':'depth','Y':'lat','X':'lon'})
                ds_raw=ds_raw.assign_coords({
                    'depth':ds_raw['Depth'],
                    'lon': ds_raw['Longitude'][0,:], 
                    'lat': ds_raw['Latitude'][:,0]}) 
        elif self.drv_type=='cfs':
            ds_raw=ds_raw.rename_dims({
                'depthBelowSea':'depth',
                'latitude':'lat','longitude':'lon'})
            ds_raw=ds_raw.assign_coords({
                'depth':ds_raw['depthBelowSea'],
                'lon': ds_raw['longitude'], 
                'lat': ds_raw['latitude']})
            ds_raw=ds_raw.drop_indexes(
                ['depthBelowSea','latitude','longitude'])
            ds_raw['pt'].values=ds_raw['pt'].values-const.K2C
            ds_raw['s'].values=ds_raw['s'].values*1000.0
        elif self.drv_type=='cpsv3':
            ds_raw=ds_raw.rename_dims({
                'st_ocean':'depth',
                'yt_ocean':'lat','xt_ocean':'lon',
                'yu_ocean':'lat','xu_ocean':'lon'})
            ds_raw=ds_raw.drop_indexes(
                ['st_ocean','xt_ocean','yt_ocean','xu_ocean','yu_ocean'])
            ds_raw=ds_raw.assign_coords({
                'depth':ds_raw['st_ocean'],
                'lon': ds_raw['xt_ocean'], 
                'lat': ds_raw['yt_ocean']})
            ds_raw['lon']=ds_raw['lon']+360.
        self.ds_raw=ds_raw
    def varname_remap(self):
        """ remap variable names to roms standard """
        self.varmap={}
        # loop the variables        
        for raw_var in self.ds_raw.data_vars:
            for roms_var in const.ROMS_VAR_NAME_MAPS:
                if raw_var in const.ROMS_VAR_NAME_MAPS[roms_var]:
                    self.varmap[roms_var]=raw_var
        # try dims then
        for raw_dim in self.ds_raw.dims:
            for roms_dim in const.ROMS_VAR_NAME_MAPS:
                if raw_dim in const.ROMS_VAR_NAME_MAPS[roms_dim]:
                    self.varmap[roms_dim]=raw_dim
        # check all included       
        for roms_var in const.ROMS_VAR_NAME_MAPS:
            if roms_var not in self.varmap:
                if self.drv_type=='cfs' and roms_var=='zeta':
                    self.varmap[roms_var]='zeta'
                    return
                utils.throw_error(
                    print_prefix+'%s not found in raw map!'% roms_var)
    
    def inter2d(self, roms_varname):
        '''
        first fill the missing value
        then interpolate to new 2d grid 
        '''

        utils.write_log(
            print_prefix+roms_varname+' 2d-interp...')
        roms_lat, roms_lon=self.lat1d, self.lon1d        
        raw_varname=self.varmap[roms_varname]
        
        var=self.ds_raw[raw_varname]
        var = var.interpolate_na(
            dim='lon', method="nearest", 
            fill_value="extrapolate") 
        
        self.ds_smp[roms_varname].values = var.interp(
            lon=roms_lon, lat=roms_lat,
            method='linear').values
    def inter3d(self, roms_varname):
        
        utils.write_log(
            print_prefix+roms_varname+' 3d-interp step1: raw data fill...')
        
        roms_lat, roms_lon=self.lat1d, self.lon1d        
        raw_varname=self.varmap[roms_varname]
        raw_var=self.ds_raw[raw_varname]

        data_template=self.roms_3dtemplate.copy(deep=True) 
        nt, nz, ny, nx=data_template.shape 

        NZ_DP=len(const.DENSE_DP)

        # fill missing value
        raw_var = raw_var.interpolate_na(
            dim="depth", method="nearest", fill_value="extrapolate")
        raw_var = raw_var.interpolate_na(
            dim="lon", method="nearest", fill_value="extrapolate")
        # horizontal interpolate before vertical interpolate
        raw_var = raw_var.interp(
            lat=roms_lat,lon=roms_lon, method='linear')

        # vertical interpolate to DENSE_DP
        raw_var=raw_var.interp(
            depth=const.DENSE_DP, method='linear',
            kwargs={"fill_value": "extrapolate"})
        #print(raw_var.isnull().sum())
         
        utils.write_log(
            print_prefix+roms_varname+' 3d-interp step2: idx-assign interp...')
        if self.drv_type in ['cpsv3', 'hycom']:
            data_template[:,-1,:,:]=raw_var.values[:,0,:,:]
            for iz in range(0,nz-1):
                idx2d=self.dp_idx[iz,:,:]
                idx3d=np.broadcast_to(idx2d,(nt,NZ_DP,ny,nx))
                data_template[:,iz,:,:]=np.take_along_axis(
                    raw_var.values,idx3d,axis=1)[:,0,:,:]
        else:
            data_template[:,-1,:,:]=raw_var.values[0,:,:]
            for iz in range(0,nz-1):
                idx2d=self.dp_idx[iz,:,:]
                idx3d=np.broadcast_to(idx2d,(NZ_DP,ny,nx))
                data_template[:,iz,:,:]=np.take_along_axis(
                    raw_var.values,idx3d,axis=0)[0,:,:]
       # deal with uv
        if roms_varname =='u':
            dz_4d=np.broadcast_to(self.dz,(nt,nz,ny,nx))
            ubar_rho=data_template[:,0,:,:]
            ubar_rho = (data_template*dz_4d).sum(dim='sc_r')/self.h
            self.ds_smp['u'].values= data_template.interp(
                lat=self.lat_u,lon=self.lon_u,method='linear').values
            
            self.ds_smp['ubar'].values = ubar_rho.interp(
                lat=self.lat_u,lon=self.lon_u,method='linear')
            

        elif roms_varname=='v':
            dz_4d=np.broadcast_to(self.dz,(nt,nz,ny,nx))
            vbar_rho=data_template[:,0,:,:]
            vbar_rho = (data_template*dz_4d).sum(dim='sc_r')/self.h
            self.ds_smp['v'].values = data_template.interp(
                lat=self.lat_v,lon=self.lon_v,method='linear')
        
            self.ds_smp['vbar'].values = vbar_rho.interp(
                lat=self.lat_v,lon=self.lon_v,method='linear')
        else: 
            self.ds_smp[roms_varname].values=data_template.values
        
    def build_clm(self, time_frm):
        '''Build climatology file'''
        utils.write_log(
            print_prefix+'build clm file@%s...' % time_frm.strftime('%Y%m%d'))
        clmsmpfn=os.path.join(self.domdb_root, self.uranus.nml_temp,'roms_d01_clmsmp.nc')
        ds_clm=xr.load_dataset(clmsmpfn)
        # loop the variables to assign values
        for varname in const.ROMS_VAR:
            ds_clm[varname].values=self.ds_smp[varname].values
 
        # deal with time vars 
        time_offset=time_frm - const.BASE_TIME
        for var in const.CLM_TIME_VAR:
            var_time=ds_clm[var+'_time']
            # pkg time
            var_time.values[:]=int(time_offset.total_seconds())*const.S2NS+const.HALF_DAY
            #ds_clm=ds_clm.assign_coords({var:var_time})
        
        clmfn=os.path.join(
            self.proj_root,'coawst_clm_%s.nc'% time_frm.strftime('%Y%m%d'))
        ds_clm.to_netcdf(clmfn)

    def build_bdy(self, time_frm):
        '''Build bdy file'''
        utils.write_log(
            print_prefix+'build bdy file@%s...' % time_frm.strftime('%Y%m%d'))
        bdysmpfn=os.path.join(
            self.domdb_root, self.uranus.nml_temp,'roms_d01_bdysmp.nc')
        ds_bdy=xr.load_dataset(bdysmpfn)
        
        for varname in ['zeta','ubar','vbar']:
            var_bdy=varname+'_south'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,0,:]
            var_bdy=varname+'_north'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,-1,:]
            var_bdy=varname+'_west'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,:,0]
            var_bdy=varname+'_east'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,:,-1]
        
        for varname in ['u','v','temp','salt']:
            var_bdy=varname+'_south'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,:,0,:]
            var_bdy=varname+'_north'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,:,-1,:]
            var_bdy=varname+'_west'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,:,:,0]
            var_bdy=varname+'_east'
            ds_bdy[var_bdy].values=self.ds_smp[varname].values[:,:,:,-1]
        # deal with time vars 
        time_offset=time_frm - const.BASE_TIME
        for var in const.BDY_TIME_VAR:
            var_time=ds_bdy[var+'_time']
            # pkg time
            var_time.values[:]=int(time_offset.total_seconds())*const.S2NS+const.HALF_DAY
            #ds_bdy=ds_bdy.assign_coords({var:var_time})

        bdyfn=os.path.join(
            self.proj_root,'coawst_bdy_%s.nc'% time_frm.strftime('%Y%m%d'))
        ds_bdy.to_netcdf(bdyfn)


  