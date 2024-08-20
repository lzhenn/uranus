#/usr/bin/env python3
"""Build ROMS workflow"""
import os, datetime
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
        try:
            self.version=self.cfg['URANUS']['cwst_version']
        except:
            self.version='3.8'
            utils.write_log(print_prefix+'no version in cfg, use 3.8')
        romscfg=self.cfg['ROMS']
        self.gen_icbc=romscfg.getboolean('gen_roms_icbc')
        self.down_drv=romscfg.getboolean('down_drv_data')
        self.restart_run=romscfg.getboolean('restart_run')
        self.strt_time,self.end_time=uranus.sim_strt_time,uranus.sim_end_time
        self.domid=uranus.domain_lv 
        if 'pseudo_init_ts' in romscfg:
            self.pseudo_init=True
            self.pseudo_init_time=utils.parse_init_time(romscfg['pseudo_init_ts'])
            self.pseudo_delta=self.pseudo_init_time-self.strt_time
            self.pseudo_end_time=self.pseudo_init_time+datetime.timedelta(hours=uranus.run_hours)
        else:
            self.pseudo_init=False
            self.pseudo_delta=datetime.timedelta(0)
        time_offset=self.strt_time - const.BASE_TIME
        self.dstart=time_offset.total_seconds()/const.DAY_IN_SEC

        self.drv_type=romscfg['drv_type']
        self.gen_surf=romscfg.getboolean('gen_surf')
        self.surf_type=romscfg['surf_type']
        self.surf_frq=romscfg['surf_frq']
        self.surf_wrfdom=romscfg['wrf_match']
        self.surf_root=utils.parse_fmt_timepath(self.strt_time, romscfg['surf_root'])
        
        if self.cfg.has_section(uranus.mode):
            if self.cfg.has_option(uranus.mode, 'tc_intense_factor'):
                self.tc_intense_factor=self.cfg.getfloat(uranus.mode, 'tc_intense_factor')
                utils.write_log(f'{print_prefix}TC INTENSIFICATION TURNED ON, factor={self.tc_intense_factor:.2f}')
                try:
                    self.sim_trck=pd.read_csv(
                        f'{self.surf_root}/tc_track.csv', parse_dates=['time'],index_col=['time'])
                    utils.write_log(f'{print_prefix}TC INTENSIFICATION TRACK LOADED')
                except FileNotFoundError:
                    utils.throw_error(f'{self.surf_root}/tc_track.csv not found')
        self.drv_dic=const.DRV_DIC[self.drv_type+'_roms']
        if self.drv_type=='roms':
            self.drv_dic['ocn_nfrq']=romscfg['drv_frq']
            self.drv_dic['ocn_file_nfrq']=romscfg['drv_frq']
        self.drv_root=utils.parse_fmt_timepath(self.strt_time, romscfg['drv_root'])
        self.proj_root=uranus.proj_root
        self.domdb_root=uranus.domdb_root
        self.dt=float(romscfg['dt'])
        self.nrst=int(romscfg['nrst'])
        self.nhis=int(romscfg['nhis'])
        utils.write_log(f'{print_prefix}ROMSMaker Initiation Done.')
    
    def build_meta(self):
        # build meta
        self.frm_time_series=pd.date_range(
            start=self.strt_time, end=self.end_time, freq=self.drv_dic['ocn_nfrq'])
        if len(self.frm_time_series)>1:
            self.ocn_nfrq=int((self.frm_time_series[1]-self.frm_time_series[0]).total_seconds())
        if self.drv_type=='roms':
            if self.domid=='d02':
                bdydomid='d01'
            elif self.domid=='d03':
                bdydomid='d01'
            self.ocnfn_lst, self.file_time_series=io.gen_roms_his_fnlst(
                self.drv_root, self.drv_dic, bdydomid, self.frm_time_series, self.ocn_nfrq)
        else:
            # paeudo init for cases before the reanalysis covered
            if self.pseudo_init:
                self.ocnfn_lst, self.file_time_series=io.gen_patternfn_lst(
                    self.drv_root, self.drv_dic, 
                    self.pseudo_init_time, self.pseudo_end_time, kw='ocn')
            else:
                self.ocnfn_lst, self.file_time_series=io.gen_patternfn_lst(
                    self.drv_root, self.drv_dic, self.strt_time, self.end_time, kw='ocn')
        # sanity check
        stat=io.check_filelist(self.ocnfn_lst)
        
        if not(stat):
            if self.uranus.fcst_flag:
                for buf_day in range(1,3):
                    utils.write_log(
                        f'{print_prefix}{self.drv_root} missing ocn file(s), try {buf_day} days back', 30)
                    buf_time=self.strt_time+datetime.timedelta(days=-buf_day)
                    self.drv_root=utils.parse_fmt_timepath(buf_time, self.cfg['ROMS']['drv_root'])
                    self.ocnfn_lst, self.file_time_series=io.gen_patternfn_lst(
                        self.drv_root, self.drv_dic, self.strt_time, self.end_time, kw='ocn')
                    stat=io.check_filelist(self.ocnfn_lst)
                    if stat:    
                        break
                if not(stat):
                    utils.throw_error(f'{print_prefix}{self.drv_root} missing ocn file(s)')
            else:
                utils.throw_error(f'{print_prefix}{self.drv_root} missing ocn file(s)')
    def rock(self):
        if self.gen_icbc:
        # down driven data 
            if self.down_drv:
                self.downdrv()
            self.build_meta()
            self.build_forc() 
            self.build_icbc() 
    def downdrv(self):
        if self.down_drv:
            self.uranus.crawler.down_ocn()
    def prepare_cplrock(self):
        nml_temp=self.uranus.nml_temp
        
 
        domfn=os.path.join(
                self.domdb_root, nml_temp, '*omp.nc')
        io.symlink_files(domfn, self.proj_root)
            
        roms_in=os.path.join(self.proj_root, f'roms_{self.domid}.in')
        run_seconds=int((self.end_time-self.strt_time).total_seconds())
        ntimes=run_seconds//self.dt
        nrst=self.nrst*const.DAY_IN_SEC//self.dt
        nhis=self.nhis*const.MIN_IN_SEC//self.dt 
        # only use half day
        ndefhis=nhis
        #ndefhis=const.DAY_IN_SEC//self.dt//24
        
        # use relative path for infile modification 
        proj_reldir=os.path.join('Projects', nml_temp)
        grdname=os.path.join(proj_reldir, f'roms_{self.domid}_omp.nc')
        if self.restart_run:
            ininame= f'roms_rst_{self.domid}.nc'
            nrrec='-1'
        else:
            ininame=os.path.join(proj_reldir, f'roms_ini_{self.domid}.nc')
            nrrec='0'
        
        '''
        for tf in self.frm_time_series[:-1]: 
            tf_str=tf.strftime('%Y%m%d%H%M')
            clmname+=os.path.join(proj_reldir, f'roms_clm_{self.domid}_{tf_str}.nc |\n')
            bryname+=os.path.join(proj_reldir, f'roms_bdy_{self.domid}_{tf_str}.nc |\n')
        tf_close=self.frm_time_series[-1]
        tf_str=tf_close.strftime('%Y%m%d%H%M')
        if self.version=='3.5': # bug fix for 3.5
            clmname+=os.path.join(proj_reldir, f'roms_clm_{self.domid}_{tf_str}.nc |\n')
            bryname+=os.path.join(proj_reldir, f'roms_bdy_{self.domid}_{tf_str}.nc |\n')
            tf_close=tf_close+datetime.timedelta(seconds=self.ocn_nfrq)
            tf_str=tf_close.strftime('%Y%m%d%H%M')
            clmname+=os.path.join(proj_reldir, f'roms_clm_{self.domid}_{tf_str}.nc')
            bryname+=os.path.join(proj_reldir, f'roms_bdy_{self.domid}_{tf_str}.nc')
        else:    
            clmname+=os.path.join(proj_reldir, f'roms_clm_{self.domid}_{tf_str}.nc')
            bryname+=os.path.join(proj_reldir, f'roms_bdy_{self.domid}_{tf_str}.nc')
        ''' 
        
        clmname=os.path.join(proj_reldir, f'roms_clm_{self.domid}.nc')
        bryname=os.path.join(proj_reldir, f'roms_bdy_{self.domid}.nc')
        frcname=os.path.join(proj_reldir, f'roms_forc_{self.domid}.nc') 
        
        #tf_start=self.strt_time       
        #tf_str=tf_start.strftime('%Y%m%d')
        sed_dic={
            'NtileI':self.uranus.ntasks_iocn, 'NtileJ':self.uranus.ntasks_jocn,
            'NTIMES':ntimes, 'DT':f"{self.dt:.1f}d0", 'NRST':nrst, 'NRREC':nrrec,
            'NHIS':nhis, 
            'NDEFHIS':ndefhis, 'NDIA':nhis, #'NAVG':nhis, 
            'GRDNAME':grdname, 'ININAME':ininame, 'CLMNAME':clmname, 'BRYNAME':bryname,
            'DSTART':f'{self.dstart:.3f}d0',
            #'TIME_REF':f"{tf_str}.0d0"
        }
        if self.uranus.active_comp[0]==0:
            sed_dic['FRCNAME']=frcname
        for key, itm in sed_dic.items():
            utils.sedline(key, f'{key} == {itm}', roms_in, count=1)

    
    def load_domain(self):
        """ load domain file """
        utils.write_log(print_prefix+'Load domain file...')
        domfn=os.path.join(
            self.domdb_root, self.uranus.nml_temp, f'roms_{self.domid}_omp.nc')
        self.ds_static=xr.load_dataset(domfn)
        ds_static=self.ds_static
        self.mask=ds_static['mask_rho'].values
        self.h=ds_static['h'].values
        self.lat1d,self.lon1d=ds_static['lat_rho'][:,0].values,ds_static['lon_rho'][0,:].values
        self.lat_u, self.lon_u=ds_static['lat_u'][:,0].values, ds_static['lon_u'][0,:].values
        self.lat_v, self.lon_v=ds_static['lat_v'][:,0].values, ds_static['lon_v'][0,:].values
       
        # load sample file
        sampfn=os.path.join(
            self.uranus.domdb_root, self.uranus.nml_temp, f'roms_{self.domid}_inismp.nc')
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
    
    def build_forc(self):
        """ build forcing for ROMS """
        if not(self.gen_surf):
            return
        utils.write_log(
            print_prefix+'build forcing from %s to %s...'%(
                self.strt_time.strftime('%Y%m%d%H'),self.end_time.strftime('%Y%m%d%H')))
        
        surf_dir=self.surf_root
        
        if not os.path.exists(surf_dir):
            utils.throw_error(print_prefix+'Surface forcing directory does not exist!')
        
        forc_time_series=pd.date_range(
            start=self.strt_time, end=self.end_time, freq=self.surf_frq)
        trange=[(it-const.BASE_TIME).total_seconds()/const.DAY_IN_SEC for it in forc_time_series]
        curr_time=self.strt_time
        
        if self.surf_type=='wrf':
            import netCDF4 as nc4
            import wrf
            surf_file=io.get_wrf_fn(curr_time, f'{self.surf_wrfdom}')
            surf_file=os.path.join(surf_dir,surf_file)
            wrf_hdl=nc4.Dataset(surf_file)
            XLAT=wrf.getvar(wrf_hdl,'XLAT')
            XLONG=wrf.getvar(wrf_hdl,'XLONG')
            ds_forc=io.gen_roms_forc_template(trange, XLAT, XLONG)
            wrf_hdl.close()
            # iter timeframes 
            for idx,curr_time in enumerate(forc_time_series):
                surf_file=io.get_wrf_fn(curr_time, f'{self.surf_wrfdom}')
                utils.write_log(print_prefix+'Read surface forcing from '+surf_file)
                surf_file=os.path.join(surf_dir,surf_file)
                wrf_hdl=nc4.Dataset(surf_file)
                
                if hasattr(self, 'sim_trck'):
                    try:
                        curr_loc=self.sim_trck.index.get_loc(curr_time) 
                        curr_rec=self.sim_trck.iloc[curr_loc]
                        utils.write_log(
                            print_prefix+f'TC_INTENSIFY: TC centered at lat:{curr_rec["lat"]:.3f}, lon:{curr_rec["lon"]:.3f}')
                    except KeyError:
                        utils.write_log(
                            print_prefix+'TC_INTENSIFY: NO TC record for '+curr_time.strftime('%Y%m%d%H'))
                        pass 
                for roms_var, wrf_var in const.ROMS_WRF_FORC_MAPS.items():
                    if roms_var=='Tair':
                        temp_var = wrf.getvar(
                            wrf_hdl, wrf_var, 
                            timeidx=wrf.ALL_TIMES, method="cat")
                        temp_var=temp_var-const.K2C
                    elif roms_var=='Pair':
                        temp_var = wrf.getvar(
                            wrf_hdl, wrf_var, 
                            timeidx=wrf.ALL_TIMES, method="cat", units='mb')
                    else:
                        temp_var = wrf.getvar(
                            wrf_hdl, wrf_var,
                            timeidx=wrf.ALL_TIMES, method="cat")
                    if ('curr_rec' in locals()) and (roms_var=='Uwind' or roms_var=='Vwind'):
                        utils.write_log(
                            print_prefix+f'TC_INTENSIFY: {roms_var} orginal range: {temp_var.min().values:.2f}, {temp_var.max().values:.2f}')
                        temp_var=mathlib.intensify_tc(
                            self.tc_intense_factor, temp_var,
                            curr_rec['lat'], curr_rec['lon'])
                        utils.write_log(
                            print_prefix+f'TC_INTENSIFY: {roms_var} intensified range: {temp_var.min().values:.2f}, {temp_var.max().values:.2f}')
                    ds_forc[roms_var].values[idx,:,:]=temp_var.values
                wrf_hdl.close()
        forc_fn=os.path.join(
            self.proj_root,f'roms_forc_{self.domid}.nc')
        ds_forc.to_netcdf(forc_fn)
    def build_icbc(self):
        """ build icbcs for ROMS """
        utils.write_log(
            print_prefix+'build icbcs from %s to %s...'%(
                self.strt_time.strftime('%Y%m%d%H'),self.end_time.strftime('%Y%m%d%H')))
        cln_lst=[ele+self.domid for ele in const.ROMS_CLEAN_LIST] 
        io.del_files(self.proj_root, cln_lst)
        # load domain file
        self.load_domain()
        
        bdysmpfn=os.path.join(
            self.domdb_root, self.uranus.nml_temp,f'roms_{self.domid}_bdysmp.nc')
        self.ds_bdy=xr.load_dataset(bdysmpfn)         
        
        clmsmpfn=os.path.join(
            self.domdb_root, self.uranus.nml_temp,f'roms_{self.domid}_clmsmp.nc')
        self.ds_clm=xr.load_dataset(clmsmpfn)
        
        self.bdy_fnlst, self.clm_fnlst=[],[]
        for it, tf in enumerate(self.frm_time_series):
            utils.write_log(
                print_prefix+'build icbcs@%s...'% tf.strftime('%Y%m%d%H%M'))
            if tf+self.pseudo_delta in self.file_time_series:
                # load raw file  
                self.load_raw(it) 
            if self.drv_type=='roms':
                self.interp_roms()
            else:
                # first deal with zeta
                roms_var='zeta'
                '''
                if self.drv_type=='cfs':
                    utils.write_log(
                        print_prefix+'cfs data set zeta=0...')
                    self.ds_smp['zeta'].values[:]=0.0
                else:
                    self.inter2d(roms_var)
                '''
                self.ds_smp['zeta'].values[:]=0.0
                # calculate vertical coordinate before 3d interpolation
                zeta = self.ds_smp['zeta'] 
                h=zeta # obtain dimensional information of zeta  
                h= self.ds_static.h.values
                z_rho = mathlib.sigma2depth(zeta, h, self.ds_smp)
                self.dz=z_rho[1:,:,:]-z_rho[:-1,:,:]
                self.dz=np.concatenate((self.dz,self.dz[-1:,:,:]),axis=0)
                self.dp_idx=mathlib.assign_depth_idx(
                    z_rho, self.mask)
                
                # Hot Spot: May Comment out 3D interp for test
                for roms_var in ['temp','salt','u','v']:
                    self.inter3d(roms_var)
                
            
            if tf==self.strt_time:
                # pkg time
                self.ds_smp['ocean_time'].values[:]=int(self.dstart*const.DAY_IN_SEC)*const.S2NS
                self.ds_smp=self.ds_smp.assign_coords({'ocean_time':self.ds_smp['ocean_time']})
                # output
                inifn=os.path.join(self.proj_root, f'roms_ini_{self.domid}.nc')
                self.ds_smp.to_netcdf(inifn)
                utils.write_log(print_prefix+'build initial conditions done!')
            if self.drv_type!='roms':
                self.build_clm(it, tf)
            self.build_bdy(it, tf)
        self.merge_icbc()
    
    def merge_icbc(self):
        '''merge all icbcs into one single file'''
        if self.drv_type!='roms':
            ds_lst=[xr.load_dataset(fn) for fn in self.clm_fnlst]
            ds_all=xr.merge(ds_lst)
            ds_all.to_netcdf(os.path.join(
                self.proj_root, 'roms_clm_%s.nc' %self.domid))
        
        ds_lst=[xr.load_dataset(fn) for fn in self.bdy_fnlst]
        ds_all=xr.merge(ds_lst)
        ds_all.to_netcdf(os.path.join(
            self.proj_root, 'roms_bdy_%s.nc' %self.domid))
 
    def load_raw(self, it):
        """ load raw GCM/ROMS files to build ICBC"""
        fn = self.ocnfn_lst[it]
        utils.write_log(print_prefix+'Load raw file: '+fn)
        self.ds_raw=xr.load_dataset(fn)
        ds_raw=self.ds_raw
        self.varname_remap()
        # deal with previous hycom
        if self.drv_type=='hycom':
            if self.varmap['lat_rho']=='Y':
                ds_raw=ds_raw.rename_dims({
                    'Depth':'depth','Y':'lat','X':'lon'})
                ds_raw=ds_raw.drop_indexes(
                    ['Depth'])
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
    def interp_roms(self):
        '''
        first fill the missing value
        then interpolate to new 2d grid 
        '''
        smpcoords={
            'lon1d_rho':self.lon1d,'lat1d_rho':self.lat1d,
            'lon1d_u':self.lon_u,'lat1d_u':self.lat_u,
            'lon1d_v':self.lon_v,'lat1d_v':self.lat_v}
        roms_vardimmap=[
            ('zeta','rho'),('temp','rho'),
            ('salt','rho'),('u','u'),('v','v'),
            ('ubar','u'),('vbar','v')]
        for (roms_varname, dimname) in roms_vardimmap:
            utils.write_log(print_prefix+roms_varname+' 2d-interp...')
            lat1d_raw=self.ds_raw[f'lat_{dimname}'][:,0].values
            lon1d_raw=self.ds_raw[f'lon_{dimname}'][0,:].values 
            raw_var=self.ds_raw[roms_varname]
            smp_var=self.ds_smp[roms_varname]
            # unify the var using 1d lat and lon
            if roms_varname in ['zeta', 'ubar','vbar']:
                da_temp= xr.DataArray(
                    raw_var.values, dims=('time', 'lat', 'lon'), 
                    coords={'time':[0],'lat': lat1d_raw, 'lon': lon1d_raw})
            else:
                ndpth=smp_var.shape[1] 
                da_temp= xr.DataArray(
                    raw_var.values[:,0:ndpth,:,:], dims=('time', 'depth', 'lat', 'lon'), 
                    coords={'time':[0],'depth': np.arange(ndpth), 'lat': lat1d_raw, 'lon': lon1d_raw})
            da_temp = da_temp.interpolate_na(
                dim='lon', method="linear", fill_value="extrapolate")
            da_temp = da_temp.interpolate_na(
                dim='lat', method="linear", fill_value="extrapolate")
            da_temp=da_temp.interp(
                {'lat':smpcoords[f'lat1d_{dimname}'],'lon':smpcoords[f'lon1d_{dimname}']},
                method='linear')
            smp_var.values =da_temp.values
            # down for debug
            #if roms_varname in ['zeta', 'ubar','vbar', 'u', 'v']:
            #    smp_var.values=smp_var.values*0.0
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
        
    def build_clm(self, idt, time_frm):
        '''Build climatology file'''
        utils.write_log(
            print_prefix+'build clm file@%s...' % time_frm.strftime('%Y%m%d%H%M'))
        ds_clm=self.ds_clm
        # loop the variables to assign values
        for varname in const.ROMS_VAR:
            ds_clm[varname].values=self.ds_smp[varname].values
 
        # deal with time vars 
        time_offset=time_frm - const.BASE_TIME
        for var in const.CLM_TIME_VAR:
            var_time=ds_clm[var+'_time']
            # pkg time
            var_time.values[:]=int(time_offset.total_seconds())*const.S2NS
            '''
            if self.version=='3.5':
                var_time.values[:]=int(time_offset.total_seconds())*const.S2NS
                #var_time.values[:]=int(time_offset.total_seconds()-1)*const.S2NS
            else:
                var_time.values[:]=int(time_offset.total_seconds())*const.S2NS
            #ds_clm=ds_clm.assign_coords({var:var_time})
            '''
        clmfn=os.path.join(
            self.proj_root,f'roms_clm_{self.domid}_%s.nc'% time_frm.strftime('%Y%m%d%H%M'))
        self.clm_fnlst.append(clmfn)
        ds_clm.to_netcdf(clmfn)
        # bug fix for 3.5
        '''
        if self.version=='3.5' and time_frm==self.frm_time_series[-1]:
            time_frm=self.frm_time_series[-1]+datetime.timedelta(
                    seconds=self.ocn_nfrq)
            for var in const.CLM_TIME_VAR:
                var_time=ds_clm[var+'_time']
                var_time.values[:]=(idt+1)*self.ocn_nfrq*const.S2NS
            clmfn=os.path.join(
                self.proj_root,'coawst_clm_%s.nc'% time_frm.strftime('%Y%m%d%H%M'))
            self.clm_fnlst.append(clmfn)
            ds_clm.to_netcdf(clmfn)
        '''
    def build_bdy(self, idt, time_frm):
        '''Build bdy file'''
        utils.write_log(
            print_prefix+'build bdy file@%s...' % time_frm.strftime('%Y%m%d%H%M'))
       
        ds_bdy=self.ds_bdy
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
        time_offset=time_frm - const.BASE_TIME
        # deal with time vars 
        for var in const.BDY_TIME_VAR:
            var_time=ds_bdy[var+'_time']
            var_time.values[:]=int(time_offset.total_seconds())*const.S2NS
            '''
            if self.version=='3.5':
                var_time.values[:]=int(time_offset.total_seconds())*const.S2NS
                #var_time.values[:]=int(time_offset.total_seconds()-1)*const.S2NS
            else:
                var_time.values[:]=int(time_offset.total_seconds())*const.S2NS
            '''
            #ds_bdy=ds_bdy.assign_coords({var:var_time})
        bdyfn=os.path.join(
            self.proj_root,f'roms_bdy_{self.domid}_%s.nc'% time_frm.strftime('%Y%m%d%H%M'))
        self.bdy_fnlst.append(bdyfn)
        ds_bdy.to_netcdf(bdyfn)
        '''
        # bug fix for 3.5
        if self.version=='3.5' and time_frm==self.frm_time_series[-1]:
            time_frm=self.frm_time_series[-1]+datetime.timedelta(
                    seconds=self.ocn_nfrq)
            for var in const.BDY_TIME_VAR:
                var_time=ds_bdy[var+'_time']
                var_time.values[:]=(idt+1)*self.ocn_nfrq*const.S2NS
            bdyfn=os.path.join(
                self.proj_root,f'roms_bdy_{self.domid}_%s.nc'% time_frm.strftime('%Y%m%d%H%M'))
            self.bdy_fnlst.append(bdyfn)
            ds_bdy.to_netcdf(bdyfn)
        ''' 
 
