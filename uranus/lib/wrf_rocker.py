#/usr/bin/env python3
"""Build WRF preprocess workflow"""
import os, shutil, subprocess
import pandas as pd
import xarray as xr
import numpy as np
from scipy.io import FortranFile
from . import utils, io, const, mathlib

# ---Module regime consts and variables---
print_prefix='lib.WRFController>>'


# ---Classes and Functions---
class WRFRocker:
    def __init__(self, uranus):
        self.uranus=uranus
        self.cfg=uranus.cfg
        self.mach_name=self.cfg['URANUS']['machine_name']
        wrfcfg=self.cfg['WRF']
        self.rewrite_geog=wrfcfg.getboolean('rewrite_geo_em')
        self.rewrite_namelist=wrfcfg.getboolean('rewrite_namelist')
        self.drv_type=wrfcfg['drv_type']
        self.drv_dic=const.DRV_DIC[self.drv_type+'_wrf']
        self.run_ungrib=wrfcfg.getboolean('run_ungrib')
        self.run_metgrid=wrfcfg.getboolean('run_metgrid')
        self.run_real=wrfcfg.getboolean('run_real')
        self.run_wrf=wrfcfg.getboolean('run_wrf')
        self.drv_root=utils.valid_path(
            utils.parse_fmt_timepath(self.uranus.sim_strt_time, wrfcfg['drv_root']))
        self.down_drv=wrfcfg.getboolean('down_drv_data')
        self.wps_root, self.wrf_root=utils.valid_path(wrfcfg['wps_root']), utils.valid_path(wrfcfg['wrf_root'])
        self.ntasks_wrf=uranus.ntasks_atm
        self.mach_meta=uranus.machine_dic
        # read drv meta
        resource_path = os.path.join('db', self.drv_type+'.csv')
        if not(self.drv_dic['use_ungrib']):
            self.df_meta=pd.read_csv(utils.fetch_pkgdata(resource_path))
 
        
        utils.write_log(f'{print_prefix}WRFMaker Initiation Done.')
    def rock(self):
        self.preprocess()
        if self.run_ungrib:
            self.ungrib()
        if self.run_metgrid:
            self.metgrid()
        if self.run_real:
            self.real()
        if self.run_wrf:
            self.wrf()
    def preprocess(self):
        
        self.clean_workspace()
        # down driven data 
        self.downdrv()
        
        if self.rewrite_geog:
            domfn=os.path.join(
                self.uranus.domdb_root, self.uranus.nml_temp, 'geo_em*')
            io.copy_files(domfn, self.wps_root)
        
        if self.rewrite_namelist:
            # WPS: deal with namelist.wps
            nml_src=os.path.join(
                self.uranus.cfgdb_root, self.uranus.nml_temp, 'namelist.wps')
            nml_dest=os.path.join(self.wps_root, 'namelist.wps')
            
            shutil.copy(nml_src, nml_dest)
            utils.sed_wrf_timeline('start_date',self.uranus.sim_strt_time,nml_dest)
            utils.sed_wrf_timeline('end_date',self.uranus.sim_end_time,nml_dest)
            interval=int(self.drv_dic['atm_nfrq'][0])*3600
            utils.sedline('interval_seconds',f'interval_seconds = {interval}',nml_dest) 
            prefix=self.drv_type.upper()
            utils.sedline('prefix',f"prefix = '{prefix}'",nml_dest) 
            utils.sedline('fg_name',f"fg_name = '{self.drv_dic['fg_name']}'",nml_dest) 
            
            # WRF: deal with namelist.input
            nml_src=os.path.join(
                self.uranus.cfgdb_root, self.uranus.nml_temp, 'namelist.input')
            nml_dest=os.path.join(self.wrf_root, 'namelist.input')
            shutil.copy(nml_src, nml_dest)
            start_time,end_time=self.uranus.sim_strt_time,self.uranus.sim_end_time
            sim_hrs=self.uranus.run_hours
            if self.drv_dic['use_ungrib']:
                nsoil=self.drv_dic['nsoil']
            else:
                nsoil = (self.df_meta['aim_v'].str.startswith('SM')).sum()
            nplv=int(self.drv_dic['plv'][2:])+1
            sed_dic={
                'run_hours':sim_hrs, 'interval_seconds':interval,
                'num_metgrid_soil_levels':nsoil, 'num_metgrid_levels':nplv,}
            for key, itm in sed_dic.items():
                utils.sedline(key, f'{key} = {itm}',nml_dest)
            
            utils.sed_wrf_timeline('start_year',start_time,nml_dest,fmt='%Y')
            utils.sed_wrf_timeline('start_month',start_time,nml_dest,fmt='%m')
            utils.sed_wrf_timeline('start_day',start_time,nml_dest,fmt='%d')
            utils.sed_wrf_timeline('start_hour',start_time,nml_dest,fmt='%H')
             
            utils.sed_wrf_timeline('end_year',end_time,nml_dest,fmt='%Y')
            utils.sed_wrf_timeline('end_month',end_time,nml_dest,fmt='%m')
            utils.sed_wrf_timeline('end_day',end_time,nml_dest,fmt='%d')
            utils.sed_wrf_timeline('end_hour',end_time,nml_dest,fmt='%H')
    def clean_workspace(self):
        if self.run_ungrib:
            io.del_files(self.wps_root, const.UNGRIB_CLEAN_LIST)
        
        if self.run_metgrid:
            io.del_files(self.wps_root, const.METGRID_CLEAN_LIST)
        
        if self.run_real:
            io.del_files(self.wrf_root, const.WRF_CLEAN_LIST)
            
    def downdrv(self):
        if self.down_drv:
            self.uranus.crawler.down_atm()
            
    def ungrib(self):
        drv_dic=self.drv_dic
        if drv_dic['use_ungrib']:
            nml_dest=os.path.join(self.wps_root, 'namelist.wps')
            utils.write_log(
                f'{self.drv_type} use ungrib.exe to generate wrf interim files')
            # link Vtable
            io.symlink_files(
                os.path.join(
                    self.wps_root,'ungrib','Variable_Tables',drv_dic['Vtable']), 
                os.path.join(self.wps_root,'Vtable'))
            for fp,prefix in zip(drv_dic['file_patterns'],drv_dic['ungrib_prefixes']):
                # ./link_grib.csh
                day_series=pd.date_range(
                    start=self.uranus.sim_strt_time, end=self.uranus.sim_end_time, freq='D')
                
                fp_single=day_series[0].strftime(fp)
                fnlist=[os.path.join('$DRV_ROOT',fp_single)]
                for day in day_series:
                    if day.strftime(fp)!=fp_single:
                        fp_single=day.strftime(fp)
                        fnlist.append(os.path.join('$DRV_ROOT',fp_single))
                cmd=f'cd {self.wps_root}; DRV_ROOT={self.drv_root}; ./link_grib.csh {" ".join(fnlist)};'
                cmd=f'{cmd}source {self.mach_meta["bashrc"]};./ungrib.exe'
                utils.sedline('prefix',f'prefix = {prefix}',nml_dest) 
                utils.write_log(print_prefix+'Run ungrib.exe: '+cmd)
                rcode=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        else:
            utils.write_log(
                f'{self.drv_type} use struct BYTEFLOW toolkit to generate wrf interim files')
            self._build_meta()
            self._gen_interim()
    def metgrid(self):
        mach_meta=self.mach_meta
        bashrc=mach_meta['bashrc']
        mpicmd=mach_meta['mpicmd']
        metgrid_np=mach_meta['metgrid_np']
        nml_dest=os.path.join(self.wps_root, 'namelist.wps')
        utils.sedline('fg_name',f'fg_name = {self.drv_dic["fg_name"]}',nml_dest)
        cmd=utils.build_execmd(
            self.mach_name, bashrc, self.wps_root, mpicmd, metgrid_np, 'metgrid.exe')
        utils.write_log(print_prefix+'Run metgrid.exe: '+cmd)
        rcode=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        
        # special for pird
        if self.mach_name in ['pird']:        
            jobid=rcode.stdout.decode().split()[3]
            chck_cmd=mach_meta['chck_cmd']
            io.hpc_quechck(chck_cmd, jobid)       
    def real(self):
        mach_meta=self.mach_meta
        bashrc=mach_meta['bashrc']
        mpicmd=mach_meta['mpicmd']
        real_np=mach_meta['real_np']
        io.symlink_files(os.path.join(self.wps_root,'met_em.d*'), self.wrf_root)
        cmd=utils.build_execmd(
            self.mach_name, bashrc, self.wrf_root, mpicmd, real_np, 'real.exe')
        utils.write_log(print_prefix+'Run real.exe: '+cmd)
        rcode=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        if self.mach_name in ['pird']:        
            jobid=rcode.stdout.decode().split()[3]
            chck_cmd=mach_meta['chck_cmd']
            io.hpc_quechck(chck_cmd, jobid)
           
    def wrf(self):
        mach_meta=self.mach_meta
        bashrc=mach_meta['bashrc']
        mpicmd=mach_meta['mpicmd']
        wrf_np=self.ntasks_wrf
        cmd=utils.build_execmd(
            self.mach_name, bashrc, self.wrf_root, mpicmd, wrf_np, 'wrf.exe')
        utils.write_log(print_prefix+'Run wrf.exe: '+cmd)
        subprocess.run(cmd, shell=True)
    
    
    def _build_meta(self):
        '''
        Build metadata for drv types which use struct BYTEFLOW
        '''
        uranus=self.uranus
        drv_dic=self.drv_dic 
        # specific configs for drv types
        
       
        # build timefrms
        init_time=uranus.sim_strt_time
        end_time=uranus.sim_end_time
        self.frm_time_series=pd.date_range(
            start=init_time, end=end_time, freq=drv_dic['atm_nfrq'])
        
        # build file names
        special=''
        if self.drv_type=='bcmm':
            special=uranus.cfg['bcmm']['scenario_name']
        
       
        # atm    
        self.atmfn_lst, self.file_time_series=io.gen_patternfn_lst(
            self.drv_root, drv_dic, init_time, end_time, 
            kw='atm',special=special)
        # lnd 
        if self.drv_type in ['cpsv3','bcmm']:
            self.lnd_root=utils.parse_fmt_timepath(self.uranus.sim_strt_time, self.cfg[self.drv_type]['lnd_root'])
            self.lndfn_lst, _=io.gen_patternfn_lst(
                self.lnd_root, drv_dic, init_time, end_time, 
                kw='lnd',special=special)
    
   
    # Below for BYTEFLOW toolkits to generate interim files
    def _gen_interim(self):
        
        self.out_slab=io.gen_wrf_mid_template(self.drv_type+'_wrf')

        for it, tf in enumerate(self.frm_time_series):
            if tf in self.file_time_series:
                self._load_raw(tf)
            self._parse_raw(it)
            self._org_wrfinterm(tf)
            
    def _parse_raw(self,it):
        tf = self.frm_time_series[it]
        df_meta=self.df_meta
        ds,drv_dic=self.ds,self.drv_dic
        LATS,LONS,PNAME=self.drv_dic['lats'],self.drv_dic['lons'],self.drv_dic['plv']
        PLVS=const.PLV_DIC[PNAME]
        self.plvs=PLVS 
        # frm in file
        itf=it % drv_dic['frm_per_file']
        if drv_dic['vcoord']=='sigmap':
            ps=io.sel_frm(ds[drv_dic['psname']], tf, itf)
        for idy, itm in df_meta.iterrows():
            src_v, aim_v=itm['src_v'], itm['aim_v']
            # skip comments
            if src_v.startswith('#'):
                continue
            lvltype=itm['type']
            lvlmark=itm['lvlmark']
            utils.write_log(
                print_prefix+'Parsing '+src_v+',lvltype='+lvltype+',lvlmark='+lvlmark)
            if lvltype == '2d-soil':
                da=io.sel_frm(self.ds_lnd[src_v], tf, itf)
                if aim_v.startswith('SM'):
                    if  itm['units']=='kg/m-3':
                        da.values=da.values*1e-2 # m^3/m-3
                    da=accum_soil_moist(
                        da, aim_v, self.drv_dic['soillv'], self.drv_dic['soil_dim_name'])
                elif aim_v.startswith('ST'):
                    da=avg_soil_temp(
                        da, aim_v, self.drv_dic['soillv'], self.drv_dic['soil_dim_name'])
                    da=xr.where(da>0, da, np.nan)
                #else: # land sea mask
                #    da=da[0,:,:]
                #    da=xr.where(da>0, 1, 0)
            else:
                da=io.sel_frm(self.ds[src_v], tf, itf)
            if lvltype=='3d':
                #continue # for test
                if 'lev' in da.coords:
                    da = da.rename({'lev': 'plev'})
 
                if lvlmark == 'Lev' and drv_dic['vcoord']=='sigmap':
                    # interpolate from hybrid to pressure level 
                    da=mathlib.hybrid2pressure(da,self.ap,self.b, ps, PLVS)
                try: 
                    self.outfrm[aim_v+'3D']=da.interp(lat=LATS, lon=LONS, plev=PLVS,
                        method='linear',kwargs={"fill_value": "extrapolate"})
                except ValueError:
                    self.outfrm[aim_v+'3D']=da.interp(lat=LATS, lon=LONS, lev=PLVS,
                        method='linear',kwargs={"fill_value": "extrapolate"})
                    
            elif lvltype in ['2d', '2d-soil']:
                da=da.interpolate_na(
                    dim="lon", method="nearest",fill_value="extrapolate")    
                self.outfrm[aim_v]=da.interp(lat=LATS, lon=LONS,
                    method='linear',kwargs={"fill_value": "extrapolate"})
                
        # for test
        #for itm in self.outfrm.keys():
        #    self.outfrm[itm].to_netcdf(itm+'.nc')
    def _load_raw(self,tf):
        '''
        Load raw netCDF data from certain drv type 
        ''' 
        # init empty raw container dict
        self.ds, self.outfrm={},{}
        fn_lst=self.atmfn_lst
        idx=self.file_time_series.get_loc(tf)
            
        utils.write_log(print_prefix+'Loading '+fn_lst[idx])
        self.ds=xr.open_dataset(fn_lst[idx]) 
        if idx==0 and self.drv_dic['vcoord']=='sigmap':
            aname,bname=self.drv_dic['acoef'],self.drv_dic['bcoef']
            acoef=self.ds[aname].values
            p0=self.ds[self.drv_dic['p0name']].values
            self.ap=acoef*p0
            self.b=self.ds[bname].values
        if self.drv_type in ['cpsv3','bcmm']:
            self.ds_lnd=xr.open_dataset(self.lndfn_lst[idx])
        
    
    def _org_wrfinterm(self, tf, tgt='main'):
        ymdH=tf.strftime('%Y-%m-%d_%H')
        df_meta=self.df_meta
        PNAME=self.drv_dic['plv']
        PLVS=const.PLV_DIC[PNAME]
        
        if tgt=='main':
            out_fn=os.path.join(self.wps_root,f'{self.drv_type.upper()}:{ymdH}')
        if tgt=='sst':
            out_fn=os.path.join(self.wps_root,f'{self.drv_type.upper()}_SST:{ymdH}')
        
        utils.write_log(print_prefix+'Writing '+out_fn)
        
        # dtype='>u4' for header (big-endian, unsigned int)
        wrf_mid = FortranFile(out_fn, 'w', header_dtype=np.dtype('>u4'))
        
        out_dic=self.out_slab
        out_dic['HDATE']=tf.strftime('%Y-%m-%d_%H:%M:%S:0000')
        
            
        for idy, itm in df_meta.iterrows():
            src_v, aim_v=itm['src_v'], itm['aim_v']
            if src_v.startswith('#'):
                continue
            lvltype=itm['type']
            out_dic['FIELD']=aim_v
            out_dic['UNIT']=itm['units']
            out_dic['DESC']=itm['desc']
            out_dic['XLVL']=200100.0
            if lvltype=='3d':
               for lvl in PLVS:
                    out_dic['XLVL']=lvl
                    out_dic['SLAB']=self.outfrm[aim_v+'3D'].sel(plev=lvl).values
                    io.write_record(wrf_mid, out_dic)
            elif lvltype=='2d' or lvltype=='2d-soil':
                out_dic['SLAB']=self.outfrm[aim_v].values
            io.write_record(wrf_mid, out_dic)
        wrf_mid.close()        

         
def accum_soil_moist(da, aim_v, model_name, lvname):
    strt_dp, end_dp=utils.decode_depth(aim_v)
    SOI_LVS=const.SOILLV_DIC[model_name]
    DP_SOI_LVS=np.diff(SOI_LVS,prepend=0.0)
    ids, ide, s_res, e_res=mathlib.find_indices(strt_dp, end_dp, SOI_LVS)
    if model_name=='cpsv3':
        # from kg m-2 to m3 m-2
        # accum
        da_r=da[ids:ide,:,:].sum(dim=lvname)    
        da_r=da_r-da[ids]*s_res/DP_SOI_LVS[ids]-da[ide]*e_res/DP_SOI_LVS[ide]
        da_r=(da_r/const.RHO_WATER)/(end_dp-strt_dp)
    elif model_name=='bcmm':
        da_r=da[ids,:,:]
    return da_r

def avg_soil_temp(da, aim_v, model_name, lvname):
    strt_dp, end_dp=utils.decode_depth(aim_v)
    SOI_LVS=const.SOILLV_DIC[model_name]
    DP_SOI_LVS=np.diff(SOI_LVS) # DP_SOI_LVS[0]=SOI_LVS[1]-SOI_LVS[0]
    ids, ide, s_res, e_res=mathlib.find_indices(strt_dp, end_dp, SOI_LVS)
    da_r=da[ids:ide,:,:]
    DP_SOI_LVS[ids], DP_SOI_LVS[ide-1]=DP_SOI_LVS[ids]-s_res, DP_SOI_LVS[ide-1]-e_res
    for idx in range(ids, ide):
        da_r[idx-ids,:,:]=da_r[idx-ids,:,:]*DP_SOI_LVS[idx]/(end_dp-strt_dp)
    da_r=da_r.sum(dim=lvname)
    return da_r
