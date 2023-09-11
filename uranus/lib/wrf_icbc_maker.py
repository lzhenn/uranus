#/usr/bin/env python3
"""Build WRF preprocess workflow"""
import os, shutil, cftime
import pandas as pd
import xarray as xr
import numpy as np
from . import utils, io, const, mathlib

# ---Module regime consts and variables---
print_prefix='lib.WRFMaker>>'

class WRFMaker:
    def __init__(self, uranus):
        self.uranus=uranus
        self.cfg=uranus.cfg
        wrfcfg=self.cfg['WRF']
        self.run_maker=wrfcfg.getboolean('preprocess_wrf')
        self.drv_type=wrfcfg['drv_type']
        self.drv_dic=const.DRV_DIC[self.drv_type]
        self.run_ungrib=wrfcfg.getboolean('run_ungrib')
        self.run_metgrid=wrfcfg.getboolean('run_metgrid')
        self.run_real=wrfcfg.getboolean('run_real')
        self.drv_root=utils.valid_path(wrfcfg['drv_root'])
        self.wps_root, self.wrf_root=utils.valid_path(wrfcfg['wps_root']), utils.valid_path(wrfcfg['wrf_root'])
        utils.write_log(f'{print_prefix}WRFMaker Initiation Done.')
    
    def make_icbc(self):
        if self.run_maker:
            self.preprocess()
            if self.run_ungrib:
                self.ungrib()
            if self.run_metgrid:
                self.metgrid()
            if self.run_real:
                self.real()
        else:
            utils.write_log('run_maker is False, No need to make wps')
            return
        
    # ---Classes and Functions---
    def preprocess(self):
        self.clean_workspace()
    
    def clean_workspace(self):
        if self.run_ungrib or self.run_metgrid:
            io.del_files(self.wps_root, const.WPS_CLEAN_LIST)
        
        if self.run_real:
            io.del_files(self.wrf_root, const.WRF_CLEAN_LIST)
            
    
    def ungrib(self):
        drv_dic=self.drv_dic
        if drv_dic['use_ungrib']:
            pass
        else:
            utils.write_log(f'{self.drv_type} use struct BYTEFLOW toolkit to generate wrf interim files')
            self._build_meta()
            self._gen_interim()
        for time_frm in self.frm_time_series:
            pass 
    def metgrid(self):
        pass
    
    def real(self):
        pass
    
    
    
    def _build_meta(self):
        '''
        Build metadata for drv types which use struct BYTEFLOW
        '''
        uranus=self.uranus
        drv_dic=self.drv_dic 
        # specific configs for drv types
        
        # read drv meta
        resource_path = os.path.join('db', self.drv_type+'.csv')
        self.df_meta=pd.read_csv(utils.fetch_pkgdata(resource_path))
        
        # build timefrms and file names
        init_time=uranus.sim_strt_time
        end_time=uranus.sim_end_time
        self.frm_time_series=pd.date_range(
            start=init_time, end=end_time, freq=drv_dic['atm_nfrq'])
       
        self.file_time_series=pd.date_range(
            start=init_time, end=end_time, freq=drv_dic['atm_file_nfrq'])
        
        self.atmfn_lst=io.gen_patternfn_lst(
            self.drv_root, drv_dic, init_time, end_time)
        if self.drv_type=='cpsv3':
            self.lnd_root=self.cfg['cpsv3']['lnd_root']
            self.lndfn_lst=io.gen_patternfn_lst(
                self.lnd_root, drv_dic, init_time, end_time, kw='lnd')
    def _gen_interim(self):
        for it, tf in enumerate(self.frm_time_series):
            if tf in self.file_time_series:
                self._load_raw(tf)
            self._parse_raw(it)
            
    def _parse_raw(self,it):
        tf = self.frm_time_series[it]
        df_meta=self.df_meta
        ds,drv_dic=self.ds,self.drv_dic
        LATS,LONS,PNAME=self.drv_dic['lats'],self.drv_dic['lons'],self.drv_dic['plv']
        PLVS=const.PLV_DIC[PNAME]
        
        # frm in file
        itf=it % drv_dic['frm_per_file']
        if drv_dic['vcoord']=='sigmap':
            ps=io.sel_frm(ds[drv_dic['psname']], tf, itf)
        for idy, itm in df_meta.iterrows():
            src_v, aim_v=itm['src_v'], itm['aim_v']
            lvltype=itm['type']
            lvlmark=itm['lvlmark']
            utils.write_log(
                print_prefix+'Parsing '+src_v+',lvltype='+lvltype+',lvlmark='+lvlmark)
            if lvltype == '2d-soil' and self.drv_type=='cpsv3':
                da=io.sel_frm(self.ds_lnd[src_v], tf, itf)
                if aim_v.startswith('SM'):
                    da=accum_soil_moist(
                        da, aim_v, self.drv_dic['soillv'], self.drv_dic['soil_dim_name'])
                elif aim_v.startswith('ST'):
                    da=avg_soil_temp(da, aim_v, self.drv_dic['soillv'], self.drv_dic['soil_dim_name'])
                else: # land sea mask
                    da=da[0,:,:]
                    da=xr.where(da>0, 1, 0)
            else:
                da=io.sel_frm(self.ds[src_v], tf, itf)
            if lvltype=='3d':
                #continue
                if lvlmark == 'Lev' and drv_dic['vcoord']=='sigmap':
                    # interpolate from hybrid to pressure level 
                    da=mathlib.hybrid2pressure(da,self.ap,self.b, ps, PLVS)
                self.outfrm[src_v]=da.interp(lat=LATS, lon=LONS, plev=PLVS,
                        method='linear',kwargs={"fill_value": "extrapolate"})
            
            elif lvltype in ['2d', '2d-soil']:
                #da=da.interpolate_na(
                #    dim="lon", method="nearest",fill_value="extrapolate")    
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
            self.ap=self.ds[aname].values
            self.b=self.ds[bname].values
        
        if self.drv_type=='cpsv3':
            self.ds_lnd=xr.open_dataset(self.lndfn_lst[idx])
            
            
def accum_soil_moist(da, aim_v, model_name, lvname):
    strt_dp, end_dp=utils.decode_depth(aim_v)
    SOI_LVS=const.SOILLV_DIC[model_name]
    DP_SOI_LVS=np.diff(SOI_LVS)
    ids, ide, s_res, e_res=mathlib.find_indices(strt_dp, end_dp, SOI_LVS)
    # from kg m-2 to m3 m-2
    # accum
    da_r=da[ids:ide,:,:].sum(dim=lvname)    
    da_r=da_r-da[ids]*s_res/DP_SOI_LVS[ids]-da[ide]*e_res/DP_SOI_LVS[ide]
    da_r=(da_r/const.RHO_WATER)/(end_dp-strt_dp)
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
