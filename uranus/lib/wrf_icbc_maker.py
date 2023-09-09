#/usr/bin/env python3
"""Build WRF preprocess workflow"""
import os, shutil
import pandas as pd
from . import utils, io, const
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
        # specific configs for drv types
        
        # read drv meta
        resource_path = os.path.join('db', self.drv_type+'.csv')
        self.df_meta=pd.read_csv(utils.fetch_pkgdata(resource_path))
        
        # build timefrms and file names
        init_time=uranus.sim_strt_time
        end_time=uranus.sim_end_time
        self.frm_time_series=pd.date_range(
            start=init_time, end=end_time, freq=self.drv_dic['atm_nfrq'])
        
        self.file_time_series=pd.date_range(
            start=inittime, end=endtime, freq=drv_dic['atm_file_nfrq'])
        
        self.atmfn_lst=io.gen_patternfn_lst(
            self.drv_root, self.drv_dic, init_time, end_time)
        if self.drv_type=='cpsv3':
            self.lnd_root=self.cfg['cpsv3']['lnd_root']
            self.lndfn_lst=io.gen_patternfn_lst(
                self.lnd_root, self.drv_dic, init_time, end_time, kw='lnd')
    def _gen_interim(self):
        for tf in self.frm_time_series:
            if tf in self.file_time_series:
                self._load_raw(tf)
            self._parse_raw(tf)
    def _parse_raw(self,tf):
        pass
    def _load_raw(self,tf):
        '''
        Load raw netCDF data from certain drv type 
        ''' 
        # init empty raw container dict
        self.ds, self.outfrm={},{}
        idf=0 # file index for CMIP6 files
        for idx, irow in self.meta_rows.iterrows(): 
            # CMIP6 regular
                utils.write_log(print_prefix+'Loading '+self.fn_lst[idx])
                ds=xr.open_dataset(self.fn_lst[idf])
                idf=idf+1
                
            # need interpolation coefficients
            if (lvlmark=='Lev' and itm['type']=='3d'):
                if not(hasattr(self,'ap')):
                    self.ap=ds['ap'].values
                    self.b=ds['b'].values
                    self.ps=ds['ps'].sel(
                    time=slice(self.etl_strt_time,self.etl_end_time))
