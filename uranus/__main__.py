#!/usr/bin/env python3
'''
Date: Sep 08, 2023
uranus is the main driver  

This is the main script to drive the model

History:
Sep 08, 2023 --- Kick off the project 

L_Zealot
'''
import sys, os, subprocess
import logging, logging.config
import shutil
import datetime

from .lib import cfgparser, utils, const, io

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

print_prefix='Uranus>>'
# path to the top-level handler
CWD=sys.path[0]

# path to this module
#MWD=os.path.split(os.path.realpath(__file__))[0]
   
class Uranus:
    '''
    Uranus is a class to drive the system
    '''

    def __init__(self, cfgfn='config.case.ini'):
        self._setup_logging()

        if not(os.path.exists(os.path.join(CWD,cfgfn))):
            utils.write_log('config file not exist, copy from pkg...')
            copy_cfg(os.path.join(CWD,cfgfn))
        self.cfg=cfgparser.read_cfg(os.path.join(CWD,cfgfn))
        cfg=self.cfg['URANUS']
        
        self.rock_flg=cfg.getboolean('rock_cpl')
        self.rock_wrf=self.cfg['WRF'].getboolean('rock_wrf')
        self.rock_roms=self.cfg['ROMS'].getboolean('rock_roms')
        self.mode=cfg['uranus_mode']
        self.nml_temp=cfg['nml_temp']
        self.machine_name=cfg['machine_name']
        
        
        self.machine_dic=const.MACHINE_DIC[self.machine_name]
        self.bashrc=self.machine_dic['bashrc']
        self.mpicmd=self.machine_dic['mpicmd']
 
        self.cfgdb_root=self.machine_dic['cfgdb_root']
        self.domdb_root=self.machine_dic['domdb_root']
        if self.mode=='shu':
            self.cplexe_root=self.cfg['WRF']['wrf_root']
        else:
            self.cplexe_root=self.machine_dic[f'{self.mode}_root']
        self.proj_root=os.path.join(
            self.cplexe_root,'Projects',self.nml_temp)
        
        self.sim_strt_time=datetime.datetime.strptime(
            cfg['model_init_ts'],'%Y%m%d%H')
        self.run_days=int(cfg['model_run_days'])
        self.sim_end_time=self.sim_strt_time+datetime.timedelta(days=self.run_days)
        
        self.arch_flag=cfg.getboolean('archive_flag')
        self.arch_root=utils.parse_fmt_timepath(self.sim_strt_time, cfg['arch_root'])
        utils.write_log('Uranus Initiation Done.')

    def waterfall(self):
        self.makewrf()
        self.makeroms()
        self.cplrock()
    
    def makewrf(self):
        from .lib import wrf_rocker
        if self.rock_wrf:
            self.wrfmaker=wrf_rocker.WRFRocker(self)
            self.wrfmaker.rock() 
    
    def makeroms(self):
        from.lib import roms_rocker
        if self.rock_roms:
            self.romsmaker=roms_rocker.ROMSRocker(self)
            self.romsmaker.build_icbc() 
    
    def cplrock(self):
        # rock the coupled model!
        if self.rock_flg:
            cfg,mode=self.cfg, self.mode
            self.ntasks_atm=int(cfg[mode]['ntasks_atm'])
            self.ntasks_iocn=int(cfg[mode]['ntasks_iocn'])
            self.ntasks_jocn=int(cfg[mode]['ntasks_jocn'])
            self.ntasks_ocn=self.ntasks_iocn*self.ntasks_jocn
            self.ntasks_all=self.ntasks_atm+self.ntasks_ocn 
            # copy files
            cfgfn=os.path.join(
                self.cfgdb_root, self.nml_temp, '*in')
            io.copy_files(cfgfn, self.proj_root)
            domfn=os.path.join(
                self.domdb_root, self.nml_temp, 'scrip*')
            io.symlink_files(domfn, self.proj_root)
            
            wrf_nml=os.path.join(
                cfg['WRF']['wrf_root'],'namelist.input')
            io.copy_files(wrf_nml, self.cplexe_root)
            
            wrf_drvfn=os.path.join(
                cfg['WRF']['wrf_root'],'wrf[l,f,i,b]*')
            io.symlink_files(wrf_drvfn, self.cplexe_root) 
 
            # sed cfg files
            cpl_in=os.path.join(self.proj_root, 'coupling.in')
            ocn_name=os.path.join('Projects', self.nml_temp, 'roms_d01.in')
            scrip_name=os.path.join('Projects', self.nml_temp, 'scrip.nc')
            sed_dic={
                'NnodesATM':self.ntasks_atm, 'NnodesOCN':self.ntasks_ocn,
                'OCN_name':ocn_name, 'SCRIP_COAWST_NAME':scrip_name
            }
            for key, itm in sed_dic.items():
                utils.sedline(key, f'{key} = {itm}', cpl_in)
            
            # roms flow
            self.romsmaker.prepare_rock()
        
            # run coawstM
            cmd=f'source {self.bashrc}; cd {self.cplexe_root};'
            cmd+=f'{self.mpicmd} -np {self.ntasks_all} ./coawstM {cpl_in}'
            cmd+=f' >& {self.cplexe_root}/coawstM.log'
            utils.write_log(print_prefix+'Run coawstM: '+cmd)
            subprocess.run(cmd, shell=True)
        
        if self.arch_flag:
            self.archive_data()
    
    def archive_data(self):
        if not(os.path.exists(self.arch_root)):
            utils.write_log(print_prefix+'mkdir '+self.arch_root)
            os.mkdir(self.arch_root)
            file_patterns=['wrf[o,r,x]*','*nc']
            for itm in file_patterns:
                io.move_files(os.path.join(self.cplexe_root,itm), self.arch_root)
    def _setup_logging(self):
        """
        Configures the logging module using the 
        'config.logging.ini' file in the installed package.
        """
        resource_path = os.path.join('conf','config.logging.ini')
        config_file=utils.fetch_pkgdata(resource_path)
        logging.config.fileConfig(
            config_file, disable_existing_loggers=False)


def copy_cfg(dest_path):
    """
    Copies a configuration file from the installed package to the destination path.
        dest_path (str): The path of the destination configuration file.
    """
    resource_path = os.path.join('conf','config.case.ini')
    src_path=utils.fetch_pkgdata(resource_path)
    shutil.copy2(src_path, dest_path)