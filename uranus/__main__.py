#!/usr/bin/env python3
'''
Date: Sep 08, 2023
uranus is the main driver  

This is the main script to drive the model

History:
Sep 08, 2023 --- Kick off the project 

L_Zealot
'''
import sys, os
import logging, logging.config
import shutil
import datetime

from .lib import cfgparser, utils

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


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
        
        self.mode=cfg['uranus_mode']
        self.nml_temp=cfg['nml_temp']
        self.machine_name=cfg['machine_name']
        self.cfgdb_root=cfg['cfgdb_root']
        self.domdb_root=cfg['domdb_root']
        
        self.sim_strt_time=datetime.datetime.strptime(
            cfg['model_init_ts'],'%Y%m%d%H')
        self.run_days=int(cfg['model_run_days'])
        self.sim_end_time=self.sim_strt_time+datetime.timedelta(days=self.run_days)
        
        utils.write_log('Uranus Initiation Done.')

    def waterfall(self):
        self.makewrf()
        
    def makewrf(self):
        from .lib import wrf_rocker
        self.wrfmaker=wrf_rocker.WRFRocker(self)
        self.wrfmaker.make_icbc() 
    
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