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
import time
from .lib import cfgparser, utils, const, io, crawler

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

    def __init__(self, cfgfn='config.case.ini',cfg=None):
        self._setup_logging()

        if cfg is None:
            if not(os.path.exists(os.path.join(CWD,cfgfn))):
                utils.write_log('config file not exist, copy from pkg...')
                copy_cfg(os.path.join(CWD,cfgfn))
            self.cfg=cfgparser.read_cfg(os.path.join(CWD,cfgfn))
        elif cfg=='rewrite':
            utils.write_log('forced rewrite config file, copy from pkg...')
            copy_cfg(os.path.join(CWD,cfgfn))
            self.cfg=cfgparser.read_cfg(os.path.join(CWD,cfgfn))
        else:
            utils.write_log('use passed config...')
            self.cfg=cfg
        cfg=self.cfg['URANUS']
        
        self.rock_flg=cfg.getboolean('rock_cpl')
        self.rock_wrf=self.cfg['WRF'].getboolean('rock_wrf')
        self.rock_roms=self.cfg['ROMS'].getboolean('rock_roms')
        self.rock_swan=self.cfg['SWAN'].getboolean('rock_swan')
        self.mode=cfg['uranus_mode']
        # atm, ocn, wav
        self.active_comp=const.MODE_FLAG_DIC[self.mode]
        self.nml_temp=cfg['nml_temp']
        self.domain_lv=cfg['domain_lv']
        self.sim_strt_time=utils.parse_init_time(cfg['model_init_ts'])
        if cfg['model_init_ts'].startswith('T'):
            self.fcst_flag=True
        else:
            self.fcst_flag=False
        self.run_hours=utils.parse_runspan(cfg['model_run_span'])
        self.sim_end_time=self.sim_strt_time+datetime.timedelta(hours=self.run_hours)


        self.machine_name=cfg['machine_name']
        self.machine_dic=const.MACHINE_DIC[self.machine_name]
        self.update_machine_dic()
        self.bashrc=self.machine_dic['bashrc']
        self.mpicmd=self.machine_dic['mpicmd']
        self.cfgdb_root=self.machine_dic['cfgdb_root']
        self.domdb_root=self.machine_dic['domdb_root']
        
            
        self.ntasks_atm, self.ntasks_wav,\
        self.ntasks_ocn, self.ntasks_iocn,self.ntasks_jocn,\
        self.ntasks_all=utils.get_ntasks(self.cfg, self.mode) 
        
        if self.mode=='shu':
            self.cplexe_root=self.cfg['WRF']['wrf_root']
        else:
            self.cplexe_root=self.machine_dic[f'{self.mode}_root']
        
        self.proj_root=os.path.join(
            self.cplexe_root,'Projects',self.nml_temp)
        io.check_mkdir(self.proj_root)
        
        if (
            self.cfg['WRF'].getboolean('down_drv_data') or\
            self.cfg['ROMS'].getboolean('down_drv_data') or\
            self.cfg['SWAN'].getboolean('down_drv_data')):
            self.crawler=crawler.Crawler(self)
        
        self.arch_flag=cfg.getboolean('archive_flag')
        self.arch_root=utils.parse_fmt_timepath(self.sim_strt_time, cfg['arch_root'])
        # copy config file 
        cfgfn=os.path.join(
        self.cfgdb_root, self.nml_temp, '*in')
        io.copy_files(cfgfn, self.proj_root)
        utils.write_log('Uranus Initiation Done.')

    def spinup_prep(self):
        pass
    def waterfall(self):
        '''
        waterfall the main stream 
        '''
        self.makewrf()
        self.makeroms()
        self.makeswan()
        self.cplrock()
 
    def update_machine_dic(self):
        try:
            spec=const.NML_SPECIFIC[self.nml_temp]
            for k, v in spec.items():
                self.machine_dic[k]=v
        except KeyError:
            pass
        try:
            self.machine_dic[f'{self.mode}_root']=self.cfg['cwst_path']
        except KeyError:
            pass
 
   
    def makewrf(self):
        from .lib import wrf_rocker
        if self.rock_wrf:
            self.wrfmaker=wrf_rocker.WRFRocker(self)
            self.wrfmaker.rock() 
    
    def makeroms(self):
        from.lib import roms_rocker
        if self.rock_roms:
            self.romsmaker=roms_rocker.ROMSRocker(self)
            self.romsmaker.rock() 
    def makeswan(self):
        from .lib import swan_rocker
        if self.rock_swan:
            self.swanmaker=swan_rocker.SWANRocker(self)
            self.swanmaker.rock()
    def cplrock(self):
        # rock the coupled model!
        if self.rock_flg:
            cfg=self.cfg
            domlv=self.domain_lv
            # copy files
          
            domfn=os.path.join(
                self.domdb_root, self.nml_temp, 'scrip*')
            io.symlink_files(domfn, self.proj_root)
            if self.ntasks_atm>0: 
                wrf_nml=os.path.join(
                    cfg['WRF']['wrf_root'],'namelist.input')
                io.copy_files(wrf_nml, self.cplexe_root)
                
                wrf_drvfn=os.path.join(
                    cfg['WRF']['wrf_root'],'wrf[l,f,i,b]*')
                io.symlink_files(wrf_drvfn, self.cplexe_root) 
    
            # sed cfg files
            cpl_in=os.path.join(self.proj_root, 'coupling.in')
            ocn_name=os.path.join('Projects', self.nml_temp, f'roms_{domlv}.in')
            wav_name=os.path.join('Projects', self.nml_temp, f'swan_{domlv}.in')
            scrip_name=os.path.join('Projects', self.nml_temp, f'scrip.{domlv}.nc')
            sed_dic={
                'NnodesATM':self.ntasks_atm, 'NnodesOCN':self.ntasks_ocn, 'NnodesWAV':self.ntasks_wav, 
                'WAV_name':wav_name,'OCN_name':ocn_name, 
                'SCRIP_COAWST_NAME':scrip_name
            }
            for key, itm in sed_dic.items():
                utils.sedline(key, f'{key} = {itm}', cpl_in, count=1)
            
            # roms flow
            self.romsmaker.prepare_cplrock()
            # swan flow
            self.swanmaker.prepare_cplrock()
            # run coawstM
            # special for cmme
            if self.machine_name == 'pird':
                cmd=f'mpirun -n {self.ntasks_all} ./coawstM'
            else:
                cmd=utils.build_execmd(
                    self.machine_name, self.bashrc, self.cplexe_root, 
                    self.mpicmd, self.ntasks_all, 'coawstM')
 
            relative_cpl_in=os.path.join(
                '.','Projects',self.nml_temp,'coupling.in')
            cmd+=f' {relative_cpl_in}'
            cmd+=f' >& {self.cplexe_root}/coawstM.{datetime.datetime.now().strftime("%y%m%d%H%M%S")}.log' 
            
            if self.machine_name == 'pird':
                cwst_sbatch=f'{self.cplexe_root}/coawstM.sh'
                utils.sedline('#SBATCH -n',f'#SBATCH -n {self.ntasks_all}', cwst_sbatch) 
                utils.sedline('time mpirun', cmd, cwst_sbatch) 
                batchcmd=utils.build_execmd(
                self.machine_name, self.bashrc, self.cplexe_root, 
                    self.mpicmd, self.ntasks_all, 'coawstM')
                utils.write_log(print_prefix+'Run coawstM: '+batchcmd)
                rcode=subprocess.run(batchcmd, shell=True, stdout=subprocess.PIPE)
                jobid=rcode.stdout.decode().split()[3]
                # special for pird
                chck_cmd='squeue | grep cmme'
                rcode=subprocess.run(chck_cmd, shell=True, stdout=subprocess.PIPE)
                stdout=rcode.stdout.decode()
                timer=180
                while jobid in stdout:
                    utils.write_log(f'{print_prefix}({timer}s check) On: {stdout}')
                    time.sleep(timer)
                    rcode=subprocess.run(chck_cmd, shell=True, stdout=subprocess.PIPE)
                    stdout=rcode.stdout.decode()
            else:
                utils.write_log(print_prefix+'Run coawstM: '+cmd)
                rcode=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
                
        if self.arch_flag:
            self.archive_data()
    
    def archive_data(self):
        io.check_mkdir(self.arch_root)
        # mv files
        file_patterns=['wrf[o,r,x]*','roms_his*','swan_his*']
        for itm in file_patterns:
            cmd=f'mv {os.path.join(self.cplexe_root,itm)} {self.arch_root}'
            utils.write_log(print_prefix+'Archive: '+cmd)
            rcode=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        # cp files
        file_patterns=['namelist.input']
        for itm in file_patterns:
            wrf_nml=os.path.join(self.cplexe_root,itm)
            io.copy_files(wrf_nml, self.arch_root)
            
        if self.active_comp[1]==1:
            io.zip_roms_his(self.arch_root) 
        #for itm in file_patterns:
        #    io.move_files(os.path.join(self.cplexe_root,itm), self.arch_root)
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
