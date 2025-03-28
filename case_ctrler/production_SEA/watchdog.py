'''
Watchdog for SEA project production run
Preprocessing, runtime, and postprocessing marches together.
the watchdog monitors the wrfout to decide driven data replacement.

    1. preprocess wrf[i,l,b,f]* and move to $runtime_temp
    2. link to $runtime_wrfdir
    3. run wrf.exe (single long run)
    4. monitor wrfout in $runtime_wrfdir 
    5. if wrfout reaches 1-day buffer, 
       5.1 redo step 1 and 2 with 1-day buffer
       5.2 postprocess wrfout and wrfrst
'''

import uranus
from uranus.lib import io, utils
import configparser
from datetime import datetime
import time
import os
import pandas as pd
import time
from threading import Thread



    
def runwrf():
    agent=uranus.Uranus(cfgfn='runtime.ini')
    agent.waterfall()


# config opts
precfg_file='preprocessor.ini'
runcfg_file='runtime.ini'
prenode='hqlx204'
init_ts='1969010100'
end_ts= '1969010200'
epoch_span='6H'
restart_run=False



src_files='/home/lzhenn/array204/WRF-4.1.5/run/wrf[i,l,b,f]*'
runtime_root='/home/lzhenn/hqnfs'
runtime_temp=f'{runtime_root}/temp'
runtime_wrfdir=f'{runtime_root}/WRF-4.1.5_org_P2/run/'
postdir=f'{runtime_root}/post'
dest_path=os.path.join(runtime_temp, 'prepare/')
runtime_path=os.path.join(runtime_temp, 'online')


span_hours = int(epoch_span[:-1])
restart_interval=span_hours*60
preprocess_span=f'{span_hours}H' # add 3-hr for wrffdda file
#preprocess_span=f'{span_hours+3}H' # add 3-hr for wrffdda file
# modify precfg
precfg=configparser.ConfigParser()
precfg.read(precfg_file)
precfg['URANUS']['machine_name']=prenode
precfg['URANUS']['model_init_ts']=init_ts
precfg['URANUS']['model_run_span']=preprocess_span
pre_agent = uranus.Uranus(cfg=precfg)

# first-time modify runcfg
runtime_cfg=configparser.ConfigParser()
runtime_cfg.read(runcfg_file)
runtime_cfg['URANUS']['model_init_ts']=init_ts
runtime_cfg['URANUS']['model_run_span']=epoch_span
runtime_cfg['WRF']['wrf_root']=runtime_wrfdir
# for restart
if restart_run:
    runtime_cfg['WRF']['nml_modification']=f'restart:.true.|restart_interval:{restart_interval:d}'
else:
    runtime_cfg['WRF']['nml_modification']=f'restart:.false.|restart_interval:{restart_interval:d}'

with open(runcfg_file, 'w') as configfile:
    runtime_cfg.write(configfile)



init_ts_obj=datetime.strptime(init_ts, '%Y%m%d%H')
end_ts_obj=datetime.strptime(end_ts,'%Y%m%d%H')
date_range = pd.date_range(
    start=init_ts_obj,
    end=end_ts_obj,
    freq=epoch_span)


nrst=len(date_range)-1    

# the first prepare
pre_agent.waterfall()

io.copy_files(src_files, dest_path)
io.del_files(runtime_path, 'wrf[i,l,b,f]*')
io.move_files(
    os.path.join(dest_path,'wrf[i,l,b,f]*'), runtime_path)  

print(date_range)

for i,idate in enumerate(date_range[:-1]):
    
    utils.write_log(
        f'{i+1:04d}/{nrst:04d} Loop: from {idate} to {date_range[i+1]}')
     
    # runtime
    wrf_thread = Thread(target=runwrf) 
    start_wrf_time=time.time()
    wrf_thread.start()

    time.sleep(30)  # Wait for 30 seconds 
    curr_ts=date_range[i+1].strftime('%Y%m%d%H') 
    
    # modify precfg and execute the next preparation
    utils.write_log(
        f'{i+1:04d}/{nrst:04d} Loop: prepare icbc for {curr_ts} run...')
    precfg['URANUS']['model_init_ts']=curr_ts
    pre_agent = uranus.Uranus(cfg=precfg)
    pre_agent.waterfall()
    io.copy_files(src_files, dest_path)
    
    
    utils.write_log(
        'wait for wrf.exe to complete before modifying runcfg...') 
    
    wrf_thread.join()
    end_wrf_time=time.time()
    delta_time=end_wrf_time - start_wrf_time
    
    io.copy_files(
        os.path.join(runtime_wrfdir,'rsl.out.0000'),
        os.path.join(runtime_wrfdir,f'rsl.out.0000.{i+1:04d}'))

    if  delta_time<300:
        utils.throw_error(f'wrf.exe run for {delta_time:d}s, not right.')
    utils.write_log('wrf.exe completed!') 
    
    io.del_files(runtime_path, 'wrf[i,l,b,f]*')
    io.move_files(
        os.path.join(dest_path,'wrf[i,l,b,f]*'), runtime_path)  
    
    utils.write_log('modifying runcfg...')
    runtime_cfg=configparser.ConfigParser()
    runtime_cfg.read(runcfg_file)
    runtime_cfg['URANUS']['model_init_ts']=curr_ts
    runtime_cfg['WRF']['nml_modification']=f'restart:.true.|restart_interval:{restart_interval:d}'
    with open(runcfg_file, 'w') as configfile:
        runtime_cfg.write(configfile)

