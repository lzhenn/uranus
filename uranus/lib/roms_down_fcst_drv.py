#!/usr/bin/python
# -*- coding: UTF-8 -*-
'''
Fetch fcst drv data  
    -- hycom
    -- cfs

Usage:
   python3 roms_down_fcst_drv.py $drv_type $yyyymmddhh $arch_path $ndays [$buffer]
e.g.
   # download hycom data initiate from 20240430 12Z to ./hycom_drv folder, 3-day data with 0 day buffer
   python3 roms_down_fcst_drv.py hycom 2024043012 ./hycom_drv 3 0

URL Flavors:
    CFS
        https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-9-month-forecast/6-hourly-ocean/2024/202404/20240426
    HYCOM
        https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/runs/GLBy0.08_930_FMRC_RUN_
        2020-12-29T12:00:00Z?var=salinity
        &north=23.0000&west=110.0000&east=120.9200&south=10.0000
        &disableProjSubset=on&horizStride=1
        &time=2021-01-06T00%3A00%3A00Z&vertCoord=&addLatLon=true&accept=netcdf4

                L_Zealot
                Apr 30, 2024

'''

import os, sys, time
import requests
import numpy as np
import datetime


# HYCOM
HYCOM_BASE0='https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/runs/GLBy0.08_930_FMRC_RUN_'
# HYCOM_VAR_RANGE='var=surf_el&north=30.0000&west=100.0000&east=130&south=10.0000&horizStride=1'
HYCOM_VAR_RANGE='var=surf_el&var=salinity&var=water_temp&var=water_u&var=water_v&north=30.0000&west=100.0000&east=130&south=10.0000&horizStride=1'
HYCOM_TAIL='T12%3A00%3A00Z&vertCoord=&addLatLon=true&accept=netcdf4'
#HYCOM_FSIZE0=33
HYCOM_FSIZE0=33000000
#----------------------------------------------------
# User Defined Part
#----------------------------------------------------
def fetch(cfg_dic={}):
    try_time=10
    sleep=120

    try:
        # drv type, hycom or cfs
        drv_type=cfg_dic['drv_type']
        # init time
        g_init_time=cfg_dic['g_init_time']
        # output dir
        fout_dir=cfg_dic['drv_dir']
        # simulation days
        sim_ndays=cfg_dic['sim_ndays']
        # buffer days in lead (in case of no data for the current day init)
        buf_days=cfg_dic['buf_days']
    except KeyError:
        # arguments in
        args=sys.argv
        drv_type=args[1]
        g_init_time=args[2]
        fout_dir=args[3]
        sim_ndays=int(args[4])
        buf_days=int(args[5])
    
    # tframes
    tfrms=np.arange(buf_days,sim_ndays+buf_days+1)
    # parser
    int_time_obj = datetime.datetime.strptime(g_init_time, '%Y%m%d%H')
    int_time_obj = int_time_obj+datetime.timedelta(days=-buf_days)
    
    if drv_type=='hycom':
        print('>>>>ROMS: HYCOM with '+str(buf_days)+' day(s) buffer: fetch from'+int_time_obj.strftime('%Y-%m-%d_%HZ'))
        # URL base str for forecast data
        url_base=HYCOM_BASE0+int_time_obj.strftime('%Y-%m-%d')+'T'+int_time_obj.strftime('%H')+':00:00Z?'
        
        for ifrm in tfrms:
            curr_filetime=int_time_obj+datetime.timedelta(days=int(ifrm))
            fn='hycom_'+curr_filetime.strftime('%Y%m%d%H')+'.nc'
            
            # sanity check
            if os.path.exists(fout_dir+'/'+fn):
                print('>>>> %s exist, test size...' % fn)
                fsize = os.path.getsize(fout_dir+'/'+fn)
                if (fsize>HYCOM_FSIZE0):
                    print('>>>> %s sanity check passed, skip this file...' % fn)
                    continue

            print('>>>>ROMS: Download %s...' % (curr_filetime.strftime('%Y%m%d')), end='')
            url_time='&time='+curr_filetime.strftime('%Y-%m-%d')
            url=url_base+HYCOM_VAR_RANGE+url_time+HYCOM_TAIL
            print(url)
            
            while try_time>0: 
                try:
                    rqst=requests.get(url)
                except:
                    try_time=try_time-1
                    print(fn+' fetch failed, sleep %ds to try again...' % sleep)
                    time.sleep(sleep)

                
                if (rqst.status_code == 200):
                    f = open(fout_dir+'/'+fn, 'wb')
                    f.write(rqst.content)
                    f.close()
                    print(fn+' status:'+str(rqst.status_code))

                    fsize = os.path.getsize(fout_dir+'/'+fn)
                    if (fsize>HYCOM_FSIZE0):
                        print('>>>> %s sanity check passed, next file...' % fn)
                        break
                try_time=try_time-1
                print(fn+' fetch failed, sleep %ds to try again...' % sleep)
                time.sleep(sleep)
            
            if try_time ==0:
                print('Maximum try time reached. Exit.')
        print('>>>>ROMS: HYCOM DATA ARCHIVED SUCCESSFULLY!!!')
    elif drv_type=='cfs':
        pass 
if __name__ == "__main__":
    fetch()



