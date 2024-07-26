#/usr/bin/env python3
"""Driven Data Crawler"""
from . import const, utils, io, roms_down_fcst_drv
import datetime, os, requests
print_prefix='lib.Crawler>>'
class Crawler():
    def __init__(self, uranus):
        utils.write_log(print_prefix+'Crawler Initiation.')
        self.cfg=uranus.cfg
        self.start_time=uranus.sim_strt_time
        self.end_time=uranus.sim_end_time
        self.uranus=uranus
    def down_ocn(self):
        utils.write_log(print_prefix+'Download driven data for ocean...')
        romsmaker=self.uranus.romsmaker
        curr_time_obj = self.start_time 
        drv_root=romsmaker.drv_root
        drv_type=romsmaker.drv_type
        io.check_mkdir(drv_root) 
        if self.uranus.fcst_flag:
            buf_day,sim_ndays=0,7
            g_init_time=curr_time_obj-datetime.timedelta(days=buf_day)
            cfg_dic={
                'drv_type':drv_type,'g_init_time':g_init_time.strftime('%Y%m%d%H'),
                'drv_dir':drv_root,'sim_ndays':sim_ndays,'buf_days':buf_day}
            roms_down_fcst_drv.fetch(cfg_dic)
        else:
            # analysis data
            if drv_type=='cfsr':
                url_base='https://www.ncei.noaa.gov/oa/prod-cfs-reanalysis/6-hourly-ocean/'
                #https://www.ncei.noaa.gov/oa/prod-cfs-reanalysis/6-hourly-ocean/1979/197901/19790101/ocnh01.gdas.1979010100.grb2
                if curr_time_obj<=self.end_time: 
                    utils.write_log(
                        f'{print_prefix}download CFSR Ocean @ {curr_time_obj.strftime("%Y%m%d%H")} to {drv_root}...')  
                    fn='ocnh01.gdas.'+curr_time_obj.strftime('%Y%m%d%H')+'.grb2'
                    # sanity check
                    if os.path.exists(drv_root+'/'+fn):
                        utils.write_log(f'{print_prefix} exists.')
                    yyyy,yyyymm=curr_time_obj.strftime('%Y'),curr_time_obj.strftime('%Y%m')
                    yyyymmdd=curr_time_obj.strftime('%Y%m%d')
                    url=url_base+yyyy+'/'+yyyymm+'/'+yyyymmdd+'/'+fn
                    try:
                        rqst=requests.get(url)
                    except:
                        write_log(f'{print_prefix}CFSR Ocean fetch failed, exit...')
                        exit()
                    if (rqst.status_code == 200):
                        f = open(drv_root+'/'+fn, 'wb')
                        f.write(rqst.content)
                        f.close()
                        utils.write_log(f'{print_prefix} {fn} done.')
                    curr_time_obj=curr_time_obj+datetime.timedelta(days=1)
            elif drv_type=='hycom':
                # parser
                print('>>>>ROMS: HYCOM fetch from '+int_time_obj.strftime('%Y-%m-%d_%HZ'))
                df_exp_info   =  pd.read_csv(CWD+'/hycom_list.txt', sep='\s+')

                # domain group binding
                df_dom_info=pd.read_csv(CWD+'/../../db/hycom_db.csv',index_col='dom_grp')
                range=df_dom_info.loc[dom_group]
                range_seg='&north='+str(range['north'])
                range_seg=range_seg+'&west='+str(range['west'])
                range_seg=range_seg+'&east='+str(range['east'])
                range_seg=range_seg+'&south='+str(range['south'])

                # find exp binding
                for idx, itm in df_exp_info.iterrows():
                    bind_time_strt=datetime.datetime.strptime(itm['date_strt'],'%Y-%m-%d')
                    if int_time_obj>= bind_time_strt:
                        exp_series=itm['exp_series']
                        exp_name=itm['exp_name']
                        break
                    
                    for ifrm in tfrms:
                        curr_filetime=int_time_obj+datetime.timedelta(days=int(ifrm))
                        # URL base str for forecast data
                        url_base='https://ncss.hycom.org/thredds/ncss/'+exp_series+'0.08/'+exp_name
                        if curr_filetime>datetime.datetime(2016,4,18):
                            url_base=url_base+'?'
                            url_var_range='var=surf_el&var=salinity&var=water_temp&var=water_u&var=water_v'
                        else:
                            url_base=url_base+'/'+curr_filetime.strftime('%Y')+'?'
                            url_var_range='var=ssh&var=salinity&var=temperature&var=u&var=v'

                        url_var_range=url_var_range+range_seg+'&horizStride=1'
                        url_time='&time='+curr_filetime.strftime('%Y-%m-%dT%H')+'%3A00%3A00Z'
                        #url_time_start='&time_start='+curr_filetime.strftime('%Y-%m-%dT%H')+'%3A00%3A00Z'
                        #url_time_end='&time_end='+curr_filetime.strftime('%Y-%m-%dT%H')+'%3A00%3A00Z'
                        #url_tail='&timeStride=1&vertCoord=&addLatLon=true&accept=netcdf4'
                        url_tail='&vertCoord=&addLatLon=true&accept=netcdf4'
                        url=url_base+url_var_range+url_time+url_tail
                        
                        fn='hycom_%s.nc' % curr_filetime.strftime('%Y%m%d%H')
                        while try_time>0: 
                            try:
                                print('>>>>ROMS: Download %s----%s-->%s' % (
                                    curr_filetime.strftime('%Y%m%d'), url, fout_dir))
                                rqst=requests.get(url)
                            except:
                                try_time=try_time-1
                                print(fn+' fetch failed, sleep %ds to try again...' % sleep)
                                time.sleep(sleep)
                            
                            if rqst.status_code == 200:
                                break
                            else:
                                try_time=try_time-1
                                print(fn+' fetch failed, sleep %ds to try again...' % sleep)
                                time.sleep(sleep)
                        if try_time==0:
                            print('Download failed!')
                            exit(1)
                        f = open(fout_dir+'/'+fn, 'wb')
                        f.write(rqst.content)
                        f.close()
                        print(fn+' status:'+str(rqst.status_code))
                    

    def down_atm(self):
        
        curr_time_obj = self.start_time 
        drv_root=self.cfg['WRF']['drv_root']
        io.check_mkdir(drv_root) 
        
        ''' Download ERA5 pres and single lv data '''
        if self.cfg['WRF']['drv_type']=='era5':
            import cdsapi
            c = cdsapi.Client()
            file_time_delta=datetime.timedelta(days=1)
            domain=self.cfg['era5_atm']['area_nwse']
            frq=self.cfg['era5_atm']['frq']
            while curr_time_obj <= self.end_time:
                date_str=curr_time_obj.strftime('%Y%m%d')
                utils.write_log(
                    f'{print_prefix}download ERA5 single layer: {date_str}') 
                
                # single layer data retriever
                sl_fn=os.path.join(drv_root,f'{date_str}-sl.grib')
                if not os.path.exists(sl_fn):
                    c.retrieve(
                        'reanalysis-era5-single-levels',
                        {
                            'product_type':'reanalysis',
                            'format':'grib',
                            'variable':const.ERA5_CONST['SL_VARS'],
                            'date':date_str,
                            'area':domain,
                            'time':f'00/to/23/by/{frq}',
                        }, sl_fn)
                # pressure layer data retriever
                utils.write_log(
                    f'{print_prefix}download multiple layers: {date_str}') 
                pl_fn=os.path.join(drv_root,f'{date_str}-pl.grib')
                if not os.path.exists(pl_fn):
                    c.retrieve(
                        'reanalysis-era5-pressure-levels',
                        {
                            'product_type':'reanalysis',
                            'format':'grib',
                            'pressure_level':const.ERA5_CONST['PL_LAYS'],
                            'date':date_str,
                            'area':domain,
                            'time':f'00/to/23/by/{frq}',
                            'variable':const.ERA5_CONST['PL_VARS'],
                            }, pl_fn)
                    
                curr_time_obj=curr_time_obj+file_time_delta
            