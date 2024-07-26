
import numpy as np
import datetime
PKG_NAME = 'uranus'


MODE_FLAG_DIC={
    # atm, ocn, wav
    'shu':[1, 0, 0],
    'neptune':[0, 1, 0],
    'calypso':[0, 0, 1],
    'aegir':[1, 1, 0],
    'poseidon':[0, 1, 1],
    'aeolus':[1, 0, 1],
    'njord':[1, 1, 1]
}   


# WRF
UNGRIB_CLEAN_LIST=['CFS', 'ERA', 'SST', 'PFILE', 'GFS', 'CPSV3','BCMM', 'GRIBFILE']
METGRID_CLEAN_LIST=['met_em']
WRF_CLEAN_LIST=['met_em','wrfinput_','wrflowinp_','wrffda_','wrfrst_']
MACHINE_DIC={
    'hqlx74':{
        'bashrc':'/home/lzhenn/.bashrc_intel20_amd',
        'mpicmd':'mpirun', 'corespernode':128, 'nodes':1,
        'metgrid_np':16, 'real_np':16, 'wrf_np':32,
        'njord_root':'/home/lzhenn/array74/Njord_Calypso/COAWST_Njord',
        'aegir_root':'/home/lzhenn/array74/Njord_Calypso/COAWST_Aegir',
        'poseidon_root':'/home/lzhenn/array74/Njord_Calypso/COAWST_Poseidon',
        'calypso_root':'/home/lzhenn/array74/Njord_Calypso/COAWST_Calypso',
        'cfgdb_root':'/home/lzhenn/array74/workspace/uranus/uranus/nml_db/',
        'domdb_root':'/home/lzhenn/array74/Njord_Calypso/domaindb/'
    },
    'hqlx130':{
        'bashrc':'/home/lzhenn/.bashrc_intel20_amd',
        'mpicmd':'mpirun', 'corespernode':128, 'nodes':1,
        'metgrid_np':16, 'real_np':16, 'wrf_np':32,
        'njord_root':'/home/lzhenn/array130/COAWST_Njord',
        'aegir_root':'/home/lzhenn/array130/COAWST_Aegir',
        'poseidon_root':'/home/lzhenn/array130/COAWST_Poseidon',
        'calypso_root':'/home/lzhenn/array130/COAWST_Calypso',
        'cfgdb_root':'/home/lzhenn/array74/workspace/uranus/uranus/nml_db/',
        'domdb_root':'/home/lzhenn/array74/Njord_Calypso/domaindb/'
    },
    'hqlx181':{ # operational njord
        'bashrc':'/home/lzhenn/.bashrc_intel20_amd',
        'mpicmd':'mpirun', 'corespernode':192, 'nodes':1,
        'metgrid_np':32, 'real_np':48, 'wrf_np':32,
        'njord_root':'/home/lzhenn/array181/op_njord/model/COAWST_Njord/',
        'poseidon_root':'/home/lzhenn/array181/op_njord/model/COAWST_Poseidon',
        'opexe_root':'/home/lzhenn/array181/op_njord/model/COAWST_Njord/',
        'cfgdb_root':'/home/lzhenn/array74/workspace/uranus/uranus/nml_db/',
        'domdb_root':'/home/lzhenn/array74/Njord_Calypso/domaindb/'
    },
    'hqlx86':{
        'bashrc':'/home/metctm1/.bashrc_i20_cwst38',
        'mpicmd':'mpirun', 'corespernode':48, 'nodes':1,
        'metgrid_np':16, 'real_np':16, 'wrf_np':32,
        'aegir_root':'/home/metctm1/array_hq86/COAWSTv38',
        'cfgdb_root':'/disk/r086/metctm1/work/uranus/uranus/nml_db/',
        'domdb_root':'/home/lzhenn/array74/workspace/aegir-implement/domaindb/'
    },
    'th2':{
        'bashrc':'/BIGDATA2/sysu_xmhu_1/s2s_cpl/bashrc',
        'mpicmd':'yhrun', 'corespernode':24, 'nodes':5,
        'metgrid_np':16, 'real_np':16, 'wrf_np':24,
        'aegir_root':'/BIGDATA2/sysu_xmhu_1/s2s_cpl/model/COAWST_Aegir',
        'cfgdb_root':'/BIGDATA2/sysu_xmhu_1/s2s_cpl/control/uranus/uranus/nml_db/',
        'domdb_root':'/BIGDATA2/sysu_xmhu_1/s2s_cpl/data/domaindb/'
    },    
    'pird':{
        'bashrc':'/g6/cmme/COAWST-S2S/bashrc_coawst_run',
        'chck_cmd':'squeue | grep cmme',
        'mpicmd':'sbatch', 'corespernode':32, 'nodes':3,
        'metgrid_np':16, 'real_np':16, 'wrf_np':64,
        'aegir_root':'/g6/cmme/COAWST-S2S/model/COAWST_Aegir',
        'opexe_root':'/g6/cmme/COAWST-S2S/model/op_model/COAWST_op/',
        'cfgdb_root':'/g6/cmme/COAWST-S2S/control/uranus/uranus/nml_db/',
        'domdb_root':'/g6/cmme/COAWST-S2S/data/domaindb/'
    }

}
for itm in [
    '47','54','62','65','69','129',
    '133','132','111','100','174','179','182']:
    MACHINE_DIC['hqlx'+itm]=MACHINE_DIC['hqlx74']

NML_SPECIFIC={
    'SEAsia_4km':{'metgrid_np':48, 'real_np':64}
}
       
DRV_DIC={
    # WRF
    'cpsv3_wrf':{
        'use_ungrib':False,'atm_nfrq':'6H','atm_file_nfrq':'1D', 'frm_per_file':4,
        'vcoord':'sigmap', 'acoef':'hyam', 'bcoef':'hybm', 'psname':'PS', 
        'lnd_file_nfrq':'1D', 'soillv':'cpsv3', 'soil_dim_name':'levsoi',
        'atm_naming_pattern':'$I%Y%m%d%H$I.cam2.h4.$F%Y-%m-%d$F-00000.nc',
        'lnd_naming_pattern':'$I%Y%m%d%H$I.clm2.h2.$F%Y-%m-%d$F-00000.nc',
        'lats':np.linspace(-89.6559642468699, 89.6559642468699, 400), 
        'lons':np.linspace(0.0,359.55,800), 'plv':'PL18',
    },
    'bcmm_wrf':{
        'use_ungrib':False,'atm_nfrq':'6H','atm_file_nfrq':'MS', 'frm_per_file':-1,
        'vcoord':'p', 'fg_name':"'BCMM'", 
        'lnd_file_nfrq':'MS', 'soillv':'bcmm', 'soil_dim_name':'depth',
        'atm_naming_pattern':'atm_$S_$F%Y_%m$F.nc4',
        'lnd_naming_pattern':'lnd.$S.$F%Y%m$F.nc',
        'lats':np.linspace(-90, 90, 145), 
        'lons':np.linspace(0.0,358.75,288), 'plv':'PL18',
    },
    'MPI-ESM1-2-HR_wrf':{
        'use_ungrib':False,'atm_nfrq':'6H','atm_file_nfrq':'MS', 'frm_per_file':-1,
        'vcoord':'p',  
        'lnd_file_nfrq':'MS', 'soillv':'bcmm', 'soil_dim_name':'depth',
        'atm_naming_pattern':'atm_$S_$F%Y_%m$F.nc4',
        'lnd_naming_pattern':'lnd.$S.$F%Y%m$F.nc',
        'lats':np.linspace(-90, 90, 145), 
        'lons':np.linspace(0.0,358.75,288), 'plv':'PL18',
    },
    'era5_wrf':{
        'use_ungrib':True,'atm_nfrq':'3H','nsoil':4,'plv':'PL26',
        'Vtable':'Vtable.ERA-interim.pl','ungrib_prefixes':['ERASL','ERAPL'],
        'fg_name':"'ERASL','ERAPL'",
        'file_patterns':['%Y%m*-sl.grib','%Y%m*-pl.grib'],
    },
    'gfs_wrf':{
        'use_ungrib':True,'atm_nfrq':'3H','nsoil':4,'plv':'PL33',
        'Vtable':'Vtable.GFS','ungrib_prefixes':['GFS'],
        'fg_name':"'GFS'",
        'file_patterns':['gfs.t??z.pgrb2.0p25.f*'],
    },
    # ROMS
    'cpsv3_roms':{
        'ocn_nfrq':'24H','ocn_file_nfrq':'1D',
        'ocn_naming_pattern':'ocean_$F%Y_%m_%d$F.nc',
    },
    'cfs_roms':{'ocn_naming_pattern':'ocnf$F%Y%m%d%H$F.01.$I%Y%m%d%H$I.grb2'},
    'cfsr_roms':{
        'ocn_nfrq':'24H','ocn_file_nfrq':'1D',
        'ocn_naming_pattern':'ocnh01.gdas.$I%Y%m%d%H$I.grb2'},
    'hycom_roms':{
        'ocn_nfrq':'24H','ocn_file_nfrq':'1D',
        'ocn_naming_pattern':'hycom_$F%Y%m%d%H$F.nc'},
    'roms_roms':{'ocn_naming_pattern':'roms_his_domid_%05d.nc4'},
    # SWAN
    'era5_swan':{'wav_naming_pattern':'$F%Y%m%d$F-wv.grib', 'wav_nfrq':'3H', 'wav_file_nfrq':'1D','frm_per_file':8}, 
    'era5mon_swan':{'wav_naming_pattern':'$F%Y%m$F-sl.grib', 'wav_nfrq':'3H', 'wav_file_nfrq':'1M','frm_per_file':240}, 
    'gfs_swan':{'wav_naming_pattern':'gfswave.t12z.global.0p16.f$A.grib2','wav_nfrq':'3H','wav_file_nfrq':'3H'}    
}

# Lower bottom level 
SOILLV_DIC={
    'cpsv3':0.01*np.array([1.75, 4.51, 9.06, 16.55, 28.91, 49.29, 82.89, 138.28, 229.61, 380.19]),
    'bcmm':np.array([0.0, 0.1, 0.4, 1.0, 2.0])
}

PLV_DIC={
    'PL14':100.0*np.array([1000.0, 925.0, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50]),
    'PL18':100.0*np.array([1000.0, 975.0, 950.0, 925.0, 900.0, 850, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50])
}


# -----------------------ROMS--------------------------

# BASE TIME FOR ROMS
BASE_TIME=datetime.datetime(1858, 11, 17, 0, 0)
S2NS=1000000000
HALF_DAY_NS=43200000000000 # half day in nanoseconds
DAY_IN_SEC=86400
HR_IN_SEC=3600
MIN_IN_SEC=60
# !!! FIXED DEPTH: NEVER CHANGE THIS BELOW !!!
DENSE_DP=np.concatenate((
    np.arange(0,50), np.arange(50, 200,10),
    np.arange(200, 1000, 100),np.arange(1000,6000,500))).astype(int)
# !!! FIXED DEPTH: NEVER CHANGE THIS ABOVE !!!

ROMS_VAR_NAME_MAPS={
    'temp':['temp','temperature','water_temp','pt'],
    'salt':['salt','salinity','s'],
    'zeta':['zeta','ssh','sea_surface_height','surf_el','eta_t'],
    'u':['u','water_u','uoe'],
    'v':['v','water_v','von'],
    'lat_rho':['Y','Latitude','lat','yt_ocean','eta_rho'],
    'lon_rho':['X','Longitude','lon','xt_ocean','xi_rho']}

ROMS_VAR=['zeta','temp','salt','u','v','ubar','vbar']
CLM_TIME_VAR=['v2d','v3d','temp','salt', 'zeta', 'ocean']
BDY_TIME_VAR=['v2d','v3d','temp','salt', 'zeta']

ROMS_WRF_FORC_MAPS={
    'Uwind': 'U10', 'Vwind': 'V10', 
    'Pair': 'slp', 'Tair': 'T2',  
    'Qair': 'rh2', 
    #'swrad': 'GSW', 
    'swrad': 'SWDOWN', 
    #'lwrad': 'LWDNB', 
    'lwrad': 'GLW', 
    'lwrad_down': 'GLW'}
FORC_TIMEVAR_LIST=['wind','pair','tair','qair','rain','srf','lrf']
FORC_VAR_DIC={'Uwind':{'long_name': 'surface u-wind component', 
                  'units': 'meter second-1', 
                  'field': 'Uwind, scalar, series', 
                  'coordinates': 'lon lat', 
                  'time': 'wind_time'},
        'Vwind':{'long_name': 'surface v-wind component',
                'units': 'meter second-1',
                'field': 'Vwind, scalar, series',
                'coordinates': 'lon lat',
                'time': 'wind_time'},
        'Pair':{'long_name': 'surface air pressure',
                'units': 'millibar',
                'field': 'Pair, scalar, series',
                'coordinates': 'lon lat',
                'time': 'pair_time'},
        'Tair':{'long_name': 'surface air temperature',
                'units': 'celsius',
                'field': 'Tair, scalar, series',
                'coordinates': 'lon lat',
                'time': 'tair_time'},
        'Qair':{'long_name': 'surface relative humidity',
                'units': 'percentage',
                'field': 'Qair, scalar, series',
                'coordinates': 'lon lat',
                'time': 'qair_time'},
        'rain':{'long_name': 'rain fall rate',
                'units':'kilogram meter-2 second-1',
                'field': 'Rain, scalar, series',
                'coordinates': 'lon lat',
                'time': 'rain_time'},
        'swrad':{'long_name': 'net solar shortwave radiation', 
                'units': 'Watts meter-2', 
                'positive_value': 'downward flux, heating', 
                'negative_value': 'upward flux, cooling', 
                'field': 'swrad, scalar, series', 
                'coordinates': 'lon lat', 
                'time': 'srf_time'},
        'lwrad':{'long_name': 'net downward longwave radiation',
                'units': 'Watts meter-2', 
                'positive_value': 'downward flux, heating', 
                'negative_value': 'upward flux, cooling', 
                'field': 'lwrad, scalar, series', 
                'coordinates': 'lon lat', 
                'time': 'lrf_time'},
        'lwrad_down':{'long_name': 'downward solar longwave radiation', 
                'units': 'Watts meter-2', 
                'positive_value': 'downward flux, heating', 
                'negative_value': 'upward flux, cooling', 
                'field': 'lwrad_down, scalar, series', 
                'coordinates': 'lon lat', 'time': 'lrf_time'}
}
ROMS_CLEAN_LIST=['roms_bdy_','roms_ini_','roms_clm_']


# -----------------------SWAN--------------------------
SWAN_CLEAN_LIST=['swan_bdy','Errfile','PRINT','lwavp_',
                'dir_','hsig_','hswell_','tps_','per_','wind_']




# Physics
RHO_WATER=1000.0
K2C=273.15



# Dataset 
# -----------------------ERA5--------------------------
ERA5_CONST={
    'SL_VARS':[
                '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
                '2m_temperature','land_sea_mask','mean_sea_level_pressure',
                'sea_ice_cover','sea_surface_temperature','skin_temperature',
                'snow_depth','soil_temperature_level_1','soil_temperature_level_2',
                'soil_temperature_level_3','soil_temperature_level_4','surface_pressure',
                'volumetric_soil_water_layer_1','volumetric_soil_water_layer_2','volumetric_soil_water_layer_3',
                'volumetric_soil_water_layer_4'],
    'PL_LAYS':[
                '50','70','100','150','200','250',
                '300','350','400','450','500','550','600',
                '650','700','750','775','800','825',
                '850','875','900','925','950','975','1000'],
    'PL_VARS':[
                'geopotential','relative_humidity','specific_humidity',
                'temperature','u_component_of_wind','v_component_of_wind'],
}
