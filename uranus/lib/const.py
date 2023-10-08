
import numpy as np
import datetime
PKG_NAME = 'uranus'

# WRF
UNGRIB_CLEAN_LIST=['CFS', 'ERA', 'SST', 'PFILE', 'GFS', 'CPSV3']
METGRID_CLEAN_LIST=['met_em']
WRF_CLEAN_LIST=['met_em','wrfinput_','wrflowinp_','wrffda_','wrfrst_']
MACHINE_DIC={
    'hqlx74':{
        'bashrc':'/home/lzhenn/.bashrc_intel20_amd',
        'mpicmd':'mpirun', 'corespernode':128, 'nodes':1,
        'metgrid_np':16, 'real_np':16, 'wrf_np':32,
        'aegir_root':'/home/lzhenn/array74/Njord_Calypso/COAWST_Aegir',
        'cfgdb_root':'/home/lzhenn/array74/workspace/aegir-implement/db/',
        'domdb_root':'/home/lzhenn/array74/workspace/aegir-implement/domaindb/'
    }
}

DRV_DIC={
    'cpsv3_wrf':{
        'use_ungrib':False,'atm_nfrq':'6H','atm_file_nfrq':'1D', 'frm_per_file':4,
        'vcoord':'sigmap', 'acoef':'hyam', 'bcoef':'hybm', 'psname':'PS', 
        'lnd_file_nfrq':'1D', 'soillv':'cpsv3', 'soil_dim_name':'levsoi',
        'atm_naming_pattern':'$I%Y%m%d%H$I.cam2.h4.$F%Y-%m-%d$F-00000.nc',
        'lnd_naming_pattern':'$I%Y%m%d%H$I.clm2.h2.$F%Y-%m-%d$F-00000.nc',
        'lats':np.linspace(-89.6559642468699, 89.6559642468699, 400), 
        'lons':np.linspace(0.0,359.55,800), 'plv':'PL18',
    },
    'cpsv3_roms':{
        'ocn_nfrq':'24H','ocn_file_nfrq':'1D',
        'ocn_naming_pattern':'ocean_$F%Y_%m_%d$F.nc',
    },
    'cfs_roms':{'ocn_naming_pattern':'ocnf$F%Y%m%d%H$F.01.$I%Y%m%d%H$I.grb2'},
    'hycom_roms':{'ocn_naming_pattern':'hycome_$F%Y%m%d%H$F.nc',},
    
     'bcmm_wrf':{
        'use_ungrib':False,'atm_nfrq':'6H','atm_file_nfrq':'MS', 'frm_per_file':-1,
        'vcoord':'p',  
        'lnd_file_nfrq':'MS', 'soillv':'bcmm', 'soil_dim_name':'depth',
        'atm_naming_pattern':'atm_$S_$F%Y_%m$F.nc4',
        'lnd_naming_pattern':'lnd.$S.$F%Y%m$F.nc',
        'lats':np.linspace(-90, 90, 145), 
        'lons':np.linspace(0.0,358.75,288), 'plv':'PL18',
    },
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
    'lat_rho':['Y','Latitude','lat','yt_ocean'],
    'lon_rho':['X','Longitude','lon','xt_ocean']}

ROMS_VAR=['zeta','temp','salt','u','v','ubar','vbar']
CLM_TIME_VAR=['v2d','v3d','temp','salt', 'zeta', 'ocean']
BDY_TIME_VAR=['v2d','v3d','temp','salt', 'zeta']



# Physics
RHO_WATER=1000.0
K2C=273.15