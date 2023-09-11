
import numpy as np

PKG_NAME = 'uranus'

# WRF
WPS_CLEAN_LIST=['met_em', 'CFS', 'ERA', 'SST', 'PFILE', 'GFS']
WRF_CLEAN_LIST=['met_em','wrfinput_','wrflowinp_','wrffda_','wrfrst_']

DRV_DIC={
    'cpsv3':{
        'use_ungrib':False,'atm_nfrq':'6H','atm_file_nfrq':'1D', 'frm_per_file':4,
        'vcoord':'sigmap', 'acoef':'hyam', 'bcoef':'hybm', 'psname':'PS', 
        'soillv':'CLM2', 'soil_dim_name':'levsoi',
        'atm_naming_pattern':'$I%Y%m%d%H$I.cam2.h4.$F%Y-%m-%d$F-00000.nc',
        'lnd_naming_pattern':'$I%Y%m%d%H$I.clm2.h2.$F%Y-%m-%d$F-00000.nc',
        'lats':np.linspace(-89.6559642468699, 89.6559642468699, 400), 
        'lons':np.linspace(0.0,359.55,800), 'plv':'PL18',
    }
}

# Lower bottom level 
SOILLV_DIC={
    'CLM2':0.01*np.array([0, 1.75, 4.51, 9.06, 16.55, 28.91, 49.29, 82.89, 138.28, 229.61, 380.19]),
}

PLV_DIC={
    'PL14':100.0*np.array([1000.0, 925.0, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50]),
    'PL18':100.0*np.array([1000.0, 975.0, 950.0, 925.0, 900.0, 850, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50])
}

# Physics
RHO_WATER=1000.0