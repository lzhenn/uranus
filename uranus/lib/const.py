
PKG_NAME = 'uranus'

# WRF
WPS_CLEAN_LIST=['met_em', 'CFS', 'ERA', 'SST', 'PFILE', 'GFS']
WRF_CLEAN_LIST=['met_em','wrfinput_','wrflowinp_','wrffda_','wrfrst_']

DRV_DIC={
    'cpsv3':{
        'use_ungrib':False,'atm_nfrq':'6H','atm_file_nfrq':'1D',
        'atm_naming_pattern':'$I%Y%m%d%H$I.cam2.h4.$F%Y-%m-%d$F-00000.nc',
        'lnd_naming_pattern':'$I%Y%m%d%H$I.clm2.h2.$F%Y-%m-%d$F-00000.nc'
    }
}
