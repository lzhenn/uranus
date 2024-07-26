#!/bin/sh
#------------------------------------------------------
# GFS WAVE DOWNLOAD SCRIPT 
#                                    Zhenning LI
#                                   Feb 19, 2022
#------------------------------------------------------

#https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20220207/00/wave/gridded/gfswave.t00z.global.0p25.f000.grib2
#
#
# Usage:
# sh gfs_wave_down_fcst_subdomain.sh 2024043012 /home/lzhenn/array181/op_njord/data/drv/wav/2024042912 7 
## ------------Below for user-defined configurations ------------
# Start time 
STRT_YMDH=$1
#STRT_YMDH=2022021812

# Archive path
ARCH_PATH=$2
#ARCH_PATH=/home/metctm1/array/data/gfs_wave/2022021812

# How long period to fecth
FCST_DAY=$3
#FCST_DAY=1

# The interval to fetch GFS output, 3-hr preferred, 
# 1-hr minimum, and no longer than 6-hr.
FRQ=3

# Init hour
INIT_HOUR=${STRT_YMDH:8:2}
# Resolution: 0p25, 0p16
RES=0p16

LON_LEFT=100
LON_RIGHT=140
LAT_TOP=50
LAT_BOTTOM=5


# try time interval in seconds
TRY_INTERVAL=120

# ------------Upper for user-defined configurations ------------

FETCH_DAY=$(date -d "${STRT_YMDH:0:8}" +%Y%m%d)
TODAY=$(date +%Y%m%d)
TIME_DELTA=`expr $TODAY - $FETCH_DAY`
if [ ! -d $ARCH_PATH ]; then
    mkdir $ARCH_PATH
fi

# not realtime, try get from envf server
#if [ $TIME_DELTA -gt 7 ]; then
#    echo "not realtime, exit"
#    exit
#fi

TOTAL_HR=`expr $FCST_DAY \* 24`


# fetch from ncep server
BASE_URL="https://nomads.ncep.noaa.gov/cgi-bin/filter_gfswave.pl"
VAR_FILTER="&var_DIRPW=on&var_HTSGW=on&var_PERPW=on"
SUBR_FILTER="&subregion=&leftlon="${LON_LEFT}"&rightlon="${LON_RIGHT}"&toplat="${LAT_TOP}"&bottomlat="${LAT_BOTTOM}
SUFFIX="&dir=%2Fgfs."${STRT_YMDH:0:8}"%2F"${INIT_HOUR}"%2Fwave%2Fgridded"

for CURR_HR in $(seq 0 $FRQ $TOTAL_HR) 
do
    TRY_TIME=10

    TSTEP=`printf "%03d" $CURR_HR`

    FN_FILTER="?file=gfswave.t"${INIT_HOUR}"z.global."${RES}".f"${TSTEP}".grib2&all_lev=on"
    FN="gfswave.t"${INIT_HOUR}"z.global."${RES}".f"${TSTEP}".grib2"
    if [ $RES == "0p25" ]; then
        NOM_SIZE=11000000
    else
        NOM_SIZE=110000
    fi
    SRC_URL=${BASE_URL}${FN_FILTER}${VAR_FILTER}${SUBR_FILTER}${SUFFIX}

    if [ -f ${ARCH_PATH}/${FN} ]; then
        FILESIZE=$(stat -c%s "${ARCH_PATH}/${FN}")
        if (( $FILESIZE>$NOM_SIZE )); then
            echo "File downloaded, skip "${FN}
            continue
        fi
    fi

    while (( $TRY_TIME>0 ))
    do
        wget ${SRC_URL} -O ${ARCH_PATH}/${FN}
        if [[ "$?" == 0 ]]; then
            FILESIZE=$(stat -c%s "${ARCH_PATH}/${FN}")
            if (( $FILESIZE<$NOM_SIZE )); then
                echo "File size not correct, try again in "$TRY_INTERVAL"s..."
                TRY_TIME=`expr $TRY_TIME - 1`
                sleep $TRY_INTERVAL
            else
                break
            fi
        else
            echo "Error downloading file ${FN}, try again in "$TRY_INTERVAL"s..."
            TRY_TIME=`expr $TRY_TIME - 1`
            sleep $TRY_INTERVAL
        fi
    done
    if (( $TRY_TIME == 0)); then
        echo "Download "${ARCH_PATH}/${FN} "failed after maximum attampts!!!"
    fi
done 
