#!/bin/bash

ID=$1
PID=$2
YR=$3
DAYS=$4
HOST=hqlx$ID
CTRL=/home/lzhenn/array74/workspace/uranus/case_ctrler/cmip6-down

#conda activate uranus-test
cd ${CTRL}
python era5drv.py -m ${ID} -p ${PID} -y ${YR} -d ${DAYS}
