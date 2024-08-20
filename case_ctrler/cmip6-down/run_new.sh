#!/bin/bash

ID=$1
PID=$2
YR=$3
DAYS=$4
SSP=$5
HOST=hqlx$ID
CTRL=`pwd`
ssh $HOST "conda activate uranus-test;cd ${CTRL};python cmip6_wrf_new.py -m ${ID} -p ${PID} -y ${YR} -d ${DAYS} -s ${SSP}"
