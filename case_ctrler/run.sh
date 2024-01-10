#!/bin/bash

ID=$1
PID=$2
YR=$3
DAYS=$4
SSP=$5
HOST=hqlx$ID
CTRL=/home/lzhenn/array74/workspace/uranus/case_ctrler

ssh $HOST "conda activate uranus;cd ${CTRL};python cmip6_wrf.py -m ${ID} -p ${PID} -y ${YR} -d ${DAYS} -s ${SSP}"
