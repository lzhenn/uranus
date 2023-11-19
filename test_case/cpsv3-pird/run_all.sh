#!/bin/bash
#SBATCH --comment=WRF
#SBATCH -J step1_coawst_s2s 
#SBATCH -n 1
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-socket=16
#SBATCH -o %j
#SBATCH -e %j
#SBATCH -p serial 
#SBATCH -t 10:00:00


source /g6/cmme/COAWST-S2S/bashrc_coawst_run 
conda activate uranus
python3 ./cpsv3.py
