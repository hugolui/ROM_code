#!/bin/sh
#PBS -N dynamic_stall_v6
#PBS -l select=1:ngpus=1
#PBS -l walltime=12:00:00
#PBS -m ae
#PBS -M hugo.slui@gmail.com
#PBS -o nohup.dat
#PBS -e nohup_erro.dat
cd ${PBS_O_WORKDIR}
module load gcc
module load cuda-toolkit/9.0.176
module load cudnn/7.0
source $HOME/tensorflow/bin/activate 
python main.py
