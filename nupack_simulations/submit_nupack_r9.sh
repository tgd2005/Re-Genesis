#!/bin/sh
#PBS -l walltime=09:00:00
#PBS -l select=1:ncpus=16:mem=128GB

cd /rds/general/user/ag4723/home/Re-Genesis/nupack_simulations
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate nupack
runcommand="python3 $HOME/Re-Genesis/nupack_simulations/r-9_nicked_onetube.py"
timeout 8h $runcommand
