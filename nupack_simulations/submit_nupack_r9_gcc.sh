#!/bin/sh -l
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:mem=16GB
#PBS -N npk_twotubes


cd $PBS_O_WORKDIR/scripts/r-9
module purge
module load tools/prod
module load Python/3.11.5-GCCcore-13.2.0
source /rds/general/user/ag4723/home/Re-Genesis/nupack_simulations/nupack/bin/activate
python3 r-9_nicked_twotubes.py


