#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:mem=16GB

cd $PBS_O_WORKDIR/scripts/r-9

module load Python/3.12.3-GCCcore-13.3.0

python r-9_nicked_twotubes.py


