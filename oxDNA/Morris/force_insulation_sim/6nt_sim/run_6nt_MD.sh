#!/bin/sh
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=4:mem=10GB

cd /rds/general/user/tc2523/home/force_insulation_sim/6nt_sim
inputfile="6nt_input.txt"

export OMP_NUM_THREADS=4

runcommand="/rds/general/user/tc2523/home/oxDNA/oxDNA2/oxDNA/build/bin/oxDNA  $inputfile"
timeout 23h $runcommand

