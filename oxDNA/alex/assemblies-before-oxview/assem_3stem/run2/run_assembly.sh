#!/bin/sh
#PBS -lwalltime=02:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/assem_3stem/run2
inputfile="input_assembly"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 119m $runcommand

