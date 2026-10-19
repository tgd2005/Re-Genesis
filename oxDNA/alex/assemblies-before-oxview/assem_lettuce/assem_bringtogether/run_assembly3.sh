#!/bin/sh
#PBS -lwalltime=02:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/assem_lettuce/assem_bringtogether/
inputfile="input_assembly3"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 115m $runcommand

