#!/bin/sh
#PBS -lwalltime=15:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/assem_lettuce/fixed_external/
inputfile="input_assembly2"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 850m $runcommand

