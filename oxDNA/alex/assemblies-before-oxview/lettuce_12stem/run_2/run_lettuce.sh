#!/bin/sh
#PBS -lwalltime=10:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/lettuce/run_2/
inputfile="input_lettuce"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 545m $runcommand

