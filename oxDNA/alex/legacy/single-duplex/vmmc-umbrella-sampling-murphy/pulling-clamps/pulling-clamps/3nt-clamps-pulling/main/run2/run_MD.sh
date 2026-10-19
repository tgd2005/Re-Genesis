#!/bin/sh
#PBS -lwalltime=20:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/oxview_3nt/main/run2
inputfile="input_MD"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 1199m $runcommand

