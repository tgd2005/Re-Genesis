#!/bin/sh
#PBS -lwalltime=00:20:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/oxview_3nt/vmmc/run1
inputfile="input_VMMC"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 19m $runcommand

