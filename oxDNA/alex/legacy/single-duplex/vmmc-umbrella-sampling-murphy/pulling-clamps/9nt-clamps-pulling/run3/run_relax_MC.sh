#!/bin/sh
#PBS -lwalltime=10:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/oxview_9nt/run3
inputfile="input_relax_MC"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 599m $runcommand

