#!/bin/sh
#PBS -lwalltime=05:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/lettuce/assembly
inputfile="input_relax_MC"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 229m $runcommand

