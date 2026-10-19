#!/bin/sh
#PBS -lwalltime=00:20:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/ag4723/home/Re-Genesis/oxDNA/alex/assem_40rod/fixed_external/
inputfile="input_relax_MD"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 15m $runcommand

