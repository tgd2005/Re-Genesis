#!/bin/sh
#PBS -lwalltime=01:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd $HOME/test/assem_4hb_sprouting
inputfile="input_assembly"

runcommand="$HOME/oxDNA/build/bin/oxDNA  $inputfile"
timeout 55m $runcommand

