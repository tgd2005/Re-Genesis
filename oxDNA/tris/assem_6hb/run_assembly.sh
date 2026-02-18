#!/bin/sh
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd $HOME/test/assem_6hb
inputfile="input_assembly"

runcommand="$HOME/oxDNA/build/bin/oxDNA  $inputfile"
timeout 465m $runcommand

