#!/bin/sh
#PBS -lwalltime=06:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd $HOME/test/assem_pushing_3q_rod
inputfile="input_assembly"

runcommand="$HOME/oxDNA/build/bin/oxDNA  $inputfile"
timeout 345m $runcommand

