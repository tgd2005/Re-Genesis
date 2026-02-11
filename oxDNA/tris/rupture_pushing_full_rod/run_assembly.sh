#!/bin/sh
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd $HOME/test/rupture_pushing_full_rod
inputfile="input_assembly"

runcommand="$HOME/oxDNA/build/bin/oxDNA  $inputfile"
timeout 465m $runcommand

