#!/bin/sh
#PBS -lwalltime=15:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd $HOME/test/4hb_lettuce_12nt/rupture_4hb_lettuce12nt
inputfile="input_rupture"

runcommand="$HOME/oxDNA/build/bin/oxDNA  $inputfile"
timeout 895m $runcommand

