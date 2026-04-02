#!/bin/sh
#PBS -lwalltime=00:20:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd $HOME/test/4hb_lettuce_12nt/assem_4hb_lettuce12nt
inputfile="input_relax_MD"

runcommand="$HOME/oxDNA/build/bin/oxDNA  $inputfile"
timeout 15m $runcommand

