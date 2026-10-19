#!/bin/sh
#PBS -lwalltime=48:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd "$PBS_O_WORKDIR"
inputfile="input_VMMC"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 47m $runcommand

