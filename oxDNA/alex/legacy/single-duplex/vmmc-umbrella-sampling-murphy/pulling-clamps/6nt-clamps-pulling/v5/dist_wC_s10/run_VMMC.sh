#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=1:mem=10GB
cd "$PBS_O_WORKDIR"
timeout 71h $HOME/oxDNA_git/oxDNA/build/bin/oxDNA input_VMMC
