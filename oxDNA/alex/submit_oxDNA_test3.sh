#!/bin/sh
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-3
module load fix_unwritable_tmp

IND="$PBS_ARRAY_INDEX"

#MODIFY THE FOLLOWING 2 LINES TO POINT TO THE DIRECTORY WITH THE INPUT FILE, DOESNT HAVE TO BE TWO LINES WITH BASE, CAN BE ONE
BASE=/rds/general/user/ag4723/home/VMMC_MELTING_WD/test
INPUT_DIR="${BASE}/test_$IND/test_${IND}_input"

#CD INTO INPUT DIRECTORY - MANDATORY FOR OXDNA
cd "$INPUT_DIR"

#POINT TO THE INPUT FILE - MANDATORY FOR OXDNA
inputfile="${INPUT_DIR}/input"

#DEFINE THE RUNCOMMAND AND RUN WITH TIMEOUT FOR SAFE CRASHING
runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA $inputfile"
timeout 8m $runcommand

# After oxDNA finishes, collect all .dat files in a separate output directory within the input directory
outdir="../test_${IND}_output"
mkdir -p "$outdir"
for f in ./*.dat; do
     [ -e "$f" ] || break

if [ "$(basename "$f")" = "initial.dat" ]; then
continue
fi

     mv "$f" "$outdir"/
done
