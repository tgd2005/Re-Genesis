#!/usr/bin/env bash


set -euo pipefail

PBS_SCRIPT="run_array_traps_MD_long2.pbs"  # change if your file has a different name
TOTAL=120
BATCH=10

for start in $(seq 1 $BATCH $TOTAL); do
  remaining=$((TOTAL - start + 1))
  n=$(( remaining < BATCH ? remaining : BATCH ))

  echo "Submitting rows $start..$((start + n - 1)) as array 1-$n"
  qsub -v START="$start" -J 1-"$n" "$PBS_SCRIPT"
done

