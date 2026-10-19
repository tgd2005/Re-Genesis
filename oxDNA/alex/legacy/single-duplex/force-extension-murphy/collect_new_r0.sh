#!/usr/bin/env bash
# Copies obs_dist.dat files for the new r0=32-38 runs (K=0.25, seeds 1-10)
# into fextensions_results/fd_pull_full_rod_dats/ with the same naming convention
# as the existing fullrod results.
#
# Run from the force-extension/ directory AFTER array_make_obs_new_r0.sh has completed.
#
# Usage:
#   cd /path/to/force-extension
#   bash collect_new_r0.sh
set -euo pipefail

ANALYSIS="fd_array_fullrod_kinking_seeds/analysis_files"
DEST="fextensions_results/fd_pull_full_rod_dats"
mkdir -p "$DEST"

for seed in 1 2 3 4 5 6 7 8 9 10; do
  for r0 in 32 33 34 35 36 37 38; do
    src="${ANALYSIS}/Seed_${seed}/K_025/r0_${r0}/obs_dist.dat"
    dst="${DEST}/obs_dist_K_025_r0_${r0}_seed_${seed}.dat"

    if [[ -f "$src" ]]; then
      cp -n -- "$src" "$dst"
      echo "Copied: $dst"
    else
      echo "MISSING: $src" >&2
    fi
  done
done
