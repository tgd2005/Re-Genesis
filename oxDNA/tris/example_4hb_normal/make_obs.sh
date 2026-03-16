#!/usr/bin/env bash
# Usage:
#   ./make_obs_multi.sh -o obs.dat file1.dist file2.dist [file3.dist ...]
set -euo pipefail

OUT="obs.dat"
if [[ "${1:-}" == "-o" ]]; then OUT="$2"; shift 2; fi
(( $# >= 1 )) || { echo "No input files"; exit 1; }

awk '
function mag(x,y,z){ return sqrt(x*x + y*y + z*z) }
BEGIN{ OFS=" "; last_step=-1 }

# New file: set offset so steps continue from previous file
FNR==1 { first_in_file=1; offset=(last_step>=0? last_step : 0) }

# Block starts on a numeric step line: step dist
$1 ~ /^-?[0-9]+$/ {
  step=$1; dist=$2;

  # next non-comment line: force on particle 15 (fx fy fz)
  do { if (getline <= 0) exit } while ($1 ~ /^#/)
  fx=$1; fy=$2; fz=$3;

  # next non-comment line: either opposite force (3 cols) or HB line
  do { if (getline <= 0) exit } while ($1 ~ /^#/)
  if (NF==3) {
    # opposite force line present -> read one more non-comment line for HB
    if (getline <= 0) exit
    while ($1 ~ /^#/) { if (getline <= 0) exit }
  }

  hb = $1 + 0
  F = mag(fx,fy,fz)

  # drop duplicate first block of subsequent files if they start at step 0
  if (first_in_file && last_step>=0 && step==0) { first_in_file=0; next }
  first_in_file=0

  s = step + offset
  print s, dist, F, hb
  last_step = s
}
' "$@" > "$OUT"

echo "Wrote $OUT"

