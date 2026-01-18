#!/usr/bin/env bash
# Usage:
#   ./make_obs_multi.sh obs_all.dat file1.dat file2.dat file3.dat [more.dat...]
#
# Each input file is an oxDNA data_output file like posttraps_*_dist.dat:
# header line: "<step> <distance> [# ...]"
# followed by two 3-vectors (we use the first = external force on particle 15).
#
# Output columns:
#   1) step  (continuous across files)
#   2) distance(15,16)
#   3) |F_ext(15)|

set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 OUTFILE file1.dat file2.dat file3.dat [more.dat...]" >&2
  exit 1
fi

OUT="$1"; shift

awk '
  # --- helpers ---
  function isnum(x){ return (x+0 == x) }
  function almosteq(a,b,eps,    d){ d=a-b; if(d<0)d=-d; return d<=eps }

  function emit(step,dist,fx,fy,fz,    mag){
    mag = sqrt(fx*fx + fy*fy + fz*fz)
    printf "%.0f %.10g %.10f\n", step, dist, mag
    last_emitted_step = step
    last_raw_step     = step_raw
    last_raw_dist     = dist
  }

  BEGIN{
    fileidx=0
    offset=0
    last_emitted_step=0
    last_raw_step=""
    last_raw_dist=""
  }

  # New file boundary
  FNR==1 {
    fileidx++
    wantF=0
    # Steps in the new file will be offset so time continues
    if (fileidx>1) offset = last_emitted_step
  }

  /^[[:space:]]*#/ || NF==0 { next }  # skip comments/blank lines

  {
    # Header line: step distance
    if (isnum($1) && isnum($2) && (NF==2 || $3 ~ /^#/)) {
      step_raw=$1
      dist    =$2
      step    = step_raw + offset

      # Detect duplicate-at-boundary: skip the first block of a file
      # if it matches the last block of the previous file (within tol).
      skip_this_block = (fileidx>1) && (step_raw==last_raw_step) && almosteq(dist,last_raw_dist,1e-9)

      wantF=1
      next
    }

    # First 3-vector after header is the external force on particle 15
    if (wantF && NF>=3 && isnum($1) && isnum($2) && isnum($3)) {
      if (!skip_this_block) {
        emit(step, dist, $1, $2, $3)
      }
      wantF=0
    }
  }
' "$@" > "$OUT"

echo "Wrote $OUT  (cols: step  distance_15_16  |F_ext(15)|)"

