#!/usr/bin/env bash
# Usage:
#   ./make_obs_mc_md.sh posttraps_rel_MC_dist.dat posttraps_rel_MD_dist.dat obs_mc_md.dat
#
# Output columns:
#   1) step (MD steps are offset to continue after MC)
#   2) distance between 15 and 16
#   3) |F_ext(15)| magnitude

set -euo pipefail
MC="${1:-posttraps_rel_MC_dist.dat}"
MD="${2:-posttraps_rel_MD_dist.dat}"
OUT="${3:-obs_mc_md.dat}"

awk '
  function isnum(x){ return (x+0==x) }

  function emit(step,dist,fx,fy,fz,   mag,row){
    mag=sqrt(fx*fx+fy*fy+fz*fz)
    row=sprintf("%.0f %.10g %.10f", step, dist, mag)
    print row
    last_emitted_step = step
  }

  BEGIN{
    fileidx=0
    offset=0
    skipping_first_md_block=0
  }

  FNR==1 {                           # New file boundary
    fileidx++
    wantF=0
    if (fileidx==2) {
      # For the MD file: offset steps by the last MC step,
      # and skip the first MD block (duplicate of last MC block)
      offset = last_emitted_step
      skipping_first_md_block = 1
    }
  }

  /^[[:space:]]*#/ || NF==0 { next }  # ignore comments/blank

  {
    # Header line: "<step> <distance> [# ...]"
    if (isnum($1) && isnum($2) && (NF==2 || $3 ~ /^#/)) {
      step_raw=$1; dist=$2
      step = step_raw + offset
      wantF=1
      next
    }

    # First 3-number line after header = force on particle 15
    if (wantF && NF>=3 && isnum($1) && isnum($2) && isnum($3)) {
      if (fileidx==2 && skipping_first_md_block==1) {
        # consume this first MD block but DO NOT emit it
        skipping_first_md_block=0
      } else {
        emit(step, dist, $1, $2, $3)
      }
      wantF=0
    }
  }
' "$MC" "$MD" > "$OUT"

echo "Wrote $OUT  (cols: step  distance_15_16  |F_ext(15)|)"

