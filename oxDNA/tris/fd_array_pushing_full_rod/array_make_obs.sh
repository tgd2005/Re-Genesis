#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# run_analysis.sh — Build a single obs.dat from MC→MD→MD_long distance logs
#
# #### PURPOSE
#   For each (k, r0) case, parse three oxDNA outputs and stitch them in time:
#     1) posttraps_rel_MC_dist.dat
#     2) posttraps_rel_MD_dist.dat
#     3) posttraps_long_MD_dist.dat
#   The result is obs.dat with 4 columns:
#     step  distance(15–16)  |F_ext_on_15|  hbonds
#
# #### ASSUMED INPUT LAYOUT
#   $ROOT/output/K_<kTag>_r0_<r0>/
#     ├─ K_<kTag>_r0_<r0>_traps_MC/        (posttraps_rel_MC_dist.dat)
#     ├─ K_<kTag>_r0_<r0>_traps_MD/        (posttraps_rel_MD_dist.dat)
#     └─ K_<kTag>_r0_<r0>_traps_MD_long/   (posttraps_long_MD_dist.dat; obs.dat is written here)
#
# #### DIST FILE FORMAT (per time block; comments may appear and are ignored)
#   Line A: "<step> <distance>    # t = <step>"
#   Line B: "<Fx15> <Fy15> <Fz15>"           # external force on particle 15
#   Line C: either "<Fx16> <Fy16> <Fz16>" THEN next line "<hbonds>"
#            or directly "<hbonds>"
#
# #### KEY BEHAVIOUR
#   • MD files restart their step counter at 0 → we keep a running step offset
#     so the final obs.dat time is continuous across MC→MD→MD_long.
#   • The very first block of any subsequent file (where raw step == 0) is
#     dropped to avoid duplicating the MC/MD join point.
#   • |F| is computed as sqrt(Fx15^2 + Fy15^2 + Fz15^2).
#
# #### CONFIGURATION
#   ROOT       : project root directory
#   K_VALUES   : array of k values (decimals); converted to tag "XYZ" via k*100
#   R0_VALUES  : array of r0 integers
#
# #### OUTPUT (per case)
#   <ROOT>/output/K_<kTag>_r0_<r0>/K_<kTag>_r0_<r0>_traps_MD_long/obs.dat
#
# #### HOW TO SET UP THE FUNCTION FOR NEW USERS
#   1) Adjust ROOT, K_VALUES, R0_VALUES as needed.
#   2) chmod +x run_analysis.sh
#   3) ./run_analysis.sh
#
# #### EXIT CODES
#   Exits non-zero if a required directory/file is missing for a given case;
#   other cases continue to be processed.
# ------------------------------------------------------------------------------



set -euo pipefail

ROOT="$HOME/test/fd_array_pushing_full_rod"

# choose what to process
K_VALUES=(0.35)
R0_VALUES=(2 3 4 5 6 7 8 9 10 11 12 13 14 16 18 20 22 24 26 28 30)

k_to_tag() {
  python3 - <<'PY' "$1"
import sys
k=float(sys.argv[1]); print(f"{int(round(k*100)):03d}")
PY
}

echo "==> ROOT: $ROOT"
[[ -d "$ROOT" ]] || { echo "ERROR: ROOT does not exist"; exit 1; }

for k in "${K_VALUES[@]}"; do
  k_tag="$(k_to_tag "$k")"
  for r0 in "${R0_VALUES[@]}"; do
    base="K_${k_tag}_r0_${r0}"
    outdir="${ROOT}/output/${base}"
    mc_src="${outdir}/${base}_traps_MC"
    md_src="${outdir}/${base}_traps_MD"
    long_src="${outdir}/${base}_traps_MD_long"
    workdir="$long_src"

    echo
    echo "============================================================"
    echo "==> k=$k (tag=$k_tag), r0=$r0"
    echo "==> WORKDIR: $workdir"

    # Require dirs and files
    for d in "$mc_src" "$md_src" "$long_src"; do
      [[ -d "$d" ]] || { echo "SKIP: missing dir $d"; continue 2; }
    done
    cd "$workdir"

    for f in \
      "${mc_src}/posttraps_rel_MC_dist.dat" \
      "${md_src}/posttraps_rel_MD_dist.dat" \
      "posttraps_long_MD_dist.dat"
    do
      [[ -f "$f" ]] || { echo "SKIP: missing $f"; continue 2; }
    done

    # Parse blocks → obs.dat
    awk '
      function mag(x,y,z){ return sqrt(x*x+y*y+z*z) }
      BEGIN{ OFS=" "; last_step=-1; file_idx=0; }

      # New file → remember we crossed a boundary
      FNR==1 { file_idx++; first_in_file=1; offset=(file_idx>1 ? last_step : 0) }

      # Start of a block: "step distance ..."
      $1 ~ /^-?[0-9]+$/ && NF>=2 {
        raw_step = $1 + 0
        dist     = $2 + 0.0

        # Next non-comment line = force on particle 15
        do { if (getline <= 0) exit } while ($1 ~ /^#/)
        fx=$1+0; fy=$2+0; fz=$3+0

        # Next non-comment: either opposite force (three floats) OR HB (single int)
        do { if (getline <= 0) exit } while ($1 ~ /^#/)
        if (NF==3) {
          # consume opposite-force line and read HB line
          do { if (getline <= 0) exit } while ($1 ~ /^#/)
        }
        hb = $1 + 0

        # Drop duplicated first block of each subsequent file (common MC→MD join)
        if (first_in_file && file_idx>1 && raw_step==0) { first_in_file=0; next }
        first_in_file=0

        step = raw_step + offset
        F = mag(fx,fy,fz)

        print step, dist, F, hb
        last_step = step
      }
    ' "${mc_src}/posttraps_rel_MC_dist.dat" \
       "${md_src}/posttraps_rel_MD_dist.dat" \
       "posttraps_long_MD_dist.dat" > obs.dat

    echo "==> Wrote: $workdir/obs.dat"
  done
done

echo
echo "All done."

