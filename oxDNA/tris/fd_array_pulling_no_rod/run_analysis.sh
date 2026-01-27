#!/usr/bin/env bash
# run_analysis.sh
#
# Runs the analysis for all (k_tag, r0_tag) combinations.
#
# We start at:
#   /Users/murphy92/Re-Genesis/oxDNA/alex/force_dist_array_fullrod
#
# Directory layout assumed:
#   output/K_{k_tag}_r0_{r0_tag}/K_{k_tag}_r0_{r0_tag}_traps_MD_long
#   output/K_{k_tag}_r0_{r0_tag}/K_{k_tag}_r0_{r0_tag}_traps_MC
#   output/K_{k_tag}_r0_{r0_tag}/K_{k_tag}_r0_{r0_tag}_traps_MD
#
# Where folder names use tags like:
#   K_035_r0_6   (k_tag is k*100 formatted as 3 digits; r0_tag is integer)
#
# For each combo, inside *_MD_long we:
#  - copy template/add_hb_column_to_dist.sh
#  - copy posttraps_rel_MC_* from ../*_traps_MC
#  - copy posttraps_rel_MD_* from  *_traps_MD
#  - generate *_dist_with_hb.dat for rel_MC, rel_MD, long_MD
#  - run make_obs.sh to create obs_with_hb.dat

set -euo pipefail

ROOT="$HOME/test/fd_array_pulling_no_rod"
TEMPLATE="${ROOT}/template/add_hb_column_to_dist.sh"

# k values (as decimals) -> tag format like 0.35 -> 035
K_VALUES=(0.35 0.45)

# r0 values (integers)
R0_VALUES=(2 3 4 5 6 7 8 9 10 11 12 13 14 16 18 20 22 24 26 28 30)

# Convert decimal like 0.35 to "035"
k_to_tag() {
  local k="$1"
  # Use python for correct rounding and zero-padding
  python3 - <<'PY' "$k"
import sys
k=float(sys.argv[1])
tag=int(round(k*100))
print(f"{tag:03d}")
PY
}

echo "==> ROOT: $ROOT"
[[ -d "$ROOT" ]] || { echo "ERROR: ROOT does not exist"; exit 1; }
[[ -f "$TEMPLATE" ]] || { echo "ERROR: Missing template script: $TEMPLATE"; exit 1; }

for k in "${K_VALUES[@]}"; do
  k_tag="$(k_to_tag "$k")"

  for r0 in "${R0_VALUES[@]}"; do
    base="K_${k_tag}_r0_${r0}"
    outdir="${ROOT}/output/${base}"
    workdir="${outdir}/${base}_traps_MD_long"
    mc_src="${outdir}/${base}_traps_MC"
    md_src="${outdir}/${base}_traps_MD"

    echo
    echo "============================================================"
    echo "==> Processing: k=$k (tag=$k_tag), r0=$r0  =>  $base"
    echo "==> WORKDIR: $workdir"

    # Skip gracefully if a combo is missing
    if [[ ! -d "$workdir" ]]; then
      echo "SKIP: missing workdir: $workdir"
      continue
    fi
    if [[ ! -d "$mc_src" ]]; then
      echo "SKIP: missing MC source dir: $mc_src"
      continue
    fi
    if [[ ! -d "$md_src" ]]; then
      echo "SKIP: missing MD source dir: $md_src"
      continue
    fi

    cd "$workdir"

    echo "==> Copying add_hb_column_to_dist.sh into workdir"
    cp -f "$TEMPLATE" ./add_hb_column_to_dist.sh
    chmod +x ./add_hb_column_to_dist.sh

    echo "==> Copying MC posttraps_rel_MC_* into workdir"
    cp -f "${mc_src}"/posttraps_rel_MC_* . || true

    echo "==> Copying MD posttraps_rel_MD_* into workdir"
    cp -f "${md_src}"/posttraps_rel_MD_* . || true

    # Check required inputs; if missing, skip this combo
    req=(
      "posttraps_rel_MC_dist.dat"
      "posttraps_rel_MC_hb.dat"
      "posttraps_rel_MD_dist.dat"
      "posttraps_rel_MD_hb.dat"
      "posttraps_long_MD_dist.dat"
      "posttraps_long_MD_hb.dat"
      "make_obs.sh"
    )
    missing=0
    for f in "${req[@]}"; do
      if [[ ! -f "$f" ]]; then
        echo "MISSING: $workdir/$f"
        missing=1
      fi
    done
    if [[ "$missing" -eq 1 ]]; then
      echo "SKIP: missing required inputs for $base"
      continue
    fi

    echo "==> Building dist_with_hb files"
    ./add_hb_column_to_dist.sh posttraps_rel_MC_dist.dat  posttraps_rel_MC_hb.dat  posttraps_rel_MC_dist_with_hb.dat
    ./add_hb_column_to_dist.sh posttraps_rel_MD_dist.dat  posttraps_rel_MD_hb.dat  posttraps_rel_MD_dist_with_hb.dat
    ./add_hb_column_to_dist.sh posttraps_long_MD_dist.dat posttraps_long_MD_hb.dat posttraps_long_MD_dist_with_hb.dat

    echo "==> Running make_obs.sh -> obs_with_hb.dat"
    chmod +x ./make_obs.sh
    ./make_obs.sh -o obs_with_hb.dat \
      posttraps_rel_MC_dist_with_hb.dat \
      posttraps_rel_MD_dist_with_hb.dat \
      posttraps_long_MD_dist_with_hb.dat

    echo "==> OK: wrote $workdir/obs_with_hb.dat"
  done
done

echo
echo "All done."

