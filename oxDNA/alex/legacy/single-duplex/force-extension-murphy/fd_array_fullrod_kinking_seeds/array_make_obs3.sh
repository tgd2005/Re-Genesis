#!/bin/sh
# POSIX version: stitches MC→MD→MD_long dist logs, writes obs_seed_<seed>.dat,
# and collects them into $ROOT/analysis_files/.

set -eu

# --- project root (edit if needed)
ROOT="$HOME/Re-Genesis/oxDNA/alex/fd_array_fullrod_kinking_seeds"
ANALYSIS_DIR="$ROOT/analysis_files"
mkdir -p "$ANALYSIS_DIR"

# what to process
K_VALUES="0.25 0.35 0.45"
R0_VALUES="2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31"

# helper: 0.35 -> 035
k_to_tag() {
  awk -v k="$1" 'BEGIN{printf "%03d", int(k*100+0.5)}'
}

echo "==> ROOT: $ROOT"
[ -d "$ROOT" ] || { echo "ERROR: ROOT does not exist" >&2; exit 1; }

for k in $K_VALUES; do
  k_tag="$(k_to_tag "$k")"

  for r0 in $R0_VALUES; do
    r0_tag="$(r0_to_tag "$r0")"   # if you have it; otherwise just r0_tag="$r0"

    found=0

    # loop over seed directories
    for seed_dir in "$outroot"/Seed_*/; do
      [ -d "$seed_dir" ] || continue
      found=1

      seed_tag="${seed_dir%/}"
      seed_tag="${seed_tag##*/}"      # Seed_123
      seed_tag="${seed_tag#Seed_}"    # 123

      base_dir="${seed_dir%/}/K_${k_tag}/r0_${r0_tag}"

      # candidate run dirs
      mc_src="${base_dir}/traps_MC"
      md_src="${base_dir}/traps_MD"
      long_src="${base_dir}/traps_MD_long"
      echo
      echo "============================================================"
      echo "==> k=$k (tag=$k_tag), r0=$r0, seed=$seed"
      echo "==> WORKDIR: $workdir"

      [ -d "$mc_src"   ] || { echo "SKIP: missing dir $mc_src";   continue; }
      [ -d "$md_src"   ] || { echo "SKIP: missing dir $md_src";   continue; }
      [ -d "$long_src" ] || { echo "SKIP: missing dir $long_src"; continue; }

      [ -f "$mc_src/posttraps_rel_MC_dist.dat"      ] || { echo "SKIP: missing MC dist";      continue; }
      [ -f "$md_src/posttraps_rel_MD_dist.dat"      ] || { echo "SKIP: missing MD dist";      continue; }
      [ -f "$long_src/posttraps_long_MD_dist.dat"   ] || { echo "SKIP: missing MD_long dist"; continue; }

      cd "$workdir"

      # Parse blocks into obs (step, dist, |F15|, hb); keep time continuous; drop first block of later files
      awk '
        function mag(x,y,z){ return sqrt(x*x+y*y+z*z) }
        BEGIN{ OFS=" "; last_step=-1; file_idx=0; }
        FNR==1 { file_idx++; first_in_file=1; offset=(file_idx>1 ? last_step : 0) }
        $1 ~ /^-?[0-9]+$/ && NF>=2 {
          raw_step = $1 + 0
          dist     = $2 + 0.0

          # next non-comment → force on particle 15
          do { if (getline <= 0) exit } while ($1 ~ /^#/)
          fx=$1+0; fy=$2+0; fz=$3+0

          # next non-comment: either opposite force (3 cols) then HB, or directly HB
          do { if (getline <= 0) exit } while ($1 ~ /^#/)
          if (NF==3) {
            do { if (getline <= 0) exit } while ($1 ~ /^#/)
          }
          hb = $1 + 0

          # drop duplicated first block of later files
          if (first_in_file && file_idx>1 && raw_step==0) { first_in_file=0; next }
          first_in_file=0

          step = raw_step + offset
          F = mag(fx,fy,fz)
          print step, dist, F, hb
          last_step = step
        }
      ' "$mc_src/posttraps_rel_MC_dist.dat" \
         "$md_src/posttraps_rel_MD_dist.dat" \
         "$long_src/posttraps_long_MD_dist.dat" > "obs_seed_${seed}.dat"

      # collect for GitHub
      cp -f "obs_seed_${seed}.dat" \
        "$ANALYSIS_DIR/obs_K_${k_tag}_r0_${r0}_seed_${seed}.dat"
      echo "==> Collected → $ANALYSIS_DIR/obs_K_${k_tag}_r0_${r0}_seed_${seed}.dat"
    done

    [ "$found" -eq 1 ] || echo "SKIP: no seed directories found for k=${k_tag}, r0=${r0}"
  done
done

echo
echo "All done. Collected files are in: $ANALYSIS_DIR"

