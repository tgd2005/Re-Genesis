#!/usr/bin/env bash
set -euo pipefail


ROOT="$HOME/Re-Genesis/oxDNA/alex/fd_array_halfrod_kinking_seeds"
ANALYSIS_DIR="$ROOT/analysis_files"
mkdir -p "$ANALYSIS_DIR"

# what to process
K_VALUES="0.25"
R0_VALUES="2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31"

# helper: 0.35 -> 035
k_to_tag() {
  awk -v k="$1" 'BEGIN{printf "%03d", int(k*100+0.5)}'
}

# helper: r0 tag (identity for your r0 list)
r0_to_tag() { echo "$1"; }

echo "==> ROOT: $ROOT"
[ -d "$ROOT" ] || { echo "ERROR: ROOT does not exist" >&2; exit 1; }

outroot="$ROOT/output"
[ -d "$outroot" ] || { echo "ERROR: output/ not found under ROOT: $outroot" >&2; exit 1; }

for k in $K_VALUES; do
  k_tag="$(k_to_tag "$k")"

  for r0 in $R0_VALUES; do
    r0_tag="$(r0_to_tag "$r0")"

    found_seed=0

    for seed_dir in "$outroot"/Seed_*/; do
      [[ -d "$seed_dir" ]] || continue
      found_seed=1

      seed="${seed_dir%/}"
      seed="${seed##*/}"     # Seed_XXX
      seed="${seed#Seed_}"   # XXX

      base_dir="${seed_dir%/}/K_${k_tag}/r0_${r0_tag}"
      [[ -d "$base_dir" ]] || continue

      mc_src="${base_dir}/traps_MC"
      md_src="${base_dir}/traps_MD"
      long_src="${base_dir}/traps_MD_long"

      # Mirror output structure in analysis dir
      dest_dir="$ANALYSIS_DIR/Seed_${seed}/K_${k_tag}/r0_${r0_tag}"
      mkdir -p "$dest_dir"

      ############################################################
      # 1) DIST observable: step dist |F15| hb force_energy
      ############################################################
      dist_files=()
      [[ -f "$mc_src/posttraps_rel_MC_dist.dat"    ]] && dist_files+=("$mc_src/posttraps_rel_MC_dist.dat")
      [[ -f "$md_src/posttraps_rel_MD_dist.dat"    ]] && dist_files+=("$md_src/posttraps_rel_MD_dist.dat")
      [[ -f "$long_src/posttraps_long_MD_dist.dat" ]] && dist_files+=("$long_src/posttraps_long_MD_dist.dat")

      if ((${#dist_files[@]} > 0)); then
        awk '
          function mag(x,y,z){ return sqrt(x*x+y*y+z*z) }

          # Read next non-comment line into $0/$1..$NF; return 1 on success else 0
          function get_noncomment(   ok){
            while ((ok = getline) > 0) {
              if ($0 ~ /^[[:space:]]*#/) continue
              if ($0 ~ /^[[:space:]]*$/) continue
              return 1
            }
            return 0
          }

          BEGIN{
            OFS=" ";
            last_step=-1;
            file_idx=0;
          }

          FNR==1{
            file_idx++;
            first_in_file=1;
            offset = (file_idx>1 ? last_step : 0);
          }

          # timestep line: integer step + dist (allow trailing comments)
          $1 ~ /^-?[0-9]+$/ && NF>=2 {
            raw_step = $1 + 0
            dist     = $2 + 0.0

            # Force on particle 15 (next non-comment line)
            if (!get_noncomment()) exit
            fx=$1+0; fy=$2+0; fz=$3+0

            # Opposite force vector (next non-comment line), ignored
            if (!get_noncomment()) exit

            # Last line contains hb and force_energy (second value)
            if (!get_noncomment()) exit
            hb = $1 + 0
            force_energy = (NF>=2 ? $2 + 0.0 : 0.0)

            # Drop duplicated first block of later files (restart at t=0)
            if (first_in_file && file_idx>1 && raw_step==0) { first_in_file=0; next }
            first_in_file=0

            step = raw_step + offset
            F = mag(fx,fy,fz)
            print step, dist, F, hb, force_energy
            last_step = step
          }
        ' "${dist_files[@]}" > "${dest_dir}/obs_dist.dat"
      fi

      ############################################################
      # 2) PAIR ENERGY observable: step CXSTCK(0,85) CXSTCK(31,32)
      ############################################################
      pair_files=()
      [[ -f "$mc_src/posttraps_rel_MC_pair_energy.dat"    ]] && pair_files+=("$mc_src/posttraps_rel_MC_pair_energy.dat")
      [[ -f "$md_src/posttraps_rel_MD_pair_energy.dat"    ]] && pair_files+=("$md_src/posttraps_rel_MD_pair_energy.dat")
      [[ -f "$long_src/posttraps_long_MD_pair_energy.dat" ]] && pair_files+=("$long_src/posttraps_long_MD_pair_energy.dat")

      if ((${#pair_files[@]} > 0)); then
        awk '
          BEGIN{
            OFS=" ";
            file_idx=0;
            last_t=-1;
            offset=0;
            have_block=0;
          }

          function flush_block(  t_out){
            if (!have_block) return;

            # drop duplicated first block of later files if t==0
            if (first_in_file && file_idx>1 && cur_t==0) {
              first_in_file=0;
              have_block=0;
              v085="nan";
              v3132="nan";
              return;
            }
            first_in_file=0;

            t_out = cur_t + offset;
            print t_out, v085, v3132;
            last_t = t_out;

            have_block=0;
            v085="nan";
            v3132="nan";
          }

          FNR==1{
            file_idx++;
            first_in_file=1;
            offset = (file_idx>1 ? last_t : 0);
            have_block=0;
            v085="nan";
            v3132="nan";
          }

          # timestep header line: contains "t = NUMBER"
          /^[[:space:]]*#/{
            if (match($0, /t[[:space:]]*=[[:space:]]*([0-9]+)/, m)) {
              flush_block();
              cur_t = m[1] + 0;
              have_block=1;
              next;
            }
            next;
          }

          # Data lines: id1 id2 ... CXSTCK ... (CXSTCK is column 9)
          have_block && NF>=9 {
            id1=$1+0; id2=$2+0;
            cx=$9+0.0;

            if ((id1==0  && id2==85) || (id1==85 && id2==0)) v085=cx;
            if ((id1==31 && id2==32) || (id1==32 && id2==31)) v3132=cx;
            next;
          }

          END{ flush_block(); }
        ' "${pair_files[@]}" > "${dest_dir}/obs_cxstck.dat"
      fi

      echo "==> Seed_${seed} K_${k_tag} r0_${r0_tag} -> ${dest_dir}"
    done

    if [[ "$found_seed" -eq 0 ]]; then
      echo "WARNING: No Seed_* directories found under $outroot" >&2
    fi
  done
done
