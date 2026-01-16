#!/usr/bin/env zsh
set -euo pipefail

if (( $# != 2 )); then
  print "Usage: $0 <k> <r0>" >&2
  exit 1
fi

k="$1"
r0="$2"
dir="../K_${k}_r0_${r0}"

out_dir="$HOME/Re-Genesis/oxDNA/alex/force_dist_v2/force-analysis"
mkdir -p "$out_dir"

# Source ONLY the function file (not .zshrc)
source "./zsh_funcs_force.zsh"

cd "$dir"

files=(
  postass_rel_MC_dist.dat
  postass_rel_MD_dist.dat
  posttraps_rel_MC_dist.dat
  posttraps_long_MD_dist.dat
)

out="${out_dir}/force_mag_K_${k}_r0_${r0}.dat"
: > "$out"

last_global=""   # last stitched step written
first_file=1

for f in "${files[@]}"; do
  [[ -f "$f" ]] || { print "Missing: $dir/$f" >&2; exit 1; }

  tmp="$(mktemp)"
  fvec_mag "$f" > "$tmp"
  [[ -s "$tmp" ]] || { print "No data from fvec_mag for: $f" >&2; rm -f "$tmp"; exit 1; }

  first_step="$(awk 'NR==1{print $1; exit}' "$tmp")"
  last_raw="$(awk 'END{print $1}' "$tmp")"

  if (( first_file )); then
    offset=0
    skip_first=0
    first_file=0
  else
    # align: raw first step maps to last_global
    offset=$(( last_global - first_step ))
    skip_first=1
  fi

  # append, skipping first row for files 2+
  awk -v off="$offset" -v skip="$skip_first" 'NR==1 && skip {next} {print $1+off, $2}' \
    "$tmp" >> "$out"

  # update last_global based on what we just appended:
  # last_global = last_raw + offset (still valid even if we skipped first row)
  last_global=$(( last_raw + offset ))

  rm -f "$tmp"
done

print "Wrote: $out"
