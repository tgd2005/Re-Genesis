#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./collect_obs.sh /path/to/search_root
# If you omit the argument, it searches in the current directory.

root="${1:-.}"
outroot="fextension_results"

# Collect both obs_dist.dat and obs_cxstck.dat by default.
# If you ONLY want dist, change -name 'obs_*.dat' to -name 'obs_dist.dat'
find "$root" -type f -name 'obs_dist.dat' -print0 |
while IFS= read -r -d '' src; do
  # Make relative path (best-effort) so "highest level" is meaningful
  rel="${src#"$root"/}"

  # Highest-level folder = first path component
  top="${rel%%/*}"

  # Extract seed and r0 from the path
  if [[ "$rel" =~ /Seed_([0-9]+)/ ]]; then
    seed="${BASH_REMATCH[1]}"
  else
    echo "Skipping (no Seed_N in path): $src" >&2
    continue
  fi

  if [[ "$rel" =~ /r0_([0-9]+)/ ]]; then
    r0="${BASH_REMATCH[1]}"
  else
    echo "Skipping (no r0_N in path): $src" >&2
    continue
  fi

  # Extract obs type from filename: obs_dist.dat or obs_cxstck.dat
  base="$(basename "$src")"
  if [[ "$base" =~ ^obs_([^.]+)\.dat$ ]]; then
    obs="${BASH_REMATCH[1]}"
  else
    echo "Skipping (unexpected filename): $src" >&2
    continue
  fi

  # Destination directory and filename
  dest_dir="${outroot}/${top}"
  mkdir -p "$dest_dir"

  # Include obs type to avoid overwriting when both dist and cxstck exist
  dest_name="obs_${obs}_K_025_r0_${r0}_seed_${seed}.dat"
  dest="${dest_dir}/${dest_name}"

  # Copy without overwriting existing files
  cp -n -- "$src" "$dest" || true
done
