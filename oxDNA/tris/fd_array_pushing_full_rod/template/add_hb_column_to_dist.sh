#!/usr/bin/env bash
# Usage:
#   ./add_hb_column_to_dist.sh input.dist input.hb output_with_hb.dist
#
# Produces a dist-like file where each block has an extra line containing
# hb_count (number of HB pairs for that step) inserted so that make_obs.sh
# can read it as the "hb number" line.

set -euo pipefail

DIST="${1:?need input.dist}"
HB="${2:?need input.hb}"
OUT="${3:?need output_with_hb.dist}"

awk -v HBFILE="$HB" '
# ---------- First pass: read HB counts into hb_count[step] ----------
FNR==NR {
  # header like: "# step 1000"
  if ($0 ~ /^[[:space:]]*#[[:space:]]*step[[:space:]]+[0-9]+/) {
    t=$0
    sub(/^[[:space:]]*#[[:space:]]*step[[:space:]]+/, "", t)
    split(t, a, /[^0-9]+/)
    cur = a[1] + 0
    hb_count[cur] = 0
    next
  }
  # count bond rows like "2 98"
  if (cur != "" && $0 !~ /^[[:space:]]*#/ && NF>=2 && $1 ~ /^-?[0-9]+$/ && $2 ~ /^-?[0-9]+$/) {
    hb_count[cur]++
  }
  next
}

# ---------- Second pass: stream the DIST file and inject hb_count ----------
function is_stepdist_line() { return ($0 !~ /^[[:space:]]*#/ && NF>=2 && $1 ~ /^-?[0-9]+$/) }
function is_comment_line()  { return ($0 ~ /^[[:space:]]*#/) }
function is_blank_line()    { return ($0 ~ /^[[:space:]]*$/) }
function is_3float_line()   { return (NF==3) }  # good enough for your dist structure

BEGIN { OFS=" " }

{
  # Copy comments/blank lines verbatim
  if (is_comment_line() || is_blank_line()) { print; next }

  # Start of a block: "step dist"
  if (is_stepdist_line()) {
    step = $1 + 0
    print $0

    # Next non-comment line: force
    while ((getline line) > 0) {
      if (line ~ /^[[:space:]]*#/) { print line; continue }
      print line
      break
    }
    if (line == "") next  # EOF safety

    # Next non-comment line: could be opposite force OR next block header (rare) OR something else
    while ((getline line2) > 0) {
      if (line2 ~ /^[[:space:]]*#/) { print line2; continue }
      # If this looks like a new block header, we should inject HB then re-process this line.
      split(line2, tmp)
      if (line2 !~ /^[[:space:]]*#/ && tmp[1] ~ /^-?[0-9]+$/ && length(tmp[2])>0 && line2 !~ /[[:alpha:]]/) {
        # inject HB then fall through by printing line2 as start of next block
        print ( (step in hb_count) ? hb_count[step] : 0 )
        print line2
        next
      }
      # Otherwise print it (this is usually opposite force)
      print line2
      break
    }

    # Now inject hb_count as a single-number line (what make_obs.sh expects)
    print ( (step in hb_count) ? hb_count[step] : 0 )
    next
  }

  # Any other non-comment line, just print
  print
}

' "$HB" "$DIST" > "$OUT"

echo "Wrote $OUT"
