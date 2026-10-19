#!/usr/bin/env python3
"""
Convert an "old-style" oxDNA topology where lines are in 5'->3' order per strand,
but still use the 4-column format:

  strand_id  base  n3  n5   (mate may have mixed semantics, but we IGNORE n3/n5)

Output (new format):
  <N> <n_strands> 5->3
  <strand1_sequence>
  <strand2_sequence>
  ...

We assume the nucleotide lines are already ordered 5'->3' for each strand.
"""

import sys
from pathlib import Path

def convert(in_path, out_path):
    # read and keep only non-empty non-comment lines
    raw_lines = in_path.read_text().splitlines()
    lines = []
    for ln in raw_lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        lines.append(s)

    if not lines:
        raise ValueError("Empty input file.")

    header = lines[0].split()
    if len(header) < 2:
        raise ValueError("Bad header line: %r" % lines[0])

    N = int(header[0])
    # We don't trust the old strand count if their file is weird; we'll compute it.
    seq_by_strand = {}
    order = []  # first-seen order of strand IDs

    for i, ln in enumerate(lines[1:]):
        parts = ln.split()
        if len(parts) < 2:
            raise ValueError("Bad nucleotide line at %d: %r" % (i, ln))
        sid = int(parts[0])
        b = parts[1]

        if sid not in seq_by_strand:
            seq_by_strand[sid] = []
            order.append(sid)
        seq_by_strand[sid].append(b)

    if sum(len(v) for v in seq_by_strand.values()) != N:
        raise ValueError("Header says N=%d, but parsed %d nucleotide lines."
                         % (N, sum(len(v) for v in seq_by_strand.values())))

    out_lines = ["%d %d 5->3" % (N, len(order))]
    for sid in order:
        out_lines.append("".join(seq_by_strand[sid]))

    out_path.write_text("\n".join(out_lines) + "\n")

def main(argv):
    if len(argv) != 3:
        sys.stderr.write("Usage: %s old.top new.top\n" % argv[0])
        return 2
    in_path = Path(argv[1])
    out_path = Path(argv[2])
    convert(in_path, out_path)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
