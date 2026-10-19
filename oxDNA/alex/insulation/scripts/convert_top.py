#!/usr/bin/env python3
"""
Convert old oxDNA topology (per-nucleotide lines: strand base n3 n5)
into newer "5->3" sequence-per-strand .top format:

<N> <n_strands> 5->3
<seq strand 1 5->3>
<seq strand 2 5->3>
...

Assumptions (standard oxDNA meaning):
- n3 = neighbor index in the 3' direction (or -1 at 3' end)
- n5 = neighbor index in the 5' direction (or -1 at 5' end)
- nucleotide indices are 0-based, in file order (first nucleotide line is index 0)
"""

import sys
from pathlib import Path

def parse_old_topology(text):
    lines = []
    for ln in text.splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        lines.append(s)

    if not lines:
        raise ValueError("Empty topology text.")

    header = lines[0].split()
    if len(header) < 2:
        raise ValueError("Bad header line: %r" % lines[0])

    N = int(header[0])
    n_strands = int(header[1])

    strand_id = []
    base = []
    n3 = []
    n5 = []

    for i, ln in enumerate(lines[1:]):
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError("Bad nucleotide line (need 4 cols) at nucleotide %d: %r" % (i, ln))
        strand_id.append(int(parts[0]))
        base.append(parts[1])
        n3.append(int(parts[2]))
        n5.append(int(parts[3]))

    if len(base) != N:
        raise ValueError("Header says N=%d, but found %d nucleotide lines." % (N, len(base)))

    return N, n_strands, strand_id, base, n3, n5

def build_sequences_5to3(N, strand_id, base, n3, n5):
    # indices per strand
    by_strand = {}
    for i, s in enumerate(strand_id):
        by_strand.setdefault(s, []).append(i)

    seqs = {}

    for s, idxs in by_strand.items():
        # 5' end is where n5 == -1
        five_ends = [i for i in idxs if n5[i] == -1]
        if len(five_ends) != 1:
            raise ValueError("Strand %s: expected exactly 1 five-prime end (n5==-1), found %d."
                             % (s, len(five_ends)))
        start = five_ends[0]

        seen = set()
        out = []
        cur = start
        while cur != -1:
            if cur in seen:
                raise ValueError("Strand %s: loop detected while traversing from 5' end at %d."
                                 % (s, start))
            seen.add(cur)
            out.append(base[cur])
            cur = n3[cur]

        if seen != set(idxs):
            missing = set(idxs) - seen
            extra = seen - set(idxs)
            raise ValueError("Strand %s: traversal mismatch. Missing %d nts, extra %d nts."
                             % (s, len(missing), len(extra)))

        seqs[s] = "".join(out)

    return seqs

def convert_file(in_path, out_path):
    text = in_path.read_text()
    N, n_strands, strand_id, base, n3, n5 = parse_old_topology(text)
    seqs = build_sequences_5to3(N, strand_id, base, n3, n5)

    strand_order = sorted(seqs.keys())
    out_lines = ["%d %d 5->3" % (N, len(strand_order))]
    for s in strand_order:
        out_lines.append(seqs[s])

    out_path.write_text("\n".join(out_lines) + "\n")

def main(argv):
    if len(argv) != 3:
        sys.stderr.write("Usage: %s old.top new.top\n" % argv[0])
        return 2
    in_path = Path(argv[1])
    out_path = Path(argv[2])
    if not in_path.exists():
        sys.stderr.write("Input not found: %s\n" % in_path)
        return 2
    convert_file(in_path, out_path)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
