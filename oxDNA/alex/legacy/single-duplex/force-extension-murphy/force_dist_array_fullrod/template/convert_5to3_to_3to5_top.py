import argparse
import sys
import re

VALID_BASES = set("ACGT")

def parse_input(path: str):
    with open(path, "r") as f:
        raw_lines = [ln.strip() for ln in f if ln.strip()]

    if not raw_lines:
        raise ValueError("Input file is empty.")

    # Header: e.g. "170 4 5->3" or "170 4"
    header_parts = raw_lines[0].split()
    if len(header_parts) < 2:
        raise ValueError("Header must have at least two fields: <N> <num_strands> [optional direction].")

    try:
        N = int(header_parts[0])
        S = int(header_parts[1])
    except ValueError as e:
        raise ValueError("First two header fields must be integers: <N> <num_strands>.") from e

    seq_lines = raw_lines[1:]
    if len(seq_lines) < S:
        raise ValueError(f"Expected {S} strand sequence lines after the header, found {len(seq_lines)}.")

    # Take exactly S strands (common format: one line per strand)
    seqs = seq_lines[:S]

    # Clean/validate sequences
    cleaned = []
    for i, s in enumerate(seqs, start=1):
        s = re.sub(r"\s+", "", s).upper()
        if not s:
            raise ValueError(f"Strand {i} is empty.")
        bad = set(s) - VALID_BASES
        if bad:
            raise ValueError(f"Strand {i} contains invalid characters {sorted(bad)}; allowed: A,C,G,T.")
        cleaned.append(s)

    total = sum(len(s) for s in cleaned)
    if total != N:
        raise ValueError(f"Header N={N}, but total bases across {S} strands is {total}.")

    return N, S, cleaned

def write_old_topology(out_path: str, N: int, S: int, seqs_5to3):
    # Convert each strand to 3'->5' by reversing the sequence
    seqs_3to5 = [s[::-1] for s in seqs_5to3]

    with open(out_path, "w") as out:
        out.write(f"{N} {S}\n")
        idx = 0  # global nucleotide index (0-based)
        for strand_id, seq in enumerate(seqs_3to5, start=1):
            L = len(seq)
            for j, base in enumerate(seq):
                prev_idx = -1 if j == 0 else idx - 1
                next_idx = -1 if j == L - 1 else idx + 1
                out.write(f"{strand_id} {base} {prev_idx} {next_idx}\n")
                idx += 1

def main():
    ap = argparse.ArgumentParser(
        description="Convert oxDNA new (5'->3' sequences per strand) topology to old explicit-neighbor (3'->5') topology."
    )
    ap.add_argument("input", help="Input topology file (new 5'->3' format)")
    ap.add_argument("output", help="Output topology file (old explicit neighbor format)")
    args = ap.parse_args()

    try:
        N, S, seqs = parse_input(args.input)
        write_old_topology(args.output, N, S, seqs)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
