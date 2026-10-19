import argparse
from pathlib import Path

def parse_k_list(s: str):
    # Accept either "0.25,0.35,0.45" or "0.25 0.35 0.45"
    parts = [p.strip() for p in s.replace(",", " ").split() if p.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("k list is empty")
    try:
        return [float(p) for p in parts]
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Could not parse k values: {e}")

def main():
    ap = argparse.ArgumentParser(
        description="Generate params.tsv with columns: #idx, k, r0, seed"
    )
    ap.add_argument(
        "--k",
        required=True,
        type=parse_k_list,
        help='List of k values, e.g. "0.25,0.35,0.45" (commas or spaces ok)',
    )
    ap.add_argument(
        "--r0-min",
        type=float,
        required=True,
        help="Minimum r0 (inclusive), e.g. 2",
    )
    ap.add_argument(
        "--r0-max",
        type=float,
        required=True,
        help="Maximum r0 (inclusive), e.g. 31",
    )
    ap.add_argument(
        "--r0-step",
        type=float,
        default=1.0,
        help="Step for r0, default 1.0",
    )
    ap.add_argument(
        "--seeds",
        type=int,
        required=True,
        help="Number of seeds, e.g. 10 (will generate seeds 1..N)",
    )
    ap.add_argument(
        "-o",
        "--out",
        default="params.tsv",
        help='Output file path (default: "params.tsv")',
    )
    ap.add_argument(
        "--no-header",
        action="store_true",
        help="If set, do not write the header line",
    )

    args = ap.parse_args()

    if args.r0_step <= 0:
        raise SystemExit("ERROR: --r0-step must be > 0")
    if args.r0_max < args.r0_min:
        raise SystemExit("ERROR: --r0-max must be >= --r0-min")
    if args.seeds <= 0:
        raise SystemExit("ERROR: --seeds must be >= 1")

    # Build r0 values inclusive, robust to float steps
    r0_values = []
    x = args.r0_min
    # epsilon to avoid missing the endpoint due to floating point
    eps = abs(args.r0_step) * 1e-9 + 1e-12
    while x <= args.r0_max + eps:
        # Clean representation (keeps integers as "2.0" if you pass floats; that's fine)
        r0_values.append(x)
        x += args.r0_step

    out_path = Path(args.out)
    idx = 1

    with out_path.open("w", encoding="utf-8") as f:
        if not args.no_header:
            f.write("#idx\tk\tr0\tseed\n")

        # Order: for seed 1 -> all k x r0 combos, then seed 2, ...
        for seed in range(1, args.seeds + 1):
            for k in args.k:
                for r0 in r0_values:
                    f.write(f"{idx}\t{k}\t{r0}\t{seed}\n")
                    idx += 1

if __name__ == "__main__":
    main()
