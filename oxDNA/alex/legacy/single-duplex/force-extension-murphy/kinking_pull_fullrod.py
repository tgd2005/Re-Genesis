#!/usr/bin/env python3
"""
Kinking frequency vs extension for pull simulations (full rod, kinking seeds).
Kinking = coaxial unstacking at EITHER nick (pair 0-85 OR pair 31-32).
"""

import os
import re
import statistics
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── constants ─────────────────────────────────────────────────────────────────
DIST_CONVERSION = 0.8       # oxDNA units → nm
CUTOFF_BLOCKS   = 250_000   # skip first 250k timestep blocks (= step 2.5e7)
ZERO_TOL        = 1e-9      # |CXSTCK| < this → coaxially unstacked (kinked)

INPUT_FOLDER = "cxstack_pair_energy_dats"
PATTERN = re.compile(r"pair_energy_K_025_r0_(\d+)_seed_(\d+)\.dat$")

# ── matplotlib style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    "axes.linewidth":    1,
    "xtick.major.width": 1,
    "ytick.major.width": 1,
    "xtick.major.size":  4,
    "ytick.major.size":  4,
    "font.size":         9,
    "axes.labelsize":    9,
    "xtick.minor.width": 1,
    "ytick.minor.width": 1,
    "xtick.minor.size":  2,
    "ytick.minor.size":  2,
    "font.sans-serif":   ["Arial"],
})


# ── fast parser: direct Python I/O, only steady-state half ───────────────────
def parse_cxstck(filepath):
    """
    Read CXSTCK (field index 8) from non-comment lines, steady-state only.
    Each timestep block = 1 comment line + 2 data lines (pair1 then pair2).
    Skip the first CUTOFF_BLOCKS blocks, then collect CXSTCK for both pairs.
    Returns (pair1_vals, pair2_vals) as lists.
    """
    cxstck = []
    blocks_seen = 0
    collecting = False

    with open(filepath) as f:
        for line in f:
            if line[0] == '#':
                blocks_seen += 1
                if blocks_seen > CUTOFF_BLOCKS:
                    collecting = True
                continue
            if not collecting:
                continue
            # split only as far as needed (field 8, 0-indexed → 9th token)
            parts = line.split(None, 9)
            if len(parts) >= 9:
                try:
                    cxstck.append(float(parts[8]))
                except ValueError:
                    pass

    return cxstck[0::2], cxstck[1::2]   # pair1, pair2


def kinking_freq(pair1, pair2):
    """
    Fraction of timesteps where EITHER nick is unstacked (CXSTCK ≈ 0).
    max(p1, p2) ≈ 0  iff  at least one pair is at zero.
    """
    n = min(len(pair1), len(pair2))
    kinked = sum(
        1 for p1, p2 in zip(pair1[:n], pair2[:n])
        if abs(max(p1, p2)) < ZERO_TOL
    )
    return kinked / n if n > 0 else float("nan")


# ── main ──────────────────────────────────────────────────────────────────────
seed_results = defaultdict(list)   # r0 → [freq per seed]

files = sorted(f for f in os.listdir(INPUT_FOLDER) if PATTERN.match(f))
total = len(files)
print(f"Processing {total} files (K=025 only, sequential)...", flush=True)

for i, fname in enumerate(files, 1):
    m = PATTERN.match(fname)
    r0, seed = int(m.group(1)), int(m.group(2))
    print(f"[{i:3d}/{total}]  r0={r0:2d}  seed={seed}", end="  ", flush=True)

    p1, p2 = parse_cxstck(os.path.join(INPUT_FOLDER, fname))
    freq = kinking_freq(p1, p2)
    seed_results[r0].append(freq)
    print(f"kink={freq:.3f}", flush=True)

# ── aggregate: mean across seeds ─────────────────────────────────────────────
r0_sorted = sorted(seed_results)
x_nm      = [DIST_CONVERSION * r0 for r0 in r0_sorted]
y_mean    = [statistics.mean(seed_results[r0]) for r0 in r0_sorted]

# ── plot ──────────────────────────────────────────────────────────────────────
# scatter all seed points
x_scatter = [DIST_CONVERSION * r0 for r0 in r0_sorted for _ in seed_results[r0]]
y_scatter = [f for r0 in r0_sorted for f in seed_results[r0]]

plt.figure(figsize=(2, 1,7))
plt.scatter(x_scatter, y_scatter, alpha=0.7, s=8, marker='o', color="#4a85cc")
plt.plot(x_nm, y_mean, color="#2308bd", lw=1.2)
plt.xlabel("Extension (nm)")
plt.ylabel("Kinking Frequency")
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)
plt.xlim(0, 26)
plt.ylim(0, 1.05)
plt.tight_layout()
plt.savefig("kinking_pull_fullrod.pdf")
plt.savefig("kinking_pull_fullrod.svg")
print("\nSaved kinking_pull_fullrod.pdf")
