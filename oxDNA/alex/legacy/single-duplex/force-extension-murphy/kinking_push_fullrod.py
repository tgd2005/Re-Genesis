#!/usr/bin/env python3
"""
Kinking frequency vs extension for push simulations (full rod).
Data from Re-Genesis crxstack files: step  CXSTCK_pair1  CXSTCK_pair2
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
CUTOFF_STEP     = 2.5e7     # steady-state cutoff
ZERO_TOL        = 1e-9

INPUT_FOLDER = "/Users/murphy92/Y3_Group_Project/Re-Genesis/oxDNA/tris/cxstack_push_full_rod_dats"
PATTERN = re.compile(r"crxstack_K_(\d+)_r0_(\d+)_seed_(\d+)\.dat$")

# ── matplotlib style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    "axes.linewidth":    1.5,
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "xtick.major.size":  8,
    "ytick.major.size":  8,
    "font.size":         14,
    "axes.labelsize":    14,
})


def parse_kinking_freq(filepath):
    """
    Read step, CXSTCK_pair1, CXSTCK_pair2. Apply SS cutoff.
    Kinking = max(pair1, pair2) ≈ 0 (either nick unstacked).
    """
    kinked = 0
    total  = 0
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                step = float(parts[0])
                p1   = float(parts[1])
                p2   = float(parts[2])
            except ValueError:
                continue
            if step < CUTOFF_STEP:
                continue
            total += 1
            if abs(max(p1, p2)) < ZERO_TOL:
                kinked += 1
    return kinked / total if total > 0 else float("nan")


# ── main ──────────────────────────────────────────────────────────────────────
seed_results = defaultdict(list)   # r0 → [freq per seed]

files = sorted(f for f in os.listdir(INPUT_FOLDER) if PATTERN.match(f))
total = len(files)
print(f"Processing {total} push files...", flush=True)

for i, fname in enumerate(files, 1):
    m = PATTERN.match(fname)
    K_str, r0, seed = m.group(1), int(m.group(2)), int(m.group(3))
    freq = parse_kinking_freq(os.path.join(INPUT_FOLDER, fname))
    seed_results[r0].append(freq)
    print(f"[{i:3d}/{total}]  r0={r0:2d}  seed={seed}  kink={freq:.3f}", flush=True)

# ── aggregate ─────────────────────────────────────────────────────────────────
r0_sorted = sorted(seed_results)
x_nm   = [DIST_CONVERSION * r0 for r0 in r0_sorted]
y_mean = [statistics.mean(seed_results[r0]) for r0 in r0_sorted]

# ── plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 4))

for r0, x_val in zip(r0_sorted, x_nm):
    ax.scatter([x_val] * len(seed_results[r0]), seed_results[r0],
               color="red", alpha=0.4, s=18, zorder=2)

ax.plot(x_nm, y_mean,
        color="red", linewidth=1.5, marker="o", markersize=5,
        label="Mean", zorder=3)

ax.set_xlabel(r"Extension $r_0$ (nm)")
ax.set_ylabel("Kinking frequency")
ax.set_xlim(left=0, right=18)
ax.set_ylim(0, 1)
ax.legend(frameon=False)
plt.tight_layout()
plt.savefig("/Users/murphy92/Y3_Group_Project/oxDNA/force-extension/kinking_push_fullrod.pdf")
plt.savefig("/Users/murphy92/Y3_Group_Project/oxDNA/force-extension/kinking_push_fullrod.svg")
print("\nSaved kinking_push_fullrod.pdf")
