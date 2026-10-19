#!/usr/bin/env python3
"""
Kinking frequency vs mean extension for pull simulations (full rod, kinking seeds).
Mean extension computed from obs_dist.dat (col 2, oxDNA units * 0.8 -> nm).
Kinking from obs_cxstck.dat: unstacked at EITHER nick.
"""

import os
import statistics
from collections import defaultdict
from scipy import stats
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── constants ─────────────────────────────────────────────────────────────────
DIST_CONVERSION = 0.8        # oxDNA units → nm
CUTOFF_STEP     = 2.5e7      # steady-state cutoff
ZERO_TOL        = 1e-9       # |CXSTCK| < this → coaxially unstacked (kinked)

ANALYSIS_DIR = "fd_array_fullrod_kinking_seeds/analysis_files"
K_TAG = "025"
SEEDS = {1, 2, 3, 4, 5}

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


def parse_mean_ext(filepath):
    """Mean distance (nm) from obs_dist.dat after steady-state cutoff."""
    vals = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                step = float(parts[0])
                dist = float(parts[1])
            except ValueError:
                continue
            if step > CUTOFF_STEP:
                vals.append(dist * DIST_CONVERSION)
    return statistics.mean(vals) if vals else float("nan")


def parse_kinking_freq(filepath):
    """Fraction of steady-state timesteps where EITHER nick is unstacked."""
    kinked = 0
    total  = 0
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
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
            if step > CUTOFF_STEP:
                total += 1
                if abs(max(p1, p2)) < ZERO_TOL:
                    kinked += 1
    return kinked / total if total > 0 else float("nan")


# ── collect data ──────────────────────────────────────────────────────────────
# per_r0[r0] = list of (mean_ext_nm, kink_freq) tuples, one per seed
per_r0 = defaultdict(list)

seed_dirs = sorted(
    d for d in os.listdir(ANALYSIS_DIR)
    if d.startswith("Seed_") and os.path.isdir(os.path.join(ANALYSIS_DIR, d))
    and int(d.split("_")[1]) in SEEDS
)
total_seeds = len(seed_dirs)
print(f"Found {total_seeds} seed directories", flush=True)

for seed_dir in seed_dirs:
    k_dir = os.path.join(ANALYSIS_DIR, seed_dir, f"K_{K_TAG}")
    if not os.path.isdir(k_dir):
        continue

    r0_dirs = sorted(
        (d for d in os.listdir(k_dir)
         if d.startswith("r0_") and os.path.isdir(os.path.join(k_dir, d))),
        key=lambda d: int(d.split("_")[1])
    )

    for r0_dir in r0_dirs:
        r0 = int(r0_dir.split("_")[1])
        base = os.path.join(k_dir, r0_dir)

        dist_file   = os.path.join(base, "obs_dist.dat")
        cxstck_file = os.path.join(base, "obs_cxstck.dat")

        if not os.path.isfile(dist_file) or not os.path.isfile(cxstck_file):
            print(f"  MISSING files for {seed_dir} r0={r0}", flush=True)
            continue

        mean_ext  = parse_mean_ext(dist_file)
        kink_freq = parse_kinking_freq(cxstck_file)
        per_r0[r0].append((mean_ext, kink_freq))
        print(f"  {seed_dir}  r0={r0:2d}  ext={mean_ext:.3f} nm  kink={kink_freq:.3f}", flush=True)

# ── aggregate ─────────────────────────────────────────────────────────────────
x_scatter, y_scatter = [], []

for r0 in sorted(per_r0):
    for ext, kink in per_r0[r0]:
        x_scatter.append(ext)
        y_scatter.append(kink)

x_arr = np.array(x_scatter)
y_arr = np.array(y_scatter)

# linear regression on all scatter points
result = stats.linregress(x_arr, y_arr)
slope, intercept, r, p, se_slope = result.slope, result.intercept, result.rvalue, result.pvalue, result.stderr
se_intercept = result.intercept_stderr
r2 = r**2

print(f"\nLinear regression:")
print(f"  slope     = {slope:.6f} ± {se_slope:.6f} nm⁻¹")
print(f"  intercept = {intercept:.6f} ± {se_intercept:.6f}")
print(f"  R²        = {r2:.4f}")
print(f"  p-value   = {p:.3e}")

# zero-crossing of pull F-E curve for 54nt full rod (interpolated from FD data)
XMAX = 25.3

x_line = np.array([0, XMAX])
y_line = slope * x_line + intercept

FIG_W, FIG_H = 2.25, 1.7
SIZE_TAG = f"{FIG_W}x{str(FIG_H).replace('.', 'p')}"


def make_axes():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, 1.00)
    tick_step = 5
    ax.set_xticks(np.arange(0, XMAX + tick_step, tick_step))
    ax.set_yticks([0.00, 0.25, 0.50, 0.75, 1.00])
    ax.set_xlabel("Mean extension (nm)")
    ax.set_ylabel("Kinking frequency")
    return fig, ax


def save_fig(fig, stem):
    for ext in ("pdf", "svg"):
        path = f"{stem}_{SIZE_TAG}.{ext}"
        fig.savefig(path)
        print(f"Saved {path}")
    plt.close(fig)


# ── plot: scatter only ────────────────────────────────────────────────────────
fig, ax = make_axes()
ax.scatter(x_scatter, y_scatter, alpha=0.7, s=8, marker='o', color="#4a85cc")
plt.tight_layout()
save_fig(fig, "kinking_pull_fullrod_vs_ext_scatter")

# ── plot: scatter + regression ────────────────────────────────────────────────
fig, ax = make_axes()
ax.scatter(x_scatter, y_scatter, alpha=0.7, s=8, marker='o', color="#4a85cc")
ax.plot(x_line, y_line, color="#2308bd", lw=1.5, zorder=5)
plt.tight_layout()
save_fig(fig, "kinking_pull_fullrod_vs_ext_reg")
