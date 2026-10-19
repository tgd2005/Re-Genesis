#!/usr/bin/env python
"""
free_energy_kcal.py  —  Macrostate ΔF: bonded vs ruptured aptamer stem.
                        Values reported in kcal mol⁻¹.

Identical to free_energy.py except ΔF is converted from kT to kcal mol⁻¹
using T = 300 K (oxDNA default simulation temperature).

  1 kT  =  R·T / 1000  =  1.9872 × 300 / 1000  =  0.5961 kcal mol⁻¹

Uses v8 last_hist.dat files (which incorporate v7 counts via init_hist stitching).
The two umbrella sampling windows are joined at their shared state (bond=0, mindist_bin=0).

Method
------
  Z_bonded   = Σ_{bond=1..6}  U_HB(bond, bin=0)
  Z_ruptured = (U_HB(0,0) / U_dist(0,0))  ×  Σ_{bin=0..6}  U_dist(0, bin)
  ΔF         = −ln(Z_ruptured / Z_bonded)  ×  kT_TO_KCAL   [kcal mol⁻¹]

Usage
-----
    python free_energy_kcal.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

# ── Configuration ──────────────────────────────────────────────────────────────

ROOT = os.path.dirname(os.path.abspath(__file__))

SYSTEMS = {
    'pulling': {
        'v8_dir':      'pulling-clamps/6nt-clamps-pulling/v8',
        'hb_weight':   'wB',
        'dist_weight': 'wA',
    },
    'pushing': {
        'v8_dir':      'pushing-clamps/6nt-pushing/v8',
        'hb_weight':   'wA',
        'dist_weight': 'wB',
    },
}

WEIGHT_LABEL = '_'.join(
    f'{s}-HB-{c["hb_weight"]}_dist-{c["dist_weight"]}'
    for s, c in SYSTEMS.items()
)

COLORS = {'pulling': '#E05A4E', 'pushing': '#4E7BE0'}

N_BOOTSTRAP    = 20000
BOOTSTRAP_SEED = 42

# Conversion factor: kT → kcal mol⁻¹  at T = 300 K
R_KCAL       = 1.98720425864e-3   # kcal mol⁻¹ K⁻¹
T_SIM        = 293.0              # K  (20 °C, experimental temperature)
KT_TO_KCAL   = R_KCAL * T_SIM    # = 0.59616 kcal mol⁻¹

# ── Style ──────────────────────────────────────────────────────────────────────

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
    "font.family":       "sans-serif",
})

# ── I/O ────────────────────────────────────────────────────────────────────────

def read_last_hist(path):
    """
    Parse last_hist.dat.
    Returns {(bond, mindist_bin): unbiased_count}  (column 4 of each data row).
    """
    counts = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            counts[(int(parts[0]), int(parts[1]))] = float(parts[3])
    return counts


def find_seeds(base_dir, weight_prefix, is_dist=False):
    """
    Return sorted list of (seed_label, last_hist_path).
    HB   dirs:  wB_s1, wB_s2, …
    Dist dirs:  dist_wA_s1, dist_wA_s2, …
    """
    prefix = f'dist_{weight_prefix}_s' if is_dist else f'{weight_prefix}_s'
    seeds = []
    for entry in sorted(os.listdir(base_dir)):
        if entry.startswith(prefix):
            hp = os.path.join(base_dir, entry, 'last_hist.dat')
            if os.path.isfile(hp):
                seeds.append((entry, hp))
    return seeds

# ── Core computation ───────────────────────────────────────────────────────────

def compute_dF_per_seed(hb_seeds, dist_seeds):
    """
    Pair HB and dist seeds by index, compute ΔF for each pair.
    Returns list of result dicts with dF in both kT and kcal mol⁻¹.
    Skips pairs with zero counts.
    """
    results = []
    for (hb_name, hb_path), (dist_name, dist_path) in zip(hb_seeds, dist_seeds):
        hb   = read_last_hist(hb_path)
        dist = read_last_hist(dist_path)

        z_bonded  = sum(hb.get((b, 0), 0.0) for b in range(1, 7))
        u_hb_00   = hb.get((0, 0), 0.0)
        u_dist_00 = dist.get((0, 0), 0.0)

        if z_bonded == 0 or u_hb_00 == 0 or u_dist_00 == 0:
            print(f'  SKIP {hb_name}+{dist_name}: '
                  f'z_bonded={z_bonded:.3g}  u_hb(0,0)={u_hb_00:.3g}  '
                  f'u_dist(0,0)={u_dist_00:.3g}')
            continue

        scale      = u_hb_00 / u_dist_00
        z_ruptured = scale * sum(dist.get((0, b), 0.0) for b in range(0, 7))
        dF_kT      = -np.log(z_ruptured / z_bonded)
        dF_kcal    = dF_kT * KT_TO_KCAL

        results.append({
            'hb_seed':    hb_name,
            'dist_seed':  dist_name,
            'z_bonded':   z_bonded,
            'z_ruptured': z_ruptured,
            'scale':      scale,
            'dF_kT':      dF_kT,
            'dF':         dF_kcal,
            '_dist_bins': [dist.get((0, b), 0.0) for b in range(7)],
        })
        print(f'  {hb_name} + {dist_name}:  '
              f'Z_bonded={z_bonded:.3e}  Z_ruptured={z_ruptured:.3e}  '
              f'scale={scale:.3e}  ΔF={dF_kT:+.4f} kT  ({dF_kcal:+.4f} kcal/mol)')

    return results

# ── Jitter helpers ─────────────────────────────────────────────────────────────

def _auto_bin_points(y, bins='auto', min_bins=3, max_bins=12):
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n == 0:
        return np.array([], dtype=int), np.array([0.0, 1.0]), {}
    if n == 1:
        return np.array([0], dtype=int), np.array([y[0]-0.5, y[0]+0.5]), {0: np.array([0])}

    edges  = np.histogram_bin_edges(y, bins=bins)
    n_bins = len(edges) - 1
    if n_bins < min_bins or n_bins > max_bins:
        n_bins = int(np.clip(n_bins, min_bins, max_bins))
        edges  = np.linspace(y.min(), y.max(), n_bins + 1)

    bin_idx        = np.digitize(y, edges[1:-1], right=False)
    bin_to_indices = {b: np.where(bin_idx == b)[0] for b in range(n_bins)}
    return bin_idx, edges, bin_to_indices


def _jitter_points(y, bin_to_indices, jitter_range=0.25):
    y         = np.asarray(y, dtype=float)
    x_offsets = np.zeros(len(y), dtype=float)
    if len(y) == 0:
        return x_offsets

    largest_len = max((len(idx) for idx in bin_to_indices.values()), default=1)

    for b, idx in bin_to_indices.items():
        count = len(idx)
        if count == 0:
            continue
        if count == 1:
            x_offsets[idx[0]] = 0.0
            continue

        cur_odd = count % 2 == 1
        lrg_odd = largest_len % 2 == 1
        if lrg_odd and cur_odd:
            denom, cand = largest_len, largest_len
        elif lrg_odd and not cur_odd:
            denom, cand = largest_len + 1, largest_len + 1
        elif not lrg_odd and not cur_odd:
            denom, cand = largest_len, largest_len
        else:
            denom, cand = largest_len + 1, largest_len + 1

        step       = jitter_range / max(denom, 1)
        candidates = (np.arange(cand) - (cand - 1) / 2.0) * step
        chosen     = candidates[np.argsort(np.abs(candidates))[:count]]
        chosen     = np.sort(chosen)
        idx_sorted = idx[np.argsort(y[idx])]
        x_offsets[idx_sorted] = chosen

    return x_offsets

# ── Plots ──────────────────────────────────────────────────────────────────────

def plot_dF(all_results, out_dir):
    """
    Strip plot with smart jitter, mean ± SD per system.  ΔF in kcal mol⁻¹.
    """
    group_order  = ['pushing', 'pulling']
    group_labels = {'pulling': 'Pull', 'pushing': 'Push'}
    group_labels = [group_labels[s] for s in group_order]
    group_values = {s: np.array([r['dF'] for r in all_results[s]]) for s in group_order}

    fig, ax0 = plt.subplots(1, 1, figsize=(3.0, 1.5))

    for xi, sys_name in enumerate(group_order):
        values = group_values[sys_name]
        color  = COLORS[sys_name]

        _, _, bin_to_indices = _auto_bin_points(values, bins='fd', min_bins=3, max_bins=10)
        x_offsets = _jitter_points(values, bin_to_indices, jitter_range=0.40)

        ax0.scatter(xi + x_offsets, values, color=color, alpha=0.7, s=15, zorder=3)

        mean = values.mean() if len(values) > 0 else np.nan
        std  = (values.std(ddof=2) if len(values) > 2
                else (values.std(ddof=1) if len(values) > 1 else 0.0))

        ax0.errorbar(xi, mean, yerr=std, fmt='none', color='black',
                     capsize=8, elinewidth=1, alpha=0.9, zorder=4)
        ax0.hlines(mean, xi - 0.25, xi + 0.25, color='black', linewidth=1, zorder=5)

    ax0.set_xlim(-0.5, len(group_order) - 0.5)
    ax0.set_xticks(range(len(group_order)))
    ax0.set_xticklabels(group_labels)
    ax0.set_ylabel(r'$\Delta F$  (kcal mol$^{-1}$)')
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)

    fig.subplots_adjust(left=0.20, right=0.95, bottom=0.15, top=0.95)
    for ext in ('pdf', 'svg'):
        path = os.path.join(out_dir, f'dF_macrostate_kcal.{ext}')
        fig.savefig(path, bbox_inches='tight', dpi=150)
        print(f'  Saved: {path}')
    plt.close(fig)

    summary = {
        s: (group_values[s].mean(), group_values[s].std(ddof=1), len(group_values[s]))
        for s in group_order if len(group_values[s]) > 0
    }
    return summary


def plot_partition_functions(all_results, out_dir):
    """Log-scale view of Z_bonded and Z_ruptured per seed — sanity check."""
    n_sys = sum(1 for r in all_results.values() if r)
    if n_sys == 0:
        return

    fig, axes = plt.subplots(1, n_sys, figsize=(3.0 * n_sys, 3.2), sharey=True)
    if n_sys == 1:
        axes = [axes]
    np.random.seed(42)

    all_z = [v for res in all_results.values()
              for r in res for v in (r['z_bonded'], r['z_ruptured']) if v > 0]
    if all_z:
        shared_ylim = (min(all_z) * 0.3, max(all_z) * 3.0)

    ax_idx = 0
    for sys_name, results in all_results.items():
        if not results:
            continue
        ax  = axes[ax_idx]; ax_idx += 1
        zb  = np.array([r['z_bonded']   for r in results])
        zr  = np.array([r['z_ruptured'] for r in results])
        xs  = np.arange(len(zb))
        jit = np.random.uniform(-0.1, 0.1, len(zb))

        ax.scatter(xs + jit - 0.15, zb, color='steelblue', s=18, alpha=0.75,
                   label='Z_bonded',   zorder=3)
        ax.scatter(xs + jit + 0.15, zr, color='firebrick', s=18, alpha=0.75,
                   label='Z_ruptured', zorder=3)

        ax.set_yscale('log')
        if all_z:
            ax.set_ylim(shared_ylim)
        ax.set_xticks(xs)
        ax.set_xticklabels([r['hb_seed'].split('_s')[-1] for r in results], fontsize=7)
        ax.set_xlabel('Seed index')
        ax.set_title(sys_name.capitalize(), fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if ax_idx == 1:
            ax.set_ylabel('Unbiased count  (log scale)')
            ax.legend(fontsize=8, frameon=False)

    fig.suptitle('Partition function components per seed  '
                 '(Z_ruptured scaled to HB window at overlap state)',
                 fontsize=8, y=1.02)
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.15, top=0.88, wspace=0.08)
    for ext in ('pdf', 'svg'):
        path = os.path.join(out_dir, f'partition_functions.{ext}')
        fig.savefig(path, bbox_inches='tight', dpi=150)
        print(f'  Saved: {path}')
    plt.close(fig)


def plot_dist_bins(all_results, out_dir):
    """Stacked bar showing scaled unbiased count per dist bin (bond=0) per seed."""
    n_sys = sum(1 for r in all_results.values() if r)
    if n_sys == 0:
        return

    bin_colors  = plt.cm.viridis_r(np.linspace(0.1, 0.9, 7))
    bin_labels  = [f'bin {b}' for b in range(7)]

    fig, axes = plt.subplots(1, n_sys, figsize=(3.5 * n_sys, 3.2), sharey=False)
    if n_sys == 1:
        axes = [axes]

    ax_idx = 0
    for sys_name, results in all_results.items():
        if not results:
            continue
        ax     = axes[ax_idx]; ax_idx += 1
        xs     = np.arange(len(results))
        bottom = np.zeros(len(results))

        for b in range(7):
            heights = np.array([r['scale'] * r['_dist_bins'][b] for r in results])
            ax.bar(xs, heights, bottom=bottom, color=bin_colors[b],
                   label=bin_labels[b], width=0.6)
            bottom += heights

        ax.set_xticks(xs)
        ax.set_xticklabels([r['hb_seed'].split('_s')[-1] for r in results], fontsize=7)
        ax.set_xlabel('Seed index')
        ax.set_title(sys_name.capitalize(), fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if ax_idx == 1:
            ax.set_ylabel('Scaled unbiased count')

    fig.legend(
        handles=[axes[0].containers[b] for b in range(7)],
        labels=bin_labels,
        loc='upper center',
        ncol=7,
        fontsize=7,
        frameon=False,
        title='mindist bin',
        title_fontsize=7,
        bbox_to_anchor=(0.5, 1.08),
    )
    fig.suptitle('Z_ruptured breakdown by mindist bin  (bond=0)',
                 fontsize=8, y=1.17)
    fig.subplots_adjust(top=0.82)
    for ext in ('pdf', 'svg'):
        path = os.path.join(out_dir, f'dist_bin_breakdown.{ext}')
        fig.savefig(path, bbox_inches='tight', dpi=150)
        print(f'  Saved: {path}')
    plt.close(fig)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir   = os.path.join(ROOT, 'free_energy_results', f'{timestamp}_{WEIGHT_LABEL}_kcal')
    os.makedirs(out_dir, exist_ok=True)
    print(f'Output directory: {out_dir}')
    print(f'Conversion: 1 kT = {KT_TO_KCAL:.5f} kcal/mol  (T = {T_SIM} K)\n')

    all_results = {}

    for sys_name, cfg in SYSTEMS.items():
        base = os.path.join(ROOT, cfg['v8_dir'])
        print(f'── {sys_name.upper()} ───────────────────────────────────────')

        hb_seeds   = find_seeds(base, cfg['hb_weight'],  is_dist=False)
        dist_seeds = find_seeds(base, cfg['dist_weight'], is_dist=True)
        print(f'  HB seeds   ({cfg["hb_weight"]}): {len(hb_seeds)}')
        print(f'  Dist seeds ({cfg["dist_weight"]}): {len(dist_seeds)}')

        all_results[sys_name] = compute_dF_per_seed(hb_seeds, dist_seeds)
        print()

    print('Generating plots...')
    summary = plot_dF(all_results, out_dir)
    plot_partition_functions(all_results, out_dir)
    plot_dist_bins(all_results, out_dir)

    print('\n── Summary ─────────────────────────────────────────────')
    print(f'{"System":<12} {"N":>4}  {"mean ΔF (kcal/mol)":>20}  {"SD":>9}')
    print('─' * 52)
    for sys_name, (mean, sd, n) in summary.items():
        print(f'{sys_name:<12} {n:>4}  {mean:>+20.4f}  {sd:>9.4f}')

    csv_path = os.path.join(out_dir, 'dF_results_kcal.csv')
    with open(csv_path, 'w') as fh:
        fh.write('system,hb_seed,dist_seed,z_bonded,z_ruptured,scale,dF_kT,dF_kcal_per_mol\n')
        for sys_name, results in all_results.items():
            for r in results:
                fh.write(
                    f'{sys_name},{r["hb_seed"]},{r["dist_seed"]},'
                    f'{r["z_bonded"]:.6e},{r["z_ruptured"]:.6e},'
                    f'{r["scale"]:.6e},{r["dF_kT"]:.6f},{r["dF"]:.6f}\n'
                )
    print(f'  Saved: {csv_path}')


if __name__ == '__main__':
    main()
