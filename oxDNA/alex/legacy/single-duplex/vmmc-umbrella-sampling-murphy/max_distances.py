#!/usr/bin/env python3
"""
max_distances.py — scan all dist_obs.dat files across v7/v8/v9/v10
and report the maximum mindistance reached per system, version, combo, and seed.
"""
import numpy as np
import glob, os, re

ROOT = os.path.dirname(os.path.abspath(__file__))

SYS_DIRS = {
    'pulling':  'pulling-clamps/6nt-clamps-pulling',
    'noclamps': 'noclamps/6nt-noclamps',
    'pushing':  'pushing-clamps/6nt-pushing',
}

VERSIONS = ['v7', 'v8', 'v9', 'v10']

results = []  # (system, version, combo, seed, n_points, max_mindist, p95)

for system, sdir in SYS_DIRS.items():
    for version in VERSIONS:
        pattern = os.path.join(ROOT, sdir, version, 'dist_*', 'dist_obs.dat')
        for fpath in sorted(glob.glob(pattern)):
            run = os.path.basename(os.path.dirname(fpath))  # e.g. dist_wA_s1
            m = re.match(r'dist_(w\w+)_(s\d+)', run)
            combo, seed = (m.group(1), m.group(2)) if m else (run, '?')
            try:
                data = np.loadtxt(fpath)
                if data.ndim < 2 or data.shape[1] < 2:
                    continue
                mindist = data[:, 1:].min(axis=1)
                results.append((system, version, combo, seed,
                                 len(mindist), mindist.max(),
                                 np.percentile(mindist, 95)))
            except Exception as e:
                print(f"  SKIP {fpath}: {e}")

# ── Per-system summary ────────────────────────────────────────────────────────
print("\n" + "="*70)
print("MAX MINDISTANCE SUMMARY  (oxDNA units;  1 unit ≈ 0.85 nm)")
print("="*70)

for system in ['pulling', 'noclamps', 'pushing']:
    rows = [r for r in results if r[0] == system]
    if not rows:
        print(f"\n{system.upper()}  — no data"); continue

    overall_max = max(r[5] for r in rows)
    best = max(rows, key=lambda r: r[5])

    print(f"\n{system.upper()}")
    print(f"  Overall max across all versions/seeds: {overall_max:.2f}  "
          f"(best: {best[1]}/{best[2]}/{best[3]})")

    # per-version breakdown
    for v in VERSIONS:
        vrows = [r for r in rows if r[1] == v]
        if not vrows:
            continue
        vmax  = max(r[5] for r in vrows)
        vp95  = max(r[6] for r in vrows)
        vbest = max(vrows, key=lambda r: r[5])
        print(f"  {v}  n_seeds={len(vrows):3d}  max={vmax:6.2f}  "
              f"best_p95={vp95:6.2f}  best_seed={vbest[2]}/{vbest[3]}")

# ── Full per-seed table ───────────────────────────────────────────────────────
print("\n" + "="*70)
print("PER-SEED TABLE")
print(f"{'system':<10} {'ver':<4} {'combo':<6} {'seed':<5} {'n_pts':>7} "
      f"{'max':>7} {'p95':>7}")
print("-"*70)
for r in sorted(results, key=lambda x: (x[0], x[1], x[2], x[3])):
    system, ver, combo, seed, n, mx, p95 = r
    print(f"{system:<10} {ver:<4} {combo:<6} {seed:<5} {n:>7d} {mx:>7.2f} {p95:>7.2f}")

