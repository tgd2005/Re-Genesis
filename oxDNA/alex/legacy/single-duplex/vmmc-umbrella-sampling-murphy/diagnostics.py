#!/usr/bin/env python
"""
diagnostics.py — VMMC umbrella sampling simulation quality assessment
======================================================================

Reads last_hist.dat for every run found under the given directories.
Works with any version (v4, v5, ...) and both HB and distance windows.

Usage
-----
    python3 diagnostics.py <dir1> [dir2] ...

Each argument is a directory that is scanned for run subdirectories
matching the pattern:  (dist_)?(w[ABC])_(s\d+)
System (pulling / noclamps / pushing) and version (v4 / v5 / ...)
are inferred automatically from the path string.

Example — compare v4 and v5 for all three systems:
    python3 diagnostics.py \\
        pulling-clamps/6nt-clamps-pulling/v4  pulling-clamps/6nt-clamps-pulling/v5 \\
        noclamps/6nt-noclamps/v4              noclamps/6nt-noclamps/v5 \\
        pushing-clamps/6nt-pushing/v4         pushing-clamps/6nt-pushing/v5

Outputs
-------
    diagnostics.csv          — per-seed metrics table
    diagnostics_sigma.pdf    — σ_naive strip plots (lower error bound, coloured)
    diagnostics_dg.pdf       — ΔG per seed coloured by σ_naive
"""

import os, sys, re, csv
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from collections import defaultdict

# ── Physical constants ──────────────────────────────────────────────────────
kB  = 0.001987   # kcal / (mol · K)
T_K = 293.15     # 20 °C
kT  = kB * T_K   # ~0.582 kcal / mol

# ── Run directory naming convention ────────────────────────────────────────
RUN_RE = re.compile(r'^(dist_)?(w[A-Z])_(s\d+)$')

SYSTEMS = ['pulling', 'noclamps', 'pushing']
WINDOWS = ['HB', 'dist']
COMBOS  = ['wA', 'wB', 'wC', 'wD', 'wE', 'wF', 'wG', 'wH']

# Target bin for the dist window ΔF calculation.
# Bin 6 (mindistance > 12 oxDNA units) is the true fully-separated state;
# simulations confirm the geometric max is ~16.5 units, well inside bin 6.
DIST_TARGET_BIN = 6


# ══════════════════════════════════════════════════════════════════════════════
# I/O
# ══════════════════════════════════════════════════════════════════════════════

def parse_last_hist(path):
    """
    Parse last_hist.dat → dict: (bond, mdbin) → unbiased_count.

    File format (4 data columns per row):
        bond   mdbin   biased_count   unbiased_count

    We keep column 4 (unbiased_count = biased_count / weight).
    This is accumulated at every simulation step (not just print steps),
    so it has much higher effective resolution than energy.dat snapshots.
    Rows where unbiased_count == 0 are stored as 0.0.
    """
    data = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 4:
                continue
            data[(int(p[0]), int(p[1]))] = float(p[3])
    return data


def parse_steps_from_hist(path):
    """
    Read the simulation time from the header of last_hist.dat.

    The first line has the format:
        #t = 4930000; extr. Ts:

    Returns the step count as an int, or 0 if not parseable.
    """
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('#t'):
                    m = re.search(r'#t\s*=\s*([\d.eE+]+)', line)
                    if m:
                        return int(float(m.group(1)))
    except (OSError, ValueError):
        pass
    return 0


def _fmt_steps(n):
    """Format a step count compactly for plot labels: 4930000 → '4.9M'."""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f}M"
    if n >= 1_000:
        return f"{n / 1_000:.0f}k"
    return str(n)


def parse_traj_hist(path):
    """
    Parse traj_hist.dat → list of (t, hist_dict) tuples.

    traj_hist.dat is a sequence of cumulative histogram snapshots, one per
    print_energy_every steps.  Each snapshot looks identical to last_hist.dat:

        #t = 10000; extr. Ts:
        0 0 0 0
        1 0 0 0
        ...

    The histograms are CUMULATIVE: snapshot at t=N contains counts from
    step 0 to step N.  To get counts for a sub-interval [t1, t2], subtract:
        counts_interval = hist[t2][state] - hist[t1][state]

    Returns a list of (t: int, hist: dict) sorted by t.
    Returns an empty list if the file does not exist or is unreadable.
    """
    snapshots = []
    try:
        with open(path) as f:
            t_cur  = None
            h_cur  = {}
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#t'):
                    # Save previous snapshot before starting a new one
                    if t_cur is not None:
                        snapshots.append((t_cur, h_cur))
                    m = re.search(r'#t\s*=\s*([\d.eE+]+)', line)
                    t_cur = int(float(m.group(1))) if m else 0
                    h_cur = {}
                elif line.startswith('#'):
                    continue
                else:
                    p = line.split()
                    if len(p) >= 4:
                        h_cur[(int(p[0]), int(p[1]))] = float(p[3])
            if t_cur is not None:
                snapshots.append((t_cur, h_cur))
    except (OSError, ValueError):
        pass
    return snapshots


def _dg_from_hist(hist, ref_state, tgt_state):
    """
    Compute ΔG = -ln(u_target / u_ref) [kT] from a histogram dict.
    Returns nan if either state has zero or negative counts.
    """
    u_ref = hist.get(ref_state, 0.0)
    u_tgt = hist.get(tgt_state, 0.0)
    if u_ref > 0.0 and u_tgt > 0.0:
        return -np.log(u_tgt / u_ref)
    return float('nan')


def parse_hb_obs(path):
    """
    Parse hb_obs.dat (HB window observable, 100-step resolution).
    Format: step  hb_count
    Returns (times, hb_counts) as numpy arrays, or two empty arrays if missing.
    """
    try:
        data = np.loadtxt(path, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] >= 2:
            return data[:, 0], data[:, 1]
    except (OSError, ValueError):
        pass
    return np.array([]), np.array([])


def parse_dist_obs(path):
    """
    Parse dist_obs.dat (dist window observable, 100-step resolution).
    Format: step  d1  d2  d3  d4  d5  d6  (7 columns)
    Returns (times, mindist) where mindist = min(d1..d6) per row.
    The min reproduces the OP mindistance observable used for binning.
    Returns two empty arrays if file is missing or has wrong format.
    """
    try:
        data = np.loadtxt(path, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] >= 7:
            return data[:, 0], np.min(data[:, 1:7], axis=1)
    except (OSError, ValueError):
        pass
    return np.array([]), np.array([])


def count_roundtrips(signal, lo_thresh, hi_thresh):
    """
    Count full round trips in a 1D observable signal.

    A half-trip is a one-way crossing from one extreme to the other:
      lo → hi  or  hi → lo  (ignoring intermediate values).
    A full round trip = 2 half-trips.

    Parameters
    ----------
    signal     : 1D array of observable values
    lo_thresh  : values ≤ this are "at low endpoint"
    hi_thresh  : values ≥ this are "at high endpoint"

    Returns
    -------
    int : number of full round trips (half-trips // 2)
    """
    if len(signal) == 0:
        return 0
    half_trips = 0
    last = None
    for val in signal:
        if val <= lo_thresh:
            cur = 'lo'
        elif val >= hi_thresh:
            cur = 'hi'
        else:
            continue
        if last is not None and cur != last:
            half_trips += 1
        last = cur
    return half_trips // 2


def compute_running_dg(snapshots, ref_state, tgt_state):
    """
    Compute ΔG at each snapshot from the cumulative histograms.

    Parameters
    ----------
    snapshots  : list of (t, hist_dict)  from parse_traj_hist
    ref_state  : (bond, mdbin) tuple  — reference / stable state
    tgt_state  : (bond, mdbin) tuple  — rare / target state

    Returns
    -------
    ts  : np.ndarray of step values
    dgs : np.ndarray of ΔG [kT], NaN for snapshots where either state = 0
    """
    ts  = np.array([t for t, _ in snapshots], dtype=float)
    dgs = np.array([_dg_from_hist(h, ref_state, tgt_state)
                    for _, h in snapshots], dtype=float)
    return ts, dgs


def block_dg_from_traj(snapshots, ref_state, tgt_state):
    """
    Split trajectory in two equal halves and compute ΔG for each.

    Because traj_hist histograms are cumulative, the second-half counts are:
        h2[state] = h_final[state] - h_mid[state]

    This avoids any dependence on energy.dat resolution.

    Returns
    -------
    block1_dg  ΔG from first half [kT]
    block2_dg  ΔG from second half [kT]
    block_diff |block2_dg - block1_dg| [kT] — convergence indicator;
               close to 0 means the two halves agree (good convergence)
    """
    if len(snapshots) < 2:
        return float('nan'), float('nan'), float('nan')

    mid = len(snapshots) // 2
    _, h_mid   = snapshots[mid - 1]   # end of first half
    _, h_final = snapshots[-1]         # end of second half

    # Second-half histogram = cumulative[final] - cumulative[mid]
    h2 = {}
    for state in h_final:
        diff = h_final[state] - h_mid.get(state, 0.0)
        if diff > 0.0:
            h2[state] = diff

    dg1 = _dg_from_hist(h_mid,   ref_state, tgt_state)
    dg2 = _dg_from_hist(h2,      ref_state, tgt_state)

    if _valid(dg1) and _valid(dg2):
        diff = abs(dg2 - dg1)
    else:
        diff = float('nan')

    return dg1, dg2, diff


# ══════════════════════════════════════════════════════════════════════════════
# Metadata inference from path
# ══════════════════════════════════════════════════════════════════════════════

def infer_system(path):
    """Extract pulling / noclamps / pushing from a directory path string."""
    p = path.lower()
    if 'pulling' in p:   return 'pulling'
    if 'noclamps' in p:  return 'noclamps'
    if 'pushing' in p:   return 'pushing'
    return 'unknown'


def infer_version(path):
    """Extract the vN label from a path string (e.g. 'v5')."""
    m = re.search(r'[\\/](v\d+)[\\/]?', path)
    if m:
        return m.group(1)
    # fallback: last path component if it matches vN
    m = re.search(r'\b(v\d+)\b', path)
    return m.group(1) if m else 'unknown'


# ══════════════════════════════════════════════════════════════════════════════
# Core metric computation  (from last_hist.dat only)
# ══════════════════════════════════════════════════════════════════════════════

def compute_metrics(hist, is_dist):
    """
    Compute all quality metrics from a last_hist.dat dictionary.

    HB window (is_dist=False)
    ─────────────────────────
    Active states:  bond = 0..6, mdbin = 0
    Reference:      bond = 6, mdbin = 0  (fully bonded — highest u, most stable)
    Target:         bond = 0, mdbin = 0  (fully melted — rarest state)

    Dist window (is_dist=True)
    ──────────────────────────
    Active states:  bond = 0, mdbin = 0..6
    Reference:      bond = 0, mdbin = 0              (unbound, strands close together)
    Target:         bond = 0, mdbin = DIST_TARGET_BIN (= 6, mindist > 12 oxDNA units)
    Bin 6 is the true fully-separated state; geometric max confirmed at ~16.5
    oxDNA units across all systems. Seeds that never visit bin 6 get NaN ΔF.

    Returns
    ───────
    u_ref, u_target   unbiased counts at the two endpoint states
    u_rare            minimum non-zero u across all active states — the worst-
                      sampled state (bottleneck for free-energy convergence)
    u_max             maximum u (almost always the reference state)
    flatness_ratio    u_max / u_rare  →  1 = perfect flat histogram
    bins_visited      comma-separated list of visited state indices
    n_bins_visited    number of visited active states
    target_bin        which bin index was used as target (informational)
    dg_kT             ΔG = −ln(u_target / u_ref)  [kT units]
                        HB:   G(bond=0) − G(bond=6)   should be > 0
                        dist: G(mdbin=K) − G(mdbin=0)
    sigma_naive_kT    sqrt(1/u_target + 1/u_ref)  [kT units]
                        This is the LOWER BOUND on std(ΔG) assuming every
                        accumulated count is an independent sample.  The true
                        std is ≥ this because consecutive counts are correlated.
                        Use cross-seed std for the actual uncertainty estimate.
                        Useful here as a relative ranking: lower = better sampled.
    """
    if is_dist:
        # Collect every (bond=0, mdbin) entry present in the file
        all_mdbins = sorted(mb for (b, mb) in hist if b == 0)
        visited    = [(0, mb) for mb in all_mdbins if hist.get((0, mb), 0.0) > 0.0]
        ref_state  = (0, 0)
        # Fixed target: bin 6 (mindist > 12 oxDNA units) = true fully-separated state.
        # Geometric max confirmed at ~16.5 units across all systems/versions.
        tgt_state  = (0, DIST_TARGET_BIN)
        bins_label = [mb for (_, mb) in visited]
        n_expected = max(all_mdbins) + 1 if all_mdbins else 0
    else:
        visited    = [(b, 0) for b in range(7) if hist.get((b, 0), 0.0) > 0.0]
        ref_state  = (6, 0)
        tgt_state  = (0, 0)
        bins_label = [b for (b, _) in visited]
        n_expected = 7

    u_vals   = [hist[s] for s in visited]
    u_rare   = min(u_vals) if u_vals else 0.0
    u_max    = max(u_vals) if u_vals else 0.0
    flatness = u_max / u_rare if u_rare > 0.0 else float('nan')

    u_ref    = hist.get(ref_state,  0.0)
    u_target = hist.get(tgt_state,  0.0)

    if u_ref > 0.0 and u_target > 0.0:
        dg_kT          = -np.log(u_target / u_ref)
        sigma_naive_kT =  np.sqrt(1.0 / u_target + 1.0 / u_ref)
    else:
        dg_kT          = float('nan')
        sigma_naive_kT = float('nan')

    return {
        'u_ref':           u_ref,
        'u_target':        u_target,
        'u_rare':          u_rare,
        'u_max':           u_max,
        'flatness_ratio':  flatness,
        'bins_visited':    ','.join(map(str, bins_label)),
        'n_bins_visited':  len(visited),
        'n_bins_expected': n_expected,
        'target_bin':      tgt_state[1] if is_dist else tgt_state[0],
        'dg_kT':           dg_kT,
        'sigma_naive_kT':  sigma_naive_kT,
    }


# ══════════════════════════════════════════════════════════════════════════════
# Level 1 gate-keeping
# ══════════════════════════════════════════════════════════════════════════════

def level1_status(n_visited, n_expected, is_dist):
    """
    PASS  — all expected bins visited
    FAIL  — HB window: at least one bond state (0–6) never visited
              → ΔG is undefined; seed must be excluded
    FLAG  — dist window: some expected bins not visited
              → soft warning; may be an interface range issue, not a bad sim
    """
    if n_visited >= n_expected:
        return 'PASS'
    return 'FLAG' if is_dist else 'FAIL'


# ══════════════════════════════════════════════════════════════════════════════
# Run discovery and analysis
# ══════════════════════════════════════════════════════════════════════════════

def discover_runs(base_dirs):
    """
    Walk each given directory, match subdirectory names against the run pattern,
    return a list of metadata dicts — one per run found.
    """
    runs = []
    for base in base_dirs:
        base = os.path.realpath(base)
        if not os.path.isdir(base):
            print(f"WARNING: {base} not found, skipping", file=sys.stderr)
            continue
        system  = infer_system(base)
        version = infer_version(base)
        for entry in sorted(os.listdir(base)):
            m = RUN_RE.match(entry)
            if not m:
                continue
            runs.append({
                'run_dir': os.path.join(base, entry),
                'source':  base,
                'system':  system,
                'version': version,
                'is_dist': m.group(1) is not None,
                'window':  'dist' if m.group(1) else 'HB',
                'combo':   m.group(2),
                'seed':    m.group(3),
            })
    return runs


def analyze_run(meta):
    """
    Analyze one run directory.  Returns a flat dict ready for CSV output.
    Merges metadata fields with computed metrics.
    """
    result = dict(meta)   # copy metadata
    hist_path = os.path.join(meta['run_dir'], 'last_hist.dat')

    if not os.path.exists(hist_path):
        result.update({
            'status': 'NO_DATA',
            'dg_kT': float('nan'), 'dg_kcal': float('nan'),
            'sigma_naive_kT': float('nan'),
            'n_bins_visited': 0, 'n_bins_expected': 0,
            'bins_visited': '', 'target_bin': float('nan'),
            'u_rare': float('nan'), 'u_max': float('nan'),
            'flatness_ratio': float('nan'),
            'u_target': float('nan'), 'u_ref': float('nan'),
        })
        return result

    hist    = parse_last_hist(hist_path)
    metrics = compute_metrics(hist, meta['is_dist'])
    result.update(metrics)
    result['status'] = level1_status(
        metrics['n_bins_visited'], metrics['n_bins_expected'], meta['is_dist']
    )
    result['dg_kcal']   = metrics['dg_kT'] * kT if _valid(metrics['dg_kT']) else float('nan')
    result['n_steps']   = parse_steps_from_hist(hist_path)
    init_hist_path      = os.path.join(meta['run_dir'], 'init_hist.dat')
    result['init_steps'] = (parse_steps_from_hist(init_hist_path)
                            if os.path.exists(init_hist_path) else 0)

    # ── traj_hist.dat block convergence ──────────────────────────────────────
    traj_path = os.path.join(meta['run_dir'], 'traj_hist.dat')
    if meta['is_dist']:
        ref_state = (0, 0)
        tgt_state = (0, DIST_TARGET_BIN)
    else:
        ref_state = (6, 0)
        tgt_state = (0, 0)

    result['ref_state'] = ref_state
    result['tgt_state'] = tgt_state

    if os.path.exists(traj_path):
        snapshots = parse_traj_hist(traj_path)
        result['traj_n_snapshots'] = len(snapshots)
        b1, b2, bdiff = block_dg_from_traj(snapshots, ref_state, tgt_state)
        result['traj_block1_dg_kT']  = b1
        result['traj_block2_dg_kT']  = b2
        result['traj_block_diff_kT'] = bdiff
        # Store running ΔG time series for convergence plot
        ts, dgs = compute_running_dg(snapshots, ref_state, tgt_state)
        result['_traj_ts']  = ts
        result['_traj_dgs'] = dgs
    else:
        result['traj_n_snapshots']   = 0
        result['traj_block1_dg_kT']  = float('nan')
        result['traj_block2_dg_kT']  = float('nan')
        result['traj_block_diff_kT'] = float('nan')
        result['_traj_ts']  = np.array([])
        result['_traj_dgs'] = np.array([])

    # ── Observable time series (hb_obs / dist_obs) ───────────────────────────
    if meta['is_dist']:
        obs_path = os.path.join(meta['run_dir'], 'dist_obs.dat')
        obs_ts, obs_signal = parse_dist_obs(obs_path)
        lo = meta.get('rt_dist_lo', 2.0)
        hi = meta.get('rt_dist_hi', 10.0)
        rt = count_roundtrips(obs_signal, lo_thresh=lo, hi_thresh=hi) if len(obs_signal) > 0 else 0
    else:
        obs_path = os.path.join(meta['run_dir'], 'hb_obs.dat')
        obs_ts, obs_signal = parse_hb_obs(obs_path)
        lo = meta.get('rt_hb_lo', 0.0)
        hi = meta.get('rt_hb_hi', 6.0)
        rt = count_roundtrips(obs_signal, lo_thresh=lo, hi_thresh=hi) if len(obs_signal) > 0 else 0

    result['_obs_ts']     = obs_ts
    result['_obs_signal'] = obs_signal
    result['n_roundtrips']  = rt
    result['obs_n_points']  = len(obs_signal)

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Continuation stitching  (chain v4→v7→v8 etc. into one trace)
# ══════════════════════════════════════════════════════════════════════════════

def stitch_continuations(rows):
    """
    For continuation runs (those with an init_hist.dat), stitch the observable
    time series and traj_hist convergence data together with their source runs so
    the full simulation history is visible as one continuous trace.

    Chain reconstruction: a row is matched to a source row in the same
    (system, window, combo) group.  We do NOT require exact matching of
    init_steps to source.n_steps — in practice the source run often continues
    past the snapshot point, so we pick the source with n_steps closest to
    (and ideally ≥) init_steps.  Seed is ignored in matching because all v8
    seeds of a combo share the same init_hist from a single v7 source.

    Stitching details:
      * source data is truncated to t ≤ init_steps (drops any part of the
        source run that happened AFTER the branch point — that data is real
        but belongs to a different trajectory fork).
      * link's own data is offset by (source.offset + init_steps), so
        absolute timeline = branch point of source + link's own progress.

    After stitching, each row gets:
        total_steps   : offset + n_steps in absolute (stitched) timeline
        _traj_ts      : absolute-step traj_hist timestamps
        _traj_dgs     : corresponding ΔG values
        _obs_ts       : absolute-step observable timestamps
        _obs_signal   : corresponding observable values
        n_roundtrips  : re-counted from the full stitched obs signal
    """
    # Group candidates by (system, window, combo) — seed-agnostic.
    by_combo = defaultdict(list)
    for r in rows:
        if r.get('n_steps', 0) > 0:
            key = (r['system'], r['window'], r['combo'])
            by_combo[key].append(r)

    def find_source(r):
        init_s = r.get('init_steps', 0)
        if init_s <= 0:
            return None
        key = (r['system'], r['window'], r['combo'])
        candidates = [c for c in by_combo.get(key, []) if c is not r]
        if not candidates:
            return None
        # Prefer sources whose n_steps ≥ init_steps (they actually contain the
        # branch point in their own trajectory). Among those, pick the one
        # closest to init_steps.
        above = [c for c in candidates if c.get('n_steps', 0) >= init_s]
        if above:
            return min(above, key=lambda c: c['n_steps'] - init_s)
        # Fallback: the best we have — pick the one with the largest n_steps
        return max(candidates, key=lambda c: c.get('n_steps', 0))

    def build_chain(r):
        """Return [oldest_ancestor, ..., r] following source links."""
        chain, seen = [], set()
        cur = r
        while cur is not None:
            rid = id(cur)
            if rid in seen:
                break  # guard against unexpected cycles
            seen.add(rid)
            chain.insert(0, cur)
            cur = find_source(cur)
        return chain

    n_stitched = 0
    for r in rows:
        chain = build_chain(r)
        if len(chain) <= 1:
            r.setdefault('total_steps', r.get('n_steps', 0))
            continue
        n_stitched += 1

        # Compute per-link absolute offsets.
        # Root: offset = 0.
        # Link i (i > 0): offset = offsets[i-1] + chain[i].init_steps
        #   — i.e. source's absolute start + branch point within source's own timeline.
        offsets = [0]
        for i in range(1, len(chain)):
            offsets.append(offsets[i - 1] + chain[i].get('init_steps', 0))
        r['total_steps'] = offsets[-1] + r.get('n_steps', 0)

        # Per-link truncation: a non-last link is truncated at the init_steps
        # of the NEXT link (the branch point) so source data past that point
        # — which belongs to a different trajectory fork — is dropped.
        def truncate_at(i):
            if i < len(chain) - 1:
                return chain[i + 1].get('init_steps', float('inf'))
            return float('inf')

        # Stitch traj_hist convergence traces
        tts_parts, tdg_parts = [], []
        for i, link in enumerate(chain):
            off = offsets[i]
            cut = truncate_at(i)
            ts = link.get('_traj_ts', np.array([]))
            dg = link.get('_traj_dgs', np.array([]))
            if len(ts) > 0:
                mask = ts <= cut
                tts_parts.append(ts[mask] + off)
                tdg_parts.append(dg[mask])
        if tts_parts:
            r['_traj_ts']  = np.concatenate(tts_parts)
            r['_traj_dgs'] = np.concatenate(tdg_parts)

        # Stitch observable time series
        ots_parts, osig_parts = [], []
        for i, link in enumerate(chain):
            off = offsets[i]
            cut = truncate_at(i)
            ts  = link.get('_obs_ts',     np.array([]))
            sig = link.get('_obs_signal', np.array([]))
            if len(ts) > 0:
                mask = ts <= cut
                ots_parts.append(ts[mask] + off)
                osig_parts.append(sig[mask])
        if ots_parts:
            r['_obs_ts']     = np.concatenate(ots_parts)
            r['_obs_signal'] = np.concatenate(osig_parts)
            r['obs_n_points'] = len(r['_obs_signal'])
            lo = (r.get('rt_dist_lo', 2.0) if r['is_dist']
                  else r.get('rt_hb_lo', 0.0))
            hi = (r.get('rt_dist_hi', 12.0) if r['is_dist']
                  else r.get('rt_hb_hi', 6.0))
            r['n_roundtrips'] = count_roundtrips(
                r['_obs_signal'], lo_thresh=lo, hi_thresh=hi)

    if n_stitched > 0:
        print(f"Stitched {n_stitched} continuation run(s) to their source runs.")


# ══════════════════════════════════════════════════════════════════════════════
# Cross-seed statistics  (added after all runs are analysed)
# ══════════════════════════════════════════════════════════════════════════════

def add_cross_seed_stats(rows, sigma_max=float('inf')):
    """
    For each group (source, system, window, combo), compute cross-seed stats
    in two passes:

    Pass 1 — using all non-FAIL seeds with valid ΔG:
        seed_n        number of usable seeds
        seed_mean_dg  mean ΔG
        seed_std_dg   std ΔG  ← real uncertainty; seeds are independent so this
                      automatically captures autocorrelation without needing
                      round-trip counts or block analysis
        seed_outlier  True if |ΔG − group mean| > 2 × seed_std_dg

    Pass 2 — filtered estimate:
        Seeds are kept only if BOTH conditions hold:
          (a) sigma_naive_kT ≤ sigma_max  — rare state was sampled enough
          (b) not a 2σ outlier from the Pass-1 group mean

        sigma_max filters on per-seed data quality independently of other seeds.
        Outlier detection then removes seeds inconsistent with their peers.
        When sigma is large for ALL seeds (common for HB window), the 2σ rule
        alone does nothing — sigma_max is what actually bites in that case.

        filtered_n        number of seeds remaining after both filters
        filtered_mean_dg  mean ΔG of kept seeds  ← use this for reporting
        filtered_std_dg   std ΔG of kept seeds
        include_in_filter True if this seed contributes to the filtered estimate
    """
    groups = defaultdict(list)
    for r in rows:
        key = (r['source'], r['system'], r['window'], r['combo'])
        groups[key].append(r)

    for key, group in groups.items():
        # Pass 1: all non-FAIL seeds with valid ΔG (unfiltered baseline)
        usable    = [r for r in group if r['status'] != 'FAIL' and _valid(r['dg_kT'])]
        usable_dg = [r['dg_kT'] for r in usable]
        n    = len(usable_dg)
        mean = float(np.mean(usable_dg))        if n >= 1 else float('nan')
        std  = float(np.std(usable_dg, ddof=1)) if n >= 2 else float('nan')

        for r in group:
            r['seed_n']       = n
            r['seed_mean_dg'] = mean
            r['seed_std_dg']  = std
            dg  = r['dg_kT']
            sig = r.get('sigma_naive_kT', float('nan'))
            if _valid(dg) and _valid(std) and std > 0:
                r['seed_outlier'] = bool(abs(dg - mean) > 2 * std)
            else:
                r['seed_outlier'] = False
            # Mark seeds that fail the sigma threshold (for plot transparency)
            r['sigma_filtered'] = _valid(sig) and sig > sigma_max

        # Pass 2: apply sigma threshold AND outlier filter
        kept = [
            r for r in usable
            if not r['seed_outlier']
            and _valid(r.get('sigma_naive_kT', float('nan')))
            and r.get('sigma_naive_kT') <= sigma_max
        ]
        kept_dg = [r['dg_kT'] for r in kept]
        fn    = len(kept_dg)
        fmean = float(np.mean(kept_dg))        if fn >= 1 else float('nan')
        fstd  = float(np.std(kept_dg, ddof=1)) if fn >= 2 else float('nan')
        kept_ids = set(id(r) for r in kept)

        for r in group:
            r['filtered_n']        = fn
            r['filtered_mean_dg']  = fmean
            r['filtered_std_dg']   = fstd
            r['include_in_filter'] = id(r) in kept_ids


# ══════════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════════

# σ_naive colour map: green (good, low σ) → yellow → red (bad, high σ)
# Log scale because the range spans several orders of magnitude.
SIGMA_CMAP = plt.cm.RdYlGn_r
SIGMA_NORM = mcolors.LogNorm(vmin=1e-3, vmax=1.0)

# Marker shapes encode Level 1 status
MARKER_SHAPE = {'PASS': 'o', 'FLAG': 's', 'FAIL': 'X', 'NO_DATA': '.'}
MARKER_EDGE  = {'PASS': 'none', 'FLAG': 'k', 'FAIL': 'darkred', 'NO_DATA': 'grey'}


def _sigma_color(sigma):
    """Map σ_naive value to an RGBA colour via the log-normalised colourmap."""
    if not _valid(sigma):
        return (0.7, 0.7, 0.7, 1.0)   # grey for nan
    return SIGMA_CMAP(SIGMA_NORM(np.clip(sigma, 1e-3, 1.0)))


def _sigma_alpha(sigma, faded=False):
    """
    Map σ_naive to a line alpha for trace plots where colour encodes seed identity.
    faded=True (outlier/FAIL) overrides to a low fixed alpha.
    """
    if faded:
        return 0.20
    if not _valid(sigma):
        return 0.45
    # sigma < 0.05 → opaque; sigma > 0.5 → dim but still visible
    return float(np.clip(0.9 - 0.6 * np.log10(max(sigma, 1e-3) / 0.01) / np.log10(0.5 / 0.01), 0.25, 0.9))


def _seed_color_map(subset):
    """
    Assign a distinct categorical colour (tab10) to each unique (combo, seed) pair
    in the given subset.  Returns a dict: seed_label → RGBA colour.
    Seed label = 'wA_s1' etc.
    """
    labels = sorted(set(f"{r['combo']}_{r['seed']}" for r in subset))
    cmap   = plt.cm.tab10
    return {lab: cmap(i % 10) for i, lab in enumerate(labels)}


def _xtick_labels(rows):
    """
    Build sorted x-tick labels of the form "version-combo\nt=<median steps>"
    from the versions present in a subset of rows.

    The second line shows the median step count across seeds in that group,
    formatted compactly (e.g. 4.9M).  This lets you instantly see which
    groups ran longer without leaving the plot.
    """
    versions = sorted(set(r['version'] for r in rows))
    labels = []
    for v in versions:
        for c in COMBOS:
            group_rows = [r for r in rows if r['version'] == v and r['combo'] == c]
            if not group_rows:
                continue  # skip combos with no data in this subset
            group_steps = [r.get('total_steps', r.get('n_steps', 0))
                           for r in group_rows
                           if r.get('total_steps', r.get('n_steps', 0)) > 0]
            base = f"{v}-{c}"
            if group_steps:
                med_t = int(np.median(group_steps))
                labels.append(f"{base}\nt={_fmt_steps(med_t)}")
            else:
                labels.append(base)
    return labels


def _make_axes_grid(nrows, ncols, figsize, suptitle):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    fig.suptitle(suptitle, fontsize=11)
    return fig, axes


def _add_sigma_colorbar(fig):
    """Add a small horizontal sigma colour bar at the top-right, outside the subplot grid."""
    sm = ScalarMappable(cmap=SIGMA_CMAP, norm=SIGMA_NORM)
    sm.set_array([])
    cax = fig.add_axes([0.78, 0.945, 0.18, 0.013])
    cb = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cb.set_label('σ_naive  (kT)', fontsize=8)
    cax.tick_params(labelsize=7)


def make_sigma_plot(rows, out_path):
    """
    Figure 1 — σ_naive strip plot.

    Layout:  2 rows (HB / dist)  ×  3 columns (pulling / noclamps / pushing)

    Each panel shows one seed per point:
      x-axis  version × combo  (e.g. v4-wA, v5-wA, v5-wB, v5-wC)
      y-axis  σ_naive  [kT units, log scale]
      colour  green → red via σ_naive  (green = good, red = bad)
      shape   ○ PASS   □ FLAG (dist missing bins)   ✕ FAIL (HB missing bonds)

    Reference lines at 0.5 kT and 0.1 kT help orient the eye.  No hard threshold
    is applied — the distribution itself tells you where "good" sits for your data.
    """
    np.random.seed(42)   # reproducible jitter

    fig, axes = _make_axes_grid(
        2, 3, figsize=(14, 8),
        suptitle=(
            'σ_naive per seed  [lower bound on ΔG error, kT units]\n'
            'Assumes all histogram counts are independent — true σ ≥ this value.\n'
            'Use cross-seed std for the real uncertainty estimate.  '
            '▲ = σ > axis max (target barely visited).  ▼ = σ < axis min (very well sampled).  Grey = NaN.'
        )
    )

    for ri, win in enumerate(WINDOWS):
        for ci, sys_ in enumerate(SYSTEMS):
            ax     = axes[ri][ci]
            subset = [r for r in rows if r['window'] == win and r['system'] == sys_]

            if not subset:
                ax.set_visible(False)
                continue

            labels   = _xtick_labels(subset)
            pos_map  = {lb.split('\n')[0]: i for i, lb in enumerate(labels)}

            for r in subset:
                sig   = r.get('sigma_naive_kT', float('nan'))
                label = f"{r['version']}-{r['combo']}"
                xi    = pos_map.get(label, 0) + np.random.uniform(-0.18, 0.18)
                color  = _sigma_color(sig)
                marker = MARKER_SHAPE.get(r['status'], 'o')
                edge   = MARKER_EDGE.get(r['status'],  'none')
                SIGMA_YMIN = 1e-4
                SIGMA_YMAX = 10.0
                if _valid(sig):
                    if sig < SIGMA_YMIN:
                        # σ below axis range — pin to bottom with downward triangle
                        ax.scatter(xi, SIGMA_YMIN, c=[color], marker='v',
                                   edgecolors=edge, linewidths=0.8, s=55, zorder=3,
                                   clip_on=False)
                    elif sig <= SIGMA_YMAX:
                        ax.scatter(xi, sig, c=[color], marker=marker,
                                   edgecolors=edge, linewidths=0.8, s=55, zorder=3)
                    else:
                        # σ above axis range — pin to top with upward triangle
                        ax.scatter(xi, SIGMA_YMAX, c=[color], marker='^',
                                   edgecolors=edge, linewidths=0.8, s=55, zorder=3,
                                   clip_on=False)
                else:
                    # NaN σ (FAIL / NO_DATA) — plot a grey marker at bottom of axis
                    ax.scatter(xi, SIGMA_YMIN, c=['lightgrey'], marker=marker,
                               edgecolors=edge, linewidths=0.8, s=55, zorder=3,
                               clip_on=False)

            ax.set_yscale('log')
            ax.set_ylim(SIGMA_YMIN, SIGMA_YMAX)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
            ax.set_title(f'{sys_}  —  {win}', fontsize=10)
            if ci == 0:
                ax.set_ylabel('σ_naive  (kT)', fontsize=9)
            ax.axhline(0.5, color='dimgrey', ls='--', lw=0.8, alpha=0.7)
            ax.axhline(0.1, color='dimgrey', ls=':',  lw=0.8, alpha=0.7)
            ax.text(0.99, 0.52, '0.5 kT', ha='right', va='bottom',
                    fontsize=7, color='dimgrey', transform=ax.get_yaxis_transform())
            ax.text(0.99, 0.103, '0.1 kT', ha='right', va='bottom',
                    fontsize=7, color='dimgrey', transform=ax.get_yaxis_transform())
            ax.grid(axis='y', alpha=0.25)

    # Shared colour bar
    _add_sigma_colorbar(fig)

    # Marker legend
    from matplotlib.lines import Line2D
    leg = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor='grey',
               markersize=9, label='PASS  (all bins visited)'),
        Line2D([0],[0], marker='s', color='w', markerfacecolor='grey',
               markeredgecolor='k', markersize=9,
               label='FLAG  (dist: some bins missing — check interface range)'),
        Line2D([0],[0], marker='X', color='w', markerfacecolor='grey',
               markeredgecolor='darkred', markersize=9,
               label='FAIL  (HB: bond state missing — seed excluded)'),
    ]
    fig.legend(handles=leg, loc='lower center', ncol=3,
               bbox_to_anchor=(0.47, 0.01), fontsize=8)

    fig.tight_layout(rect=[0, 0.10, 1.0, 0.91])
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def make_dg_plot(rows, out_path, min_roundtrips=0):
    """
    Figure 2 — ΔG per seed, coloured by σ_naive.

    Same 2 × 3 panel layout as Figure 1.

    Each point is one seed's ΔG estimate, coloured by σ_naive so you can see
    at a glance which seeds to trust.

    Marker encoding:
      ○  normal seed included in the filtered estimate
      ✕  FAIL seed (missing bond state — never used)
      ◇  outlier seed (|ΔG − group mean| > 2 × seed_std) — excluded from bars

    Error bar encoding (per combo):
      solid black bar   filtered mean ± filtered_std  (outliers removed)
      dashed grey bar   unfiltered mean ± seed_std    (all non-FAIL seeds)

    Showing both lets you see how much the outliers were pulling the estimate.
    The solid bar is what should be used for reporting.
    """
    np.random.seed(42)

    fig, axes = _make_axes_grid(
        2, 3, figsize=(14, 8),
        suptitle='ΔG per seed  [kT units],  coloured by σ_naive\n'
                 'Solid bar = filtered mean ± std (outliers removed, n≥5)  '
                 '|  Dashed bar = unfiltered mean ± std (all non-FAIL, n≥5)\n'
                 'Number above each dot = round-trip count.  ◇ = outlier excluded from filtered bar'
    )

    for ri, win in enumerate(WINDOWS):
        for ci, sys_ in enumerate(SYSTEMS):
            ax     = axes[ri][ci]
            subset = [r for r in rows
                      if r['window'] == win and r['system'] == sys_
                      and r['status'] not in ('NO_DATA',)]

            if not subset:
                ax.set_visible(False)
                continue

            labels  = _xtick_labels(subset)
            pos_map = {lb.split('\n')[0]: i for i, lb in enumerate(labels)}

            # Individual seeds
            for r in subset:
                dg  = r.get('dg_kT', float('nan'))
                sig = r.get('sigma_naive_kT', float('nan'))
                if not _valid(dg):
                    continue
                label  = f"{r['version']}-{r['combo']}"
                xi_pos = pos_map.get(label)
                if xi_pos is None:
                    continue
                xi     = xi_pos + np.random.uniform(-0.18, 0.18)
                color  = _sigma_color(sig)
                is_out = r.get('seed_outlier', False)
                is_fail = r.get('status') == 'FAIL'
                marker = 'X' if is_fail else ('D' if is_out else 'o')
                edge   = 'darkred' if is_fail else ('black' if is_out else 'none')
                alpha  = 0.35 if (is_out or is_fail) else 1.0
                ax.scatter(xi, dg, c=[color], marker=marker, alpha=alpha,
                           edgecolors=edge, linewidths=0.9, s=55, zorder=3)
                rt = r.get('n_roundtrips', 0)
                ax.annotate(str(rt), (xi, dg),
                            textcoords='offset points', xytext=(0, 5),
                            fontsize=5, ha='center', color='#333333')

            # Collect per-combo filtered and unfiltered stats from the rows
            # (already computed by add_cross_seed_stats — just read them off)
            seen_labels = set()
            for r in subset:
                label = f"{r['version']}-{r['combo']}"
                if label in seen_labels:
                    continue
                seen_labels.add(label)
                xi = pos_map.get(label)
                if xi is None:
                    continue

                if min_roundtrips == 0:
                    # Unfiltered bar (dashed grey) — only meaningful without RT filter
                    umean = r.get('seed_mean_dg', float('nan'))
                    ustd  = r.get('seed_std_dg',  float('nan'))
                    if _valid(umean) and _valid(ustd) and r.get('seed_n', 0) >= 5:
                        ax.errorbar(xi - 0.12, umean, yerr=ustd, fmt='_',
                                    color='grey', capsize=4, lw=1.2,
                                    linestyle='--', zorder=4, alpha=0.7)

                # Filtered bar (solid black) — always shown regardless of RT filter
                fmean = r.get('filtered_mean_dg', float('nan'))
                fstd  = r.get('filtered_std_dg',  float('nan'))
                if _valid(fmean) and _valid(fstd) and r.get('filtered_n', 0) >= 5:
                    xi_f = xi + (0.12 if min_roundtrips == 0 else 0)
                    ax.errorbar(xi_f, fmean, yerr=fstd, fmt='_',
                                color='black', capsize=4, lw=1.8, zorder=5)
                elif _valid(fmean):
                    xi_f = xi + (0.12 if min_roundtrips == 0 else 0)
                    ax.scatter(xi_f, fmean, marker='_',
                               color='black', s=200, zorder=5)

            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
            ax.set_title(f'{sys_}  —  {win}', fontsize=10)
            if ci == 0:
                ax.set_ylabel('ΔG  (kT)', fontsize=9)
            ax.grid(axis='y', alpha=0.25)

    # Shared colour bar
    _add_sigma_colorbar(fig)

    from matplotlib.lines import Line2D
    leg = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor='grey',
               markersize=9, label='Included seed'),
        Line2D([0],[0], marker='D', color='w', markerfacecolor='grey',
               markeredgecolor='black', markersize=9, label='Outlier (excluded, faded)'),
        Line2D([0],[0], marker='X', color='w', markerfacecolor='grey',
               markeredgecolor='darkred', markersize=9, label='FAIL (excluded, faded)'),
        Line2D([0],[0], marker='_', color='black',
               markersize=14, lw=2, label='Filtered mean ± std  (outliers removed, solid)'),
        Line2D([0],[0], marker='_', color='grey', alpha=0.7,
               markersize=14, lw=1.5, label='Unfiltered mean ± std  (all non-FAIL, dashed)'),
    ]
    fig.legend(handles=leg, loc='lower center', ncol=5,
               bbox_to_anchor=(0.47, 0.01), fontsize=8)

    fig.tight_layout(rect=[0, 0.10, 1.0, 0.91])
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def make_convergence_plot(rows, out_path):
    """
    Figure 3 — ΔG(t) convergence time series.

    Layout:  2 rows (HB / dist)  ×  3 columns (pulling / noclamps / pushing)

    Each panel shows one line per seed:
      x-axis  simulation steps
      y-axis  ΔG(t)  [kT], computed from cumulative traj_hist.dat counts
      colour  seed identity (tab10 categorical — same seed = same colour within panel)
      alpha   σ_naive (more opaque = better sampled; faded = outlier/FAIL)

    A converged seed has ΔG(t) that plateaus by the end of the run.
    Dashed vertical line = trajectory midpoint used for block convergence check.
    Seeds with no traj_hist.dat data are silently skipped.
    Note: lines start only once both endpoint states have been visited (ΔG defined).
    """
    fig, axes = _make_axes_grid(
        2, 3, figsize=(14, 9),
        suptitle=(
            'ΔG(t) convergence  [kT units]  —  from cumulative traj_hist.dat counts\n'
            'Colour = seed identity (tab10).  Alpha = σ_naive (faint = poorly sampled or outlier).\n'
            'Lines begin once both endpoint states are first visited.  '
            'Dashed grey = trajectory midpoint.'
        )
    )

    for ri, win in enumerate(WINDOWS):
        for ci, sys_ in enumerate(SYSTEMS):
            ax     = axes[ri][ci]
            subset = [r for r in rows
                      if r['window'] == win and r['system'] == sys_
                      and len(r.get('_traj_ts', [])) > 0]

            if not subset:
                ax.set_visible(False)
                continue

            seed_colors = _seed_color_map(subset)
            legend_handles = []

            for r in subset:
                ts  = r['_traj_ts']
                dgs = r['_traj_dgs']
                if len(ts) == 0:
                    continue
                faded = r.get('seed_outlier', False) or r['status'] == 'FAIL'
                rt    = r.get('n_roundtrips', 0)
                label = f"{r['combo']}_{r['seed']}  ({rt}RT)"
                color = seed_colors[f"{r['combo']}_{r['seed']}"]
                alpha = _sigma_alpha(r.get('sigma_naive_kT', float('nan')), faded=faded)
                lw    = 0.8 if faded else 1.3
                line, = ax.plot(ts / 1e6, dgs, color=color, alpha=alpha, lw=lw,
                                label=label)
                legend_handles.append(line)

                mid_idx = len(ts) // 2
                if mid_idx > 0:
                    ax.axvline(ts[mid_idx - 1] / 1e6,
                               color='dimgrey', lw=0.5, ls='--', alpha=0.25)

            ax.legend(handles=legend_handles, fontsize=6, loc='upper left',
                      framealpha=0.6, handlelength=1.5)
            ax.set_xlabel('steps  (×10⁶)', fontsize=8)
            ax.set_title(f'{sys_}  —  {win}', fontsize=10)
            if ci == 0:
                ax.set_ylabel('ΔG  (kT)', fontsize=9)
            ax.grid(alpha=0.2)

    fig.tight_layout(rect=[0, 0.01, 1.0, 1.0])
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def make_obs_plot(rows, out_path):
    """
    Figure 4 — Observable time traces at 100-step resolution.

    Layout: 2 rows (HB / dist) × 3 columns (pulling / noclamps / pushing)

    HB row:   bond count vs steps, from hb_obs.dat
              y-axis 0–6; dashed grey lines at each integer bond level.
    Dist row: min(d1..d6) vs steps, from dist_obs.dat
              y-axis 0–~14; dashed grey lines at bin interfaces (2,4,6,8,10,12).

    Colour = seed identity (tab10 categorical) so individual seeds are trackable.
    Alpha  = σ_naive (more opaque = better sampled; faded = outlier/FAIL).
    Legend lists each seed with its round-trip count.

    Round trips are defined as:
      HB:   alternating visits to hb ≤ 1  and  hb ≥ 5
      Dist: alternating visits to mindist ≤ 2  and  mindist ≥ 10
    """
    first = next((r for r in rows if r.get('obs_n_points', 0) > 0), {})
    hb_lo  = first.get('rt_hb_lo',   0.0)
    hb_hi  = first.get('rt_hb_hi',   6.0)
    d_lo   = first.get('rt_dist_lo', 2.0)
    d_hi   = first.get('rt_dist_hi', 10.0)

    fig, axes = _make_axes_grid(
        2, 3, figsize=(14, 9),
        suptitle=(
            'Observable time traces  (100-step resolution)\n'
            'HB: bond count from hb_obs.dat  |  '
            'Dist: min(6 pair distances) from dist_obs.dat\n'
            f'RT thresholds — HB: hb ≤ {hb_lo} / hb ≥ {hb_hi}  |  '
            f'Dist: mindist ≤ {d_lo} / mindist ≥ {d_hi}\n'
            'Colour = seed identity (tab10).  Faded = outlier/FAIL.  '
            'Legend shows round-trip count per seed.'
        )
    )

    OBS_YLABELS = {'HB': 'bond count', 'dist': 'min pair dist  (oxDNA units)'}
    HB_REFS   = [0, 1, 2, 3, 4, 5, 6]
    DIST_REFS = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]

    for ri, win in enumerate(WINDOWS):
        for ci, sys_ in enumerate(SYSTEMS):
            ax     = axes[ri][ci]
            subset = [r for r in rows
                      if r['window'] == win and r['system'] == sys_
                      and len(r.get('_obs_ts', [])) > 0]

            if not subset:
                ax.set_visible(False)
                continue

            seed_colors   = _seed_color_map(subset)
            legend_handles = []

            for r in subset:
                ts  = r['_obs_ts']
                sig = r['_obs_signal']
                if len(ts) == 0:
                    continue
                faded = r.get('seed_outlier', False) or r['status'] == 'FAIL'
                label = f"{r['combo']}_{r['seed']}"
                color = seed_colors[label]
                alpha = 0.15 if faded else 0.35
                lw    = 0.4
                rt    = r.get('n_roundtrips', 0)
                leg_label = f"{label}  ({rt}RT)"
                line, = ax.plot(ts / 1e6, sig, color=color, alpha=alpha, lw=lw,
                                label=leg_label)
                legend_handles.append(line)

            refs = HB_REFS if win == 'HB' else DIST_REFS
            for ref in refs:
                ax.axhline(ref, color='dimgrey', lw=0.35, ls=':', alpha=0.5)

            ax.legend(handles=legend_handles, fontsize=6, loc='best',
                      framealpha=0.6, handlelength=1.5)

            if win == 'HB':
                ax.set_ylim(-0.3, 6.5)
            else:
                ax.set_ylim(-0.3, 14.0)

            ax.set_xlabel('steps  (×10⁶)', fontsize=8)
            ax.set_ylabel(OBS_YLABELS[win], fontsize=8)
            ax.set_title(f'{sys_}  —  {win}', fontsize=10)
            ax.grid(alpha=0.15)

    fig.tight_layout(rect=[0, 0.01, 1.0, 1.0])
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# CSV output helpers
# ══════════════════════════════════════════════════════════════════════════════

CSV_FIELDS = [
    'source', 'system', 'version', 'window', 'combo', 'seed', 'status',
    'bins_visited', 'n_bins_visited', 'n_bins_expected', 'target_bin',
    'u_rare', 'u_max', 'flatness_ratio', 'u_target', 'u_ref',
    'dg_kT', 'dg_kcal', 'sigma_naive_kT', 'n_steps',
    'traj_n_snapshots', 'traj_block1_dg_kT', 'traj_block2_dg_kT', 'traj_block_diff_kT',
    'n_roundtrips', 'obs_n_points',
    'seed_n', 'seed_mean_dg', 'seed_std_dg', 'seed_outlier', 'sigma_filtered',
    'filtered_n', 'filtered_mean_dg', 'filtered_std_dg', 'include_in_filter',
]


def _valid(v):
    """Return True if v is a finite float (not nan or inf)."""
    try:
        return np.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _fmt(v):
    """Format a value for CSV output."""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, float):
        return f'{v:.5g}' if _valid(v) else 'nan'
    return str(v) if v is not None else ''


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='VMMC umbrella sampling diagnostics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('dirs', nargs='+', metavar='DIR',
                        help='Version directories to scan (e.g. v4/ v5/)')
    parser.add_argument('--sigma-max', type=float, default=float('inf'),
                        metavar='FLOAT',
                        help=(
                            'Exclude seeds with sigma_naive > this value [kT] '
                            'from the filtered ΔG estimate. '
                            'Look at diagnostics_sigma.pdf first to choose a '
                            'sensible threshold (e.g. --sigma-max 0.2). '
                            'Default: no threshold (inf).'
                        ))
    parser.add_argument('--rt-hb-lo', type=float, default=0.0, metavar='FLOAT',
                        help=(
                            'HB round-trip lower threshold: hb_count ≤ this value '
                            'counts as "at bond=0 endpoint". Default: 0 (exact).'
                        ))
    parser.add_argument('--rt-hb-hi', type=float, default=6.0, metavar='FLOAT',
                        help=(
                            'HB round-trip upper threshold: hb_count ≥ this value '
                            'counts as "at bond=6 endpoint". Default: 6 (exact).'
                        ))
    parser.add_argument('--rt-dist-lo', type=float, default=2.0, metavar='FLOAT',
                        help=(
                            'Dist round-trip lower threshold: mindist ≤ this value '
                            'counts as "at bin-0 endpoint". Default: 2.0 (bin 0 upper edge).'
                        ))
    parser.add_argument('--rt-dist-hi', type=float, default=12.0, metavar='FLOAT',
                        help=(
                            'Dist round-trip upper threshold: mindist ≥ this value '
                            'counts as "at bin-6 endpoint". Default: 12.0 (bin 6 lower edge).'
                        ))
    parser.add_argument('--prefix', type=str, default='', metavar='STR',
                        help=(
                            'Prefix for all output files, e.g. "v8_" produces '
                            'v8_diagnostics.csv, v8_diagnostics_sigma.pdf, etc. '
                            'Default: no prefix.'
                        ))
    parser.add_argument('--min-roundtrips', type=int, default=0, metavar='INT',
                        help=(
                            'Only include seeds with at least this many round trips '
                            'in all plots and stats. Applied to all plot types. '
                            'Default: 0 (no filter). Use 1 to show only seeds that '
                            'completed at least one full round trip.'
                        ))
    args = parser.parse_args()

    sigma_max = args.sigma_max
    if sigma_max < float('inf'):
        print(f"sigma_max = {sigma_max} kT  "
              f"(seeds with σ_naive > {sigma_max} excluded from filtered estimate)")
    else:
        print("sigma_max = inf  (no sigma threshold — use --sigma-max to filter)")

    print(f"RT thresholds — HB: hb ≤ {args.rt_hb_lo} / hb ≥ {args.rt_hb_hi}  |  "
          f"Dist: mindist ≤ {args.rt_dist_lo} / mindist ≥ {args.rt_dist_hi}")

    print(f"Scanning {len(args.dirs)} path(s)...")

    run_metas = discover_runs(args.dirs)
    if not run_metas:
        print("No run directories found. Check paths and naming convention.")
        sys.exit(1)
    print(f"Found {len(run_metas)} run directories.")

    rt_cfg = {
        'rt_hb_lo':   args.rt_hb_lo,
        'rt_hb_hi':   args.rt_hb_hi,
        'rt_dist_lo': args.rt_dist_lo,
        'rt_dist_hi': args.rt_dist_hi,
    }
    for m in run_metas:
        m.update(rt_cfg)

    rows_all = [analyze_run(m) for m in run_metas]
    stitch_continuations(rows_all)

    # Sigma plot always shows every seed — it's a data-quality diagnostic and
    # filtering by RT would hide exactly the seeds you most need to inspect.
    # All other plots (dG, convergence, obs, CSV) respect --min-roundtrips.
    min_rt = args.min_roundtrips
    if min_rt > 0:
        rows = [r for r in rows_all if r.get('n_roundtrips', 0) >= min_rt]
        print(f"RT filter (≥{min_rt} round trips): {len(rows)}/{len(rows_all)} seeds kept")
    else:
        rows = rows_all

    add_cross_seed_stats(rows_all, sigma_max=sigma_max)   # sigma uses all seeds
    add_cross_seed_stats(rows,     sigma_max=sigma_max)   # dG/CSV use filtered

    # ── CSV ──────────────────────────────────────────────────────────────────
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pfx        = args.prefix.rstrip('_')   # e.g. "v8"
    sig_str    = f"sig{args.sigma_max:.3g}" if args.sigma_max < float('inf') else "sigall"
    rt_str     = f"rt{args.min_roundtrips}"
    ts_str     = datetime.now().strftime('%m%d_%H%M')
    subdir_name = '_'.join(filter(None, [pfx, rt_str, sig_str, ts_str]))
    out_dir    = os.path.join(script_dir, 'diagnostics_output', subdir_name)
    os.makedirs(out_dir, exist_ok=True)
    csv_path  = os.path.join(out_dir, f'{subdir_name}.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k, float('nan'))) for k in CSV_FIELDS})
    print(f"Written: {csv_path}  ({len(rows)} rows)\n")

    # ── Plots ─────────────────────────────────────────────────────────────────
    make_sigma_plot(rows_all,   os.path.join(out_dir, f'{subdir_name}_sigma.pdf'))
    make_dg_plot(rows,          os.path.join(out_dir, f'{subdir_name}_dg.pdf'),
                 min_roundtrips=min_rt)
    make_convergence_plot(rows, os.path.join(out_dir, f'{subdir_name}_convergence.pdf'))
    make_obs_plot(rows,         os.path.join(out_dir, f'{subdir_name}_obs.pdf'))

    # ── Console summary table ─────────────────────────────────────────────────
    print()
    hdr  = f"{'System':<10} {'Win':<5} {'Ver':<5} "
    hdr += f"{'PASS':>5} {'FLAG':>5} {'FAIL':>5} {'NO_DATA':>8}  "
    hdr += f"{'outliers':>8}  {'median σ':>9}  {'med RT':>7}  {'filtered ΔG±std (kT)':>22}"
    print(hdr)
    print('─' * len(hdr))

    versions = sorted(set(r['version'] for r in rows))
    for sys_ in SYSTEMS:
        for win in WINDOWS:
            for ver in versions:
                sub = [r for r in rows
                       if r['system'] == sys_ and r['window'] == win
                       and r['version'] == ver]
                if not sub:
                    continue
                counts   = defaultdict(int)
                n_out    = 0
                for r in sub:
                    counts[r['status']] += 1
                    if r.get('seed_outlier', False):
                        n_out += 1
                sigs = [r['sigma_naive_kT'] for r in sub if _valid(r.get('sigma_naive_kT'))]
                med_s = f"{np.median(sigs):.3f}" if sigs else '—'

                rts = [r['n_roundtrips'] for r in sub if r.get('obs_n_points', 0) > 0]
                med_rt = f"{int(np.median(rts))}" if rts else '—'

                # Filtered ΔG: read from first row that has it
                # (same value across the group — just pick the first non-nan)
                fmeans = [r['filtered_mean_dg'] for r in sub if _valid(r.get('filtered_mean_dg'))]
                fstds  = [r['filtered_std_dg']  for r in sub if _valid(r.get('filtered_std_dg'))]
                fns    = [r['filtered_n']        for r in sub if r.get('filtered_n', 0) > 0]
                # Group by combo to show per-combo filtered ΔG
                combo_lines = []
                for combo in COMBOS:
                    csub = [r for r in sub if r['combo'] == combo]
                    if not csub:
                        continue
                    fm = csub[0].get('filtered_mean_dg', float('nan'))
                    fs = csub[0].get('filtered_std_dg',  float('nan'))
                    fn = csub[0].get('filtered_n', 0)
                    if _valid(fm) and _valid(fs) and fn >= 2:
                        combo_lines.append(f"{combo}:{fm:.2f}±{fs:.2f}(n={fn})")
                    elif _valid(fm) and fn >= 1:
                        combo_lines.append(f"{combo}:{fm:.2f}(n={fn})")
                dg_s = '  '.join(combo_lines) if combo_lines else '—'

                print(f"{sys_:<10} {win:<5} {ver:<5} "
                      f"{counts['PASS']:>5} {counts['FLAG']:>5} {counts['FAIL']:>5} "
                      f"{counts['NO_DATA']:>8}  {n_out:>8}  {med_s:>9}  {med_rt:>7}  {dg_s}")

    print("""
Column guide (CSV)
──────────────────
flatness_ratio    u_max / u_rare — ratio of highest to lowest unbiased count
                  across visited states. Pure diagnostic, not used for filtering.
                  1 = perfectly flat histogram (ideal weights). Large = some
                  states sampled much less than others, weights need updating.

sigma_naive_kT    sqrt(1/u_target + 1/u_ref) — lower bound on std(ΔG) [kT].
                  Assumes all histogram counts are independent. True std ≥ this.
                  Useful for relative ranking across seeds, not as an absolute
                  error bar. Use filtered_std_dg for reporting.

seed_std_dg       Std of ΔG across ALL non-FAIL seeds in the same combo group.
                  Includes outliers — shows raw spread before filtering.

seed_outlier      True if |ΔG − seed_mean| > 2 × seed_std.
                  These seeds are excluded from the filtered estimate.

filtered_mean_dg  Mean ΔG after removing outlier seeds. Use this for reporting.
filtered_std_dg   Std of ΔG across non-outlier seeds. The real uncertainty.
filtered_n        Number of seeds contributing to the filtered estimate.
include_in_filter True if this seed is included in the filtered estimate.

traj_block1_dg_kT ΔG computed from the first half of the trajectory (cumulative
                  counts at t_mid).  Based on traj_hist.dat snapshots.
traj_block2_dg_kT ΔG computed from the second half (counts[final] - counts[mid]).
traj_block_diff   |block2 - block1| in kT.  Close to 0 means the two halves agree.
                  A large diff (>1 kT) signals the histogram hasn't converged yet.
                  Use this to rank seeds for continuation: prefer small block_diff.

n_roundtrips      Full round trips counted from the high-resolution observable file
                  (hb_obs.dat or dist_obs.dat, 100-step resolution).
                  HB:   trips between hb ≤ 1 and hb ≥ 5.
                  Dist: trips between mindist ≤ 2 (bin 0) and mindist ≥ 10 (bin 5+).
                  Each count = one complete low→high→low or high→low→high cycle.
                  More round trips = more independent sampling of the free energy.
obs_n_points      Number of rows in the observable file (= n_steps / 100).
""")


if __name__ == '__main__':
    main()


24
1637
761
1493
1864
