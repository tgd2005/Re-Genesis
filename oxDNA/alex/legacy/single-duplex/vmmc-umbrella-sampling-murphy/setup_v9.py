#!/usr/bin/env python3
"""
setup_v9.py  —  create v9 run directories with updated weights
==============================================================

Iterative weight update derived from pooled v7 data. Fresh start (no init_hist).
10 seeds per combo (wA_s1..s10, wB_s1..s10, wC_s1..s10) per system per window.

HB window (30 per system, 90 total):
  - New system-specific weights from v7 pooled unbiased counts
  - A: ideal  B: bins 0,1,2 ×2  C: bins 0,1,2 ×0.5
  - NO init_hist_file
  - Starting conf: fully-bonded relaxed structure (bond=6)

Dist window (30 per system, 90 total):
  - New system-specific dist weights from v7 pooled unbiased counts
  - Same A/B/C variant scheme
  - NO init_hist_file
  - Starting conf: bond=0 conf (strands separated, close in box)

Weight computation:
  - Pool all v7 seeds' last_hist.dat for each system+window
  - HB:   w_ideal(b) = u_pooled(b=6) / u_pooled(b)   [reference = bond=6]
  - Dist: w_ideal(mb) = u_max / u_pooled(mb)           [reference = most-visited bin]
  - Zero-count fallback: 10 × max(w_ideal for nonzero bins)

Wfiles written to links/6nt/:
  - wfile_v9_hb_{pull/nocl/push}_{A/B/C}.txt
  - wfile_v9_dist_{pull/nocl/push}_{A/B/C}.txt

Usage:
    python3 setup_v9.py           # compute weights, create all v9 dirs
    python3 setup_v9.py --dry-run # preview weights and what would be created
"""

import os, shutil, random, argparse, glob

ROOT = os.path.dirname(os.path.abspath(__file__))

# ── Static lookup tables ─────────────────────────────────────────────────────

SYS_DIR = {
    'pulling':  'pulling-clamps/6nt-clamps-pulling',
    'noclamps': 'noclamps/6nt-noclamps',
    'pushing':  'pushing-clamps/6nt-pushing',
}

SYS_TOP = {
    'pulling':  '6nt-pulling-relaxed-NEW',
    'noclamps': '6nt-noclamps-relaxed-NEW',
    'pushing':  '6nt-pushing-relaxed-NEW',
}

ALL_PAIRS = {
    'pulling':  [(15,56), (16,55), (17,54), (18,53), (19,52), (20,51)],
    'noclamps': [(54,95), (55,94), (56,93), (57,92), (58,91), (59,90)],
    'pushing':  [(99,140),(100,139),(101,138),(102,137),(103,136),(104,135)],
}

DIST_START_CONF = {
    'pulling':  'conf_1390000.dat',
    'noclamps': 'conf_6240000.dat',
    'pushing':  'conf_4180000.dat',
}

# v7 HB seeds per system (combo, seed_id)
V7_HB_SEEDS = {
    'pulling':  [('wA','s1'), ('wB','s2'), ('wB','s4'), ('wC','s2'), ('wC','s5')],
    'noclamps': [('wA','s1'), ('wA','s3'), ('wB','s2'), ('wB','s3'), ('wC','s1')],
    'pushing':  [('wA','s1'), ('wA','s3'), ('wB','s3'), ('wB','s5'), ('wC','s3')],
}

# v7 dist seeds (same for all systems): dist_{seed_str}
V7_DIST_SEEDS = ['wA_s1', 'wA_s2', 'wB_s1', 'wB_s2', 'wC_s1']

SYS_ABBREV = {'pulling': 'pull', 'noclamps': 'nocl', 'pushing': 'push'}

COMBOS  = ['wA', 'wB', 'wC']
N_SEEDS = 10

# ── File content templates ───────────────────────────────────────────────────

RUN_VMMC_SH = """\
#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd "$PBS_O_WORKDIR"
inputfile="input_VMMC"

runcommand="$HOME/oxDNA_git/oxDNA/build/bin/oxDNA  $inputfile"
timeout 71h $runcommand
"""


def _obs_hb():
    """HB window observable: hb_list count at 100-step resolution."""
    return (
        "data_output_1 = {\n"
        "  name = hb_obs.dat\n"
        "  print_every = 100\n"
        "  col_1 = {\n"
        "    type = step\n"
        "  }\n"
        "  col_2 = {\n"
        "    type = hb_list\n"
        "    only_count = true\n"
        "    order_parameters_file = op-setA.txt\n"
        "  }\n"
        "}"
    )


def _obs_dist(system):
    """Dist window observable: step + distance for all 6 stem pairs.

    Outputs 7 columns: step, d(pair1), d(pair2), ..., d(pair6).
    Post-process with min(col2..col7) to get the true mindistance,
    matching what the mindistance OP computes.
    """
    pairs = ALL_PAIRS[system]
    lines = [
        "data_output_1 = {",
        "  name = dist_obs.dat",
        "  print_every = 100",
        "  col_1 = {",
        "    type = step",
        "  }",
    ]
    for i, (p1, p2) in enumerate(pairs, start=2):
        lines += [
            f"  col_{i} = {{",
            f"    type = distance",
            f"    particle_1 = {p1}",
            f"    particle_2 = {p2}",
            f"    PBC = true",
            f"  }}",
        ]
    lines.append("}")
    return "\n".join(lines)


def _input_vmmc(top_base, op_file, is_dist, seed, system, use_init_hist):
    obs = _obs_dist(system) if is_dist else _obs_hb()
    init_hist_line = "init_hist_file = init_hist.dat\n" if use_init_hist else ""
    return (
        f"####  PROGRAM PARAMETERS  ####\n"
        f"backend = CPU\n"
        f"seed = {seed}\n"
        f"interaction_type = DNA2\n"
        f"####    SIM PARAMETERS    ####\n"
        f"sim_type = VMMC\n"
        f"ensemble = NVT\n"
        f"steps = 1e10\n"
        f"maxclust = 50\n"
        f"delta_translation = 0.25\n"
        f"delta_rotation = 0.15\n"
        f"T = 20C\n"
        f"verlet_skin = 1.25\n"
        f"umbrella_sampling = 1\n"
        f"op_file = {op_file}\n"
        f"weights_file = wfile.txt\n"
        f"salt_concentration = 0.5\n"
        f"max_backbone_force = 5\n"
        f"max_backbone_force_far = 10\n"
        f"### Seq. Dependence ##########\n"
        f"use_average_seq = 1\n"
        f"small_system = 1\n"
        f"####    INPUT / OUTPUT    ####\n"
        f"topology = {top_base}.top\n"
        f"conf_file = last_conf.dat\n"
        f"trajectory_file = trajectory.dat\n"
        f"no_stdout_energy = 0\n"
        f"restart_step_counter = 1\n"
        f"energy_file = energy.dat\n"
        f"print_conf_interval = 1e4\n"
        f"print_energy_every = 1e4\n"
        f"time_scale = linear\n"
        f"default_weight = 0.\n"
        f"safe_weights = 0\n"
        f"{init_hist_line}"
        f"{obs}\n"
    )


# ── Core helpers ─────────────────────────────────────────────────────────────

def _symlink(link_path, target):
    if os.path.lexists(link_path):
        os.remove(link_path)
    os.symlink(target, link_path)


def _find_topology(sys_path, top_base):
    """Find topology file in v7, v5, or v4 (any run dir)."""
    for ver in ['v7', 'v5', 'v4']:
        pattern = os.path.join(sys_path, ver, '*', f'{top_base}.top')
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
    return None


def _ensure_op_file(vN_dir, sys_path, op_fname):
    """Ensure real op file exists in vN_dir, copy from v7/v5/v4 if needed."""
    vN_op = os.path.join(vN_dir, op_fname)
    if not os.path.exists(vN_op):
        for ver in ['v7', 'v5', 'v4']:
            candidate = os.path.realpath(
                os.path.join(sys_path, ver, op_fname))
            if os.path.isfile(candidate):
                shutil.copy2(candidate, vN_op)
                return True
        print(f"  WARNING: could not find {op_fname} to copy into {os.path.basename(vN_dir)}/")
        return False
    return True


def read_hist(path):
    """Parse last_hist.dat → dict {(bond, mdbin): unbiased_count}."""
    data = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                p = line.split()
                if len(p) < 4:
                    continue
                data[(int(p[0]), int(p[1]))] = float(p[3])
    except (OSError, IOError):
        pass
    return data


def find_relaxed_conf(system):
    """Find the original relaxed bonded conf (bond=6 starting structure).

    Checks:
      - pulling:  v4/wA_s1/{top_base}.dat
      - noclamps: sys_path/{top_base}.dat
      - pushing:  sys_path/{top_base}.dat
    Falls back to any v4 run dir containing the dat file.
    """
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    dat_fname = f'{top_base}.dat'

    # Primary locations (mirrors v7 HB_START_CONF logic)
    if system == 'pulling':
        primary = os.path.join(sys_path, 'v4', 'wA_s1', dat_fname)
    else:
        primary = os.path.join(sys_path, dat_fname)

    if os.path.exists(primary):
        return primary

    # Fallback: search v7, v5, v4 run dirs
    for ver in ['v7', 'v5', 'v4']:
        pattern = os.path.join(sys_path, ver, '*', dat_fname)
        matches = glob.glob(pattern)
        if matches:
            return matches[0]

    return None


# ── Weight computation ────────────────────────────────────────────────────────

def pool_v7_hists(system, window):
    """Pool all v7 seeds' last_hist.dat for system+window.

    window = 'HB' or 'dist'.
    Returns dict {(bond, mdbin): summed_unbiased_count}.
    """
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    pooled = {}

    if window == 'HB':
        seeds = [(combo, sid) for (combo, sid) in V7_HB_SEEDS[system]]
        run_dirs = [os.path.join(sys_path, 'v7', f'{c}_{s}')
                    for (c, s) in seeds]
    else:
        run_dirs = [os.path.join(sys_path, 'v7', f'dist_{s}')
                    for s in V7_DIST_SEEDS]

    for run_dir in run_dirs:
        hist_path = os.path.join(run_dir, 'last_hist.dat')
        hist = read_hist(hist_path)
        for key, u in hist.items():
            pooled[key] = pooled.get(key, 0.0) + u

    return pooled


def compute_hb_weights(system):
    """Compute ideal HB weights for system from pooled v7 HB data.

    w(b) = u_pooled(b=6, mdbin=0) / u_pooled(b, mdbin=0)
    Returns list of 7 floats [w0..w6].
    """
    pooled = pool_v7_hists(system, 'HB')
    u_ref = pooled.get((6, 0), 0.0)
    if u_ref == 0.0:
        print(f"  WARNING: {system} HB pooled u at bond=6 is zero; using 1.0 as reference")
        u_ref = 1.0

    weights_raw = []
    for b in range(7):
        u = pooled.get((b, 0), 0.0)
        if u > 0:
            weights_raw.append(u_ref / u)
        else:
            weights_raw.append(None)  # placeholder

    # Zero-count fallback: 10 × max of nonzero ideal weights
    nonzero = [w for w in weights_raw if w is not None]
    fallback = 10.0 * max(nonzero) if nonzero else 1.0
    weights = [w if w is not None else fallback for w in weights_raw]
    return weights


def compute_dist_weights(system):
    """Compute ideal dist weights for system from pooled v7 dist data.

    w(mb) = u_max / u_pooled(bond=0, mdbin=mb)
    where u_max = max(u_pooled for mb in 0..6 with u>0).
    Returns list of 7 floats [w0..w6].
    """
    pooled = pool_v7_hists(system, 'dist')
    u_vals = [pooled.get((0, mb), 0.0) for mb in range(7)]
    u_max = max(u for u in u_vals if u > 0) if any(u > 0 for u in u_vals) else 1.0

    weights_raw = []
    for mb in range(7):
        u = u_vals[mb]
        if u > 0:
            weights_raw.append(u_max / u)
        else:
            weights_raw.append(None)

    nonzero = [w for w in weights_raw if w is not None]
    fallback = 10.0 * max(nonzero) if nonzero else 1.0
    weights = [w if w is not None else fallback for w in weights_raw]
    return weights


def make_variants(weights_a):
    """Return A/B/C variants: A=ideal, B=bins0-2 ×2, C=bins0-2 ×0.5."""
    w_b = list(weights_a)
    w_c = list(weights_a)
    for i in range(3):
        w_b[i] = weights_a[i] * 2.0
        w_c[i] = weights_a[i] * 0.5
    return {'wA': weights_a, 'wB': w_b, 'wC': w_c}


def write_hb_wfile(path, weights):
    """Write HB wfile: 7 rows, bond=0..6, mdbin=0."""
    with open(path, 'w') as f:
        for b, w in enumerate(weights):
            f.write(f"{b} 0  {w:.6g}\n")


def write_dist_wfile(path, weights):
    """Write dist wfile: 7 rows, bond=0, mdbin=0..6."""
    with open(path, 'w') as f:
        for mb, w in enumerate(weights):
            f.write(f"0 {mb}  {w:.6g}\n")


def write_all_wfiles(dry_run=False):
    """Compute and write all v9 wfiles to links/6nt/. Returns nested dict of weights."""
    links_dir = os.path.join(ROOT, 'links', '6nt')
    all_weights = {}  # system -> {'hb': {X: [...]}, 'dist': {X: [...]}}

    print("=" * 60)
    print("WEIGHT COMPUTATION  (pooled from v7 data)")
    print("=" * 60)

    for system in ['pulling', 'noclamps', 'pushing']:
        abbrev = SYS_ABBREV[system]
        all_weights[system] = {'hb': {}, 'dist': {}}

        # HB weights
        hb_a = compute_hb_weights(system)
        hb_variants = make_variants(hb_a)
        all_weights[system]['hb'] = hb_variants

        print(f"\n{system.upper()} HB weights:")
        print(f"  {'Bin':>4}  {'A':>12}  {'B':>12}  {'C':>12}")
        for b in range(7):
            print(f"  b={b}   {hb_variants['wA'][b]:>12.4g}  "
                  f"{hb_variants['wB'][b]:>12.4g}  "
                  f"{hb_variants['wC'][b]:>12.4g}")

        # Dist weights
        dist_a = compute_dist_weights(system)
        dist_variants = make_variants(dist_a)
        all_weights[system]['dist'] = dist_variants

        print(f"\n{system.upper()} dist weights:")
        print(f"  {'Bin':>4}  {'A':>12}  {'B':>12}  {'C':>12}")
        for mb in range(7):
            print(f"  mb={mb}  {dist_variants['wA'][mb]:>12.4g}  "
                  f"{dist_variants['wB'][mb]:>12.4g}  "
                  f"{dist_variants['wC'][mb]:>12.4g}")

        if not dry_run:
            os.makedirs(links_dir, exist_ok=True)
            for X in ['wA', 'wB', 'wC']:
                letter = X[1]  # 'A', 'B', or 'C'
                hb_fname  = f'wfile_v9_hb_{abbrev}_{letter}.txt'
                dist_fname = f'wfile_v9_dist_{abbrev}_{letter}.txt'
                write_hb_wfile(os.path.join(links_dir, hb_fname),
                               hb_variants[X])
                write_dist_wfile(os.path.join(links_dir, dist_fname),
                                 dist_variants[X])
                print(f"  Wrote {hb_fname}")
                print(f"  Wrote {dist_fname}")

    return all_weights


# ── Setup functions ──────────────────────────────────────────────────────────

def setup_hb(system, combo, seed_id, dry_run=False):
    """Create one v9 HB run directory — fresh start, new weights."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    abbrev   = SYS_ABBREV[system]
    letter   = combo[1]  # 'A', 'B', 'C'
    run_name = f"{combo}_{seed_id}"
    dst_dir  = os.path.join(sys_path, 'v9', run_name)
    v9_dir   = os.path.join(sys_path, 'v9')

    relaxed_conf = find_relaxed_conf(system)
    top_path     = _find_topology(sys_path, top_base)

    missing = []
    if not relaxed_conf:
        missing.append(f'{top_base}.dat (relaxed bonded conf)')
    if not top_path:
        missing.append(f'{top_base}.top')
    if missing:
        print(f"  SKIP  {run_name}: missing {', '.join(missing)}")
        return False

    wfile_name = f'wfile_v9_hb_{abbrev}_{letter}.txt'
    print(f"  HB  {SYS_DIR[system]}/v9/{run_name}  [fresh, conf: {os.path.basename(relaxed_conf)}]")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(top_path, os.path.join(dst_dir, f'{top_base}.top'))
    shutil.copy2(relaxed_conf, os.path.join(dst_dir, 'last_conf.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{wfile_name}')

    _ensure_op_file(v9_dir, sys_path, 'op-setA.txt')
    _symlink(os.path.join(dst_dir, 'op-setA.txt'), '../op-setA.txt')

    sh_path = os.path.join(dst_dir, 'run_VMMC.sh')
    with open(sh_path, 'w') as f:
        f.write(RUN_VMMC_SH)
    os.chmod(sh_path, 0o755)

    seed = random.randint(10000, 99999)
    with open(os.path.join(dst_dir, 'input_VMMC'), 'w') as f:
        f.write(_input_vmmc(top_base, 'op-setA.txt',
                            is_dist=False, seed=seed,
                            system=system, use_init_hist=False))
    return True


def setup_dist(system, combo, seed_id, dry_run=False):
    """Create one v9 dist run directory — fresh start, new weights."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    abbrev   = SYS_ABBREV[system]
    letter   = combo[1]
    run_name = f"dist_{combo}_{seed_id}"
    dst_dir  = os.path.join(sys_path, 'v9', run_name)
    v9_dir   = os.path.join(sys_path, 'v9')

    conf_fname       = DIST_START_CONF[system]
    start_conf_path  = os.path.join(sys_path, conf_fname)
    top_path         = _find_topology(sys_path, top_base)

    missing = []
    if not os.path.exists(start_conf_path):
        missing.append(conf_fname)
    if not top_path:
        missing.append(f'{top_base}.top')
    if missing:
        print(f"  SKIP  {run_name}: missing {', '.join(missing)}")
        return False

    wfile_name = f'wfile_v9_dist_{abbrev}_{letter}.txt'
    print(f"  dist {SYS_DIR[system]}/v9/{run_name}  [fresh, conf: {conf_fname}]")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(top_path, os.path.join(dst_dir, f'{top_base}.top'))
    shutil.copy2(start_conf_path, os.path.join(dst_dir, 'last_conf.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{wfile_name}')

    _ensure_op_file(v9_dir, sys_path, 'op-setA_dist.txt')
    _symlink(os.path.join(dst_dir, 'op-setA_dist.txt'), '../op-setA_dist.txt')

    sh_path = os.path.join(dst_dir, 'run_VMMC.sh')
    with open(sh_path, 'w') as f:
        f.write(RUN_VMMC_SH)
    os.chmod(sh_path, 0o755)

    seed = random.randint(10000, 99999)
    with open(os.path.join(dst_dir, 'input_VMMC'), 'w') as f:
        f.write(_input_vmmc(top_base, 'op-setA_dist.txt',
                            is_dist=True, seed=seed,
                            system=system, use_init_hist=False))
    return True


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Set up v9 runs (fresh start, updated weights from v7)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--dry-run', action='store_true',
                    help='Preview weights and what would be created, write nothing')
    args = ap.parse_args()

    random.seed()

    if args.dry_run:
        print("DRY RUN — no files will be written\n")

    # Compute and (optionally) write weight files
    write_all_wfiles(dry_run=args.dry_run)

    n_ok = n_skip = 0

    print()
    print("=" * 60)
    print("HB WINDOW  (fresh start, new v9 weights, no init_hist)")
    print("=" * 60)
    for system in ['pulling', 'noclamps', 'pushing']:
        print(f"\n{system.upper()}")
        for combo in COMBOS:
            for i in range(1, N_SEEDS + 1):
                ok = setup_hb(system, combo, f's{i}', dry_run=args.dry_run)
                n_ok += ok; n_skip += (not ok)

    print()
    print("=" * 60)
    print("DIST WINDOW  (fresh start, new v9 weights, no init_hist)")
    print("=" * 60)
    for system in ['pulling', 'noclamps', 'pushing']:
        print(f"\n{system.upper()}")
        for combo in COMBOS:
            for i in range(1, N_SEEDS + 1):
                ok = setup_dist(system, combo, f's{i}', dry_run=args.dry_run)
                n_ok += ok; n_skip += (not ok)

    print(f"\n{'Would create' if args.dry_run else 'Created'}: {n_ok}   "
          f"Skipped: {n_skip}\n")

    if not args.dry_run and n_ok > 0:
        print("Next steps:")
        print("  1. Spot-check a few input_VMMC and links/6nt/wfile_v9_*.txt")
        print("  2. Push to cluster:")
        print("       pushnew oxview-assemblies/pulling-clamps/6nt-clamps-pulling/v9/")
        print("       pushnew oxview-assemblies/noclamps/6nt-noclamps/v9/")
        print("       pushnew oxview-assemblies/pushing-clamps/6nt-pushing/v9/")
        print("       pushnew oxview-assemblies/links/6nt/")
        print("  3. Submit on cluster:")
        print("       for d in */*/v9/*_s*/; do")
        print("         cd \"$REMOTE_SYNC_ROOT/$d\" && qsub run_VMMC.sh && cd \"$REMOTE_SYNC_ROOT\"")
        print("       done")


if __name__ == '__main__':
    main()
