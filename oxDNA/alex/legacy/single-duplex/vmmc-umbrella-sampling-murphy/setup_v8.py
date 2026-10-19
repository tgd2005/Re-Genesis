#!/usr/bin/env python3
"""
setup_v8.py  —  create v8 run directories
==========================================

HB window (15 runs per system, 45 total):
  - Continues best v7 HB seed per combo (highest unbiased count at bond=0, mdbin=0)
  - 10 seeds per combo (wA_s1..s10, wB_s1..s10, wC_s1..s10)
  - init_hist.dat = best v7 seed's last_hist.dat
  - Starting conf: best v7 seed's last_conf.dat
  - wfile symlink: same as v7 (wfile_A/B/C.txt from links/6nt/)

Dist window (15 runs per system, 45 total):
  - Continues best v7 dist seed per combo (explored highest mindist bin)
  - 10 seeds per combo (dist_wA_s1..s10, dist_wB_s1..s10, dist_wC_s1..s10)
  - init_hist.dat = best v7 dist seed's last_hist.dat
  - Starting conf: best v7 dist seed's last_conf.dat
  - wfile symlink: same v7 dist wfiles (wfile_v7_dist_{sys}_{X}.txt)

Total v8: 30 HB + 30 dist = 60 runs per system, 180 total.

Usage:
    python3 setup_v8.py           # create all v8 dirs
    python3 setup_v8.py --dry-run # preview only, no files written
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

# HB weight files — same as v7 (frozen by init_hist_file)
HB_WFILE = {
    'wA': 'wfile_A.txt',
    'wB': 'wfile_B.txt',
    'wC': 'wfile_C.txt',
}

# Dist weight files — same as v7
DIST_WFILE = {
    ('pulling',  'wA'): 'wfile_v7_dist_pull_A.txt',
    ('pulling',  'wB'): 'wfile_v7_dist_pull_B.txt',
    ('pulling',  'wC'): 'wfile_v7_dist_pull_C.txt',
    ('noclamps', 'wA'): 'wfile_v7_dist_nocl_A.txt',
    ('noclamps', 'wB'): 'wfile_v7_dist_nocl_B.txt',
    ('noclamps', 'wC'): 'wfile_v7_dist_nocl_C.txt',
    ('pushing',  'wA'): 'wfile_v7_dist_push_A.txt',
    ('pushing',  'wB'): 'wfile_v7_dist_push_B.txt',
    ('pushing',  'wC'): 'wfile_v7_dist_push_C.txt',
}

# v7 HB seeds that were created (system → list of (combo, seed_id))
V7_HB_SEEDS = {
    'pulling':  [('wA','s1'), ('wB','s2'), ('wB','s4'), ('wC','s2'), ('wC','s5')],
    'noclamps': [('wA','s1'), ('wA','s3'), ('wB','s2'), ('wB','s3'), ('wC','s1')],
    'pushing':  [('wA','s1'), ('wA','s3'), ('wB','s3'), ('wB','s5'), ('wC','s3')],
}

# v7 dist seeds that were created (same for all systems)
V7_DIST_SEEDS = ['wA_s1', 'wA_s2', 'wB_s1', 'wB_s2', 'wC_s1']

COMBOS = ['wA', 'wB', 'wC']
N_SEEDS = 10

# Dist starting conf (fallback if no v7 dist seed has usable last_conf)
DIST_START_CONF = {
    'pulling':  'conf_1390000.dat',
    'noclamps': 'conf_6240000.dat',
    'pushing':  'conf_4180000.dat',
}

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


def best_hb_seed(system, combo):
    """Find v7 HB seed with highest unbiased count at (bond=0, mdbin=0).

    Returns (seed_id, hist_dict) or None if no v7 seeds found for this combo.
    """
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    candidates = [(c, s) for (c, s) in V7_HB_SEEDS[system] if c == combo]
    if not candidates:
        return None

    best_id = None
    best_count = -1.0
    best_hist = {}
    for (_, seed_id) in candidates:
        run_dir = os.path.join(sys_path, 'v7', f'{combo}_{seed_id}')
        hist_path = os.path.join(run_dir, 'last_hist.dat')
        hist = read_hist(hist_path)
        count = hist.get((0, 0), 0.0)
        if count > best_count:
            best_count = count
            best_id = seed_id
            best_hist = hist
    return (best_id, best_hist) if best_id is not None else None


def best_dist_seed(system, combo):
    """Find v7 dist seed for this combo that explored the highest mindist bin.

    Criterion: max mdbin with u>0 in last_hist.dat (bond=0 rows only).
    Ties broken by highest u at that bin.
    Returns (seed_id, hist_dict) or None.
    """
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    # v7 dist seeds are: dist_wA_s1, dist_wA_s2, dist_wB_s1, dist_wB_s2, dist_wC_s1
    # Filter to this combo
    candidates = [s for s in V7_DIST_SEEDS if s.startswith(f'{combo}_')]
    if not candidates:
        return None

    best_id = None
    best_max_bin = -1
    best_u_at_max = -1.0
    best_hist = {}
    for seed_str in candidates:
        run_dir = os.path.join(sys_path, 'v7', f'dist_{seed_str}')
        hist_path = os.path.join(run_dir, 'last_hist.dat')
        hist = read_hist(hist_path)
        # Find max mdbin with u > 0 for bond=0
        max_bin = -1
        u_at_max = 0.0
        for (b, mb), u in hist.items():
            if b == 0 and u > 0 and mb > max_bin:
                max_bin = mb
                u_at_max = u
            elif b == 0 and mb == max_bin and u > u_at_max:
                u_at_max = u
        if (max_bin > best_max_bin) or (max_bin == best_max_bin and u_at_max > best_u_at_max):
            best_max_bin = max_bin
            best_u_at_max = u_at_max
            best_id = seed_str
            best_hist = hist
    return (best_id, best_hist) if best_id is not None else None


# ── Setup functions ──────────────────────────────────────────────────────────

def setup_hb(system, combo, new_seed_id, src_seed_id, dry_run=False):
    """Create one v8 HB run directory, continuing from best v7 HB seed."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    run_name = f"{combo}_{new_seed_id}"
    src_dir  = os.path.join(sys_path, 'v7', f'{combo}_{src_seed_id}')
    dst_dir  = os.path.join(sys_path, 'v8', run_name)
    v8_dir   = os.path.join(sys_path, 'v8')

    top_path = _find_topology(sys_path, top_base)

    # Validate sources
    missing = []
    if not os.path.exists(os.path.join(src_dir, 'last_hist.dat')):
        missing.append(f'v7/{combo}_{src_seed_id}/last_hist.dat')
    if not os.path.exists(os.path.join(src_dir, 'last_conf.dat')):
        missing.append(f'v7/{combo}_{src_seed_id}/last_conf.dat')
    if not top_path:
        missing.append(f'{top_base}.top')
    if missing:
        print(f"  SKIP  {run_name}: missing {', '.join(missing)}")
        return False

    print(f"  HB  {SYS_DIR[system]}/v8/{run_name}  ← v7/{combo}_{src_seed_id}")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(top_path, os.path.join(dst_dir, f'{top_base}.top'))
    shutil.copy2(os.path.join(src_dir, 'last_conf.dat'),
                 os.path.join(dst_dir, 'last_conf.dat'))
    shutil.copy2(os.path.join(src_dir, 'last_hist.dat'),
                 os.path.join(dst_dir, 'init_hist.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{HB_WFILE[combo]}')

    _ensure_op_file(v8_dir, sys_path, 'op-setA.txt')
    _symlink(os.path.join(dst_dir, 'op-setA.txt'), '../op-setA.txt')

    sh_path = os.path.join(dst_dir, 'run_VMMC.sh')
    with open(sh_path, 'w') as f:
        f.write(RUN_VMMC_SH)
    os.chmod(sh_path, 0o755)

    seed = random.randint(10000, 99999)
    with open(os.path.join(dst_dir, 'input_VMMC'), 'w') as f:
        f.write(_input_vmmc(top_base, 'op-setA.txt',
                            is_dist=False, seed=seed,
                            system=system, use_init_hist=True))
    return True


def setup_dist(system, combo, new_seed_id, src_seed_str, dry_run=False):
    """Create one v8 dist run directory, continuing from best v7 dist seed."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    run_name = f"dist_{combo}_{new_seed_id}"
    src_dir  = os.path.join(sys_path, 'v7', f'dist_{src_seed_str}')
    dst_dir  = os.path.join(sys_path, 'v8', run_name)
    v8_dir   = os.path.join(sys_path, 'v8')

    top_path = _find_topology(sys_path, top_base)

    # Validate sources
    missing = []
    if not os.path.exists(os.path.join(src_dir, 'last_hist.dat')):
        missing.append(f'v7/dist_{src_seed_str}/last_hist.dat')
    if not os.path.exists(os.path.join(src_dir, 'last_conf.dat')):
        missing.append(f'v7/dist_{src_seed_str}/last_conf.dat')
    if not top_path:
        missing.append(f'{top_base}.top')
    if missing:
        print(f"  SKIP  {run_name}: missing {', '.join(missing)}")
        return False

    print(f"  dist {SYS_DIR[system]}/v8/{run_name}  ← v7/dist_{src_seed_str}")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(top_path, os.path.join(dst_dir, f'{top_base}.top'))
    shutil.copy2(os.path.join(src_dir, 'last_conf.dat'),
                 os.path.join(dst_dir, 'last_conf.dat'))
    shutil.copy2(os.path.join(src_dir, 'last_hist.dat'),
                 os.path.join(dst_dir, 'init_hist.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{DIST_WFILE[(system, combo)]}')

    _ensure_op_file(v8_dir, sys_path, 'op-setA_dist.txt')
    _symlink(os.path.join(dst_dir, 'op-setA_dist.txt'), '../op-setA_dist.txt')

    sh_path = os.path.join(dst_dir, 'run_VMMC.sh')
    with open(sh_path, 'w') as f:
        f.write(RUN_VMMC_SH)
    os.chmod(sh_path, 0o755)

    seed = random.randint(10000, 99999)
    with open(os.path.join(dst_dir, 'input_VMMC'), 'w') as f:
        f.write(_input_vmmc(top_base, 'op-setA_dist.txt',
                            is_dist=True, seed=seed,
                            system=system, use_init_hist=True))
    return True


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Set up v8 runs (HB + dist continuation from best v7 seeds)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--dry-run', action='store_true',
                    help='Preview what would be created, write nothing')
    args = ap.parse_args()

    random.seed()

    if args.dry_run:
        print("DRY RUN — no files will be written\n")

    n_ok = n_skip = 0

    print("=" * 60)
    print("HB WINDOW  (init_hist + last_conf from best v7 HB seed)")
    print("=" * 60)
    for system in ['pulling', 'noclamps', 'pushing']:
        print(f"\n{system.upper()}")
        for combo in COMBOS:
            result = best_hb_seed(system, combo)
            if result is None:
                print(f"  SKIP  {combo}: no v7 HB seeds found for this combo")
                n_skip += N_SEEDS
                continue
            src_seed_id, hist = result
            u_target = hist.get((0, 0), 0.0)
            print(f"  Best v7 seed: {combo}_{src_seed_id}  (u at bond=0,mdbin=0 = {u_target:.1f})")
            for i in range(1, N_SEEDS + 1):
                ok = setup_hb(system, combo, f's{i}', src_seed_id,
                              dry_run=args.dry_run)
                n_ok += ok; n_skip += (not ok)

    print()
    print("=" * 60)
    print("DIST WINDOW  (init_hist + last_conf from best v7 dist seed)")
    print("=" * 60)
    for system in ['pulling', 'noclamps', 'pushing']:
        print(f"\n{system.upper()}")
        for combo in COMBOS:
            result = best_dist_seed(system, combo)
            if result is None:
                print(f"  SKIP  {combo}: no v7 dist seeds found for this combo")
                n_skip += N_SEEDS
                continue
            src_seed_str, hist = result
            max_bin = max((mb for (b, mb), u in hist.items()
                           if b == 0 and u > 0), default=-1)
            print(f"  Best v7 seed: dist_{src_seed_str}  (max mdbin explored = {max_bin})")
            for i in range(1, N_SEEDS + 1):
                ok = setup_dist(system, combo, f's{i}', src_seed_str,
                                dry_run=args.dry_run)
                n_ok += ok; n_skip += (not ok)

    print(f"\n{'Would create' if args.dry_run else 'Created'}: {n_ok}   "
          f"Skipped: {n_skip}\n")

    if not args.dry_run and n_ok > 0:
        print("Next steps:")
        print("  1. Spot-check a few input_VMMC and init_hist.dat files")
        print("  2. Push to cluster:")
        print("       pushnew oxview-assemblies/pulling-clamps/6nt-clamps-pulling/v8/")
        print("       pushnew oxview-assemblies/noclamps/6nt-noclamps/v8/")
        print("       pushnew oxview-assemblies/pushing-clamps/6nt-pushing/v8/")
        print("       pushnew oxview-assemblies/links/6nt/")
        print("  3. Submit on cluster:")
        print("       for d in */*/v8/*_s*/; do")
        print("         cd \"$REMOTE_SYNC_ROOT/$d\" && qsub run_VMMC.sh && cd \"$REMOTE_SYNC_ROOT\"")
        print("       done")


if __name__ == '__main__':
    main()
