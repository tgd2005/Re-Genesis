#!/usr/bin/env python3
"""
setup_v7.py  —  create v7 run directories
==========================================

HB window (15 runs per system, 45 total):
  - Same 5 best v4 seeds per system as v6
  - init_hist_file carried from v4 (weights UNCHANGED)
  - Starting conf: fully-bonded relaxed dat (bond=6)
  - Observable: hb_list count at 100-step resolution

Dist window (15 runs per system, 45 total):
  - Fresh start — NO init_hist_file (new weight files, bin 6 now allowed)
  - New wfile_v7_dist_*.txt: per-bin weights derived from v5 data, bin 6 = 1
    (bin 6 previously blocked by default_weight=0 — now allowed so system
     can visit mindist>12 if geometrically accessible; expect it to stay ~0)
  - Starting conf: bond=0, dist_bin=0 conf (strands close, fully separated)
  - Observable: distance between ALL 6 stem pairs (7 columns) so post-processing
    can compute the true mindistance (min over all 6), not just one pair

Weight file design:
  - interfaces = 2., 4., 6., 8., 10., 12. → 7 bins (0-2, 2-4, ..., ≥12)
  - Bin 5 (10-12) = reference state (weight=1); bin 6 (≥12) also weight=1
  - A/B/C variants: bins 0,1,2 × 1 / ×2 / ×0.5 (bins 3,4,5,6 fixed)

Usage:
    python3 setup_v7.py           # create all v7 dirs
    python3 setup_v7.py --dry-run # preview only, no files written
"""

import os, shutil, random, argparse, glob

ROOT = os.path.dirname(os.path.abspath(__file__))

# ── HB seeds: same best v4 seeds as v6 ───────────────────────────────────────
# (system, 'HB', 'v4', combo, seed_id)
HB_SEEDS = [
    ('pulling',  'wA', 's1'),
    ('pulling',  'wB', 's2'),
    ('pulling',  'wB', 's4'),
    ('pulling',  'wC', 's2'),
    ('pulling',  'wC', 's5'),
    ('noclamps', 'wB', 's3'),
    ('noclamps', 'wB', 's2'),
    ('noclamps', 'wA', 's1'),
    ('noclamps', 'wA', 's3'),
    ('noclamps', 'wC', 's1'),
    ('pushing',  'wA', 's1'),
    ('pushing',  'wA', 's3'),
    ('pushing',  'wB', 's3'),
    ('pushing',  'wB', 's5'),
    ('pushing',  'wC', 's3'),
]

# ── Dist seeds: fresh runs, 5 per system, distributed across A/B/C ───────────
# (system, combo, seed_id)  — no source version, always fresh
DIST_SEEDS = [
    ('pulling',  'wA', 's1'),
    ('pulling',  'wA', 's2'),
    ('pulling',  'wB', 's1'),
    ('pulling',  'wB', 's2'),
    ('pulling',  'wC', 's1'),
    ('noclamps', 'wA', 's1'),
    ('noclamps', 'wA', 's2'),
    ('noclamps', 'wB', 's1'),
    ('noclamps', 'wB', 's2'),
    ('noclamps', 'wC', 's1'),
    ('pushing',  'wA', 's1'),
    ('pushing',  'wA', 's2'),
    ('pushing',  'wB', 's1'),
    ('pushing',  'wB', 's2'),
    ('pushing',  'wC', 's1'),
]

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

# HB weight files — must match v4 source (frozen by init_hist_file)
HB_WFILE = {
    'wA': 'wfile_A.txt',
    'wB': 'wfile_B.txt',
    'wC': 'wfile_C.txt',
}

# Dist weight files — new v7 files with per-bin weights + bin 6 allowed
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

# All 6 stem pairs per system (v4/v5 topology indices, from CLAUDE.md)
ALL_PAIRS = {
    'pulling':  [(15,56), (16,55), (17,54), (18,53), (19,52), (20,51)],
    'noclamps': [(54,95), (55,94), (56,93), (57,92), (58,91), (59,90)],
    'pushing':  [(99,140),(100,139),(101,138),(102,137),(103,136),(104,135)],
}

# HB starting conf: fully-bonded relaxed structure
HB_START_CONF = {
    'pulling':  ('v4/wA_s1', '6nt-pulling-relaxed-NEW.dat'),
    'noclamps': ('.',         '6nt-noclamps-relaxed-NEW.dat'),
    'pushing':  ('.',         '6nt-pushing-relaxed-NEW.dat'),
}

# Dist starting conf: bond=0, mindist~0 (strands fully separated but close in box)
DIST_START_CONF = {
    'pulling':  ('.',  'conf_1390000.dat'),
    'noclamps': ('.',  'conf_6240000.dat'),
    'pushing':  ('.',  'conf_4180000.dat'),
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


# ── Core setup functions ─────────────────────────────────────────────────────

def _symlink(link_path, target):
    if os.path.lexists(link_path):
        os.remove(link_path)
    os.symlink(target, link_path)


def _find_topology(sys_path, top_base):
    """Find topology file in v5 or v4 (any run dir)."""
    for ver in ['v5', 'v4']:
        pattern = os.path.join(sys_path, ver, '*', f'{top_base}.top')
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
    return None


def _ensure_op_file(v7_dir, sys_path, op_fname):
    """Ensure real op file exists in v7/, copy from v5/v4 if needed."""
    v7_op = os.path.join(v7_dir, op_fname)
    if not os.path.exists(v7_op):
        for ver in ['v5', 'v4']:
            candidate = os.path.realpath(
                os.path.join(sys_path, ver, op_fname))
            if os.path.isfile(candidate):
                shutil.copy2(candidate, v7_op)
                return True
        print(f"  WARNING: could not find {op_fname} to copy into v7/")
        return False
    return True


def setup_hb(system, combo, seed_id, dry_run=False):
    """Create one v7 HB run directory, continuing from v4."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    run_name = f"{combo}_{seed_id}"
    src_dir  = os.path.join(sys_path, 'v4', run_name)
    dst_dir  = os.path.join(sys_path, 'v7', run_name)
    v7_dir   = os.path.join(sys_path, 'v7')

    conf_subdir, conf_fname = HB_START_CONF[system]
    start_conf_path = os.path.join(sys_path, conf_subdir, conf_fname)

    # Validate sources
    for fpath, label in [
        (src_dir,                                        f'v4/{run_name}'),
        (os.path.join(src_dir, 'last_hist.dat'),         'last_hist.dat'),
        (os.path.join(src_dir, f'{top_base}.top'),       f'{top_base}.top'),
        (start_conf_path,                                conf_fname),
    ]:
        if not os.path.exists(fpath):
            print(f"  SKIP  missing {label}: {fpath}")
            return False

    print(f"  HB  {SYS_DIR[system]}/v7/{run_name}  ← v4/{run_name}  [conf: {conf_fname}]")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(os.path.join(src_dir, f'{top_base}.top'), dst_dir)
    shutil.copy2(start_conf_path, os.path.join(dst_dir, 'last_conf.dat'))
    shutil.copy2(os.path.join(src_dir, 'last_hist.dat'),
                 os.path.join(dst_dir, 'init_hist.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{HB_WFILE[combo]}')

    _ensure_op_file(v7_dir, sys_path, 'op-setA.txt')
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


def setup_dist(system, combo, seed_id, dry_run=False):
    """Create one v7 dist run directory — fresh start, no init_hist."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    run_name = f"dist_{combo}_{seed_id}"
    dst_dir  = os.path.join(sys_path, 'v7', run_name)
    v7_dir   = os.path.join(sys_path, 'v7')

    conf_subdir, conf_fname = DIST_START_CONF[system]
    start_conf_path = os.path.join(sys_path, conf_subdir, conf_fname)

    top_path = _find_topology(sys_path, top_base)

    # Validate sources
    for fpath, label in [
        (start_conf_path, conf_fname),
        (top_path or '',  f'{top_base}.top'),
    ]:
        if not fpath or not os.path.exists(fpath):
            print(f"  SKIP  missing {label}: {fpath}")
            return False

    print(f"  dist {SYS_DIR[system]}/v7/{run_name}  [fresh, conf: {conf_fname}]")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(top_path, os.path.join(dst_dir, f'{top_base}.top'))
    shutil.copy2(start_conf_path, os.path.join(dst_dir, 'last_conf.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{DIST_WFILE[(system, combo)]}')

    _ensure_op_file(v7_dir, sys_path, 'op-setA_dist.txt')
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
        description='Set up v7 runs (HB continuation + fresh dist)',
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
    cur_system = None

    print("=" * 60)
    print("HB WINDOW  (init_hist from v4, start from relaxed conf)")
    print("=" * 60)
    for system, combo, seed_id in HB_SEEDS:
        if system != cur_system:
            cur_system = system
            print(f"\n{system.upper()}")
        ok = setup_hb(system, combo, seed_id, dry_run=args.dry_run)
        n_ok += ok; n_skip += (not ok)

    print()
    print("=" * 60)
    print("DIST WINDOW  (fresh start, new v7 weights, all-pairs obs)")
    print("=" * 60)
    cur_system = None
    for system, combo, seed_id in DIST_SEEDS:
        if system != cur_system:
            cur_system = system
            print(f"\n{system.upper()}")
        ok = setup_dist(system, combo, seed_id, dry_run=args.dry_run)
        n_ok += ok; n_skip += (not ok)

    print(f"\n{'Would create' if args.dry_run else 'Created'}: {n_ok}   "
          f"Skipped: {n_skip}\n")

    if not args.dry_run and n_ok > 0:
        print("Next steps:")
        print("  1. Spot-check a few input_VMMC and links/6nt/wfile_v7_dist_*.txt")
        print("  2. Push to cluster:")
        print("       pushnew oxview-assemblies/pulling-clamps/6nt-clamps-pulling/v7/")
        print("       pushnew oxview-assemblies/noclamps/6nt-noclamps/v7/")
        print("       pushnew oxview-assemblies/pushing-clamps/6nt-pushing/v7/")
        print("       pushnew oxview-assemblies/links/6nt/")
        print("  3. Submit on cluster:")
        print("       for d in */*/v7/*_s*/; do")
        print("         cd \"$REMOTE_SYNC_ROOT/$d\" && qsub run_VMMC.sh && cd \"$REMOTE_SYNC_ROOT\"")
        print("       done")


if __name__ == '__main__':
    main()
