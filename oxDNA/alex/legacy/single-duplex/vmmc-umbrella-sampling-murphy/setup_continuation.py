#!/usr/bin/env python3
"""
setup_continuation.py  —  create v6 continuation run directories
=================================================================

HB seeds come from v4 (all-PASS, ~13-20M steps, zero FAILs across all systems,
shared wfile_A/B/C.txt weights).
Dist seeds come from v5 (v4 dist never converged; system-specific v5 weights).

Each v6 run directory contains:
  input_VMMC        new random seed + init_hist_file + inline observable
  run_VMMC.sh       identical PBS script (72h, 1CPU, 10GB)
  init_hist.dat     copy of source last_hist.dat  (init_hist_file points here)
  last_conf.dat     copy of source last_conf.dat  (new conf_file)
  <system>.top      copy of topology
  wfile.txt         symlink to same weight file as source (WEIGHTS UNCHANGED)
  op-setA.txt       symlink to ../op-setA.txt         (HB window)
  op-setA_dist.txt  symlink to ../op-setA_dist.txt    (dist window)

Observables (data_output_1 in input_VMMC, print_every=100):
  HB window   hb_obs.dat   col: step | hb_list count (via op-setA.txt pairs)
  dist window dist_obs.dat col: step | distance between first stem pair

Usage:
    python3 setup_continuation.py           # create all v6 dirs
    python3 setup_continuation.py --dry-run # preview only, no files written
"""

import os, shutil, random, argparse

ROOT = os.path.dirname(os.path.abspath(__file__))

# ── Seed selections ──────────────────────────────────────────────────────────
# (system, window, source_version, combo, seed_id)
# HB from v4, dist from v5.

SEEDS = [
    # ── pulling HB — v4 (wA-s1 has sigma=0.003; wB/wC cluster 7-9 kT) ──────
    ('pulling',  'HB',   'v4', 'wA', 's1'),
    ('pulling',  'HB',   'v4', 'wB', 's2'),
    ('pulling',  'HB',   'v4', 'wB', 's4'),
    ('pulling',  'HB',   'v4', 'wC', 's2'),
    ('pulling',  'HB',   'v4', 'wC', 's5'),
    # ── pulling dist — v5 (block_diff 0.09-0.6; ΔG -7 to -9 kT cluster) ────
    ('pulling',  'dist', 'v5', 'wC', 's9'),
    ('pulling',  'dist', 'v5', 'wA', 's3'),
    ('pulling',  'dist', 'v5', 'wB', 's1'),
    ('pulling',  'dist', 'v5', 'wA', 's5'),
    ('pulling',  'dist', 'v5', 'wC', 's1'),
    # ── noclamps HB — v4 (block_diff 0.69-0.93 for wB; σ 0.002-0.004) ──────
    ('noclamps', 'HB',   'v4', 'wB', 's3'),
    ('noclamps', 'HB',   'v4', 'wB', 's2'),
    ('noclamps', 'HB',   'v4', 'wA', 's1'),
    ('noclamps', 'HB',   'v4', 'wA', 's3'),
    ('noclamps', 'HB',   'v4', 'wC', 's1'),
    # ── noclamps dist — v5 (only -3 to -7 kT cluster; skip stuck >0 seeds) ─
    ('noclamps', 'dist', 'v5', 'wC', 's2'),
    ('noclamps', 'dist', 'v5', 'wA', 's7'),
    ('noclamps', 'dist', 'v5', 'wA', 's3'),
    ('noclamps', 'dist', 'v5', 'wC', 's7'),
    ('noclamps', 'dist', 'v5', 'wA', 's2'),
    # ── pushing HB — v4 (wA-s1 block_diff=0.24; all σ 0.002-0.003) ─────────
    ('pushing',  'HB',   'v4', 'wA', 's1'),
    ('pushing',  'HB',   'v4', 'wA', 's3'),
    ('pushing',  'HB',   'v4', 'wB', 's3'),
    ('pushing',  'HB',   'v4', 'wB', 's5'),
    ('pushing',  'HB',   'v4', 'wC', 's3'),
    # ── pushing dist — v5 (best-converged window; σ 0.001-0.003) ─────────────
    ('pushing',  'dist', 'v5', 'wC', 's9'),
    ('pushing',  'dist', 'v5', 'wC', 's10'),
    ('pushing',  'dist', 'v5', 'wA', 's9'),
    ('pushing',  'dist', 'v5', 'wA', 's3'),
    ('pushing',  'dist', 'v5', 'wB', 's9'),
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

# Weight files — must match source version (frozen by init_hist_file)
HB_WFILE = {
    'wA': 'wfile_A.txt',
    'wB': 'wfile_B.txt',
    'wC': 'wfile_C.txt',
}
DIST_WFILE = {
    ('pulling',  'wA'): 'wfile_v5_dist_pull_A.txt',
    ('pulling',  'wB'): 'wfile_v5_dist_pull_B.txt',
    ('pulling',  'wC'): 'wfile_v5_dist_pull_C.txt',
    ('noclamps', 'wA'): 'wfile_v5_dist_nocl_A.txt',
    ('noclamps', 'wB'): 'wfile_v5_dist_nocl_B.txt',
    ('noclamps', 'wC'): 'wfile_v5_dist_nocl_C.txt',
    ('pushing',  'wA'): 'wfile_v5_dist_push_A.txt',
    ('pushing',  'wB'): 'wfile_v5_dist_push_B.txt',
    ('pushing',  'wC'): 'wfile_v5_dist_push_C.txt',
}

# First stem pair per system — used as distance observable in dist window.
# Pairs from CLAUDE.md (v4/v5 topology indices).
DIST_OBS_PAIR = {
    'pulling':  (15, 56),
    'noclamps': (54, 95),
    'pushing':  (99, 140),
}

# ── Starting conformations ────────────────────────────────────────────────────
# HB window: start from fully-bonded relaxed conf (bond=6).
# Relative to ROOT/SYS_DIR[system].  All three systems already have this file.
HB_START_CONF = {
    'pulling':  ('v4/wA_s1', '6nt-pulling-relaxed-NEW.dat'),
    'noclamps': ('.',         '6nt-noclamps-relaxed-NEW.dat'),
    'pushing':  ('.',         '6nt-pushing-relaxed-NEW.dat'),
}

# Dist window: start from bond=0, dist_bin≈0 conf (strands close together).
# These were the original v5 dist starting confs.
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


def _obs_block(is_dist, system):
    """
    Inline data_output_1 block for the input file.

    HB window:   hb_obs.dat  — step | hb_list count via op-setA.txt pairs
                 Gives HB count for the specific stem pairs at 100-step resolution.
                 Use this to count bond=6→0→6 round trips without needing
                 energy.dat (which only has 1e4 resolution).

    Dist window: dist_obs.dat — step | distance between first stem pair
                 p1, p2 = first pair from CLAUDE.md for each system.
                 Use to count dist=0→max→0 round trips at 100-step resolution.
    """
    if not is_dist:
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
    else:
        p1, p2 = DIST_OBS_PAIR[system]
        return (
            f"data_output_1 = {{\n"
            f"  name = dist_obs.dat\n"
            f"  print_every = 100\n"
            f"  col_1 = {{\n"
            f"    type = step\n"
            f"  }}\n"
            f"  col_2 = {{\n"
            f"    type = distance\n"
            f"    particle_1 = {p1}\n"
            f"    particle_2 = {p2}\n"
            f"    PBC = true\n"
            f"  }}\n"
            f"}}"
        )


def _input_vmmc(top_base, op_file, is_dist, seed, system):
    obs = _obs_block(is_dist, system)
    return f"""\
####  PROGRAM PARAMETERS  ####
backend = CPU
seed = {seed}
interaction_type = DNA2
####    SIM PARAMETERS    ####
sim_type = VMMC
ensemble = NVT
steps = 1e10
maxclust = 50
delta_translation = 0.25
delta_rotation = 0.15
T = 20C
verlet_skin = 1.25
umbrella_sampling = 1
op_file = {op_file}
weights_file = wfile.txt
salt_concentration = 0.5
max_backbone_force = 5
max_backbone_force_far = 10
### Seq. Dependence ##########
use_average_seq = 1
small_system = 1
####    INPUT / OUTPUT    ####
topology = {top_base}.top
conf_file = last_conf.dat
trajectory_file = trajectory.dat
no_stdout_energy = 0
restart_step_counter = 1
energy_file = energy.dat
print_conf_interval = 1e4
print_energy_every = 1e4
time_scale = linear
default_weight = 0.
safe_weights = 0
init_hist_file = init_hist.dat
{obs}
"""


# ── Core setup function ──────────────────────────────────────────────────────

def setup_run(system, window, src_ver, combo, seed_id, dry_run=False):
    """Create one v6 run directory.  Returns True on success, False if skipped."""
    is_dist  = (window == 'dist')
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]

    run_name = f"dist_{combo}_{seed_id}" if is_dist else f"{combo}_{seed_id}"
    src_dir  = os.path.join(sys_path, src_ver, run_name)
    dst_dir  = os.path.join(sys_path, 'v6',    run_name)
    v6_dir   = os.path.join(sys_path, 'v6')

    # ── Resolve the correct starting conf ────────────────────────────────────
    if is_dist:
        conf_subdir, conf_fname = DIST_START_CONF[system]
        start_conf_path = os.path.join(sys_path, conf_subdir, conf_fname)
    else:
        conf_subdir, conf_fname = HB_START_CONF[system]
        start_conf_path = os.path.join(sys_path, conf_subdir, conf_fname)

    # ── Source validation ────────────────────────────────────────────────────
    if not os.path.isdir(src_dir):
        print(f"  SKIP  source not found:  {src_dir}")
        return False
    for fpath, label_ in [(os.path.join(src_dir, 'last_hist.dat'), 'last_hist.dat'),
                          (os.path.join(src_dir, f'{top_base}.top'), f'{top_base}.top'),
                          (start_conf_path, conf_fname)]:
        if not os.path.exists(fpath):
            print(f"  SKIP  missing {label_}:  {fpath}")
            return False

    label = (f"{SYS_DIR[system]}/v6/{run_name}  ← {src_ver}/{run_name}"
             f"  [conf: {conf_fname}]")
    print(f"  {label}")
    if dry_run:
        return True

    # ── Create destination ───────────────────────────────────────────────────
    os.makedirs(dst_dir, exist_ok=True)

    # Hard-copy topology, starting conf, and history
    shutil.copy2(os.path.join(src_dir, f'{top_base}.top'), dst_dir)
    shutil.copy2(start_conf_path, os.path.join(dst_dir, 'last_conf.dat'))
    shutil.copy2(os.path.join(src_dir, 'last_hist.dat'),
                 os.path.join(dst_dir, 'init_hist.dat'))

    # ── Weight file symlink ──────────────────────────────────────────────────
    wfile_name = HB_WFILE[combo] if not is_dist else DIST_WFILE[(system, combo)]
    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{wfile_name}')

    # ── Op file: ensure real copy exists in v6/, then symlink from run dir ──
    op_fname = 'op-setA_dist.txt' if is_dist else 'op-setA.txt'
    v6_op    = os.path.join(v6_dir, op_fname)
    if not os.path.exists(v6_op):
        # Copy real file from v5 (identical to v4 for HB; v5 is authoritative for dist)
        for src_ver_op in ['v5', 'v4']:
            candidate = os.path.realpath(os.path.join(sys_path, src_ver_op, op_fname))
            if os.path.isfile(candidate):
                shutil.copy2(candidate, v6_op)
                break
        else:
            print(f"  WARNING: could not find {op_fname} to copy into v6/")
    _symlink(os.path.join(dst_dir, op_fname), f'../{op_fname}')

    # ── run_VMMC.sh ──────────────────────────────────────────────────────────
    sh_path = os.path.join(dst_dir, 'run_VMMC.sh')
    with open(sh_path, 'w') as f:
        f.write(RUN_VMMC_SH)
    os.chmod(sh_path, 0o755)

    # ── input_VMMC ───────────────────────────────────────────────────────────
    new_seed = random.randint(10000, 99999)
    with open(os.path.join(dst_dir, 'input_VMMC'), 'w') as f:
        f.write(_input_vmmc(top_base, op_fname, is_dist, new_seed, system))

    return True


def _symlink(link_path, target):
    """Create a symlink, replacing any existing one."""
    if os.path.lexists(link_path):
        os.remove(link_path)
    os.symlink(target, link_path)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Set up v6 continuation runs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--dry-run', action='store_true',
                    help='Preview what would be created, write nothing')
    args = ap.parse_args()

    random.seed()   # use system entropy for new oxDNA seeds

    if args.dry_run:
        print("DRY RUN — no files will be written\n")

    n_ok = n_skip = 0
    cur_system = None
    for system, window, src_ver, combo, seed_id in SEEDS:
        if system != cur_system:
            cur_system = system
            print(f"\n{'─'*60}\n{system.upper()}\n{'─'*60}")
        result = setup_run(system, window, src_ver, combo, seed_id,
                           dry_run=args.dry_run)
        if result:
            n_ok += 1
        else:
            n_skip += 1

    print(f"\n{'Would create' if args.dry_run else 'Created'}: {n_ok}   "
          f"Skipped: {n_skip}\n")

    if not args.dry_run and n_ok > 0:
        print("Next steps:")
        print("  1. Spot-check a few input_VMMC files")
        print("  2. pushnew the three v6/ trees and links/6nt/ to the cluster:")
        print("       pushnew oxview-assemblies/pulling-clamps/6nt-clamps-pulling/v6/")
        print("       pushnew oxview-assemblies/noclamps/6nt-noclamps/v6/")
        print("       pushnew oxview-assemblies/pushing-clamps/6nt-pushing/v6/")
        print("  3. Submit on cluster:")
        print("       for d in */*/v6/*_s*/; do")
        print("         cd \"$REMOTE_SYNC_ROOT/$d\" && qsub run_VMMC.sh && cd \"$REMOTE_SYNC_ROOT\"")
        print("       done")


if __name__ == '__main__':
    main()
