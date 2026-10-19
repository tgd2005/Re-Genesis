#!/usr/bin/env python3
"""
setup_v10.py  —  create v10 run directories with manual intuitive weights
=========================================================================

Manual weight combos based on analysis of v7 results (fresh start, no init_hist).

HB window  — 5 combos (wD/wE/wF/wG/wH) × 10 seeds:
  Key insight: pushing HB converges well because b=0 is accidentally over-weighted
  (wA weight ≈ 22× ideal). Strategy: intentionally over-weight b=0 to force more
  b=0 sampling and accelerate round trips.

  D: v9_A but w(b0) ×10            — moderate b=0 boost (system-specific)
  E: v9_A but w(b0) ×20            — aggressive b=0 boost (system-specific)
  F: v9_A but w(b0) ×5, w(b2) ×2  — mild b=0 boost + nudge through b=2 (system-specific)
  G: pushing v7/v8 weights for ALL  — use known-good pushing weights everywhere
  H: pushing v9_A ideal for ALL     — use pushing ideal weights everywhere

Dist window — 4 combos (wD/wE/wF/wG) × 10 seeds:
  Key insight: pulling dist achieves round trips; use pulling weights elsewhere.
  Pushing dist is physically stuck (mb5/mb0 ratio 10³–10⁵); try extreme mb0 boosts.
  Noclamps has a real free energy barrier at mb3→mb4.

  D: targeted system-specific (pull mb0×4, nocl mb3×10+mb4×5, push mb0×5)
  E: aggressive system-specific (pull mb0×10, nocl flatten mb0-mb4, push mb0×20)
  F: pulling v7 wA weights for ALL  — [5130,122,14.2,4.11,1.51,1,1]
  G: pulling v9_A weights for ALL   — [1604,59.9,7.74,2.06,1.01,1,2.99]

Wfiles written to links/6nt/:
  - wfile_v10_hb_{pull/nocl/push}_{D/E/F/G/H}.txt
  - wfile_v10_dist_{pull/nocl/push}_{D/E/F/G}.txt

Total: 150 HB + 120 dist = 270 runs.

Usage:
    python3 setup_v10.py           # write wfiles, create all v10 dirs
    python3 setup_v10.py --dry-run # preview only, no files written
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

SYS_ABBREV = {'pulling': 'pull', 'noclamps': 'nocl', 'pushing': 'push'}

COMBOS_HB   = ['wD', 'wE', 'wF', 'wG', 'wH']
COMBOS_DIST = ['wD', 'wE', 'wF', 'wG']
N_SEEDS = 10

# ── Hardcoded weight values (precomputed from v7 pooled data) ─────────────────

# v9_A ideal weights (reference values)
V9_HB_A = {
    'pulling':  [4615,    1077400, 54268, 3650,  265.9,  19.2,  1],
    'noclamps': [1473,     597000, 31680, 2310,  156.3,  13.5,  1],
    'pushing':  [ 276,     139300, 12850,  982,   61.05,  4.703, 1],
}

V9_DIST_A = {
    'pulling':  [1604,   59.9,   7.74,  2.06,  1.01,  1,     2.99],
    'noclamps': [43.98,  14.61,  68.24, 18.15, 3.858, 1,     1.099],
    'pushing':  [40898, 1037,   162.0,  18.74, 3.184, 1,     1.253],
}

# Pushing v7/v8 HB weights (same wfile used in v7 wA for pushing) — used as wG
_PUSH_V7_HB = [6056.85, 1055510, 33208.79, 1958.96, 154.06, 14.23, 1]

# Pulling v7 dist weights (wfile_v7_dist_pull_A) — used as dist wF
_PULL_V7_DIST = [5130, 122, 14.2, 4.11, 1.51, 1, 1]

# v10 HB weight combos: system → combo → list of 7 floats [w_b0..w_b6]
# wD/wE/wF are system-specific (based on v9_A ideal); wG/wH same for all.
V10_HB = {
    'pulling': {
        # wD: b0 ×10
        'wD': [46150,  1077400,  54268, 3650,  265.9, 19.2,  1],
        # wE: b0 ×20
        'wE': [92300,  1077400,  54268, 3650,  265.9, 19.2,  1],
        # wF: b0 ×5, b2 ×2
        'wF': [23075,  1077400, 108536, 3650,  265.9, 19.2,  1],
        # wG: pushing v7/v8 weights for all (same as known-good pushing wA)
        'wG': _PUSH_V7_HB,
        # wH: pushing v9_A ideal for all
        'wH': [276,     139300,  12850,  982,   61.05,  4.703, 1],
    },
    'noclamps': {
        'wD': [14730,   597000,  31680, 2310,  156.3, 13.5,  1],
        'wE': [29460,   597000,  31680, 2310,  156.3, 13.5,  1],
        'wF': [ 7365,   597000,  63360, 2310,  156.3, 13.5,  1],
        'wG': _PUSH_V7_HB,
        'wH': [276,     139300,  12850,  982,   61.05,  4.703, 1],
    },
    'pushing': {
        'wD': [ 2760,   139300,  12850,  982,   61.05,  4.703, 1],
        'wE': [ 5520,   139300,  12850,  982,   61.05,  4.703, 1],
        'wF': [ 1380,   139300,  25700,  982,   61.05,  4.703, 1],
        'wG': _PUSH_V7_HB,
        'wH': [276,     139300,  12850,  982,   61.05,  4.703, 1],
    },
}

# v10 dist weight combos: system → combo → list of 7 floats [w_mb0..w_mb6]
# wD/wE are system-specific; wF/wG are pulling weights applied to all systems.
V10_DIST = {
    'pulling': {
        # wD: targeted — mb0 ×4 (from v9_A)
        'wD': [  6416,   59.9,   7.74,  2.06, 1.01, 1, 2.99],
        # wE: aggressive — mb0 ×10
        'wE': [ 16040,   59.9,   7.74,  2.06, 1.01, 1, 2.99],
        # wF: pulling v7 wA weights for all
        'wF': _PULL_V7_DIST,
        # wG: pulling v9_A for all
        'wG': [  1604,   59.9,   7.74,  2.06, 1.01, 1, 2.99],
    },
    'noclamps': {
        # wD: targeted — mb3 ×10, mb4 ×5 (from v9_A)
        'wD': [ 43.98,  14.61,  68.24, 181.5, 19.29, 1, 1.099],
        # wE: aggressive — flatten mb0-mb4 to 19.07 (semi-flat histogram)
        'wE': [ 19.07,  19.07,  19.07,  19.07, 19.07, 1, 1.099],
        'wF': _PULL_V7_DIST,
        'wG': [  1604,   59.9,   7.74,   2.06,  1.01, 1, 2.99],
    },
    'pushing': {
        # wD: targeted — mb0 ×5 (from v9_A)
        'wD': [204490,  1037,  162.0,  18.74, 3.184, 1, 1.253],
        # wE: aggressive — mb0 ×20
        'wE': [817960,  1037,  162.0,  18.74, 3.184, 1, 1.253],
        'wF': _PULL_V7_DIST,
        'wG': [  1604,   59.9,   7.74,   2.06,  1.01, 1, 2.99],
    },
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


def find_relaxed_conf(system):
    """Find the original relaxed bonded conf (bond=6 starting structure).

    Checks:
      - pulling:  v4/wA_s1/{top_base}.dat
      - noclamps: sys_path/{top_base}.dat
      - pushing:  sys_path/{top_base}.dat
    Falls back to any v4/v5/v7 run dir containing the dat file.
    """
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    dat_fname = f'{top_base}.dat'

    if system == 'pulling':
        primary = os.path.join(sys_path, 'v4', 'wA_s1', dat_fname)
    else:
        primary = os.path.join(sys_path, dat_fname)

    if os.path.exists(primary):
        return primary

    for ver in ['v7', 'v5', 'v4']:
        pattern = os.path.join(sys_path, ver, '*', dat_fname)
        matches = glob.glob(pattern)
        if matches:
            return matches[0]

    return None


# ── Wfile writing ─────────────────────────────────────────────────────────────

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
    """Write all v10 wfiles to links/6nt/ and print a summary table."""
    links_dir = os.path.join(ROOT, 'links', '6nt')

    print("=" * 60)
    print("V10 WEIGHT FILES")
    print("=" * 60)

    for system in ['pulling', 'noclamps', 'pushing']:
        abbrev = SYS_ABBREV[system]
        hb_header = "  ".join(f"{c:>12}" for c in COMBOS_HB)
        print(f"\n{system.upper()} HB weights:")
        print(f"  {'Bin':>4}  {hb_header}")
        for b in range(7):
            row = "  " + f"b={b}  " + "  ".join(
                f"{V10_HB[system][c][b]:>12.4g}" for c in COMBOS_HB)
            print(row)

        dist_header = "  ".join(f"{c:>12}" for c in COMBOS_DIST)
        print(f"\n{system.upper()} dist weights:")
        print(f"  {'Bin':>4}  {dist_header}")
        for mb in range(7):
            row = "  " + f"mb={mb} " + "  ".join(
                f"{V10_DIST[system][c][mb]:>12.4g}" for c in COMBOS_DIST)
            print(row)

        if not dry_run:
            os.makedirs(links_dir, exist_ok=True)
            for combo in COMBOS_HB:
                letter = combo[1]
                hb_fname = f'wfile_v10_hb_{abbrev}_{letter}.txt'
                write_hb_wfile(os.path.join(links_dir, hb_fname),
                               V10_HB[system][combo])
                print(f"  Wrote {hb_fname}")
            for combo in COMBOS_DIST:
                letter = combo[1]
                dist_fname = f'wfile_v10_dist_{abbrev}_{letter}.txt'
                write_dist_wfile(os.path.join(links_dir, dist_fname),
                                 V10_DIST[system][combo])
                print(f"  Wrote {dist_fname}")


# ── Setup functions ──────────────────────────────────────────────────────────

def setup_hb(system, combo, seed_id, dry_run=False):
    """Create one v10 HB run directory — fresh start, manual weights."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    abbrev   = SYS_ABBREV[system]
    letter   = combo[1]  # D/E/F/G
    run_name = f"{combo}_{seed_id}"
    dst_dir  = os.path.join(sys_path, 'v10', run_name)
    v10_dir  = os.path.join(sys_path, 'v10')

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

    wfile_name = f'wfile_v10_hb_{abbrev}_{letter}.txt'
    print(f"  HB  {SYS_DIR[system]}/v10/{run_name}  [fresh, conf: {os.path.basename(relaxed_conf)}]")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(top_path, os.path.join(dst_dir, f'{top_base}.top'))
    shutil.copy2(relaxed_conf, os.path.join(dst_dir, 'last_conf.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{wfile_name}')

    _ensure_op_file(v10_dir, sys_path, 'op-setA.txt')
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
    """Create one v10 dist run directory — fresh start, manual weights."""
    sys_path = os.path.join(ROOT, SYS_DIR[system])
    top_base = SYS_TOP[system]
    abbrev   = SYS_ABBREV[system]
    letter   = combo[1]
    run_name = f"dist_{combo}_{seed_id}"
    dst_dir  = os.path.join(sys_path, 'v10', run_name)
    v10_dir  = os.path.join(sys_path, 'v10')

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

    wfile_name = f'wfile_v10_dist_{abbrev}_{letter}.txt'
    print(f"  dist {SYS_DIR[system]}/v10/{run_name}  [fresh, conf: {conf_fname}]")
    if dry_run:
        return True

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(top_path, os.path.join(dst_dir, f'{top_base}.top'))
    shutil.copy2(start_conf_path, os.path.join(dst_dir, 'last_conf.dat'))

    _symlink(os.path.join(dst_dir, 'wfile.txt'),
             f'../../../../links/6nt/{wfile_name}')

    _ensure_op_file(v10_dir, sys_path, 'op-setA_dist.txt')
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
        description='Set up v10 runs (fresh start, manual intuitive weights)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--dry-run', action='store_true',
                    help='Preview weights and what would be created, write nothing')
    args = ap.parse_args()

    random.seed()

    if args.dry_run:
        print("DRY RUN — no files will be written\n")

    # Print (and optionally write) weight files
    write_all_wfiles(dry_run=args.dry_run)

    n_ok = n_skip = 0

    n_hb_combos   = len(COMBOS_HB)   * N_SEEDS * 3
    n_dist_combos = len(COMBOS_DIST) * N_SEEDS * 3
    print()
    print("=" * 60)
    print(f"HB WINDOW  ({len(COMBOS_HB)} combos × {N_SEEDS} seeds × 3 systems = {n_hb_combos} runs; "
          "fresh start, no init_hist)")
    print("=" * 60)
    for system in ['pulling', 'noclamps', 'pushing']:
        print(f"\n{system.upper()}")
        for combo in COMBOS_HB:
            for i in range(1, N_SEEDS + 1):
                ok = setup_hb(system, combo, f's{i}', dry_run=args.dry_run)
                n_ok += ok; n_skip += (not ok)

    print()
    print("=" * 60)
    print(f"DIST WINDOW  ({len(COMBOS_DIST)} combos × {N_SEEDS} seeds × 3 systems = {n_dist_combos} runs; "
          "fresh start, no init_hist)")
    print("=" * 60)
    for system in ['pulling', 'noclamps', 'pushing']:
        print(f"\n{system.upper()}")
        for combo in COMBOS_DIST:
            for i in range(1, N_SEEDS + 1):
                ok = setup_dist(system, combo, f's{i}', dry_run=args.dry_run)
                n_ok += ok; n_skip += (not ok)

    print(f"\n{'Would create' if args.dry_run else 'Created'}: {n_ok}   "
          f"Skipped: {n_skip}\n")

    if not args.dry_run and n_ok > 0:
        print("Next steps:")
        print("  1. Spot-check a few input_VMMC and links/6nt/wfile_v10_*.txt")
        print("  2. Push to cluster:")
        print("       pushnew oxview-assemblies/pulling-clamps/6nt-clamps-pulling/v10/")
        print("       pushnew oxview-assemblies/noclamps/6nt-noclamps/v10/")
        print("       pushnew oxview-assemblies/pushing-clamps/6nt-pushing/v10/")
        print("       pushnew oxview-assemblies/links/6nt/")
        print("  3. Submit on cluster:")
        print("       for d in */*/v10/*_s*/; do")
        print("         cd \"$REMOTE_SYNC_ROOT/$d\" && qsub run_VMMC.sh && cd \"$REMOTE_SYNC_ROOT\"")
        print("       done")


if __name__ == '__main__':
    main()
