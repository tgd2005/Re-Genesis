# oxDNA/alex: mechanical rupture of an aptamer stem by a bent-duplex clamp

Simulation work by Alexander for **Re-Genesis**, the Imperial College London Year-3 Bioengineering group project (2025–26), supervised by Dr Thomas Ouldridge and Dr Max Earle.

## Question
Can a stiff DNA rod, hybridised onto a duplex "clamp", open the stem of an aptamer mechanically, and does the clamp have to **push** or **pull** to do it? The aim is a force-driven, target-agnostic way to switch an aptamer off.

## What was done
- **Force–extension (oxDNA MD, harmonic traps).** Arrays over trap stiffness, trap separation and seeds for full-, half-, 14- and 40-nt rods, and the frequency of nick kinking (loss of coaxial stacking) against extension.
- **Stem-rupture free energy (oxDNA2 VMMC umbrella sampling).** Pulling clamp, pushing clamp and a no-clamp reference. Two windows (H-bond count; minimum pair distance) are joined at their shared state. The 6-nt stem is the production system: v8 gives the result below, and v9–v10 test re-derived and hand-tuned weights. Relaxed 3-, 9- and 12-nt stems were built but not run to production.

## Result
6-nt stem, oxDNA2, 20 °C, 0.5 M NaCl. ΔF = F(ruptured) − F(bonded), mean ± SD over 10 seeds:

| Clamp | ΔF (kT) |
|---|---|
| Pulling | −0.002 ± 0.317 |
| Pushing | **−6.835 ± 0.208** |

A pushing clamp makes the open stem favourable by about 7 kT; a pulling clamp has no effect. Plate-reader FRET assays on the same devices point the same way: pushing devices lost stem fluorescence, pulling devices matched the control, and a 12-nt stem responded less than a 6-nt one.
Data and figure: `legacy/single-duplex/vmmc-umbrella-sampling-murphy/free_energy_results/20260414_214053_*/` (`dF_results.csv`, `dF_macrostate.pdf`).

## Layout
| Path | Contents |
|---|---|
| `force_dist_array_*/`, `force_dist_v2/`, `assem_*/`, `traps*/` | First force–extension, assembly and trap runs |
| `force-extension/` | Later force–extension arrays |
| `oxview-assemblies/` | Cluster copy of the VMMC run tree (`pulling-clamps/`, `pushing-clamps/`, `noclamps/`, v1–v10) with the early daily notes |
| `assemblies-before-oxview/`, `pitch/`, `insulation/`, `attic/` | Early assemblies, helical-pitch and insulation tests, superseded runs |
| `legacy/single-duplex/vmmc-umbrella-sampling-murphy/` | Local copy of the same VMMC tree, with the shared weight files (`links/6nt/`) and the analysis: `free_energy.py` (final ΔF), `diagnostics.py` (convergence), `setup_v*.py`, results, and `METHOD_NOTES.md` (parameters, run structure, known issues) |
| `legacy/single-duplex/force-extension-murphy/` | Kinking-seed arrays, `kinking_*.py` and figures, `01_f_extensions.ipynb` |
| `legacy/single-duplex/alex-workspace/` | Earlier copies of files that also appear in the folders above |

**Kept:** inputs (oxDNA inputs, topologies, starting and last configurations, order-parameter and weight files, PBS scripts), `last_hist.dat` histograms, small observables, analysis code and figures. **Not kept:** trajectories, per-step time series above 2 MB, pair-energy dumps and PBS logs, all of which can be regenerated from the inputs. This includes about 176 GB of raw output from the directional-pulling trap runs (`attic/traps_v2`).

## Next
The work continues outside this repo as a POLNET Physics of Life summer bursary project with Dr Ouldridge: switch-strand actuation of the stem (a 2D free-energy landscape) and coaxial stacking at nicked junctions. Open items from here: check v8 against the v9–v10 weight variants, run the 3/9/12-nt stems, and extend to multi-duplex (4HB/6HB) devices.
