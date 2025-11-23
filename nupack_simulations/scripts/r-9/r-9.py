from nupack import *
import os


"""
Design two DNA strands for a target duplex using the NUPACK 4 Python API.

Requirements:
- Python 3 and the `nupack` package installed.

Inputs:
- No command-line arguments; all model, domains, and constraints are defined in this script.

Outputs:
- Prints ensemble defect for 100 design trials and the best trial's strand sequences.
- Saves the best Strand 1, Strand 2, and Strand 3 sequences to 'best_strands.txt'.
"""

#### PARAMETERS TO ADJUST ####
loop_length = 54  # length of the loop domain in strand 2
rod_length = 26  # length of the rod domain in strand 3
clamp_length = 15  # length of the clamp domains in strands 1 and 2

# NOTICE: loop_length - rod_length / 2 should be an integer.
blankregion = (loop_length - rod_length) / 2
if not blankregion.is_integer():
    raise ValueError("loop_length - rod_length must be an even number.")
blankregion = int(blankregion)


model = Model(
    material='dna',          # DNA
    ensemble='stacking',
    celsius=27,              # Temperature 27 °C
    sodium=1.0,              # [Na+] = 1 M
    magnesium=0.0            # [Mg2+] = 0 M
)

options = DesignOptions(
    f_stop=0.04,
    seed=12345,
)

d1  = Domain(f'N{clamp_length}', name='Domain 1')   # N right clamp length
d2  = Domain("GAG", name="5'Overhang") # d2-d8 are the lettuce domains.
d3  = Domain("GCTCGGGACGAC", name="5'Stem")
d4  = Domain('AGGGATCAGGGATGTAGCGGCAGTGGGTGA', name='LettuceCore')
d5  = Domain("GTCGTCCCGAGC", name="3'Stem")
d6  = Domain("TGCAATCGTAA", name="3'Overhang")
d7  = Domain(f'N{clamp_length}', name='Domain 7')    # left clamp length
d8 = Domain(f'N{blankregion}', name='Domain 8')    # loop - rod = unpaired region
d9  = Domain(f'N{rod_length}', name='Domain 9')    # N20; Rod structure, variable length in simulations.
d10 = Domain(f'N{blankregion}', name='Domain 10')    # loop - rod = unpaired region

# Strand 1 = Lettuce_S-9-1
Strand1 = TargetStrand(
    [d1, d2, d3, d4, d5, d6, d7],
    name='Lettuce_R-9-1'
)

# Strand 2 = Clamp loop structure
Strand2 = TargetStrand(
    [~d1, d8 , d9, d10, ~d7],
    name='Clamp loop strand'
)

# Strand 3 = Rod structure
Strand3 = TargetStrand(
    [~d9],
    name='Clamp rod strand'
)

lettuce_r9_structure = '.3(12.30)12.11'
# legacy structure for 2 strands: '(15(22.22(10.5)10.4)22(15+)15.54)15'
s1s2_complex_structure = f'({clamp_length}' + lettuce_r9_structure + f'({clamp_length}+){clamp_length}.{blankregion}.{rod_length}.{blankregion}){clamp_length}'
s2s3_complex_structure = f'.{clamp_length}.{blankregion}({rod_length}.{blankregion}.{clamp_length}+){rod_length}'

# complex when only lettuce and loop are interacting
C1 = TargetComplex(
    [Strand1, Strand2],
    structure=s1s2_complex_structure,
    name='Complex 1'
)

# complex when loop and rod are interacting
C2 = TargetComplex(
    [Strand2, Strand3],
    structure=s2s3_complex_structure,
    name='Complex 2'
)

#print(f"The number of nucleotides in the target structure:{C1.nt(model.alphabet)}")

tube1 = TargetTube(
    on_targets={C1: 1e-6},
    off_targets=SetSpec(max_size=2),
    name='Tube 1'
)
tube2 = TargetTube(
    on_targets={C2: 1e-6},
    off_targets=SetSpec(max_size=2),
    name='Tube 2'
)

# ---------- Design job (tube_design) ----------

# GC content constraint
sim1 = Similarity([d1], f'S{clamp_length}', limits=[0.40, 0.55])
sim2 = Similarity([d7], f'S{clamp_length}', limits=[0.40, 0.55])
sim3 = Similarity([d8], f'S{blankregion}', limits=[0.40, 0.55])
sim4 = Similarity([d9], f'S{rod_length}', limits=[0.40, 0.55])
sim5 = Similarity([d10], f'S{blankregion}', limits=[0.40, 0.55])
# Diversity constraint
div1 = Diversity(word=10, types=3, scope=[d1, d7, d8, d9, d10])
# Pattern constraint
pat1 = Pattern(['GGGG'], scope=[d1, d7, d8, d9, d10])  # No GGGG in variable domains


design = tube_design(
    tubes=[tube1, tube2],
    hard_constraints=[sim1, sim2, sim3, sim4, sim5, div1, pat1],
    soft_constraints=[],
    defect_weights=None,
    options=None,
    model=model
)

# Run 100 trials
results = design.run(trials=100)

for i, r in enumerate(results):
    ed = r.defects.ensemble_defect   # scalar
    print(f"Trial {i}: ensemble defect = {ed:.4f}")


best_idx, best_result = min(
    enumerate(results),
    key=lambda kv: kv[1].defects.ensemble_defect
)

print(f"Best trial: {best_idx}, ensemble defect = {best_result.defects.ensemble_defect:.4f}")



target_strands = [Strand1, Strand2, Strand3]

print("Best strand sequences:")
for ts in target_strands:
    designed_strand = best_result.to_analysis(ts)
    seq = str(designed_strand)
    print(f"{ts.name}: {seq}")

os.chdir('../')
os.chdir("../output") # Change directory to output folder for saving results
with open("strands3_r91.txt", "w") as f:
    f.write(f"Best trial ensemble defect: {best_result.defects.ensemble_defect:.4f}\n")
    for ts in target_strands:
        seq = str(best_result.to_analysis(ts))
        f.write(f"{ts.name}\t{seq}\n")
    print(f"Done! Saved best sequences to {f.name}")
