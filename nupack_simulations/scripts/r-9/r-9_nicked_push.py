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
rod_length = 54  # length of the rod domain in strand 3
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

# Lettuce r-9-1 sequences (abstracted here for design simplicity)
overhang5p = 'GAG'
stem1a  = 'GC'
stem1_bump = 'TC'
stem1b  = 'TCGGGACGAC'
hairpin='AGGGATCAGGGATGTAGCGGCAGTGGGTGA'
stem2  = 'GTCGTCCCGAGC'
overhang3p = 'TGCAATCGTAA'

clamp5p  = Domain(f'N{clamp_length}', name='Domain 1')   # N left clamp length
d1  = Domain(overhang5p, name='overhang5p') # d1-d7 are the lettuce domains.
d2  = Domain(stem1a, name='stem1a')
d3  = Domain(stem1_bump, name='stem1_bump')
d4 = Domain(stem1b, name='stem1b')
d5 = Domain(hairpin, name='hairpin')
d6 = Domain(stem2, name='stem2')
d7 = Domain(overhang3p, name='overhang3p')
clamp3p = Domain(f'N{clamp_length}', name='Domain 5')    # right clamp length
d8 = Domain(f'N{blankregion}', name='Domain 6')    # loop - rod = unpaired region
d9 = Domain(f'N{rod_length}', name='Domain 7')    # Legacy strand for rod
d10 = Domain(f'N{blankregion}', name='Domain 8')    # loop - rod = unpaired region
d11 = Domain(f'N{rod_length/2}', name = 'Domain 9')
d12 = Domain(f'N{rod_length/2}', name = 'Domain 10')

# Strand 1 = Lettuce_r-9-1
Strand1 = TargetStrand(
    [clamp5p, d2, d3, d4, d5, d6, clamp3p],
    name='Lettuce-r9'
)

# Strand 2 = Nicked loop strand - left side
Strand2 = TargetStrand(
    [~clamp5p, d8, d11],
    name='Nicked loop strand - left'
)

# Strand 3 = Nicked loop strand - right side
Strand3 = TargetStrand(
    [d12, d10, ~clamp3p],
    name='Nicked loop strand - right'
)

# Strand 6 = Nicked rod strand - left side
Strand4 = TargetStrand(
    [~d11],
    name='Nicked rod strand - left'
)

# Strand 7 = Nicked rod strand - right side
Strand5 = TargetStrand(
    [~d12],
    name='Nicked rod strand - right'
)

lettuce_structure = f"){clamp_length}(2.2(10.30)12({clamp_length}"
rightrod_structure = f"({rod_length/2}"
rightloop_structure = f"){rod_length/2}.{blankregion}({clamp_length}"
leftloop_structure = f"){clamp_length}.{blankregion}({rod_length/2}"
leftrod_structure = f"){rod_length/2}"

structure = rightrod_structure+"+"+rightloop_structure+"+"+lettuce_structure+"+"+leftloop_structure+"+"+leftrod_structure

# complex when only lettuce and loop are interacting
C1 = TargetComplex(
    [Strand5, Strand3, Strand1, Strand2, Strand4],
    structure=structure,
    name='Complex 1'
)

#print(f"The number of nucleotides in the target structure:{C1.nt(model.alphabet)}")

tube1 = TargetTube(
    on_targets={C1: 1e-6},
    off_targets=SetSpec(max_size=2),
    name='Tube 1'
)


# ---------- Design job (tube_design) ----------

# GC content constraint
sim1 = Similarity([clamp5p], f'S{clamp_length}', limits=[0.40, 0.55])
sim2 = Similarity([clamp3p], f'S{clamp_length}', limits=[0.40, 0.55])
sim3 = Similarity([d8], f'S{blankregion}', limits=[0.40, 0.55])
sim4 = Similarity([d10], f'S{blankregion}', limits=[0.40, 0.55])
sim5 = Similarity([d11], f'S{rod_length/2}', limits=[0.40, 0.55])
sim6 = Similarity([d12], f'S{rod_length/2}', limits=[0.40, 0.55])
# Diversity constraint
div1 = Diversity(word=10, types=3, scope=[clamp5p, clamp3p, d8, d10, d11, d12])
pat1 = Pattern(['GGGG', 'CCCC'], scope=[clamp5p, clamp3p, d8, d10, d11, d12]) # avoid G-quadruplexes


design = tube_design(
    tubes=[tube1],
    hard_constraints=[sim1, sim2, sim3, sim4, sim5, sim6, div1, pat1],
    soft_constraints=[],
    defect_weights=None,
    options=None,
    model=model
)

# Run 10 trials
results = design.run(trials=10)

for i, r in enumerate(results):
    ed = r.defects.ensemble_defect   # scalar
    print(f"Trial {i}: ensemble defect = {ed:.4f}")


best_idx, best_result = min(
    enumerate(results),
    key=lambda kv: kv[1].defects.ensemble_defect
)

print(f"Best trial: {best_idx}, ensemble defect = {best_result.defects.ensemble_defect:.4f}")


target_strands = [Strand1, Strand2, Strand3, Strand4, Strand5]

print("Best strand sequences:")
for ts in target_strands:
    designed_strand = best_result.to_analysis(ts)
    seq = str(designed_strand)
    print(f"{ts.name}: {seq}")

os.chdir('../')
os.chdir('../output/r-9')
with open("strands3_r9_nicked_push.txt", "w") as f:
    f.write(f"Best trial ensemble defect: {best_result.defects.ensemble_defect:.4f}\n")
    for ts in target_strands:
        seq = str(best_result.to_analysis(ts))
        f.write(f"{ts.name}\t{seq}\n")
    print(f"Done! Saved best sequences to {f.name}")
