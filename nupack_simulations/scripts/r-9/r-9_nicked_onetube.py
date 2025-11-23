from nupack import *


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
rod_length = 36  # length of the rod domain in strand 3
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


stem1  = 'GCTCGGGACGAC'
stem2  = 'GTCGTCCCGAGC'
hairpin='AGGGATCAGGGATGTAGCGGCAGTGGGTGA'

d1  = Domain(f'N{clamp_length}', name='Domain 1')   # N left clamp length
d2  = Domain(stem1, name='Domain 2') # d2-d4 are the lettuce domains.
d3  = Domain(hairpin, name='Domain 3')
d4  = Domain(stem2, name='Domain 4')
d5 = Domain(f'N{clamp_length}', name='Domain 5')    # right clamp length
d6 = Domain(f'N{blankregion}', name='Domain 6')    # loop - rod = unpaired region
d7 = Domain(f'N{blankregion}', name='Domain 8')    # loop - rod = unpaired region
d8 = Domain(f'N{rod_length/2}', name = 'Domain 9')
d9 = Domain(f'N{rod_length/2}', name = 'Domain 10')

# Strand 1 = Lettuce_r-9-1
Strand1 = TargetStrand(
    [d1, d2, d3, d4, d5],
    name='Lettuce-r9'
)

# Strand 4 = Nicked loop strand - left side
Strand2 = TargetStrand(
    [d8, d6, ~d1],
    name='Nicked loop strand - left'
)

# Strand 5 = Nicked loop strand - right side
Strand3 = TargetStrand(
    [~d5, d7, d9],
    name='Nicked loop strand - right'
)

# Strand 6 = Nicked rod strand - left side
Strand4 = TargetStrand(
    [~d8],
    name='Nicked rod strand - left'
)

# Strand 7 = Nicked rod strand - right side
Strand5 = TargetStrand(
    [~d9],
    name='Nicked rod strand - right'
)


# legacy structure for 2 strands: '(15(22.22(10.5)10.4)22(15+)15.54)15'
nicked_complex_structure = f'({rod_length/2}+){rod_length/2}.{blankregion}({clamp_length}+){clamp_length}({len(stem1)}.{len(hairpin)}){len(stem2)}({clamp_length}+){clamp_length}.{blankregion}({rod_length/2}+){rod_length/2}'



C1 = TargetComplex(
    [Strand4, Strand2, Strand1, Strand3, Strand5],
    structure=nicked_complex_structure,
    name='Complex 1'
)

#print(f"The number of nucleotides in the target structure:{C1.nt(model.alphabet)}")

tube1 = TargetTube(
    on_targets={C1: 1e-6},
    off_targets=SetSpec(max_size=5),
    name='Tube 1'
)


# ---------- Design job (tube_design) ----------

# GC content constraint
sim1 = Similarity([d1], f'S{clamp_length}', limits=[0.35, 0.55])
sim2 = Similarity([d5], f'S{clamp_length}', limits=[0.35, 0.55])
sim3 = Similarity([d6], f'S{blankregion}', limits=[0.35, 0.55])
sim4 = Similarity([d7], f'S{blankregion}', limits=[0.35, 0.55])
# Diversity constraint
div1 = Diversity(word=10, types=3, scope=[d1, d5, d6, d7, d8, d9])
pat1 = Pattern(['GGGG', 'CCCC', 'AAAA', 'TTTT'], scope=[d1, d5, d6, d7, d8, d9]) # avoid G-quadruplexes


design = tube_design(
    tubes=[tube1],
    hard_constraints=[sim1, sim2, sim3, sim4, sim4, div1, pat1],
    soft_constraints=[],
    defect_weights=None,
    options=None,
    model=model
)

# Run 100 trials
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

with open("strands3_r9_nicked.txt", "w") as f:
    f.write(f"Best trial ensemble defect: {best_result.defects.ensemble_defect:.4f}\n")
    for ts in target_strands:
        seq = str(best_result.to_analysis(ts))
        f.write(f"{ts.name}\t{seq}\n")
    print(f"Done! Saved best sequences to {f.name}")
