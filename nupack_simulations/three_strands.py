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
rod_length = 20  # length of the rod domain in strand 3
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
d2  = Domain('AACGTGCTCAAGGTGTCTTAGT', name='Domain 2') # d2-d8 are the lettuce domains.
d3  = Domain('AGGGATGATGCGGCAGTGGGCT', name='Domain 3')
d4  = Domain('TCGCAGAAGA', name='Domain 4')
d5  = Domain('GTGTT', name='Domain 5')
d6  = Domain('TCTTCTGCGA', name='Domain 6')
d7  = Domain('GGGG', name='Domain 7')
d8  = Domain('ACTAAGACACCTTGAGCACGTT', name='Domain 8')
d9  = Domain(f'N{clamp_length}', name='Domain 9')    # left clamp length
d10 = Domain(f'N{blankregion}', name='Domain 10')    # loop - rod = unpaired region
d11  = Domain(f'N{rod_length}', name='Domain 11')    # N20; Rod structure, variable length in simulations.
d12 = Domain(f'N{blankregion}', name='Domain 12')    # loop - rod = unpaired region

# Strand 1 = Lettuce_S-9-1
Strand1 = TargetStrand(
    [d1, d2, d3, d4, d5, d6, d7, d8, d9],
    name='Lettuce_S-9-1'
)

# Strand 2 = Clamp loop structure
Strand2 = TargetStrand(
    [~d1, d10 , d11, d12, ~d9],
    name='Clamp loop strand'
)

# Strand 3 = Rod structure
Strand3 = TargetStrand(
    [~d11],
    name='Clamp rod strand'
)

# legacy structure for 2 strands: '(15(22.22(10.5)10.4)22(15+)15.54)15'
s1s2_complex_structure = '(15(22.22(10.5)10.4)22(15+)15.17.20.17)15'
s2s3_complex_structure = '.15.17(20.17.15+)20'

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

#GC content constraint
sim1 = Similarity([d1], 'S15', limits=[0.35, 0.55])
sim2 = Similarity([d9], 'S15', limits=[0.35, 0.55])
sim3 = Similarity([d10], 'S17', limits=[0.35, 0.55])
sim4 = Similarity([d11], 'S20', limits=[0.35, 0.55])
sim5 = Similarity([d12], 'S17', limits=[0.35, 0.55])
# Diversity constraint
div1 = Diversity(word=5, types=2, scope=[d1])
div2 = Diversity(word=5, types=2, scope=[d9])
div3 = Diversity(word=5, types=2, scope=[d11])


design = tube_design(
    tubes=[tube1, tube2],
    hard_constraints=[sim1, sim2, sim3, sim4, sim5, div1, div2, div3],
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

with open("best_strands3.txt", "w") as f:
    f.write(f"Best trial ensemble defect: {best_result.defects.ensemble_defect:.4f}\n")
    for ts in target_strands:
        seq = str(best_result.to_analysis(ts))
        f.write(f"{ts.name}\t{seq}\n")
    print(f"Done! Saved best sequences to {f.name}")
