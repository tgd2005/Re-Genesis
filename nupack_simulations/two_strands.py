from nupack import *


"""
Design two DNA strands for a target duplex using the NUPACK 4 Python API.

Requirements:
- Python 3 and the `nupack` package installed.

Inputs:
- No command-line arguments; all model, domains, and constraints are defined in this script.

Outputs:
- Prints ensemble defect for 100 design trials and the best trial’s strand sequences.
- Saves the best Strand 1 and Strand 2 sequences to 'best_strands.txt'.
"""


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

d1  = Domain('N15', name='Domain 1')   # N15
d2  = Domain('AACGTGCTCAAGGTGTCTTAGT', name='Domain 2')
d3  = Domain('AGGGATGATGCGGCAGTGGGCT', name='Domain 3')
d4  = Domain('TCGCAGAAGA', name='Domain 4')
d5  = Domain('GTGTT', name='Domain 5')
d6  = Domain('TCTTCTGCGA', name='Domain 6')
d7  = Domain('GGGG', name='Domain 7')
d8  = Domain('ACTAAGACACCTTGAGCACGTT', name='Domain 8')
d9  = Domain('N15', name='Domain 9')    # N15
d10  = Domain('N54', name='Domain 10')    # N54


Strand1 = TargetStrand(
    [d1, d2, d3, d4, d5, d6, d7, d8, d9],
    name='Strand 1'
)

Strand2 = TargetStrand(
    [~d1, d10, ~d9],
    name='Strand 2'
)

complex_structure = '(15(22.22(10.5)10.4)22(15+)15.54)15'

C1 = TargetComplex(
    [Strand1, Strand2],
    structure=complex_structure,
    name='Complex 1'
)

#print(f"The number of nucleotides in the target structure:{C1.nt(model.alphabet)}")

tube1 = TargetTube(
    on_targets={C1: 1e-6},
    off_targets=SetSpec(max_size=2),
    name='Tube 1'
)

# ---------- Design job (tube_design) ----------

#GC content constraint
sim1 = Similarity([d1], 'S15', limits=[0.45, 0.55])
sim2 = Similarity([d9], 'S15', limits=[0.45, 0.55])
sim3 = Similarity([d10], 'S54', limits=[0.45, 0.55])


design = tube_design(
    tubes=[tube1],
    hard_constraints=[sim1, sim2, sim3],
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



target_strands = [Strand1, Strand2]

print("Best strand sequences:")
for ts in target_strands:
    designed_strand = best_result.to_analysis(ts)
    seq = str(designed_strand)
    print(f"{ts.name}: {seq}")

with open("best_strands.txt", "w") as f:
    for ts in target_strands:
        seq = str(best_result.to_analysis(ts))
        f.write(f"{ts.name}\t{seq}\n")
    print(f"Done! Saved best sequences to {f.name}")
