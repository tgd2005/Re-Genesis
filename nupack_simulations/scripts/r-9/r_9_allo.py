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
allo_length = 35  # length of the allosteric region in strand 2
linker_length = 5 # length of the all
clamp_length = 15  # length of the clamp domains in strands 1 and 2

# NOTICE: loop_length - rod_length / 2 should be an integer.
# blankregion = (loop_length - rod_length) / 2
# if not blankregion.is_integer():
    #raise ValueError("loop_length - rod_length must be an even number.")
#blankregion = int(blankregion)


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

# clamp region is temperarily removed for full hybridization
d1  = Domain(f'N{clamp_length}', name='Domain 1')   # N right clamp length
d2  = Domain("GAG", name="5'Overhang") # d2-d8 are the lettuce domains.
d3  = Domain("GCTCGGGACGAC", name="5'Stem")
d4  = Domain('AGGGATCAGGGATGTAGCGGCAGTGGGTGA', name='LettuceCore')
d5  = Domain("GTCGTCCCGAGC", name="3'Stem")
d6  = Domain("TGCAATCGTAA", name="3'Overhang")
d7  = Domain(f'N{clamp_length}', name='Domain 7')    # left clamp length
#d8 = Domain(f'N{blankregion}', name='Domain 8')    # loop - rod = unpaired region
#d9  = Domain(f'N{rod_length}', name='Domain 9')    # N20; Rod structure, variable length in simulations.
d10 = Domain(f'N{linker_length}', name='Domain 10')    # allosteric region 
# Strand 1 = Lettuce_S-9-1
Strand1 = TargetStrand(
    [d1, d2, d3, d4, d5, d6, d7],
    name='Lettuce_R-9-1'
)

# Strand 2 = Clamp loop structure
Strand2 = TargetStrand(
    [~d7,  d10, ~d1],
    name='allosteric strand'
)



lettuce_r9_structure = '.3(12.30)12.11'
# legacy structure for 2 strands: '(15(22.22(10.5)10.4)22(15+)15.54)15'
s1s2_complex_structure = (
    f'({clamp_length}'       # Strand 1  5' 
    + lettuce_r9_structure + # Lettuce
    f'({clamp_length}'       # Strand 1  3' 
    + '+' +                  # --- 链的分隔 ---
    f'){clamp_length}'       # Strand 2  5' 
    f'.{linker_length}'      # Linker
    f'){clamp_length}'       # Strand 2  3'
)

# complex when only lettuce and loop are interacting
C1 = TargetComplex(
    [Strand1, Strand2],
    structure=s1s2_complex_structure,
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
sim1 = Similarity([d1], f'S{clamp_length}', limits=[0.40, 0.55])
sim2 = Similarity([d7], f'S{clamp_length}', limits=[0.40, 0.55])
#sim3 = Similarity([d8], f'S{blankregion}', limits=[0.40, 0.55])
#sim4 = Similarity([d9], f'S{rod_length}', limits=[0.40, 0.55])
sim5 = Similarity([d10], f'S{linker_length}', limits=[0.40, 0.55])
# Diversity constraint
#div1 = Diversity(word=10, types=3, scope=[d1, d7, d10])
# Pattern constraint
pat1 = Pattern(['GGGG'], scope=[d1, d7, d10])  # No GGGG in variable domains


design = tube_design(
    tubes=[tube1],
    hard_constraints=[sim1, sim2, sim5, pat1],
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

output_dir = "output"
os.makedirs(output_dir, exist_ok=True) 


output_file = os.path.join(output_dir, "strands3_r91.txt")


try:
    with open(output_file, "w") as f:
        
        f.write("# NUPACK Design Result\n")
        f.write(f"# Target Structure: {s1s2_complex_structure}\n")
        f.write(f"# Best Trial Index: {best_idx}\n")
        f.write(f"# Ensemble Defect:  {best_result.defects.ensemble_defect:.6f}\n")
        f.write("-" * 40 + "\n\n")

     
        target_strands = [Strand1, Strand2]
        
        for ts in target_strands:
            
            seq = str(best_result.to_analysis(ts))
            
        
            f.write(f">{ts.name}\n")
            
            f.write(f"{seq}\n")
            
            
            print(f"[{ts.name}] saved.")
            print(f"Seq: {seq}\n")

    print(f"\n Success! Results saved to: {os.path.abspath(output_file)}")

except IOError as e:
    print(f"\n Error saving file: {e}")