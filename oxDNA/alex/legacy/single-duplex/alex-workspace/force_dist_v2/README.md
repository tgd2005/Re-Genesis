## How to navigate the force_dist repository

A little overview of the end goal of this set of simulations:
to characterise the force-distance behaviour of our DNA nanodevice, as a function of both stiffness and equilibrium length.
DNA assemblies and simulations in oxDNA are sensitive to failure because the process of using traps results in highly stressed energetic configurations of the system.
Moreover, this system requires many traps to assemble which in itself may result extreme configurations which also fail.

To resolve this problem, we set up the system in such a way that traps are incorporated/removed gradually, and allowing the system to relax at each metastable state.

Because of the large amount of files generated in these pipelines, I want to emphasize adhering to a specific naming convention that might not always make sense
but is crucial for keeping everything simple and easy to follow.

### Protocol for running simulations

Each step always consists of running an MC simulation followed by an MD simulation for relaxation purposes as per Lorenzo's documentation on relaxing initial DNA
configurations.

For external force, external.conf is only the new traps whilst external2.conf is the new traps + assembly traps. The only reason they are named like this is
because of the order in which I made them.

The steps given in order are:
1) Run input_postass_rel_MC, using the INITIAL_CONF_FILE and external2.conf -> runs a relaxation MC producing observables with the prefix "postass_relax_MC" or "postass_rel_MC"
2) Run input_postass_rel_MD, using postass_relax_MC_last.conf and external2.conf -> runs a relaxation MD producing observables with the prefix "postass_relax_MD" or "postass_rel_MD"
3) Run input_only_traps_MC, using postass_relax_MD_last.conf and external.conf -> observables have prefix "posttraps_relax_MC" or "posttraps_rel_MC"
4) Run input_only_traps_MD, using posttraps_relax_MD_last.conf and external.conf -> observables have prefix "posttraps_relax_MD" or "posttraps_rel_MD"

What hasn't been done yet, and would probably be wise to do so, is run one final MD simulation same as (4) but for a longer time span and dt for scientific robustness.
I propose that for this we use input file name "input_only_traps_long_MD" observables prefix "posttraps_long_MD" to avoid further confusion.

### Optional Observables

There are two types of optional observables:
PREFIX_dist.dat contains the following columns: {"step"  "distance between nt 15 and 16"  "mutual trap force vectors for nt 15 and 16"}
PREFIX_hb.dat contains information on whether nucleotides are hybridised in the system. Useful for evaluating system failures without oxview.

## For Sasha

Since you are going to be working with variable stiffness I propose you name your folders following "K_VALUE_r0_VALUE".
Also, when it comes to initial conditions, name your configuration files as so: "ass_relaxed_K_VALUE_r0_PREVVALUE". Follow my example.
