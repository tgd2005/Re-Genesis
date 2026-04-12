#!/bin/sh
#PBS -N oxDNA_Array_12nt
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=10GB
#PBS -J 1-10


MAIN_DIR="/rds/general/user/tc2523/home/Re-Genesis/Re-Genesis/oxDNA/Morris/force_insulation_sim/9nt_sim_0.0.1"
cd $MAIN_DIR


ID=${PBS_ARRAY_INDEX}


RUN_DIR="run_seed_${ID}"
mkdir -p $RUN_DIR


cp 9nt_input.txt 9nt-insulation-RELAXED.top 9nt-insulation-RELAXED.dat $RUN_DIR/
cp 9nt-insulation-RELAXED_force.txt 9_nt_hb_region.txt $RUN_DIR/


cd $RUN_DIR

# 环境变量设置
export OMP_NUM_THREADS=8


inputfile="9nt_input.txt"
oxdna_bin="/rds/general/user/tc2523/home/oxDNA/oxDNA2/oxDNA/build/bin/oxDNA"


runcommand="$oxdna_bin $inputfile \
    seed=$((ID + 10000)) \
    energy_file=energy_s${ID}.dat \
    trajectory_file=traj_s${ID}.dat \
    lastconf_file=last_s${ID}.conf \
    data_output_1.name=hb_list_s${ID}.dat \
    data_output_2.name=hb_distance_s${ID}.dat \
    data_output_3.name=hb_energy_s${ID}.dat "

# 7. 执行模拟
timeout 23h $runcommand