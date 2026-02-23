#!/bin/bash
#SBATCH --job-name=discrete_global_highseg_96gb
#SBATCH --output=logs/discrete_global_highseg_96gb_%A_%a.out
#SBATCH --error=logs/discrete_global_highseg_96gb_%A_%a.err
#SBATCH --array=0,1,2
#SBATCH --time=02:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id ${SLURM_ARRAY_TASK_ID} \
    --condition high_seg \
    --mode dapi \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_global/high_seg/dapi \
    --global-solve \
    --global-time-limit 600 \
    --global-mip-gap 0.05 \
    --max-em-iterations 10
