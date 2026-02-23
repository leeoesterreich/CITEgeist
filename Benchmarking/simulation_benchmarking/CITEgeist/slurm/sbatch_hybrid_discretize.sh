#!/bin/bash
#SBATCH --job-name=hybrid_discretize
#SBATCH --output=logs/hybrid_discretize_%A_%a.out
#SBATCH --error=logs/hybrid_discretize_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Hybrid approach:
# 1. Run continuous model to get optimal proportions (~0.90)
# 2. Discretize using nuclei counts: c_t = round(p_t * N)
# 3. Measure how much correlation we lose from discretization
python Benchmarking/simulation_benchmarking/CITEgeist/src/test_discretize_continuous.py \
    --replicate-id ${SLURM_ARRAY_TASK_ID} \
    --condition mixed \
    --mode dapi \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist/output_hybrid_discretize/mixed/dapi
