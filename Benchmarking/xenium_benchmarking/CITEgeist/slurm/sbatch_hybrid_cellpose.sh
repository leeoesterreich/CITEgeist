#!/bin/bash
#SBATCH --job-name=hybrid_xenium
#SBATCH --output=logs/hybrid_xenium_%A_%a.out
#SBATCH --error=logs/hybrid_xenium_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=04:00:00
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
# 1. Run continuous model (CLR preprocessing, QP optimization)
# 2. Discretize using nuclei counts from Cellpose
# 3. Achieves ~94% of continuous model quality with integer cell counts

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_cellpose \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0
