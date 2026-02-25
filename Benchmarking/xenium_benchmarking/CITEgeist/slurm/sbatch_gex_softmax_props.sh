#!/bin/bash
#SBATCH --job-name=gex_softmax_props
#SBATCH --output=logs/gex_softmax_props_%A_%a.out
#SBATCH --error=logs/gex_softmax_props_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# SOFTMAX PROPS variant: Apply softmax to proportions (sharpen allocation)
# WITHOUT KL penalty - just L2 regularization
# This is the user's simpler idea: concentrate on dominant cell types

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/gex_softmax_props \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0 \
    --kl-temperature 0.3
