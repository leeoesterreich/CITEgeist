#!/bin/bash
#SBATCH --job-name=gex_no_l2
#SBATCH --output=logs/gex_no_l2_%A_%a.out
#SBATCH --error=logs/gex_no_l2_%A_%a.err
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

# Test with L2 regularization DISABLED (lambda_gex_reg=0.0)
# This should give concentrated gene allocation to most enriched cell types
# Instead of spreading uniformly as with lambda_gex_reg=0.01

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/gex_no_l2 \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0 \
    --lambda-gex-reg 0.0
