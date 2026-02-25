#!/bin/bash
#SBATCH --job-name=gex_improved
#SBATCH --output=logs/gex_improved_%A_%a.out
#SBATCH --error=logs/gex_improved_%A_%a.err
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

# GEX deconvolution improvements test:
# 1. 75th percentile threshold (instead of 50th)
# 2. No smoothing (instead of 80/20)
# 3. Enrichment only (no proportion weighting)
#
# Uses the hybrid approach (continuous → discretize) as baseline
# but with improved GEX deconvolution

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/gex_improved \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0
