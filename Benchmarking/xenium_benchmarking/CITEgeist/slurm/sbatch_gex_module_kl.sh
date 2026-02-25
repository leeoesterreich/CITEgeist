#!/bin/bash
#SBATCH --job-name=gex_module_kl
#SBATCH --output=logs/gex_module_kl_%A_%a.out
#SBATCH --error=logs/gex_module_kl_%A_%a.err
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

# FULL variant: Module-aware enrichment + KL regularization
# Expected: Best GEX deconvolution performance

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/gex_module_kl \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0 \
    --use-module-enrichment \
    --module-weight 0.5 \
    --use-kl-regularization \
    --kl-temperature 0.3 \
    --lambda-kl 0.1
