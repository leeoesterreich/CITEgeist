#!/bin/bash
#SBATCH --job-name=gex_kl_only
#SBATCH --output=logs/gex_kl_only_%A_%a.out
#SBATCH --error=logs/gex_kl_only_%A_%a.err
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

# KL-ONLY variant: KL regularization without module enrichment
# Isolates the effect of replacing L2 with KL regularization

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/gex_kl_only \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0 \
    --use-kl-regularization \
    --kl-temperature 0.3 \
    --lambda-kl 0.1
