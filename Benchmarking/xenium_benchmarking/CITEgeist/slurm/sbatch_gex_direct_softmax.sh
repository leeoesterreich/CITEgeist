#!/bin/bash
#SBATCH --job-name=gex_direct_softmax
#SBATCH --output=logs/gex_direct_softmax_%A_%a.out
#SBATCH --error=logs/gex_direct_softmax_%A_%a.err
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

# DIRECT SOFTMAX: Skip Gurobi entirely, use direct allocation
# This is what worked in the previous analysis:
# - per_gene_r: 0.288 (+5.5% vs baseline)
# - variance_ratio: 1.53 (-49% vs baseline's 3.02)
#
# Key: kl_temperature < 1.0 AND use_kl_regularization=False
# triggers direct softmax mode (no Gurobi, no L2)

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/gex_direct_softmax \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0 \
    --kl-temperature 0.3
