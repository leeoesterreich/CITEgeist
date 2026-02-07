#!/bin/bash
#SBATCH --job-name=compare_ach7_scres
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --output=slurm_log/compare_ach7_scresolve.out
#SBATCH --error=slurm_log/compare_ach7_scresolve.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

echo "=== Comparing All Methods on ACHIEVABLE-7 Including scResolve ==="
echo "Date: $(date)"
echo ""

# Activate environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation

python src/compare_achievable_7_with_scresolve.py

echo ""
echo "=== Comparison complete ==="
