#!/bin/bash
#SBATCH --job-name=eval_compare
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/compare_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/compare_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

echo "=============================================="
echo "Full 5-Method Comparison (Protein GT)"
echo "=============================================="
echo "Start time: $(date)"

CITEGEIST_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"
export PATH="${CITEGEIST_ENV}/bin:$PATH"
source activate "${CITEGEIST_ENV}"

SRC_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src"

python "${SRC_DIR}/compare_all_methods.py"

echo ""
echo "End time: $(date)"
