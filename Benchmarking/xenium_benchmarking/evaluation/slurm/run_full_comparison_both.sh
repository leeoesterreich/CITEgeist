#!/bin/bash
#SBATCH --job-name=eval_both
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/compare_both_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/compare_both_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

echo "=============================================="
echo "Full Comparison: Proportions + GEX"
echo "=============================================="
echo "Start time: $(date)"

CITEGEIST_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"
export PATH="${CITEGEIST_ENV}/bin:$PATH"
source activate "${CITEGEIST_ENV}"

SRC_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src"

echo ""
echo "=========================================="
echo "PART 1: PROPORTIONS (5 methods)"
echo "=========================================="
python "${SRC_DIR}/compare_all_methods.py"

echo ""
echo "=========================================="
echo "PART 2: GEX DECONVOLUTION (3 methods)"
echo "=========================================="
python "${SRC_DIR}/compare_all_methods_gex.py"

echo ""
echo "End time: $(date)"
