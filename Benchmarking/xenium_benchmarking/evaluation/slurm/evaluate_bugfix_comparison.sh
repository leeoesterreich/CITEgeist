#!/bin/bash
#SBATCH --job-name=eval_bugfix
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/eval_bugfix_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/eval_bugfix_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Compare pre-fix vs post-fix CITEgeist benchmark results
# ============================================================================

set -e

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BASE_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
EVAL_SRC="${BASE_DIR}/evaluation/src"
GT_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_granular_gt"
RESULTS_DIR="${BASE_DIR}/evaluation/results_bugfix_comparison"

mkdir -p "${RESULTS_DIR}"

# Environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "=============================================="
echo "Bug Fix Comparison: Pre-fix vs Post-fix"
echo "=============================================="
echo ""

# Evaluate pre-fix results (output_achievable_7)
echo "--- Evaluating PRE-FIX results (output_achievable_7) ---"
python "${EVAL_SRC}/evaluate_benchmark.py" \
    --gt-dir "${GT_DIR}" \
    --pred-dir "${BASE_DIR}/CITEgeist/output_achievable_7" \
    --n-regions 5 \
    --prefix "Xenium" \
    --output "${RESULTS_DIR}/pre_fix_results.json"

echo ""
echo "--- Evaluating POST-FIX results (output/manual) ---"
python "${EVAL_SRC}/evaluate_benchmark.py" \
    --gt-dir "${GT_DIR}" \
    --pred-dir "${BASE_DIR}/CITEgeist/output/manual" \
    --n-regions 5 \
    --prefix "Xenium" \
    --output "${RESULTS_DIR}/post_fix_results.json"

echo ""
echo "=============================================="
echo "Comparison Complete!"
echo "=============================================="
echo "Pre-fix results: ${RESULTS_DIR}/pre_fix_results.json"
echo "Post-fix results: ${RESULTS_DIR}/post_fix_results.json"
