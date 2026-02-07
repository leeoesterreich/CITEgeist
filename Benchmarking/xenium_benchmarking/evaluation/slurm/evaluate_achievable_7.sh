#!/bin/bash
#SBATCH --job-name=eval_achievable_7
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/eval_achievable_7_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/eval_achievable_7_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Unified Evaluation: Compare all methods on ACHIEVABLE-7 benchmark
# ============================================================================
# This evaluates all methods against the 7-type achievable ground truth.
# GT types are collapsed: 10 -> 7 (Myofibroblasts+Stromal -> Fibroblasts)
#
# Methods evaluated:
# - CITEgeist_achievable_7: Manual 7-type profiles
# - CITEgeist_autodiscovery_achievable_7: Autodiscovered profiles mapped to 7
# - Cell2Location_achievable_7: Collapse 10-type predictions to 7
# - RCTD_achievable_7: Collapse 10-type predictions to 7
# - Tangram_achievable_7: Collapse 10-type predictions to 7
# - Seurat_achievable_7: Collapse 10-type predictions to 7
# ============================================================================

set -e

# Directories (use absolute paths - SCRIPT_DIR resolves incorrectly in SLURM)
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
EVAL_DIR="${BASE_DIR}/evaluation"
SRC_DIR="${EVAL_DIR}/src"
GT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt"
OUTPUT_DIR="${EVAL_DIR}/results_achievable_7"

echo "=============================================="
echo "ACHIEVABLE-7 Benchmark Evaluation"
echo "=============================================="
echo "Start time: $(date)"
echo ""
echo "Evaluating with 7-type achievable ground truth:"
echo "  - Myofibroblasts + Stromal -> Fibroblasts (both VIM+)"
echo "  - Vascular Stromal -> Endothelial (CD31+)"
echo "  - CD8+ T + Proliferating T -> CD8+ T (both CD3E+)"
echo "  - Mixed Immune -> CD4+ T (HLA-DR+ T)"
echo ""
echo "Methods:"
echo "  - CITEgeist_achievable_7 (manual 7-type profiles)"
echo "  - CITEgeist_autodiscovery_achievable_7 (autodiscovered -> 7 mapping)"
echo "  - Cell2Location_achievable_7 (10-type -> 7-type collapse)"
echo "  - RCTD_achievable_7 (10-type -> 7-type collapse)"
echo "  - Tangram_achievable_7 (10-type -> 7-type collapse)"
echo "  - Seurat_achievable_7 (10-type -> 7-type collapse)"
echo ""
echo "Base directory: ${BASE_DIR}"
echo "Ground truth: ${GT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${EVAL_DIR}/slurm/slurm_log"

# Use CITEgeist_env (has scipy, pandas, matplotlib)
ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"
MINICONDA_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/miniconda3"
export PATH="${ENV_PATH}/bin:${MINICONDA_PATH}/bin:$PATH"
source activate "${ENV_PATH}"

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Check which method outputs are available
echo "Checking available results..."
for method in CITEgeist Cell2Location Tangram RCTD Seurat; do
    if [ "${method}" = "CITEgeist" ]; then
        # Check both achievable_7 and autodiscovery outputs
        achievable_dir="${BASE_DIR}/${method}/output_achievable_7"
        autodiscovery_dir="${BASE_DIR}/${method}/output_autodiscovery"
        if [ -d "${achievable_dir}" ] && [ "$(ls -A ${achievable_dir} 2>/dev/null)" ]; then
            echo "  ${method}_achievable_7: AVAILABLE"
        else
            echo "  ${method}_achievable_7: NOT AVAILABLE (run xenium_benchmark_achievable_7.sh first)"
        fi
        if [ -d "${autodiscovery_dir}" ] && [ "$(ls -A ${autodiscovery_dir} 2>/dev/null)" ]; then
            echo "  ${method}_autodiscovery_achievable_7: AVAILABLE"
        else
            echo "  ${method}_autodiscovery_achievable_7: NOT AVAILABLE"
        fi
    else
        method_dir="${BASE_DIR}/${method}/output_granular"
        if [ -d "${method_dir}" ] && [ "$(ls -A ${method_dir} 2>/dev/null)" ]; then
            echo "  ${method}_achievable_7: AVAILABLE (will collapse 10->7)"
        else
            echo "  ${method}_achievable_7: NOT AVAILABLE"
        fi
    fi
done
echo ""

# Run evaluation for 7-type achievable methods
python "${SRC_DIR}/evaluate_all_methods.py" \
    --base-dir "${BASE_DIR}" \
    --gt-dir "${GT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-regions 14 \
    --prefix "Xenium" \
    --methods \
        CITEgeist_achievable_7 \
        CITEgeist_autodiscovery_achievable_7 \
        Cell2Location_achievable_7 \
        RCTD_achievable_7 \
        Tangram_achievable_7 \
        Seurat_achievable_7

echo ""
echo "=============================================="
echo "ACHIEVABLE-7 Evaluation Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/"
echo ""
echo "Key files:"
echo "  - ${OUTPUT_DIR}/full_results.json (complete results)"
echo "  - ${OUTPUT_DIR}/method_summary.csv (summary by method)"
echo "  - ${OUTPUT_DIR}/comparison_table.csv (all results)"
echo "  - ${OUTPUT_DIR}/method_comparison.png (comparison figure)"
