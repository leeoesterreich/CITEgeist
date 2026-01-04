#!/bin/bash
#SBATCH --job-name=eval_all_methods
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/eval_all_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/eval_all_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Unified Evaluation: Compare all deconvolution methods on Xenium benchmark
# Methods: CITEgeist, Cell2Location, Tangram, RCTD, Seurat

set -e

# Directories (use absolute paths - SCRIPT_DIR resolves incorrectly in SLURM)
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
EVAL_DIR="${BASE_DIR}/evaluation"
SRC_DIR="${EVAL_DIR}/src"
GT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_rna_gt"
OUTPUT_DIR="${EVAL_DIR}/results"

echo "=============================================="
echo "Unified Xenium Benchmark Evaluation"
echo "=============================================="
echo "Start time: $(date)"
echo ""
echo "Comparing methods:"
echo "  - CITEgeist (proportions + GEX)"
echo "  - Cell2Location (proportions + GEX)"
echo "  - Tangram (proportions + GEX)"
echo "  - RCTD (proportions only)"
echo "  - Seurat (proportions only)"
echo ""
echo "Base directory: ${BASE_DIR}"
echo "Ground truth: ${GT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

mkdir -p "${OUTPUT_DIR}"

# Use CITEgeist_env (has scipy, pandas, matplotlib)
ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"
MINICONDA_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/miniconda3"
export PATH="${ENV_PATH}/bin:${MINICONDA_PATH}/bin:$PATH"
source activate "${ENV_PATH}"

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Check which methods have results available
echo "Checking available results..."
for method in CITEgeist Cell2Location Tangram RCTD Seurat; do
    method_dir="${BASE_DIR}/${method}/output_rna_gt"
    if [ -d "${method_dir}" ] && [ "$(ls -A ${method_dir} 2>/dev/null)" ]; then
        echo "  ${method}: AVAILABLE"
    else
        echo "  ${method}: NOT AVAILABLE"
    fi
done
echo ""

# Run unified evaluation
python "${SRC_DIR}/evaluate_all_methods.py" \
    --base-dir "${BASE_DIR}" \
    --gt-dir "${GT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-regions 5 \
    --prefix "Xenium"

echo ""
echo "=============================================="
echo "Evaluation Complete!"
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
