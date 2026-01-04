#!/bin/bash
#SBATCH --job-name=eval_rna_gt             # Job name
#SBATCH --output=slurm_log/%x.out          # Standard output
#SBATCH --error=slurm_log/%x.err           # Standard error
#SBATCH --ntasks=1                         # Number of tasks
#SBATCH --cpus-per-task=1                  # CPU cores per task
#SBATCH --mem=16G                          # Memory per node
#SBATCH --time=01:00:00                    # Time limit (extended for GEX)
#SBATCH --cluster=htc                      # Cluster
#SBATCH --partition=htc                    # Partition

# ============================================================================
# Evaluate CITEgeist Benchmark Results (RNA-based Ground Truth)
# ============================================================================
# This script evaluates CITEgeist predictions against RNA-based ground truth.
#
# Evaluates BOTH:
# 1. Cell type proportions (JSD, RMSE, MAE, Pearson)
# 2. Gene expression deconvolution (RMSE, NRMSE, MAE on log1p data)
#
# Reference:
#   Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
#   Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
#   https://doi.org/10.1186/s12859-025-06044-0
# ============================================================================

# Color codes for output
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
RED="\033[1;31m"
RESET="\033[0m"

# Directories
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
EVALUATION_DIR="${BENCH_DIR}/evaluation"
GT_DIR="${PSEUDOVISIUM_DIR}/data_rna_gt"
PRED_DIR="${CITEGEIST_DIR}/output_rna_gt"
PROP_OUTPUT_JSON="${BENCH_DIR}/results/rna_gt_proportions_evaluation.json"
GEX_OUTPUT_JSON="${BENCH_DIR}/results/rna_gt_gex_evaluation.json"
SLURM_LOG_DIR="${EVALUATION_DIR}/slurm/slurm_log"

# Create directories
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "$(dirname ${PROP_OUTPUT_JSON})"

# Record start time
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Evaluation started"

# ============================================================================
# Environment Setup
# ============================================================================

echo -e "${BLUE}Setting up environment...${RESET}"

# Activate CITEgeist conda environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to activate CITEgeist environment!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}CITEgeist environment activated.${RESET}"

cd "${REPO_ROOT}"

# ============================================================================
# Part 1: Evaluate Cell Type Proportions
# ============================================================================

echo -e "\n${BLUE}============================================${RESET}"
echo -e "${BLUE}Part 1: Cell Type Proportion Evaluation${RESET}"
echo -e "${BLUE}============================================${RESET}"
echo -e "${BLUE}Ground truth: ${GT_DIR}/ground_truth${RESET}"
echo -e "${BLUE}Predictions:  ${PRED_DIR}${RESET}"
echo -e "${BLUE}Output:       ${PROP_OUTPUT_JSON}${RESET}"

python "${EVALUATION_DIR}/src/evaluate_benchmark.py" \
    --gt-dir "${GT_DIR}" \
    --pred-dir "${PRED_DIR}" \
    --n-regions 5 \
    --output "${PROP_OUTPUT_JSON}" \
    --prefix "Xenium"

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Proportion evaluation failed!${RESET}" >&2
    # Continue to GEX evaluation anyway
else
    echo -e "${GREEN}Proportion evaluation completed${RESET}"
fi

# ============================================================================
# Part 2: Evaluate Gene Expression Deconvolution
# ============================================================================

echo -e "\n${BLUE}============================================${RESET}"
echo -e "${BLUE}Part 2: Gene Expression Evaluation${RESET}"
echo -e "${BLUE}============================================${RESET}"
echo -e "${BLUE}Ground truth: ${GT_DIR}/ground_truth_gex${RESET}"
echo -e "${BLUE}Predictions:  ${PRED_DIR}/*_pass1/layers${RESET}"
echo -e "${BLUE}Output:       ${GEX_OUTPUT_JSON}${RESET}"

# Check if GEX ground truth exists
if [ -d "${GT_DIR}/ground_truth_gex" ]; then
    python "${EVALUATION_DIR}/src/evaluate_gex.py" \
        --gt-dir "${GT_DIR}" \
        --pred-dir "${PRED_DIR}" \
        --n-regions 5 \
        --output "${GEX_OUTPUT_JSON}" \
        --prefix "Xenium"

    if [ $? -ne 0 ]; then
        echo -e "${YELLOW}Warning: GEX evaluation failed (may not have GEX predictions yet)${RESET}" >&2
    else
        echo -e "${GREEN}GEX evaluation completed${RESET}"
    fi
else
    echo -e "${YELLOW}Warning: GEX ground truth not found at ${GT_DIR}/ground_truth_gex${RESET}"
    echo -e "${YELLOW}Run generate_gex_gt.sh first to generate GEX ground truth${RESET}"
fi

# ============================================================================
# Done
# ============================================================================

END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "\n[${YELLOW}${END_TIME}${RESET}] All evaluations completed"
echo -e "${GREEN}Proportion results: ${PROP_OUTPUT_JSON}${RESET}"
echo -e "${GREEN}GEX results:        ${GEX_OUTPUT_JSON}${RESET}"
