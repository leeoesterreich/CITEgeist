#!/bin/bash
#SBATCH --job-name=eval_rna_gt             # Job name
#SBATCH --output=slurm_log/%x.out          # Standard output
#SBATCH --error=slurm_log/%x.err           # Standard error
#SBATCH --ntasks=1                         # Number of tasks
#SBATCH --cpus-per-task=1                  # CPU cores per task
#SBATCH --mem=16G                          # Memory per node
#SBATCH --time=00:30:00                    # Time limit
#SBATCH --cluster=htc                      # Cluster
#SBATCH --partition=htc                    # Partition

# ============================================================================
# Evaluate CITEgeist Benchmark Results (RNA-based Ground Truth)
# ============================================================================
# This script evaluates CITEgeist predictions against RNA-based ground truth.
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
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
GT_DIR="${BENCH_DIR}/data_rna_gt"
PRED_DIR="${BENCH_DIR}/output_rna_gt"
OUTPUT_JSON="${BENCH_DIR}/results/rna_gt_evaluation.json"
SLURM_LOG_DIR="${BENCH_DIR}/slurm/slurm_log"

# Create directories
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "$(dirname ${OUTPUT_JSON})"

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

# ============================================================================
# Run Evaluation
# ============================================================================

echo -e "${BLUE}Evaluating benchmark results...${RESET}"
echo -e "${BLUE}Ground truth: ${GT_DIR}${RESET}"
echo -e "${BLUE}Predictions:  ${PRED_DIR}${RESET}"
echo -e "${BLUE}Output:       ${OUTPUT_JSON}${RESET}"

cd "${REPO_ROOT}"

python "${BENCH_DIR}/src/evaluate_benchmark.py" \
    --gt-dir "${GT_DIR}" \
    --pred-dir "${PRED_DIR}" \
    --n-regions 5 \
    --output "${OUTPUT_JSON}" \
    --prefix "Xenium"

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Evaluation failed!${RESET}" >&2
    exit 1
fi

# ============================================================================
# Done
# ============================================================================

END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${END_TIME}${RESET}] Evaluation completed"
echo -e "${GREEN}Results saved to: ${OUTPUT_JSON}${RESET}"
