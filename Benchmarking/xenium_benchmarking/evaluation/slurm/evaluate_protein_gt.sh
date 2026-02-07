#!/bin/bash
#SBATCH --job-name=eval_protein_gt
#SBATCH --output=slurm_log/%x.out
#SBATCH --error=slurm_log/%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Evaluate CITEgeist on Protein-Based Ground Truth
# ============================================================================
# This script evaluates CITEgeist (with negative marker redistribution)
# against protein-based single-cell ground truth.
#
# Ground Truth: Single-cell protein gating (8 cell types)
# Method: CITEgeist_protein_negmarkers
# ============================================================================

# Color codes
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
RED="\033[1;31m"
RESET="\033[0m"

# Directories
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
EVAL_DIR="${BENCH_DIR}/evaluation"
GT_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt"
OUTPUT_DIR="${EVAL_DIR}/results_protein_gt"
SLURM_LOG_DIR="${EVAL_DIR}/slurm/slurm_log"

mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

echo -e "${BLUE}Setting up environment...${RESET}"

# Activate environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to activate environment!${RESET}" >&2
    exit 1
fi

cd "${REPO_ROOT}"

echo -e "${YELLOW}Evaluating CITEgeist_protein_negmarkers against protein GT...${RESET}"
echo -e "${BLUE}Ground truth: ${GT_DIR}${RESET}"
echo -e "${BLUE}Output: ${OUTPUT_DIR}${RESET}"

# Run evaluation
python "${EVAL_DIR}/src/evaluate_all_methods.py" \
    --gt-dir "${GT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --methods CITEgeist_protein_negmarkers \
    --prefix Xenium \
    --n-regions 5

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Evaluation failed!${RESET}" >&2
    exit 1
fi

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Evaluation complete!${RESET}"
echo -e "${GREEN}Results: ${OUTPUT_DIR}${RESET}"
echo -e "${GREEN}======================================${RESET}"
