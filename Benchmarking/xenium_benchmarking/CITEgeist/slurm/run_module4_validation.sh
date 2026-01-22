#!/bin/bash
#SBATCH --job-name=module4_validation
#SBATCH --output=slurm_log/%x_%j.out
#SBATCH --error=slurm_log/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

# ============================================================================
# Module 4/5 Validation Pipeline
# ============================================================================
# Validates Module 4 (Anchored Program Discovery) and Module 5 (Cross-Sample
# Integration) on three data sources:
# 1. Single-cell Xenium RCC data
# 2. Manual deconvolved data (output_achievable_7)
# 3. Autodiscovery deconvolved data (output_autodiscovery)
#
# Prerequisites:
# - Achievable_7 and Autodiscovery benchmarks must be complete
# - Xenium single-cell data must be accessible
#
# Output:
# - Module 4/5 results for each source
# - Concordance analysis comparing all sources
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
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
OUTPUT_DIR="${CITEGEIST_DIR}/output_module4_validation"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

# Create directories
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Record start
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Module 4/5 Validation started"

# ============================================================================
# Environment Setup
# ============================================================================

echo -e "${BLUE}Setting up environment...${RESET}"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

# Load Gurobi (needed for some Module 4 operations)
module load gurobi/12.0.3
if [ $? -ne 0 ]; then
    echo -e "${RED}Warning: Failed to load Gurobi module${RESET}" >&2
fi

# ============================================================================
# Check Prerequisites
# ============================================================================

echo -e "${BLUE}Checking prerequisites...${RESET}"

# Check achievable_7 output
ACHIEVABLE_DIR="${CITEGEIST_DIR}/output_achievable_7"
if [ ! -d "${ACHIEVABLE_DIR}" ] || [ -z "$(ls -A ${ACHIEVABLE_DIR}/Xenium_region_*_pass1 2>/dev/null)" ]; then
    echo -e "${RED}Error: Achievable_7 output not found or incomplete${RESET}"
    echo -e "${YELLOW}Run xenium_benchmark_achievable_7.sh first${RESET}"
    exit 1
fi
echo -e "${GREEN}  Achievable_7 output: OK${RESET}"

# Check autodiscovery output
AUTODISCOVERY_DIR="${CITEGEIST_DIR}/output_autodiscovery"
if [ ! -d "${AUTODISCOVERY_DIR}" ] || [ -z "$(ls -A ${AUTODISCOVERY_DIR}/Xenium_region_*_pass1 2>/dev/null)" ]; then
    echo -e "${RED}Error: Autodiscovery output not found or incomplete${RESET}"
    echo -e "${YELLOW}Run run_autodiscovery_benchmark.sh first${RESET}"
    exit 1
fi
echo -e "${GREEN}  Autodiscovery output: OK${RESET}"

# Check Xenium single-cell data
XENIUM_DATA="/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"
if [ ! -f "${XENIUM_DATA}/cell_feature_matrix.h5" ]; then
    echo -e "${RED}Error: Xenium single-cell data not found${RESET}"
    exit 1
fi
echo -e "${GREEN}  Xenium single-cell data: OK${RESET}"

# ============================================================================
# Run Module 4/5 Validation
# ============================================================================

cd "${REPO_ROOT}"

echo -e "${YELLOW}Running Module 4/5 validation on all sources...${RESET}"

$PYTHON "${CITEGEIST_DIR}/src/run_module4_validation.py" \
    --source all \
    --all-regions \
    --K-programs 5 \
    --max-cells 50000 \
    --output-dir "${OUTPUT_DIR}"

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Module 4/5 validation failed!${RESET}" >&2
    exit 1
fi

echo -e "${GREEN}Module 4/5 validation completed${RESET}"

# ============================================================================
# Run Concordance Analysis
# ============================================================================

echo -e "${YELLOW}Running concordance analysis...${RESET}"

$PYTHON "${CITEGEIST_DIR}/src/analyze_concordance.py" \
    --validation-dir "${OUTPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}/concordance"

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Concordance analysis failed!${RESET}" >&2
    exit 1
fi

echo -e "${GREEN}Concordance analysis completed${RESET}"

# ============================================================================
# Summary
# ============================================================================

END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${END_TIME}${RESET}] Pipeline completed"

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Module 4/5 Validation Complete${RESET}"
echo -e "${GREEN}======================================${RESET}"
echo -e "Results saved to: ${OUTPUT_DIR}"
echo -e "Concordance results: ${OUTPUT_DIR}/concordance"
echo -e ""
echo -e "Key files:"
echo -e "  - validation_summary.json"
echo -e "  - concordance/concordance_summary.json"
echo -e "  - concordance/figures/*.png"
