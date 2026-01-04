#!/bin/bash
#SBATCH --job-name=compare_benchmarks
#SBATCH --output=slurm_log/%x_%j.out
#SBATCH --error=slurm_log/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

# ============================================================================
# Mode C: Compare Benchmark Results
# ============================================================================
# Compares CITEgeist performance across benchmarking modes:
# - Manual Profiles: Hand-crafted cell_profile_dict
# - Auto-Discovery: Module 1+2 discovered profiles
#
# Generates comparison figures and summary statistics.
# ============================================================================

# Color codes for output
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
RESET="\033[0m"

# Directories
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"

MANUAL_DIR="${CITEGEIST_DIR}/output_granular"
AUTODISCOVERY_DIR="${CITEGEIST_DIR}/output_autodiscovery"
OUTPUT_DIR="${CITEGEIST_DIR}/figures"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${SLURM_LOG_DIR}"

echo -e "${BLUE}Comparing benchmark results...${RESET}"
echo -e "${BLUE}Manual results: ${MANUAL_DIR}${RESET}"
echo -e "${BLUE}Auto-discovery results: ${AUTODISCOVERY_DIR}${RESET}"
echo -e "${BLUE}Output directory: ${OUTPUT_DIR}${RESET}"

# Set up environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

cd "${REPO_ROOT}"

$PYTHON "${CITEGEIST_DIR}/src/compare_benchmark_modes.py" \
    --manual-dir "${MANUAL_DIR}" \
    --autodiscovery-dir "${AUTODISCOVERY_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-regions 5

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Comparison failed!${RESET}" >&2
    exit 1
fi

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Benchmark Comparison Complete${RESET}"
echo -e "${GREEN}Figures saved to: ${OUTPUT_DIR}${RESET}"
echo -e "${GREEN}======================================${RESET}"
