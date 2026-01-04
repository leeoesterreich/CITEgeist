#!/bin/bash
#SBATCH --job-name=autodiscovery_bench
#SBATCH --output=slurm_log/%x_%a.out
#SBATCH --error=slurm_log/%x_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-4

# ============================================================================
# Mode B: Full Auto-Discovery Benchmark
# ============================================================================
# Runs the complete CITEgeist pipeline with automatic profile discovery:
# 1. Module 1: Identify interesting markers
# 2. Module 2a: Pairwise colocalization analysis
# 3. Module 2b: Profile discovery via hierarchical clustering
# 4. Module 2c: Profile selection by spatial variance
# 5. Module 3: Cell proportion optimization
# 6. Evaluation against RNA-based ground truth (10 cell types)
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
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
INPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_granular_gt"
OUTPUT_DIR="${CITEGEIST_DIR}/output_autodiscovery"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

# Create directories if needed
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Auto-discovery benchmark started for region $SLURM_ARRAY_TASK_ID" | tee -a "${SLURM_LOG_DIR}/autodiscovery_runtime.log"

# ============================================================================
# Environment Setup
# ============================================================================

echo -e "${BLUE}Setting up environment...${RESET}"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

# Load Gurobi module for optimization
module load gurobi/12.0.3
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to load Gurobi module!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}Gurobi module loaded.${RESET}"

# ============================================================================
# Run Auto-Discovery Benchmark (Mode B)
# ============================================================================

REGION_ID=$SLURM_ARRAY_TASK_ID
echo -e "${YELLOW}Processing region ${REGION_ID} with AUTO-DISCOVERY pipeline${RESET}"
echo -e "${BLUE}Input directory: ${INPUT_DIR}${RESET}"
echo -e "${BLUE}Output directory: ${OUTPUT_DIR}${RESET}"

cd "${REPO_ROOT}"

$PYTHON "${CITEGEIST_DIR}/src/run_benchmark_autodiscovery.py" \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 4.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --fdr-alpha 0.05 \
    --top-k 3 \
    --min-marginal-gain 0.001

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Auto-discovery benchmark failed on region ${REGION_ID}!${RESET}" >&2
    exit 1
fi

echo -e "${GREEN}Auto-discovery benchmark completed for region ${REGION_ID}${RESET}"

# ============================================================================
# Timing
# ============================================================================

END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo -e "[${YELLOW}${END_TIME}${RESET}] Job completed for region ${REGION_ID}" | tee -a "${SLURM_LOG_DIR}/autodiscovery_runtime.log"
echo -e "[${GREEN}AUTODISCOVERY_RUNTIME${RESET}]: Region ${REGION_ID} took ${RUNTIME} seconds (${RUNTIME_MINUTES} minutes)" | tee -a "${SLURM_LOG_DIR}/autodiscovery_runtime.log"

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Region ${REGION_ID} Auto-Discovery Benchmark Complete${RESET}"
echo -e "${GREEN}======================================${RESET}"
