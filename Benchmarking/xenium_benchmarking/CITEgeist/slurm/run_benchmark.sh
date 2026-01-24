#!/bin/bash
#SBATCH --job-name=citegeist_bench              # Job name
#SBATCH --output=slurm_log/%x_%a.out            # Standard output
#SBATCH --error=slurm_log/%x_%a.err             # Standard error
#SBATCH --ntasks=1                              # Number of tasks
#SBATCH --cpus-per-task=4                       # CPU cores per task
#SBATCH --mem=64G                               # Memory per node
#SBATCH --time=24:00:00                         # Time limit
#SBATCH --cluster=htc                           # Cluster
#SBATCH --partition=htc                         # Partition
#SBATCH --array=0-4                             # 5 regions (adjust for more)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Unified CITEgeist Xenium Benchmark
# ============================================================================
#
# Two modes:
# - manual: Achievable-7 hand-crafted profiles (fair benchmark baseline)
# - hierarchical: Autodiscovery with hierarchy learning (flagship feature)
#
# Usage:
#   sbatch run_benchmark.sh                    # Default: manual mode
#   MODE=hierarchical sbatch run_benchmark.sh  # Hierarchical autodiscovery
#
# ============================================================================

# Get mode from environment variable (default: manual)
MODE=${MODE:-manual}

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
INPUT_DIR="${PSEUDOVISIUM_DIR}/data_granular_gt"
OUTPUT_DIR="${CITEGEIST_DIR}/output/${MODE}"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

# Create directories if needed
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Job started for region $SLURM_ARRAY_TASK_ID (${MODE} mode)" | tee -a "${SLURM_LOG_DIR}/citegeist_bench_runtime.log"

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

# Load Gurobi module
module load gurobi/12.0.3
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to load Gurobi module!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}Gurobi module loaded.${RESET}"

# ============================================================================
# Run CITEgeist Benchmark
# ============================================================================

REGION_ID=$SLURM_ARRAY_TASK_ID
echo -e "${YELLOW}Processing region ${REGION_ID} with ${MODE} mode${RESET}"
echo -e "${BLUE}Input directory: ${INPUT_DIR}${RESET}"
echo -e "${BLUE}Output directory: ${OUTPUT_DIR}${RESET}"

# Change to repo root
cd "${REPO_ROOT}"

# Run benchmark with specified mode
echo -e "${BLUE}Running CITEgeist on region ${REGION_ID}...${RESET}"
python "${CITEGEIST_DIR}/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --mode ${MODE} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 4.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --run-gex

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: CITEgeist failed on region ${REGION_ID}!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}CITEgeist completed successfully for region ${REGION_ID}${RESET}"

# ============================================================================
# Cleanup and Timing
# ============================================================================

END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo -e "[${YELLOW}${END_TIME}${RESET}] Job completed for region ${REGION_ID}" | tee -a "${SLURM_LOG_DIR}/citegeist_bench_runtime.log"
echo -e "[${GREEN}RUNTIME${RESET}]: Region ${REGION_ID} (${MODE}) took ${RUNTIME} seconds (${RUNTIME_MINUTES} minutes)" | tee -a "${SLURM_LOG_DIR}/citegeist_bench_runtime.log"

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Region ${REGION_ID} ${MODE} benchmark complete${RESET}"
echo -e "${GREEN}======================================${RESET}"
