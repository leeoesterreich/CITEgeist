#!/bin/bash
#SBATCH --job-name=profile_disc_eval
#SBATCH --output=slurm_log/%x_%a.out
#SBATCH --error=slurm_log/%x_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-4

# ============================================================================
# Mode A: Profile Discovery Evaluation
# ============================================================================
# Tests whether CITEgeist's Module 1+2 correctly discovers and groups expected
# marker combinations from Xenium protein data.
#
# Evaluates:
# - Profile Recovery Rate: How many expected groupings were discovered?
# - Marker Grouping Accuracy: Are expected markers grouped together?
# - Profile Purity: Are discovered profiles free of forbidden combinations?
# - Coverage: What fraction of expected markers appear in profiles?
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
OUTPUT_DIR="${CITEGEIST_DIR}/output_profile_discovery"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

# Create directories if needed
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Profile discovery evaluation started for region $SLURM_ARRAY_TASK_ID" | tee -a "${SLURM_LOG_DIR}/profile_disc_runtime.log"

# ============================================================================
# Environment Setup
# ============================================================================

echo -e "${BLUE}Setting up environment...${RESET}"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

# ============================================================================
# Run Profile Discovery Evaluation (Mode A)
# ============================================================================

REGION_ID=$SLURM_ARRAY_TASK_ID
echo -e "${YELLOW}Processing region ${REGION_ID} for profile discovery evaluation${RESET}"
echo -e "${BLUE}Input directory: ${INPUT_DIR}${RESET}"
echo -e "${BLUE}Output directory: ${OUTPUT_DIR}${RESET}"

cd "${REPO_ROOT}"

$PYTHON "${CITEGEIST_DIR}/src/evaluate_profile_discovery.py" \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --fdr-alpha 0.05 \
    --top-k 3 \
    --n-permutations 199

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Profile discovery evaluation failed on region ${REGION_ID}!${RESET}" >&2
    exit 1
fi

echo -e "${GREEN}Profile discovery evaluation completed for region ${REGION_ID}${RESET}"

# ============================================================================
# Timing
# ============================================================================

END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo -e "[${YELLOW}${END_TIME}${RESET}] Job completed for region ${REGION_ID}" | tee -a "${SLURM_LOG_DIR}/profile_disc_runtime.log"
echo -e "[${GREEN}PROFILE_DISC_RUNTIME${RESET}]: Region ${REGION_ID} took ${RUNTIME} seconds (${RUNTIME_MINUTES} minutes)" | tee -a "${SLURM_LOG_DIR}/profile_disc_runtime.log"

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Region ${REGION_ID} Profile Discovery Evaluation Complete${RESET}"
echo -e "${GREEN}======================================${RESET}"
