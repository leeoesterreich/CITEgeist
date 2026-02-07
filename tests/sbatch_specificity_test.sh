#!/bin/bash
#SBATCH --job-name=specificity_test
#SBATCH --output=slurm_log/specificity_%x_%a.out
#SBATCH --error=slurm_log/specificity_%x_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-9

# ============================================================================
# Test specificity-weighted profile discovery on all datasets
# ============================================================================
# Array mapping:
#   0-4:  simulated_high_seg replicate 0-4
#   5-9:  simulated_mixed replicate 0-4
#   10-14: xenium region 0-4
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

# Create log directory
mkdir -p "${REPO_ROOT}/tests/slurm_log"

# Set up environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

cd "${REPO_ROOT}"

# Determine dataset and replicate from array index
TASK_ID=$SLURM_ARRAY_TASK_ID

if [ $TASK_ID -lt 5 ]; then
    DATASET="simulated_high_seg"
    REPLICATE=$TASK_ID
elif [ $TASK_ID -lt 10 ]; then
    DATASET="simulated_mixed"
    REPLICATE=$((TASK_ID - 5))
else
    DATASET="xenium"
    REPLICATE=$((TASK_ID - 10))
fi

echo "Running specificity test on ${DATASET} replicate/region ${REPLICATE}"
echo "Task ID: ${TASK_ID}"

OUTPUT_DIR="${REPO_ROOT}/test_results/specificity_weighting"
mkdir -p "${OUTPUT_DIR}"

$PYTHON tests/test_specificity_weighting.py \
    --dataset "${DATASET}" \
    --replicate ${REPLICATE} \
    --output-dir "${OUTPUT_DIR}"

if [ $? -eq 0 ]; then
    echo "Test completed successfully for ${DATASET} ${REPLICATE}"
else
    echo "Test FAILED for ${DATASET} ${REPLICATE}" >&2
    exit 1
fi
