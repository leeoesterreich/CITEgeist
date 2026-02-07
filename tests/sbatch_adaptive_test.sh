#!/bin/bash
#SBATCH --job-name=adaptive_coloc
#SBATCH --output=slurm_log/adaptive_%x_%a.out
#SBATCH --error=slurm_log/adaptive_%x_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-39

# ============================================================================
# Test adaptive colocalization enhancements on simulated data
# ============================================================================
# Array mapping (40 jobs = 4 configs × 2 datasets × 5 replicates):
#
#   Config 0: baseline (no enhancements)
#   Config 1: expression_fallback only
#   Config 2: multiscale only
#   Config 3: both enhancements
#
#   Dataset 0: high_seg (strong spatial signal)
#   Dataset 1: mixed (weak spatial signal)
#
#   Jobs 0-4:   baseline, high_seg, rep 0-4
#   Jobs 5-9:   baseline, mixed, rep 0-4
#   Jobs 10-14: expr_fallback, high_seg, rep 0-4
#   Jobs 15-19: expr_fallback, mixed, rep 0-4
#   Jobs 20-24: multiscale, high_seg, rep 0-4
#   Jobs 25-29: multiscale, mixed, rep 0-4
#   Jobs 30-34: both, high_seg, rep 0-4
#   Jobs 35-39: both, mixed, rep 0-4
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

# Create log directory
mkdir -p "${REPO_ROOT}/tests/slurm_log"

# Set up environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

cd "${REPO_ROOT}"

# Parse array index
TASK_ID=$SLURM_ARRAY_TASK_ID

# Determine config (0-3), dataset (0-1), replicate (0-4)
# Each config has 10 jobs (2 datasets × 5 replicates)
CONFIG=$((TASK_ID / 10))
REMAINDER=$((TASK_ID % 10))
DATASET_ID=$((REMAINDER / 5))
REPLICATE=$((REMAINDER % 5))

# Config names
case $CONFIG in
    0)
        CONFIG_NAME="baseline"
        USE_EXPR_FALLBACK="false"
        USE_MULTISCALE="false"
        ;;
    1)
        CONFIG_NAME="expr_fallback"
        USE_EXPR_FALLBACK="true"
        USE_MULTISCALE="false"
        ;;
    2)
        CONFIG_NAME="multiscale"
        USE_EXPR_FALLBACK="false"
        USE_MULTISCALE="true"
        ;;
    3)
        CONFIG_NAME="both"
        USE_EXPR_FALLBACK="true"
        USE_MULTISCALE="true"
        ;;
esac

# Dataset names
if [ $DATASET_ID -eq 0 ]; then
    DATASET="high_seg"
else
    DATASET="mixed"
fi

echo "============================================"
echo "Adaptive Colocalization Test"
echo "============================================"
echo "Task ID: ${TASK_ID}"
echo "Config: ${CONFIG_NAME} (${CONFIG})"
echo "Dataset: ${DATASET} (${DATASET_ID})"
echo "Replicate: ${REPLICATE}"
echo "Expression Fallback: ${USE_EXPR_FALLBACK}"
echo "Multi-scale: ${USE_MULTISCALE}"
echo "============================================"

OUTPUT_DIR="${REPO_ROOT}/test_results/adaptive_colocalization/${CONFIG_NAME}"
mkdir -p "${OUTPUT_DIR}"

# Run the test
$PYTHON tests/test_adaptive_benchmark.py \
    --dataset "${DATASET}" \
    --replicate ${REPLICATE} \
    --config "${CONFIG_NAME}" \
    --use-expr-fallback "${USE_EXPR_FALLBACK}" \
    --use-multiscale "${USE_MULTISCALE}" \
    --output-dir "${OUTPUT_DIR}"

if [ $? -eq 0 ]; then
    echo "Test completed successfully"
else
    echo "Test FAILED" >&2
    exit 1
fi
