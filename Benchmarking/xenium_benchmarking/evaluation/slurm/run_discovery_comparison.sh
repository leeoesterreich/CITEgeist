#!/bin/bash
#SBATCH --job-name=discovery_cmp
#SBATCH --output=slurm_log/%x_%a.out
#SBATCH --error=slurm_log/%x_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-9
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Module 1-2 vs Leiden Discovery Comparison
# ============================================================================
#
# Array indices 0-4: spot resolution (regions 0-4)
# Array indices 5-9: single-cell resolution (regions 0-4)
#
# Usage:
#   sbatch run_discovery_comparison.sh           # Run all
#   SWEEP=1 sbatch run_discovery_comparison.sh   # Also run top_k sweep
# ============================================================================

SWEEP=${SWEEP:-0}

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
EVAL_SRC="${REPO_ROOT}/Benchmarking/xenium_benchmarking/evaluation/src"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/evaluation/results/discovery_comparison"

# Activate environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"
mkdir -p "${OUTPUT_DIR}"

# Determine region and resolution from array index
TASK_ID=${SLURM_ARRAY_TASK_ID}
if [ ${TASK_ID} -lt 5 ]; then
    REGION=${TASK_ID}
    RES_LEVEL="spot"
else
    REGION=$((TASK_ID - 5))
    RES_LEVEL="cell"
fi

echo "Region: ${REGION}, Resolution: ${RES_LEVEL}"

# Run Leiden baseline
echo "Running Leiden baseline..."
python "${EVAL_SRC}/leiden_baseline_comparison.py" \
    --region ${REGION} \
    --resolution-level ${RES_LEVEL} \
    --output-dir "${OUTPUT_DIR}"

# Run Module 1-2
echo "Running Module 1-2..."
python "${EVAL_SRC}/module12_discovery_runner.py" \
    --region ${REGION} \
    --resolution-level ${RES_LEVEL} \
    --output-dir "${OUTPUT_DIR}"

# Optional: top_k sweep (spot only to save compute)
if [ "${SWEEP}" = "1" ] && [ "${RES_LEVEL}" = "spot" ]; then
    echo "Running top_k sweep..."
    python "${EVAL_SRC}/module12_discovery_runner.py" \
        --region ${REGION} \
        --resolution-level ${RES_LEVEL} \
        --top-k-sweep \
        --output-dir "${OUTPUT_DIR}"
fi

echo "Done: Region ${REGION}, ${RES_LEVEL}"
