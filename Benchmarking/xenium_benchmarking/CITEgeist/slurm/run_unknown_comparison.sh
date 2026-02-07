#!/bin/bash
#SBATCH --job-name=citegeist_unknown_cmp
#SBATCH --output=slurm_log/%x_%a.out
#SBATCH --error=slurm_log/%x_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-9
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Unknown Fraction Comparison: Xenium Benchmark
# ============================================================================
# Array 0-4: WITH unknown (5 regions)
# Array 5-9: WITHOUT unknown (5 regions)
# ============================================================================

# Determine condition and region from array task ID
if [ $SLURM_ARRAY_TASK_ID -lt 5 ]; then
    REGION_ID=$SLURM_ARRAY_TASK_ID
    CONDITION="with_unknown"
    NO_UNKNOWN_FLAG=""
else
    REGION_ID=$((SLURM_ARRAY_TASK_ID - 5))
    CONDITION="no_unknown"
    NO_UNKNOWN_FLAG="--no-unknown"
fi

# Directories
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
INPUT_DIR="${PSEUDOVISIUM_DIR}/data_granular_gt"
OUTPUT_DIR="${CITEGEIST_DIR}/output/unknown_comparison/${CONDITION}"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

START_TIMESTAMP=$(date +%s)
echo "[$(date)] Region ${REGION_ID} (${CONDITION}) started"

# Environment Setup
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

echo "Running CITEgeist region ${REGION_ID} (${CONDITION})..."
python "${CITEGEIST_DIR}/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --mode manual \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 250.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --run-gex \
    ${NO_UNKNOWN_FLAG}

EXIT_CODE=$?
END_TIMESTAMP=$(date +%s)
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))

echo "[$(date)] Region ${REGION_ID} (${CONDITION}) completed (exit=${EXIT_CODE}, ${RUNTIME}s)"
exit $EXIT_CODE
