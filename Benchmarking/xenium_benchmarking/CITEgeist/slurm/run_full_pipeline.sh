#!/bin/bash
#SBATCH --job-name=xenium_full           # Job name
#SBATCH --output=slurm_log/%x_%a.out     # Standard output
#SBATCH --error=slurm_log/%x_%a.err      # Standard error
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=8                # CPU cores per task
#SBATCH --mem=128G                       # Memory per node
#SBATCH --time=12:00:00                  # Time limit
#SBATCH --cluster=htc                    # Cluster
#SBATCH --partition=htc                  # Partition
#SBATCH --array=0-4                      # 5 regions
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Full CITEgeist Pipeline (Modules 1-5) on Xenium Data
# ============================================================================
# Runs:
#   - Module 1-2: Preprocessing
#   - Module 3a: Cell proportion estimation
#   - Module 3b: Gene expression deconvolution
#   - Benchmarking: Proportion and GEX metrics
#   - Module 4: Anchored program discovery
# ============================================================================

# Directories
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
INPUT_DIR="${PSEUDOVISIUM_DIR}/data"
OUTPUT_DIR="${CITEGEIST_DIR}/output"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

# Create directories
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[${START_TIME}] Job started for region $SLURM_ARRAY_TASK_ID"

# ============================================================================
# Environment Setup
# ============================================================================

echo "Setting up environment..."

# Activate CITEgeist conda environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate CITEgeist environment!"
    exit 1
fi
echo "CITEgeist environment activated."

# Load Gurobi module
module load gurobi/12.0.3
if [ $? -ne 0 ]; then
    echo "Error: Failed to load Gurobi module!"
    exit 1
fi
echo "Gurobi module loaded."

# ============================================================================
# Run Full Pipeline
# ============================================================================

REGION_ID=$SLURM_ARRAY_TASK_ID
echo ""
echo "============================================"
echo "Processing Region ${REGION_ID}"
echo "============================================"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

cd "${REPO_ROOT}"

# Run full pipeline (Modules 1-4)
python "${CITEGEIST_DIR}/src/run_full_pipeline.py" \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 4.0

EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "Error: Pipeline failed on region ${REGION_ID} with exit code ${EXIT_CODE}"
    exit $EXIT_CODE
fi

# ============================================================================
# Timing
# ============================================================================

END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo ""
echo "============================================"
echo "[${END_TIME}] Region ${REGION_ID} completed"
echo "Runtime: ${RUNTIME} seconds (${RUNTIME_MINUTES} minutes)"
echo "============================================"

# Log runtime
echo "[${END_TIME}] Region ${REGION_ID}: ${RUNTIME}s (${RUNTIME_MINUTES}min)" >> "${SLURM_LOG_DIR}/runtime.log"
