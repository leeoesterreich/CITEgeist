#!/bin/bash
#SBATCH --job-name=xenium_granular          # Job name
#SBATCH --output=slurm_log/%x_%a.out        # Standard output
#SBATCH --error=slurm_log/%x_%a.err         # Standard error
#SBATCH --ntasks=1                          # Number of tasks
#SBATCH --cpus-per-task=4                   # CPU cores per task
#SBATCH --mem=64G                           # Memory per node
#SBATCH --time=24:00:00                     # Time limit (extended for GEX deconvolution)
#SBATCH --cluster=htc                       # Cluster
#SBATCH --partition=htc                     # Partition
#SBATCH --array=0-4                         # 5 regions
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Xenium Benchmarking - CITEgeist with GRANULAR 10-Cell-Type Ground Truth
# ============================================================================
# This script runs CITEgeist on pseudo-Visium spots using the UNSIMPLIFIED
# 10-cell-type RNA clustering for ground truth.
#
# This provides maximum granularity to highlight CITEgeist's proteomic
# integration advantage for distinguishing:
# - CD8+ T cells vs Proliferating T vs Mixed Immune
# - Myofibroblasts vs Stromal vs Vascular Stromal
# - Endothelial vs Vascular Stromal
#
# Cell types (10):
# 1. CD8+ T cells     - CD3E=378, CD8A=210
# 2. Macrophages      - CD68=430, CD163=88
# 3. Mixed Immune     - CD3E=142, CD8A=118, HLA-DR=142
# 4. Epithelial       - PanCK=39
# 5. Myofibroblasts   - alphaSMA=108, Vimentin=374
# 6. Stromal          - Vimentin+, mixed low
# 7. Endothelial      - CD31=168
# 8. B cells          - CD20=293, CD45RA=398
# 9. Proliferating T  - CD3E=679, PCNA=83
# 10. Vascular Stromal - CD31=53, Vimentin=209
# ============================================================================

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
INPUT_DIR="${PSEUDOVISIUM_DIR}/data_granular_gt"    # Granular 10-cell-type data
OUTPUT_DIR="${CITEGEIST_DIR}/output_granular"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

# Create directories if needed
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Job started for region $SLURM_ARRAY_TASK_ID" | tee -a "${SLURM_LOG_DIR}/xenium_granular_runtime.log"

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
# Run CITEgeist with GRANULAR 10-cell-type profiles
# ============================================================================

REGION_ID=$SLURM_ARRAY_TASK_ID
echo -e "${YELLOW}Processing region ${REGION_ID} with GRANULAR 10-cell-type ground truth${RESET}"
echo -e "${BLUE}Input directory: ${INPUT_DIR}${RESET}"
echo -e "${BLUE}Output directory: ${OUTPUT_DIR}${RESET}"

# Change to repo root
cd "${REPO_ROOT}"

# Run CITEgeist benchmark with GRANULAR cell profiles (10 cell types)
# Includes both cell proportion estimation AND gene expression deconvolution
echo -e "${BLUE}Running CITEgeist on region ${REGION_ID}...${RESET}"
python "${CITEGEIST_DIR}/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 4.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --use-granular-profiles \
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

echo -e "[${YELLOW}${END_TIME}${RESET}] Job completed for region ${REGION_ID}" | tee -a "${SLURM_LOG_DIR}/xenium_granular_runtime.log"
echo -e "[${GREEN}XENIUM_GRANULAR_RUNTIME${RESET}]: Region ${REGION_ID} took ${RUNTIME} seconds (${RUNTIME_MINUTES} minutes)" | tee -a "${SLURM_LOG_DIR}/xenium_granular_runtime.log"

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Region ${REGION_ID} GRANULAR benchmark complete${RESET}"
echo -e "${GREEN}======================================${RESET}"
