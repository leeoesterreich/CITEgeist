#!/bin/bash
#SBATCH --job-name=xenium_protein_negmarkers    # Job name
#SBATCH --output=slurm_log/%x_%a.out            # Standard output
#SBATCH --error=slurm_log/%x_%a.err             # Standard error
#SBATCH --ntasks=1                              # Number of tasks
#SBATCH --cpus-per-task=4                       # CPU cores per task
#SBATCH --mem=64G                               # Memory per node
#SBATCH --time=4:00:00                          # Time limit
#SBATCH --cluster=htc                           # Cluster
#SBATCH --partition=htc                         # Partition
#SBATCH --array=0-4                             # 5 regions
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Xenium Benchmarking - CITEgeist with Protein-Based GT + Negative Markers
# ============================================================================
# This script runs CITEgeist on pseudo-Visium spots using:
# 1. Protein-based ground truth (from create_protein_gt.py)
# 2. Negative marker redistribution (post-hoc correction)
#
# Cell types (8):
# 1. B cells         - CD20+
# 2. CD4+ T cells    - CD3E+ CD4+ CD8A-
# 3. CD8+ T cells    - CD3E+ CD8A+
# 4. Macrophages     - CD68+ CD3E-
# 5. Endothelial     - CD31+ CD68- CD3E-
# 6. Epithelial      - PanCK+ or E-Cadherin high
# 7. Myofibroblasts  - alphaSMA high
# 8. Stromal         - Vimentin+ (remaining)
#
# Negative marker competitions fix Stromal over-prediction:
# - Stromal loses to Endothelial when CD31 is high
# - Stromal loses to Macrophages when CD68 is high
# - Stromal loses to T cells when CD3E is high
# - Stromal loses to Myofibroblasts when alphaSMA is high
# - Stromal loses to Epithelial when PanCK is high
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
INPUT_DIR="${PSEUDOVISIUM_DIR}/data_protein_gt"    # Protein-based GT
OUTPUT_DIR="${CITEGEIST_DIR}/output_protein_negmarkers"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

# Create directories if needed
mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Job started for region $SLURM_ARRAY_TASK_ID" | tee -a "${SLURM_LOG_DIR}/xenium_protein_negmarkers_runtime.log"

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
# Run CITEgeist with Protein GT + Negative Markers
# ============================================================================

REGION_ID=$SLURM_ARRAY_TASK_ID
echo -e "${YELLOW}Processing region ${REGION_ID} with protein GT + negative markers${RESET}"
echo -e "${BLUE}Input directory: ${INPUT_DIR}${RESET}"
echo -e "${BLUE}Output directory: ${OUTPUT_DIR}${RESET}"

# Change to repo root
cd "${REPO_ROOT}"

# Run CITEgeist with achievable profiles and negative marker redistribution
echo -e "${BLUE}Running CITEgeist on region ${REGION_ID}...${RESET}"
python "${CITEGEIST_DIR}/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 250.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --lambda-laplacian 0.0 \
    --no-unknown \
    --use-achievable-profiles \
    --use-negative-markers \
    --transfer-fraction 0.6

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

echo -e "[${YELLOW}${END_TIME}${RESET}] Job completed for region ${REGION_ID}" | tee -a "${SLURM_LOG_DIR}/xenium_protein_negmarkers_runtime.log"
echo -e "[${GREEN}PROTEIN_NEGMARKERS_RUNTIME${RESET}]: Region ${REGION_ID} took ${RUNTIME} seconds (${RUNTIME_MINUTES} minutes)" | tee -a "${SLURM_LOG_DIR}/xenium_protein_negmarkers_runtime.log"

echo -e "${GREEN}======================================${RESET}"
echo -e "${GREEN}Region ${REGION_ID} protein GT + negmarkers complete${RESET}"
echo -e "${GREEN}======================================${RESET}"
