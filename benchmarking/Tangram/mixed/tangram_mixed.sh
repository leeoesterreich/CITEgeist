#!/bin/bash
#SBATCH --job-name=tangram_mixed  # Job name
#SBATCH --output=slurm_log/%x_%a.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=slurm_log/%x_%a.err    # Standard error (%x = job name, %j = job ID)
#SBATCH --ntasks=1                # Number of tasks (1 per job)
#SBATCH --cpus-per-task=1         # Number of CPU cores per task
#SBATCH --mem=32G                 # Memory per node
#SBATCH --time=03:00:00           # Time limit (hh:mm:ss)
#SBATCH -M gpu                    # Partition/queue name
#SBATCH --gres=gpu:1              # Request one GPU
#SBATCH --array=0-4               # Array of jobs (one per replicate)

# Color codes
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
RED="\033[1;31m"
RESET="\033[0m"

# Load Miniconda environment
export PATH=/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/miniconda3/bin:$PATH
echo -e "${BLUE}Loading Miniconda environment...${RESET}"
source /bgfs/alee/LO_LAB/Personal/Brent_Schlegel/miniconda3/bin/activate /bgfs/alee/LO_LAB/Personal/Brent_Schlegel/miniconda3/envs/tangram-env
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to activate environment!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}Miniconda environment loaded successfully.${RESET}"

# Change to project directory
echo -e "${BLUE}Changing to project directory...${RESET}"
cd /ihome/acillo/bts76/Simulations/scCube_12k/Tangram/mixed
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to change directory!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}Directory changed successfully.${RESET}"

# Replicate names
echo -e "${YELLOW}Defining replicate files...${RESET}"
REPLICATES=(
    "Wu_rep_0_GEX.h5ad"
    "Wu_rep_1_GEX.h5ad"
    "Wu_rep_2_GEX.h5ad"
    "Wu_rep_3_GEX.h5ad"
    "Wu_rep_4_GEX.h5ad"
)

# Print current task info
REPLICATE_NAME="${REPLICATES[$SLURM_ARRAY_TASK_ID]}"
echo -e "${YELLOW}Starting job ${SLURM_ARRAY_TASK_ID} for replicate: ${REPLICATE_NAME}${RESET}"

# Run Python script
echo -e "${BLUE}Running Tangram.py script${RESET}"
python Tangram.py \
    --spatial_path "/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/replicates/mixed/h5ad_objects/${REPLICATE_NAME}" \
    --output_dir "/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/Tangram/mixed/tangram_map_${SLURM_ARRAY_TASK_ID}" \
    --wu_reference_path "/bgfs/alee/LO_LAB/General/Public_Data/BC-Datasets/Wu_2021/Wu_etal_2021_BRCA_scRNASeq" \
    --device "cuda:0"

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Tangram.py script failed on replicate ${REPLICATE_NAME}!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}Tangram.py script completed successfully for replicate: ${REPLICATE_NAME}${RESET}"
