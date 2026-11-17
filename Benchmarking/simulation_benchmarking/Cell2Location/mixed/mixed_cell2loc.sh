#!/bin/bash
#SBATCH --job-name=mx_cell2location  # Job name
#SBATCH --output=logs/%x_%a.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=logs/%x_%a.err    # Standard error (%x = job name, %j = job ID)
#SBATCH --ntasks=1                # Number of tasks (1 per job)
#SBATCH --cpus-per-task=2         # Number of CPU cores per task
#SBATCH --mem=32G                 # Memory per node
#SBATCH --time=05:00:00           # Time limit (hh:mm:ss)
# gpu-race-submit.sh sets --cluster and --partition
#SBATCH --gres=gpu:1              # Request one GPU
#SBATCH --array=0-4               # Array of jobs (one per replicate)

# Color codes
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
RED="\033[1;31m"
RESET="\033[0m"

# Record the start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "[${YELLOW}${START_TIME}${RESET}] Job started for replicate $SLURM_ARRAY_TASK_ID" | tee -a logs/mx_cell2location_runtime.log

# Load Miniconda environment
eval "$(conda shell.bash hook 2>/dev/null)"
echo -e "${BLUE}Loading cell2location environment...${RESET}"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/cell2location_env
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to activate environment!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}Miniconda environment loaded successfully.${RESET}"

# Change to project directory
echo -e "${BLUE}Changing to project directory...${RESET}"
cd ${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/Cell2Location/mixed
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
echo -e "${BLUE}Running Python script: mixed_cell2loc.py${RESET}"
python mixed_cell2loc.py "${REPLICATE_NAME}"

# Check for errors
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Python script failed on replicate ${REPLICATE_NAME}!${RESET}" >&2
    exit 1
fi
echo -e "${GREEN}Python script completed successfully for replicate: ${REPLICATE_NAME}${RESET}"

# Record the end time
END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')

# Calculate total runtime
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo -e "[${YELLOW}${END_TIME}${RESET}] Job completed for replicate $SLURM_ARRAY_TASK_ID" | tee -a logs/mx_cell2location_runtime.log
echo -e "[${GREEN}MX_CELL2LOCATION_TOTAL_RUNTIME${RESET}]: Replicate $SLURM_ARRAY_TASK_ID took $RUNTIME seconds ($RUNTIME_MINUTES minutes)" | tee -a logs/mx_cell2location_runtime.log

