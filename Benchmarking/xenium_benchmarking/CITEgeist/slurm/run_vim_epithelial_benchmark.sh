#!/bin/bash
#SBATCH --job-name=vim_epi_bench
#SBATCH --output=slurm_log/%x_%a.out
#SBATCH --error=slurm_log/%x_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Experiment: Vimentin added to Epithelial profile
# ============================================================================
# Tests hypothesis that fibroblasts absorb proportion from VIM+ epithelial
# cells because only fibroblasts claim Vimentin in the baseline profiles.
#
# Change: Epithelial ["PanCK"] → ["PanCK", "Vimentin"]
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
CITEGEIST_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
INPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_granular_gt"
OUTPUT_DIR="${CITEGEIST_DIR}/output/vim_epithelial"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[${START_TIME}] VIM-Epithelial experiment: region $SLURM_ARRAY_TASK_ID"

# Environment setup
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

module load gurobi/12.0.3

cd "${REPO_ROOT}"

python "${CITEGEIST_DIR}/src/run_benchmark_vim_epithelial.py" \
    --region-id ${SLURM_ARRAY_TASK_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 4.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --run-gex

if [ $? -ne 0 ]; then
    echo "ERROR: Failed on region ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[${END_TIME}] Completed region ${SLURM_ARRAY_TASK_ID}"
