#!/bin/bash
#SBATCH --job-name=qc_12patient
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/slurm/slurm_log/%x_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/slurm/slurm_log/%x_%A_%a.err
#SBATCH --cluster=htc
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Self-consistency QC for all 12 HCC22-088 patient samples.
# Each task: load per-cell h5ad (~800-6000 cells x ~9000 genes), run QC, save PDFs.

set -eo pipefail

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
SRC_SCRIPT="${REPO_ROOT}/CITEgeist/examples/run_qc_12patient.py"
LOG_DIR="${REPO_ROOT}/CITEgeist/slurm/slurm_log"

mkdir -p "${LOG_DIR}"

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

NUMBA_CACHE_DIR=/tmp/numba_cache \
MPLBACKEND=Agg \
python "${SRC_SCRIPT}" --array-idx "${SLURM_ARRAY_TASK_ID}"
