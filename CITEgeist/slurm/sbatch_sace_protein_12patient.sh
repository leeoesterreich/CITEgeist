#!/bin/bash
#SBATCH --job-name=sace_protein_12pt
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/slurm/slurm_log/%x_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/slurm/slurm_log/%x_%A_%a.err
#SBATCH --cluster=htc
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# SACE per-cell protein annotation (v2) for all 12 HCC22-088 patient samples.
# Reads proportions from module3_cuopt_qp (Epithelial profile), runs SACE GEX
# + SACE protein (M3.5), and overwrites single_cell.h5ad with updated outputs.
# Run AFTER sbatch_qp_12patient.sh completes.  No GPU required.

set -eo pipefail

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
SRC_SCRIPT="${REPO_ROOT}/CITEgeist/examples/run_sace_protein_12patient.py"
LOG_DIR="${REPO_ROOT}/CITEgeist/slurm/slurm_log"

mkdir -p "${LOG_DIR}"

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}/CITEgeist"

NUMBA_CACHE_DIR=/tmp/numba_cache \
python "${SRC_SCRIPT}" \
  --array-index "${SLURM_ARRAY_TASK_ID}"

echo "SACE protein complete for array index ${SLURM_ARRAY_TASK_ID}"
