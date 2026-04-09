#!/bin/bash
#SBATCH --job-name=mdk_expand
#SBATCH --partition=htc
#SBATCH --array=0-11
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/mdk_expand_%A_%a.out
#SBATCH --error=logs/mdk_expand_%A_%a.err

set -eo pipefail

echo "=== SACE Refilter + MDK Analysis ==="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
set -u

export HDF5_USE_FILE_LOCKING=FALSE

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
mkdir -p logs

python CITEgeist/examples/run_sace_refilter_mdk.py ${SLURM_ARRAY_TASK_ID}

echo "=== Done: $(date) ==="
