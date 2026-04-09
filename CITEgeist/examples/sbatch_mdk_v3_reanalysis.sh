#!/bin/bash
#SBATCH --job-name=mdk_v3
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/mdk_v3_%j.out
#SBATCH --error=logs/mdk_v3_%j.err

set -eo pipefail

echo "=== MDK v3 Single-Cell Reanalysis ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"

# Activate environment (nounset disabled — conda scripts have unbound vars)
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
set -u

# Disable HDF5 file locking for NFS
export HDF5_USE_FILE_LOCKING=FALSE

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
mkdir -p logs

python CITEgeist/examples/run_mdk_v3_reanalysis.py

echo "=== Done: $(date) ==="
