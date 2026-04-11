#!/bin/bash
#SBATCH --job-name=biv-morans-v3
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/logs/biv_morans_v3_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

export HDF5_USE_FILE_LOCKING=FALSE

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

echo "=== Bivariate Moran's I v3 recomputation ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

python CITEgeist/examples/run_bivariate_morans_v3.py

echo "=== Done: $(date) ==="
