#!/bin/bash
#SBATCH --job-name=cell_resolution_validation
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_log/cell_resolution_%A_%a.out
#SBATCH --error=slurm_log/cell_resolution_%A_%a.err
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

module load python/3.10
conda activate CITEgeist_env

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd $REPO_ROOT

python Benchmarking/xenium_benchmarking/CITEgeist/src/validate_cell_resolution.py \
    --region $SLURM_ARRAY_TASK_ID \
    --max-cells 50000
