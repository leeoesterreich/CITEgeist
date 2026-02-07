#!/bin/bash
#SBATCH --job-name=gating_validation
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_log/gating_validation_%A_%a.out
#SBATCH --error=slurm_log/gating_validation_%A_%a.err
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate CITEgeist environment!" >&2
    exit 1
fi

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd $REPO_ROOT

# Create slurm_log directory if needed
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log

python Benchmarking/xenium_benchmarking/CITEgeist/src/validate_gating_classification.py \
    --region $SLURM_ARRAY_TASK_ID \
    --max-cells 10000
