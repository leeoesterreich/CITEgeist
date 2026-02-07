#!/bin/bash
#SBATCH --job-name=panel_check
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=slurm_log/panel_check_%j.out
#SBATCH --error=slurm_log/panel_check_%j.err
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

python Benchmarking/xenium_benchmarking/CITEgeist/src/check_panel_discrimination.py
