#!/bin/bash
#SBATCH --job-name=pass1_diag
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_log/pass1_diag_%j.out
#SBATCH --error=slurm_log/pass1_diag_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate CITEgeist environment!" >&2
    exit 1
fi

module load gurobi/12.0.3

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd $REPO_ROOT

python Benchmarking/xenium_benchmarking/CITEgeist/src/diagnose_pass1.py
