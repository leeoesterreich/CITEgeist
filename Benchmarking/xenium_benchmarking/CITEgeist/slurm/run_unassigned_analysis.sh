#!/bin/bash
#SBATCH --job-name=unassigned_analysis
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_log/unassigned_analysis_%j.out
#SBATCH --error=slurm_log/unassigned_analysis_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

echo Job started

ENV=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
export PATH=$ENV/bin:$PATH
source activate $ENV

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
python Benchmarking/xenium_benchmarking/CITEgeist/src/analyze_unassigned_doublets.py

echo Job finished
