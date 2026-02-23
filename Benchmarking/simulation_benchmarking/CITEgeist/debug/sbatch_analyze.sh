#!/bin/bash
#SBATCH --job-name=analyze_mixed_highseg
#SBATCH --output=analyze_mixed_highseg_%j.out
#SBATCH --error=analyze_mixed_highseg_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Load conda
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CITEgeist/debug

python analyze_mixed_vs_highseg.py
