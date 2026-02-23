#!/bin/bash
#SBATCH --job-name=unknown_rna
#SBATCH --output=unknown_rna_%j.out
#SBATCH --error=unknown_rna_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
python Benchmarking/xenium_pseudovisium/src/analyze_unknown_rna_distribution.py
