#!/bin/bash
#SBATCH --job-name=create_rna_gt
#SBATCH --output=create_rna_gt_%j.out
#SBATCH --error=create_rna_gt_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
python Benchmarking/xenium_pseudovisium/src/create_rna_gt.py
