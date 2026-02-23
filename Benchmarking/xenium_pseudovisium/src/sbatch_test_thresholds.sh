#!/bin/bash
#SBATCH --job-name=test_thresh
#SBATCH --output=test_thresh_%j.out
#SBATCH --error=test_thresh_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
python Benchmarking/xenium_pseudovisium/src/test_lower_thresholds.py
