#!/bin/bash
#SBATCH --job-name=discovery_cmp
#SBATCH --output=tests/slurm_log/discovery_cmp_%j.out
#SBATCH --error=tests/slurm_log/discovery_cmp_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python tests/test_discovery_comparison.py
