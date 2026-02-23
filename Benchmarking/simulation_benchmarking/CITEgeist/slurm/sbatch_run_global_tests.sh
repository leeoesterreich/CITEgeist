#!/bin/bash
#SBATCH --job-name=global_iqp_tests
#SBATCH --output=logs/global_iqp_tests_%j.out
#SBATCH --error=logs/global_iqp_tests_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python -m pytest CITEgeist/tests/test_discrete_global_solve.py -v
