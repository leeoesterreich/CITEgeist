#!/bin/bash
#SBATCH --job-name=discrete_slack_only
#SBATCH --output=logs/discrete_slack_only_%A_%a.out
#SBATCH --error=logs/discrete_slack_only_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Test ONLY constraint slack (±1 cell) with CLR preprocessing
# NO regularization (lambda_reg=0), NO marker weighting
# This isolates the effect of flexible nuclei constraints
python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id ${SLURM_ARRAY_TASK_ID} \
    --condition mixed \
    --mode dapi \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_slack_only/mixed/dapi \
    --use-clr-preprocessing \
    --constraint-slack 1 \
    --lambda-reg 0.0 \
    --global-solve \
    --global-time-limit 600 \
    --global-mip-gap 0.05 \
    --max-em-iterations 10
