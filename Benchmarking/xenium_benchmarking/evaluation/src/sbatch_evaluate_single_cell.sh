#!/bin/bash
#SBATCH --job-name=eval_sc_res
#SBATCH --output=logs/eval_sc_res_%j.out
#SBATCH --error=logs/eval_sc_res_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Evaluate single-cell resolution benchmark results

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/evaluation/src/evaluate_single_cell_resolution.py \
    --sc-dir Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_resolution \
    --gt-dir Benchmarking/xenium_pseudovisium/data_protein_gt \
    --regions 0,1,2,3,4 \
    --output Benchmarking/xenium_benchmarking/evaluation/results/single_cell_resolution_results.json
