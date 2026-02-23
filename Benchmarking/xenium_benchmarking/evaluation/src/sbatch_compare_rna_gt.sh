#!/bin/bash
#SBATCH --job-name=compare_rna_gt
#SBATCH --output=compare_rna_gt_%j.out
#SBATCH --error=compare_rna_gt_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Compare all methods against RNA GT with T cells combined
python Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py \
    --gt-dir Benchmarking/xenium_pseudovisium/data_rna_gt \
    --combine-t-cells \
    --output-suffix "_rna_gt"
