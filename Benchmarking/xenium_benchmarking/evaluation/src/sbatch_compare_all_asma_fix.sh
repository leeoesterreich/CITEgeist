#!/bin/bash
#SBATCH --job-name=compare_all_asma
#SBATCH --output=compare_all_asma_%j.out
#SBATCH --error=compare_all_asma_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Re-evaluate all methods against αSMA-only fibroblast GT

echo "=== Comparing All Methods Against αSMA-Only Fibroblast GT ==="
echo "Start time: $(date)"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py

echo "=== Done ==="
echo "End time: $(date)"
