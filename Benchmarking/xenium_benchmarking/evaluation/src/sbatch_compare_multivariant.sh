#!/bin/bash
#SBATCH --job-name=compare_multivariant
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/compare_multivariant_%j.out
#SBATCH --error=logs/compare_multivariant_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run proportion comparison
echo "=== Running Proportion Comparison ==="
python Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py

# Run GEX comparison
echo "=== Running GEX Comparison ==="
python Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py

echo "=== Complete ==="
