#!/bin/bash
#SBATCH --job-name=regen_gex_gt
#SBATCH --output=regen_gex_gt_%j.out
#SBATCH --error=regen_gex_gt_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Regenerate GEX ground truth with αSMA-only fibroblast assignments

echo "=== Regenerating GEX Ground Truth ==="
echo "Start time: $(date)"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_pseudovisium/src/generate_protein_gt_gex.py

echo "=== Done ==="
echo "End time: $(date)"
