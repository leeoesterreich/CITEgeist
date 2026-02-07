#!/bin/bash
#SBATCH --job-name=protein_gt_gex
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/slurm/slurm_log/protein_gt_gex_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/slurm/slurm_log/protein_gt_gex_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

echo "=============================================="
echo "Generate GEX Ground Truth (Protein-Gated)"
echo "=============================================="
echo "Start time: $(date)"

CITEGEIST_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"
export PATH="${CITEGEIST_ENV}/bin:$PATH"
source activate "${CITEGEIST_ENV}"

SRC_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/src"

python "${SRC_DIR}/generate_protein_gt_gex.py"

echo ""
echo "Verifying output..."
ls -la /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth_gex/
for i in 0 1 2 3 4; do
    echo "Region ${i}:"
    ls /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth_gex/Xenium_region_${i}/
done

echo ""
echo "End time: $(date)"
