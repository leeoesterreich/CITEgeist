#!/bin/bash
#SBATCH --job-name=gen_granular_gt
#SBATCH --output=slurm_log/gen_granular_gt_%j.out
#SBATCH --error=slurm_log/gen_granular_gt_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc

# ============================================================================
# Generate Granular 10-Cell-Type Ground Truth Dataset
# ============================================================================
# Creates pseudo-Visium spots with UNSIMPLIFIED RNA clustering (10 cell types)
# to provide maximum granularity for benchmarking CITEgeist.
#
# Cell types: CD8+ T cells, Macrophages, Mixed Immune, Epithelial,
#             Myofibroblasts, Stromal, Endothelial, B cells,
#             Proliferating T, Vascular Stromal
# ============================================================================

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Set up environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

echo "============================================================"
echo "Generating Granular 10-Cell-Type Ground Truth Dataset"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Time: $(date)"
echo "Python: $PYTHON"
echo ""

$PYTHON Benchmarking/xenium_pseudovisium/src/generate_granular_dataset.py

echo ""
echo "Completed at: $(date)"
