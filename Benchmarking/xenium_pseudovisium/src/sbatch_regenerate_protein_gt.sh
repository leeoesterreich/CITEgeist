#!/bin/bash
#SBATCH --job-name=regen_protein_gt
#SBATCH --output=regen_protein_gt_%j.out
#SBATCH --error=regen_protein_gt_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Regenerate protein-based ground truth with αSMA-only fibroblast definition
# (Vimentin gate removed due to lack of specificity in RCC)

echo "=== Regenerating Protein GT (αSMA-only fibroblasts) ==="
echo "Start time: $(date)"

# Activate CITEgeist conda environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Backup old ground truth
BACKUP_DIR="Benchmarking/xenium_pseudovisium/data_protein_gt_backup_$(date +%Y%m%d_%H%M%S)"
echo "Backing up old ground truth to: ${BACKUP_DIR}"
cp -r Benchmarking/xenium_pseudovisium/data_protein_gt "${BACKUP_DIR}"

# Regenerate ground truth
python Benchmarking/xenium_pseudovisium/src/create_protein_gt.py

echo "=== Done ==="
echo "End time: $(date)"
