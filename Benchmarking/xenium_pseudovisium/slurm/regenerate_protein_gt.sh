#!/bin/bash
#SBATCH --job-name=regen_protein_gt
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/slurm/slurm_log/regen_protein_gt_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/slurm/slurm_log/regen_protein_gt_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"

echo "Re-generating protein GT with achievable-7 types..."
echo "Start time: $(date)"

python "${REPO}/Benchmarking/xenium_pseudovisium/src/create_protein_gt.py"

echo ""
echo "Verifying output..."
python -c "
import pandas as pd
gt = pd.read_csv('${REPO}/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth/Xenium_region_0_prop.csv', index_col=0)
print('Columns:', list(gt.columns))
print('Shape:', gt.shape)
ct_cols = [c for c in gt.columns if c not in ['n_cells','spot_x','spot_y']]
print(f'Cell type columns ({len(ct_cols)}):', ct_cols)
for c in ct_cols:
    print(f'  {c}: mean={gt[c].mean():.4f}')
"

echo ""
echo "Done! End time: $(date)"
