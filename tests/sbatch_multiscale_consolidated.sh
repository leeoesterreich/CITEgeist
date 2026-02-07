#!/bin/bash
#SBATCH --job-name=multiscale_bench
#SBATCH --output=slurm_log/multiscale_%A_%a.out
#SBATCH --error=slurm_log/multiscale_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=HTC
#SBATCH --array=0-9  # 2 test sets × 5 replicates

# Consolidated benchmark for multi-scale neighborhoods
# Tests profile discovery with larger scales [6, 12, 24, 48, 64]

# Activate conda environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory: $(pwd)"

# Test sets and replicates
TEST_SETS=("high_seg" "mixed")
N_REPLICATES=5

# Calculate indices
test_set_index=$((SLURM_ARRAY_TASK_ID / N_REPLICATES))
replicate=$((SLURM_ARRAY_TASK_ID % N_REPLICATES))  # 0-indexed

# Get parameter values
TEST_SET=${TEST_SETS[$test_set_index]}

echo "======================================"
echo "Testing MULTI-SCALE NEIGHBORHOODS on ${TEST_SET} replicate ${replicate}"
echo "======================================"

# Create output directory
mkdir -p slurm_log

# Run Python benchmark
python -c "
import sys
import os
import numpy as np
import scanpy as sc

sys.path.append('CITEgeist')
from model.marker_interest import identify_interesting_markers
from model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,
)

# Load data
test_set = '${TEST_SET}'
rep = ${replicate}

data_dir = f'replicates/{test_set}/h5ad_objects/'
adata = sc.read_h5ad(f'{data_dir}/Wu_rep_{rep}_CITE.h5ad')
print(f'Loaded {test_set} rep {rep}: {adata.shape}')

# Get expression matrix
X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
coords = adata.obsm['spatial']
marker_names = list(adata.var_names)

# Module 1: Identify interesting markers
interest_result = identify_interesting_markers(
    X, coords, marker_names,
    kurtosis_threshold=2.0,
    morans_threshold=0.05,
    gmm_snr_threshold=1.0,
    verbose=True,
)
print(f'Module 1: {len(interest_result.interesting_markers)} interesting markers')

# Module 2a: Analyze colocalization with multi-scale (default: [6, 12, 24, 48, 64])
coloc_result = analyze_marker_colocalization(
    X, coords, marker_names,
    markers_to_analyze=interest_result.interesting_markers,
    neighbor_k=6,
    smooth_k=6,
    n_permutations=199,
    verbose=True,
)

# Check multi-scale results
pairs_with_multiscale = sum(1 for p in coloc_result.pairs if p.bivariate_morans_per_scale is not None)
print(f'Pairs with multi-scale: {pairs_with_multiscale}/{len(coloc_result.pairs)}')

if pairs_with_multiscale > 0:
    # Distribution of best scales
    scale_counts = {}
    for p in coloc_result.pairs:
        if p.bivariate_morans_best_scale is not None:
            scale = p.bivariate_morans_best_scale
            scale_counts[scale] = scale_counts.get(scale, 0) + 1
    print(f'Best scale distribution: {scale_counts}')

# Module 2b: Discover profiles
discovery_result = discover_profiles(
    coloc_result,
    fdr_alpha=0.05,
    top_k=5,
    verbose=True,
)
print(f'Module 2b: {len(discovery_result.profiles)} profiles')

# Module 2c: Select profiles
selection_result = select_profiles(
    X, coords, marker_names,
    profiles=discovery_result.profiles,
    interesting_markers=interest_result.interesting_markers,
    colocalization_result=coloc_result,
    min_spatial_explained=0.90,
    min_protein_explained=0.90,
    max_profiles=15,
    min_profiles=5,
    verbose=True,
)
print(f'Module 2c: {selection_result.optimal_n} profiles selected')

# Ground truth cell types (9 types, 2 markers each)
GROUND_TRUTH_TYPES = [
    'B-cells', 'CAFs', 'Cancer_Epithelial', 'Endothelial', 'Myeloid',
    'Normal_Epithelial', 'PVL', 'Plasmablasts', 'T-cells'
]

# Count cell types found
found_types = set()
for profile in selection_result.selected_profiles:
    for marker in profile:
        # Extract cell type from marker name (e.g., 'B-cells_Protein_1' -> 'B-cells')
        for ct in GROUND_TRUTH_TYPES:
            if ct.replace('_', ' ') in marker or ct in marker:
                found_types.add(ct)

print()
print('=== RESULTS ===')
print(f'Cell types found: {len(found_types)}/{len(GROUND_TRUTH_TYPES)}')
print(f'Found: {sorted(found_types)}')
print(f'Missing: {sorted(set(GROUND_TRUTH_TYPES) - found_types)}')
print('Done')
"

echo ""
echo "Job completed at $(date)"
