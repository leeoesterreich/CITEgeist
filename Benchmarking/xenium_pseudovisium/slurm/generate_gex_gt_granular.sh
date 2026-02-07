#!/bin/bash
#SBATCH --job-name=gex_gt_granular
#SBATCH --output=slurm_log/gex_gt_granular_%j.out
#SBATCH --error=slurm_log/gex_gt_granular_%j.err
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=HTC

# Generate GEX ground truth for granular 10-cell-type classification
# Uses RNA-based clustering WITHOUT simplification to preserve all 10 cell types

# Activate conda environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium

# Create slurm_log directory if needed
mkdir -p slurm_log

echo "Generating GEX ground truth for granular 10 cell types..."
echo "Output directory: data_granular_gt/ground_truth_gex"

# Run the generate script with granular (unsimplified) cell types
# Need to modify the script call to use simplify=False
python -c "
import sys
sys.path.insert(0, 'src')

import os
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from load_xenium import load_xenium_data, split_gex_protein
from create_pseudo_spots import create_pseudo_visium_spots
from rna_cell_types import load_rna_clusters, XENIUM_KIDNEY_RNA_CLUSTER_MAP, GRANULAR_CELLTYPE_MAP
from split_regions import split_tissue_regions

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants
DATA_DIR = '/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma'
OUTPUT_BASE = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt')

logger.info('=' * 60)
logger.info('Generating GRANULAR (10 cell type) GEX Ground Truth')
logger.info('=' * 60)

# Load data
logger.info('Loading Xenium data...')
adata = load_xenium_data(DATA_DIR)
adata_gex_cells, adata_protein_cells = split_gex_protein(adata)
logger.info(f'  Cells: {len(adata)}, Genes: {adata_gex_cells.shape[1]}')

# Create pseudo-spots
logger.info('Creating pseudo-Visium spots...')
adata_spots, cell_to_spot = create_pseudo_visium_spots(
    adata,
    spot_diameter=55.0,
    center_spacing=100.0,
    min_cells=3,
)
spot_names = list(adata_spots.obs_names)
logger.info(f'  Created {len(spot_names)} spots')

# Load RNA clusters (no simplification - keep all 10 types)
logger.info('Loading RNA clusters with GRANULAR 10 cell types (no simplification)...')
clusters_df = load_rna_clusters(DATA_DIR, n_clusters=10)

# Apply granular mapping (keeps all 10 types)
cell_types = clusters_df['Cluster'].map(XENIUM_KIDNEY_RNA_CLUSTER_MAP)
# Note: GRANULAR_CELLTYPE_MAP just returns the same names, so we use the raw cluster map
cell_types = cell_types.fillna('Unknown')

# Align with adata
common_cells = cell_types.index.intersection(adata.obs_names)
cell_types = cell_types.loc[common_cells]
cell_types = cell_types.reindex(adata.obs_names).fillna('Unknown')

logger.info('Cell type distribution:')
for ct, count in cell_types.value_counts().items():
    logger.info(f'  {ct}: {count} ({100*count/len(cell_types):.1f}%)')

# Calculate GEX ground truth
logger.info('Calculating GEX ground truth for all spots...')

unique_types = sorted(cell_types.unique())
n_spots = len(spot_names)
gene_names = list(adata_gex_cells.var_names)

logger.info(f'  Cell types: {len(unique_types)}, Spots: {n_spots}, Genes: {len(gene_names)}')

# Extract aligned data
cell_spot_ids = cell_to_spot.reindex(adata_gex_cells.obs_names)['spot_id'].values
cell_types_aligned = cell_types.reindex(adata_gex_cells.obs_names).values

gex_by_celltype = {}

for ct in unique_types:
    if ct == 'Unknown':
        continue
    logger.info(f'  Processing {ct}...')

    # Initialize (genes x spots) matrix
    ct_matrix = pd.DataFrame(0.0, index=gene_names, columns=spot_names)

    ct_mask = cell_types_aligned == ct

    for spot_idx, spot_name in enumerate(spot_names):
        spot_mask = cell_spot_ids == spot_name
        combined_mask = ct_mask & spot_mask

        if combined_mask.sum() > 0:
            if hasattr(adata_gex_cells.X, 'toarray'):
                cell_gex = adata_gex_cells.X[combined_mask, :].toarray()
            else:
                cell_gex = adata_gex_cells.X[combined_mask, :]

            summed_gex = cell_gex.sum(axis=0).flatten()
            ct_matrix.iloc[:, spot_idx] = summed_gex

    n_active_spots = (ct_matrix.sum(axis=0) > 0).sum()
    gex_by_celltype[ct] = ct_matrix
    logger.info(f'    {ct}: {n_active_spots} spots with expression')

# Split into 5 regions
logger.info('Splitting into 5 regions...')
region_masks = split_tissue_regions(adata_spots, n_regions=5)

# Save GEX ground truth per region
logger.info('Saving GEX ground truth per region...')
for region_id, mask in enumerate(region_masks):
    region_gt_dir = OUTPUT_BASE / 'ground_truth_gex' / f'Xenium_region_{region_id}'
    region_gt_dir.mkdir(parents=True, exist_ok=True)

    region_spot_names = [spot_names[i] for i in range(len(spot_names)) if mask[i]]

    for ct, df in gex_by_celltype.items():
        region_df = df[region_spot_names]

        # Safe filename
        ct_safe = ct.replace(' ', '_').replace('+', 'pos')
        filename = f'{ct_safe}_GT.csv'
        filepath = region_gt_dir / filename
        region_df.to_csv(filepath)

    logger.info(f'Region {region_id} ({mask.sum()} spots): Saved {len(gex_by_celltype)} cell type files')

logger.info('=' * 60)
logger.info('GRANULAR GEX ground truth generation complete!')
logger.info(f'Output: {OUTPUT_BASE / \"ground_truth_gex\"}')
logger.info('=' * 60)
"

echo "GEX ground truth generation completed."
