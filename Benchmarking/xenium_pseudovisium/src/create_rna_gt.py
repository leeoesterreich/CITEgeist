#!/usr/bin/env python
"""
Create RNA-based ground truth using Xenium's gene expression k-means clustering.

This script uses RNA clustering (not protein gating) to assign cell types,
then aggregates to pseudo-Visium spot-level proportions.

Key differences from protein GT:
- Uses RNA k-means clusters (analysis.tar.gz) for cell type assignment
- ALL cells get assigned (no Unknown category)
- Proportions sum to 1.0 per spot
- 6 cell types (T cells combined, RNA can't distinguish CD4/CD8)

Reference:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
"""

import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from rna_cell_types import load_rna_cell_types, SIMPLIFIED_CELLTYPE_MAP

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 6 cell types (T cells combined since RNA can't distinguish CD4/CD8)
RNA_CELL_TYPE_ORDER = [
    "B cells",
    "T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def calculate_spot_proportions(
    cell_types: pd.Series,
    cell_to_spot: pd.DataFrame,
    spot_coords: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate cell type proportions per spot from RNA-based cell types.

    Unlike protein GT, RNA assigns ALL cells so proportions sum to 1.0.
    """
    spot_ids = spot_coords.index.tolist()
    proportions = pd.DataFrame(0.0, index=spot_ids, columns=RNA_CELL_TYPE_ORDER)
    n_cells = pd.Series(0, index=spot_ids)

    # Align cell types with cell_to_spot mapping
    common_cells = cell_types.index.intersection(cell_to_spot.index)
    cell_types_aligned = cell_types.loc[common_cells]
    cell_to_spot_aligned = cell_to_spot.loc[common_cells]

    logger.info(f"Aligned {len(common_cells)} cells with spot mapping")

    # Calculate proportions per spot
    for spot_id in spot_ids:
        cells_in_spot = cell_to_spot_aligned[cell_to_spot_aligned['spot_id'] == spot_id].index

        if len(cells_in_spot) == 0:
            continue

        spot_cell_types = cell_types_aligned.loc[cells_in_spot]
        n_cells[spot_id] = len(cells_in_spot)

        # Count each type - ALL cells are assigned so total = n_cells
        type_counts = spot_cell_types.value_counts()
        total = len(spot_cell_types)

        for ct in RNA_CELL_TYPE_ORDER:
            if ct in type_counts.index:
                proportions.loc[spot_id, ct] = type_counts[ct] / total

    # Add metadata
    proportions['n_cells'] = n_cells
    proportions['spot_x'] = spot_coords['x']
    proportions['spot_y'] = spot_coords['y']

    return proportions


def main():
    """Create RNA-based ground truth for all regions."""

    # Paths
    xenium_dir = Path('/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma')
    pseudovisium_dir = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium')

    # Use existing cell_to_spot mapping from protein GT
    protein_gt_dir = pseudovisium_dir / 'data_protein_gt'

    # Output to new directory
    output_dir = pseudovisium_dir / 'data_rna_gt'
    output_dir.mkdir(exist_ok=True)
    (output_dir / 'ground_truth').mkdir(exist_ok=True)

    logger.info("=" * 60)
    logger.info("Creating RNA-Based Ground Truth")
    logger.info("=" * 60)

    # Load RNA-based cell types using existing module
    logger.info("Loading RNA k-means clusters...")
    cell_types = load_rna_cell_types(
        xenium_dir,
        n_clusters=10,
        simplify=True,  # Combines T cells, etc.
    )
    logger.info(f"Loaded {len(cell_types)} cell type assignments")

    # Show distribution
    logger.info("\nRNA-based cell type distribution:")
    for ct, count in cell_types.value_counts().items():
        pct = 100 * count / len(cell_types)
        logger.info(f"  {ct}: {count:,} ({pct:.1f}%)")

    # Verify no Unknown
    n_unknown = (cell_types == "Unknown").sum()
    if n_unknown > 0:
        logger.warning(f"Found {n_unknown} Unknown cells - these will be excluded")
        cell_types = cell_types[cell_types != "Unknown"]

    # Save cell-level classifications
    cell_types.to_csv(output_dir / 'cell_type_assignments.csv', header=['cell_type'])

    # Load cell-to-spot mapping (reuse from protein GT)
    logger.info("\nLoading cell-to-spot mapping...")
    cell_to_spot = pd.read_csv(
        protein_gt_dir / 'cell_to_spot_mapping.csv',
        index_col=0
    )
    cell_to_spot.index = cell_to_spot.index.astype(str)
    logger.info(f"Loaded mapping for {len(cell_to_spot)} cells")

    # Copy cell_to_spot mapping to RNA GT dir for completeness
    cell_to_spot.to_csv(output_dir / 'cell_to_spot_mapping.csv')

    # Process each region
    logger.info("\nCalculating spot-level proportions per region...")

    all_stats = []

    for region_id in range(5):
        logger.info(f"\n--- Region {region_id} ---")

        # Load region's spot coordinates from h5ad (GEX has spatial coords)
        h5ad_path = protein_gt_dir / 'h5ad_objects' / f'Xenium_region_{region_id}_GEX.h5ad'
        if not h5ad_path.exists():
            logger.warning(f"h5ad not found: {h5ad_path}, skipping region")
            continue

        import scanpy as sc
        adata = sc.read_h5ad(h5ad_path)

        # Get spot coordinates
        spot_coords = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y']
        )

        # Filter cell_to_spot to this region's spots
        region_spots = set(spot_coords.index)
        region_cell_to_spot = cell_to_spot[cell_to_spot['spot_id'].isin(region_spots)]

        logger.info(f"Region {region_id}: {len(region_spots)} spots, {len(region_cell_to_spot)} cells")

        # Calculate proportions
        proportions = calculate_spot_proportions(
            cell_types,
            region_cell_to_spot,
            spot_coords,
        )

        # Filter to spots with cells
        proportions = proportions[proportions['n_cells'] > 0]

        # Verify proportions sum to ~1.0
        prop_sums = proportions[RNA_CELL_TYPE_ORDER].sum(axis=1)
        logger.info(f"Proportion sums: mean={prop_sums.mean():.4f}, min={prop_sums.min():.4f}, max={prop_sums.max():.4f}")

        # Save
        output_path = output_dir / 'ground_truth' / f'Xenium_region_{region_id}_prop.csv'
        proportions.to_csv(output_path)
        logger.info(f"Saved: {output_path}")

        # Stats
        stats = {
            'region_id': region_id,
            'n_spots': len(proportions),
            'n_cells': int(proportions['n_cells'].sum()),
            'mean_cells_per_spot': float(proportions['n_cells'].mean()),
            'prop_sum_mean': float(prop_sums.mean()),
        }
        for ct in RNA_CELL_TYPE_ORDER:
            stats[f'{ct}_mean'] = float(proportions[ct].mean())
        all_stats.append(stats)

    # Save summary
    summary = {
        'description': 'RNA-based ground truth using k-means clustering',
        'cell_types': RNA_CELL_TYPE_ORDER,
        'n_cell_types': len(RNA_CELL_TYPE_ORDER),
        'total_cells': len(cell_types),
        'unknown_rate': 0.0,  # RNA assigns all cells
        'regions': all_stats,
    }

    with open(output_dir / 'dataset_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info("\n" + "=" * 60)
    logger.info("RNA Ground Truth Creation Complete")
    logger.info("=" * 60)
    logger.info(f"Output: {output_dir}")
    logger.info(f"Cell types: {RNA_CELL_TYPE_ORDER}")
    logger.info(f"Total cells: {len(cell_types):,}")
    logger.info(f"Unknown rate: 0% (RNA assigns all cells)")


if __name__ == '__main__':
    main()
