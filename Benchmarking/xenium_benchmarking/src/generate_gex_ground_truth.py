"""
Generate gene expression ground truth for Xenium benchmarking.

This script creates per-cell-type gene expression matrices from single-cell data,
which serve as ground truth for GEX deconvolution benchmarking.
"""

import os
import sys
import logging
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from load_xenium import load_xenium_data, split_gex_protein
from create_pseudo_spots import create_pseudo_visium_spots
from define_cell_types import XENIUM_CELL_PROFILE_DICT, classify_cells_by_protein
from split_regions import split_tissue_regions

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def calculate_spot_gex_by_celltype(
    adata_gex_cells: sc.AnnData,
    cell_to_spot: np.ndarray,
    cell_types: pd.Series,
    spot_names: List[str],
) -> Dict[str, pd.DataFrame]:
    """
    Calculate per-cell-type gene expression sums for each spot.

    Args:
        adata_gex_cells: AnnData with gene expression at single-cell level
        cell_to_spot: Array mapping each cell to a spot index (-1 = unassigned)
        cell_types: Series with cell type labels for each cell
        spot_names: List of spot names

    Returns:
        Dict mapping cell type names to DataFrames (genes x spots)
    """
    unique_types = sorted(cell_types.unique())
    n_spots = len(spot_names)
    gene_names = list(adata_gex_cells.var_names)

    logger.info(f"Calculating GEX ground truth for {len(unique_types)} cell types, {n_spots} spots")

    gex_by_celltype = {}

    for ct in unique_types:
        logger.info(f"  Processing {ct}...")

        # Initialize (genes x spots) matrix for this cell type
        ct_matrix = pd.DataFrame(
            0.0,
            index=gene_names,
            columns=spot_names,
        )

        # Find cells of this type
        ct_mask = cell_types == ct

        for spot_idx in range(n_spots):
            # Find cells of this type in this spot
            spot_mask = cell_to_spot == spot_idx
            combined_mask = ct_mask.values & spot_mask

            if combined_mask.sum() > 0:
                # Sum gene expression across all cells of this type in this spot
                if hasattr(adata_gex_cells.X, 'toarray'):
                    cell_gex = adata_gex_cells.X[combined_mask, :].toarray()
                else:
                    cell_gex = adata_gex_cells.X[combined_mask, :]

                summed_gex = cell_gex.sum(axis=0).flatten()
                ct_matrix.iloc[:, spot_idx] = summed_gex

        n_active_spots = (ct_matrix.sum(axis=0) > 0).sum()
        gex_by_celltype[ct] = ct_matrix
        logger.info(f"    {ct}: {n_active_spots} spots with expression")

    return gex_by_celltype


def save_region_gex_ground_truth(
    gex_by_celltype: Dict[str, pd.DataFrame],
    region_mask: np.ndarray,
    spot_names: List[str],
    output_dir: str,
) -> None:
    """
    Save GEX ground truth for a single region.

    Args:
        gex_by_celltype: Dict mapping cell type to (genes x spots) DataFrames
        region_mask: Boolean mask for spots in this region
        spot_names: Full list of spot names
        output_dir: Directory to save files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get spots in this region
    region_spot_names = [spot_names[i] for i in range(len(spot_names)) if region_mask[i]]

    for ct, df in gex_by_celltype.items():
        # Subset to region spots
        region_df = df[region_spot_names]

        # Replace spaces/special chars for filenames
        ct_safe = ct.replace(" ", "_").replace("+", "pos")
        filename = f"{ct_safe}_GT.csv"
        filepath = output_dir / filename
        region_df.to_csv(filepath)

    logger.info(f"  Saved {len(gex_by_celltype)} cell type GEX files to {output_dir}")


def main():
    """Generate GEX ground truth for all Xenium regions."""
    DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"
    OUTPUT_BASE = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/data")

    # Load original data
    logger.info("Loading Xenium data...")
    adata = load_xenium_data(DATA_DIR)
    adata_gex_cells, adata_protein_cells = split_gex_protein(adata)

    logger.info(f"  Cells: {len(adata)}")
    logger.info(f"  Genes: {adata_gex_cells.shape[1]}")
    logger.info(f"  Proteins: {adata_protein_cells.shape[1]}")

    # Create pseudo-spots (to get cell_to_spot mapping)
    logger.info("\nCreating pseudo-Visium spots...")
    adata_spots, cell_to_spot = create_pseudo_visium_spots(
        adata,
        spot_diameter=55.0,
        center_spacing=100.0,
        min_cells=3,
    )
    spot_names = list(adata_spots.obs_names)
    logger.info(f"  Created {len(spot_names)} spots")

    # Classify cells
    logger.info("\nClassifying cells by protein expression...")
    cell_types = classify_cells_by_protein(
        adata_protein_cells,
        XENIUM_CELL_PROFILE_DICT,
        threshold_method="percentile",
        percentile=50.0,
    )

    # Rename Unassigned to Unknown to match CITEgeist naming
    cell_types = cell_types.replace("Unassigned", "Unknown")

    type_counts = cell_types.value_counts()
    logger.info("  Cell type distribution:")
    for ct, count in type_counts.items():
        logger.info(f"    {ct}: {count} ({100*count/len(cell_types):.1f}%)")

    # Calculate GEX ground truth for all spots
    logger.info("\nCalculating GEX ground truth for all spots...")
    gex_by_celltype = calculate_spot_gex_by_celltype(
        adata_gex_cells,
        cell_to_spot,
        cell_types,
        spot_names,
    )

    # Split into regions
    logger.info("\nSplitting into 5 regions...")
    region_masks = split_tissue_regions(adata_spots, n_regions=5)

    # Save GEX ground truth for each region
    logger.info("\nSaving GEX ground truth per region...")
    for region_id, mask in enumerate(region_masks):
        region_gt_dir = OUTPUT_BASE / "ground_truth_gex" / f"Xenium_region_{region_id}"
        logger.info(f"Region {region_id} ({mask.sum()} spots):")
        save_region_gex_ground_truth(
            gex_by_celltype,
            mask,
            spot_names,
            str(region_gt_dir),
        )

    logger.info("\nGEX ground truth generation complete!")


if __name__ == "__main__":
    main()
