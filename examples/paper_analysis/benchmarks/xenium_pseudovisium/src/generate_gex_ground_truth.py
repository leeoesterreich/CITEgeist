"""
Generate gene expression ground truth for Xenium benchmarking.

This script creates per-cell-type gene expression matrices from single-cell data,
which serve as ground truth for GEX deconvolution benchmarking.

Supports two cell type classification methods:
1. RNA-based (recommended): Uses RNA clustering from Xenium analysis output
2. Protein-based: Uses protein marker thresholding

Reference for RNA-based approach:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0
"""

import argparse
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
from rna_cell_types import load_rna_cell_types
from split_regions import split_tissue_regions

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def calculate_spot_gex_by_celltype(
    adata_gex_cells: sc.AnnData,
    cell_to_spot_df: pd.DataFrame,
    cell_types: pd.Series,
    spot_names: List[str],
) -> Dict[str, pd.DataFrame]:
    """
    Calculate per-cell-type gene expression sums for each spot.

    Args:
        adata_gex_cells: AnnData with gene expression at single-cell level
        cell_to_spot_df: DataFrame mapping cells to spots (columns: spot_idx, spot_id)
        cell_types: Series with cell type labels for each cell
        spot_names: List of spot names

    Returns:
        Dict mapping cell type names to DataFrames (genes x spots)
    """
    unique_types = sorted(cell_types.unique())
    n_spots = len(spot_names)
    gene_names = list(adata_gex_cells.var_names)

    logger.info(f"Calculating GEX ground truth for {len(unique_types)} cell types, {n_spots} spots")

    # Extract spot_id column as array aligned with cells
    # Ensure alignment with adata_gex_cells
    cell_spot_ids = cell_to_spot_df.reindex(adata_gex_cells.obs_names)["spot_id"].values

    # Convert cell_types to array aligned with adata
    cell_types_aligned = cell_types.reindex(adata_gex_cells.obs_names).values

    gex_by_celltype = {}

    for ct in unique_types:
        logger.info(f"  Processing {ct}...")

        # Initialize (genes x spots) matrix for this cell type
        ct_matrix = pd.DataFrame(
            0.0,
            index=gene_names,
            columns=spot_names,
        )

        # Find cells of this type (as boolean array)
        ct_mask = cell_types_aligned == ct

        for spot_idx, spot_name in enumerate(spot_names):
            # Find cells of this type in this spot
            spot_mask = cell_spot_ids == spot_name
            combined_mask = ct_mask & spot_mask

            if combined_mask.sum() > 0:
                # Sum gene expression across all cells of this type in this spot
                if hasattr(adata_gex_cells.X, "toarray"):
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
    parser = argparse.ArgumentParser(description="Generate gene expression ground truth for Xenium benchmarking")
    parser.add_argument(
        "--use-rna-types",
        action="store_true",
        default=True,
        help="Use RNA-based cell types (recommended, default: True)",
    )
    parser.add_argument(
        "--use-protein-types",
        action="store_true",
        help="Use protein-based cell types (overrides --use-rna-types)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: auto-determined based on cell type method)",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="/path/to/Xenium_RCC",
        help="Xenium data directory",
    )
    args = parser.parse_args()

    # Determine cell type method
    use_rna = not args.use_protein_types  # RNA is default unless protein is specified

    # Set output directory based on method
    REPO_ROOT = Path("/path/to/CITEgeist_analysis")
    if args.output_dir:
        OUTPUT_BASE = Path(args.output_dir)
    elif use_rna:
        OUTPUT_BASE = REPO_ROOT / "benchmarks" / "xenium_pseudovisium" / "data_rna_gt"
    else:
        OUTPUT_BASE = REPO_ROOT / "benchmarks" / "xenium_pseudovisium" / "data"

    logger.info("=" * 60)
    logger.info("Gene Expression Ground Truth Generation")
    logger.info("=" * 60)
    logger.info(f"Cell type method: {'RNA-based' if use_rna else 'Protein-based'}")
    logger.info(f"Output directory: {OUTPUT_BASE}")
    logger.info("=" * 60)

    # Load original data
    logger.info("\nLoading Xenium data...")
    adata = load_xenium_data(args.data_dir)
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

    # Classify cells based on method
    if use_rna:
        logger.info("\nClassifying cells using RNA clustering (recommended)...")
        cell_types = load_rna_cell_types(
            args.data_dir,
            n_clusters=10,
            simplify=True,  # Use 6 simplified cell types
        )
        # Align cell types with adata index
        common_cells = cell_types.index.intersection(adata.obs_names)
        cell_types = cell_types.loc[common_cells]
        # Reindex to match adata order
        cell_types = cell_types.reindex(adata.obs_names).fillna("Unknown")
    else:
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

    logger.info("\n" + "=" * 60)
    logger.info("GEX ground truth generation complete!")
    logger.info(f"Output saved to: {OUTPUT_BASE / 'ground_truth_gex'}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
