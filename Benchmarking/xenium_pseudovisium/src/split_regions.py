"""
Split tissue into non-overlapping regions for pseudo-replicates.

This module provides functions to divide the pseudo-Visium spots
into spatial regions for creating independent benchmark replicates.

IMPORTANT: RNA-based clustering is now the default for ground truth.
---------------------------------------------------------------------
As of 2025, the recommended approach for cell type annotation in Xenium
is to use gene expression clustering rather than protein thresholding.
This avoids circular logic in benchmarking.

Reference:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0
"""

import logging
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


def split_tissue_regions(
    adata_spots: sc.AnnData,
    n_regions: int = 5,
    method: str = "horizontal_strips",
) -> List[np.ndarray]:
    """
    Split spots into non-overlapping spatial regions.

    Args:
        adata_spots: AnnData with pseudo-spot data (must have obsm['spatial'])
        n_regions: Number of regions to create
        method: Splitting method:
            - "horizontal_strips": divide by x-coordinate
            - "vertical_strips": divide by y-coordinate

    Returns:
        List of boolean mask arrays, one per region
    """
    spatial = adata_spots.obsm["spatial"]
    n_spots = len(adata_spots)

    if method == "horizontal_strips":
        # Divide by x-coordinate (horizontal strips along y-axis)
        coord = spatial[:, 0]  # x
    elif method == "vertical_strips":
        # Divide by y-coordinate (vertical strips along x-axis)
        coord = spatial[:, 1]  # y
    else:
        raise ValueError(f"Unknown method: {method}")

    # Calculate boundaries
    coord_min, coord_max = coord.min(), coord.max()
    boundaries = np.linspace(coord_min, coord_max, n_regions + 1)

    # Assign spots to regions
    region_masks = []
    for i in range(n_regions):
        if i == n_regions - 1:
            # Last region includes the boundary
            mask = (coord >= boundaries[i]) & (coord <= boundaries[i + 1])
        else:
            mask = (coord >= boundaries[i]) & (coord < boundaries[i + 1])
        region_masks.append(mask)

    # Log region sizes
    for i, mask in enumerate(region_masks):
        logger.info(f"Region {i}: {mask.sum()} spots ({100*mask.sum()/n_spots:.1f}%)")

    return region_masks


def save_region_datasets(
    adata_gex: sc.AnnData,
    adata_protein: sc.AnnData,
    proportions_df: pd.DataFrame,
    region_masks: List[np.ndarray],
    output_dir: str,
    prefix: str = "Xenium",
) -> List[dict]:
    """
    Save separate h5ad files and ground truth for each region.

    Args:
        adata_gex: AnnData with gene expression (spot-level)
        adata_protein: AnnData with protein expression (spot-level)
        proportions_df: DataFrame with ground truth proportions
        region_masks: List of boolean masks for each region
        output_dir: Base output directory
        prefix: Prefix for filenames

    Returns:
        List of dicts with paths to saved files per region
    """
    output_dir = Path(output_dir)
    h5ad_dir = output_dir / "h5ad_objects"
    gt_dir = output_dir / "ground_truth"

    h5ad_dir.mkdir(parents=True, exist_ok=True)
    gt_dir.mkdir(parents=True, exist_ok=True)

    saved_files = []

    for region_id, mask in enumerate(region_masks):
        region_spots = adata_gex.obs_names[mask]

        # Subset data
        adata_gex_region = adata_gex[region_spots, :].copy()
        adata_protein_region = adata_protein[region_spots, :].copy()
        props_region = proportions_df.loc[region_spots].copy()

        # File paths
        gex_path = h5ad_dir / f"{prefix}_region_{region_id}_GEX.h5ad"
        protein_path = h5ad_dir / f"{prefix}_region_{region_id}_CITE.h5ad"
        gt_path = gt_dir / f"{prefix}_region_{region_id}_prop.csv"

        # Save files
        adata_gex_region.write_h5ad(gex_path)
        adata_protein_region.write_h5ad(protein_path)
        props_region.to_csv(gt_path)

        logger.info(
            f"Region {region_id}: saved {len(region_spots)} spots "
            f"to {gex_path.name}, {protein_path.name}, {gt_path.name}"
        )

        saved_files.append({
            "region_id": region_id,
            "n_spots": len(region_spots),
            "gex_path": str(gex_path),
            "protein_path": str(protein_path),
            "gt_path": str(gt_path),
        })

    return saved_files


def prepare_citegeist_input(
    adata_gex: sc.AnnData,
    adata_protein: sc.AnnData,
) -> Tuple[sc.AnnData, sc.AnnData]:
    """
    Prepare AnnData objects for CITEgeist input.

    Ensures proper formatting:
    - feature_types column set correctly
    - var_names made unique if needed

    Args:
        adata_gex: Gene expression AnnData
        adata_protein: Protein expression AnnData

    Returns:
        Tuple of (adata_gex_prepared, adata_protein_prepared)
    """
    adata_gex = adata_gex.copy()
    adata_protein = adata_protein.copy()

    # Ensure feature_types are set correctly
    adata_gex.var["feature_types"] = "Gene Expression"
    adata_protein.var["feature_types"] = "Antibody Capture"

    # Make var_names unique if needed
    adata_gex.var_names_make_unique()
    adata_protein.var_names_make_unique()

    return adata_gex, adata_protein


def create_all_region_datasets(
    data_dir: str,
    output_dir: str,
    n_regions: int = 5,
    spot_diameter: float = 55.0,
    center_spacing: float = 100.0,
    min_cells: int = 3,
    use_rna_clusters: bool = True,
    n_clusters: int = 10,
    simplify_cell_types: bool = True,
    threshold_method: str = "percentile",
    percentile: float = 50.0,
) -> dict:
    """
    Complete pipeline to create region datasets from Xenium data.

    By default, uses RNA-based clustering for cell type annotation (RECOMMENDED).
    This avoids circular logic where protein markers are used both for ground truth
    and for deconvolution.

    Reference:
        Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
        Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
        https://doi.org/10.1186/s12859-025-06044-0

    Args:
        data_dir: Path to Xenium data directory
        output_dir: Path to output directory
        n_regions: Number of regions to create
        spot_diameter: Diameter of pseudo-spots (µm)
        center_spacing: Center-to-center spacing (µm)
        min_cells: Minimum cells per spot
        use_rna_clusters: If True (default), use RNA-based clustering for ground truth.
                          If False, use protein thresholding (NOT RECOMMENDED).
        n_clusters: Number of RNA clusters to use (default 10)
        simplify_cell_types: If True, merge related cell types (e.g., Endothelial_2 -> Endothelial)
        threshold_method: Method for protein thresholding (only used if use_rna_clusters=False)
        percentile: Percentile for thresholding (only used if use_rna_clusters=False)

    Returns:
        Dict with summary and paths to saved files
    """
    try:
        from .load_xenium import load_xenium_data, split_gex_protein
        from .create_pseudo_spots import create_pseudo_visium_spots, save_cell_to_spot_mapping
        from .generate_ground_truth import calculate_spot_proportions, add_spatial_coordinates
        from .rna_cell_types import load_rna_cell_types, calculate_spot_proportions_rna
    except ImportError:
        from load_xenium import load_xenium_data, split_gex_protein
        from create_pseudo_spots import create_pseudo_visium_spots, save_cell_to_spot_mapping
        from generate_ground_truth import calculate_spot_proportions, add_spatial_coordinates
        from rna_cell_types import load_rna_cell_types, calculate_spot_proportions_rna

    # Only import protein-based classification if needed
    if not use_rna_clusters:
        try:
            from .define_cell_types import XENIUM_CELL_PROFILE_DICT, classify_cells_by_protein
        except ImportError:
            from define_cell_types import XENIUM_CELL_PROFILE_DICT, classify_cells_by_protein

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading Xenium data...")
    adata = load_xenium_data(data_dir)
    adata_gex_cells, adata_protein_cells = split_gex_protein(adata)

    # Create pseudo-spots
    logger.info("Creating pseudo-Visium spots...")
    adata_spots, cell_to_spot_df = create_pseudo_visium_spots(
        adata,
        spot_diameter=spot_diameter,
        center_spacing=center_spacing,
        min_cells=min_cells,
    )

    # Split spots by feature type
    gex_mask = adata_spots.var["feature_types"] == "Gene Expression"
    protein_mask = adata_spots.var["feature_types"] == "Protein Expression"
    adata_gex_spots = adata_spots[:, gex_mask].copy()
    adata_protein_spots = adata_spots[:, protein_mask].copy()

    # Prepare for CITEgeist
    adata_gex_spots, adata_protein_spots = prepare_citegeist_input(
        adata_gex_spots, adata_protein_spots
    )

    # Classify cells and calculate proportions
    if use_rna_clusters:
        # RECOMMENDED: RNA-based clustering avoids circular logic
        # Reference: Zhao et al. (2025). BMC Bioinformatics, 26(1), 25.
        logger.info(f"Using RNA-based clustering for ground truth (n_clusters={n_clusters})...")
        cell_types = load_rna_cell_types(
            data_dir,
            n_clusters=n_clusters,
            simplify=simplify_cell_types,
        )

        proportions_df = calculate_spot_proportions_rna(
            cell_to_spot_df,
            cell_types,
            list(adata_spots.obs_names),
        )
    else:
        # NOT RECOMMENDED: Protein thresholding creates circular logic
        logger.warning(
            "Using protein thresholding for ground truth. "
            "This is NOT RECOMMENDED as it creates circular logic with protein-based deconvolution. "
            "See: Zhao et al. (2025). BMC Bioinformatics, 26(1), 25."
        )
        cell_types = classify_cells_by_protein(
            adata_protein_cells,
            XENIUM_CELL_PROFILE_DICT,
            threshold_method=threshold_method,
            percentile=percentile,
        )

        proportions_df = calculate_spot_proportions(
            cell_to_spot_df,
            cell_types,
            list(adata_spots.obs_names),
            include_unassigned=True,
        )

    proportions_df = add_spatial_coordinates(proportions_df, adata_spots)

    # Split into regions
    logger.info(f"Splitting into {n_regions} regions...")
    region_masks = split_tissue_regions(adata_spots, n_regions=n_regions)

    # Add region_id to cell_to_spot mapping
    cell_to_spot_df["region_id"] = None
    for region_id, mask in enumerate(region_masks):
        region_spot_names = set(adata_spots.obs_names[mask])
        cell_mask = cell_to_spot_df["spot_id"].isin(region_spot_names)
        cell_to_spot_df.loc[cell_mask, "region_id"] = region_id

    # Save cell-to-spot mapping
    save_cell_to_spot_mapping(cell_to_spot_df, output_dir)

    # Save region datasets
    logger.info("Saving region datasets...")
    saved_files = save_region_datasets(
        adata_gex_spots,
        adata_protein_spots,
        proportions_df,
        region_masks,
        output_dir,
    )

    # Summary
    # Get cell types from proportions columns (excluding metadata columns)
    metadata_cols = ["n_cells", "x", "y", "spot_x", "spot_y"]
    cell_type_cols = [c for c in proportions_df.columns if c not in metadata_cols]

    summary = {
        "data_dir": str(data_dir),
        "output_dir": str(output_dir),
        "n_regions": n_regions,
        "total_spots": len(adata_spots),
        "total_cells": len(adata),
        "n_genes": adata_gex_spots.shape[1],
        "n_proteins": adata_protein_spots.shape[1],
        "ground_truth_method": "rna_clustering" if use_rna_clusters else "protein_thresholding",
        "ground_truth_reference": "Zhao et al. (2025). BMC Bioinformatics, 26(1), 25. https://doi.org/10.1186/s12859-025-06044-0",
        "cell_types": cell_type_cols,
        "regions": saved_files,
    }

    # Save summary
    summary_path = output_dir / "dataset_summary.json"
    import json
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Dataset creation complete. Summary saved to {summary_path}")

    return summary


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"
    OUTPUT_DIR = "Benchmarking/xenium_benchmarking/data"

    sys.path.insert(0, "Benchmarking/xenium_benchmarking/src")

    # Using RNA-based clustering by default (RECOMMENDED)
    # Reference: Zhao et al. (2025). BMC Bioinformatics, 26(1), 25.
    # https://doi.org/10.1186/s12859-025-06044-0
    print("Creating all region datasets with RNA-based ground truth...")
    summary = create_all_region_datasets(
        data_dir=DATA_DIR,
        output_dir=OUTPUT_DIR,
        n_regions=5,
        spot_diameter=55.0,
        center_spacing=100.0,
        min_cells=3,
        use_rna_clusters=True,  # Default: RNA-based clustering (recommended)
        n_clusters=10,
        simplify_cell_types=True,
    )

    print("\nSummary:")
    print(f"  Ground truth method: {summary['ground_truth_method']}")
    print(f"  Total spots: {summary['total_spots']}")
    print(f"  Total cells: {summary['total_cells']}")
    print(f"  Genes: {summary['n_genes']}")
    print(f"  Proteins: {summary['n_proteins']}")
    print(f"  Cell types: {summary['cell_types']}")
    print(f"  Regions: {summary['n_regions']}")
    for region in summary["regions"]:
        print(f"    Region {region['region_id']}: {region['n_spots']} spots")
