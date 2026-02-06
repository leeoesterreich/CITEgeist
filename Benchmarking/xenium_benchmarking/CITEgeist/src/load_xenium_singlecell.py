"""
Load Xenium single-cell data for CITEgeist analysis.

This module provides functions to load Xenium single-cell data (not aggregated
into pseudo-Visium spots) for direct application of CITEgeist modules.

Key difference from pseudo-Visium:
- Each observation is a SINGLE CELL, not a spot containing multiple cells
- Modules 1-2 work directly on single cells (resolution-independent)
- Module 3 becomes "soft profile assignment" (no deconvolution needed)
- Module 4 discovers programs from actual cell-type-specific expression
"""

import json
import logging
import sys
from pathlib import Path
from typing import Tuple, Optional, List, Dict

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# Add xenium_pseudovisium src to path for existing loaders
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "xenium_pseudovisium" / "src"))

logger = logging.getLogger(__name__)

# Default data directory
XENIUM_DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"


def load_xenium_singlecell(
    data_dir: str = XENIUM_DATA_DIR,
    region_id: Optional[int] = None,
    region_bounds: Optional[Tuple[float, float, float, float]] = None,
    max_cells: Optional[int] = None,
    seed: int = 42,
    min_counts: int = 0,
) -> Tuple[sc.AnnData, sc.AnnData]:
    """
    Load Xenium single-cell data with GEX and protein split.

    This function loads the full single-cell resolution data from Xenium,
    optionally subsetting to a specific region or random sample.

    Args:
        data_dir: Path to Xenium output directory.
        region_id: If provided, load cells from this region (0-4) using
                   pre-defined region boundaries matching pseudo-Visium regions.
        region_bounds: If provided, custom bounds as (x_min, x_max, y_min, y_max).
                       Overrides region_id if both provided.
        max_cells: Maximum cells to load (randomly sampled if data exceeds).
                   Useful for memory management with large datasets.
        seed: Random seed for subsampling.
        min_counts: Minimum total counts per cell to keep.

    Returns:
        Tuple of (adata_gex, adata_protein) with aligned observations.
        - adata_gex: Gene expression (n_cells × n_genes)
        - adata_protein: Protein expression (n_cells × n_proteins)
        - Both have obsm['spatial'] with cell coordinates
    """
    try:
        from load_xenium import load_xenium_data, split_gex_protein
    except ImportError:
        from Benchmarking.xenium_pseudovisium.src.load_xenium import load_xenium_data, split_gex_protein

    logger.info(f"Loading Xenium single-cell data from {data_dir}")

    # Load full single-cell data
    adata = load_xenium_data(data_dir, min_counts=min_counts)
    logger.info(f"Loaded {adata.shape[0]} cells × {adata.shape[1]} features")

    # Get spatial coordinates
    coords = adata.obsm["spatial"]

    # Determine region bounds
    if region_bounds is not None:
        x_min, x_max, y_min, y_max = region_bounds
        logger.info(f"Using custom region bounds: x=[{x_min:.1f}, {x_max:.1f}], y=[{y_min:.1f}, {y_max:.1f}]")
    elif region_id is not None:
        # Calculate region bounds matching pseudo-Visium horizontal strips
        x_coords = coords[:, 0]
        x_full_min, x_full_max = x_coords.min(), x_coords.max()
        n_regions = 5
        boundaries = np.linspace(x_full_min, x_full_max, n_regions + 1)
        x_min = boundaries[region_id]
        x_max = boundaries[region_id + 1]
        y_min = coords[:, 1].min()
        y_max = coords[:, 1].max()
        logger.info(f"Region {region_id} bounds: x=[{x_min:.1f}, {x_max:.1f}], y=[{y_min:.1f}, {y_max:.1f}]")
    else:
        x_min, x_max = None, None
        y_min, y_max = None, None

    # Subset by region if specified
    if x_min is not None:
        x = coords[:, 0]
        y = coords[:, 1]
        region_mask = (x >= x_min) & (x < x_max) & (y >= y_min) & (y <= y_max)

        # Handle last region including boundary
        if region_id is not None and region_id == 4:
            region_mask = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)

        adata = adata[region_mask, :].copy()
        logger.info(f"After region subsetting: {adata.shape[0]} cells")

    # Subsample if max_cells specified
    if max_cells is not None and adata.shape[0] > max_cells:
        rng = np.random.default_rng(seed)
        indices = rng.choice(adata.shape[0], size=max_cells, replace=False)
        indices = np.sort(indices)
        adata = adata[indices, :].copy()
        logger.info(f"After subsampling: {adata.shape[0]} cells")

    # Split into GEX and protein
    adata_gex, adata_protein = split_gex_protein(adata)

    logger.info(f"GEX: {adata_gex.shape[0]} cells × {adata_gex.shape[1]} genes")
    logger.info(f"Protein: {adata_protein.shape[0]} cells × {adata_protein.shape[1]} proteins")

    return adata_gex, adata_protein


def get_region_bounds(
    data_dir: str = XENIUM_DATA_DIR,
    n_regions: int = 5,
) -> List[Tuple[float, float, float, float]]:
    """
    Get region boundaries for the Xenium dataset.

    Returns boundaries matching the pseudo-Visium horizontal strip regions.

    Args:
        data_dir: Path to Xenium output directory.
        n_regions: Number of regions (default: 5).

    Returns:
        List of (x_min, x_max, y_min, y_max) tuples for each region.
    """
    try:
        from load_xenium import load_xenium_data
    except ImportError:
        from Benchmarking.xenium_pseudovisium.src.load_xenium import load_xenium_data

    # Load just coordinates (minimal processing)
    adata = load_xenium_data(data_dir, min_counts=0)
    coords = adata.obsm["spatial"]

    x_coords = coords[:, 0]
    y_coords = coords[:, 1]

    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()

    # Horizontal strip boundaries (matching pseudo-Visium)
    boundaries = np.linspace(x_min, x_max, n_regions + 1)

    region_bounds = []
    for i in range(n_regions):
        region_bounds.append((
            boundaries[i],
            boundaries[i + 1],
            y_min,
            y_max,
        ))

    return region_bounds


def get_quadrant_bounds(
    data_dir: str = XENIUM_DATA_DIR,
) -> List[Tuple[float, float, float, float]]:
    """
    Get quadrant boundaries for the Xenium dataset.

    Quadrants are defined by X/Y midpoint:
        Q0: bottom-left, Q1: bottom-right, Q2: top-left, Q3: top-right

    Args:
        data_dir: Path to Xenium output directory.

    Returns:
        List of (x_min, x_max, y_min, y_max) tuples for each quadrant (Q0-Q3).
    """
    try:
        from load_xenium import load_xenium_data
    except ImportError:
        from Benchmarking.xenium_pseudovisium.src.load_xenium import load_xenium_data

    adata = load_xenium_data(data_dir, min_counts=0)
    coords = adata.obsm["spatial"]

    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    x_mid = (x_min + x_max) / 2
    y_mid = (y_min + y_max) / 2

    return [
        (x_min, x_mid, y_min, y_mid),  # Q0: bottom-left
        (x_mid, x_max, y_min, y_mid),  # Q1: bottom-right
        (x_min, x_mid, y_mid, y_max),  # Q2: top-left
        (x_mid, x_max, y_mid, y_max),  # Q3: top-right
    ]


def load_all_regions(
    data_dir: str = XENIUM_DATA_DIR,
    max_cells_per_region: Optional[int] = None,
    seed: int = 42,
) -> Dict[int, Tuple[sc.AnnData, sc.AnnData]]:
    """
    Load single-cell data for all 5 regions.

    Args:
        data_dir: Path to Xenium output directory.
        max_cells_per_region: Maximum cells per region (for memory management).
        seed: Random seed for subsampling.

    Returns:
        Dict mapping region_id (0-4) to (adata_gex, adata_protein) tuples.
    """
    results = {}
    for region_id in range(5):
        logger.info(f"\n{'='*50}")
        logger.info(f"Loading region {region_id}")
        logger.info(f"{'='*50}")
        adata_gex, adata_protein = load_xenium_singlecell(
            data_dir=data_dir,
            region_id=region_id,
            max_cells=max_cells_per_region,
            seed=seed,
        )
        results[region_id] = (adata_gex, adata_protein)
        logger.info(f"Region {region_id}: {adata_gex.shape[0]} cells")

    return results


def get_dataset_summary(data_dir: str = XENIUM_DATA_DIR) -> Dict:
    """
    Get summary statistics for the Xenium dataset.

    Returns:
        Dict with dataset statistics.
    """
    try:
        from load_xenium import load_xenium_data, split_gex_protein
    except ImportError:
        from Benchmarking.xenium_pseudovisium.src.load_xenium import load_xenium_data, split_gex_protein

    adata = load_xenium_data(data_dir, min_counts=0)
    coords = adata.obsm["spatial"]

    adata_gex, adata_protein = split_gex_protein(adata)

    summary = {
        "total_cells": adata.shape[0],
        "n_genes": adata_gex.shape[1],
        "n_proteins": adata_protein.shape[1],
        "protein_names": list(adata_protein.var_names),
        "spatial_extent": {
            "x_min": float(coords[:, 0].min()),
            "x_max": float(coords[:, 0].max()),
            "y_min": float(coords[:, 1].min()),
            "y_max": float(coords[:, 1].max()),
        },
        "region_bounds": [
            {
                "region_id": i,
                "x_min": b[0],
                "x_max": b[1],
                "y_min": b[2],
                "y_max": b[3],
            }
            for i, b in enumerate(get_region_bounds(data_dir))
        ],
    }

    return summary


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    print("=" * 60)
    print("Xenium Single-Cell Data Loader Test")
    print("=" * 60)

    # Get dataset summary
    print("\nDataset summary:")
    summary = get_dataset_summary()
    print(f"  Total cells: {summary['total_cells']:,}")
    print(f"  Genes: {summary['n_genes']}")
    print(f"  Proteins: {summary['n_proteins']}")
    print(f"  Proteins: {summary['protein_names']}")
    print(f"  Spatial extent: x=[{summary['spatial_extent']['x_min']:.1f}, {summary['spatial_extent']['x_max']:.1f}]")
    print(f"                  y=[{summary['spatial_extent']['y_min']:.1f}, {summary['spatial_extent']['y_max']:.1f}]")

    # Load a single region (region 0) for testing
    print("\n" + "=" * 60)
    print("Loading region 0 (subset of 50,000 cells for testing)")
    print("=" * 60)

    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=0,
        max_cells=50000,
        seed=42,
    )

    print(f"\nGEX shape: {adata_gex.shape}")
    print(f"Protein shape: {adata_protein.shape}")
    print(f"Proteins: {list(adata_protein.var_names)}")
    print(f"Spatial coordinates shape: {adata_gex.obsm['spatial'].shape}")

    # Quick validation: check coordinates are within expected range
    coords = adata_gex.obsm["spatial"]
    print(f"\nSpatial coordinate ranges:")
    print(f"  X: [{coords[:, 0].min():.1f}, {coords[:, 0].max():.1f}]")
    print(f"  Y: [{coords[:, 1].min():.1f}, {coords[:, 1].max():.1f}]")

    print("\nData loader test complete!")
