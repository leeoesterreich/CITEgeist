"""
Create pseudo-Visium spots from Xenium single-cell data.

This module provides functions to create a hexagonal grid matching
Visium geometry and aggregate single cells into pseudo-spots.
"""

import logging
from pathlib import Path
from typing import Tuple, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)

# Visium geometry constants (in micrometers)
VISIUM_SPOT_DIAMETER = 55.0  # µm
VISIUM_CENTER_SPACING = 100.0  # µm center-to-center


def create_hexagonal_grid(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    spot_diameter: float = VISIUM_SPOT_DIAMETER,
    center_spacing: float = VISIUM_CENTER_SPACING,
) -> np.ndarray:
    """
    Generate hexagonal grid coordinates matching Visium geometry.

    Creates a hexagonal (honeycomb) pattern where:
    - Spots are separated by center_spacing in x direction
    - Odd rows are offset by center_spacing/2 in x
    - Row spacing in y is center_spacing * sqrt(3)/2 for regular hexagons

    Args:
        x_min: Minimum x coordinate (µm)
        x_max: Maximum x coordinate (µm)
        y_min: Minimum y coordinate (µm)
        y_max: Maximum y coordinate (µm)
        spot_diameter: Diameter of each spot (µm), default 55
        center_spacing: Center-to-center distance (µm), default 100

    Returns:
        Array of (x, y) spot center coordinates, shape (n_spots, 2)
    """
    # Row spacing for hexagonal pattern
    row_spacing = center_spacing * np.sqrt(3) / 2

    # Generate grid
    spots = []
    row_idx = 0
    y = y_min + spot_diameter / 2  # Start with half-diameter margin

    while y <= y_max - spot_diameter / 2:
        # Offset odd rows by half the center spacing
        x_offset = (center_spacing / 2) if (row_idx % 2 == 1) else 0
        x = x_min + spot_diameter / 2 + x_offset

        while x <= x_max - spot_diameter / 2:
            spots.append([x, y])
            x += center_spacing

        y += row_spacing
        row_idx += 1

    spot_centers = np.array(spots)
    logger.info(
        f"Created hexagonal grid: {len(spot_centers)} spots "
        f"({spot_diameter}µm diameter, {center_spacing}µm spacing)"
    )

    return spot_centers


def assign_cells_to_spots(
    cell_coords: np.ndarray,
    spot_centers: np.ndarray,
    spot_radius: float = VISIUM_SPOT_DIAMETER / 2,
) -> np.ndarray:
    """
    Assign each cell to a spot based on distance to spot center.

    Cells within spot_radius of a spot center are assigned to that spot.
    Cells outside all spots are assigned -1.

    Args:
        cell_coords: (n_cells, 2) array of cell (x, y) coordinates
        spot_centers: (n_spots, 2) array of spot center coordinates
        spot_radius: Radius of each spot (half of diameter)

    Returns:
        Array of spot indices for each cell, -1 if not in any spot
    """
    # Build KD-tree for efficient nearest-neighbor search
    tree = cKDTree(spot_centers)

    # Find nearest spot for each cell
    distances, spot_indices = tree.query(cell_coords, k=1)

    # Mark cells outside spot radius as unassigned (-1)
    assignments = np.where(distances <= spot_radius, spot_indices, -1)

    n_assigned = np.sum(assignments >= 0)
    n_unassigned = np.sum(assignments < 0)
    logger.info(
        f"Cell assignment: {n_assigned} assigned ({100*n_assigned/len(assignments):.1f}%), "
        f"{n_unassigned} unassigned"
    )

    return assignments


def aggregate_counts_per_spot(
    adata_cells: sc.AnnData,
    cell_to_spot_idx: np.ndarray,
    spot_centers: np.ndarray,
    min_cells: int = 1,
) -> Tuple[sc.AnnData, pd.DataFrame]:
    """
    Aggregate gene/protein counts from cells to spots.

    Sums counts from all cells assigned to each spot.

    Args:
        adata_cells: AnnData with single-cell data
        cell_to_spot_idx: Array mapping each cell to a spot index (-1 = unassigned)
        spot_centers: (n_spots, 2) array of spot center coordinates
        min_cells: Minimum cells per spot to keep (default 1)

    Returns:
        Tuple of:
        - AnnData with aggregated spot-level data
        - DataFrame mapping cell_id -> spot_id (with spot_x, spot_y)
    """
    n_spots = len(spot_centers)
    n_features = adata_cells.shape[1]

    # Convert to sparse if not already
    if sparse.issparse(adata_cells.X):
        X_cells = adata_cells.X.tocsr()
    else:
        X_cells = sparse.csr_matrix(adata_cells.X)

    # Aggregate counts per spot
    spot_counts = []
    cells_per_spot = []

    for spot_idx in range(n_spots):
        cell_mask = cell_to_spot_idx == spot_idx
        n_cells = cell_mask.sum()
        cells_per_spot.append(n_cells)

        if n_cells > 0:
            # Sum counts from cells in this spot
            spot_count = np.asarray(X_cells[cell_mask].sum(axis=0)).flatten()
        else:
            spot_count = np.zeros(n_features)

        spot_counts.append(spot_count)

    # Create spot matrix
    X_spots = np.vstack(spot_counts)
    cells_per_spot = np.array(cells_per_spot)

    # Create spot AnnData (with all spots including empty ones)
    all_spot_names = [f"spot_{i}" for i in range(n_spots)]
    adata_spots = sc.AnnData(
        X=sparse.csr_matrix(X_spots),
        obs=pd.DataFrame(
            {
                "n_cells": cells_per_spot,
                "x": spot_centers[:, 0],
                "y": spot_centers[:, 1],
            },
            index=all_spot_names,
        ),
        var=adata_cells.var.copy(),
    )
    adata_spots.obsm["spatial"] = spot_centers.copy()

    # Create cell-to-spot mapping DataFrame BEFORE filtering
    # This maps cell_id -> spot_name (not index!)
    cell_to_spot_df = pd.DataFrame({
        "spot_idx": cell_to_spot_idx,
        "spot_id": [f"spot_{i}" if i >= 0 else None for i in cell_to_spot_idx],
        "spot_x": [spot_centers[i, 0] if i >= 0 else np.nan for i in cell_to_spot_idx],
        "spot_y": [spot_centers[i, 1] if i >= 0 else np.nan for i in cell_to_spot_idx],
    }, index=adata_cells.obs_names)

    # Filter spots with too few cells
    if min_cells > 0:
        keep_mask = adata_spots.obs["n_cells"] >= min_cells
        n_before = len(adata_spots)
        kept_spot_names = set(adata_spots.obs_names[keep_mask])
        adata_spots = adata_spots[keep_mask, :].copy()

        # Update cell_to_spot_df: set spot_id to None for cells in filtered spots
        cell_to_spot_df.loc[~cell_to_spot_df["spot_id"].isin(kept_spot_names), "spot_id"] = None

        logger.info(
            f"Filtered spots: {n_before} -> {len(adata_spots)} "
            f"(min_cells={min_cells})"
        )

    # Summary statistics
    logger.info(
        f"Spot statistics: "
        f"mean={adata_spots.obs['n_cells'].mean():.1f} cells/spot, "
        f"median={adata_spots.obs['n_cells'].median():.1f}, "
        f"range=[{adata_spots.obs['n_cells'].min()}, {adata_spots.obs['n_cells'].max()}]"
    )

    return adata_spots, cell_to_spot_df


def create_pseudo_visium_spots(
    adata_cells: sc.AnnData,
    spot_diameter: float = VISIUM_SPOT_DIAMETER,
    center_spacing: float = VISIUM_CENTER_SPACING,
    min_cells: int = 3,
    margin: float = 0.0,
) -> Tuple[sc.AnnData, pd.DataFrame]:
    """
    High-level function to create pseudo-Visium spots from single-cell data.

    Args:
        adata_cells: AnnData with single-cell data (must have obsm['spatial'])
        spot_diameter: Diameter of each spot (µm)
        center_spacing: Center-to-center distance (µm)
        min_cells: Minimum cells per spot to keep
        margin: Additional margin around tissue boundary

    Returns:
        Tuple of:
        - adata_spots: AnnData with pseudo-spot data
        - cell_to_spot: DataFrame mapping cell_id -> spot_id (with coordinates)
    """
    # Get tissue boundaries from cell coordinates
    spatial = adata_cells.obsm["spatial"]
    x_min, x_max = spatial[:, 0].min() - margin, spatial[:, 0].max() + margin
    y_min, y_max = spatial[:, 1].min() - margin, spatial[:, 1].max() + margin

    logger.info(f"Tissue boundaries: x=[{x_min:.1f}, {x_max:.1f}], y=[{y_min:.1f}, {y_max:.1f}]")

    # Create hexagonal grid
    spot_centers = create_hexagonal_grid(
        x_min, x_max, y_min, y_max,
        spot_diameter=spot_diameter,
        center_spacing=center_spacing,
    )

    # Assign cells to spots (returns array of spot indices)
    spot_radius = spot_diameter / 2
    cell_to_spot_idx = assign_cells_to_spots(spatial, spot_centers, spot_radius)

    # Aggregate counts (returns AnnData and DataFrame mapping)
    adata_spots, cell_to_spot_df = aggregate_counts_per_spot(
        adata_cells, cell_to_spot_idx, spot_centers, min_cells=min_cells
    )

    return adata_spots, cell_to_spot_df


def get_cell_type_counts_per_spot(
    cell_to_spot: Union[np.ndarray, pd.DataFrame],
    cell_types: pd.Series,
    n_spots: Optional[int] = None,
    spot_names: Optional[list] = None,
) -> pd.DataFrame:
    """
    Count cells of each type per spot.

    Args:
        cell_to_spot: Array mapping each cell to a spot index, OR
                      DataFrame with 'spot_id' column mapping cell_id -> spot_name
        cell_types: Series with cell type labels (same length as cell_to_spot)
        n_spots: Total number of spots (required if cell_to_spot is array)
        spot_names: List of spot names to include (if None, infer from mapping)

    Returns:
        DataFrame with cell type counts per spot
    """
    # Get unique cell types
    unique_types = sorted(cell_types.unique())

    # Handle DataFrame format (new)
    if isinstance(cell_to_spot, pd.DataFrame):
        if spot_names is None:
            spot_names = sorted(cell_to_spot["spot_id"].dropna().unique())

        counts = pd.DataFrame(0, index=spot_names, columns=unique_types)

        for spot_name in spot_names:
            cell_mask = cell_to_spot["spot_id"] == spot_name
            if cell_mask.sum() > 0:
                # Get cell types for cells in this spot
                spot_cell_ids = cell_to_spot.index[cell_mask]
                spot_types = cell_types.loc[cell_types.index.isin(spot_cell_ids)]
                type_counts = spot_types.value_counts()
                for ct, count in type_counts.items():
                    if ct in counts.columns:
                        counts.loc[spot_name, ct] = count

    # Handle array format (legacy)
    else:
        if n_spots is None:
            raise ValueError("n_spots required when cell_to_spot is an array")

        counts = pd.DataFrame(
            0,
            index=[f"spot_{i}" for i in range(n_spots)],
            columns=unique_types,
        )

        for spot_idx in range(n_spots):
            cell_mask = cell_to_spot == spot_idx
            if cell_mask.sum() > 0:
                spot_types = cell_types.iloc[cell_mask]
                type_counts = spot_types.value_counts()
                for ct, count in type_counts.items():
                    counts.loc[f"spot_{spot_idx}", ct] = count

    return counts


def save_cell_to_spot_mapping(
    cell_to_spot_df: pd.DataFrame,
    output_dir: Union[str, Path],
    filename: str = "cell_to_spot_mapping.csv",
) -> Path:
    """
    Save cell-to-spot mapping to CSV file.

    Args:
        cell_to_spot_df: DataFrame mapping cell_id -> spot_id
        output_dir: Output directory
        filename: Filename for the mapping file

    Returns:
        Path to saved file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filepath = output_dir / filename
    cell_to_spot_df.to_csv(filepath)
    logger.info(f"Saved cell-to-spot mapping: {filepath}")

    return filepath


def load_cell_to_spot_mapping(
    data_dir: Union[str, Path],
    filename: str = "cell_to_spot_mapping.csv",
) -> pd.DataFrame:
    """
    Load cell-to-spot mapping from CSV file.

    Args:
        data_dir: Directory containing the mapping file
        filename: Filename for the mapping file

    Returns:
        DataFrame with cell_id index and spot_id, spot_x, spot_y columns
    """
    filepath = Path(data_dir) / filename

    if not filepath.exists():
        raise FileNotFoundError(f"Cell-to-spot mapping not found: {filepath}")

    df = pd.read_csv(filepath, index_col=0)
    logger.info(f"Loaded cell-to-spot mapping: {len(df)} cells from {filepath}")

    return df


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    sys.path.insert(0, "Benchmarking/xenium_benchmarking/src")
    from load_xenium import load_xenium_data

    DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"

    print("Loading Xenium data...")
    adata = load_xenium_data(DATA_DIR)

    print("\nCreating pseudo-Visium spots...")
    adata_spots, cell_to_spot = create_pseudo_visium_spots(
        adata,
        spot_diameter=55.0,
        center_spacing=100.0,
        min_cells=3,
    )

    print(f"\nResult: {adata_spots.shape[0]} spots x {adata_spots.shape[1]} features")
    print(f"Cells per spot: mean={adata_spots.obs['n_cells'].mean():.1f}, "
          f"median={adata_spots.obs['n_cells'].median():.1f}")
