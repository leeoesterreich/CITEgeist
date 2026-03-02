"""Create pseudo-Visium spots from Visium HD single-cell data.

Aggregates single-cell data into Visium-like spots for benchmarking
morphology-based cell type assignment.
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from typing import Tuple


def create_hex_grid(
    bounds: Tuple[float, float, float, float],
    spacing_um: float,
    pixel_size: float,
) -> pd.DataFrame:
    """Create hexagonal grid of pseudo-Visium spots.

    Args:
        bounds: (x_min, x_max, y_min, y_max) in pixels
        spacing_um: Center-to-center spacing in microns
        pixel_size: Microns per pixel

    Returns:
        DataFrame with spot_id, x, y columns
    """
    x_min, x_max, y_min, y_max = bounds
    spacing_px = spacing_um / pixel_size

    # Hex grid: offset every other row by half spacing
    spots = []
    spot_id = 0

    y = y_min + spacing_px / 2
    row = 0
    while y < y_max:
        x_offset = (spacing_px / 2) if row % 2 else 0
        x = x_min + spacing_px / 2 + x_offset
        while x < x_max:
            spots.append({'spot_id': spot_id, 'x': x, 'y': y})
            spot_id += 1
            x += spacing_px
        y += spacing_px * np.sqrt(3) / 2  # Hex row spacing
        row += 1

    return pd.DataFrame(spots)


def assign_cells_to_spots(
    cells: pd.DataFrame,
    spots: pd.DataFrame,
    radius_um: float,
    pixel_size: float,
) -> pd.DataFrame:
    """Assign cells to nearest spot within radius.

    Args:
        cells: DataFrame with cell_id, x, y, cell_type columns
        spots: DataFrame with spot_id, x, y columns
        radius_um: Maximum distance from spot center in microns
        pixel_size: Microns per pixel

    Returns:
        cells DataFrame with added spot_id column
    """
    radius_px = radius_um / pixel_size

    # Build KD-tree for spots
    spot_coords = spots[['x', 'y']].values
    tree = cKDTree(spot_coords)

    # Query nearest spot for each cell
    cell_coords = cells[['x', 'y']].values
    distances, indices = tree.query(cell_coords)

    # Assign to spot if within radius
    result = cells.copy()
    result['spot_id'] = -1  # Default: unassigned
    mask = distances <= radius_px
    result.loc[mask, 'spot_id'] = spots.iloc[indices[mask]]['spot_id'].values

    # Remove unassigned cells
    result = result[result['spot_id'] >= 0]

    return result


def compute_spot_proportions(
    mapping: pd.DataFrame,
    min_cells: int = 5,
) -> pd.DataFrame:
    """Compute cell type proportions per spot.

    Args:
        mapping: DataFrame with spot_id, cell_type columns
        min_cells: Minimum cells per spot (filter smaller spots)

    Returns:
        DataFrame with spot_id as index, cell types as columns, proportions as values
    """
    # Count cells per spot
    spot_counts = mapping.groupby('spot_id').size()
    valid_spots = spot_counts[spot_counts >= min_cells].index

    # Filter to valid spots
    filtered = mapping[mapping['spot_id'].isin(valid_spots)]

    # Compute proportions
    proportions = filtered.groupby(['spot_id', 'cell_type']).size().unstack(fill_value=0)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    return proportions
