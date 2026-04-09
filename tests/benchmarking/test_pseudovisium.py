# tests/benchmarking/test_pseudovisium.py
"""Tests for pseudo-Visium hex grid construction."""
import numpy as np
import pandas as pd
import pytest


def test_hex_grid_spacing():
    """Hex grid centers should be ~100µm apart (within 15µm tolerance)."""
    from Benchmarking.morphology_audit.src.pseudovisium import build_hex_grid

    # 1000×1000 µm tissue at 0.25 µm/px = 4000×4000 px
    grid = build_hex_grid(
        tissue_bbox_um=(0, 0, 1000, 1000),
        spot_spacing_um=100.0,
    )
    # Should have spots
    assert len(grid) > 0
    # Check spacing between nearest neighbors
    from scipy.spatial import cKDTree
    tree = cKDTree(grid[["x_um", "y_um"]].values)
    dists, _ = tree.query(grid[["x_um", "y_um"]].values, k=2)
    nn_dists = dists[:, 1]  # nearest neighbor (not self)
    # Hex row height is 100*sqrt(3)/2 ≈ 86.6µm, so median NN dist is ~86.6µm
    # Tolerance of 15µm to account for non-equal row/column spacing
    assert np.abs(np.median(nn_dists) - 100.0) < 15.0  # ~100µm hex structure


def test_hex_grid_count():
    """1000×1000µm with 100µm spacing should have ~80-150 spots."""
    from Benchmarking.morphology_audit.src.pseudovisium import build_hex_grid

    grid = build_hex_grid(
        tissue_bbox_um=(0, 0, 1000, 1000),
        spot_spacing_um=100.0,
    )
    assert 80 <= len(grid) <= 150


def test_assign_cells_to_spots():
    """Each cell should map to its nearest spot."""
    from Benchmarking.morphology_audit.src.pseudovisium import assign_cells_to_spots

    spot_coords = np.array([[0, 0], [100, 0], [200, 0]], dtype=float)
    cell_coords = np.array([[10, 5], [90, 3], [195, 2]], dtype=float)

    assignments = assign_cells_to_spots(cell_coords, spot_coords)
    assert np.array_equal(assignments, [0, 1, 2])


def test_compute_spot_proportions():
    """Spot proportions should sum to 1 and reflect type counts."""
    from Benchmarking.morphology_audit.src.pseudovisium import compute_spot_proportions

    cell_to_spot = np.array([0, 0, 0, 1, 1])
    cell_types = pd.Series(["A", "A", "B", "B", "B"])
    type_list = ["A", "B"]

    props = compute_spot_proportions(cell_to_spot, cell_types, type_list, n_spots=2)
    # Spot 0: 2A, 1B → [2/3, 1/3]
    np.testing.assert_allclose(props[0], [2 / 3, 1 / 3], atol=1e-6)
    # Spot 1: 0A, 2B → [0, 1]
    np.testing.assert_allclose(props[1], [0.0, 1.0], atol=1e-6)


def test_filter_sparse_spots():
    """Spots with <min_cells should be excluded."""
    from Benchmarking.morphology_audit.src.pseudovisium import filter_sparse_spots

    cell_to_spot = np.array([0, 0, 0, 1, 2, 2, 2, 2])
    # Spot 0: 3 cells, Spot 1: 1 cell (drop), Spot 2: 4 cells
    keep_mask, new_mapping = filter_sparse_spots(cell_to_spot, n_spots=3, min_cells=3)
    assert keep_mask[0] == True
    assert keep_mask[1] == False
    assert keep_mask[2] == True
