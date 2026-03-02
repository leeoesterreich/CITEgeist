"""Tests for pseudo-Visium spot creation."""
import numpy as np
import pandas as pd
import pytest
from create_pseudo_visium import (
    create_hex_grid,
    assign_cells_to_spots,
    compute_spot_proportions,
)


def test_create_hex_grid_basic():
    """Test hex grid generation covers bounding box."""
    bounds = (0, 1000, 0, 1000)  # x_min, x_max, y_min, y_max
    spacing_um = 100
    pixel_size = 0.5  # um/pixel

    spots = create_hex_grid(bounds, spacing_um, pixel_size)

    assert isinstance(spots, pd.DataFrame)
    assert 'spot_id' in spots.columns
    assert 'x' in spots.columns
    assert 'y' in spots.columns
    assert len(spots) > 0
    # spacing_px = 100um / 0.5(um/px) = 200px
    # Roughly (1000/200)^2 = 25 spots for square grid; hex is similar
    assert 15 < len(spots) < 50


def test_assign_cells_to_spots():
    """Test cell-to-spot assignment."""
    # Create mock cells
    cells = pd.DataFrame({
        'cell_id': [1, 2, 3, 4],
        'x': [50, 55, 150, 155],
        'y': [50, 55, 50, 55],
        'cell_type': ['A', 'A', 'B', 'B'],
    })

    # Create mock spots
    spots = pd.DataFrame({
        'spot_id': [0, 1],
        'x': [50, 150],
        'y': [50, 50],
    })

    mapping = assign_cells_to_spots(cells, spots, radius_um=30, pixel_size=1.0)

    assert 'spot_id' in mapping.columns
    # Cells 1,2 should be in spot 0; cells 3,4 in spot 1
    assert mapping[mapping['cell_id'] == 1]['spot_id'].iloc[0] == 0
    assert mapping[mapping['cell_id'] == 3]['spot_id'].iloc[0] == 1


def test_compute_spot_proportions():
    """Test proportion computation."""
    mapping = pd.DataFrame({
        'cell_id': [1, 2, 3, 4, 5],
        'spot_id': [0, 0, 0, 1, 1],
        'cell_type': ['A', 'A', 'B', 'B', 'B'],
    })

    proportions = compute_spot_proportions(mapping, min_cells=2)

    assert 0 in proportions.index
    assert 1 in proportions.index
    # Spot 0: 2A, 1B -> A=0.67, B=0.33
    assert abs(proportions.loc[0, 'A'] - 0.667) < 0.01
    assert abs(proportions.loc[0, 'B'] - 0.333) < 0.01
