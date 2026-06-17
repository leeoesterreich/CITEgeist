"""Tests for nuclear morphology feature extraction."""
import numpy as np
import pytest
from CITEgeist.model.morphology.morphology_features import extract_nucleus_features, largest_remainder_discretize


def test_extract_features_single_nucleus():
    """Test feature extraction from a single circular nucleus."""
    # Create 50x50 mask with single circular nucleus (label=1)
    mask = np.zeros((50, 50), dtype=np.int32)
    center = (25, 25)
    radius = 10
    y, x = np.ogrid[:50, :50]
    circle = (x - center[1])**2 + (y - center[0])**2 <= radius**2
    mask[circle] = 1

    features_df = extract_nucleus_features(mask)

    assert len(features_df) == 1
    assert 'nucleus_id' in features_df.columns
    assert 'area' in features_df.columns
    assert 'circularity' in features_df.columns
    assert 'eccentricity' in features_df.columns
    assert 'solidity' in features_df.columns
    assert 'centroid_x' in features_df.columns
    assert 'centroid_y' in features_df.columns
    # Circle should have high circularity (close to 1)
    assert features_df['circularity'].iloc[0] > 0.9


def test_extract_features_multiple_nuclei():
    """Test feature extraction from multiple nuclei."""
    mask = np.zeros((100, 100), dtype=np.int32)
    # Nucleus 1: circle at (25, 25)
    y, x = np.ogrid[:100, :100]
    mask[((x - 25)**2 + (y - 25)**2) <= 64] = 1
    # Nucleus 2: circle at (75, 75)
    mask[((x - 75)**2 + (y - 75)**2) <= 100] = 2

    features_df = extract_nucleus_features(mask)

    assert len(features_df) == 2
    assert set(features_df['nucleus_id']) == {1, 2}
    # Nucleus 2 should have larger area (r=10 vs r=8)
    area_1 = features_df[features_df['nucleus_id'] == 1]['area'].iloc[0]
    area_2 = features_df[features_df['nucleus_id'] == 2]['area'].iloc[0]
    assert area_2 > area_1


def test_extract_features_empty_mask():
    """Test handling of empty mask."""
    mask = np.zeros((50, 50), dtype=np.int32)
    features_df = extract_nucleus_features(mask)
    assert len(features_df) == 0


def test_extract_features_elongated_nucleus():
    """Test that elongated nucleus has lower circularity and higher eccentricity."""
    mask = np.zeros((50, 100), dtype=np.int32)
    # Elongated ellipse
    mask[20:30, 20:80] = 1

    features_df = extract_nucleus_features(mask)

    assert len(features_df) == 1
    # Elongated shape should have lower circularity
    assert features_df['circularity'].iloc[0] < 0.7
    # And higher eccentricity
    assert features_df['eccentricity'].iloc[0] > 0.8


# --- Largest Remainder Discretization Tests ---

def test_largest_remainder_exact():
    """Test discretization when proportions divide evenly."""
    proportions = np.array([0.5, 0.5, 0.0])
    n_total = 4
    counts = largest_remainder_discretize(proportions, n_total)
    assert list(counts) == [2, 2, 0]
    assert counts.sum() == n_total


def test_largest_remainder_rounding():
    """Test discretization with remainders."""
    proportions = np.array([0.4, 0.35, 0.25])
    n_total = 5
    counts = largest_remainder_discretize(proportions, n_total)
    # 0.4*5=2.0, 0.35*5=1.75, 0.25*5=1.25
    # Floor: [2, 1, 1] = 4, need 1 more
    # Remainders: [0.0, 0.75, 0.25] -> give to index 1
    assert list(counts) == [2, 2, 1]
    assert counts.sum() == n_total


def test_largest_remainder_single_cell():
    """Test with single cell - assigns to highest proportion."""
    proportions = np.array([0.4, 0.35, 0.25])
    n_total = 1
    counts = largest_remainder_discretize(proportions, n_total)
    assert list(counts) == [1, 0, 0]
    assert counts.sum() == n_total


def test_largest_remainder_zero_total():
    """Test with zero cells."""
    proportions = np.array([0.5, 0.3, 0.2])
    n_total = 0
    counts = largest_remainder_discretize(proportions, n_total)
    assert list(counts) == [0, 0, 0]
    assert counts.sum() == 0


# --- Cell Morphology Feature Tests ---

def test_extract_cell_features_single_cell():
    """Test cell feature extraction from nucleus + cell masks."""
    from CITEgeist.model.morphology.morphology_features import extract_cell_features

    # Create nucleus (small circle) and cell (larger circle) masks
    nucleus_mask = np.zeros((50, 50), dtype=np.int32)
    cell_mask = np.zeros((50, 50), dtype=np.int32)
    y, x = np.ogrid[:50, :50]

    # Nucleus: r=5 at center
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 25] = 1
    # Cell: r=15 at center (same label)
    cell_mask[((x - 25)**2 + (y - 25)**2) <= 225] = 1

    features_df = extract_cell_features(nucleus_mask, cell_mask)

    assert len(features_df) == 1
    assert 'cell_id' in features_df.columns
    # Nuclear features
    assert 'nucleus_area' in features_df.columns
    assert 'nucleus_circularity' in features_df.columns
    # Cell features
    assert 'cell_area' in features_df.columns
    assert 'cell_circularity' in features_df.columns
    # Ratio features
    assert 'nc_ratio' in features_df.columns
    assert 'cytoplasm_area' in features_df.columns

    # Cell area should be larger than nucleus
    row = features_df.iloc[0]
    assert row['cell_area'] > row['nucleus_area']
    # N:C ratio should be < 1
    assert 0 < row['nc_ratio'] < 1
    # Cytoplasm = cell - nucleus
    assert row['cytoplasm_area'] == row['cell_area'] - row['nucleus_area']


def test_extract_cell_features_multiple_cells():
    """Test cell features with multiple cells of different shapes."""
    from CITEgeist.model.morphology.morphology_features import extract_cell_features

    nucleus_mask = np.zeros((100, 100), dtype=np.int32)
    cell_mask = np.zeros((100, 100), dtype=np.int32)
    y, x = np.ogrid[:100, :100]

    # Cell 1: Small nucleus, large cell (low N:C ratio - like macrophage)
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 16] = 1  # r=4
    cell_mask[((x - 25)**2 + (y - 25)**2) <= 225] = 1    # r=15

    # Cell 2: Large nucleus, small cell (high N:C ratio - like lymphocyte)
    nucleus_mask[((x - 75)**2 + (y - 75)**2) <= 36] = 2  # r=6
    cell_mask[((x - 75)**2 + (y - 75)**2) <= 64] = 2     # r=8

    features_df = extract_cell_features(nucleus_mask, cell_mask)

    assert len(features_df) == 2

    cell1 = features_df[features_df['cell_id'] == 1].iloc[0]
    cell2 = features_df[features_df['cell_id'] == 2].iloc[0]

    # Cell 1 should have lower N:C ratio (more cytoplasm)
    assert cell1['nc_ratio'] < cell2['nc_ratio']
    # Cell 1 should have larger cytoplasm area
    assert cell1['cytoplasm_area'] > cell2['cytoplasm_area']
