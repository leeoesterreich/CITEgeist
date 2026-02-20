"""Tests for nuclear morphology feature extraction."""
import numpy as np
import pytest
from CITEgeist.model.morphology_features import extract_nucleus_features


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
