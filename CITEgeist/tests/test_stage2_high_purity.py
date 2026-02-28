# CITEgeist/tests/test_stage2_high_purity.py
"""Tests for high-purity spot detection."""
import sys
from pathlib import Path

import pytest
import pandas as pd
import numpy as np
import torch

# Add model directory to path for direct imports (avoid heavy __init__.py)
MODEL_DIR = Path(__file__).parent.parent / "model"
sys.path.insert(0, str(MODEL_DIR))


def test_find_high_purity_spots():
    """Identify spots where one type dominates."""
    from stage2_high_purity import find_high_purity_spots

    # Create mock proportions
    proportions = pd.DataFrame({
        'B cells': [0.8, 0.3, 0.1, 0.75],
        'T cells': [0.1, 0.6, 0.2, 0.15],
        'Macrophages': [0.1, 0.1, 0.7, 0.10],
    }, index=['spot_0', 'spot_1', 'spot_2', 'spot_3'])

    high_purity = find_high_purity_spots(proportions, threshold=0.70)

    assert len(high_purity) == 3  # spots 0, 2, 3 have >70% dominance
    assert high_purity['spot_0'] == 'B cells'
    assert high_purity['spot_2'] == 'Macrophages'
    assert high_purity['spot_3'] == 'B cells'
    assert 'spot_1' not in high_purity  # max is 0.6 < 0.70


def test_collect_embeddings_by_type():
    """Collect VAE embeddings grouped by dominant type."""
    from stage2_high_purity import collect_embeddings_by_type

    high_purity = {'spot_0': 'B cells', 'spot_1': 'T cells'}

    # Mock patch loader
    def mock_load_patches(spot_id):
        n = 5 if spot_id == 'spot_0' else 3
        return torch.randn(n, 2, 96, 96)

    # Mock encoder
    class MockEncoder:
        def __call__(self, x):
            return torch.randn(x.shape[0], 128), None
        def eval(self):
            pass

    embeddings = collect_embeddings_by_type(
        high_purity,
        mock_load_patches,
        MockEncoder(),
        device='cpu',
    )

    assert 'B cells' in embeddings
    assert 'T cells' in embeddings
    assert embeddings['B cells'].shape == (5, 128)
    assert embeddings['T cells'].shape == (3, 128)


def test_compute_type_centroids():
    """Compute centroid embeddings for each type."""
    from stage2_high_purity import compute_type_centroids

    embeddings = {
        'B cells': torch.randn(100, 32),
        'T cells': torch.randn(50, 32),
    }

    centroids = compute_type_centroids(embeddings)

    assert centroids.shape == (2, 32)
    # Should be normalized
    norms = torch.norm(centroids, dim=1)
    assert torch.allclose(norms, torch.ones(2), atol=1e-5)
