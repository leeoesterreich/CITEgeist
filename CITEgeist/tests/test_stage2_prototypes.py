# CITEgeist/tests/test_stage2_prototypes.py
"""Tests for Stage 2 type prototypes."""
import sys
from pathlib import Path

import pytest
import torch

# Add model directory to path for direct imports (avoid heavy __init__.py)
MODEL_DIR = Path(__file__).parent.parent / "model"
sys.path.insert(0, str(MODEL_DIR))


def test_prototypes_initialization():
    """Prototypes initialize with correct shape and unit norm."""
    from stage2_prototypes import Stage2Prototypes

    protos = Stage2Prototypes(n_types=7, projection_dim=32)

    P = protos()  # Get prototype matrix
    assert P.shape == (7, 32)

    # Check unit norm
    norms = torch.norm(P, dim=1)
    assert torch.allclose(norms, torch.ones(7), atol=1e-5)


def test_prototypes_from_centroids():
    """Prototypes can be initialized from embedding centroids."""
    from stage2_prototypes import Stage2Prototypes

    # Simulate centroids from high-purity spots
    centroids = torch.randn(7, 32)
    centroids = torch.nn.functional.normalize(centroids, dim=1)

    protos = Stage2Prototypes(n_types=7, projection_dim=32)
    protos.init_from_centroids(centroids)

    P = protos()
    assert torch.allclose(P, centroids, atol=1e-5)


def test_cosine_distances():
    """Compute cosine distances from embeddings to prototypes."""
    from stage2_prototypes import Stage2Prototypes

    protos = Stage2Prototypes(n_types=7, projection_dim=32)
    embeddings = torch.randn(16, 32)
    embeddings = torch.nn.functional.normalize(embeddings, dim=1)

    distances = protos.cosine_distances(embeddings)

    assert distances.shape == (16, 7)
    # Cosine distance range is [0, 2]
    assert distances.min() >= 0
    assert distances.max() <= 2
