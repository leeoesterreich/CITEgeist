# CITEgeist/tests/test_stage2_projection.py
"""Tests for Stage 2 projection head."""
import sys
from pathlib import Path

import pytest
import torch

# Add model directory to path for direct imports (avoid heavy __init__.py)
MODEL_DIR = Path(__file__).parent.parent / "model"
sys.path.insert(0, str(MODEL_DIR))


def test_projection_head_output_shape():
    """Projection head outputs correct dimensions with L2 normalization."""
    from stage2_projection import Stage2ProjectionHead

    head = Stage2ProjectionHead(latent_dim=128, projection_dim=32)
    z = torch.randn(16, 128)  # batch of 16 embeddings

    p = head(z)

    assert p.shape == (16, 32)
    # Check L2 normalized (unit norm)
    norms = torch.norm(p, dim=1)
    assert torch.allclose(norms, torch.ones(16), atol=1e-5)


def test_projection_head_gradient_flow():
    """Gradients flow through projection head."""
    from stage2_projection import Stage2ProjectionHead

    head = Stage2ProjectionHead(latent_dim=128, projection_dim=32)
    z = torch.randn(8, 128, requires_grad=True)

    p = head(z)
    loss = p.sum()
    loss.backward()

    assert z.grad is not None
    assert head.layers[0].weight.grad is not None
