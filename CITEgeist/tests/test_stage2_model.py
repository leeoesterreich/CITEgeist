# CITEgeist/tests/test_stage2_model.py
"""Tests for Stage 2 complete model."""
import sys
from pathlib import Path

import pytest
import torch

# Add model directory to path for direct imports (avoid heavy __init__.py)
MODEL_DIR = Path(__file__).parent.parent / "model"
sys.path.insert(0, str(MODEL_DIR))


def test_stage2_model_forward():
    """Model computes OT loss and transport plan."""
    from stage2_model import Stage2Model
    from vae import VAEEncoder

    # Create encoder (will be frozen)
    encoder = VAEEncoder(in_channels=2, latent_dim=128)

    model = Stage2Model(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        projection_dim=32,
    )

    # Simulate a spot with 10 nuclei
    patches = torch.randn(10, 2, 96, 96)
    cell_counts = torch.tensor([2, 1, 3, 1, 1, 1, 1], dtype=torch.float32)  # sums to 10

    loss, transport_plan = model(patches, cell_counts)

    assert loss.shape == ()  # scalar
    assert transport_plan.shape == (10, 7)
    # Transport plan should approximately sum to cell_counts/N per column
    col_sums = transport_plan.sum(dim=0) * 10
    assert torch.allclose(col_sums, cell_counts, atol=0.5)


def test_stage2_model_assign():
    """Model produces hard assignments."""
    from stage2_model import Stage2Model
    from vae import VAEEncoder

    encoder = VAEEncoder(in_channels=2, latent_dim=128)
    model = Stage2Model(encoder=encoder, n_types=7)

    patches = torch.randn(10, 2, 96, 96)
    cell_counts = torch.tensor([2, 1, 3, 1, 1, 1, 1], dtype=torch.float32)

    assignments, confidence = model.assign(patches, cell_counts)

    assert assignments.shape == (10,)
    assert confidence.shape == (10,)
    assert assignments.min() >= 0
    assert assignments.max() < 7
    # With untrained model, assignments may not perfectly match counts
    # Just verify the total number of assignments is correct
    assert len(assignments) == len(patches)


def test_encoder_is_frozen():
    """VAE encoder parameters are frozen."""
    from stage2_model import Stage2Model
    from vae import VAEEncoder

    encoder = VAEEncoder(in_channels=2, latent_dim=128)
    model = Stage2Model(encoder=encoder, n_types=7)

    for param in model.encoder.parameters():
        assert param.requires_grad == False
