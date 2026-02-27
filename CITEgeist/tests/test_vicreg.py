# CITEgeist/tests/test_vicreg.py
"""Tests for VICReg loss module."""
import torch
import pytest


def test_vicreg_loss_returns_components():
    """VICReg should return total loss and components."""
    from CITEgeist.model.vicreg import vicreg_loss

    z = torch.randn(32, 128)  # Batch of 32, 128-dim latents
    z_aug = torch.randn(32, 128)  # Augmented view

    loss, components = vicreg_loss(z, z_aug)

    assert isinstance(loss, torch.Tensor)
    assert loss.dim() == 0  # Scalar
    assert "invariance" in components
    assert "variance" in components
    assert "covariance" in components


def test_vicreg_invariance_zero_for_identical():
    """Invariance loss should be 0 for identical views."""
    from CITEgeist.model.vicreg import vicreg_loss

    z = torch.randn(32, 128)

    _, components = vicreg_loss(z, z)  # Same input twice

    assert components["invariance"].item() < 1e-6


def test_vicreg_variance_penalizes_collapse():
    """Variance loss should be high when all embeddings identical."""
    from CITEgeist.model.vicreg import vicreg_loss

    # Collapsed embeddings: all same vector
    z = torch.ones(32, 128)
    z_aug = torch.ones(32, 128) + torch.randn(32, 128) * 0.01

    _, components = vicreg_loss(z, z_aug)

    # Variance loss should be high (std < 1 in all dims)
    assert components["variance"].item() > 100  # 128 dims * ~1 each


def test_vicreg_covariance_penalizes_correlation():
    """Covariance loss should penalize correlated dimensions."""
    from CITEgeist.model.vicreg import vicreg_loss

    # Create embeddings with correlated dimensions
    base = torch.randn(32, 1)
    z = base.expand(32, 128).clone()  # All dims identical = highly correlated
    z_aug = z + torch.randn(32, 128) * 0.1

    _, components = vicreg_loss(z, z_aug)

    # Covariance loss should be high
    assert components["covariance"].item() > 10


def test_vicreg_batch_size_one_raises():
    """VICReg should raise error for batch_size=1."""
    from CITEgeist.model.vicreg import vicreg_loss

    z = torch.randn(1, 128)  # Single sample
    z_aug = torch.randn(1, 128)

    try:
        vicreg_loss(z, z_aug)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "batch_size >= 2" in str(e)
