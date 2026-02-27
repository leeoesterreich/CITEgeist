# CITEgeist/tests/test_train_vae_vicreg.py
"""Tests for VICReg-enhanced VAE training."""
import torch
import pytest
import tempfile
import numpy as np
from pathlib import Path


def test_train_vae_with_vicreg_flag():
    """Training with --use-vicreg should include VICReg loss."""
    from CITEgeist.model.train_vae import train_vae

    # Create temporary patches
    with tempfile.TemporaryDirectory() as tmpdir:
        patches_dir = Path(tmpdir) / "patches"
        patches_dir.mkdir()
        output_dir = Path(tmpdir) / "output"

        # Create dummy patches
        patches = np.random.randn(100, 2, 96, 96).astype(np.float32)
        np.save(patches_dir / "test_patches.npy", patches)

        # Train with VICReg
        history = train_vae(
            patches_dir=str(patches_dir),
            output_dir=str(output_dir),
            in_channels=2,
            latent_dim=32,
            epochs=2,
            batch_size=16,
            device="cpu",
            use_vicreg=True,
        )

        # Should have VICReg loss components in history
        assert "vicreg_loss" in history
        assert "vicreg_invariance" in history
        assert "vicreg_variance" in history
        assert "vicreg_covariance" in history


def test_train_vae_without_vicreg_unchanged():
    """Training without --use-vicreg should work as before."""
    from CITEgeist.model.train_vae import train_vae

    with tempfile.TemporaryDirectory() as tmpdir:
        patches_dir = Path(tmpdir) / "patches"
        patches_dir.mkdir()
        output_dir = Path(tmpdir) / "output"

        patches = np.random.randn(100, 2, 96, 96).astype(np.float32)
        np.save(patches_dir / "test_patches.npy", patches)

        history = train_vae(
            patches_dir=str(patches_dir),
            output_dir=str(output_dir),
            in_channels=2,
            latent_dim=32,
            epochs=2,
            batch_size=16,
            device="cpu",
            use_vicreg=False,
        )

        # Should NOT have VICReg components
        assert "vicreg_loss" not in history
        assert "loss" in history
        assert "recon_loss" in history
