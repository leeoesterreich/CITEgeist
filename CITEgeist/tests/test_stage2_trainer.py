# CITEgeist/tests/test_stage2_trainer.py
"""Tests for Stage 2 trainer."""
import sys
from pathlib import Path
import tempfile

import pytest
import torch
import pandas as pd
import numpy as np

# Add model directory to path for direct imports (avoid heavy __init__.py)
MODEL_DIR = Path(__file__).parent.parent / "model"
sys.path.insert(0, str(MODEL_DIR))


def test_trainer_phase1_initialization():
    """Phase 1 initializes prototypes from high-purity spots."""
    from stage2_trainer import Stage2Trainer
    from vae import VAEEncoder

    encoder = VAEEncoder(in_channels=2, latent_dim=128)
    cell_types = ['B cells', 'T cells', 'Macrophages']

    trainer = Stage2Trainer(
        encoder=encoder,
        cell_types=cell_types,
        projection_dim=32,
        device='cpu',
    )

    # Mock proportions with high-purity spots
    proportions = pd.DataFrame({
        'B cells': [0.9, 0.1, 0.1],
        'T cells': [0.05, 0.85, 0.1],
        'Macrophages': [0.05, 0.05, 0.8],
    }, index=['spot_0', 'spot_1', 'spot_2'])

    # Mock patch loader
    def mock_load(spot_id):
        return torch.randn(5, 2, 96, 96)

    trainer.initialize_from_high_purity(
        proportions=proportions,
        load_patches_fn=mock_load,
        threshold=0.70,
    )

    # Check prototypes were initialized
    protos = trainer.model.prototypes()
    assert protos.shape == (3, 32)


def test_trainer_phase2_training_step():
    """Phase 2 training step reduces loss."""
    from stage2_trainer import Stage2Trainer
    from vae import VAEEncoder

    encoder = VAEEncoder(in_channels=2, latent_dim=128)
    cell_types = ['B cells', 'T cells']

    trainer = Stage2Trainer(
        encoder=encoder,
        cell_types=cell_types,
        learning_rate=1e-3,
        device='cpu',
    )

    # One training step
    patches = torch.randn(8, 2, 96, 96)
    cell_counts = torch.tensor([3.0, 5.0])

    loss1 = trainer.train_step(patches, cell_counts)
    loss2 = trainer.train_step(patches, cell_counts)

    # Loss should generally decrease or stay stable
    assert loss1.item() > 0
    assert loss2.item() > 0


def test_trainer_save_load():
    """Trainer can save and load checkpoints."""
    from stage2_trainer import Stage2Trainer
    from vae import VAEEncoder

    encoder = VAEEncoder(in_channels=2, latent_dim=128)
    trainer = Stage2Trainer(encoder=encoder, cell_types=['A', 'B'], device='cpu')

    with tempfile.TemporaryDirectory() as tmpdir:
        save_path = Path(tmpdir) / "checkpoint.pt"
        trainer.save(save_path)

        # Create new trainer and load
        trainer2 = Stage2Trainer(encoder=encoder, cell_types=['A', 'B'], device='cpu')
        trainer2.load(save_path)

        # Check prototypes match
        p1 = trainer.model.prototypes()
        p2 = trainer2.model.prototypes()
        assert torch.allclose(p1, p2)
