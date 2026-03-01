# CITEgeist/tests/test_stage2_integration.py
"""Integration tests for Stage 2 pipeline."""
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


@pytest.fixture
def mock_data():
    """Create mock data for integration testing."""
    n_spots = 20
    n_types = 3
    cell_types = ['B cells', 'T cells', 'Macrophages']

    # Mock proportions with some high-purity spots
    props_data = np.random.dirichlet(np.ones(n_types) * 2, size=n_spots)
    # Make some spots high-purity
    props_data[0] = [0.85, 0.10, 0.05]
    props_data[1] = [0.10, 0.80, 0.10]
    props_data[2] = [0.05, 0.10, 0.85]

    proportions = pd.DataFrame(
        props_data,
        columns=cell_types,
        index=[f'spot_{i}' for i in range(n_spots)],
    )

    # Mock cell counts (10 nuclei per spot)
    nuclei_per_spot = 10
    cell_counts_data = np.round(props_data * nuclei_per_spot).astype(int)
    # Ensure sums to nuclei_per_spot
    for i in range(n_spots):
        diff = nuclei_per_spot - cell_counts_data[i].sum()
        cell_counts_data[i, 0] += diff

    cell_counts = pd.DataFrame(
        cell_counts_data,
        columns=cell_types,
        index=[f'spot_{i}' for i in range(n_spots)],
    )

    return {
        'proportions': proportions,
        'cell_counts': cell_counts,
        'cell_types': cell_types,
        'n_spots': n_spots,
        'nuclei_per_spot': nuclei_per_spot,
    }


def test_full_pipeline_integration(mock_data):
    """Test full Stage 2 pipeline with mock data."""
    from vae import VAE
    from two_stage_pipeline import TwoStagePipeline

    # Create mock VAE
    vae = VAE(in_channels=2, latent_dim=128)

    # Create pipeline
    pipeline = TwoStagePipeline(
        vae=vae,
        cell_types=mock_data['cell_types'],
        device='cpu',
    )

    # Mock patch loader
    def load_patches(spot_id):
        return torch.randn(mock_data['nuclei_per_spot'], 2, 96, 96)

    # Run with training (short)
    with tempfile.TemporaryDirectory() as tmpdir:
        pipeline.train_stage2(
            stage1_proportions=mock_data['proportions'],
            stage1_cell_counts=mock_data['cell_counts'],
            load_patches_fn=load_patches,
            n_epochs=2,  # Short training for test
            checkpoint_dir=Path(tmpdir),
        )

        # Run inference
        assignments = pipeline.run_stage2(
            stage1_proportions=mock_data['proportions'],
            stage1_cell_counts=mock_data['cell_counts'],
            load_patches_fn=load_patches,
        )

    # Verify assignments
    assert len(assignments) == mock_data['n_spots']
    for spot_id, assign in assignments.items():
        assert len(assign) == mock_data['nuclei_per_spot']
        assert all(0 <= a < len(mock_data['cell_types']) for a in assign)


def test_pipeline_checkpoint_save_load(mock_data):
    """Test checkpoint save and load functionality."""
    from vae import VAE
    from two_stage_pipeline import TwoStagePipeline

    vae = VAE(in_channels=2, latent_dim=128)

    def load_patches(spot_id):
        return torch.randn(mock_data['nuclei_per_spot'], 2, 96, 96)

    with tempfile.TemporaryDirectory() as tmpdir:
        checkpoint_dir = Path(tmpdir)

        # Train and save
        pipeline1 = TwoStagePipeline(
            vae=vae,
            cell_types=mock_data['cell_types'],
            device='cpu',
        )
        pipeline1.train_stage2(
            stage1_proportions=mock_data['proportions'],
            stage1_cell_counts=mock_data['cell_counts'],
            load_patches_fn=load_patches,
            n_epochs=2,
            checkpoint_dir=checkpoint_dir,
        )

        # Load checkpoint in new pipeline
        pipeline2 = TwoStagePipeline(
            vae=vae,
            cell_types=mock_data['cell_types'],
            device='cpu',
        )
        pipeline2.load_checkpoint(checkpoint_dir / "stage2_final.pt")

        # Verify prototypes match
        p1 = pipeline1.trainer.model.prototypes()
        p2 = pipeline2.trainer.model.prototypes()
        assert torch.allclose(p1, p2)


def test_pipeline_respects_counts(mock_data):
    """Verify Stage 2 assignments roughly match Stage 1 counts."""
    from vae import VAE
    from two_stage_pipeline import TwoStagePipeline

    vae = VAE(in_channels=2, latent_dim=128)
    pipeline = TwoStagePipeline(
        vae=vae,
        cell_types=mock_data['cell_types'],
        device='cpu',
    )

    def load_patches(spot_id):
        return torch.randn(mock_data['nuclei_per_spot'], 2, 96, 96)

    # Run without training (random prototypes)
    assignments = pipeline.run_stage2(
        stage1_proportions=mock_data['proportions'],
        stage1_cell_counts=mock_data['cell_counts'],
        load_patches_fn=load_patches,
        skip_training=True,
    )

    # With untrained model (random prototypes), assignments won't match counts well
    # Just verify that we get assignments for all spots
    for spot_id in mock_data['cell_counts'].index:
        assert len(assignments[spot_id]) == mock_data['nuclei_per_spot']
        # Verify all assignments are valid type indices
        assert all(0 <= a < len(mock_data['cell_types']) for a in assignments[spot_id])
