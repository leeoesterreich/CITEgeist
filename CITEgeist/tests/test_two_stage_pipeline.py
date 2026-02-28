# CITEgeist/tests/test_two_stage_pipeline.py
"""Tests for two-stage pipeline orchestration."""
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


def test_pipeline_end_to_end():
    """Pipeline runs Stage 1 + Stage 2 and produces assignments."""
    from two_stage_pipeline import TwoStagePipeline
    from vae import VAE

    # Create mock VAE
    vae = VAE(in_channels=2, latent_dim=128)
    cell_types = ['B cells', 'T cells', 'Macrophages']

    # Mock Stage 1 results (from hybrid method)
    stage1_proportions = pd.DataFrame({
        'B cells': [0.4, 0.2, 0.3],
        'T cells': [0.3, 0.5, 0.2],
        'Macrophages': [0.3, 0.3, 0.5],
    }, index=['spot_0', 'spot_1', 'spot_2'])

    stage1_cell_counts = pd.DataFrame({
        'B cells': [4, 2, 3],
        'T cells': [3, 5, 2],
        'Macrophages': [3, 3, 5],
    }, index=['spot_0', 'spot_1', 'spot_2'])

    nuclei_counts = pd.Series([10, 10, 10], index=['spot_0', 'spot_1', 'spot_2'])

    pipeline = TwoStagePipeline(
        vae=vae,
        cell_types=cell_types,
        device='cpu',
    )

    # Mock patch loader
    def load_patches(spot_id):
        n = int(nuclei_counts[spot_id])
        return torch.randn(n, 2, 96, 96)

    # Run Stage 2 assignment
    assignments = pipeline.run_stage2(
        stage1_proportions=stage1_proportions,
        stage1_cell_counts=stage1_cell_counts,
        load_patches_fn=load_patches,
        skip_training=True,  # Skip training for quick test
    )

    assert 'spot_0' in assignments
    assert len(assignments['spot_0']) == 10
    assert all(0 <= a < 3 for a in assignments['spot_0'])
