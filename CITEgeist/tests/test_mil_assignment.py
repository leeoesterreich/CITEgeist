"""Tests for MIL-based nucleus assignment in Module 3b."""
import numpy as np
import pandas as pd
import pytest
import torch

from CITEgeist.model.module3b_nucleus_assignment import (
    run_nucleus_assignment_mil,
    NucleusAssignmentResult,
)


@pytest.fixture
def fake_spot_data():
    """Create minimal fake data for 2 spots, 3 types, 10 nuclei."""
    n_types = 3
    cell_types = ["TypeA", "TypeB", "TypeC"]

    # 2 spots: spot_0 has 4 nuclei, spot_1 has 6
    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': range(10),
        'spot_id': ['spot_0'] * 4 + ['spot_1'] * 6,
    })

    proportions = pd.DataFrame({
        'spot_id': ['spot_0', 'spot_1'],
        'TypeA': [0.5, 0.33],
        'TypeB': [0.25, 0.33],
        'TypeC': [0.25, 0.34],
    })

    nuclei_counts = pd.Series({'spot_0': 4, 'spot_1': 6})

    # Pre-computed embeddings: dict of spot_id -> (N, 384) array
    embeddings = {
        'spot_0': np.random.randn(4, 384).astype(np.float32),
        'spot_1': np.random.randn(6, 384).astype(np.float32),
    }

    return nuclei_spot_map, proportions, nuclei_counts, cell_types, embeddings


def test_mil_assignment_returns_result(fake_spot_data):
    """MIL assignment should return NucleusAssignmentResult."""
    nuclei_spot_map, proportions, nuclei_counts, cell_types, embeddings = fake_spot_data

    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        n_epochs=5,  # Short training for test
        lambda_prior=1.0,
        device="cpu",
    )

    assert isinstance(result, NucleusAssignmentResult)
    assert result.method == "mil"
    assert len(result.assignments) == 10
    # All assigned types should be valid
    for nid, ctype in result.assignments.items():
        assert ctype in cell_types


def test_mil_assignment_respects_counts(fake_spot_data):
    """Assignments per spot should respect discretized count constraints."""
    nuclei_spot_map, proportions, nuclei_counts, cell_types, embeddings = fake_spot_data

    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        n_epochs=5,
        lambda_prior=0.0,
        device="cpu",
    )

    # Check count constraints per spot
    for spot_id in ['spot_0', 'spot_1']:
        spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
        n_nuclei = int(nuclei_counts[spot_id])
        assigned_count = sum(
            1 for nid in spot_nuclei['nucleus_id']
            if nid in result.assignments
        )
        assert assigned_count == len(spot_nuclei)
