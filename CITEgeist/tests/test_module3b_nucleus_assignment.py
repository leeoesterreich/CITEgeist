"""Tests for Module 3b: Per-Nucleus Assignment."""
import numpy as np
import pandas as pd
import pytest
from CITEgeist.model.module3b_nucleus_assignment import (
    NucleusAssignmentResult,
    run_nucleus_assignment,
)


def test_run_nucleus_assignment_basic():
    """Test end-to-end nucleus assignment."""
    np.random.seed(42)

    # Create mock mask with 50 nuclei (5 per spot, 10 spots)
    mask = np.zeros((100, 100), dtype=np.int32)
    nucleus_id = 1
    for i in range(10):
        for j in range(5):
            cx, cy = 10 + i*9, 10 + j*18
            y, x = np.ogrid[:100, :100]
            circle = ((x - cx)**2 + (y - cy)**2) <= 16
            mask[circle] = nucleus_id
            nucleus_id += 1

    # Mock nuclei-to-spot mapping
    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': range(1, 51),
        'spot_id': [f'spot_{i // 5}' for i in range(50)],
    })

    # Mock proportions (10 spots x 3 types)
    proportions = pd.DataFrame({
        'spot_id': [f'spot_{i}' for i in range(10)],
        'type_0': [0.5] * 10,
        'type_1': [0.3] * 10,
        'type_2': [0.2] * 10,
    })

    # Mock nuclei counts
    nuclei_counts = pd.Series(
        [5] * 10,
        index=[f'spot_{i}' for i in range(10)]
    )

    cell_types = ['type_0', 'type_1', 'type_2']

    result = run_nucleus_assignment(
        mask=mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
    )

    assert isinstance(result, NucleusAssignmentResult)
    assert len(result.assignments) == 50  # all nuclei assigned
    assert result.classifier is not None
    assert result.morphology_features is not None
    assert len(result.morphology_features) == 50


def test_result_assignments_valid_types():
    """Test that all assignments are valid cell types."""
    np.random.seed(123)

    # Minimal setup
    mask = np.zeros((50, 50), dtype=np.int32)
    mask[10:15, 10:15] = 1
    mask[10:15, 30:35] = 2

    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2],
        'spot_id': ['spot_0', 'spot_0'],
    })

    proportions = pd.DataFrame({
        'spot_id': ['spot_0'],
        'Cancer': [0.5],
        'Immune': [0.5],
    })

    nuclei_counts = pd.Series([2], index=['spot_0'])

    result = run_nucleus_assignment(
        mask=mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=['Cancer', 'Immune'],
    )

    # All assigned types should be valid
    for nid, cell_type in result.assignments.items():
        assert cell_type in ['Cancer', 'Immune']
