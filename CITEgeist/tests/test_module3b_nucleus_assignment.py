"""Tests for Module 3b: Per-Nucleus Assignment."""
import numpy as np
import pandas as pd
import pytest
from CITEgeist.model.module3b_nucleus_assignment import (
    NucleusAssignmentResult,
    run_nucleus_assignment,
    random_assign_nuclei_to_types,
)


def test_random_assign_nuclei_basic():
    """Test random assignment respects count constraints."""
    nucleus_ids = np.array([1, 2, 3, 4, 5])
    counts = np.array([3, 2])  # 3 of type 0, 2 of type 1

    rng = np.random.default_rng(42)
    assignments = random_assign_nuclei_to_types(nucleus_ids, counts, rng)

    assert len(assignments) == 5
    type_counts = [0, 0]
    for nid, type_idx in assignments.items():
        type_counts[type_idx] += 1

    assert type_counts == [3, 2]


def test_run_nucleus_assignment_random_default():
    """Test that random assignment is the default."""
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

    assert isinstance(result, NucleusAssignmentResult)
    assert result.method == "random"
    assert result.classifier is None
    assert result.morphology_features is None
    assert result.assignment_probs is None
    assert len(result.assignments) == 2


def test_run_nucleus_assignment_morphology():
    """Test end-to-end nucleus assignment with morphology guidance."""
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
        use_morphology=True,  # Enable morphology guidance
    )

    assert isinstance(result, NucleusAssignmentResult)
    assert result.method == "morphology"
    assert len(result.assignments) == 50  # all nuclei assigned
    assert result.classifier is not None
    assert result.morphology_features is not None
    assert len(result.morphology_features) == 50


def test_result_assignments_valid_types():
    """Test that all assignments are valid cell types (random mode)."""
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
        random_seed=123,  # For reproducibility
    )

    # All assigned types should be valid
    for nid, cell_type in result.assignments.items():
        assert cell_type in ['Cancer', 'Immune']


def test_run_nucleus_assignment_with_cell_features():
    """Test assignment using cell morphology features."""
    # Create synthetic nucleus and cell masks
    nucleus_mask = np.zeros((100, 100), dtype=np.int32)
    cell_mask = np.zeros((100, 100), dtype=np.int32)
    y, x = np.ogrid[:100, :100]

    # Two nuclei/cells
    nucleus_mask[((x - 30)**2 + (y - 30)**2) <= 25] = 1
    cell_mask[((x - 30)**2 + (y - 30)**2) <= 100] = 1
    nucleus_mask[((x - 70)**2 + (y - 70)**2) <= 25] = 2
    cell_mask[((x - 70)**2 + (y - 70)**2) <= 100] = 2

    # Both in same spot
    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2],
        'spot_id': ['spot_A', 'spot_A'],
    })

    proportions = pd.DataFrame({
        'spot_id': ['spot_A'],
        'TypeA': [0.5],
        'TypeB': [0.5],
    })

    nuclei_counts = pd.Series({'spot_A': 2})
    cell_types = ['TypeA', 'TypeB']

    result = run_nucleus_assignment(
        mask=nucleus_mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        cell_mask=cell_mask,
        use_morphology=True,  # Enable morphology to test cell features
    )

    # Should have assignments for both nuclei
    assert len(result.assignments) == 2
    assert result.method == "morphology"
    # Each should be assigned to one of the types
    assert set(result.assignments.values()) <= {'TypeA', 'TypeB'}
    # Morphology features should include cell features
    assert 'cell_area' in result.morphology_features.columns
    assert 'nc_ratio' in result.morphology_features.columns
