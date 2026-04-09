"""Tests for Module 3b: Per-Nucleus Assignment."""
import numpy as np
import pandas as pd
import pytest
from CITEgeist.model.assignment.module3b_nucleus_assignment import (
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


def test_single_cell_adata_with_none_morphology():
    """Test create_single_cell_adata handles morphology_features=None."""
    from CITEgeist.model.assignment.single_cell_output import create_single_cell_adata

    cell_gex = pd.DataFrame(
        np.random.rand(3, 5),
        index=[1, 2, 3],
        columns=[f'gene_{i}' for i in range(5)],
    )
    assignments = {1: 'TypeA', 2: 'TypeB', 3: 'TypeA'}

    adata = create_single_cell_adata(
        cell_gex=cell_gex,
        morphology_features=None,
        assignments=assignments,
        sample_name='test_sample',
    )

    assert adata.n_obs == 3
    assert adata.n_vars == 5
    assert 'cell_type' in adata.obs.columns
    assert 'spatial' not in adata.obsm  # No coords without morphology
    assert adata.uns['source_sample'] == 'test_sample'


def test_single_cell_adata_with_morphology():
    """Test create_single_cell_adata stores morphology in obsm."""
    from CITEgeist.model.assignment.single_cell_output import create_single_cell_adata

    cell_gex = pd.DataFrame(
        np.random.rand(3, 5),
        index=[1, 2, 3],
        columns=[f'gene_{i}' for i in range(5)],
    )
    assignments = {1: 'TypeA', 2: 'TypeB', 3: 'TypeA'}
    morph = pd.DataFrame({
        'nucleus_id': [1, 2, 3],
        'spot_id': ['s0', 's0', 's1'],
        'centroid_x': [10.0, 20.0, 30.0],
        'centroid_y': [15.0, 25.0, 35.0],
        'area': [100, 150, 120],
        'circularity': [0.8, 0.9, 0.85],
        'eccentricity': [0.3, 0.2, 0.4],
        'solidity': [0.95, 0.9, 0.92],
        'major_axis_length': [12.0, 14.0, 13.0],
        'minor_axis_length': [10.0, 12.0, 11.0],
        'aspect_ratio': [1.2, 1.17, 1.18],
    })

    adata = create_single_cell_adata(
        cell_gex=cell_gex,
        morphology_features=morph,
        assignments=assignments,
        sample_name='test_sample',
    )

    assert adata.n_obs == 3
    assert 'spatial' in adata.obsm
    assert 'morphology' in adata.obsm
    assert adata.obsm['morphology'].shape[1] == 7  # 7 mask features present
    assert 'morphology_feature_names' in adata.uns


