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


def test_constrained_hungarian_dapi(tmp_path):
    """Test constrained Hungarian assignment with synthetic DAPI patches."""
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_constrained

    np.random.seed(42)
    cell_types = ['TypeA', 'TypeB', 'TypeC']
    n_spots = 20
    n_nuclei_per_spot = 5

    # Create synthetic patches (2-channel, 48x48) per spot
    nuclei_records = []
    for i in range(n_spots):
        spot_id = f"spot_{i}"
        patches = np.random.rand(n_nuclei_per_spot, 2, 48, 48).astype(np.float32)
        # Make TypeA-dominant spots brighter in DAPI
        if i < 10:
            patches[:, 0, :, :] += 2.0  # brighter DAPI
        np.save(tmp_path / f"{spot_id}_patches.npy", patches)
        nuc_ids = np.arange(i * n_nuclei_per_spot, (i + 1) * n_nuclei_per_spot)
        np.save(tmp_path / f"{spot_id}_nucleus_ids.npy", nuc_ids)

        for nid in nuc_ids:
            nuclei_records.append({'nucleus_id': nid, 'spot_id': spot_id})

    nuclei_spot_map = pd.DataFrame(nuclei_records)

    # Proportions: first 10 spots are TypeA-dominant, rest are TypeB-dominant
    prop_data = []
    for i in range(n_spots):
        if i < 10:
            prop_data.append({'spot_id': f"spot_{i}", 'TypeA': 0.6, 'TypeB': 0.2, 'TypeC': 0.2})
        else:
            prop_data.append({'spot_id': f"spot_{i}", 'TypeA': 0.2, 'TypeB': 0.6, 'TypeC': 0.2})
    proportions = pd.DataFrame(prop_data)

    nuclei_counts = pd.Series(
        [n_nuclei_per_spot] * n_spots,
        index=[f"spot_{i}" for i in range(n_spots)]
    )

    result = run_nucleus_assignment_constrained(
        patches_dir=str(tmp_path),
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        purity_threshold=0.5,
    )

    assert isinstance(result, NucleusAssignmentResult)
    assert result.method == "constrained_hungarian"
    assert len(result.assignments) > 0
    # All assignments should be valid cell types
    for nid, ctype in result.assignments.items():
        assert ctype in cell_types


def test_constrained_hungarian_no_patches(tmp_path):
    """Test constrained Hungarian with empty patches dir raises ValueError."""
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_constrained

    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2],
        'spot_id': ['spot_0', 'spot_0'],
    })
    proportions = pd.DataFrame({
        'spot_id': ['spot_0'],
        'TypeA': [0.5],
        'TypeB': [0.5],
    })
    nuclei_counts = pd.Series([2], index=['spot_0'])

    with pytest.raises(ValueError, match="No features collected"):
        run_nucleus_assignment_constrained(
            patches_dir=str(tmp_path),
            nuclei_spot_map=nuclei_spot_map,
            proportions=proportions,
            nuclei_counts=nuclei_counts,
            cell_types=['TypeA', 'TypeB'],
        )


def test_single_cell_adata_with_none_morphology():
    """Test create_single_cell_adata handles morphology_features=None."""
    from CITEgeist.model.single_cell_output import create_single_cell_adata

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
    from CITEgeist.model.single_cell_output import create_single_cell_adata

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


def test_vae_assignment_runs():
    """VAE assignment requires trained models - placeholder test."""
    pytest.skip("Integration test - requires trained models")


def test_vae_assignment_import():
    """Test that VAE assignment function can be imported."""
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_vae
    assert callable(run_nucleus_assignment_vae)


def test_vae_assignment_missing_checkpoint():
    """Test that missing checkpoint raises FileNotFoundError."""
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_vae

    # Create minimal inputs
    image = np.random.rand(3, 100, 100).astype(np.float32)
    nuclei_df = pd.DataFrame({
        'nucleus_id': [1, 2],
        'spot_id': ['spot_0', 'spot_0'],
        'bbox_x_min': [10, 50],
        'bbox_y_min': [10, 50],
        'bbox_x_max': [30, 70],
        'bbox_y_max': [30, 70],
    })
    proportions = pd.DataFrame({
        'spot_id': ['spot_0'],
        'TypeA': [0.5],
        'TypeB': [0.5],
    })

    with pytest.raises(FileNotFoundError):
        run_nucleus_assignment_vae(
            image=image,
            nuclei_df=nuclei_df,
            proportions=proportions,
            cell_types=['TypeA', 'TypeB'],
            vae_checkpoint='nonexistent_vae.pt',
            prototype_checkpoint='nonexistent_proto.pt',
            device='cpu',
        )


def test_vae_assignment_missing_columns():
    """Test that missing required columns raises ValueError."""
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_vae

    image = np.random.rand(3, 100, 100).astype(np.float32)
    # Missing bbox columns
    nuclei_df = pd.DataFrame({
        'nucleus_id': [1, 2],
        'spot_id': ['spot_0', 'spot_0'],
    })
    proportions = pd.DataFrame({
        'spot_id': ['spot_0'],
        'TypeA': [0.5],
        'TypeB': [0.5],
    })

    with pytest.raises(ValueError, match="missing required columns"):
        run_nucleus_assignment_vae(
            image=image,
            nuclei_df=nuclei_df,
            proportions=proportions,
            cell_types=['TypeA', 'TypeB'],
            vae_checkpoint='vae.pt',
            prototype_checkpoint='proto.pt',
            device='cpu',
        )
