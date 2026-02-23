"""Tests for single-cell AnnData output."""
import numpy as np
import pandas as pd
import pytest
import anndata as ad
from CITEgeist.model.single_cell_output import create_single_cell_adata


def test_create_adata_basic():
    """Test basic AnnData creation."""
    cell_gex = pd.DataFrame(
        np.random.rand(10, 5),
        index=range(10),
        columns=[f'gene{i}' for i in range(5)]
    )

    morphology = pd.DataFrame({
        'nucleus_id': range(10),
        'spot_id': [f'spot{i//2}' for i in range(10)],
        'centroid_x': np.random.rand(10) * 100,
        'centroid_y': np.random.rand(10) * 100,
        'area': np.random.rand(10) * 500,
        'circularity': np.random.rand(10),
    })

    assignments = {i: f'type{i % 3}' for i in range(10)}

    adata = create_single_cell_adata(
        cell_gex=cell_gex,
        morphology_features=morphology,
        assignments=assignments,
        sample_name='test_sample',
    )

    assert isinstance(adata, ad.AnnData)
    assert adata.n_obs == 10
    assert adata.n_vars == 5
    assert 'cell_type' in adata.obs.columns
    assert 'spot_id' in adata.obs.columns
    assert 'x' in adata.obs.columns
    assert 'y' in adata.obs.columns


def test_adata_spatial_coords():
    """Test that spatial coordinates are properly stored."""
    cell_gex = pd.DataFrame(
        np.ones((3, 2)),
        index=[1, 2, 3],
        columns=['g1', 'g2']
    )

    morphology = pd.DataFrame({
        'nucleus_id': [1, 2, 3],
        'spot_id': ['s1', 's1', 's2'],
        'centroid_x': [10.0, 20.0, 30.0],
        'centroid_y': [15.0, 25.0, 35.0],
        'area': [100, 200, 300],
    })

    assignments = {1: 'A', 2: 'B', 3: 'A'}

    adata = create_single_cell_adata(cell_gex, morphology, assignments, 'test')

    assert list(adata.obs['x']) == [10.0, 20.0, 30.0]
    assert list(adata.obs['y']) == [15.0, 25.0, 35.0]
    assert 'spatial' in adata.obsm
    assert adata.obsm['spatial'].shape == (3, 2)
