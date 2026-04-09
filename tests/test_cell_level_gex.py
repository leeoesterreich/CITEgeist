"""Tests for cell-level GEX distribution."""
import numpy as np
import pandas as pd
import pytest
from CITEgeist.model.gex.cell_level_gex import distribute_gex_to_cells


def test_distribute_gex_equal_split():
    """Test equal GEX distribution among cells of same type."""
    # 2 spots, 3 genes
    # Spot A: 2 macrophages
    # Spot B: 1 macrophage, 1 fibroblast
    deconvolved_gex = pd.DataFrame({
        'gene1': [100.0, 50.0, 30.0],
        'gene2': [200.0, 100.0, 60.0],
        'gene3': [300.0, 150.0, 90.0],
    }, index=['spotA:::Macrophage', 'spotB:::Macrophage', 'spotB:::Fibroblast'])

    assignments = {
        1: 'Macrophage',  # spot A
        2: 'Macrophage',  # spot A
        3: 'Macrophage',  # spot B
        4: 'Fibroblast',  # spot B
    }

    nucleus_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2, 3, 4],
        'spot_id': ['spotA', 'spotA', 'spotB', 'spotB'],
    })

    cell_gex = distribute_gex_to_cells(deconvolved_gex, assignments, nucleus_spot_map)

    # Check shape
    assert cell_gex.shape == (4, 3)

    # Spot A: 2 macs split 100 -> 50 each
    assert cell_gex.loc[1, 'gene1'] == 50.0
    assert cell_gex.loc[2, 'gene1'] == 50.0

    # Spot B: 1 mac gets full 50
    assert cell_gex.loc[3, 'gene1'] == 50.0

    # Spot B: 1 fib gets full 30
    assert cell_gex.loc[4, 'gene1'] == 30.0


def test_distribute_gex_preserves_total():
    """Test that total expression is preserved after distribution."""
    np.random.seed(42)
    deconvolved_gex = pd.DataFrame(
        np.random.rand(3, 10) * 100,
        index=['spot1:::TypeA', 'spot1:::TypeB', 'spot2:::TypeA'],
        columns=[f'gene{i}' for i in range(10)]
    )

    assignments = {1: 'TypeA', 2: 'TypeB', 3: 'TypeA', 4: 'TypeA'}
    nucleus_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2, 3, 4],
        'spot_id': ['spot1', 'spot1', 'spot2', 'spot2'],
    })

    cell_gex = distribute_gex_to_cells(deconvolved_gex, assignments, nucleus_spot_map)

    # Total per spot-type should be preserved
    # spot1:::TypeA -> nucleus 1 only (nucleus 2 is TypeB)
    spot1_typeA_original = deconvolved_gex.loc['spot1:::TypeA'].sum()
    spot1_typeA_cells = cell_gex.loc[1].sum()
    assert np.isclose(spot1_typeA_cells, spot1_typeA_original)


def test_allocate_gex_type_reference_preserves_spot_totals():
    """Test that type-specific reference allocation preserves spot totals."""
    from CITEgeist.model.cell_level_gex import allocate_gex_type_reference
    spot_gex = pd.DataFrame(
        {'gene_A': [10.0, 20.0], 'gene_B': [5.0, 15.0]},
        index=['spot_0', 'spot_1']
    )
    proportions = pd.DataFrame(
        {'TypeX': [0.6, 0.3], 'TypeY': [0.4, 0.7]},
        index=['spot_0', 'spot_1']
    )
    hard_labels = np.array(['TypeX', 'TypeX', 'TypeY', 'TypeY', 'TypeX'])
    scores = np.array([[0.8, 0.2], [0.7, 0.3], [0.3, 0.7], [0.2, 0.8], [0.6, 0.4]])
    type_names = ['TypeX', 'TypeY']
    barcodes = np.array(['spot_0', 'spot_0', 'spot_0', 'spot_1', 'spot_1'])
    nucleus_ids = np.arange(5)

    cell_gex = allocate_gex_type_reference(
        hard_labels=hard_labels, scores=scores, type_names=type_names,
        barcodes=barcodes, nucleus_ids=nucleus_ids,
        proportions=proportions, spot_gex=spot_gex,
    )
    assert cell_gex.shape == (5, 2)
    for spot in ['spot_0', 'spot_1']:
        mask = barcodes == spot
        spot_sum = cell_gex.values[mask].sum(axis=0)
        expected = spot_gex.loc[spot].values
        np.testing.assert_allclose(spot_sum, expected, rtol=0.01)


def test_allocate_gex_handles_missing_spots():
    """Test that cells in missing spots get zero GEX."""
    from CITEgeist.model.cell_level_gex import allocate_gex_type_reference
    spot_gex = pd.DataFrame({'g1': [10.0]}, index=['spot_0'])
    proportions = pd.DataFrame({'T': [1.0]}, index=['spot_0'])
    hard_labels = np.array(['T', 'T'])
    scores = np.array([[1.0], [1.0]])
    barcodes = np.array(['spot_0', 'spot_MISSING'])
    nucleus_ids = np.arange(2)
    cell_gex = allocate_gex_type_reference(
        hard_labels=hard_labels, scores=scores, type_names=['T'],
        barcodes=barcodes, nucleus_ids=nucleus_ids,
        proportions=proportions, spot_gex=spot_gex,
    )
    assert cell_gex.shape == (2, 1)
    assert cell_gex.iloc[0, 0] == 10.0  # spot_0 cell gets all GEX
    assert cell_gex.iloc[1, 0] == 0.0   # missing spot gets 0
