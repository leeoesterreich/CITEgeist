"""Tests for single_cell_qc.py."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from CITEgeist.model.qc import QCResult
from CITEgeist.model.qc.single_cell_qc import (
    check_compartment_emptiness,
    compute_cell_metrics,
    detect_empty_cells,
    run_single_cell_qc,
)


@pytest.mark.unit
class TestComputeCellMetrics:
    """Tests for per-cell metric computation."""

    def test_returns_dataframe_with_expected_columns(self, mock_sc_adata):
        df = compute_cell_metrics(mock_sc_adata)
        assert isinstance(df, pd.DataFrame)
        for col in ["total_umi", "n_genes", "pct_mt", "pct_ribo"]:
            assert col in df.columns

    def test_umi_counts_nonnegative(self, mock_sc_adata):
        df = compute_cell_metrics(mock_sc_adata)
        assert (df["total_umi"] >= 0).all()

    def test_pct_mt_between_0_and_100(self, mock_sc_adata):
        df = compute_cell_metrics(mock_sc_adata)
        assert (df["pct_mt"] >= 0).all()
        assert (df["pct_mt"] <= 100).all()

    def test_n_genes_leq_total_genes(self, mock_sc_adata):
        df = compute_cell_metrics(mock_sc_adata)
        assert (df["n_genes"] <= mock_sc_adata.n_vars).all()


@pytest.mark.unit
class TestDetectEmptyCells:
    """Tests for empty cell detection."""

    def test_returns_boolean_mask(self, mock_sc_adata):
        metrics = compute_cell_metrics(mock_sc_adata)
        mask = detect_empty_cells(metrics)
        assert mask.dtype == bool
        assert len(mask) == len(mock_sc_adata)

    def test_high_threshold_flags_more_cells(self, mock_sc_adata):
        metrics = compute_cell_metrics(mock_sc_adata)
        mask_low = detect_empty_cells(metrics, umi_threshold=10, genes_threshold=5)
        mask_high = detect_empty_cells(metrics, umi_threshold=500, genes_threshold=100)
        assert mask_high.sum() >= mask_low.sum()

    def test_zero_threshold_flags_nothing(self, mock_sc_adata):
        metrics = compute_cell_metrics(mock_sc_adata)
        mask = detect_empty_cells(metrics, umi_threshold=0, genes_threshold=0)
        assert mask.sum() == 0


@pytest.mark.unit
class TestCompartmentCheck:
    """Tests for compartment-level empty cell flagging."""

    def test_flags_majority_empty_compartments(self):
        """If >50% of a type's cells are empty, flag it."""
        cell_types = pd.Series(["A", "A", "A", "A", "B", "B"])
        is_empty = np.array([True, True, True, False, False, False])
        flags = check_compartment_emptiness(cell_types, is_empty)
        assert any("A" in f for f in flags)
        assert not any("B" in f for f in flags)


@pytest.mark.unit
class TestRunSingleCellQC:
    """Tests for the main entry point."""

    def test_returns_qcresult(self, mock_sc_adata):
        result = run_single_cell_qc(mock_sc_adata)
        assert isinstance(result, QCResult)

    def test_metrics_has_summary_table(self, mock_sc_adata):
        result = run_single_cell_qc(mock_sc_adata)
        assert "summary_table" in result.metrics
        assert isinstance(result.metrics["summary_table"], pd.DataFrame)

    def test_metrics_has_is_empty_mask(self, mock_sc_adata):
        result = run_single_cell_qc(mock_sc_adata)
        assert "is_empty" in result.metrics
        assert len(result.metrics["is_empty"]) == len(mock_sc_adata)

    def test_produces_figures(self, mock_sc_adata):
        result = run_single_cell_qc(mock_sc_adata)
        assert "violin_qc" in result.figures
        assert "empty_heatmap" in result.figures
