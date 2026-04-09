"""Tests for report.py orchestrator."""

import json
import os

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from CITEgeist.model.qc import QCResult
from CITEgeist.model.qc.report import run_qc


@pytest.mark.unit
class TestRunQCSelfConsistency:
    """Tests for self_consistency mode."""

    def test_returns_dict_of_qcresults(self, mock_sc_adata, mock_proportions):
        results = run_qc(mock_sc_adata, mock_proportions, mode="self_consistency")
        assert isinstance(results, dict)
        assert "single_cell" in results
        assert "marker_enrichment" in results
        assert isinstance(results["single_cell"], QCResult)

    def test_skips_proportion_and_gex_qc(self, mock_sc_adata, mock_proportions):
        results = run_qc(mock_sc_adata, mock_proportions, mode="self_consistency")
        assert "proportion" not in results
        assert "gex" not in results


@pytest.mark.unit
class TestRunQCBenchmark:
    """Tests for benchmark mode."""

    def test_requires_gt_proportions(self, mock_sc_adata, mock_proportions):
        with pytest.raises(ValueError, match="gt_proportions"):
            run_qc(mock_sc_adata, mock_proportions, mode="benchmark")

    def test_includes_proportion_and_gex(
        self, mock_sc_adata, mock_proportions, mock_gt_proportions, mock_gt_gex_layers
    ):
        results = run_qc(
            mock_sc_adata, mock_proportions, mode="benchmark",
            gt_proportions=mock_gt_proportions,
            gt_gex_layers=mock_gt_gex_layers,
        )
        assert "proportion" in results
        assert "gex" in results


@pytest.mark.unit
class TestSaveFigures:
    """Tests for figure saving."""

    def test_saves_pdfs(
        self, mock_sc_adata, mock_proportions, mock_gt_proportions, tmp_path,
    ):
        run_qc(
            mock_sc_adata, mock_proportions, mode="benchmark",
            gt_proportions=mock_gt_proportions,
            output_dir=str(tmp_path),
        )
        assert os.path.exists(os.path.join(tmp_path, "qc_summary.json"))
        pdfs = [f for f in os.listdir(str(tmp_path)) if f.endswith(".pdf")]
        assert len(pdfs) > 0, "Expected at least one PDF figure to be saved"

    def test_summary_json_readable(
        self, mock_sc_adata, mock_proportions, mock_gt_proportions, tmp_path,
    ):
        run_qc(
            mock_sc_adata, mock_proportions, mode="benchmark",
            gt_proportions=mock_gt_proportions,
            output_dir=str(tmp_path),
        )
        with open(os.path.join(tmp_path, "qc_summary.json")) as f:
            data = json.load(f)
        assert "flags" in data
        assert "mode" in data
        assert "n_cells_total" in data
        assert "n_cells_empty" in data
        assert "n_cells_analyzed" in data
