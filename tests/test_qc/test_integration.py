"""Integration test: full QC pipeline end-to-end."""

import json
import os

import numpy as np
import pandas as pd
import pytest

from CITEgeist.model.qc import QCResult
from CITEgeist.model.qc.report import run_qc


@pytest.mark.integration
class TestFullBenchmarkQC:
    """End-to-end benchmark mode QC."""

    def test_benchmark_produces_all_outputs(
        self, mock_sc_adata, mock_proportions, mock_gt_proportions,
        mock_gt_gex_layers, tmp_path,
    ):
        results = run_qc(
            mock_sc_adata, mock_proportions,
            mode="benchmark",
            gt_proportions=mock_gt_proportions,
            gt_gex_layers=mock_gt_gex_layers,
            output_dir=str(tmp_path),
        )

        # All modules ran
        assert "single_cell" in results
        assert "proportion" in results
        assert "marker_enrichment" in results

        # Summary JSON written
        summary_path = os.path.join(tmp_path, "qc_summary.json")
        assert os.path.exists(summary_path)
        with open(summary_path) as f:
            summary = json.load(f)
        assert summary["mode"] == "benchmark"
        assert summary["n_cells_total"] == 500

        # At least some PDF figures saved
        pdfs = [f for f in os.listdir(tmp_path) if f.endswith(".pdf")]
        assert len(pdfs) > 0


@pytest.mark.integration
class TestFullSelfConsistencyQC:
    """End-to-end self-consistency mode QC."""

    def test_self_consistency_produces_outputs(
        self, mock_sc_adata, mock_proportions, tmp_path,
    ):
        results = run_qc(
            mock_sc_adata, mock_proportions,
            mode="self_consistency",
            output_dir=str(tmp_path),
        )

        assert "single_cell" in results
        assert "marker_enrichment" in results
        assert "proportion" not in results

        summary_path = os.path.join(tmp_path, "qc_summary.json")
        assert os.path.exists(summary_path)
