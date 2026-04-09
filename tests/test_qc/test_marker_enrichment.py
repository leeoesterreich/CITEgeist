"""Tests for marker_enrichment.py."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from CITEgeist.model.qc import QCResult
from CITEgeist.model.qc.marker_enrichment import (
    check_cross_patient_consistency,
    check_internal_coherence,
    compute_marker_enrichment,
    run_marker_enrichment,
)


def _make_enriched_adata():
    """Create AnnData where Macrophage cells have high CD68/CD163."""
    rng = np.random.default_rng(42)
    n_cells = 200
    n_genes = 10
    gene_names = ["CD68", "CD163", "CSF1R", "CD8A", "CD8B",
                  "EPCAM", "PECAM1", "COL1A1", "GENE_X", "GENE_Y"]

    types = (["Macrophages"] * 100) + (["Epithelial"] * 100)
    X = rng.random((n_cells, n_genes)) * 2

    # Boost Macrophage markers in Macrophage cells
    X[:100, 0] += 10  # CD68
    X[:100, 1] += 8   # CD163
    X[:100, 2] += 6   # CSF1R

    # Boost Epithelial markers in Epithelial cells
    X[100:, 5] += 10  # EPCAM

    return AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame({
            "cell_type": types,
            "spot_id": [f"s{i % 50}" for i in range(n_cells)],
        }, index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=gene_names),
    )


@pytest.mark.unit
class TestComputeEnrichment:
    """Tests for per-type marker enrichment."""

    def test_macrophage_markers_enriched(self):
        adata = _make_enriched_adata()
        result = compute_marker_enrichment(adata)
        mac = result[result["cell_type"] == "Macrophages"]
        cd68 = mac[mac["marker"] == "CD68"]
        assert len(cd68) == 1
        assert cd68.iloc[0]["log2fc"] > 1.0
        assert cd68.iloc[0]["qvalue"] < 0.05

    def test_returns_expected_columns(self):
        adata = _make_enriched_adata()
        result = compute_marker_enrichment(adata)
        for col in ["cell_type", "marker", "log2fc", "pvalue", "qvalue", "auc"]:
            assert col in result.columns


@pytest.mark.unit
class TestCrossPatientConsistency:
    """Tests for cross-patient consistency checking."""

    def test_consistent_markers_not_flagged(self):
        data = pd.DataFrame({
            "cell_type": ["Mac"] * 12,
            "marker": ["CD68"] * 12,
            "log2fc": [2.0] * 12,
            "qvalue": [0.01] * 12,
            "patient_id": [f"P{i}" for i in range(12)],
        })
        flags = check_cross_patient_consistency(data)
        assert len(flags) == 0

    def test_inconsistent_markers_flagged(self):
        log2fcs = [2.0] * 5 + [-1.0] * 7
        qvals = [0.01] * 5 + [0.5] * 7
        data = pd.DataFrame({
            "cell_type": ["Mac"] * 12,
            "marker": ["CD68"] * 12,
            "log2fc": log2fcs,
            "qvalue": qvals,
            "patient_id": [f"P{i}" for i in range(12)],
        })
        flags = check_cross_patient_consistency(data)
        assert len(flags) > 0


@pytest.mark.unit
class TestInternalCoherence:
    """Tests for proportion ↔ GEX coherence."""

    def test_coherent_data_high_concordance(self):
        proportions = pd.DataFrame({
            "Macrophages": [0.8, 0.1, 0.9, 0.05],
            "Epithelial": [0.2, 0.9, 0.1, 0.95],
        }, index=["s0", "s1", "s2", "s3"])

        gex_layers = {
            "Macrophages": pd.DataFrame(
                {
                    "CD68": [10.0, 1.0, 12.0, 0.5],
                    "CD163": [8.0, 0.5, 10.0, 0.3],
                },
                index=["s0", "s1", "s2", "s3"],
            ),
            "Epithelial": pd.DataFrame(
                {
                    "EPCAM": [1.0, 10.0, 0.5, 12.0],
                    "KRT18": [0.5, 9.0, 0.3, 11.0],
                },
                index=["s0", "s1", "s2", "s3"],
            ),
        }

        result = check_internal_coherence(proportions, gex_layers)
        # Macrophages is dominant (>0.3) in s0, s2
        # Epithelial is dominant (>0.3) in s1, s3
        assert "Macrophages" in result
        assert "Epithelial" in result
        # Both should have high concordance if marker expression matches dominance
        assert result["Macrophages"]["n_dominant_spots"] == 2
        assert result["Epithelial"]["n_dominant_spots"] == 2

    def test_missing_markers_graceful(self):
        proportions = pd.DataFrame({
            "Macrophages": [0.8, 0.1],
            "Epithelial": [0.2, 0.9],
        }, index=["s0", "s1"])

        # GEX layers with no canonical markers
        gex_layers = {
            "Macrophages": pd.DataFrame(
                {"GENE_X": [5.0, 2.0]}, index=["s0", "s1"]
            ),
            "Epithelial": pd.DataFrame(
                {"GENE_Y": [1.0, 8.0]}, index=["s0", "s1"]
            ),
        }

        result = check_internal_coherence(proportions, gex_layers)
        # Should return empty dict since no canonical markers available
        assert len(result) == 0


@pytest.mark.unit
class TestRunMarkerEnrichment:
    """Tests for the main entry point."""

    def test_returns_qcresult(self):
        adata = _make_enriched_adata()
        rng = np.random.default_rng(42)
        proportions = pd.DataFrame(
            rng.dirichlet([1, 1], size=50),
            columns=["Macrophages", "Epithelial"],
            index=[f"s{i}" for i in range(50)],
        )
        result = run_marker_enrichment(adata, proportions)
        assert isinstance(result, QCResult)
        assert "enrichment_heatmap" in result.figures
