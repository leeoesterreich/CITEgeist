"""Tests for gt_rebuild/rna_limits.py (synthetic)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import rna_limits  # noqa: E402

pytestmark = pytest.mark.unit


def test_coverage_lists_absent_genes():
    out = rna_limits.panel_coverage(["CD8B", "ICOS", "OTHER"])
    assert out["present"] == ["CD8B", "ICOS"]
    assert out["absent"] == ["TCF7", "LEF1", "CD40LG", "GZMH"]


def test_detection_fraction_hand_counted():
    counts = pd.DataFrame({"CD4": [1, 0, 2, 5, 9]}, index=list("abcde"))
    assert rna_limits.detection_fraction(counts, list("abcd"), "CD4") == pytest.approx(0.75)


def test_ambiguous_share_is_cell_weighted():
    tsub = pd.DataFrame({"subcluster": ["0", "1"], "n_cells": [10, 30], "ambiguous": [True, False]})
    assert rna_limits.ambiguous_share(tsub) == pytest.approx(0.25)


def test_rna_limits_end_to_end():
    ids = [f"c{i}" for i in range(8)]
    counts = pd.DataFrame({"CD4": [1, 0, 1, 1, 0, 0, 1, 0]}, index=ids)
    labels7 = pd.Series(["CD4+ T cells"] * 4 + ["Macrophages"] * 4, index=ids)
    level1 = pd.Series(["0"] * 4 + ["1"] * 2 + ["2"] * 2, index=ids)
    names = pd.DataFrame({"cluster": ["0", "1", "2"], "label": ["T cells", "Macrophages", "Macrophages"]})
    tsub = pd.DataFrame({"subcluster": ["0"], "n_cells": [4], "ambiguous": [False]})
    out = rna_limits.rna_limits(counts, labels7, level1, names, tsub, cd4_gene="CD4")
    assert out["cd4_detected_frac_tcells"] == pytest.approx(0.75)
    assert out["cd4_detected_frac_by_macrophage_cluster"] == {"1": 0.0, "2": 0.5}
    assert out["ambiguous_tcell_share"] == 0.0
    assert out["subcluster_table"] == tsub.to_dict("records")


def test_dotplot_writes_file(tmp_path):
    ad = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    rng = np.random.RandomState(0)
    genes = ["CD8A", "CD8B", "CD4", "IL7R"]
    a = ad.AnnData(rng.poisson(1.0, (20, 4)).astype(np.float32))
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(20)]
    tsub = pd.Series(["0", "1"] * 10, index=a.obs_names)
    p = tmp_path / "dot.png"
    rna_limits.dotplot(a, tsub, p)
    assert p.exists() and p.stat().st_size > 0
