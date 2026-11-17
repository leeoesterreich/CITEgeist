"""Tests for Benchmarking/xenium_benchmarking/gt_rebuild/tcell_split.py (synthetic data only)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy")
pytest.importorskip("leidenalg")
pytest.importorskip("igraph")

import anndata as ad  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import gt_config  # noqa: E402
from gt_rebuild import tcell_split as ts  # noqa: E402

pytestmark = pytest.mark.unit

CD8 = ["A", "B"]
CD4 = ["C", "D"]
SUB = pd.Series(["s1", "s1", "s2", "s2"], index=list("wxyz"), name="subcluster")


def _x(**cols):
    return pd.DataFrame(cols, index=list("wxyz"), dtype=float)


def _row(table, sub):
    return table.set_index("subcluster").loc[sub]


def test_zscore_zero_variance_gene_is_zeros_not_nan():
    z = ts.zscore_genes(_x(A=[3, 3, 1, 1], K=[5, 5, 5, 5]))
    assert z["A"].tolist() == [1.0, 1.0, -1.0, -1.0]  # mean 2, population sd 1
    assert z["K"].tolist() == [0.0, 0.0, 0.0, 0.0]
    assert not z.isna().any().any()


def test_score_is_cd8_mean_minus_cd4_mean_of_z():
    # CD8 genes z = [1,1,-1,-1]; CD4 genes constant -> z 0. s1 = 1 - 0 = +1, s2 = -1 - 0 = -1.
    x = _x(A=[3, 3, 1, 1], B=[3, 3, 1, 1], C=[2, 2, 2, 2], D=[2, 2, 2, 2])
    table, contrib = ts.subcluster_scores(ts.zscore_genes(x), SUB, CD8, CD4)
    assert _row(table, "s1")["score"] == pytest.approx(1.0)
    assert _row(table, "s2")["score"] == pytest.approx(-1.0)
    assert _row(table, "s1")["cd8_z"] == pytest.approx(1.0)
    assert _row(table, "s1")["cd4_z"] == pytest.approx(0.0)
    assert _row(table, "s1")["n_cells"] == 2
    c = contrib.set_index(["subcluster", "gene"])
    assert c.loc[("s1", "A"), "mean_z"] == pytest.approx(1.0)
    assert c.loc[("s1", "A"), "side"] == "CD8"
    assert c.loc[("s2", "C"), "side"] == "CD4"


def test_z_scoring_stops_high_count_gene_dominating():
    # CD4 gene C has 100x the counts: raw means say CD4 in s1 (2 - 100 < 0) but s1 is CD8-high by z.
    x = _x(A=[2, 2, 0, 0], C=[100, 100, 300, 300])
    raw_s1 = x.loc[["w", "x"], "A"].mean() - x.loc[["w", "x"], "C"].mean()
    assert raw_s1 < 0
    table, _ = ts.subcluster_scores(ts.zscore_genes(x), SUB, ["A"], ["C"])
    assert _row(table, "s1")["score"] == pytest.approx(2.0)  # 1 - (-1)
    assert _row(table, "s2")["score"] == pytest.approx(-2.0)


def test_score_zero_is_cd4():
    assert ts.split_label(0.0) == "CD4+ T cells"
    assert ts.split_label(1e-9) == "CD8+ T cells"
    assert ts.split_label(-1.0) == "CD4+ T cells"


def test_ambiguous_flag_below_0_7_still_gets_sign_label():
    # CD8 genes: A z=[1,1,-1,-1], B constant -> mean 0.5 in s1; CD4 genes constant. Score s1 = +0.5.
    x = _x(A=[3, 3, 1, 1], B=[2, 2, 2, 2], C=[2, 2, 2, 2], D=[2, 2, 2, 2])
    table, _ = ts.subcluster_scores(ts.zscore_genes(x), SUB, CD8, CD4)
    r = _row(table, "s1")
    assert r["score"] == pytest.approx(0.5)
    assert abs(r["score"]) < gt_config.AMBIGUOUS_ABS_SCORE
    assert bool(r["ambiguous"]) is True
    assert ts.split_label(r["score"]) == "CD8+ T cells"


def test_absent_genes_dropped_and_listed():
    x = _x(A=[3, 3, 1, 1], C=[1, 1, 3, 3])
    table, contrib = ts.subcluster_scores(ts.zscore_genes(x), SUB, ["A", "CD8B"], ["C"])
    assert table.attrs["absent_genes"] == ["CD8B"]
    assert _row(table, "s1")["genes_used"] == "A;C"
    assert "CD8B" not in set(contrib["gene"])
    # score uses only present genes: A z=+1, C z=-1 in s1 -> 1 - (-1) = 2
    assert _row(table, "s1")["score"] == pytest.approx(2.0)
    with pytest.raises(ValueError):
        ts.subcluster_scores(ts.zscore_genes(x), SUB, ["CD8B"], ["C"])


# ---- run_split on a synthetic panel -----------------------------------------------------------

BG = [f"G{i}" for i in range(40)]
PANEL = BG + ["CD8A", "GZMA", "CD4", "IL7R"]


def _adata(n_per=50, seed=0):
    rng = np.random.default_rng(seed)
    blocks, names = [], []
    # cluster "0": CD8-high T; "1": CD4-high T; "2": non-T (different background program)
    for k, (hi_genes, hi_idx) in enumerate([(["CD8A", "GZMA"], None), (["CD4", "IL7R"], None), (None, range(0, 10))]):
        mu = np.full(len(PANEL), 0.3)
        if hi_genes:
            mu[[PANEL.index(g) for g in hi_genes]] = 8.0
        else:
            mu[list(hi_idx)] = 8.0
        blocks.append(rng.poisson(mu, size=(n_per, len(PANEL))))
        names += [f"c{k}_{i}" for i in range(n_per)]
    a = ad.AnnData(csr_matrix(np.vstack(blocks).astype(np.float32)))
    a.obs_names = names
    a.var_names = PANEL
    level1 = pd.Series([n.split("_")[0][1:] for n in names], index=names, name="cluster")
    cn = pd.DataFrame({"cluster": ["0", "1", "2"], "label": ["T cells", "T cells", "Epithelial"]})
    return a, level1, cn


def test_every_t_cell_gets_exactly_one_label():
    a, level1, cn = _adata()
    out = ts.run_split(a, level1, cn, seed=0)
    assert list(out.columns) == ["cell_id", "subcluster", "label7"]
    assert out["cell_id"].is_unique
    assert out["label7"].isin(["CD8+ T cells", "CD4+ T cells"]).all()
    assert not out.isna().any().any()
    assert set(out["cell_id"]) == set(level1.index[level1.isin(["0", "1"])])
    lab = out.set_index("cell_id")["label7"]
    assert (lab[[c for c in lab.index if c.startswith("c0_")]] == "CD8+ T cells").all()
    assert (lab[[c for c in lab.index if c.startswith("c1_")]] == "CD4+ T cells").all()


def test_non_t_clusters_untouched():
    a, level1, cn = _adata()
    x_before = a.X.copy().toarray()
    obs_before = list(a.obs_names)
    l1_before = level1.copy()
    out = ts.run_split(a, level1, cn, seed=0)
    assert not any(c.startswith("c2_") for c in out["cell_id"])
    assert np.array_equal(a.X.toarray(), x_before)
    assert list(a.obs_names) == obs_before
    pd.testing.assert_series_equal(level1, l1_before)


def test_run_split_full_matches_run_split_and_returns_tables():
    a, level1, cn = _adata()
    labels, table, contrib, absent = ts.run_split_full(a, level1, cn, seed=0)
    pd.testing.assert_frame_equal(labels, ts.run_split(a, level1, cn, seed=0))
    assert list(table.columns) == [
        "subcluster",
        "n_cells",
        "cd8_z",
        "cd4_z",
        "score",
        "ambiguous",
        "genes_used",
    ]
    assert list(contrib.columns) == ["subcluster", "gene", "side", "mean_z"]
    assert len(table) > 0 and len(contrib) > 0
    assert set(table["subcluster"]) == set(labels["subcluster"])
    assert isinstance(absent, list)
    assert "GZMB" in absent  # in the recipe list but not in the synthetic panel


def test_level1_cell_missing_from_adata_raises_clear_error():
    a, level1, cn = _adata()
    level1 = pd.concat([level1, pd.Series(["0"], index=["ghost"], name="cluster")])
    with pytest.raises(ValueError, match="not in adata"):
        ts.run_split(a, level1, cn, seed=0)
