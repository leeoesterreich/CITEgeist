"""Tests for gt_rebuild/marker_names.py (synthetic, hand-computed)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import cluster_names, gt_config  # noqa: E402
from gt_rebuild import marker_names as mn  # noqa: E402

pytestmark = pytest.mark.unit

SETS = {"A": ["a1", "a2"], "B": ["b1", "b2"], "C": ["c1"]}


def _fixture():
    """4 clusters x 2 cells; cluster means (rows 0..3) are set directly so z-scores can be computed by hand."""
    means = pd.DataFrame(
        {
            "a1": [4.0, 0.0, 0.0, 2.0],
            "a2": [4.0, 0.0, 0.0, 2.0],
            "b1": [0.0, 4.0, 0.0, 2.0],
            "b2": [0.0, 4.0, 0.0, 2.0],
            "c1": [0.0, 0.0, 5.0, 5.0],
        },
        index=["0", "1", "2", "3"],
    )
    rows, ids = [], []
    for c in means.index:
        for k in range(2):  # two cells: mean-preserving +-1 split
            rows.append(means.loc[c].values + (1 if k else -1))
            ids.append(f"{c}_{k}")
    logx = pd.DataFrame(rows, index=ids, columns=means.columns)
    clusters = pd.Series([i.split("_")[0] for i in ids], index=ids)
    return logx, clusters


def test_scores_match_hand_computation_and_argmax():
    logx, cl = _fixture()
    out = mn.score_clusters(logx, cl, SETS).set_index("cluster")
    # z over clusters with the population sd (ddof=0), computed independently with numpy
    a = np.array([4.0, 0, 0, 2])
    za = (a - a.mean()) / a.std()
    b = np.array([0.0, 4, 0, 2])
    zb = (b - b.mean()) / b.std()
    c = np.array([0.0, 0, 5, 5])
    zc = (c - c.mean()) / c.std()
    assert out.loc["0", "score_A"] == pytest.approx(za[0])
    assert out.loc["1", "score_B"] == pytest.approx(zb[1])
    assert out.loc["3", "score_C"] == pytest.approx(zc[3])
    assert list(out["label"]) == ["A", "B", "C", "C"]
    assert out.loc["2", "top_score"] == pytest.approx(zc[2])
    assert out.loc["0", "n_cells"] == 2
    assert out.loc["0", "margin"] == pytest.approx(out.loc["0", "top_score"] - out.loc["0", "second_score"])


def test_small_margin_is_ambiguous_but_labelled():
    logx, cl = _fixture()
    out = mn.score_clusters(logx, cl, SETS, margin=0.25).set_index("cluster")
    # cluster 3 has A and B and C all positive-ish; label is its argmax and ambiguity is margin < 0.25
    for c, r in out.iterrows():
        assert bool(r["ambiguous"]) == (r["margin"] < 0.25)
    wide = mn.score_clusters(logx, cl, SETS, margin=100.0)
    assert wide["ambiguous"].all() and wide["label"].notna().all()


def test_zero_sd_gene_gives_zero_not_nan():
    logx, cl = _fixture()
    logx["a2"] = 3.0  # constant across clusters
    out = mn.score_clusters(logx, cl, SETS)
    assert not out.filter(like="score_").isna().any().any()
    z = out.attrs["gene_z"]
    assert (z["a2"] == 0).all() and not z.isna().any().any()


def test_missing_marker_gene_raises():
    logx, cl = _fixture()
    with pytest.raises(ValueError, match="b2"):
        mn.score_clusters(logx.drop(columns=["b2"]), cl, SETS)


def test_every_cluster_gets_one_real_type_with_real_sets():
    genes = sorted({g for v in gt_config.MARKER_SETS.values() for g in v})
    rng = np.random.RandomState(0)
    ids = [f"c{i}" for i in range(30)]
    logx = pd.DataFrame(rng.rand(30, len(genes)), index=ids, columns=genes)
    cl = pd.Series([str(i % 5) for i in range(30)], index=ids)
    out = mn.score_clusters(logx, cl, gt_config.MARKER_SETS)
    assert len(out) == 5 and set(out["label"]) <= set(gt_config.CT6)
    assert sorted(c for c in out.columns if c.startswith("score_")) == sorted(f"score_{t}" for t in gt_config.CT6)


def test_recipe_carries_marker_rules():
    assert gt_config.RECIPE["marker_sets"] == gt_config.MARKER_SETS
    assert gt_config.RECIPE["marker_ambiguous_margin"] == 0.25
    assert set(gt_config.MARKER_SETS) == set(gt_config.CT6)
    assert gt_config.MARKER_SETS["T cells"] == ["CD3D", "CD3E", "CD3G", "TRAC", "CD2"]
    assert sum(len(v) for v in gt_config.MARKER_SETS.values()) == 40


def test_config_resolution_and_variants():
    r = gt_config.RECIPE
    assert r["level1"]["leiden_resolution"] == 2.0
    assert r["stability_variants"] == {
        "leiden_1.5": {"resolution": 1.5, "seed": 0},
        "leiden_2.5": {"resolution": 2.5, "seed": 0},
        "seed_1": {"resolution": 2.0, "seed": 1},
    }
    assert r["tcell"] == {"n_pcs": 20, "n_neighbors": 15, "leiden_resolution": 1.0}


def test_written_files_read_back_by_reader(tmp_path):
    logx, cl = _fixture()
    sets = {t: [f"g{i}"] for i, t in enumerate(gt_config.CT6)}
    rng = np.random.RandomState(1)
    ids = [f"c{i}" for i in range(24)]
    big = pd.DataFrame(rng.rand(24, 6), index=ids, columns=[f"g{i}" for i in range(6)])
    cls = pd.Series([str(i % 4) for i in range(24)], index=ids)
    mn.write_outputs(mn.score_clusters(big, cls, sets), tmp_path)
    names = cluster_names.read_cluster_names(tmp_path / "cluster_names.csv")
    assert len(names) == 4 and "delta_next" in names.columns
    z = pd.read_csv(tmp_path / "marker_gene_z.csv", index_col=0)
    assert z.shape == (4, 6) and not z.isna().any().any()


def test_stage_name_writes_files_from_synthetic_adata(tmp_path, monkeypatch):
    anndata = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    from gt_rebuild import run_stage

    genes = sorted({g for v in gt_config.MARKER_SETS.values() for g in v})
    rng = np.random.RandomState(0)
    ids = [f"c{i}" for i in range(40)]
    a = anndata.AnnData(rng.poisson(2, (40, len(genes))).astype(np.float32))
    a.obs_names, a.var_names = ids, genes
    monkeypatch.setattr(run_stage, "_load_adata", lambda: (a, None))
    pd.DataFrame({"cell_id": ids, "cluster": [str(i % 4) for i in range(40)]}).to_csv(tmp_path / "l1.csv", index=False)
    run_stage.stage_name(tmp_path / "l1.csv", tmp_path / "out")
    names = cluster_names.read_cluster_names(tmp_path / "out" / "cluster_names.csv")
    assert len(names) == 4 and set(names["label"]) <= set(gt_config.CT6)
    assert (tmp_path / "out" / "marker_gene_z.csv").exists()
