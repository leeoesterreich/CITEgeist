"""Unit tests for the SACE proportion-source ablation glue: proportion alignment, SPOTS arm, matched-gene restriction."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

D = Path(__file__).resolve().parent.parent / "Benchmarking" / "xenium_benchmarking" / "spots_baseline"
sys.path.insert(0, str(D))
import eval_sace_arms as E  # noqa: E402
import sace_props_ablation as S  # noqa: E402

TYPES = ["A", "B", "C"]


def test_align_props_order_missing_rows_and_normalization():
    df = pd.DataFrame({"C": [1.0, 0.0], "A": [1.0, 0.0], "B": [2.0, 0.0], "Z": [9.0, 9.0]}, index=["s1", "s2"])
    out = S.align_props(df, ["s2", "s1", "s3"], TYPES)
    assert list(out.columns) == TYPES and list(out.index) == ["s2", "s1", "s3"]  # column order = types, row order = spots
    assert np.allclose(out.loc["s1"], [0.25, 0.5, 0.25])  # normalized; extra column Z dropped before normalization
    assert (out.loc["s2"] == 0).all() and (out.loc["s3"] == 0).all()  # all-zero / absent spot stays zero (no NaN, no uniform)
    assert not out.isna().any().any()


def test_spots_arm_props_uses_raw_cite_counts_and_aligns_to_sace_spots():
    profile = {"A": {"Major": ["a1"], "Minor": []}, "B": {"Major": ["b1"], "Minor": []}, "C": {"Major": ["c1"], "Minor": []}}
    X = np.array([[3.0, 0, 4.0], [0.0, 0.0, 0.0], [1.0, 0, 0]])  # s1 mixed, s2 empty -> uniform, s3 pure A
    cite = ad.AnnData(X=X, obs=pd.DataFrame(index=["s1", "s2", "s3"]), var=pd.DataFrame(index=["a1", "b1", "c1"]))
    out = S.spots_arm_props(cite, profile, ["s3", "s1"])  # SACE dropped s2; order differs from CITE
    assert list(out.index) == ["s3", "s1"]
    assert np.allclose(out.loc["s3"], [1, 0, 0])
    assert np.allclose(out.loc["s1"], [0.6 / 1.4, 0, 0.8 / 1.4])
    assert np.allclose(out.sum(axis=1), 1.0)


def test_matched_gene_restriction(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(40)]
    spots = [f"spot_{i}" for i in range(30)]
    gt = pd.DataFrame(rng.poisson(3, (30, 40)).astype(float), index=spots, columns=genes)
    pred = gt.copy()
    pred.iloc[:, 20:] = rng.poisson(3, (30, 20)).astype(float)  # genes 20+ are noise; genes 0-19 are perfect
    (tmp_path / "Xenium_region_0").mkdir()
    gt.to_csv(tmp_path / "Xenium_region_0" / "B_cells_GT.csv")
    ld = tmp_path / "layers"
    ld.mkdir()
    pred.to_csv(ld / "B_cells_layer.csv")
    monkeypatch.setattr(E.C, "GT_ROOT", tmp_path)
    full = E.block_metrics(ld, 0, "B_cells")
    sub = E.block_metrics(ld, 0, "B_cells", set(genes[:20]))
    assert sub["n_genes"] == 20 and full["n_genes"] == 40
    assert sub["flattened_r"] > 0.999 and full["flattened_r"] < sub["flattened_r"]  # the subset really restricts the genes scored


def test_matched_pooled_r_restricts_to_gene_set_and_pools_blocks(tmp_path, monkeypatch):
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(30)]
    spots = [f"spot_{i}" for i in range(25)]
    gt_root, lay = tmp_path / "gt", tmp_path / "layers"
    for reg in range(5):
        (gt_root / f"Xenium_region_{reg}").mkdir(parents=True)
        (lay / f"Xenium_region_{reg}" / "layers").mkdir(parents=True)
        gt = pd.DataFrame(rng.poisson(4, (25, 30)).astype(float), index=spots, columns=genes)
        pred = gt.copy()
        pred.iloc[:, 15:] = rng.poisson(4, (25, 15)).astype(float)  # genes 15+ are noise
        gt.to_csv(gt_root / f"Xenium_region_{reg}" / "T_GT.csv")
        pred.to_csv(lay / f"Xenium_region_{reg}" / "layers" / "T_layer.csv")
    monkeypatch.setattr(E.C, "GT_ROOT", gt_root)
    monkeypatch.setattr(E.C, "CELL_TYPES", ["T"])
    good = E.matched_pooled_r(lay, set(genes[:15]))
    allg = E.matched_pooled_r(lay, set(genes))
    assert good > 0.999 and allg < good  # restriction to the supplied gene set is honoured, case-insensitively upper-cased
    assert abs(E.matched_pooled_r(lay, {g.lower().upper() for g in genes[:15]}) - good) < 1e-12


def test_qp_saved_props_reads_run_dir_and_drops_recon_error(tmp_path, monkeypatch):
    import json

    (tmp_path / "res").mkdir()
    (tmp_path / "res" / "r.json").write_text(json.dumps({"_inputs": {"run_dirs": {"2": "run2"}}}))
    (tmp_path / "run2" / "predictions").mkdir(parents=True)
    pd.DataFrame({"B": [1.0, 3.0], "A": [1.0, 1.0], "C": [2.0, 0.0], "recon_error": [9.0, 9.0]}, index=["s1", "s2"]).to_csv(
        tmp_path / "run2" / "predictions" / "baseline.csv.gz"
    )
    monkeypatch.setattr(S, "REPO", tmp_path)
    monkeypatch.setattr(S, "XB", tmp_path)
    monkeypatch.setattr(S.Path, "read_text", lambda self, *a, **k: open(tmp_path / "res" / "r.json").read() if self.name == "CITEgeist_7type_results.json" else open(self).read())
    out = S.qp_saved_props(2, ["s2", "s1"], TYPES)
    assert list(out.columns) == TYPES and list(out.index) == ["s2", "s1"]
    assert np.allclose(out.loc["s1"], [0.25, 0.25, 0.5]) and np.allclose(out.loc["s2"], [0.25, 0.75, 0.0])


def test_paired_summary_uses_only_blocks_common_to_all_arms():
    a = pd.DataFrame({"region": [0, 0, 1], "cell_type": ["x", "y", "x"], "flat_r": [0.1, 0.9, 0.5], "rmse": [1.0, 1.0, 1.0]})
    b = pd.DataFrame({"region": [0, 1], "cell_type": ["x", "x"], "flat_r": [0.3, 0.7], "rmse": [2.0, 2.0]})  # lacks block (0, y)
    out = E.paired_summary({"a": a, "b": b})
    assert out["n_common_blocks"] == 2
    assert abs(out["a"]["flat_r"] - 0.3) < 1e-9 and abs(out["b"]["flat_r"] - 0.5) < 1e-9  # a's 0.9 block is excluded


@pytest.mark.unit
def test_marker_share_arm_reads_written_preds(tmp_path, monkeypatch):
    d = tmp_path / "Benchmarking/xenium_benchmarking/CITEgeist/output_marker_share_7type/Xenium_region_3/predictions"
    d.mkdir(parents=True)
    pd.DataFrame({"B": [1.0, 3.0], "A": [1.0, 1.0], "C": [2.0, 0.0]}, index=["s1", "s2"]).to_csv(d / "baseline.csv.gz")
    monkeypatch.setenv("MAIN", str(tmp_path))
    out = S.marker_share_props(3, ["s2", "s1", "s9"], TYPES)
    assert list(out.columns) == TYPES and list(out.index) == ["s2", "s1", "s9"]
    assert np.allclose(out.loc["s1"], [0.25, 0.25, 0.5]) and (out.loc["s9"] == 0).all()
    assert "marker_share" in S.ARMS
    assert S.arm_region_dir("out", "marker_share", "R") == Path("out/R")  # flat: matches compare_xenium_gex_methods layout
    assert S.arm_region_dir("out", "spots", "R") == Path("out/spots/R")
