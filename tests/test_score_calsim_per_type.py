"""score_calsim_per_type.py: per-type mean r on synthetic sim CSVs, the 6-type merge, and the GT-GEX sigma rule."""

from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
SIMB = REPO / "Benchmarking" / "simulation_benchmarking"

pytestmark = pytest.mark.unit


def _mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# stage2_gate.py is a script that imports its sibling `score_calsim_per_type` (its own directory is on sys.path when run); make that true under pytest too,
# otherwise the gate test passes or fails depending on which other test file put this directory on sys.path first
if str(SIMB) not in sys.path:
    sys.path.insert(0, str(SIMB))
S = _mod("score_calsim_per_type_t", SIMB / "score_calsim_per_type.py")
GATE = _mod("stage2_gate_t", SIMB / "stage2_gate.py")


def _props(rng, n, noise=0.0, base=None):
    x = rng.dirichlet(np.ones(9), size=n) if base is None else base + noise * rng.random((n, 9))
    x = x / x.sum(axis=1, keepdims=True)
    return pd.DataFrame(x, index=[f"s{i}" for i in range(n)], columns=S.TYPES_9)


def make_sim(tmp, methods=S.METHODS, n_spots=40, seed=0, gex=True):
    """GT + canonical props for every method (method i = GT + increasing noise), and small GT GEX layers."""
    rng = np.random.default_rng(seed)
    root, canon = tmp / "replicates_calsim_v2", tmp / "regen" / "canonical_calsim_v2"
    for c in S.CONDITIONS:
        for r in range(S.N_REPS):
            gt = _props(rng, n_spots)
            gt["spot_x"], gt["spot_y"] = 1.0, 2.0
            (root / c / "ST_sim").mkdir(parents=True, exist_ok=True)
            gt.to_csv(root / c / "ST_sim" / f"Wu_ST_{r}_prop.csv")
            gt = gt.drop(columns=["spot_x", "spot_y"])
            for i, m in enumerate(methods):
                pr = _props(rng, n_spots, noise=0.05 * (i + 1), base=gt.to_numpy())
                pr = pr.iloc[::-1]  # different row order: the scorer must align on the GT index
                d = canon / m / c / f"Wu_rep_{r}"
                d.mkdir(parents=True, exist_ok=True)
                pr.to_csv(d / "props.csv")
            if gex:
                lay = root / c / "ST_GEX_sim" / f"sample_{r}" / "layers"
                lay.mkdir(parents=True, exist_ok=True)
                for t in S.TYPES_9:
                    genes = ["g1", "g2", "g3", "zero"] + (["extra"] if t == "PVL" else [])
                    v = rng.poisson(3, size=(n_spots, len(genes))).astype(float)
                    v[:, genes.index("zero")] = 0
                    pd.DataFrame(v, index=gt.index, columns=genes).to_csv(lay / f"{t}_GT.csv")
    return canon, root


def test_sim_type_map_equals_generator():
    tree = ast.parse((REPO / "Benchmarking/simulation_generation/generate_protein_sim.py").read_text())
    node = next(
        n
        for n in tree.body
        if isinstance(n, ast.Assign) and any(getattr(t, "id", "") == "SIM_TYPE_MAP" for t in n.targets)
    )
    assert ast.literal_eval(node.value) == S.SIM_TYPE_MAP
    assert len(S.TYPES_9) == 9 and len(S.TYPES_6) == 6
    assert S.TYPES_6 == sorted(["B cells", "Endothelial", "Epithelial", "Fibroblasts", "Macrophages", "T cells"])


def test_per_type_mean_r_matches_stage2_gate(tmp_path):
    canon, root = make_sim(tmp_path, gex=False)
    res = S.score(canon, root)
    n = 0
    for m in S.METHODS:
        for c in S.CONDITIONS:
            cell = res["methods"][m][c]
            assert cell["n_reps"] == 5 and len(cell["per_type_r_9"]["per_rep"]) == 5
            for r in range(5):
                want = GATE._per_type_mean_r(
                    canon / m / c / f"Wu_rep_{r}" / "props.csv", root / c / "ST_sim" / f"Wu_ST_{r}_prop.csv"
                )
                assert cell["per_type_r_9"]["per_rep"][r] == pytest.approx(want, abs=1e-12)
                n += 1
            v = np.array(cell["per_type_r_9"]["per_rep"])
            # headline is type-first: mean of the per-type scores, which is what the by_type field holds
            assert cell["per_type_r_9"]["mean"] == pytest.approx(np.mean(list(cell["per_type_9_by_type"].values())))
            assert cell["per_type_r_6"]["mean"] == pytest.approx(np.mean(list(cell["per_type_6_by_type"].values())))
            # error bar: population SD (ddof=0) of the per-replicate means / sqrt(n), the Xenium rule
            assert cell["per_type_r_9"]["sem"] == pytest.approx(v.std(ddof=0) / np.sqrt(5))
    assert n == 6 * 2 * 5
    # more noise -> lower r (methods were generated with increasing noise)
    means = [res["methods"][m]["high_seg"]["per_type_r_9"]["mean"] for m in S.METHODS]
    assert means[0] > means[-1]


def test_six_type_merge_sums_partner_types(tmp_path):
    df = pd.DataFrame(np.arange(18, dtype=float).reshape(2, 9), columns=S.TYPES_9)
    m6 = S.merge_6(df)
    assert list(m6.columns) == S.TYPES_6
    row = df.iloc[0]
    assert m6.iloc[0]["Epithelial"] == row["Cancer Epithelial"] + row["Normal Epithelial"]
    assert m6.iloc[0]["Fibroblasts"] == row["CAFs"] + row["PVL"]
    assert m6.iloc[0]["B cells"] == row["B-cells"] + row["Plasmablasts"]
    assert m6.iloc[0]["Macrophages"] == row["Myeloid"]
    assert m6.to_numpy().sum() == df.to_numpy().sum()


def test_merge_can_separate_from_9type(tmp_path):
    """A method that swaps partner types is poor on 9 types but perfect on the 6-type merge."""
    rng = np.random.default_rng(3)
    gt = _props(rng, 50)
    swapped = gt.rename(
        columns={
            "CAFs": "PVL",
            "PVL": "CAFs",
            "Cancer Epithelial": "Normal Epithelial",
            "Normal Epithelial": "Cancer Epithelial",
            "B-cells": "Plasmablasts",
            "Plasmablasts": "B-cells",
        }
    )
    r9 = S._mean_defined(S.per_type_r(swapped[S.TYPES_9], gt))
    r6 = S._mean_defined(S.per_type_r(S.merge_6(swapped[S.TYPES_9]), S.merge_6(gt)))
    assert r6 == pytest.approx(1.0) and r9 < 0.9


def test_constant_type_is_left_out_and_counted():
    gt = pd.DataFrame({"a": [0.1, 0.2, 0.3], "b": [0.5, 0.5, 0.5]})
    pr = pd.DataFrame({"a": [0.1, 0.25, 0.3], "b": [0.4, 0.5, 0.6]})
    rs = S.per_type_r(pr, gt)
    assert np.isnan(rs["b"]) and np.isfinite(rs["a"])
    assert S._mean_defined(rs) == rs["a"]


# ---- hand-computed: aggregation order and the two constant rules (numbers worked out on paper, not by the code) ----
# 2 replicates, 3 spots, 3 types. r of GT [1,2,3] against pred [1,2,3] = 1, [3,2,1] = -1, [1,3,2] = 0.5.
#   rep1: A pred=GT -> 1;  B pred [1,3,2] -> 0.5;  C pred constant -> 0 (constant-prediction rule)
#   rep2: A pred [3,2,1] -> -1;  B GT constant -> undefined;  C pred=GT -> 1
# type-first, const pred = 0 (the Xenium rule): A (1-1)/2=0, B 0.5 (rep1 only), C (0+1)/2=0.5 -> mean 1/3
# per-rep means over defined types: rep1 (1+0.5+0)/3 = 0.5, rep2 (-1+1)/2 = 0 -> SD (ddof=0) 0.25, s.e.m. 0.25/sqrt(2)
# the three wrong readings each give a different number:
#   type-first but constant prediction DROPPED: C = 1 (rep2 only) -> (0+0.5+1)/3 = 0.5
#   replicate-first, const pred = 0:  (0.5 + 0)/2 = 0.25
#   replicate-first, const pred dropped: rep1 (1+0.5)/2 = 0.75, rep2 0 -> 0.375
def _hand_reps():
    idx = ["s1", "s2", "s3"]
    g, flat = [1.0, 2.0, 3.0], [2.0, 2.0, 2.0]  # (scaled by 0.1 below: r is scale-free, and the evaluator renormalizes rows summing > 1.5)
    gts = [
        pd.DataFrame({"A": g, "B": g, "C": g}, index=idx),
        pd.DataFrame({"A": g, "B": flat, "C": g}, index=idx),
    ]
    prs = [
        pd.DataFrame({"A": g, "B": [1.0, 3.0, 2.0], "C": flat}, index=idx),
        pd.DataFrame({"A": [3.0, 2.0, 1.0], "B": g, "C": g}, index=idx),
    ]
    return [x * 0.1 for x in gts], [x * 0.1 for x in prs]


def test_hand_computed_order_and_constant_rules():
    gts, prs = _hand_reps()
    by_rep = [S.per_type_r(p, g) for p, g in zip(prs, gts)]
    assert by_rep[0]["A"] == pytest.approx(1.0) and by_rep[0]["B"] == pytest.approx(0.5) and by_rep[0]["C"] == 0.0
    assert by_rep[1]["A"] == pytest.approx(-1.0) and np.isnan(by_rep[1]["B"]) and by_rep[1]["C"] == pytest.approx(1.0)
    agg = S.aggregate_per_type_r(by_rep)
    assert agg["by_type"] == pytest.approx({"A": 0.0, "B": 0.5, "C": 0.5})
    assert agg["mean"] == pytest.approx(1 / 3)
    assert agg["per_rep"] == pytest.approx([0.5, 0.0])
    assert agg["sd"] == pytest.approx(0.25) and agg["sem"] == pytest.approx(0.25 / np.sqrt(2))
    assert agg["n_undefined_type_r"] == 1 and agg["n_types_unscored"] == 0
    # none of the wrong readings gives the right number
    for wrong in (0.5, 0.25, 0.375):
        assert agg["mean"] != pytest.approx(wrong)


def test_type_undefined_in_every_replicate_is_dropped_and_counted():
    gts, prs = _hand_reps()
    for g in gts:
        g["B"] = 2.0  # B constant in the GT of both replicates
    agg = S.aggregate_per_type_r([S.per_type_r(p, g) for p, g in zip(prs, gts)])
    assert "B" not in agg["by_type"] and agg["n_types_unscored"] == 1
    assert agg["mean"] == pytest.approx(np.mean([0.0, 0.5]))  # A 0, C 0.5


def test_evaluator_agrees_with_the_simulation_scorer_on_the_hand_example(tmp_path):
    """The Xenium strict evaluator (frozen, shipped) must give the same 1/3 on the same data: one definition."""
    sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking" / "evaluation" / "src"))
    import evaluate_benchmark as EB

    gts, prs = _hand_reps()
    (tmp_path / "gt" / "ground_truth").mkdir(parents=True)
    (tmp_path / "pred").mkdir()
    for i, (g, p) in enumerate(zip(gts, prs)):
        g.to_csv(tmp_path / "gt" / "ground_truth" / f"Xenium_region_{i}_prop.csv")
        p.to_csv(tmp_path / "pred" / f"Xenium_region_{i}_cell_prop_global_results.csv")
    types = ["A", "B", "C"]
    ev = EB.evaluability_table(tmp_path / "gt", types, n_regions=2)
    masks = EB.common_spot_masks(tmp_path / "gt", {"m": str(tmp_path / "pred")}, ev, n_regions=2)
    res = EB.evaluate_strict(tmp_path / "gt", str(tmp_path / "pred"), types, masks, ev, n_regions=2)
    assert res["summary"]["per_type_mean_r"] == pytest.approx(1 / 3)
    assert res["summary"]["n_zero_const_pred"] == 1


def test_missing_type_or_spot_raises(tmp_path):
    canon, root = make_sim(tmp_path, methods=["RCTD"], gex=False)
    f = canon / "RCTD" / "mixed" / "Wu_rep_2" / "props.csv"
    pd.read_csv(f, index_col=0).drop(columns=["PVL"]).to_csv(f)
    with pytest.raises(ValueError, match="9 sim types"):
        S.score(canon, root, methods=["RCTD"])
    canon, root = make_sim(tmp_path / "b", methods=["RCTD"], gex=False)
    f = canon / "RCTD" / "high_seg" / "Wu_rep_0" / "props.csv"
    pd.read_csv(f, index_col=0).iloc[1:].to_csv(f)
    with pytest.raises(ValueError, match="GT spots missing"):
        S.score(canon, root, methods=["RCTD"])


def test_undefined_replicate_raises():
    nan = float("nan")
    with pytest.raises(ValueError, match="undefined"):  # a replicate with no defined type has no per-replicate mean
        S.aggregate_per_type_r([{"a": 0.5, "b": 0.2}, {"a": nan, "b": nan}])
    with pytest.raises(ValueError):
        S.aggregate_per_type_r([{"a": nan}, {"a": nan}])


def test_gene_sigma_is_one_sd_per_gene_per_block(tmp_path):
    _canon, root = make_sim(tmp_path, methods=["RCTD"])
    sg = S.gene_sigma(root)
    files = sorted(root.glob("*/ST_GEX_sim/sample_*/layers/*_GT.csv"))
    assert len(files) == 2 * 5 * 9
    vals = []
    for f in files:  # independent loop: one SD per gene inside each file (block), genes expressed in that block, 4-gene intersection
        raw = pd.read_csv(f, index_col=0)[["g1", "g2", "g3", "zero"]].to_numpy()
        for k in range(raw.shape[1]):
            if raw[:, k].sum() > 0:
                vals.append(np.log1p(raw[:, k]).std(ddof=1))
    vals = np.array(vals)
    assert sg["n_genes"] == 4 and sg["n_expressed"] == 3 and sg["n_blocks"] == 90
    assert sg["n_obs"] == len(vals) == 270  # 'zero' is never expressed; 'extra' (PVL only) is outside the intersection
    assert sg["mean_expressed"] == pytest.approx(vals.mean())
    assert sg["p90_expressed"] == pytest.approx(np.percentile(vals, 90))
    # the old rule (one SD over rows pooled across all files) gives a different number on these data
    stack = np.vstack([np.log1p(pd.read_csv(f, index_col=0)[["g1", "g2", "g3"]].to_numpy()) for f in files])
    assert sg["mean_expressed"] != pytest.approx(stack.std(axis=0, ddof=1).mean(), abs=1e-4)


def test_gene_sigma_hand_computed_literal(tmp_path):
    lay = tmp_path / "c" / "ST_GEX_sim" / "sample_0" / "layers"
    lay.mkdir(parents=True)
    # block A: g1 counts 0,1,3 -> log1p 0, ln2, ln4 -> SD = ln 2 = 0.693147; g2 all zero -> not expressed, left out
    pd.DataFrame({"g1": [0, 1, 3], "g2": [0, 0, 0]}).to_csv(lay / "A_GT.csv")
    # block B: g1 2,2,2 -> SD 0; g2 1,3,5 -> log1p ln2, ln4, ln6 -> SD 0.555549
    pd.DataFrame({"g1": [2, 2, 2], "g2": [1, 3, 5]}).to_csv(lay / "B_GT.csv")
    sg = S.gene_sigma(tmp_path, conditions=("c",), n_reps=1)
    assert sg["n_obs"] == 3 and sg["n_blocks"] == 2 and sg["n_expressed"] == 2
    assert sg["mean_expressed"] == pytest.approx((0.693147 + 0.0 + 0.555549) / 3, abs=1e-5)  # 0.416232
    assert sg["median_expressed"] == pytest.approx(0.555549, abs=1e-5)


def test_main_writes_stamped_new_file_and_refuses_overwrite(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "scripts"))
    import provenance

    monkeypatch.setattr(
        provenance, "_code_fields", lambda p: {"_produced_by": p, "_code_commit": "c" * 40, "_code_dirty": False}
    )
    canon, root = make_sim(tmp_path)
    out = tmp_path / "out.json"
    S.main(["--canonical", str(canon), "--sim-root", str(root), "--out", str(out)])
    ok, _msg = provenance.verify(out)
    assert ok
    import json

    d = json.loads(out.read_text())
    assert d["_metric"] == S.METRIC and d["_produced_by"] == S.PRODUCER and d["_code_dirty"] is False
    assert set(d["methods"]) == set(S.METHODS) and d["gt_gex_sigma"]["n_expressed"] == 3
    assert d["merge_6"]["Fibroblasts"] == ["CAFs", "PVL"]
    with pytest.raises(SystemExit, match="exists"):
        S.main(["--canonical", str(canon), "--sim-root", str(root), "--out", str(out)])


def test_reference_only_scores_no_citegeist(tmp_path, monkeypatch):
    """v41 Figure 3: --reference-only scores the five reference-based methods only (CITEgeist dir may be absent)."""
    sys.path.insert(0, str(REPO / "scripts"))
    import json

    import provenance

    monkeypatch.setattr(
        provenance, "_code_fields", lambda p: {"_produced_by": p, "_code_commit": "c" * 40, "_code_dirty": False}
    )
    assert S.REFERENCE_METHODS == ["Cell2Location", "RCTD", "Seurat", "CARD", "Tangram"]
    canon, root = make_sim(tmp_path, methods=S.REFERENCE_METHODS)
    out = tmp_path / "ref.json"
    S.main(["--canonical", str(canon), "--sim-root", str(root), "--out", str(out), "--reference-only"])
    d = json.loads(out.read_text())
    assert set(d["methods"]) == set(S.REFERENCE_METHODS) == set(d["props_md5"])
    assert d["methods_scored"] == S.REFERENCE_METHODS
    with pytest.raises(FileNotFoundError):  # without the flag the missing CITEgeist dir is an error, not a skip
        S.main(["--canonical", str(canon), "--sim-root", str(root), "--out", str(tmp_path / "all.json")])
