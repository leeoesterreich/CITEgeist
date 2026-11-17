"""Tests for gt_rebuild/composition_arms.py (synthetic data only)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XB))

import benchmark_paths as bp  # noqa: E402
from gt_rebuild import composition_arms as ca  # noqa: E402

pytestmark = pytest.mark.unit

CT = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]
NATURAL = dict(zip(CT, [10, 30, 30, 10, 10, 10]))
GT = dict(zip(CT, [0.05, 0.30, 0.25, 0.15, 0.15, 0.10]))
OPTION1 = ("natural", "gt_matched", "uniform", "fibro_skew")
OPTION2 = ("natural", "gt_matched", "fibro_skew")


@pytest.fixture
def option1(monkeypatch):
    monkeypatch.setattr(ca, "COMPOSITION_ARM_SET", OPTION1)


def test_fibro_skew_vector(option1):
    v = ca.arm_vectors(NATURAL, GT)["fibro_skew"]
    assert v["Fibroblasts"] == pytest.approx(0.5)
    for t in CT:
        if t != "Fibroblasts":
            assert v[t] == pytest.approx(NATURAL[t] / 100 * 0.5 / 0.9)
    assert abs(sum(v.values()) - 1) < 1e-12


def test_uniform_vector_is_one_sixth(option1):
    v = ca.arm_vectors(NATURAL, GT)["uniform"]
    assert all(x == pytest.approx(1 / 6) for x in v.values())


def test_gt_matched_vector_is_normalised_gt(option1):
    v = ca.arm_vectors(NATURAL, GT)["gt_matched"]
    assert v == pytest.approx(GT)


def test_arm_set_option2_drops_uniform_and_raises_n(monkeypatch):
    avail = {t: 1000 for t in CT}
    avail["B cells"] = 60
    monkeypatch.setattr(ca, "COMPOSITION_ARM_SET", OPTION1)
    v1 = ca.arm_vectors(NATURAL, GT)
    n1 = ca.largest_common_n(v1, avail)
    monkeypatch.setattr(ca, "COMPOSITION_ARM_SET", OPTION2)
    v2 = ca.arm_vectors(NATURAL, GT)
    assert "uniform" not in v2 and tuple(v2) == OPTION2
    assert ca.largest_common_n(v2, avail) > n1


def test_default_is_option2():
    assert ca.COMPOSITION_ARM_SET == OPTION2


def test_no_hardcoded_arm_names_outside_arm_set():
    assert bp.ARMS == ca.COMPOSITION_ARM_SET
    pkg = XB / "gt_rebuild"
    for f in pkg.glob("*.py"):
        if f.name == "composition_arms.py":
            continue
        assert "fibro_skew" not in f.read_text(), f.name
    body = (XB / "benchmark_paths.py").read_text()
    assert "fibro_skew" not in body and "gt_matched" not in body


def test_comp_helpers_validate_arm():
    for fn, args in (
        (bp.comp_ref_dir, ("rctd", "bogus", 0)),
        (bp.comp_c2l_model, ("bogus", 0)),
        (bp.comp_pred_dir, ("RCTD", "bogus", 0)),
    ):
        with pytest.raises(ValueError):
            fn(*args)
    assert bp.comp_c2l_model(bp.ARMS[0], 0).name == f"{bp.ARMS[0]}_s0"


def test_largest_common_n_is_bound_by_scarcest_type():
    vecs = {"a": {t: 1 / 6 for t in CT}}
    avail = {t: 1000 for t in CT}
    avail["B cells"] = 50
    # 50 B cells at share 1/6 -> N=300 gives exactly 50; N=301 gives 50 (50.17 -> 50 + remainder
    # ties go by type order, B first) so scan the exact boundary by hand:
    n = ca.largest_common_n(vecs, avail)
    assert ca.integer_counts(vecs["a"], n)["B cells"] <= 50
    assert ca.integer_counts(vecs["a"], n + 1)["B cells"] > 50 or n + 1 > 6 * 50 + 5
    assert 300 <= n <= 305


def test_integer_counts_sum_to_n_and_largest_remainder():
    v = {
        "B cells": 0.25,
        "T cells": 0.25,
        "Macrophages": 0.25,
        "Fibroblasts": 0.25,
        "Epithelial": 0.0,
        "Endothelial": 0.0,
    }
    c = ca.integer_counts(v, 10)
    assert sum(c.values()) == 10 and c["Epithelial"] == 0
    assert sorted(c[t] for t in CT[:4]) == [2, 2, 3, 3]
    w = dict(zip(CT, [0.3, 0.3, 0.2, 0.1, 0.05, 0.05]))
    c = ca.integer_counts(w, 7)
    assert sum(c.values()) == 7
    assert (c["B cells"], c["T cells"], c["Macrophages"], c["Fibroblasts"]) == (2, 2, 2, 1)


def test_all_arms_equal_n(option1):
    vecs = ca.arm_vectors(NATURAL, GT)
    avail = {t: 500 for t in CT}
    n = ca.largest_common_n(vecs, avail)
    for v in vecs.values():
        assert sum(ca.integer_counts(v, n).values()) == n


def _obs():
    rows = []
    rng = np.random.default_rng(0)
    for t in CT:
        for s, k in (("S1", 40), ("S2", 25), ("S3", 5)):
            rows += [(t, s)] * int(k + rng.integers(0, 3))
    df = pd.DataFrame(rows, columns=["cell_type", "sample"])
    df.index = [f"c{i}" for i in range(len(df))]
    return df


def test_stratified_sample_no_stratum_overdrawn_and_no_replacement():
    obs = _obs()
    counts = {t: 30 for t in CT}
    counts["B cells"] = int((obs.cell_type == "B cells").sum())  # take every B cell
    idx = ca.stratified_sample(obs, counts, 0)
    assert idx.is_unique and idx.isin(obs.index).all()
    sub = obs.loc[idx]
    for t in CT:
        assert (sub.cell_type == t).sum() == counts[t]
    avail = obs.groupby(["cell_type", "sample"]).size()
    got = sub.groupby(["cell_type", "sample"]).size()
    assert (got <= avail.reindex(got.index)).all()
    # proportional to stratum availability
    tm = obs[obs.cell_type == "T cells"]
    exp = tm["sample"].value_counts(normalize=True) * 30
    gt = sub[sub.cell_type == "T cells"]["sample"].value_counts()
    assert (gt.reindex(exp.index).fillna(0) - exp).abs().max() < 1


def test_stratified_sample_same_seed_identical_different_seed_differs():
    obs = _obs()
    counts = {t: 20 for t in CT}
    a = ca.stratified_sample(obs, counts, 0)
    b = ca.stratified_sample(obs, counts, 0)
    c = ca.stratified_sample(obs, counts, 1)
    assert list(a) == list(b)
    assert set(a) != set(c)


def test_stratified_sample_rejects_overdraw():
    obs = _obs()
    with pytest.raises(ValueError):
        ca.stratified_sample(obs, {"B cells": 10**6}, 0)


def test_map_reference_labels_collapses_t_subtypes_and_reports_unmapped():
    s = pd.Series(["CD4+ T cells", "CD8+ T cells", "T cells", "B cells", "Mystery"])
    mapped, unmapped = ca.map_reference_labels(s)
    assert list(mapped[:4]) == ["T cells", "T cells", "T cells", "B cells"]
    assert unmapped == {"Mystery": 1}


def test_manifest_refuses_overwrite(tmp_path, option1):
    ref = tmp_path / "ref.csv"
    ref.write_text("x\n1\n")
    gt = tmp_path / "gt.csv"
    gt.write_text("y\n2\n")
    vecs = ca.arm_vectors(NATURAL, GT)
    n = 60
    p = tmp_path / "arms_manifest.json"
    kw = dict(
        vectors=vecs,
        n=n,
        natural_counts=NATURAL,
        gt_overall=GT,
        reference_csv=ref,
        gt_labels_csv=gt,
        produced_by="Benchmarking/xenium_benchmarking/gt_rebuild/composition_arms.py",
    )
    ca.write_arms_manifest(p, **kw)
    import json

    m = json.loads(p.read_text())
    assert m["seeds"] == [0, 1, 2] and m["n"] == n
    assert m["arm_set"] == list(OPTION1)
    assert all(sum(c.values()) == n for c in m["integer_counts"].values())
    assert "_code_commit" in m
    with pytest.raises(FileExistsError):
        ca.write_arms_manifest(p, **kw)


def test_freeze_refuses_without_signoff(tmp_path, monkeypatch):
    from gt_rebuild import manifest

    def boom(*a, **k):
        raise manifest.SignoffMissing("no signoff")

    monkeypatch.setattr(manifest, "require_rnagt_signoff", boom)
    out = tmp_path / "arms_manifest.json"
    with pytest.raises(manifest.SignoffMissing):
        ca.main(["freeze", "--out", str(out)])
    assert not out.exists()


def _brute_largest_n(vectors, available, limit=3000):
    return max(
        n
        for n in range(1, limit)
        if all(ca.integer_counts(v, n)[t] <= available[t] for v in vectors.values() for t in CT)
    )


def test_largest_common_n_matches_brute_force_when_old_bound_was_too_low():
    # rare binding type (v=1/130): the old bound floor(avail/v)+6 was far below the true maximum
    v = dict(zip(CT, [x / 130 for x in (1, 10, 30, 20, 20, 49)]))
    avail = {t: 10000 for t in CT}
    avail["B cells"] = 2
    old_bound = min(int(avail[t] / v[t]) for t in CT if v[t] > 0) + len(CT)
    brute = _brute_largest_n({"a": v}, avail)
    assert brute > old_bound  # the case that exposed the unsound bound
    assert ca.largest_common_n({"a": v}, avail) == brute


def test_largest_common_n_matches_brute_force_small_grid():
    for w in ([2, 20, 20, 20, 19, 19], [3, 17, 20, 20, 20, 20], [1, 1, 1, 1, 1, 5]):
        v = {t: x / sum(w) for t, x in zip(CT, w)}
        for a in (1, 3, 7):
            avail = {t: 500 for t in CT}
            avail["B cells"] = a
            assert ca.largest_common_n({"a": v}, avail) == _brute_largest_n({"a": v}, avail)


def test_real_manifest_n_unchanged():
    import json

    m = json.loads((XB / "composition_arms" / "arms_manifest.json").read_text())
    assert ca.largest_common_n(m["vectors"], m["natural_counts"]) == m["n"] == 3912
