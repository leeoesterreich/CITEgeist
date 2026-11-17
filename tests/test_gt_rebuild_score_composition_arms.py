"""Task 25: gt_rebuild/score_composition_arms.py (strict scoring of the composition-sensitivity arms)."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
XEN = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XEN))
sys.path.insert(0, str(XEN / "evaluation" / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from gt_rebuild import manifest as mf  # noqa: E402
from gt_rebuild import score_composition_arms as sca  # noqa: E402

pytestmark = pytest.mark.unit

SPOTS = [f"s{i}" for i in range(200)]
TYPES = ["A", "B"]
ARMS = ("a1", "a2")
SEEDS = [0, 1]
METHODS = ["M1", "M2"]


def _write_pred(d, r, df):
    f = Path(d) / f"Xenium_region_{r}" / f"Xenium_region_{r}_cell_prop_finetuned_results.csv"
    f.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(f)


@pytest.fixture
def world(tmp_path):
    rng = np.random.default_rng(1)
    gt = tmp_path / "gt"
    (gt / "ground_truth").mkdir(parents=True)
    gts = {}
    for r in range(5):
        a = rng.uniform(0.1, 0.9, len(SPOTS))
        df = pd.DataFrame({"A": a, "B": 1 - a}, index=SPOTS)
        gts[r] = df.copy()
        df["n_cells"], df["spot_x"], df["spot_y"] = 5, 0.0, 0.0
        df.to_csv(gt / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_label="spot_id")
    full = tmp_path / "full"
    for r in range(5):
        _write_pred(
            full, r, pd.DataFrame({"A": gts[r]["A"] + rng.normal(0, 0.05, len(SPOTS)), "B": gts[r]["B"]}, index=SPOTS)
        )
    runs = tmp_path / "runs"
    for m in METHODS:
        for arm in ARMS:
            for s in SEEDS:
                for r in range(5):
                    p = pd.DataFrame({"A": gts[r]["A"] + rng.normal(0, 0.1, len(SPOTS)), "B": gts[r]["B"]}, index=SPOTS)
                    _write_pred(runs / m / f"{arm}_s{s}", r, p)
    return {"tmp": tmp_path, "gt": gt, "full": full, "runs": runs}


def _core(w, **kw):
    args = dict(
        gt_dir=w["gt"],
        mask_pred_dirs={"full": w["full"]},
        run_dir=lambda m, arm, s: w["runs"] / m / f"{arm}_s{s}",
        arm_set=ARMS,
        seeds=SEEDS,
        methods=METHODS,
        types=TYPES,
    )
    args.update(kw)
    return sca.score_core(**args)


def test_score_arms_refuses_without_signoff(tmp_path, monkeypatch):
    (tmp_path / "manifest.json").write_text("{}")
    monkeypatch.setattr(mf.gt_config, "OUT_ROOT", tmp_path)
    monkeypatch.setattr(sca, "score_core", lambda *a, **k: pytest.fail("scored without signoff"))
    with pytest.raises(mf.SignoffMissing):
        sca.score_arms()


def test_expected_rows_derived_from_arm_set(world):
    df, shas, _ = _core(world)
    assert len(df) == len(ARMS) * len(SEEDS) * len(METHODS)
    assert set(df["arm"]) == set(ARMS) and set(df["seed"]) == set(SEEDS)
    assert len(shas) == len(df) * 5  # one sha per prediction file read
    with pytest.raises(ValueError, match="expected"):
        sca.check_row_count(df.iloc[:-1], ARMS, SEEDS, METHODS)


def test_mask_not_silently_shrunk_when_arm_lacks_spots(world):
    f = world["runs"] / "M1" / "a2_s1" / "Xenium_region_2" / "Xenium_region_2_cell_prop_finetuned_results.csv"
    pd.read_csv(f, index_col=0).iloc[:-1].to_csv(f)
    with pytest.raises(ValueError, match="absent"):
        _core(world)


def test_arm_scores_use_full_size_mask(world):
    f = world["full"] / "Xenium_region_0" / "Xenium_region_0_cell_prop_finetuned_results.csv"
    pd.read_csv(f, index_col=0).iloc[:-1].to_csv(f)  # full-size run has 199/200 spots
    with pytest.raises(sca.MaskMismatch):
        _core(world, expected_mask_spots={0: 200, 1: 200, 2: 200, 3: 200, 4: 200})
    df, _, mr = _core(world, expected_mask_spots={0: 199, 1: 200, 2: 200, 3: 200, 4: 200})
    assert mr[0] == 199


def test_summary_has_spread_and_no_pvalue_or_slope(world):
    df, _, _ = _core(world)
    s = sca.summarise(df)
    assert len(s) == len(METHODS) * len(ARMS)
    for c in ("per_type_mean_r_mean", "per_type_mean_r_min", "per_type_mean_r_max"):
        assert c in s.columns
    assert (s["per_type_mean_r_min"] <= s["per_type_mean_r_mean"]).all()
    assert (s["per_type_mean_r_mean"] <= s["per_type_mean_r_max"]).all()
    for c in s.columns:
        assert "p_" not in c and "slope" not in c and "sem" not in c, c


def test_per_type_r_columns_present(world):
    df, _, _ = _core(world)
    assert {"r_A", "r_B", "per_type_mean_r", "flat_r", "mae", "jsd"} <= set(df.columns)
    assert np.allclose(df["per_type_mean_r"], df[["r_A", "r_B"]].mean(axis=1))


def test_full_size_reported_as_endpoint_not_arm(tmp_path):
    res = tmp_path / "res"
    d = res / "method_comparison_singler_6type"
    d.mkdir(parents=True)
    summ = {"per_type_mean_r": 0.7, "overall_mean_pearson_r": 0.5, "overall_mean_mae": 0.1, "overall_mean_jsd": 0.3}
    for m, f in {"RCTD": "RCTD_6type_results.json", "CARD": "CARD_results.json"}.items():
        (d / f).write_text(json.dumps({"summary": dict(summ, **{"A_mean_r": 0.6})}))
    for f in ("CITEgeist_6type_results.json", "CITEgeist_v41p_6type_results.json"):
        (d / f).write_text(json.dumps({"summary": dict(summ, **{"A_mean_r": 0.6})}))
    ref = sca.reference_rows(res, ["RCTD", "CARD"], ["A"])
    kinds = {(r["method"], r["kind"]) for r in ref}
    assert ("RCTD", "full_size_reference") in kinds and ("CARD", "full_size_reference") in kinds
    assert ("CITEgeist", "reference_line") in kinds and ("CITEgeist_v41p", "reference_line") in kinds
    assert all("arm" not in r for r in ref)
    assert ref[0]["sha256"]


def test_write_refuses_overwrite_and_stamps(world, tmp_path):
    df, shas, _ = _core(world)
    out = tmp_path / "out"
    out.mkdir()
    payload = {"arms": df.to_dict(orient="records"), "note": "x"}
    sca.write_outputs(out, df, payload)
    assert (out / "composition_sensitivity.csv").is_file()
    assert (out / "composition_sensitivity.csv.provenance.json").is_file()
    stamped = json.loads((out / "composition_sensitivity.json").read_text())
    assert "_payload_sha256" in stamped
    with pytest.raises(FileExistsError):
        sca.write_outputs(out, df, payload)


def test_note_about_gt_matched_fibroblasts():
    assert "495" in sca.GT_MATCHED_NOTE and "lower bound" in sca.GT_MATCHED_NOTE
