"""Task 18: gt_rebuild/score_methods.py (strict scoring of all methods, provenance-gated)."""

import hashlib
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
from gt_rebuild import score_methods as sm  # noqa: E402

pytestmark = pytest.mark.unit

A = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
SPOTS = [f"s{i}" for i in range(6)]
METHODS = ["CITEgeist", "RCTD", "Seurat", "CARD", "Cell2Location", "Tangram"]


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


@pytest.fixture
def world(tmp_path):
    gt = tmp_path / "gt"
    (gt / "ground_truth").mkdir(parents=True)
    rng = np.random.default_rng(0)
    for r in range(5):
        a = rng.uniform(0.1, 0.9, 6)
        df = pd.DataFrame({"A": a, "B": 1 - a}, index=SPOTS)
        df["n_cells"], df["spot_x"], df["spot_y"] = 5, 0.0, 0.0
        df.to_csv(gt / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_label="spot_id")
    preds, entries = {}, []
    for m in METHODS:
        d = tmp_path / m
        for r in range(5):
            gdf = pd.read_csv(gt / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_col=0)
            p = pd.DataFrame({"A": gdf["A"] + rng.normal(0, 0.05, 6), "B": gdf["B"]}, index=SPOTS)
            if m == "CITEgeist":
                f = d / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"
            else:
                f = d / f"Xenium_region_{r}" / f"Xenium_region_{r}_cell_prop_finetuned_results.csv"
            f.parent.mkdir(parents=True, exist_ok=True)
            p.to_csv(f)
            entries.append({"relpath": str(f.relative_to(tmp_path)), "sha256": _sha(f)})
        preds[m] = d
    inputs = tmp_path / "inputs.json"
    inputs.write_text(json.dumps({"inputs": {"groups": {"all": entries}}}))
    return {"tmp": tmp_path, "gt": gt, "preds": preds, "inputs": inputs, "out": tmp_path / "out"}


def _run(w, **kw):
    return sm.score(
        n_types=6,
        gt_dir=w["gt"],
        pred_dirs=w["preds"],
        out_dir=w["out"],
        inputs_json=w["inputs"],
        base=w["tmp"],
        types=["A", "B"],
        **kw,
    )


def test_score_methods_refuses_without_signoff(tmp_path, monkeypatch):
    (tmp_path / "manifest.json").write_text("{}")
    monkeypatch.setenv("XEN_LABELSET", "rnagt")
    monkeypatch.setattr(mf.gt_config, "OUT_ROOT", tmp_path)
    monkeypatch.setattr(sm, "score", lambda *a, **k: pytest.fail("scored without signoff"))
    with pytest.raises(mf.SignoffMissing):
        sm.main(6)
    assert not (tmp_path / "out").exists()


def test_main_refuses_outside_rnagt(monkeypatch):
    monkeypatch.delenv("XEN_LABELSET", raising=False)
    with pytest.raises(RuntimeError, match="rnagt"):
        sm.main(6)


def test_refuses_without_inputs_json(world):
    world["inputs"].unlink()
    with pytest.raises(FileNotFoundError, match="inputs.json"):
        _run(world)
    assert not world["out"].exists()


def test_refuses_on_sha_mismatch(world):
    f = world["preds"]["RCTD"] / "Xenium_region_2" / "Xenium_region_2_cell_prop_finetuned_results.csv"
    f.write_text(f.read_text() + "\n")  # content changed after hashing
    with pytest.raises(sm.ProvenanceError, match="sha256"):
        _run(world)
    assert not world["out"].exists()


def test_refuses_on_absent_entry(world):
    d = json.loads(world["inputs"].read_text())
    ents = d["inputs"]["groups"]["all"]
    d["inputs"]["groups"]["all"] = [e for e in ents if "Seurat" not in e["relpath"] or "region_3" not in e["relpath"]]
    world["inputs"].write_text(json.dumps(d))
    with pytest.raises(sm.ProvenanceError, match="no entry"):
        _run(world)
    assert not world["out"].exists()


def test_results_json_records_exact_baseline_path_and_sha256(world):
    _run(world)
    res = json.loads((world["out"] / "method_comparison_singler_6type" / "CITEgeist_6type_results.json").read_text())
    files = res["_inputs"]["pred_files"]
    assert len(files) == 5
    for r, e in enumerate(sorted(files, key=lambda e: e["path"])):
        p = world["preds"]["CITEgeist"] / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"
        assert e["path"] == str(p)
        assert e["sha256"] == _sha(p)
    assert len(res["_inputs"]["evaluable_table_sha256"]) == 64
    assert res["summary"]["n_types_scored"] == 2


def test_outputs_named_and_complete(world):
    _run(world)
    d = world["out"] / "method_comparison_singler_6type"
    names = {p.name for p in d.iterdir()}
    assert names == {
        "CITEgeist_6type_results.json",
        "RCTD_6type_results.json",
        "Seurat_6type_results.json",
        "CARD_results.json",
        "Cell2Location_6type_results.json",
        "Tangram_6type_results.json",
        "region_type_table_6type.csv",
        "mask_report_6type.json",
        "evaluability_6type.json",
    }
    tab = pd.read_csv(d / "region_type_table_6type.csv")
    assert set(tab["method"]) == set(METHODS) and len(tab) == 6 * 5 * 2
    assert "reason" in tab.columns
    rep = json.loads((d / "mask_report_6type.json").read_text())
    assert rep["CITEgeist"]["0"]["spots_dropped"] == 0


def test_n_types_scored_mismatch_raises_naming_type(world, monkeypatch):
    gt = world["gt"] / "ground_truth"
    for r in range(5):  # B constant in every region -> never scorable
        df = pd.read_csv(gt / f"Xenium_region_{r}_prop.csv", index_col=0)
        df["B"] = 0.5
        df.to_csv(gt / f"Xenium_region_{r}_prop.csv", index_label="spot_id")
    with pytest.raises(RuntimeError, match="'B'"):
        _run(world)


def test_types_exclude_metadata_columns(world):
    assert sm.gt_types(world["gt"], 2) == ["A", "B"]
