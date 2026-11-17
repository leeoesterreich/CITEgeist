"""SPOTS regeneration for the v41p rerun: saved predictions, shared perturbation, strict zero-marker mode."""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import run_sensitivity as RS  # noqa: E402
from protein_sensitivity import spots_sensitivity as S  # noqa: E402
from protein_sensitivity import strict_score as SS  # noqa: E402
from protein_sensitivity import variant as V  # noqa: E402

pytestmark = pytest.mark.unit

PROFILE = {"A": {"Major": ["A1", "A2"]}, "B": {"Major": ["B1"]}, "C": {"Major": ["C1"]}}
VARS = ["A1", "A2", "B1", "C1", "X1"]
ZERO = {"axis": "ablation", "panel": "zero", "profile": {}, "sentinel": ""}
THIN = {"axis": "thinning", "level": 0.5, "seed": 1, "sentinel": ""}


@pytest.fixture
def cite(monkeypatch):
    monkeypatch.setitem(B.PROFILE, 6, PROFILE)
    rng = np.random.default_rng(42)
    X = rng.integers(1, 60, size=(12, 5)).astype(float)
    return ad.AnnData(X=X, obs=pd.DataFrame(index=[f"s{i}" for i in range(12)]), var=pd.DataFrame(index=VARS))


def _metrics(pred, u):
    return {"pearson_r": 0.5, "rmse": 0.1, "mae": 0.1, "jsd": 0.1, "per_type": {}}


def test_zero_markers_fail_mode_gives_failed_no_markers_and_no_prediction(cite):
    saved, seen = [], []

    def metrics(pred, u):
        seen.append(pred)
        return _metrics(pred, u)

    rows = S.run_unit(
        6, "0", [ZERO], cite=cite.copy(), metrics_fn=metrics, save_pred=lambda c, p: saved.append(c), zero_markers="fail"
    )
    assert rows[0]["status"] == "failed_no_markers"
    assert np.isnan(rows[0]["pearson_r"])
    assert saved == [] and seen == []


def test_legacy_default_still_fills_uniform(cite):
    store = []

    def metrics(pred, u):
        store.append(pred)
        return _metrics(pred, u)

    rows = S.run_unit(6, "0", [ZERO], cite=cite.copy(), metrics_fn=metrics)
    assert rows[0]["status"] == "ok"
    assert np.allclose(store[0].to_numpy(), 1 / 3)


def test_unknown_zero_markers_mode_is_refused(cite):
    with pytest.raises(ValueError):
        S.run_unit(6, "0", [ZERO], cite=cite.copy(), metrics_fn=_metrics, zero_markers="bogus")


def test_spots_saves_one_prediction_per_ok_condition(cite):
    saved = {}
    conds = [THIN, {**THIN, "seed": 2}, ZERO]
    rows = S.run_unit(
        6,
        "0",
        conds,
        cite=cite.copy(),
        metrics_fn=_metrics,
        save_pred=lambda c, p: saved.__setitem__(c, p),
        zero_markers="fail",
    )
    assert [r["status"] for r in rows] == ["ok", "ok", "failed_no_markers"]
    assert sorted(saved) == ["thinning_0.5_s1", "thinning_0.5_s2"]


def _perturbed_hashes(cite, qp_fn):
    hashes = {}
    conds = [THIN, {**THIN, "seed": 2}]
    RS.evaluate_unit(
        "xenium",
        "0",
        6,
        conds,
        loader=lambda u: (S._NoGex(), cite.copy()),
        qp_fn=qp_fn,
        sace_fn=None,
        prop_metrics_fn=_metrics,
        gex_metrics_fn=None,
        profile=PROFILE,
        on_perturbed=lambda cid, c2: hashes.__setitem__(cid, V.perturbed_sha256(c2)),
    )
    return hashes


def test_spots_and_citegeist_wiring_hash_same_perturbation(cite):
    def other_qp(g, c2, prof, u, *extra):
        return pd.DataFrame(1.0, index=c2.obs_names, columns=list(prof))

    def spots_qp(g, c2, prof, u, *extra):
        df, _ = S.spots_predict(np.asarray(c2.X), c2.var_names, prof, list(prof))
        df.index = c2.obs_names
        return df

    a, b = _perturbed_hashes(cite, other_qp), _perturbed_hashes(cite, spots_qp)
    assert a == b and len(a) == 2
    assert a["thinning_0.5_s1"] != a["thinning_0.5_s2"]


def _write_unit(d, hashes, statuses=None):
    d.mkdir(parents=True)
    with open(d / "row_meta.jsonl", "w") as fh:
        for cid, h in hashes.items():
            row = {"condition_id": cid, "perturbed_sha256": h}
            if statuses:
                row["status"] = statuses[cid]
            fh.write(json.dumps(row) + "\n")
    return d


def test_check_shared_perturbation_detects_one_changed_hash(tmp_path):
    h = {"baseline": "a", "thin_s0": "b", "thin_s1": "c"}
    a = _write_unit(tmp_path / "a", h)
    assert SS.check_shared_perturbation(a, _write_unit(tmp_path / "b", dict(h))) == []
    last = list(h)[-1]  # mutate the LAST row
    bad = SS.check_shared_perturbation(a, _write_unit(tmp_path / "c", {**h, last: "zzz"}))
    assert len(bad) == 1 and last in bad[0]


def test_null_hash_on_ok_row_is_an_error_even_when_both_null(tmp_path):
    h = {"baseline": "a", "thin_s0": "b", "thin_s1": "c"}
    st = {k: "ok" for k in h}
    a = _write_unit(tmp_path / "a", h, st)
    last = list(h)[-1]  # LAST row null on both sides: None == None must not pass
    bad = SS.check_shared_perturbation(
        _write_unit(tmp_path / "n1", {**h, last: None}, st), _write_unit(tmp_path / "n2", {**h, last: None}, st)
    )
    assert len(bad) == 2 and all(last in b and "null" in b for b in bad)
    assert SS.check_shared_perturbation(a, _write_unit(tmp_path / "ok", dict(h), st)) == []


def test_null_hash_allowed_on_failed_no_markers_row(tmp_path):
    h = {"baseline": "a", "zero": None}
    st = {"baseline": "ok", "zero": "failed_no_markers"}
    assert SS.check_shared_perturbation(_write_unit(tmp_path / "a", h, st), _write_unit(tmp_path / "b", dict(h), st)) == []
    bad = SS.check_shared_perturbation(_write_unit(tmp_path / "c", h, st), _write_unit(tmp_path / "d", {**h, "zero": "x"}, st))
    assert len(bad) == 1


def test_check_shared_perturbation_raises_on_empty_or_mismatched_ids(tmp_path):
    h = {"baseline": "a", "thin_s0": "b"}
    a = _write_unit(tmp_path / "a", h)
    with pytest.raises(ValueError):
        SS.check_shared_perturbation(a, _write_unit(tmp_path / "empty", {}))
    with pytest.raises(ValueError):
        SS.check_shared_perturbation(_write_unit(tmp_path / "empty2", {}), a)
    with pytest.raises(ValueError):
        SS.check_shared_perturbation(a, _write_unit(tmp_path / "other", {"baseline": "a", "different": "b"}))


def test_spots_variant_profile_has_cd138():
    v41p = V.RunVariant(profile="ms_v41", brightness=False)
    new, old = V.profile_for(6, v41p), V.profile_for(6, None)
    names = sorted({m for p in new.values() for k in ("Major",) for m in p[k]} | {m for p in old.values() for m in p["Major"]})
    assert "CD138" in names
    rng = np.random.default_rng(0)
    X = rng.integers(1, 50, size=(8, len(names))).astype(float)
    X2 = X.copy()
    X2[:, names.index("CD138")] *= 7  # change ONLY the CD138 channel
    col = "B cells"
    n1, _ = S.spots_predict(X, names, new, list(new))
    n2, _ = S.spots_predict(X2, names, new, list(new))
    o1, _ = S.spots_predict(X, names, old, list(old))
    o2, _ = S.spots_predict(X2, names, old, list(old))
    assert not np.allclose(n1[col], n2[col])
    assert np.allclose(o1[col], o2[col])


def test_spots_cli_writes_row_meta_and_run_variant(cite, tmp_path, monkeypatch):
    monkeypatch.setattr(V, "profile_for", lambda n, variant: PROFILE)
    cfile = tmp_path / "conds.json"
    cfile.write_text(json.dumps([{"axis": "baseline", "sentinel": ""}, THIN, ZERO]))
    out = tmp_path / "unit"
    S.run_unit_cli(
        [str(out), "6", "0", "--conditions", str(cfile), "--profile", "ms_v41"], cite=cite.copy(), metrics_fn=_metrics
    )
    meta = json.loads((out / "run_meta.json").read_text())
    assert meta["run_variant"]["profile"] == "ms_v41" and meta["run_variant"]["brightness"] is False
    rm = [json.loads(x) for x in (out / "row_meta.jsonl").read_text().splitlines()]
    assert [r["condition_id"] for r in rm] == ["baseline", "thinning_0.5_s1", "ablation_zero"]
    assert rm[1]["perturbed_sha256"] and rm[2]["perturbed_sha256"] is None and rm[2]["status"] == "failed_no_markers"
    assert sorted(p.name for p in (out / "predictions").iterdir()) == ["baseline.csv.gz", "thinning_0.5_s1.csv.gz"]
    m = pd.read_csv(out / "metrics.csv")
    assert list(m.columns) == RS.ROW_KEYS and list(m["status"]) == ["ok", "ok", "failed_no_markers"]
