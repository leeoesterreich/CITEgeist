# tests/benchmarking/test_protein_conditions.py
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import perturb as P  # noqa: E402
from protein_sensitivity.conditions import build  # noqa: E402

pytestmark = pytest.mark.unit

PROTO = json.loads((Path(__file__).resolve().parents[2] / "Benchmarking/protein_sensitivity/protocol.json").read_text())


def _expected(profile, key):
    """Independent (axis, level, seed, panel, sentinel) set from protocol.json + perturb."""
    sent = {}
    for tag, items in PROTO["sentinels"].items():
        for it in items:
            panel = it.get("panel")
            if isinstance(panel, dict):
                panel = panel[key]
            sent[(it["axis"], it.get("level"), panel)] = tag
    exp = {("baseline", None, None, None, "baseline")}
    for ax in ("thinning", "dropout", "noise"):
        for lvl in PROTO["axes"][ax]:
            for seed in PROTO["seeds"]:
                exp.add((ax, float(lvl), seed, None, sent.get((ax, lvl, None), "")))
    for ax, panels in (
        ("ablation", P.ablation_panels(profile) + P.leave_one_type_out(profile)),
        ("type_removal", P.type_removal_series(profile)),
    ):
        for name, _ in panels:
            exp.add((ax, None, None, name, sent.get((ax, None, name), "")))
    return exp


def _got(conds):
    return {(c["axis"], c.get("level"), c.get("seed"), c.get("panel"), c["sentinel"]) for c in conds}


@pytest.mark.parametrize("dataset,n,key", [("xenium", 6, "xenium6"), ("xenium", 7, "xenium7"), ("sim", 6, "sim")])
def test_sweep_matches_independent_expectation(dataset, n, key):
    profile = B.sim_profile() if dataset == "sim" else B.PROFILE[n]
    conds = json.loads(build("sweep", dataset, n))
    assert _got(conds) == _expected(profile, key)
    assert len(conds) == len(_got(conds))  # no duplicates


def test_xenium6_counts():
    profile = B.PROFILE[6]
    conds = json.loads(build("sweep", "xenium", 6))
    stoch = [c for c in conds if c["axis"] in ("thinning", "dropout", "noise")]
    assert len(stoch) == (4 + 3 + 4) * 3 == 33
    det = [c for c in conds if c["axis"] in ("ablation", "type_removal")]
    assert len(det) == len(P.ablation_panels(profile)) + len(P.leave_one_type_out(profile)) + len(
        P.type_removal_series(profile)
    )
    assert sum(c["axis"] == "baseline" for c in conds) == 1
    # sentinels exist (22 per unit) but SACE is not run on Xenium 6-type
    assert sum(c["sentinel"] != "" for c in conds) == 22
    assert sum(c["sace"] for c in conds) == 0


@pytest.mark.parametrize("dataset,n", [("xenium", 7), ("sim", 6)])
def test_sace_counts(dataset, n):
    conds = json.loads(build("sweep", dataset, n))
    assert sum(c["sace"] for c in conds) == 22
    assert all(c["sace"] == (c["sentinel"] != "") for c in conds)


def test_severe_panel_is_last_nonzero_series_entry():
    for key, profile in (("xenium6", B.PROFILE[6]), ("xenium7", B.PROFILE[7]), ("sim", B.sim_profile())):
        names = [n for n, _ in P.type_removal_series(profile)]
        assert PROTO["sentinels"]["severe"][3]["panel"][key] == names[-2]
        assert names[-1] == "zero"


def test_gate_sace_flags():
    assert json.loads(build("gate", "xenium", 6)) == [{"axis": "baseline", "sentinel": "baseline", "sace": False}]
    assert json.loads(build("gate", "xenium", 7))[0]["sace"] is True
    assert json.loads(build("gate", "sim", 6))[0]["sace"] is True
