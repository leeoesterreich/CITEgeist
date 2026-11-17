"""A12.3: synthetic-gate machinery (small runs) and fail-closed verdict functions with mutations."""
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

import gate_verdicts as V
import spatial_utils as S
import synthetic_g6 as G6
import synthetic_g7 as G7

pytestmark = pytest.mark.unit

REMEDIATION = Path(__file__).resolve().parents[2] / "Benchmarking" / "midkine_remediation"
HAVE_DATA = (S.DATA_ROOT / "HCC22-088-P4-S2_1i_rep" / "outs" / "spatial" / "tissue_positions.csv").exists()
needs_data = pytest.mark.skipif(not HAVE_DATA, reason="Visium coordinates not on this filesystem")


@pytest.fixture(scope="module")
def lat():
    return G6.make_lattice(bws=(4.0, 6.0, 10.0))


def test_wilson_interval():
    lo, hi = G6.wilson(50, 100)
    assert lo == pytest.approx(0.4038, abs=1e-3) and hi == pytest.approx(0.5962, abs=1e-3)
    assert G6.wilson(0, 10)[0] == 0.0
    assert np.isnan(G6.wilson(0, 0)[0])


@needs_data
def test_lattice_sizes_match_protocol_qp_recount():
    proto = json.loads((REMEDIATION / "protocol.DRAFT.json").read_text())
    lit = proto["specimens"]["expected_cancer_compartments_stored_prop"]["ge_0.10"]
    for spec in ("HCC22-088-P4-S2_1i_rep", "HCC22-088-P1-S2"):
        xy, q = S.load_lattice(spec, 0.10, qp_dir=S.REPO_ROOT / proto["inputs"]["synthetic_lattice_dir"])
        assert len(xy) == lit[spec] == len(q)             # the lattice is the independent QP recount, not a SACE object


@needs_data
def test_simulate_properties(lat):
    sim = G6.simulate("null_autocorrelated", lat, np.random.default_rng(0))
    n = len(lat.xy)
    assert sim["group"].shape == (n,) and abs(sim["group"].mean() - G6.HIGH_FRACTION) < 0.01
    sc = sim["scores"][0]
    assert (sc == 0).mean() > 0.05 and (sc > 0).mean() > 0.3          # zero-inflated hurdle outcome
    # labels are spatially clustered: neighbours share the label far more often than chance
    W = S.knn_weights(lat.xy, 6)
    same = float((W @ sim["group"] * sim["group"]).sum() / sim["group"].sum())
    assert same > 0.7 > G6.HIGH_FRACTION
    # injection multiplies only the high group's scores
    inj = G6.simulate("null_autocorrelated", lat, np.random.default_rng(0), inject_fold=3.0)
    ratio = inj["scores"][0][sim["group"] == 1] / np.where(sc[sim["group"] == 1] == 0, np.nan, sc[sim["group"] == 1])
    assert np.nanmax(np.abs(ratio - 3.0)) < 1e-9
    assert (inj["scores"][0][sim["group"] == 0] == sc[sim["group"] == 0]).all()


@needs_data
def test_type1_and_injection_run_end_to_end_small(lat):
    r = G6.type1("null_autocorrelated", n_datasets=4, B=49, lat=lat)
    for part in ("magnitude", "participation"):
        assert r[part]["n_ok"] + r[part]["n_not_estimable"] == 4 and 0 <= r[part]["type1"] <= 1
    inj = G6.injection(n_draws=2, B=49, lat=lat)
    assert set(inj) >= {"recovery_rate", "gmr_mean", "frac_q_lt_0.05"}


def _fake_g6(tmp, prefix, rates, n_ok=1000):
    for c, (m, p) in rates.items():
        payload = {"run": "x", "result": {"n_datasets": 1000, "magnitude": {"n_ok": n_ok, "type1": m}, "participation": {"n_ok": n_ok, "type1": p}}}
        import provenance

        provenance.write_stamped_json(tmp / f"g6_{prefix}_{c}.json", payload, "test")


def test_g6_verdict_fail_closed_and_mutations(tmp_path):
    good = {c: (0.05, 0.06) for c in V.G6_NULL_CASES}
    _fake_g6(tmp_path, "t1", good)
    v = V.g6_type1("t1", synth=tmp_path)
    assert all(x["pass"] for x in v.values()) and len(v) == 3
    # one case inflated -> that case fails, others still pass
    _fake_g6(tmp_path, "t1", {"null_noisy_nst": (0.05, 0.12)})
    v = V.g6_type1("t1", synth=tmp_path)
    assert not v["null_noisy_nst"]["pass"] and v["null_autocorrelated"]["pass"]
    # delete the LAST required case file -> raises (does not silently pass)
    (tmp_path / f"g6_t1_{V.G6_NULL_CASES[-1]}.json").unlink()
    with pytest.raises(FileNotFoundError):
        V.g6_type1("t1", synth=tmp_path)
    # matcher broken: no cases -> raises
    with pytest.raises(ValueError):
        V.g6_type1("t1", cases=(), synth=tmp_path)
    # no estimable datasets -> raises
    _fake_g6(tmp_path, "z", {c: (0.0, 0.0) for c in V.G6_NULL_CASES}, n_ok=0)
    with pytest.raises(ValueError):
        V.g6_type1("z", synth=tmp_path)
    # tampered payload -> stamp verification fails
    p = tmp_path / "g6_t1_null_autocorrelated.json"
    d = json.loads(p.read_text())
    d["result"]["magnitude"]["type1"] = 0.0
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError):
        V.g6_type1("t1", cases=("null_autocorrelated",), synth=tmp_path)


def test_g7_summarize_fail_closed_and_detects_inflation():
    def row(rho, msr, naive):
        base = {"type1_two_sided": 0.0, "ci95": [0, 0], "type1_upper": 0.0}
        d = {"lattice": "L", "n": 10, "rho": rho, "x": "continuous", "n_pairs": 1000}
        for k in ("msr_raw", "msr_rank", "naive_raw", "naive_rank"):
            d[k] = dict(base)
        d["msr_raw"]["type1_two_sided"], d["naive_raw"]["type1_two_sided"] = msr, naive
        return d

    ok = {"settings": [row(0.5, 0.05, 0.1), row(0.95, 0.06, 0.6)]}
    s = G7.summarize(ok)
    assert s["msr_gate_pass"] and s["naive_inflated_at_high_rho"]
    bad = {"settings": ok["settings"] + [row(0.8, 0.12, 0.3)]}
    assert not G7.summarize(bad)["msr_gate_pass"]
    # mutation: naive null NOT inflated at high rho -> the gate's own failure check reports it
    flat = {"settings": [row(0.5, 0.05, 0.05), row(0.95, 0.05, 0.05)]}
    assert not G7.summarize(flat)["naive_inflated_at_high_rho"]
    with pytest.raises(ValueError):
        G7.summarize({"settings": []})                 # empty scan never passes
