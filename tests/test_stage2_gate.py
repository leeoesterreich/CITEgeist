"""Stage 2 gate (spec: median across replicates per condition; Spearman(b, b*) >= 0.5 AND props r >= plain - 0.01)."""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "simulation_benchmarking"))
import stage2_gate as sg  # noqa: E402

pytestmark = pytest.mark.unit

PLAIN = 0.70


def _df(sp, r, plain=PLAIN, n=5, conds=("high_seg", "mixed")):
    rows = []
    for c in conds:
        for i in range(n):
            s = sp[i] if isinstance(sp, list) else sp
            q = r[i] if isinstance(r, list) else r
            rows.append({"condition": c, "replicate": i, "spearman_b_bstar": s, "props_r": q, "props_r_plain": plain})
    return pd.DataFrame(rows)


def test_gate_boundaries():
    assert sg.stage2_gate(_df(0.5, PLAIN))["outcome"] == "pass"
    assert sg.stage2_gate(_df(0.4999, PLAIN))["outcome"] == "drop"
    assert sg.stage2_gate(_df(0.9, PLAIN - 0.01))["outcome"] == "pass"
    assert sg.stage2_gate(_df(0.9, PLAIN - 0.0101))["outcome"] == "drop"


def test_gate_uses_median_across_replicates():
    # 2 of 5 replicates fail badly: median still passes; 3 of 5 fail: median fails
    assert sg.stage2_gate(_df([0.9, 0.9, 0.9, 0.0, 0.0], PLAIN))["outcome"] == "pass"
    assert sg.stage2_gate(_df([0.9, 0.9, 0.0, 0.0, 0.0], PLAIN))["outcome"] == "drop"
    assert sg.stage2_gate(_df(0.9, [0.8, 0.8, 0.8, 0.0, 0.0]))["outcome"] == "pass"
    assert sg.stage2_gate(_df(0.9, [0.8, 0.8, 0.0, 0.0, 0.0]))["outcome"] == "drop"


def test_failure_means_drop():
    df = _df(0.9, PLAIN)
    df.loc[df.condition == "mixed", "spearman_b_bstar"] = 0.1  # one condition fails => whole gate drops
    out = sg.stage2_gate(df)
    assert out["outcome"] == "drop"
    assert out["per_condition"]["high_seg"]["pass"] is True and out["per_condition"]["mixed"]["pass"] is False


def test_undefined_spearman_fails():
    df = _df(float("nan"), PLAIN)
    assert sg.stage2_gate(df)["outcome"] == "drop"


def test_empty_input_raises():
    with pytest.raises(ValueError):
        sg.stage2_gate(pd.DataFrame(columns=["condition", "replicate", "spearman_b_bstar", "props_r", "props_r_plain"]))


def test_missing_column_raises():
    with pytest.raises(ValueError):
        sg.stage2_gate(_df(0.9, PLAIN).drop(columns=["props_r_plain"]))


def test_plain_only_path_after_drop_emits_no_gate(tmp_path):
    rep = sg.build_report(None, chosen_method="")
    assert rep["plain_only"] is True
    assert "stage2_gate" not in rep and "outcome" not in rep
    assert "no chosen" in rep["note"].lower()
    p = tmp_path / "r.json"
    sg.write_report(p, rep)
    assert json.loads(p.read_text())["plain_only"] is True


def test_two_arm_report_embeds_gate():
    rep = sg.build_report(_df(0.9, PLAIN), chosen_method="rna_slope_adj")
    assert rep["plain_only"] is False and rep["stage2_gate"]["outcome"] == "pass"
    assert "recovery" in rep["note"].lower()  # recovery/stress test, not independent validation


def test_two_arm_without_table_raises():
    with pytest.raises(ValueError):
        sg.build_report(None, chosen_method="rna_slope_adj")


def test_nan_plain_fails():
    out = sg.stage2_gate(_df(0.9, 0.8, plain=float("nan")))
    assert out["outcome"] == "drop"
    assert any("plain" in r for r in out["per_condition"]["mixed"]["fail_reasons"])


def test_nan_candidate_fails():
    out = sg.stage2_gate(_df(0.9, float("nan")))
    assert out["outcome"] == "drop"
    assert any("candidate" in r for r in out["per_condition"]["mixed"]["fail_reasons"])


def test_both_nan_fails():
    assert sg.stage2_gate(_df(0.9, float("nan"), plain=float("nan")))["outcome"] == "drop"


def _write_arm(d, cond, rep, props, brightness):
    u = d / f"sim_6_{cond}_{rep}"
    (u / cond / f"Wu_rep_{rep}").mkdir(parents=True)
    (u / "diagnostics").mkdir()
    props.to_csv(u / cond / f"Wu_rep_{rep}" / "props.csv")
    (u / "diagnostics" / "baseline.json").write_text(json.dumps({"brightness": brightness}))


def test_per_rep_table_fixture(tmp_path):
    import numpy as np

    rng = np.random.default_rng(0)
    idx = [f"s{i}" for i in range(30)]
    gt = pd.DataFrame(rng.random((30, 3)), index=idx, columns=["A", "B", "C"])
    sim, chosen, plain = tmp_path / "sim", tmp_path / "chosen", tmp_path / "plain"
    markers = ["A_P1", "B_P1", "C_P1", "A_P2"]
    bstar = {"A_P1": 1.0, "B_P1": 2.0, "C_P1": 3.0, "A_P2": 4.0}
    (sim / "calsim_reports").mkdir(parents=True)
    for cond in ("high_seg", "mixed"):
        (sim / cond / "ST_sim").mkdir(parents=True)
        for rep in range(5):
            gt.to_csv(sim / cond / "ST_sim" / f"Wu_ST_{rep}_prop.csv")
            (sim / "calsim_reports" / f"bstar_{cond}_{rep}.json").write_text(
                json.dumps({"report": {"sims": {f"{cond}_{rep}": {"bstar": bstar}}}}))
            # chosen: perfectly ranked b via used_brightness (dict form) and plain numeric form for one marker set
            b = {m: {"used_brightness": bstar[m] * 10} for m in markers}
            _write_arm(chosen, cond, rep, gt, b)  # props == GT => r = 1
            _write_arm(plain, cond, rep, gt * 0.5 + 0.1, {})
    t = sg.per_rep_table(sim, chosen, plain)
    assert len(t) == 10 and list(t.columns) == sg.COLUMNS
    assert (t["spearman_b_bstar"] == 1.0).all()
    assert np.allclose(t["props_r"], 1.0) and np.allclose(t["props_r_plain"], 1.0)
    # reversed used_brightness => Spearman -1 (proves the key is actually read)
    for cond in ("high_seg", "mixed"):
        for rep in range(5):
            b = {m: {"used_brightness": -bstar[m]} for m in markers}
            (chosen / f"sim_6_{cond}_{rep}" / "diagnostics" / "baseline.json").write_text(json.dumps({"brightness": b}))
    assert (sg.per_rep_table(sim, chosen, plain)["spearman_b_bstar"] == -1.0).all()
