"""A12.3 G10 (protocol rev 6, descriptive regime, P4-S2 headline scope): the H1 claim tree is exhaustive and disjoint over all 48 combinations (mutation-verified)."""
import json
from pathlib import Path

import pytest

import claim_tree as C

pytestmark = pytest.mark.unit

PROTO = Path(__file__).resolve().parents[2] / "Benchmarking/midkine_remediation/protocol.DRAFT.json"


def test_enumeration_size_is_48_and_distinct():
    allc = list(C.enumerate_inputs())
    assert len(allc) == 48 == 2**4 * 3
    assert len(set(allc)) == 48
    assert {i.d for i in allc} == set(C.D_VALUES)


def test_inputs_contain_no_p_q_interval_or_control_dependency():
    # rev 6: no hurdle p/q-value, interval property (b / has_ci), BH q, below_gate (fixed True in scope) or controls gate is an input
    assert C.Inputs._fields == ("est", "s", "c", "d", "f")
    assert C.SCOPE_SPECIMEN == "HCC22-088-P4-S2_1i_rep"


def test_exactly_one_condition_for_every_combination():
    n_by_branch = {b: 0 for b in C.BRANCH_IDS}
    for i in C.enumerate_inputs():
        hits = C.evaluate(i)
        assert len(hits) == 1, (i, hits)
        n_by_branch[hits[0]] += 1
    assert all(v > 0 for v in n_by_branch.values())  # no dead branch
    assert sum(n_by_branch.values()) == 48
    # independent hand counts (not derived from the tree)
    assert n_by_branch == {
        "not_estimable": 24,                 # est False: 2^3 (s, c, f) * 3 (d)
        "not_higher": 12,                    # est, not s: 2^2 (c, f) * 3
        "direction_not_robust": 9,           # est, s, not (c and f): 3 of the 4 (c, f) combinations * 3
        "positive_between_spot_agrees": 1,   # est, s, c, f, d == agree
        "positive_native_only": 2,           # est, s, c, f, d in {disagree, not_estimable}
    }


def test_condition_ids_and_inputs_match_protocol():
    proto = json.loads(PROTO.read_text())
    block = proto["decision_rules"]["hurdle_branches_for_H1_per_part"]
    assert tuple(c["id"] for c in block["conditions_disjoint_and_exhaustive"]) == C.BRANCH_IDS
    assert set(block["inputs"]) == set(C.Inputs._fields) | {"robust"}
    assert block["scope"]["specimen"] == C.SCOPE_SPECIMEN
    assert "48" in block["exhaustiveness_test"] and "384" not in json.dumps(proto["decision_rules"]) and "1536" not in json.dumps(proto["decision_rules"])


def _mutate(name, fn):
    conds = dict(C.CONDITIONS)
    conds[name] = fn
    return conds


def test_weakening_one_condition_makes_overlap():
    # drop `s` from direction_not_robust -> overlaps not_higher
    m = _mutate("direction_not_robust", lambda i: i.est and not (i.c and i.f))
    r = C.check_exhaustive(m)
    assert r["multi"] > 0 and r["exactly_one"] < 48


def test_overlapping_two_conditions_is_caught():
    # positive_native_only ignores d -> overlaps positive_between_spot_agrees
    m = _mutate("positive_native_only", lambda i: i.est and i.s and i.c and i.f)
    assert C.check_exhaustive(m)["multi"] > 0
    with pytest.raises(AssertionError):
        for i in C.enumerate_inputs():
            C.decide(i, m)


def test_dropping_f_from_a_positive_branch_is_caught():
    m = _mutate("positive_between_spot_agrees", lambda i: i.est and i.s and i.c and i.d == "agree")
    r = C.check_exhaustive(m)
    assert r["multi"] > 0  # c and not f now satisfies both this branch and direction_not_robust


def test_gap_is_caught():
    m = _mutate("not_higher", lambda i: False)
    r = C.check_exhaustive(m)
    assert r["zero"] > 0


def test_unmutated_check_reports_clean():
    assert C.check_exhaustive() == {"n": 48, "zero": 0, "multi": 0, "exactly_one": 48}


def test_hand_cases_and_reasons():
    base = dict(est=True, s=True, c=True, d="agree", f=True)
    I = lambda **k: C.Inputs(**{**base, **k})  # noqa: E731
    assert C.decide(I())["branch"] == "positive_between_spot_agrees"
    assert C.decide(I(d="disagree"))["branch"] == "positive_native_only"
    assert C.decide(I(d="not_estimable"))["branch"] == "positive_native_only"
    assert C.decide(I(est=False))["branch"] == "not_estimable"
    assert C.decide(I(s=False))["branch"] == "not_higher"
    # not_higher does not depend on robustness inputs
    assert C.decide(I(s=False, c=False, f=False))["branch"] == "not_higher"
    # nuclei-free sign disagrees (P4-S2 is below gate) -> direction not robust, reason segmentation-dependent
    assert C.decide(I(f=False)) == {"branch": "direction_not_robust", "reason": "segmentation-dependent"}
    # normalisation sign disagreement (c False) -> other reason, also when f fails
    assert C.decide(I(c=False)) == {"branch": "direction_not_robust", "reason": "not robust to normalisation"}
    assert C.decide(I(c=False, f=False))["reason"] == "not robust to normalisation"
