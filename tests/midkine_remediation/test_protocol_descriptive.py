"""Protocol rev 6 consistency: the hurdle sections are internally consistent with DESCRIPTIVE-ONLY inference."""
import json
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

PROTO = Path(__file__).resolve().parents[2] / "Benchmarking/midkine_remediation/protocol.DRAFT.json"
# keys that may DESCRIBE superseded permutation/BH machinery or negate it
HISTORICAL = {"superseded_permutation_tests", "inference_status", "g6_fix_round_result", "g6_gate_status", "open_limitation",
              "case_null_autocorrelated", "multiple_testing", "regime", "required_statements_with_the_headline",
              "justification", "mutation_checks_note", "pre_freeze_deviations", "H2", "H3", "strongest_pair_phrase", "wording",
              "negative_control", "status", "estimands", "bootstrap", "hurdle_only_score_injection", "gates"}
PVAL = re.compile(r"\b[pq]\s*[<=]\s*0?\.\d|block[- ]permutation|\bBH\b|Benjamini|\bq-value|smallest (attainable|one-sided|possible) p|controls_pass", re.I)


@pytest.fixture(scope="module")
def proto():
    return json.loads(PROTO.read_text())


def _walk(o, path=()):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from _walk(v, path + (k,))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from _walk(v, path + (str(i),))
    else:
        yield path, str(o)


def test_no_live_p_or_q_dependency_in_hurdle_and_claim_tree(proto):
    hits = []
    for sec in ("hurdle", "decision_rules"):
        for path, s in _walk(proto[sec], (sec,)):
            if sec == "decision_rules" and path[1] in ("morans_primary", "morans_comparators", "trajectory", "nmf"):
                continue
            if HISTORICAL & set(path):
                continue
            if PVAL.search(s):
                hits.append((".".join(path), s[:100]))
    assert not hits, hits
    # non-vacuity: the scan sees the strings it is supposed to scan, and the regex fires on a known live phrase
    assert sum(1 for _ in _walk(proto["hurdle"])) > 100
    assert PVAL.search("block-permutation q<0.05") and PVAL.search("BH over 10 tests")


def test_claim_tree_inputs_have_no_hurdle_pvalue_interval_or_control_input(proto):
    ins = proto["decision_rules"]["hurdle_branches_for_H1_per_part"]["inputs"]
    assert set(ins) == {"est", "s", "c", "d", "f", "robust"}
    assert not ({"a", "b", "e", "has_ci", "below_gate", "controls_pass"} & set(ins))
    assert proto["hurdle"]["validation_controls"]["real_data_post_freeze"]["controls_in_claim_tree"] is False


def test_revision_g6_and_checklist_say_failed_descriptive(proto):
    assert proto["revision"].startswith("9 ") and proto["protocol_id"] == "A12_midkine_remediation_v9"
    g6 = [x for x in proto["gates_fail_closed"]["list"] if x.startswith("G6")]
    assert len(g6) == 1 and "FAILED" in g6[0] and "DESCRIPTIVE" in g6[0]
    a123 = [x for x in proto["pre_freeze_checklist"] if x.startswith("A12.3")]
    assert len(a123) == 1 and "FAILED" in a123[0] and "not a pass" in a123[0]
    assert proto["hurdle"]["validation_controls"]["synthetic_pre_freeze"]["g6_gate_status"].startswith("FAILED")


def test_stale_items_removed(proto):
    s = json.dumps(proto)
    assert "+1.06" not in s and "-0.38" not in s
    assert "incl. log(n_st_eff)" not in s
    assert "drop_log_n_st_eff_covariate" not in proto["hurdle"]["validation_controls"]["mutation_checks_required"]
    assert "reintroduce_log_n_st_eff_covariate" in proto["hurdle"]["validation_controls"]["mutation_checks_required"]
    assert ">=90 %" not in s and "block-permutation q<0.05 in" not in s


def test_bootstrap_deviation_user_override_and_k0_fallback_recorded(proto):
    dev = proto["hurdle"]["pre_freeze_deviations"]
    assert len(dev) >= 2 and "500 resamples" in dev[0]["what"] and "2000" in dev[0]["what"] and "20261014" in dev[0]["what"]
    ov = dev[-1]
    assert ov["date"] == "2026-10-15" and "USER OVERRIDE" in ov["what"] and "g6_fix_round_selection_rule.md" in ov["what"] and "not trustworthy" in ov["what"]
    assert ov["decided_by"].startswith("user")
    b = proto["hurdle"]["inference"]["bootstrap"]
    assert b["n"] == 2000 and b["seed"] == 20261014 and b["min_valid_fraction"] == 0.95
    assert "never widened" in b["interval_policy"] and "NOMINAL" in b["interval_policy"]
    g = proto["genotype_groups"]
    assert g["high_threshold_k0"] == 0.01
    assert "NOT independent of the synthetic generator" in g["k0_justification"] and "2/202" in g["k0_justification"]
    assert "EXCLUDING the tested outcome MDK" in g["label_variable"]
    fb = g["k0_fallback"]
    assert fb["trigger"].startswith("ONLY if either group") and "fewer than 30" in fb["trigger"] and "NOT a trigger" in fb["trigger"]
    assert fb["algorithm_frozen"]["numpy_quantile_method"] == "linear" and fb["algorithm_frozen"]["q"] == 0.70 and "ties" in fb["fallback"] and fb["outcome_blind"] is True
    assert "NO further fallback" in fb["after_selection"]
    assert "k_st" not in g and "historical_label_k_st" in g


def test_participation_has_no_interval_h2_ranks_gmr_only_and_h3_is_descriptive(proto):
    q = proto["hurdle"]["questions"]
    assert "POINT ESTIMATE ONLY" in q["H1"]["estimands"]["participation"]
    assert "AMONG ACTIVE SENDERS" in q["H1"]["estimands"]["magnitude"] and "CONDITIONAL" in q["H1"]["estimands"]["magnitude"]
    assert "GMR ONLY" in q["H2"]["effect"] and "participation" in q["H2"]["effect"].lower()
    assert "ABSOLUTE" in q["H2"]["population"] and "contextual" in q["H2"]["population"] and "competition ranking" in q["H2"]["tie_rule"]
    assert q["H2"]["significant_wording_allowed"] is False and q["H2"]["stands_out_wording_allowed"] is False
    assert "DESCRIPTIVE SIDE-BY-SIDE" in q["H3"]["what"] and "No contrast statistic" in q["H3"]["inference"]
    assert q["H1"]["headline_for"] == ["HCC22-088-P4-S2_1i_rep"]
    assert "proxy for ESR1-D538G" in q["H1"]["wording_guard"]
    assert not any("0.874-0.900" in s for _, s in _walk(proto["hurdle"]["questions"]["H1"]["required_statements_with_the_headline"]))
    assert "0.874, 0.896, 0.898 and 0.900" in " ".join(q["H1"]["required_statements_with_the_headline"])


def test_open_question_is_labelled_non_blocking(proto):
    oq = proto["open_questions_pending_user"]
    assert len(oq) == 1 and oq[0]["id"] == "genotype_citable_record" and "NON-BLOCKING FOR COMPUTATION" in oq[0]["status"]


def test_no_stale_per_cell_label_or_confirmatory_wording_outside_history(proto):
    s = json.dumps({k: v for k, v in proto.items() if k not in ("genotype_groups",)})
    assert "k_st" not in s.replace("k_st_", "")  # no live k_st label anywhere else in the protocol
    assert "confirmatory" not in json.dumps(proto).lower()
    assert "envs/COMMOT" in proto["commot"]["env"] and "REQUIRED" in proto["commot"]["env"]
    assert proto["commot"]["run_matrix"]["analysis_side_variants_no_commot_rerun"]["A03"].startswith("k0 0.005 on k_share")
    assert "same k_share labels" in proto["commot"]["run_matrix"]["runs"]["R11"]
