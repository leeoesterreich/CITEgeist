"""B22 / B221: robustness of the MDK-secretory bivariate Moran's I (Definition C, B9) and its matched controls (B17) to sequencing depth, shared denominators, neighbourhood size
and the RNA-per-cell covariate. Descriptive only; nothing here changes the frozen k=6 primary or is called a cause.

Stages: B22 (run v2b22) ran checks 9, 10, 11, 13, 14, 15. B221 (run v2b221, first rerun of B22, same pattern as X5 -> X51) runs the same six checks and adds check 17 (``d5_fixed_moran``:
forms none, D5, D5x, D5xc); v2b22 stays registered and is no longer read. Since the protocol deviation 2026-11-06 (``PROTOCOL_DEVIATION_2026-11-06_rna_total.md``) the stage builds every
``hurdle.covariate_matrix`` with ``rna="total"`` (RNA content = the raw allocated-count total): check 9's adjustment, check 11's D3, check 13's D3-D7, and D5 of checks 13 and 17. The pre-deviation
per-cell D5 is the check 13 arm ``D5_percell``. Run ``--stage B22`` only to reproduce the code path; v2b22 itself was made before the deviation (per-cell covariate).

Promoted from the exploratory checks 9, 10, 11, 13, 14, 15 and 17 (``exploratory_v42/check{9,10,11,13,14,15,17}_*.py``, git history). The analysis code lives in the modules
``depth_thin`` (check 9), ``depth_moran`` (10), ``depth_removed_moran`` (11), ``split_d3`` (13), ``share_moran`` (14), ``k_sensitivity`` (15) and ``d5_fixed_moran`` (17), unchanged; this file is the stage
entry that runs the checks IN THIS ORDER in ONE job and writes ONE artifact, ``moran_robustness_results.json``: per check the summary fields of the exploratory
``check*_summary.json`` (gates, cohort lines, per-specimen records), the exploratory CSV tables as records (``tables``), the cohort table (``cohort_table``, one row per
analysis arm / form / k, derived by ``cohort_row`` from the per-specimen results), the sha256 of the csv / markdown files written beside it, and a ``config`` block.

Gates (fail closed; ``gate_failure.json``, no artifact, no stamp, exit 2)
  * check 9 Part 1: the primary and the adjusted rho recomputed here equal the stored check 7 / check 8 values (1e-6 relative), read from the stamped B21 artifact
    (the exploratory gate read the check 7 / 8 summary files).
  * checks 10, 11, 14, 15: the unthinned / unprojected recompute equals the stored B9 Definition C I (1e-9) and p (exactly); checks 11, 14, 15 also reproduce the stored B17
    control blocks (1e-9). The stored B9 / B17 / B12 values come from the stamped B9, B17 and B12 stage directories (inputs).
  * check 13 gates D1 and D3 against the check 11 result computed earlier IN THE SAME RUN (the exploratory gate read the check 11 summary file); check 15 gates its k=6 blocks
    against the check 14 and check 11 results computed earlier in the same run (it read their summary files). Same tolerance (1e-9), p exactly.
  * check 14: N_P, N_S and N_A recomputed from E and obs equal the stored layers (1e-4); check 15: the stored cohort mean I equals the fixed S1 size.
  * check 17: form none reproduces the stored B9 I (1e-9), p (exactly) and the B17 control blocks; form D5 reproduces check 13 of THIS run (status, I, p, control blocks, 1e-9); the
    recomputed rna_per_cell_equiv equals the stored obs column (relative 1e-4).
  * B221 only (``prior_arm_problems``): the stored v2b22 ``cohort_table`` arms that do NOT use ``hurdle.covariate_matrix`` (``IDENTICAL_PRIOR_ARMS``: check 9 and 10 equal-depth arms,
    check 11 none / D1 / D2, check 13 none / D1 / D4 / D6, every check 14 form, every check 15 form and k) are reproduced by this run's table within rtol = atol = 1e-12. The arms that use it
    (``MOVING_PRIOR_ARMS``: check 11 D3, check 13 D3 / D5 / D7) may move, but must exist. The diagnostic arm ``check13|full_counts_D5_percell`` (the pre-deviation D5) must equal the stored
    ``check13|full_counts_D5`` within 1e-9. The rows allowed beyond the stored ones are the check 17 rows and ``D5_percell``. A stored arm in neither list, a missing arm or any other
    difference fails closed (gate_failure.json, no artifact, no stamp).
Seeds, draw counts, surrogate counts, targets, forms and scenarios are the exploratory ones, written into ``config``. Verdict labels printed by the exploratory scripts
(``reading``, ``label``) are written as values, not gates.

Inputs come from the CLI (set by ``.commands/submit-v2-moran.sh`` through ``chain_preflight.py``): the B5 run (``--compartments-run``), the B9 run (``--morans-run``), the B17 run
(``--controls-run``), the B12 run (``--keratin-run``), the B21 run (``--e2-run``) and ``--out-root``. Code dependencies (not inputs): the committed ``regulon_inputs/`` tables
and ``_e2_scoring.py`` (through keratin_e2 / er_regulon), the committed ``protocol.json``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import chain_guard as PF  # noqa: E402
import depth_moran as DM  # noqa: E402
import depth_removed_moran as DR  # noqa: E402
import d5_fixed_moran as D5  # noqa: E402
import depth_thin as DT  # noqa: E402
import hurdle as H  # noqa: E402
import k_sensitivity as KS  # noqa: E402
import promoted_common as CM  # noqa: E402
import protocol_io as PIO  # noqa: E402
import share_moran as SM  # noqa: E402
import split_d3 as SD  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/moran_robustness.py"
STAGE = "B22"           # the first run (v2b22); B221 (v2b221) is the rerun with check 17 (same code, same artifact name)
STAGES = ("B22", "B221")
PRIOR_STAGE = {"B221": "B22"}   # the stage whose stored cohort table a rerun must reproduce
ARTIFACT = "moran_robustness_results.json"
RNA_COVARIATE = "total"   # protocol deviation 2026-11-06; every module below builds hurdle.covariate_matrix with the same mode (checked here, fail closed at import)
if not (DT.RNA_COVARIATE == DR.RNA_COVARIATE == SD.RNA_COVARIATE == RNA_COVARIATE):
    raise ImportError("depth_thin / depth_removed_moran / split_d3 / moran_robustness disagree on RNA_COVARIATE")
HEADLINE = DT.HEADLINE
CHECKS = ("check9", "check10", "check11", "check13", "check14", "check15", "check17")  # the order of the run
PRIOR_RTOL = PRIOR_ATOL = 1e-12
DIAG_RTOL = DIAG_ATOL = 1e-9       # the diagnostic per-cell D5 arm against the stored D5 (same tolerance as check 13's own gates)
# B221 prior-arm gate, explicit. Arms that do not use hurdle.covariate_matrix must reproduce the stored v2b22 row (1e-12); arms that do may move.
IDENTICAL_PRIOR_ARMS = tuple(
    [("check9", f"equal_depth_200_{tn}") for tn in DT.TESTS]
    + [("check10", f"equal_depth_200_{tn}") for tn in ("m1", "d_mdk", "d_sec", "k3")]
    + [("check11", f"full_counts_{fm}") for fm in DR.FORMS if fm != "D3"]
    + [("check13", f"full_counts_{fm}") for fm in ("none", "D1", "D4", "D6")]       # D4 = (depth, tumour content), D6 = tumour content: column subsets WITHOUT the RNA column (columns are standardised one by one), same seed labels
    + [("check14", f"share_{fm}") for fm in SM.FORMS14]
    + [("check15", f"{fm}_k{k}") for fm in KS.FORMS15 for k in KS.KS])
MOVING_PRIOR_ARMS = tuple([("check11", "full_counts_D3")] + [("check13", f"full_counts_{fm}") for fm in ("D3", "D5", "D7")])
DIAGNOSTIC_ARM = ("check13", "full_counts_D5_percell")      # the pre-deviation D5 (depth + RNA per cell) ...
DIAGNOSTIC_PRIOR = ("check13", "full_counts_D5")             # ... must equal the stored v2b22 D5
RULE_FILES = {c: f"Benchmarking/midkine_remediation/exploratory_v42/{c}_DECISION_RULE.md" for c in CHECKS}
READING = {
    "common": "Descriptive only. Per-specimen numbers are not corrected for multiple looks except the named BH families. One ESR1-D538G patient (P4-S2): no claim about the mutation or about cause. "
              "k=6 stays the primary and is frozen in the protocol.",
    "check9": "Take the specimens where the primary and the adjusted rho have OPPOSITE signs (discordant). If there are < 4 for a comparison: 'not decidable by sign'. Otherwise, if the equal-depth mean "
              "rho has the primary's sign in >= 2/3 of the discordant specimens: the adjusted form over-corrects; if it has the adjusted's sign in >= 2/3: the adjusted form is supported; else mixed. "
              "Magnitudes are NOT compared across forms (equal-depth rho is attenuated by sparsity): only signs, the cohort count of BH q < 0.05 and the control percentile. 'Under-corrects' is flagged "
              "only if the equal-depth cohort mean is closer to zero than the adjusted cohort mean and the equal-depth BH q < 0.05 count is below the adjusted one.",
    "check10": "R1: for each B9 q < 0.05 specimen, 'remains q<0.05 at equal depth' if its equal-depth BH q < 0.05, 'lost' if estimable and q >= 0.05, else 'not estimable'. R2: mean I remains above controls iff the "
               "cohort mean equal-depth real I exceeds the cohort mean of the control-a means AND of the control-b means. R3: number of estimable specimens with real above p95 in control a and in control b. "
               "R4: a loss is 'informative' if the S1 planted association is recovered (combined p < 0.05) in that specimen, else 'not informative (thinning could not see an association of that size)'; S2, S3 "
               "reported beside. R5: counts of BH q < 0.05 and means of I(log_spot_umi, x) unthinned and thinned. R6: B12 keratin-share vs secretory Spearman at equal depth beside the stored forms. "
               "A zero result at equal depth is read together with R4 and detection.",
    "check11": "Per form D1, D2, D3: R1 count of the six B9 q<0.05 specimens remaining q<0.05; R2 mean I stays above controls iff the cohort mean real I exceeds both control-mean cohort means; R3 number of "
               "specimens with real above p95 of control a and of control b; R4 a loss is informative only if the S1 plant is recovered under that form (combined p < 0.05), else not informative; S0 false "
               "positives and S1, S2, S3, A, S1A recovery counts; R5 counts of BH q<0.05 and mean I for I(depth, x_form) and I(depth, y_form) with Pearson r. Nothing is interpreted past these outcomes.",
    "check13": "removal(F) = 1 - M(F)/M(D1), M = cohort mean real I over the specimens estimable in F and D1. BOTH: removal(D4) >= 2/3 and removal(D5) >= 2/3. TUMOUR_CONTENT: removal(D4) >= 2/3 and "
               "removal(D5) < 1/3. RNA_PER_CELL: removal(D5) >= 2/3 and removal(D4) < 1/3. NEITHER_ALONE: both < 2/3. UNCLASSIFIED: one >= 2/3 and the other between 1/3 and 2/3. D6 and D7 get the numbers "
               "but no label. Protocol deviation 2026-11-06: D3-D7 use RNA content (log rna_total) where the first run used RNA per cell, and the label name RNA_PER_CELL now means 'RNA content'; the "
               "pre-deviation per-cell D5 is the diagnostic form D5_percell (no label).",
    "check14": "Cohort = the B9-estimable specimens; 'mean control' = cohort mean of the per-specimen control mean, a and b separately. DENOMINATOR_BASELINE: under BOTH S_A and S_B, BOTH mean control means are "
               "<= 0.034/3. SIGNAL_SURVIVES: under S_B, (cohort mean real I) minus (mean control a) >= 0.007 AND minus (mean control b) >= 0.007, AND at least 3 of the 6 B9 q<0.05 specimens remain q<0.05 "
               "under S_B. Both labels can be true; if neither holds: NEITHER.",
    "check17": "SELF_INCLUSION: under D5x, (cohort mean real I) minus (mean control a) >= 0.007 AND minus (mean control b) >= 0.007, AND at least 3 of the 6 B9 q<0.05 specimens remain q<0.05 under D5x "
               "(BH within D5x). COUPLING: the same three conditions under D5xc. NEITHER: both false. Both can be true; the label is the '+'-join of the true flags, or NEITHER. D5x = D5 with the RNA-per-cell "
               "covariate recomputed per pair without the pair's genes (same denominator); D5xc = the covariate is the log raw allocated-count total without the pair's genes (no shared denominator). "
               "Protocol deviation 2026-11-06: D5 itself is now depth + RNA content (the raw total), so D5xc differs from D5 only by the pair's own genes; the difference is reported (d5xc_minus_d5).",
    "check15": "Per form (none, S_B) and k: gt_a, gt_b = number of specimens whose real I exceeds the p95 of control a, b; ratio_a, ratio_b = (cohort mean real I)/(cohort mean control a, b mean), defined only "
               "when that control mean is > 0. WIDER_HELPS if some k in {15, 18} satisfies, under BOTH none and S_B, relative to k=6 of the same form: gt_a(k) >= gt_a(6)+2, gt_b(k) >= gt_b(6)+2, "
               "ratio_a(k) > ratio_a(6), ratio_b(k) > ratio_b(6). WIDER_NO_GAIN otherwise. D1 is reported but does not enter the label.",
}


def config() -> dict:
    """Every seed, draw count, surrogate count, threshold, label and reading rule of the stage (written into the artifact)."""
    return {
        "seeds": {"SEED9": DT.SEED9, "SEED10": DM.SEED10, "SEED11": DR.SEED11, "SEED14": SM.SEED14, "SEED15": KS.SEED15},
        "seed_rule": "SeedSequence([seed, crc32(label)]) per test / draw label (labels contain the specimen, test, draw index and target); control draws: control_pairs.make_draws",
        "k": DT.K.K, "tolerances": {"gate_equal_default": {"rtol": 1e-9, "atol": 1e-12}, "check11_TOL": DR.TOL, "check9_part1_gate": {"rtol": 1e-6, "atol": 1e-9}, "check14_reconstruction": SM.RECON_TOL,
                                   "check15_s1_size": {"rtol": 0, "atol": 1e-15}, "check17_check13_gate": {"rtol": 1e-9, "atol": 1e-9}, "prior_run_arms": {"rtol": PRIOR_RTOL, "atol": PRIOR_ATOL}},
        "rna_covariate": RNA_COVARIATE, "covariates": list(DT.COVS), "protocol_deviation": H.deviation_record(),
        "check9": {"target": DT.TARGET, "n_draws": DT.N_DRAWS, "min_valid_draws": DT.MIN_VALID, "n_surrogates": DT.B, "n_controls": DT.N_CTRL, "tests": {k: list(v) for k, v in DT.TESTS.items()},
                   "comparisons": {k: list(v) for k, v in DT.COMPARISONS.items()}, "covariates": list(DT.COVS), "min_eligible": DT.MIN_EL, "min_group": DT.MIN_GROUP, "min_mdk_nonzero": DT.MIN_NZ},
        "check10": {"target": DM.TARGET, "n_draws": DM.N_DRAWS, "min_valid_draws": DM.MIN_VALID, "b1": DM.B1, "b3": DM.B3, "n_ctrl": DM.N_CTRL, "min_keratin_nonzero": DM.MIN_KNZ,
                    "scenarios": list(DM.SCENARIOS), "b_max": DM.B_MAX, "b_steps": DM.B_STEPS},
        "check11": {"n_real": DR.N_REAL, "n_ctrl": DR.N_CTRL, "b_power": DR.R.B, "forms": list(DR.FORMS), "scenarios": list(DR.SCEN), "b_max": DR.B_MAX, "b_steps": DR.B_STEPS},
        "check13": {"forms": list(SD.FORMS13), "design_columns_of_covariate_matrix": {k: list(v) for k, v in SD.DESIGN.items()},
                    "percell_forms": dict(SD.PERCELL_FORMS), "rna_covariate_of_the_other_forms": SD.RNA_COVARIATE, "s1_size": SD.S1_SIZE},
        "check14": {"forms": list(SM.FORMS14), "scenarios": list(SM.SCEN14), "s1_size": SM.S1_SIZE, "ctrl_base": SM.CTRL_BASE, "excess_min": SM.EXCESS_MIN},
        "check15": {"ks": list(KS.KS), "forms": list(KS.FORMS15), "label_forms": list(KS.LABEL_FORMS), "scenarios": list(KS.SCEN15), "s1_size": KS.S1_SIZE},
        "check17": {"forms": list(D5.FORMS17), "scenarios": list(DR.SCEN), "excess_min": D5.EXCESS_MIN, "covariate_gate_rel_tol": D5.REL_TOL_COV},
        "reading_rules": READING, "rule_files": RULE_FILES,
        "labels": {"check13": ["BOTH", "TUMOUR_CONTENT", "RNA_PER_CELL", "NEITHER_ALONE", "UNCLASSIFIED"], "check14": ["DENOMINATOR_BASELINE", "SIGNAL_SURVIVES", "NEITHER"],
                   "check15": ["WIDER_HELPS", "WIDER_NO_GAIN"], "check17": ["SELF_INCLUSION", "COUPLING", "SELF_INCLUSION+COUPLING", "NEITHER"], "check9": ["adjusted over-corrects", "adjusted supported", "mixed", "not decidable by sign (< 4 discordant specimens)"],
                   "verdicts": "written as values, not gates"},
        "not_reproduced": "the 500-count sensitivity runs of checks 9 and 10 (the exploratory jobs with TARGET=500) are not part of this stage; the rule files fix the target at 200",
        "descriptive_only": True}


# --------------------------------------------------------------------------- cohort table (one small function, tested on a toy example)
def cohort_row(check: str, arm: str, statistic: str, records: dict) -> dict:
    """One row of the cohort table from per-specimen records ``{specimen: {"I", "a_mean", "b_mean", "gt_a", "gt_b", "q", "power"}}`` (only the specimens estimable in the arm).
    ``I`` = the arm's statistic (Moran's I, or Spearman rho for checks 9 and 10 Part 3); ``a_mean`` / ``b_mean`` = the control a / b mean of that specimen; ``gt_a`` / ``gt_b`` = real above the
    control p95; ``q`` = BH q within the arm; ``power`` = ``{scenario: combined p of the planted-association test, None if not estimable}``. A field no record has (None) gives None.
    Control means, ratios and above-p95 counts use the specimens that have controls (``n_with_controls``; ``mean_real_with_controls`` is the mean real I over those, the
    exploratory check 10 R2 line); ``mean_real`` and the q / power counts use every specimen estimable in the arm. Ratio = (mean real over specimens with controls) /
    (cohort mean control), None unless that control mean is > 0 (check 15 rule)."""
    vals = list(records.values())

    def mean(key):
        return float(np.mean([v[key] for v in vals])) if vals and all(v.get(key) is not None for v in vals) else None

    def count(key):
        vs = [v[key] for v in vals if v.get(key) is not None]
        return int(sum(bool(x) for x in vs)) if vs else None

    wc = [v for v in vals if v.get("a_mean") is not None]            # the specimens that have matched controls (check 10: controls are not estimable in 3 of 11)
    mi = mean("I")
    mic = float(np.mean([v["I"] for v in wc])) if wc else None
    ma = float(np.mean([v["a_mean"] for v in wc])) if wc else None
    mb = float(np.mean([v["b_mean"] for v in wc if v.get("b_mean") is not None])) if wc and all(v.get("b_mean") is not None for v in wc) else None

    def ratio(m):
        return mic / m if mic is not None and m is not None and m > 0 else None

    power = {}
    for sc in sorted({sc for v in vals for sc in (v.get("power") or {})}):
        ps = [(v.get("power") or {}).get(sc) for v in vals]
        power[sc] = {"n_estimable": int(sum(p is not None for p in ps)), "n_combined_p_lt_0.05": int(sum(p is not None and p < 0.05 for p in ps))}
    qs = [v.get("q") for v in vals]
    p = records.get(HEADLINE)
    return {"check": check, "arm": arm, "statistic": statistic, "n_specimens": len(vals), "mean_real": mi, "n_with_controls": len(wc), "mean_real_with_controls": mic, "mean_ctrl_a": ma, "mean_ctrl_b": mb, "ratio_a": ratio(ma), "ratio_b": ratio(mb),
            "n_above_p95_a": count("gt_a"), "n_above_p95_b": count("gt_b"),
            "n_q_lt_0.05": int(sum(q is not None and q < 0.05 for q in qs)) if any(q is not None for q in qs) else None, "power": power,
            "p4s2": ({"specimen": HEADLINE, "real": p.get("I"), "ctrl_a_mean": p.get("a_mean"), "ctrl_b_mean": p.get("b_mean"), "above_p95_a": p.get("gt_a"), "above_p95_b": p.get("gt_b"),
                      "q": p.get("q"), "power": dict(p.get("power") or {})} if p else None)}


def _pc(v):
    """Combined p of a power / test aggregate, None if it is not estimable."""
    return v["p_combined"] if v and v.get("status") == "ok" else None


def _form_record(v: dict, power: dict) -> dict:
    c = v["controls"]
    return {"I": v["I"], "a_mean": c["a"]["mean"], "b_mean": c["b"]["mean"], "gt_a": c["a"]["gt95"], "gt_b": c["b"]["gt95"], "q": v.get("q_bh"), "power": power}


def records_forms11(per: dict, fm: str) -> dict:
    """Checks 11 and 13 (power under ``power.scenarios[sc].forms[fm]``)."""
    return {s: _form_record(r["forms"][fm], {sc: _pc(b.get("forms", {}).get(fm)) for sc, b in r["power"]["scenarios"].items()})
            for s, r in per.items() if r.get("status") == "ok" and r["forms"][fm].get("status") == "ok"}


def records_forms14(per: dict, fm: str) -> dict:
    return {s: _form_record(r["forms"][fm], {sc: _pc(v) for sc, v in r["power"][fm].items()}) for s, r in per.items() if r.get("status") == "ok" and r["forms"][fm].get("status") == "ok"}


def records_forms15(per: dict, fm: str, k: int) -> dict:
    return {s: _form_record(r["k"][str(k)]["forms"][fm], {sc: _pc(r["k"][str(k)]["power"][fm][sc]) for sc in KS.SCEN15})
            for s, r in per.items() if r.get("status") == "ok" and r["k"][str(k)]["forms"][fm].get("status") == "ok"}


def records_check10(per: dict, tn: str) -> dict:
    """Check 10 equal-depth test ``tn`` (m1 = the MDK-secretory I with its controls and power; d_mdk, d_sec, k3 have no controls / power)."""
    out = {}
    for s, r in per.items():
        v = r["tests"][tn] if r.get("status") == "ok" else None
        if not v or v["status"] != "ok":
            continue
        rec = {"I": v["mean"], "q": v.get("q_bh")}
        if tn == "m1":
            c = r["controls"]
            if c["status"] == "ok":
                rec.update({"a_mean": c["a"]["mean"], "b_mean": c["b"]["mean"], "gt_a": c["a"]["specimen_gt_p95"], "gt_b": c["b"]["specimen_gt_p95"]})
            rec["power"] = {sc: _pc(r["power"][sc]) for sc in DM.SCENARIOS}
        out[s] = rec
    return out


def records_check9(part2: dict, tn: str) -> dict:
    out = {}
    for s, r in part2.items():
        v = r["tests"][tn] if r.get("status") == "ok" else None
        if not v or v["status"] != "ok":
            continue
        out[s] = {"I": v["rho_mean"], "q": v.get("q_bh"), "gt_a": (v["n_draws_gt_p95"] > v["n_valid"] / 2) if "n_draws_gt_p95" in v else None}
    return out


def cohort_table(checks: dict) -> list:
    """The cohort table: one row per analysis arm / form / k / thinning, from the stage's own per-specimen results."""
    rows = []
    for tn in DT.TESTS:
        rows.append(cohort_row("check9", f"equal_depth_200_{tn}", "spearman_rho", records_check9(checks["check9"]["part2"], tn)))
    p10 = checks["check10"]["per_specimen"]
    for tn, stat in (("m1", "morans_I"), ("d_mdk", "morans_I_depth_vs_MDK"), ("d_sec", "morans_I_depth_vs_secretory"), ("k3", "spearman_rho_keratin_share_vs_secretory")):
        rows.append(cohort_row("check10", f"equal_depth_200_{tn}", stat, records_check10(p10, tn)))
    for fm in DR.FORMS:
        rows.append(cohort_row("check11", f"full_counts_{fm}", "morans_I", records_forms11(checks["check11"]["per_specimen"], fm)))
    for fm in SD.FORMS13:
        rows.append(cohort_row("check13", f"full_counts_{fm}", "morans_I", records_forms11(checks["check13"]["per_specimen"], fm)))
    for fm in SM.FORMS14:
        rows.append(cohort_row("check14", f"share_{fm}", "morans_I", records_forms14(checks["check14"]["per_specimen"], fm)))
    for fm in KS.FORMS15:
        for k in KS.KS:
            rows.append(cohort_row("check15", f"{fm}_k{k}", "morans_I", records_forms15(checks["check15"]["per_specimen"], fm, k)))
    for fm in D5.FORMS17:
        rows.append(cohort_row("check17", f"full_counts_{fm}", "morans_I", records_forms11(checks["check17"]["per_specimen"], fm)))
    return rows


# --------------------------------------------------------------------------- B221 gate: every arm of the stored v2b22 table is reproduced
def _same(a, b, path: str, out: list, rtol: float = PRIOR_RTOL, atol: float = PRIOR_ATOL) -> None:
    """Append the differences between ``a`` (this run) and ``b`` (stored): numbers rtol / atol (default 1e-12; NaN equals NaN), everything else exactly."""
    if isinstance(b, dict):
        if not isinstance(a, dict) or a.keys() != b.keys():
            out.append(f"{path}: keys differ ({sorted(a) if isinstance(a, dict) else type(a).__name__} vs {sorted(b)})")
            return
        for k in b:
            _same(a[k], b[k], f"{path}/{k}", out, rtol, atol)
    elif isinstance(b, list):
        if not isinstance(a, list) or len(a) != len(b):
            out.append(f"{path}: list length differs")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            _same(x, y, f"{path}[{i}]", out, rtol, atol)
    elif isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str) or isinstance(b, str):
        if a != b or type(a) is not type(b):
            out.append(f"{path}: {a!r} vs {b!r}")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not ((math.isnan(a) and math.isnan(b)) or math.isclose(a, b, rel_tol=rtol, abs_tol=atol)):
            out.append(f"{path}: {a!r} vs {b!r}")
    elif a != b:
        out.append(f"{path}: {a!r} vs {b!r}")


def prior_arm_problems(rows: list, prior_rows: list) -> list:
    """Problems found comparing this run's cohort table ``rows`` with the stored one ``prior_rows`` (the v2b22 artifact's ``cohort_table``); [] = the gate passes.
    * every stored arm must be in ``IDENTICAL_PRIOR_ARMS`` (reproduced within 1e-12, fields and all) or in ``MOVING_PRIOR_ARMS`` (must exist here, may move); a stored arm in neither list fails;
    * every arm of ``IDENTICAL_PRIOR_ARMS`` must be stored (a list entry the stored table lacks fails: the lists and the table must agree);
    * the diagnostic arm ``DIAGNOSTIC_ARM`` (the pre-deviation per-cell D5) must equal the stored ``DIAGNOSTIC_PRIOR`` (D5) within 1e-9;
    * the only arms allowed beyond the stored ones are the check 17 rows and ``DIAGNOSTIC_ARM``.
    Both tables must be non-empty and (check, arm) keys unique."""
    key = lambda r: (r["check"], r["arm"])  # noqa: E731
    new, old = {key(r): r for r in rows}, {key(r): r for r in prior_rows}
    out = []
    if not prior_rows or not rows:
        out.append("empty cohort table")
    if len(new) != len(rows) or len(old) != len(prior_rows):
        out.append("cohort table: (check, arm) keys are not unique")
    ident, moving = set(IDENTICAL_PRIOR_ARMS), set(MOVING_PRIOR_ARMS)
    out += [f"arm {k}: stored in the prior run but in neither IDENTICAL_PRIOR_ARMS nor MOVING_PRIOR_ARMS (unclassified)" for k in sorted(old.keys() - ident - moving)]
    out += [f"arm {k}: listed as identical but not stored in the prior run" for k in sorted(ident - old.keys())]
    out += [f"arm {k}: stored in the prior run, missing here" for k in sorted(old.keys() - new.keys())]
    out += [f"arm {k}: not in the prior run and neither a check 17 arm nor the diagnostic arm" for k in sorted(new.keys() - old.keys()) if k[0] != "check17" and k != DIAGNOSTIC_ARM]
    for k in sorted((old.keys() & new.keys()) & ident):
        _same(new[k], old[k], f"arm {k[0]}|{k[1]}", out)
    if DIAGNOSTIC_ARM not in new:
        out.append(f"diagnostic arm {DIAGNOSTIC_ARM} missing")
    elif DIAGNOSTIC_PRIOR in old:
        _same({**new[DIAGNOSTIC_ARM], "arm": DIAGNOSTIC_PRIOR[1]}, old[DIAGNOSTIC_PRIOR], f"diagnostic {DIAGNOSTIC_ARM[0]}|{DIAGNOSTIC_ARM[1]} vs stored {DIAGNOSTIC_PRIOR[1]}", out, DIAG_RTOL, DIAG_ATOL)
    return out


def prior_run_dir(stage: str, out_root, protocol_path=None) -> Path:
    """The stored run a rerun is gated against: the registered run id of ``PRIOR_STAGE[stage]``, beside ``out_root`` (never a literal id)."""
    import chain_preflight as CPF

    return Path(out_root).parent / CPF.run_id_for(PRIOR_STAGE[stage], CPF.protocol_sha12(protocol_path))


# --------------------------------------------------------------------------- tables, hashes, json
def _coerce(v):
    if v is None or v == "":
        return None
    if v in ("True", "False"):
        return v == "True"
    if re.fullmatch(r"-?\d+", v):
        return int(v)
    try:
        return float(v)
    except ValueError:
        return v


def read_table(path) -> list:
    """A csv written beside the artifact as records (empty -> None, numbers -> int / float)."""
    with open(path, newline="") as fh:
        return [{k: _coerce(v) for k, v in r.items()} for r in csv.DictReader(fh)]


def sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError(f"not serialisable: {type(o)}")


def specimens_for(protocol: dict, b9: dict, subset=None) -> tuple:
    """(specimens of checks 9 and 10 = the protocol list, P4-S2 first; specimens of checks 11, 13, 14, 15 = the B9-estimable ones in B9 order). ``subset`` restricts both (tests only)."""
    s_all = sorted(protocol["specimens"]["ordered"], key=lambda s: (s != HEADLINE, s))
    s_b9 = [s for s, r in b9["per_specimen"].items() if r.get("status") == "ok"]
    if subset:
        s_all, s_b9 = [s for s in s_all if s in subset], [s for s in s_b9 if s in subset]
    return s_all, s_b9


def compute(comp_run, b9: dict, b17: dict, b12: dict, e2: dict, protocol: dict, out, meta_common: dict, specimens=None) -> dict:
    """Run the seven checks in order in one process. ``b9`` / ``b17`` / ``b12`` / ``e2`` = the stamped B9 / B17 / B12 / B21 results. Returns ``{"checks", "files"}``; raises ``GateFailureAt``
    (a ``GateFailure`` that names the checks that had finished)."""
    s_all, s_b9 = specimens_for(protocol, b9, specimens)
    stored7, stored8 = e2["check7"]["per_specimen"], e2["check8"]["per_specimen"]
    checks, files = {}, {}
    steps = (
        ("check9", lambda: DT.compute(comp_run, stored7, stored8, s_all, out, meta_common)),
        ("check10", lambda: DM.compute(comp_run, b9, b17, b12, s_all, out, meta_common)),
        ("check11", lambda: DR.compute(comp_run, b9, b17, s_b9, out, meta_common)),
        ("check13", lambda: SD.compute(comp_run, b9, b17, checks["check11"], s_b9, out, meta_common)),  # gates D1, D3 against the check 11 result of THIS run
        ("check14", lambda: SM.compute(comp_run, b9, b17, s_b9, out, meta_common)),
        ("check15", lambda: KS.compute(comp_run, b9, b17, checks["check14"], checks["check11"], s_b9, out, meta_common)),  # gates k=6 against the check 14 / check 11 results of THIS run
        ("check17", lambda: D5.compute(comp_run, b9, b17, checks["check13"], s_b9, out, meta_common)),  # gates D5 against the check 13 result of THIS run
    )
    for name, step in steps:
        t = time.time()
        try:
            checks[name], files[name] = step()
        except CM.GateFailure as e:
            raise GateFailureAt(str(e), list(checks)) from e
        print(f"{name} done {time.time() - t:.0f}s", flush=True)
    return {"checks": checks, "files": files}


class GateFailureAt(CM.GateFailure):
    """A gate failure that remembers which checks had finished."""

    def __init__(self, msg: str, done: list):
        super().__init__(msg)
        self.done = done


def run(comp_run, morans_run, controls_run, keratin_run, e2_run, out_root, run_id: str = "", specimens=None, protocol_path=None, stage: str = STAGE, prior_run=None) -> int:
    """Compute, write ``moran_robustness_results.json`` + csv / markdown tables + ``_COMPLETE.json``. 0 ok, 2 gate failure (``gate_failure.json``, nothing else of the result).
    ``prior_run`` = the stored v2b22 run dir: every arm of its cohort table must be reproduced (B221; required for it, ``main`` resolves it)."""
    protocol = PIO.load_protocol(protocol_path or HERE / "protocol.json")
    out = Path(out_root)
    out.mkdir(parents=True, exist_ok=True)
    rd = lambda d, name: json.loads((Path(d) / name).read_text())  # noqa: E731
    b9, b17, b12, e2 = rd(morans_run, "morans_results.json"), rd(controls_run, "control_pairs_results.json"), rd(keratin_run, "keratin_secretory_summary.json"), rd(e2_run, "e2_covariation_results.json")
    code = HERE.parents[1]
    head = subprocess.run(["git", "-C", str(code), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    prior = json.loads((Path(prior_run) / ARTIFACT).read_text()) if prior_run else None
    meta = {"stage": stage, "run_id": run_id, "produced_by": PRODUCED_BY, "descriptive_only": True, "config": config(),
            "inputs": {"compartments_run": str(comp_run), "morans_results": str(Path(morans_run) / "morans_results.json"), "control_pairs_results": str(Path(controls_run) / "control_pairs_results.json"),
                       "keratin_secretory_summary": str(Path(keratin_run) / "keratin_secretory_summary.json"), "e2_covariation_results": str(Path(e2_run) / "e2_covariation_results.json")},
            "blas_threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, "git_commit": head, "specimens_subset": sorted(specimens) if specimens else None}
    if prior is not None:
        meta["inputs"]["prior_run_result"] = {"path": str(Path(prior_run) / ARTIFACT), "sha256": sha256(Path(prior_run) / ARTIFACT)}
    try:
        res = compute(comp_run, b9, b17, b12, e2, protocol, out, {"git_commit": head, "out": str(out)}, specimens)
        if prior is not None:
            bad = prior_arm_problems(json.loads(json.dumps(cohort_table(res["checks"]), default=_json_default)), prior["cohort_table"])  # JSON round trip: the stored side is JSON
            if bad:
                raise GateFailureAt(f"GATE FAILED {len(bad)} difference(s) from the stored {prior.get('run_id')} cohort table (rtol = atol = {PRIOR_RTOL}): " + "; ".join(bad[:10]), list(res["checks"]))
            meta["prior_arms_gate"] = {"prior_run_id": prior.get("run_id"), "n_arms_stored": len(prior["cohort_table"]), "n_identical_arms_compared": len(IDENTICAL_PRIOR_ARMS),
                                       "identical_arms": [list(k) for k in IDENTICAL_PRIOR_ARMS], "moving_arms": [list(k) for k in MOVING_PRIOR_ARMS],
                                       "diagnostic": {"arm": list(DIAGNOSTIC_ARM), "equals_stored": list(DIAGNOSTIC_PRIOR), "rtol": DIAG_RTOL, "atol": DIAG_ATOL},
                                       "rtol": PRIOR_RTOL, "atol": PRIOR_ATOL, "passed": True}
    except GateFailureAt as e:
        PIO.atomic_write_json(out / "gate_failure.json", {**meta, "gate_failure": str(e), "checks_done": e.done})
        print(str(e), file=sys.stderr, flush=True)
        return 2
    checks, files = res["checks"], res["files"]
    tables = {Path(n).stem: read_table(out / n) for fl in files.values() for n in fl if n.endswith(".csv")}
    meta["files_sha256"] = {n: sha256(out / n) for fl in files.values() for n in fl}
    body = {"checks": checks, "tables": tables, "cohort_table": cohort_table(checks)}
    PIO.atomic_write_json(out / ARTIFACT, json.loads(json.dumps({**meta, **body}, default=_json_default)))
    PIO.atomic_write_json(out / "_COMPLETE.json", {"path": str(out), "jobs": [os.environ.get("SLURM_JOB_ID") or "no-slurm"], "sha": head,
                                                  "completed": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "produced_by": PRODUCED_BY})
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compartments-run", type=Path, required=True, help="the B5 run dir (id from v42_run_ids.yaml)")
    ap.add_argument("--morans-run", type=Path, required=True, help="the B9 run dir; its Definition C values are the gate")
    ap.add_argument("--controls-run", type=Path, required=True, help="the B17 run dir; its control blocks are the gate")
    ap.add_argument("--keratin-run", type=Path, required=True, help="the B12 run dir; the stored keratin-secretory Spearman values")
    ap.add_argument("--e2-run", type=Path, required=True, help="the B21 run dir; the stored check 7 / 8 values gate check 9 Part 1")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<B22 id>; never overwritten")
    ap.add_argument("--run-id", required=True, help="the id of --stage in v42_run_ids.yaml")
    ap.add_argument("--stage", choices=STAGES, default=STAGE, help="B22 = the first run (no check 17); B221 = the rerun with check 17, gated against the registered B22 run")
    ap.add_argument("--prior-run", type=Path, default=None, help="B221 only: the stored B22 run dir; default = the registered B22 id beside --out-root")
    ap.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry(a.stage, a.protocol, run_id=a.run_id, out_root=a.out_root, comp_run=a.compartments_run, morans_run=a.morans_run, controls_run=a.controls_run, keratin_run=a.keratin_run, e2_run=a.e2_run)
    if a.out_root.exists():
        raise SystemExit(f"refusing: {a.out_root} exists (a run dir is never overwritten)")
    prior = (a.prior_run or prior_run_dir(a.stage, a.out_root, a.protocol)) if a.stage in PRIOR_STAGE else None
    if prior is not None and not (Path(prior) / ARTIFACT).is_file():
        raise SystemExit(f"refusing: the stored run to gate against has no {ARTIFACT}: {prior}")
    stamp = PF.stage_stamp_begin(a.stage, a.protocol, a.out_root / ARTIFACT, Path(__file__).resolve(), a.run_id)  # job-start stamp gate, before any compute
    rc = run(a.compartments_run, a.morans_run, a.controls_run, a.keratin_run, a.e2_run, a.out_root, a.run_id, protocol_path=a.protocol, stage=a.stage, prior_run=prior)
    if rc == 0 and stamp:
        stamp.finish(comp_run=a.compartments_run, morans_run=a.morans_run, controls_run=a.controls_run, keratin_run=a.keratin_run, e2_run=a.e2_run)
    print("done" if rc == 0 else "gate failed", a.run_id, a.out_root, flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
