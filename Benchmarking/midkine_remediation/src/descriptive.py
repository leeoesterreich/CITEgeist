"""DESCRIPTIVE H1 estimator and runner (protocol rev 7; spec 9.4.7-9.4.8). No p-value or q-value exists here.

Per specimen, on the eligible cancer compartments of one run:
  * magnitude  = geometric-mean ratio GMR (high / low) of the MDK-NCL sender score AMONG ACTIVE SENDERS (score > 0),
                 exp(group coefficient) of  log(score | active) ~ group + covariates + block fixed effects, with a NOMINAL
                 95% block-bootstrap percentile interval (2000 resamples of 6-pitch tiles, seed 20261014). The interval is
                 never widened; its measured synthetic coverage is carried in the output. A resample that fails (rank
                 deficiency, an empty group, a non-finite fit) is invalid; the interval needs >= MIN_VALID_FRACTION valid
                 resamples, otherwise it is a typed not_estimable. Wording never depends on whether the interval excludes 1.
  * participation = covariate-adjusted odds ratio of being an active sender (logistic, group + covariates): POINT ESTIMATE
                 ONLY, no interval.
Typed failures (never NaN with status ok): group_size, active_senders, mixed_blocks, covariate_domain, nonfinite_input,
collinear, separation, nonconvergence, nonfinite_estimate, bootstrap_valid_fraction, point_estimate_not_estimable, commot_run_failed.

Production entry (rev 7): ``run_h1_specimen`` takes NO tunable setting (2000 resamples, seed 20261014, 6-pitch tiles, interval on R00
only are frozen), requires R00/R01/R02/R10/R11 for the headline specimen (a missing run raises ``MissingRunError``; a failed COMMOT run is a
typed result), and ``write_h1`` only accepts a result whose recorded settings equal the loaded FROZEN protocol. Tests that need a
smaller bootstrap call ``_run_h1_for_tests``; its provenance records the settings actually used, so it can never be stamped as production.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Dict, Optional

import numpy as np

import claim_tree as CT
import hurdle as H
import protocol_io as PIO
from compartments import select_label  # noqa: F401  (re-exported: the label is chosen outside, from k_share only)

PROTOCOL_ID = "A12_midkine_remediation_v9"
HEADLINE_SPECIMEN = CT.SCOPE_SPECIMEN
BOOT_N = 2000
BOOT_SEED = 20261014
MIN_VALID_FRACTION = 0.95      # tails hold 5% of the resamples; if more than 5% are dropped the nominal endpoints are not interpretable
LEVEL = 0.95
MEASURED_GMR_COVERAGE = (0.874, 0.896, 0.898, 0.900)   # nominal 95% interval, 500 synthetic datasets x 4 null cases (synthetic/g6_fix/boot_*_s6.json)
PRODUCTION_PATH_COVERAGE = {"null_autocorrelated": 0.884, "null_noisy_nst_kshare": 0.900, "null_depth_confounded": 0.902, "null_autocorrelated_stress": 0.904}   # synthetic/g6_fix/prod_coverage/prod_*_s6.json: SAME nominal interval through gmr_interval at 2000 resamples (rev 7)
SEPARATION_BETA = 15.0
ESTIMAND_MAGNITUDE = "geometric-mean ratio (high/low) of the MDK-NCL sender score among active senders (compartments with a nonzero score)"
ESTIMAND_PARTICIPATION = "covariate-adjusted odds ratio (high/low) of being an active sender (nonzero MDK-NCL sender score)"
RUNS_FOR_CLAIM_TREE = ("R00", "R01", "R02", "R10", "R11")
REASON_CODES = ("group_size", "active_senders", "mixed_blocks", "covariate_domain", "nonfinite_input", "collinear", "separation",
                "nonconvergence", "nonfinite_estimate", "bootstrap_valid_fraction", "point_estimate_not_estimable", "commot_run_failed")
MIN_ELIGIBLE = 2 * H.MIN_PER_GROUP        # early gate: below this no label can give 30 per group; checked before pitch, blocks or quantiles
INTERVAL_RUNS = ("R00",)
H4_TOL = 1e-12
RNA_COVARIATE = "total"        # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md): the RNA covariate of every production hurdle run


class ScopeError(ValueError):
    """The claim tree is defined for the headline specimen only."""


def _nonest(code: str, reason: str, **extra) -> dict:
    assert code in REASON_CODES, code
    return {"status": "not_estimable", "reason_code": code, "reason": reason, **extra}


# --------------------------------------------------------------------------- point estimates
def _check_inputs(score) -> Optional[dict]:
    if not np.isfinite(score).all() or (score < 0).any():
        return _nonest("nonfinite_input", "sender score has non-finite or negative entries")
    return None


def _design(score, group, cov, block):
    """Active-sender design [group, covariates, ALL block dummies] with NO intercept (equivalent to intercept + dummies minus one,
    but stays full rank in a bootstrap resample that drops the reference block). The group coefficient is column 0."""
    act = score > 0
    g = np.asarray(group).astype(float)
    _, b = np.unique(np.asarray(block)[act], return_inverse=True)
    n = int(act.sum())
    cols = [g[act]]
    if cov is not None and np.asarray(cov).shape[1]:
        cols += list(np.asarray(cov, float)[act].T)
    Dm = np.zeros((n, int(b.max()) + 1))
    Dm[np.arange(n), b] = 1.0
    return np.log(score[act]), np.column_stack(cols + [Dm]), act, b


def fit_magnitude(score, group, cov, block, gate_group=None) -> dict:
    """GMR among active senders (OLS with block fixed effects); typed not_estimable instead of any NaN.
    ``gate_group`` (default: ``group``) is the binary label the estimability gates are evaluated on (A05: ``group`` is the continuous label)."""
    score = np.asarray(score, float)
    bad = _check_inputs(score)
    if bad:
        return bad
    group = np.asarray(group)
    gate = (group if gate_group is None else np.asarray(gate_group)).astype(int)
    block = np.asarray(block)
    est = H.estimability(gate, block, "magnitude", score > 0)
    if not est["ok"]:
        return _nonest(est["reason_code"], est["reason"])
    y, X, act, _ = _design(score, group, cov, block)
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank < X.shape[1]:
        return _nonest("collinear", f"design rank {rank} < {X.shape[1]} columns (group term collinear with covariates/blocks)")
    with np.errstate(over="ignore"):
        gmr = float(np.exp(beta[0]))                         # the ACTUAL exponentiated estimate must be finite (no ok with gmr = inf)
    if not np.isfinite(beta[0]) or not math.isfinite(gmr) or not gmr > 0:
        return _nonest("nonfinite_estimate", "group coefficient or its exponential (the GMR) is not finite")
    g = gate[act]
    return {"status": "ok", "reason": "", "reason_code": "", "value": gmr, "gmr": gmr, "coef": float(beta[0]), "estimand": ESTIMAND_MAGNITUDE,
            "n_active": int(act.sum()), "n_active_high": int(g.sum()), "n_active_low": int((1 - g).sum()), "n_mixed_blocks": est["n_mixed_blocks"]}


def fit_participation(score, group, cov, block, gate_group=None) -> dict:
    """Odds ratio of being an active sender (logistic, group + covariates); point estimate only. ``gate_group`` as in fit_magnitude."""
    score = np.asarray(score, float)
    bad = _check_inputs(score)
    if bad:
        return bad
    g = np.asarray(group, float)
    gate = (g if gate_group is None else np.asarray(gate_group)).astype(int)
    block = np.asarray(block)
    est = H.estimability(gate, block, "participation")
    if not est["ok"]:
        return _nonest(est["reason_code"], est["reason"])
    a = (score > 0).astype(float)
    X = np.column_stack([np.ones(len(a)), g] + ([np.asarray(cov, float)] if cov is not None and np.asarray(cov).shape[1] else []))
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return _nonest("collinear", "logistic design is rank deficient")
    try:
        beta, mu, conv = H._irls(X, a)
    except np.linalg.LinAlgError as e:
        return _nonest("collinear", f"singular information matrix: {e}")
    if not np.isfinite(beta).all():
        return _nonest("nonfinite_estimate", "logistic coefficients are not finite")
    if np.abs(beta).max() > SEPARATION_BETA or ((mu < 1e-8) | (mu > 1 - 1e-8)).any():
        return _nonest("separation", f"|coefficient| {np.abs(beta).max():.1f} > {SEPARATION_BETA} or fitted probabilities numerically 0/1 (complete or quasi-complete separation)")
    if not conv:
        return _nonest("nonconvergence", "IRLS did not converge")
    with np.errstate(over="ignore"):
        orr = float(np.exp(beta[1]))
    if not math.isfinite(orr):
        return _nonest("nonfinite_estimate", "odds ratio is not finite")
    return {"status": "ok", "reason": "", "reason_code": "", "value": orr, "odds_ratio": orr, "coef": float(beta[1]), "estimand": ESTIMAND_PARTICIPATION,
            "n_active": int(a.sum()), "n_mixed_blocks": est["n_mixed_blocks"]}


# --------------------------------------------------------------------------- interval
def _weighted_group_coef(y, X, w) -> Optional[float]:
    """Group coefficient (column 0) of the weighted LS fit; None if the resample is invalid (empty group, rank deficiency, non-finite)."""
    keep = w > 0
    Xk, yk, wk = X[keep], y[keep], w[keep]
    if wk[Xk[:, 0] > 0.5].sum() <= 0 or wk[Xk[:, 0] < 0.5].sum() <= 0:
        return None
    Xk = Xk[:, np.abs(Xk).sum(axis=0) > 0]               # drop block columns absent from this resample (the group column always survives)
    if Xk.shape[0] <= Xk.shape[1]:
        return None
    A = Xk.T @ (wk[:, None] * Xk)
    ev = np.linalg.eigvalsh(A)
    if not (ev[0] > 1e-10 * ev[-1]):
        return None
    try:
        beta = np.linalg.solve(A, Xk.T @ (wk * yk))
    except np.linalg.LinAlgError:
        return None
    return float(beta[0]) if np.isfinite(beta[0]) else None


def gmr_interval(score, group, cov, block, B: int = BOOT_N, seed: int = BOOT_SEED, min_valid_fraction: float = MIN_VALID_FRACTION, coverage: Optional[dict] = None) -> dict:
    """Nominal 95% block-bootstrap percentile interval for the GMR; typed not_estimable if the point estimate or too many resamples fail.
    ``coverage`` = {"values": [...], "source": {...}} (the stamped B24 record, ``coverage_record``); None = the per-cell series ``MEASURED_GMR_COVERAGE``."""
    meta = {"level": LEVEL, "nominal": True, "widened": False, "method": "block bootstrap, resample 6-pitch tiles with replacement, weighted refit, percentile",
            "n_requested": int(B), "seed": int(seed), "min_valid_fraction": float(min_valid_fraction),
            "measured_synthetic_coverage": list(MEASURED_GMR_COVERAGE if coverage is None else coverage["values"]),
            "estimand": ESTIMAND_MAGNITUDE}
    score = np.asarray(score, float)
    pt = fit_magnitude(score, group, cov, block)
    if pt["status"] != "ok":
        return _nonest("point_estimate_not_estimable", f"magnitude point estimate not estimable ({pt['reason_code']})", **meta)
    group = np.asarray(group).astype(int)
    block = np.asarray(block)
    y, X, act, _ = _design(score, group, cov, block)
    ub = np.unique(block)
    pos = {u: i for i, u in enumerate(ub)}
    blk_act = np.array([pos[v] for v in block[act]])
    rng = np.random.default_rng(seed)
    coefs = []
    for _ in range(B):
        mult = np.bincount(rng.integers(0, len(ub), len(ub)), minlength=len(ub)).astype(float)
        c = _weighted_group_coef(y, X, mult[blk_act])
        if c is not None:
            coefs.append(c)
    frac = len(coefs) / B
    meta.update({"n_valid": len(coefs), "valid_fraction": float(frac)})
    if frac < min_valid_fraction:
        return _nonest("bootstrap_valid_fraction", f"valid resamples {len(coefs)}/{B} = {frac:.3f} < {min_valid_fraction}", **meta)
    with np.errstate(over="ignore"):
        lo, hi = np.quantile(np.exp(coefs), [0.025, 0.975], method="linear")
    if not (math.isfinite(lo) and math.isfinite(hi)):
        return _nonest("nonfinite_estimate", "interval endpoints are not finite", **meta)
    return {"status": "ok", "reason": "", "reason_code": "", "lo": float(lo), "hi": float(hi), **meta}


def estimate_continuous(score, k_share, gate_group, cov, block) -> dict:
    """A05: effect per 1 SD of standardised log1p(k_share), POINT estimates only (no interval); gates on the primary binary label ``gate_group``."""
    from compartments import continuous_label

    z = continuous_label(k_share)
    return {"magnitude": fit_magnitude(score, z, cov, block, gate_group), "participation": fit_participation(score, z, cov, block, gate_group),
            "unit": "per 1 SD of standardised log1p(k_share)"}


def estimate_h1(score, group, cov, block, with_interval: bool = True, B: int = BOOT_N, seed: int = BOOT_SEED, min_valid_fraction: float = MIN_VALID_FRACTION,
                coverage: Optional[dict] = None) -> dict:
    mag = fit_magnitude(score, group, cov, block)
    if with_interval:
        mag["interval"] = gmr_interval(score, group, cov, block, B, seed, min_valid_fraction, coverage)
    return {"magnitude": mag, "participation": fit_participation(score, group, cov, block)}


# --------------------------------------------------------------------------- fail-closed output check
def assert_no_nan_ok(obj, _path="") -> None:
    """Every dict with status == 'ok' holds only finite numbers and carries a value (or an interval's lo/hi); NaN under not_estimable is allowed."""
    def finite_scan(o, p):
        if isinstance(o, dict):
            for k, v in o.items():
                finite_scan(v, f"{p}/{k}")
        elif isinstance(o, (list, tuple)):
            for i, v in enumerate(o):
                finite_scan(v, f"{p}[{i}]")
        elif isinstance(o, (float, np.floating)):
            assert math.isfinite(o), f"non-finite number at {p} under a status ok result"

    if isinstance(obj, dict):
        if obj.get("status") == "ok":
            assert ("value" in obj and obj["value"] is not None) or ("lo" in obj and "hi" in obj), f"ok result without a value at {_path}"
            for k, v in obj.items():
                if isinstance(v, dict) and v.get("status") == "not_estimable":
                    continue                           # a typed failure may hold NaN-like extras
                finite_scan(v, f"{_path}/{k}")
        for k, v in obj.items():
            assert_no_nan_ok(v, f"{_path}/{k}")
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            assert_no_nan_ok(v, f"{_path}[{i}]")


# --------------------------------------------------------------------------- claim-tree inputs (headline specimen only)
def _val(runs, run, part) -> Optional[float]:
    r = runs[run][part]
    return float(r["value"]) if r["status"] == "ok" else None


def claim_inputs(part: str, runs: Dict[str, dict], specimen: str) -> CT.Inputs:
    """Map estimate_h1 outputs of R00/R01/R02/R10/R11 to the tree inputs. A missing run is a KeyError; a not_estimable run is a typed False/not_estimable."""
    if specimen != HEADLINE_SPECIMEN:
        raise ScopeError(f"the claim tree is defined for {HEADLINE_SPECIMEN} only, not {specimen}")
    v = {r: _val(runs, r, part) for r in RUNS_FOR_CLAIM_TREE}
    est = v["R00"] is not None
    above = lambda x: x is not None and x > 1.0  # noqa: E731
    d = "not_estimable" if v["R10"] is None else ("agree" if above(v["R10"]) else "disagree")
    return CT.Inputs(est=est, s=above(v["R00"]), c=all(above(v[r]) for r in ("R00", "R01", "R02")), d=d, f=above(v["R11"]))


def decide_claim(part: str, runs: Dict[str, dict], specimen: str) -> dict:
    i = claim_inputs(part, runs, specimen)
    return {**CT.decide(i), "part": part, "inputs": i._asdict()}


# --------------------------------------------------------------------------- H2 rank (supplementary): GMR only
def rank_gmr(entries: Dict[str, dict], target: str) -> dict:
    """Competition ranking (ties share the lower rank) on the GMR among specimens with an estimable magnitude AND the absolute label rule; others are contextual."""
    for k, e in entries.items():
        if "odds_ratio" in e:
            raise ValueError(f"{k}: H2 ranks on the GMR only (no participation/odds ratio)")
    ranked = {k: e["value"] for k, e in entries.items() if e["status"] == "ok" and str(e["rule"]).startswith("absolute")}
    contextual = [k for k, e in entries.items() if e["status"] == "ok" and not str(e["rule"]).startswith("absolute")]
    not_est = [k for k, e in entries.items() if e["status"] != "ok"]
    ranks = {k: 1 + sum(1 for w in ranked.values() if w > v + 1e-12) for k, v in ranked.items()}
    out = {"ranks": dict(sorted(ranks.items(), key=lambda kv: (kv[1], kv[0]))), "k": len(ranked), "contextual": contextual, "not_estimable": not_est,
           "target_rank": ranks.get(target), "tie_rule": "competition ranking; values within 1e-12 share the lower rank"}
    out["wording"] = (f"{target} had the {ranks[target]}-th largest GMR of {len(ranked)} estimable specimens using the absolute rule (one ESR1-D538G tumour; rank only)"
                      if target in ranks else f"{target} is not ranked (not estimable or not using the absolute label rule)")
    return out


# --------------------------------------------------------------------------- H4 (supplementary): two statistics, ranked separately
def h4_statistics(np_mdk) -> dict:
    """Per-specimen H4 statistics from the N-P MDK values (pooled depth reference) of the specimen's ELIGIBLE cancer compartments:
    ``median`` = median of the values, ``fraction_nonzero`` = fraction of compartments with value STRICTLY > 0."""
    v = np.asarray(np_mdk, float)
    if v.size == 0 or not np.isfinite(v).all():
        return {"status": "not_estimable", "reason_code": "nonfinite_input", "reason": "no eligible compartments or non-finite N-P MDK values"}
    return {"status": "ok", "reason": "", "reason_code": "", "median": float(np.median(v)), "fraction_nonzero": float((v > 0).mean()), "n": int(v.size)}


def rank_h4(values: Dict[str, Optional[float]], target: str, reference_ids) -> dict:
    """ONE statistic, ONE reference set (frozen algorithm; spec 9.4.6 H4).

    * reference median  m = numpy.median of values[r] for r in ``reference_ids`` with a value (target excluded from the reference);
    * ranked population = reference specimens with a value plus the target; deviation d_i = |v_i - m|;
    * direction of the target = above / below / equal (equal iff |v_target - m| <= H4_TOL);
    * rank = 1 + #{ j : d_j > d_target + H4_TOL }  (competition ranking, 1 = largest deviation; a tie within H4_TOL shares the LOWER rank
      number; a specimen exactly at the median has d = 0 and ranks last-or-tied; no jitter, no secondary key);
    * a specimen without a value is not ranked and is listed; a target without a value has rank None.
    Rank only: no p-value, no interval, no combined rank across the two statistics or the two reference sets.
    """
    if target in reference_ids:
        raise ValueError("the target must not be part of its own reference set")
    ref = {r: values.get(r) for r in reference_ids}
    ref_ok = {r: float(v) for r, v in ref.items() if v is not None}
    unranked = sorted(r for r, v in ref.items() if v is None)
    tv = values.get(target)
    if not ref_ok:
        raise ValueError("empty reference set: H4 is not estimable without reference specimens")
    m = float(np.median(list(ref_ok.values())))
    out = {"reference_median": m, "n_reference": len(ref_ok), "unranked": unranked, "tie_tolerance": H4_TOL, "target": target}
    if tv is None:
        return {**out, "status": "not_estimable", "rank": None, "direction": None, "n_ranked": len(ref_ok), "tied_with": [],
                "wording": f"{target} has no estimable value for this statistic and is not ranked"}
    tv = float(tv)
    dev = {r: abs(v - m) for r, v in ref_ok.items()}
    dt = abs(tv - m)
    rank = 1 + sum(1 for d in dev.values() if d > dt + H4_TOL)
    tied = sorted(r for r, d in dev.items() if abs(d - dt) <= H4_TOL)
    direction = "equal" if abs(tv - m) <= H4_TOL else ("above" if tv > m else "below")
    n = len(ref_ok) + 1
    return {**out, "status": "ok", "value": tv, "deviation": dt, "direction": direction, "rank": rank, "n_ranked": n, "tied_with": tied,
            "wording": f"{target} had the {rank}-th largest absolute deviation of {n} specimens from the reference median ({direction} the reference median); rank only"}


def h4_all(stats: Dict[str, dict], target: str, wt_ids, surgical_wt_ids) -> dict:
    """The four frozen H4 results: {median, fraction_nonzero} x {all WT, surgical WT}. ``stats`` = specimen -> h4_statistics output."""
    out = {}
    for stat in ("median", "fraction_nonzero"):
        vals = {k: (v[stat] if v["status"] == "ok" else None) for k, v in stats.items()}
        out[stat] = {"all_wt": rank_h4(vals, target, list(wt_ids)), "surgical_wt": rank_h4(vals, target, list(surgical_wt_ids))}
    return out


# --------------------------------------------------------------------------- A01-A07 mechanics (frozen at rev 7; the sweep runner is built after the freeze)
_A_COMMON = {"specimens": [HEADLINE_SPECIMEN], "outcome": "R00 N-P MDK-NCL native COMMOT sender score (no COMMOT rerun; scores exist for every compartment)",
             "participation": "point estimate only", "claim_tree": False}
SWEEP_MECHANICS = {
    "common": _A_COMMON,
    "A01": {"changes": "eligibility min allocated counts 50", "recompute": ["eligibility", "labels (primary rule incl. its one fallback)", "blocks (make_blocks + median pitch on the new eligible set)", "covariate standardisation (new eligible set)"],
            "magnitude_interval": "nominal 95% production bootstrap", "label_gate": "primary rule", "fallback": "primary rule only"},
    "A02": {"changes": "eligibility min allocated counts 500", "recompute": ["eligibility", "labels (primary rule incl. its one fallback)", "blocks (make_blocks + median pitch on the new eligible set)", "covariate standardisation (new eligible set)"],
            "magnitude_interval": "nominal 95% production bootstrap", "label_gate": "primary rule", "fallback": "primary rule only"},
    "A03": {"changes": "k0 0.005", "recompute": ["labels"], "magnitude_interval": "nominal 95% production bootstrap", "label_gate": "30 eligible compartments per group", "fallback": "none"},
    "A04": {"changes": "k0 0.02", "recompute": ["labels"], "magnitude_interval": "nominal 95% production bootstrap", "label_gate": "30 eligible compartments per group", "fallback": "none"},
    "A05": {"changes": "continuous standardised log1p(k_share)", "recompute": ["label"], "magnitude_interval": "none (point estimates only)", "label_gate": "primary binary label", "fallback": "none"},
    "A06": {"changes": "top vs bottom tertile of k_share, middle dropped", "recompute": ["labels", "analysis subset"], "magnitude_interval": "nominal 95% production bootstrap",
            "label_gate": "30 compartments per tertile group on the subset", "fallback": "none",
            "note": "blocks keep the primary tile assignment (not re-merged on the subset); covariates keep the primary standardisation (never per subset)"},
    "A07": {"changes": "R_B: six-keratin allocated count >= 1", "recompute": ["labels"], "magnitude_interval": "nominal 95% production bootstrap", "label_gate": "30 eligible compartments per group", "fallback": "none"},
}
SWEEP_RULES = {"A01-A04,A06,A07": "same estimator, gates (30 compartments per group, 10 active senders per group, 10 mixed blocks) and typed failures as the primary run",
               "unchanged_unless_listed": "eligibility rows, tile assignment and covariate standardisation are those of the primary run"}


# --------------------------------------------------------------------------- runner
class MissingRunError(ValueError):
    """A required COMMOT run (R00/R01/R02/R10/R11 for the headline specimen) was not supplied at all."""


def frozen_settings(rna: Optional[str] = None) -> dict:
    """The settings the production entry uses (module constants; the protocol JSON must say the same, see ``assert_settings_frozen``).

    ``rna`` = 'per_cell' (the frozen protocol's covariates; the dict is exactly the pre-deviation one) or 'total' (protocol deviation 2026-11-06: ``covariates`` carries
    ``log_rna_total`` and ``rna_covariate`` + ``protocol_deviation`` (the deviation file's path and sha256) are added). Default: ``RNA_COVARIATE``."""
    rna = RNA_COVARIATE if rna is None else rna
    d = {"bootstrap_n": BOOT_N, "bootstrap_seed": BOOT_SEED, "min_valid_fraction": MIN_VALID_FRACTION, "level": LEVEL, "interval_runs": list(INTERVAL_RUNS),
         "tile_side_pitches": H.BLOCK_SIDE_PITCHES, "block_min_n": H.BLOCK_MIN_N, "pitch": "median nearest-neighbour distance of the eligible compartments' coordinates",
         "min_per_group": H.MIN_PER_GROUP, "min_active_per_group": H.MIN_ACTIVE_PER_GROUP, "min_mixed_blocks": H.MIN_MIXED_BLOCKS,
         "covariates": list(H.primary_covariates(rna))}
    if rna == "total":
        d["rna_covariate"] = "total"
        d["protocol_deviation"] = H.deviation_record()
    return d


def coverage_record(values, source: dict) -> dict:
    """The measured synthetic coverage a result carries: four values and where they come from (the stamped B24 artifact: run id, path, sha256, case order)."""
    values = [float(v) for v in values]
    if len(values) != len(MEASURED_GMR_COVERAGE) or not all(np.isfinite(values)):
        raise ValueError(f"coverage needs {len(MEASURED_GMR_COVERAGE)} finite values, got {values}")
    return {"values": values, "source": dict(source)}


def provenance_fields(settings: Optional[dict] = None, coverage: Optional[dict] = None) -> dict:
    """Provenance of a run: the settings ACTUALLY USED. ``production`` is derived by comparing them with the frozen settings (either covariate mode)."""
    st = frozen_settings() if settings is None else settings
    out = _provenance_fields(st, coverage)
    if coverage is not None:
        out["coverage_source"] = coverage["source"]
    return out


def _provenance_fields(st: dict, coverage: Optional[dict]) -> dict:
    return {
        "protocol_id": PROTOCOL_ID,
        "production": st == frozen_settings("per_cell") or st == frozen_settings("total"),
        "settings": st,
        "bootstrap": {"n": st["bootstrap_n"], "seed": st["bootstrap_seed"], "min_valid_fraction": st["min_valid_fraction"], "level": st["level"],
                      "tile_side_pitches": st["tile_side_pitches"], "widened": False},
        "covariates": list(st["covariates"]),
        "blocks": {"side_pitches": st["tile_side_pitches"], "min_n": st["block_min_n"]},
        "gates": {"min_per_group": st["min_per_group"], "min_active_per_group": st["min_active_per_group"], "min_mixed_blocks": st["min_mixed_blocks"]},
        "measured_synthetic_gmr_coverage": list(MEASURED_GMR_COVERAGE if coverage is None else coverage["values"]),
        "hypothesis_tests": "none",
    }


def _failed_runs(code: str, reason: str) -> dict:
    fail = {"status": "not_estimable", "reason_code": code, "reason": reason}
    return {"magnitude": {**fail, "interval": {**fail}}, "participation": dict(fail)}


def _required_runs(specimen: str) -> tuple:
    return RUNS_FOR_CLAIM_TREE if specimen == HEADLINE_SPECIMEN else ("R00",)


def _run_h1(specimen: str, obs, xy, scores: Dict[str, object], settings: dict, pitch: Optional[float], coverage: Optional[dict] = None) -> dict:
    if settings.get("rna_covariate", "per_cell") != RNA_COVARIATE:
        raise ValueError(f"settings record rna_covariate {settings.get('rna_covariate', 'per_cell')!r} but this runner computes {RNA_COVARIATE!r}")
    missing = [r for r in _required_runs(specimen) if r not in scores]
    if missing:
        raise MissingRunError(f"{specimen}: required run(s) {missing} not supplied (a failed COMMOT run must be passed as a typed result, not omitted)")
    xy = np.asarray(xy, float)
    n = len(obs)
    if len(xy) != n:
        raise ValueError("coordinates and eligible compartments must align")
    prov = provenance_fields(settings, coverage)
    base = {"specimen": specimen, "scope": "headline" if specimen == HEADLINE_SPECIMEN else "context_descriptive_only", "n_eligible": int(n), "provenance_fields": prov}
    if n < MIN_ELIGIBLE:                                   # early gate: before pitch, blocks, quantiles or any label
        why = f"{n} eligible compartments < {MIN_ELIGIBLE} (= 2 x {H.MIN_PER_GROUP}); no label can reach {H.MIN_PER_GROUP} per group"
        fail = {"status": "not_estimable", "reason_code": "group_size", "reason": why, "rule": None, "fallback_triggered": False}
        runs = {r: _failed_runs("group_size", why) for r in scores}
        claim = {part: decide_claim(part, runs, specimen) for part in ("magnitude", "participation")} if specimen == HEADLINE_SPECIMEN else None
        out = {**base, "n_blocks": None, "label": fail, "runs": runs, "claim_tree": claim}
        assert_no_nan_ok(out)
        return out
    blocks = H.make_blocks(xy, pitch if pitch is not None else H.median_pitch(xy), settings["tile_side_pitches"], settings["block_min_n"])      # coordinates only
    lab = select_label(obs["k_share"].to_numpy())
    label_out = {k: v for k, v in lab.items() if k != "high"}
    runs: Dict[str, dict] = {}
    for run, sc in scores.items():
        if isinstance(sc, dict):                           # a COMMOT run that failed or was not estimable
            runs[run] = _failed_runs("commot_run_failed", f"{run}: COMMOT run {sc.get('status', 'failed')}: {sc.get('reason', '')}")
            continue
        sc = np.asarray(sc, float)
        if len(sc) != n:
            raise ValueError(f"{run}: score vector has {len(sc)} entries for {n} eligible compartments")
        if lab["status"] != "ok":
            runs[run] = _failed_runs(lab["reason_code"], lab["reason"])
            continue
        try:
            cov = H.covariate_matrix(obs, nuclei_free=(run == "R11"), rna=RNA_COVARIATE)
        except H.CovariateDomainError as e:
            runs[run] = _failed_runs("covariate_domain", str(e))
            continue
        runs[run] = estimate_h1(sc, lab["high"].astype(int), cov, blocks, with_interval=run in settings["interval_runs"], B=settings["bootstrap_n"], seed=settings["bootstrap_seed"],
                                min_valid_fraction=settings["min_valid_fraction"], coverage=coverage)
    claim = {part: decide_claim(part, runs, specimen) for part in ("magnitude", "participation")} if specimen == HEADLINE_SPECIMEN else None
    out = {**base, "n_blocks": int(len(np.unique(blocks))), "label": label_out, "runs": runs, "claim_tree": claim}
    assert_no_nan_ok(out)
    return out


def run_h1_specimen(specimen: str, obs, xy, scores: Dict[str, object], coverage: Optional[dict] = None) -> dict:
    """PRODUCTION entry. H1 for one specimen at the FROZEN settings (no tunable argument). ``obs`` = the specimen's ELIGIBLE cancer compartments
    (with k_share etc.), ``xy`` their coordinates, ``scores[run]`` = the MDK-NCL sender score of that run aligned to ``obs`` OR a typed failure
    dict ({'status': 'failed'|'not_estimable', 'reason': ...}). The headline specimen needs R00/R01/R02/R10/R11 (``MissingRunError`` otherwise).
    One label (``select_label`` on k_share, chosen before any outcome is touched) and one block assignment serve every run; R11 swaps in the
    nuclei-free ``rna_per_cell_equiv_nf`` covariate and nothing else (RNA_COVARIATE 'total' since the 2026-11-06 deviation: the raw total has no nuclei term, so R11 differs
    from R00 only in its outcome). ``coverage`` = ``coverage_record`` of the stamped B24 run (``write_h1`` refuses a 'total' result without it). The claim tree is evaluated
    only for the headline specimen."""
    return _run_h1(specimen, obs, xy, scores, frozen_settings(), None, coverage)


def _run_h1_for_tests(specimen: str, obs, xy, scores: Dict[str, object], boot_B: int = BOOT_N, boot_seed: int = BOOT_SEED, pitch: Optional[float] = None,
                      interval_runs=INTERVAL_RUNS, coverage: Optional[dict] = None) -> dict:
    """TEST ONLY: same runner with overridable bootstrap size/seed/pitch/interval runs. The provenance records the settings used, so a result
    that departs from the frozen settings carries ``production: false`` and ``write_h1`` refuses it."""
    st = {**frozen_settings(), "bootstrap_n": int(boot_B), "bootstrap_seed": int(boot_seed), "interval_runs": list(interval_runs)}
    if pitch is not None:
        st["pitch"] = float(pitch)
    return _run_h1(specimen, obs, xy, scores, st, pitch, coverage)


def assert_settings_frozen(result: dict, proto: dict) -> None:
    """The settings recorded in ``result`` must equal the loaded protocol's frozen values AND the code constants (independent sources)."""
    import re

    f = result["provenance_fields"]
    st = f["settings"]
    b = proto["hurdle"]["inference"]["bootstrap"]
    gate = proto["hurdle"]["inference"]["estimability_gate"]
    est = proto["hurdle"]["estimable"]
    side = int(re.search(r"side (\d+) spot pitches", proto["hurdle"]["spatial_blocks"]["definition"]).group(1))
    mn = int(re.search(r"<(\d+) eligible", proto["hurdle"]["spatial_blocks"]["merge_rule"]).group(1))
    want = {"bootstrap_n": b["n"], "bootstrap_seed": b["seed"], "min_valid_fraction": b["min_valid_fraction"], "tile_side_pitches": float(side), "block_min_n": mn,
            "min_per_group": est["min_compartments_per_group"], "min_active_per_group": est["min_active_senders_per_group_for_magnitude"],
            "min_mixed_blocks": gate["min_mixed_label_blocks"], "interval_runs": ["R00"], "level": LEVEL}
    bad = {k: (st.get(k), v) for k, v in want.items() if st.get(k) != v}
    if bad:
        raise PIO.GateFailure(f"result settings differ from the frozen protocol (recorded, protocol): {bad}")
    _check_covariate_settings(st, proto)
    if st != frozen_settings(st.get("rna_covariate", "per_cell")) or not f.get("production") or f.get("protocol_id") != proto["protocol_id"]:
        raise PIO.GateFailure("result was not produced at the frozen production settings (pitch rule included) / protocol id")
    if st.get("rna_covariate") == "total" and not (f.get("coverage_source") or {}).get("sha256"):
        raise PIO.GateFailure("a rna_covariate='total' result must carry the stamped B24 coverage (provenance_fields.coverage_source with a sha256), not the per-cell constant")


def _check_covariate_settings(st: dict, proto: dict) -> None:
    """The covariates recorded in ``st`` must be the protocol's (``hurdle.models.covariates``) OR exactly the recorded deviation substitution (``log_rna_per_cell_equiv`` -> ``log_rna_total``,
    ``rna_covariate`` 'total', ``protocol_deviation`` = the committed deviation file's current path + sha256); nothing else."""
    import re

    want = [re.sub(r"^(\w+)\((\w+)\)$", r"\1_\2", c) for c in proto["hurdle"]["models"]["covariates"]]        # "log(spot_umi)" -> "log_spot_umi"
    if want != list(H.PRIMARY_COVARIATES):
        raise PIO.GateFailure(f"the protocol's covariates {want} differ from the code's frozen ones {list(H.PRIMARY_COVARIATES)}")
    rna = st.get("rna_covariate", "per_cell")
    if rna == "per_cell":
        if st.get("covariates") != want or "protocol_deviation" in st:
            raise PIO.GateFailure(f"result covariates {st.get('covariates')} differ from the frozen protocol {want} and no deviation is recorded")
        return
    sub = [("log_rna_total" if c == "log_rna_per_cell_equiv" else c) for c in want]
    if rna != "total" or st.get("covariates") != sub or st.get("protocol_deviation") != H.deviation_record():
        raise PIO.GateFailure(f"result covariates {st.get('covariates')} / rna_covariate {rna!r} / protocol_deviation {st.get('protocol_deviation')!r} are neither the protocol's "
                              f"covariates nor exactly the recorded substitution {sub} with the current deviation record")


def write_h1(out_dir, result: dict, produced_by: str, protocol_path, run_id: str) -> dict:
    """Write one specimen's production result: FROZEN protocol loaded and hash-verified (``protocol_io.load_protocol``), result settings checked against it,
    output namespaced ``.../<protocol sha12>/<run_id>/...`` into an EMPTY target directory (``check_output_target``), protocol sha256 and the
    provenance stamp added, atomic write. Raises on any mismatch; never writes a non-production result."""
    proto = PIO.load_protocol(protocol_path)                  # raises unless status FROZEN and the hash matches
    if proto.get("protocol_id") != PROTOCOL_ID:
        raise PIO.GateFailure(f"protocol_id {proto.get('protocol_id')!r} != code {PROTOCOL_ID!r}")
    assert_settings_frozen(result, proto)
    assert_no_nan_ok(result)
    target = PIO.check_output_target(Path(out_dir), proto["_protocol_sha12"], run_id)
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
    import provenance

    stamped = provenance.stamp_payload(PIO.with_protocol_hash(result, proto), produced_by)
    PIO.atomic_write_json(target / f"h1_{result['specimen']}.json", stamped)
    return stamped
