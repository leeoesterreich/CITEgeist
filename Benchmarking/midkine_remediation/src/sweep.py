"""A01-A07 analysis-side sensitivity sweep (protocol ``commot.run_matrix.analysis_variant_mechanics``; ``descriptive.SWEEP_MECHANICS``).

Headline specimen only. Every variant works on the R00 (N-P) native COMMOT MDK-NCL sender score, which exists for EVERY compartment, so
there is no COMMOT rerun. Participation is a point estimate only; no claim tree. Inputs:

  * ``obs_all`` : all compartments of the specimen (any type; ``compartments.eligible_cancer`` picks the cancer ones) with the columns
                  ``cell_type, qp_prop, total_counts, k_share, keratin_counts, spot_umi, rna_per_cell_equiv, n_st_eff, depth_ratio`` (``build_compartments`` output;
                  the RNA covariate is ``hurdle.rna_total`` = rna_per_cell_equiv x n_st_eff x depth_ratio);
  * ``xy_all``  : their coordinates; ``score_all`` : the R00 sender score aligned to ``obs_all`` (or a typed failure dict).

What each variant changes (everything else is the primary run's: eligibility rows, tile assignment, covariate standardisation):
  A01/A02  allocated-count gate 50 / 500 -> eligibility, labels (primary rule incl. its one fallback), blocks (make_blocks with the median
           pitch of the NEW eligible set) and covariate standardisation are all recomputed over the new set.
  A03/A04  k0 0.005 / 0.02 on k_share: absolute rule, 30-per-group gate, NO fallback.
  A05      continuous standardised log1p(k_share): point estimates only (no interval); gates on the primary binary label.
  A06      top vs bottom tertile of k_share, middle dropped: fit on the tertile subset with the PRIMARY blocks (not re-merged) and the
           PRIMARY covariate standardisation (never recomputed per subset); gates evaluated on the subset.
  A07      R_B: six-keratin allocated count >= 1: label only, 30-per-group gate, NO fallback.
A typed failure is never replaced by another rule.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Optional

import numpy as np

import compartments as C
import descriptive as D
import hurdle as H
import protocol_io as PIO

RNA_COVARIATE = "total"                                        # protocol deviation 2026-11-06; must equal descriptive.RNA_COVARIATE (the runner refuses other settings)
VARIANTS = ("A01", "A02", "A03", "A04", "A05", "A06", "A07")
ELIGIBILITY_MIN_COUNTS = {"A01": 50.0, "A02": 500.0}          # primary: compartments.eligible_cancer default (200)
K0_VARIANT = {"A03": C.K0_SENSITIVITY[0], "A04": C.K0_SENSITIVITY[1]}
POINT_ONLY_VARIANTS = ("A05",)                                 # no magnitude interval
INTERVAL_VARIANTS = tuple(v for v in VARIANTS if v not in POINT_ONLY_VARIANTS)


def _failed(code: str, reason: str, label: Optional[dict] = None, n: Optional[int] = None) -> dict:
    out = {"label": label if label is not None else {"status": "not_estimable", "reason_code": code, "reason": reason}, "n_eligible": n,
           **D._failed_runs(code, reason)}
    if label is None:
        out["label"]["rule"] = None
    return out


def _strip(lab: dict) -> dict:
    return {k: v for k, v in lab.items() if k not in ("high", "high_mask", "low_mask")}


def _estimate(sc, group, cov, blocks, interval: bool, st: dict, coverage: Optional[dict] = None) -> dict:
    return D.estimate_h1(sc, group, cov, blocks, with_interval=interval, B=st["bootstrap_n"], seed=st["bootstrap_seed"], min_valid_fraction=st["min_valid_fraction"], coverage=coverage)


class _Primary:
    """Primary-run objects the variants inherit unless they say otherwise."""

    def __init__(self, obs_all, xy_all, score_all, st: dict, coverage: Optional[dict] = None):
        self.st = st
        self.coverage = coverage
        self.obs_all, self.xy_all = obs_all, np.asarray(xy_all, float)
        self.score_all = score_all
        if len(self.xy_all) != len(obs_all):
            raise ValueError("coordinates and compartments must align")
        if not isinstance(score_all, dict):
            self.score_all = np.asarray(score_all, float)
            if len(self.score_all) != len(obs_all):
                raise ValueError(f"score vector has {len(self.score_all)} entries for {len(obs_all)} compartments")
        el = C.eligible_cancer(obs_all).to_numpy()
        self.el = el
        self.obs, self.xy = obs_all[el], self.xy_all[el]
        self.sc = None if isinstance(score_all, dict) else self.score_all[el]
        self.n = int(el.sum())
        self.blocks = self.cov = self.lab = None
        self.cov_error: Optional[str] = None
        if self.n >= D.MIN_ELIGIBLE:
            self.blocks = H.make_blocks(self.xy, H.median_pitch(self.xy), st["tile_side_pitches"], st["block_min_n"])
            self.lab = C.select_label(self.obs["k_share"].to_numpy())
            try:
                self.cov = H.covariate_matrix(self.obs, rna=RNA_COVARIATE)
            except H.CovariateDomainError as e:
                self.cov_error = str(e)


def _run_variant(v: str, P: _Primary) -> dict:
    st = P.st
    if isinstance(P.score_all, dict):
        return _failed("commot_run_failed", f"R00 COMMOT run {P.score_all.get('status', 'failed')}: {P.score_all.get('reason', '')}")
    interval = v in INTERVAL_VARIANTS

    if v in ELIGIBILITY_MIN_COUNTS:                                        # A01 / A02: new eligible set, EVERYTHING recomputed
        el = C.eligible_cancer(P.obs_all, min_counts=ELIGIBILITY_MIN_COUNTS[v]).to_numpy()
        obs, xy, sc = P.obs_all[el], P.xy_all[el], P.score_all[el]
        n = int(el.sum())
        if n < D.MIN_ELIGIBLE:
            return _failed("group_size", f"{n} eligible compartments < {D.MIN_ELIGIBLE} (= 2 x {H.MIN_PER_GROUP}); no label can reach {H.MIN_PER_GROUP} per group", n=n)
        blocks = H.make_blocks(xy, H.median_pitch(xy), st["tile_side_pitches"], st["block_min_n"])
        lab = C.select_label(obs["k_share"].to_numpy())                    # primary rule, incl. its one fallback
        if lab["status"] != "ok":
            return _failed(lab["reason_code"], lab["reason"], _strip(lab), n)
        try:
            cov = H.covariate_matrix(obs, rna=RNA_COVARIATE)               # standardised over the NEW eligible set
        except H.CovariateDomainError as e:
            return _failed("covariate_domain", str(e), _strip(lab), n)
        return {"label": _strip(lab), "n_eligible": n, "n_blocks": int(len(np.unique(blocks))), **_estimate(sc, lab["high"].astype(int), cov, blocks, interval, st, P.coverage)}

    # A03-A07 inherit the primary eligible rows, blocks and covariates
    if P.n < D.MIN_ELIGIBLE:
        return _failed("group_size", f"{P.n} eligible compartments < {D.MIN_ELIGIBLE} (= 2 x {H.MIN_PER_GROUP}); no label can reach {H.MIN_PER_GROUP} per group", n=P.n)
    if P.cov is None:
        return _failed("covariate_domain", P.cov_error or "covariate matrix unavailable", n=P.n)
    k = P.obs["k_share"].to_numpy()
    base = {"n_eligible": P.n, "n_blocks": int(len(np.unique(P.blocks)))}

    if v in K0_VARIANT:                                                    # A03 / A04: absolute rule, gate, no fallback
        lab = C.select_label_absolute(k, K0_VARIANT[v])
        if lab["status"] != "ok":
            return _failed(lab["reason_code"], lab["reason"], _strip(lab), P.n)
        return {"label": _strip(lab), **base, **_estimate(P.sc, lab["high"].astype(int), P.cov, P.blocks, interval, st, P.coverage)}

    if v == "A05":                                                         # continuous; gates on the primary binary label
        lab = P.lab
        if lab["status"] != "ok":
            return _failed(lab["reason_code"], lab["reason"], _strip(lab), P.n)
        est = D.estimate_continuous(P.sc, k, lab["high"].astype(int), P.cov, P.blocks)
        return {"label": {**_strip(lab), "gate_label": "primary binary label", "continuous": "standardised log1p(k_share)"}, **base, **est}

    if v == "A06":                                                         # tertile subset; primary blocks and covariates
        tl = C.tertile_labels(k)
        if tl["status"] != "ok":
            return _failed(tl["reason_code"], tl["reason"], _strip(tl), P.n)
        sub = tl["high_mask"] | tl["low_mask"]
        return {"label": {**_strip(tl), "n_subset": int(sub.sum()), "value": float(tl["thresholds"][1])}, **base,
                **_estimate(P.sc[sub], tl["high_mask"][sub].astype(int), P.cov[sub], P.blocks[sub], interval, st, P.coverage)}

    if v == "A07":                                                         # R_B: keratin count >= 1, gate, no fallback
        high = C.rb_label(P.obs["keratin_counts"].to_numpy())
        lab = C._label_result(k, high, f"rb_keratin_count_ge_{C.RB_MIN_COUNT:g}", C.RB_MIN_COUNT, False)
        if lab["status"] != "ok":
            return _failed(lab["reason_code"], lab["reason"], _strip(lab), P.n)
        return {"label": _strip(lab), **base, **_estimate(P.sc, high.astype(int), P.cov, P.blocks, interval, st, P.coverage)}

    raise ValueError(f"unknown variant {v!r}")


def _run_sweep(specimen: str, obs_all, xy_all, score_all, settings: dict, variants=VARIANTS, coverage: Optional[dict] = None) -> dict:
    if settings.get("rna_covariate", "per_cell") != RNA_COVARIATE:
        raise ValueError(f"settings record rna_covariate {settings.get('rna_covariate', 'per_cell')!r} but this sweep computes {RNA_COVARIATE!r}")
    if specimen != D.HEADLINE_SPECIMEN:
        raise D.ScopeError(f"the A01-A07 sweep is defined for {D.HEADLINE_SPECIMEN} only, not {specimen}")
    bad = [v for v in variants if v not in VARIANTS]
    if bad:
        raise ValueError(f"unknown variant(s) {bad}")
    P = _Primary(obs_all, xy_all, score_all, settings, coverage)
    out = {"specimen": specimen, "scope": "headline", "run": "R00", "outcome": "MDK-NCL native COMMOT sender score (N-P), no COMMOT rerun",
           "participation": "point estimate only", "claim_tree": None, "provenance_fields": D.provenance_fields(settings, coverage),
           "variants": {v: {"changes": D.SWEEP_MECHANICS[v]["changes"], **_run_variant(v, P)} for v in variants}}
    D.assert_no_nan_ok(out)
    return out


def run_sweep(specimen: str, obs_all, xy_all, score_all, coverage: Optional[dict] = None) -> dict:
    """PRODUCTION entry: all seven variants at the frozen settings (2000 resamples, seed 20261014, 6-pitch tiles); no tunable argument other than ``coverage``
    (``descriptive.coverage_record`` of the stamped B24 run; ``write_sweep`` refuses a 'total' result without it)."""
    return _run_sweep(specimen, obs_all, xy_all, score_all, D.frozen_settings(), VARIANTS, coverage)


def _run_sweep_for_tests(specimen: str, obs_all, xy_all, score_all, boot_B: int = D.BOOT_N, boot_seed: int = D.BOOT_SEED, variants=VARIANTS, coverage: Optional[dict] = None) -> dict:
    """TEST ONLY: smaller bootstrap / a variant subset; provenance records the settings used so ``write_sweep`` refuses it."""
    st = {**D.frozen_settings(), "bootstrap_n": int(boot_B), "bootstrap_seed": int(boot_seed)}
    return _run_sweep(specimen, obs_all, xy_all, score_all, st, variants, coverage)


def write_sweep(out_dir, result: dict, produced_by: str, protocol_path, run_id: str) -> dict:
    """Same gate as ``descriptive.write_h1`` (protocol loaded and hash-verified, settings equal the frozen ones, empty namespaced target, atomic write)."""
    proto = PIO.load_protocol(protocol_path)
    if proto.get("protocol_id") != D.PROTOCOL_ID:
        raise PIO.GateFailure(f"protocol_id {proto.get('protocol_id')!r} != code {D.PROTOCOL_ID!r}")
    D.assert_settings_frozen(result, proto)
    if set(result["variants"]) != set(VARIANTS):
        raise PIO.GateFailure(f"a production sweep holds all of {VARIANTS}, got {sorted(result['variants'])}")
    D.assert_no_nan_ok(result)
    target = PIO.check_output_target(Path(out_dir), proto["_protocol_sha12"], run_id)
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
    import provenance

    stamped = provenance.stamp_payload(PIO.with_protocol_hash(result, proto), produced_by)
    PIO.atomic_write_json(target / f"sweep_{result['specimen']}.json", stamped)
    return stamped
