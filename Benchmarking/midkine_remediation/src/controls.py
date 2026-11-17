"""Descriptive real-data controls (protocol hurdle.validation_controls.real_data_post_freeze; spec 9.4.7). They never gate a branch.

Negative control (P4-S2, R00 / N-P, MDK-NCL magnitude, the 3 primary covariates, the frozen blocks, eligibility and labels):
  4000 ACCEPTED moving-block label shifts (``hurdle.shift_pool``, seed 20261015), then 999 of them selected WITHOUT replacement
  with a separately fixed seed (20261016). The magnitude GMR is refit under each shifted label (scores, covariates, blocks
  fixed). Reported beside the observed GMR: usable / failed counts (with typed reasons), min, 2.5%, 25%, 50%, 75%, 97.5%, max
  and every usable GMR. NO percentile of the observed value, p-value or threshold is computed.

Positive control: multiply the LINEAR allocated MDK counts of the ORIGINALLY defined eligible keratin-high group by 3.0
(``inject_mdk_layers``; labels, eligibility and the LR list are NOT recomputed), rerun COMMOT (envs/COMMOT) and the estimator,
and report baseline GMR, injected GMR and delta log(GMR) (``injection_report``).
"""
from __future__ import annotations

from collections import Counter
from typing import Dict

import numpy as np

import descriptive as D
import hurdle as H

NEG_POOL_N = 4000
NEG_POOL_SEED = 20261015
NEG_SELECT_N = 999
NEG_SELECT_SEED = 20261016
MIN_USABLE_SHIFTS = 900    # of 999 selected: below 90% usable the null spread is dropped (typed not_estimable). 900 keeps >= 22 shifts in each 2.5% tail of the reported
                           # quantiles (0.025 * 900 = 22.5); a stricter cut would discard runs whose spread is still interpretable, a laxer one would report tails of a few points.
INJECT_FOLD = 3.0
RNA_COVARIATE = "total"    # protocol deviation 2026-11-06: the covariate matrix handed to ``negative_control`` is built with this mode (a12_9_real / positive_control)


def negative_control(score, group, cov, block, xy, pitch: float, n_pool: int = NEG_POOL_N, n_select: int = NEG_SELECT_N,
                     pool_seed: int = NEG_POOL_SEED, select_seed: int = NEG_SELECT_SEED, pool=None, min_usable: int = MIN_USABLE_SHIFTS) -> dict:
    """Overall ``status``: ok iff the observed GMR is estimable AND at least ``min_usable`` of the ``n_select`` shifts gave a finite usable GMR;
    otherwise typed not_estimable (observed_not_estimable / too_few_usable_shifts) and the summary is withheld. Production uses the defaults
    (999 selected, 900 required); a smaller n_select / min_usable is for toy tests only."""
    if n_select > n_pool:
        raise ValueError("cannot select more shifts than the pool holds")
    group = np.asarray(group).astype(int)
    if pool is None:
        try:
            pool = H.shift_pool(xy, pitch, n_pool, pool_seed)
        except H.ShiftPoolUnfillable as e:                                      # data limitation, not a crash: typed not_estimable, no numbers
            return {"status": "not_estimable", "reason_code": "shift_pool_unfillable",
                    "reason": f"only {e.n_found} of the {e.n_requested} required acceptable label shifts exist for this specimen",
                    "n_shifts_found": e.n_found, "n_shifts_requested": e.n_requested,
                    "acceptance": e.acceptance, "seeds": {"pool": int(pool_seed), "select": int(select_seed)},
                    "n_requested": int(n_select), "gates_branch": False}
    sel = np.random.default_rng(select_seed).choice(len(pool), size=n_select, replace=False)
    gmrs, reasons = [], Counter()
    for s in sel:
        r = D.fit_magnitude(score, group[pool[s]], cov, block)
        if r["status"] == "ok" and np.isfinite(r["gmr"]) and r["gmr"] > 0:
            gmrs.append(r["gmr"])
        else:
            reasons[r["reason_code"] or "nonfinite_estimate"] += 1
    obs = D.fit_magnitude(score, group, cov, block)
    out = {
        "n_pool": int(len(pool)), "n_requested": int(n_select), "n_usable": len(gmrs), "n_failed": int(n_select - len(gmrs)),
        "failure_reasons": dict(reasons), "seeds": {"pool": int(pool_seed), "select": int(select_seed)},
        "selected_pool_indices": [int(i) for i in sel], "gmrs": [float(v) for v in gmrs], "observed": obs, "gates_branch": False,
        "frozen": {"specimen": D.HEADLINE_SPECIMEN, "run": "R00", "outcome": "MDK-NCL magnitude (GMR among active senders)",
                   "covariates": list(H.primary_covariates(RNA_COVARIATE)), "rna_covariate": RNA_COVARIATE, "blocks": "fixed", "eligibility": "fixed", "labels": "observed labels moved by shifts"},
    }
    if obs["status"] != "ok":
        out.update(status="not_estimable", reason_code="observed_not_estimable", reason=f"observed GMR not estimable ({obs.get('reason_code')})")
    elif len(gmrs) < min_usable:
        out.update(status="not_estimable", reason_code="too_few_usable_shifts", reason=f"{len(gmrs)} usable shifts < {min_usable} required of {n_select}")
    else:
        out.update(status="ok", reason_code="", reason="")
    out["min_usable_required"] = int(min_usable)
    if out["status"] == "ok":
        q = np.quantile(gmrs, [0.025, 0.25, 0.5, 0.75, 0.975], method="linear")
        out["value"] = float(q[2])                      # the null median, a descriptive anchor only (never compared with the observed value)
        out["summary"] = {"min": float(min(gmrs)), "q025": float(q[0]), "q25": float(q[1]), "median": float(q[2]), "q75": float(q[3]), "q975": float(q[4]), "max": float(max(gmrs))}
    return out


def inject_mdk_layers(layers: Dict[str, np.ndarray], genes, cancer_type: str, spot_idx, fold: float = INJECT_FOLD) -> Dict[str, np.ndarray]:
    """Copy of the SACE layers with the MDK column of ``cancer_type`` multiplied by ``fold`` in the rows ``spot_idx`` only."""
    if not fold > 1.0:
        raise ValueError("an injection needs fold > 1")
    genes = list(genes)
    if "MDK" not in genes:
        raise ValueError("no MDK gene in the layers")
    out = {t: v.copy() for t, v in layers.items()}
    out[cancer_type][np.asarray(spot_idx, int), genes.index("MDK")] *= fold
    return out


def injection_report(baseline: dict, injected: dict) -> dict:
    """Baseline GMR, injected GMR and delta log(GMR) (typed not_estimable if either run is)."""
    base = {"fold_applied": INJECT_FOLD, "gates_branch": False}
    if baseline["status"] != "ok":
        return {**base, "status": "not_estimable", "reason_code": "baseline_not_estimable", "reason": f"baseline run: {baseline.get('reason', '')}"}
    if injected["status"] != "ok":
        return {**base, "status": "not_estimable", "reason_code": "injected_not_estimable", "reason": f"injected run: {injected.get('reason', '')}"}
    b, i = float(baseline["gmr"]), float(injected["gmr"])
    return {**base, "status": "ok", "reason": "", "reason_code": "", "baseline_gmr": b, "injected_gmr": i, "delta_log_gmr": float(np.log(i) - np.log(b)), "value": float(np.log(i) - np.log(b))}
