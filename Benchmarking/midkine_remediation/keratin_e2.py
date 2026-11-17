"""B21 part 1 (check 7): keratin presence, MDK, and the estrogen programme on eligible cancer compartments. Descriptive only.

Promoted from the exploratory check 7 (``exploratory_v42/check7_keratin_presence_e2.py``, git history). The analysis functions are that code, unchanged
(``presence``, ``odds_ratio``, ``pct_below``, ``resid_cols``, ``spearman_cols``, ``bh``, ``seed_for``, ``_nonest``, ``msr_pair``, ``both_forms``, ``control_block``,
``e2_scores``, ``_family``, ``_ctrl_family``, ``_dig``, ``cohort``, ``flat_rows``, ``gate_a07`` are AST-identical; tests/midkine_remediation/test_keratin_e2.py checks it).
Changed: ``load`` and ``analyze`` take the B5 compartments run as the keyword ``comp_run`` (it was a module constant); ``a07_reference`` reads the stamped B8 run it is
given. ``main`` is gone: the stage entry is ``e2_covariation.py``. Rules: ``exploratory_v42/check7_DECISION_RULE.md`` (git history).

Reused, not rewritten: keratin_secretory_association.{msr_spearman, residualize, spearman}, control_pairs.{deciles, make_draws, prepare_pool},
compartments.{eligible_cancer, rb_label}, hurdle.covariate_matrix, a12_9_inputs.load_score, a12_10_real.{_dense, _bh, _column},
morans.SpectralNull, spatial_utils.knn_weights, _e2_scoring.{compute_directional_score, signature_union} (code dependency, unchanged).
"""
from __future__ import annotations

import json
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))

import a12_9_inputs as A  # noqa: E402
import a12_10_real as R  # noqa: E402
import compartments as C  # noqa: E402
import control_pairs as CP  # noqa: E402
import hurdle as H  # noqa: E402
import keratin_secretory_association as K  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

sys.path.insert(0, str(HERE))
import _e2_scoring as E2  # noqa: E402

HEADLINE = "HCC22-088-P4-S2_1i_rep"
SEED = 20261107
B = 9999
N_CTRL = 200
MIN_EL = 30
MIN_GROUP = 30
MIN_NZ = 10
RNA_COVARIATE = "total"   # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md)
PHASES = ("early", "late")
THR_MAIN, THR_SENS = 1.0, 0.5


# --------------------------------------------------------------------------- pure helpers (tested in tests/midkine_remediation/test_keratin_e2.py)
def presence(counts, thr: float) -> np.ndarray:
    """Present iff count >= thr (the A07 rule at thr = 1)."""
    return np.asarray(counts, float) >= float(thr)


def odds_ratio(k1: int, n1: int, k0: int, n0: int) -> dict:
    """Raw odds ratio (present-group over absent-group) of an event; Haldane-Anscombe (+0.5) and flagged if any cell is 0."""
    cells = [k1, n1 - k1, k0, n0 - k0]
    hald = any(c == 0 for c in cells)
    a, b, c, d = (x + 0.5 for x in cells) if hald else cells
    return {"odds_ratio": float((a / b) / (c / d)), "haldane_anscombe": bool(hald), "cells": [int(x) for x in cells]}


def pct_below(real: float, ctrl) -> float:
    """Percentile of ``real`` among controls = 100 * mean(control < real) (as control_pairs)."""
    return float(100.0 * np.mean(np.asarray(ctrl, float) < real))


def resid_cols(Mx, Z) -> np.ndarray:
    """OLS residuals of every column of Mx on [1, Z]."""
    Xd = np.column_stack([np.ones(len(Mx)), Z])
    beta, *_ = np.linalg.lstsq(Xd, np.asarray(Mx, float), rcond=None)
    return np.asarray(Mx, float) - Xd @ beta


def spearman_cols(fixed, Mx) -> np.ndarray:
    """Spearman rho of ``fixed`` against each column of Mx."""
    zx = rankdata(fixed)
    zx = (zx - zx.mean()) / zx.std()
    Rk = rankdata(Mx, axis=0)
    Rc = Rk - Rk.mean(axis=0)
    sd = Rk.std(axis=0)
    out = np.full(Mx.shape[1], np.nan)
    ok = sd > 0
    out[ok] = (zx @ Rc[:, ok]) / (len(fixed) * sd[ok])
    return out


def bh(pvals):
    return R._bh(list(pvals)) if len(pvals) else []


def seed_for(label: str) -> int:
    return int(np.random.SeedSequence([SEED, zlib.crc32(label.encode())]).generate_state(1)[0])


def _nonest(reason: str) -> dict:
    return {"status": "not_estimable", "reason": reason}


def msr_pair(fixed, field, sn, label: str, B: int = B) -> dict:
    fixed, field = np.asarray(fixed, float), np.asarray(field, float)
    if fixed.std() < 1e-12 or field.std() < 1e-12:
        return _nonest("constant variable")
    r = K.msr_spearman(fixed, field, sn, B=B, seed=seed_for(label))
    return {"status": "ok", "rho": r["rho"], "p_two_sided": r["p_two_sided"], "n_surrogates": r["n_surrogates"]}


def both_forms(fixed, field, sn, Z, label: str, B: int = B) -> dict:
    """Primary and covariate-adjusted (both variables residualised on Z) Spearman + MSR."""
    out = {"primary": msr_pair(fixed, field, sn, label + "|primary", B)}
    if Z is None:
        out["adjusted"] = _nonest("covariate matrix unavailable")
    else:
        rf, rd = K.residualize(fixed, Z), K.residualize(field, Z)
        out["adjusted"] = msr_pair(rf, rd, sn, label + "|adjusted", B)
    return out


def control_block(real: dict, fixed, cols, Z) -> dict:
    """Real rho vs rho of control columns against ``fixed``; primary and adjusted forms. ``real`` = both_forms output."""
    out = {}
    for form in ("primary", "adjusted"):
        if real[form]["status"] != "ok" or (form == "adjusted" and Z is None):
            out[form] = _nonest("real statistic not estimable")
            continue
        if form == "primary":
            rho = spearman_cols(fixed, cols)
        else:
            rho = spearman_cols(K.residualize(fixed, Z), resid_cols(cols, Z))
        rho = rho[np.isfinite(rho)]
        r = real[form]["rho"]
        out[form] = {"status": "ok", "n_controls": int(len(rho)), "percentile": pct_below(r, rho), "control_p5": float(np.percentile(rho, 5)),
                     "control_median": float(np.median(rho)), "control_p95": float(np.percentile(rho, 95)), "real_gt_p95": bool(r > np.percentile(rho, 95)),
                     "real_lt_p5": bool(r < np.percentile(rho, 5))}
    return out


# --------------------------------------------------------------------------- data
def load(spec: str, *, comp_run) -> dict:
    import anndata as ad

    a = ad.read_h5ad(str(Path(comp_run) / f"{spec}.h5ad"), backed="r")
    el = C.eligible_cancer(a.obs).to_numpy().astype(bool)
    obs = a.obs[el]
    genes = list(a.var_names)
    X = R._dense(a.X[el])
    Emdk = R._dense(a.layers["E"][el])[:, genes.index(C.MDK)].astype(float)
    xy = np.asarray(a.obsm["spatial"], float)[el]
    a.file.close()
    return {"obs": obs, "genes": genes, "X": X, "E_mdk": Emdk, "xy": xy, "n": int(el.sum())}


def e2_scores(obs, genes, X) -> dict:
    import anndata as ad
    import pandas as pd

    ad_ = ad.AnnData(X=np.asarray(X, np.float32), obs=pd.DataFrame({"spot_barcode": list(obs.index)}, index=list(obs.index)))
    ad_.var_names = pd.Index(genes)
    out = {}
    for ph in PHASES:
        s = E2.compute_directional_score(ad_, ph)
        out[ph] = None if s is None else s.reindex(list(obs.index)).to_numpy(float)
    return out


# --------------------------------------------------------------------------- one specimen
def analyze(spec: str, commot_root, B: int = B, *, comp_run) -> dict:
    d = load(spec, comp_run=comp_run)
    n, obs, genes, X = d["n"], d["obs"], d["genes"], d["X"]
    kc = obs["keratin_counts"].to_numpy(float)
    res = {"specimen": spec, "n_eligible": int(n)}
    kp, kp05 = presence(kc, THR_MAIN), presence(kc, THR_SENS)
    mdk_np = R._column(X, genes, C.MDK)
    mp, mp05 = presence(d["E_mdk"], THR_MAIN), presence(d["E_mdk"], THR_SENS)
    res["partA"] = {"n_present": int(kp.sum()), "pct_present": float(100 * kp.mean()) if n else None,
                    "n_present_sens05": int(kp05.sum()), "pct_present_sens05": float(100 * kp05.mean()) if n else None}
    if n < MIN_EL:
        res["partB"] = _nonest(f"{n} eligible < {MIN_EL}")
        res["partC"] = _nonest(f"{n} eligible < {MIN_EL}")
        return res
    n1, n0 = int(kp.sum()), int(n - kp.sum())
    pc = n > M.DENSE_EIGH_MAX_N
    mdk_nz = int((mdk_np > 0).sum())
    try:
        Z, zerr = H.covariate_matrix(obs, rna=RNA_COVARIATE), ""
    except H.CovariateDomainError as e:
        Z, zerr = None, str(e)
    sc = A.load_score(commot_root / spec, "R00", list(obs.index))
    sender = None if isinstance(sc, dict) else np.asarray(sc, float)
    res["covariate_error"] = zerr
    res["sender_status"] = "ok" if sender is not None else sc.get("reason", "failed")
    sn = None if pc else M.SpectralNull(knn_weights(d["xy"], K.K))

    # ---- decile-matched controls (shared by Part B and C)
    ctrl_cols, ctrl_err = None, ""
    try:
        mdk_i = genes.index(C.MDK)
        dec, pool = CP.prepare_pool([X], genes, X, mdk_i, [])
        e2set = E2.signature_union()
        pool = np.array([i for i in pool if genes[i] not in e2set])
        dr = CP.make_draws(genes, dec, pool, mdk_i, [], spec + "|check7", n_ctrl=N_CTRL)
        ctrl_cols = X[:, dr["a"]]
        res["control_pool_size"] = int(len(pool))
    except ValueError as e:
        ctrl_err = str(e)

    # ---- Part B
    if min(n1, n0) < MIN_GROUP:
        res["partB"] = _nonest(f"keratin groups {n1}/{n0} < {MIN_GROUP} each")
    elif pc or mdk_nz < MIN_NZ or mdk_np.std() == 0:
        res["partB"] = _nonest(f"n={n} > ceiling or MDK>0 in {mdk_nz} < {MIN_NZ}")
    else:
        kf = kp.astype(float)
        b = {"n_present": n1, "n_absent": n0, "mdk_nonzero": mdk_nz}
        for thr_tag, kk, mm in (("", kp, mp), ("_sens05", kp05, mp05)):
            m1, m0 = int(mm[kk].sum()), int(mm[~kk].sum())
            g1, g0 = int(kk.sum()), int((~kk).sum())
            sec = {"n_present": g1, "n_absent": g0, "mdk_present_in_present": m1, "mdk_present_in_absent": m0}
            if min(g1, g0) >= MIN_GROUP:
                sec.update({"pct_mdk_present_keratin_present": 100.0 * m1 / g1, "pct_mdk_present_keratin_absent": 100.0 * m0 / g0, **odds_ratio(m1, g1, m0, g0)})
            else:
                sec["note"] = f"groups {g1}/{g0} < {MIN_GROUP}"
            if thr_tag == "":
                b["mdk_presence"] = sec
            else:
                b["sens05"] = sec
                if min(g1, g0) >= MIN_GROUP:
                    sec["spearman_raw_and_msr"] = msr_pair(kk.astype(float), mdk_np, sn, spec + "|B|sens05", B)
        b["mdk_np_present_group"] = {"median": float(np.median(mdk_np[kp])), "mean": float(mdk_np[kp].mean())}
        b["mdk_np_absent_group"] = {"median": float(np.median(mdk_np[~kp])), "mean": float(mdk_np[~kp].mean())}
        b["spearman_keratin_presence_vs_mdk_np"] = both_forms(kf, mdk_np, sn, Z, spec + "|B", B)
        if ctrl_cols is not None:
            b["matched_gene_control"] = control_block(b["spearman_keratin_presence_vs_mdk_np"], kf, ctrl_cols, Z)
        else:
            b["matched_gene_control"] = _nonest(ctrl_err)
        if sender is not None:
            act = sender != 0
            b["signalling_mdk_ncl"] = {
                "pct_active_keratin_present": float(100 * act[kp].mean()), "pct_active_keratin_absent": float(100 * act[~kp].mean()),
                "n_active_present": int(act[kp].sum()), "n_active_absent": int(act[~kp].sum()),
                "median_active_score_present": float(np.median(sender[kp & act])) if (kp & act).any() else None,
                "median_active_score_absent": float(np.median(sender[~kp & act])) if (~kp & act).any() else None}
        else:
            b["signalling_mdk_ncl"] = _nonest(res["sender_status"])
        res["partB"] = {"status": "ok", **b}

    # ---- Part C
    if pc or mdk_nz < MIN_NZ or mdk_np.std() == 0:
        res["partC"] = _nonest(f"n={n} > ceiling or MDK>0 in {mdk_nz} < {MIN_NZ}")
        return res
    scores = e2_scores(obs, genes, X)
    c = {"mdk_nonzero": mdk_nz}
    for ph in PHASES:
        s = scores[ph]
        if s is None:
            c[ph] = _nonest("E2 phase not scorable")
            continue
        blk = {"score_sd": float(s.std()), "mdk_np": both_forms(mdk_np, s, sn, Z, f"{spec}|C|{ph}|mdk", B)}
        blk["matched_gene_control"] = control_block(blk["mdk_np"], s, ctrl_cols, Z) if ctrl_cols is not None else _nonest(ctrl_err)
        blk["sender_mdk_ncl"] = both_forms(sender, s, sn, Z, f"{spec}|C|{ph}|sender", B) if sender is not None else _nonest(res["sender_status"])
        c[ph] = blk
    res["partC"] = {"status": "ok", **c}
    return res


# --------------------------------------------------------------------------- cohort
def _family(per_spec: dict, getter) -> dict:
    vals = {s: getter(r) for s, r in per_spec.items()}
    ok = {s: v for s, v in vals.items() if v and v.get("status") == "ok"}
    q = bh([v["p_two_sided"] for v in ok.values()])
    for (s, v), qi in zip(ok.items(), q):
        v["q_bh"] = float(qi)
    rhos = [v["rho"] for v in ok.values()]
    ranked = sorted(ok, key=lambda s: -ok[s]["rho"])
    return {"n_estimable": len(ok), "n_positive": sum(r > 0 for r in rhos), "n_q_lt_0.05": sum(v["q_bh"] < 0.05 for v in ok.values()),
            "mean_rho": float(np.mean(rhos)) if rhos else None, "median_rho": float(np.median(rhos)) if rhos else None,
            "rank_by_rho_desc": {s: i + 1 for i, s in enumerate(ranked)}}


def _ctrl_family(per_spec: dict, getter) -> dict:
    cs = [getter(r) for r in per_spec.values()]
    ok = [c for c in cs if c and c.get("status") == "ok"]
    return {"n_with_control": len(ok), "pct_specimens_real_gt_p95": float(100 * np.mean([c["real_gt_p95"] for c in ok])) if ok else None,
            "n_real_gt_p95": int(sum(c["real_gt_p95"] for c in ok)), "median_percentile": float(np.median([c["percentile"] for c in ok])) if ok else None}


def _dig(d, *keys):
    """Nested lookup; None if a key is missing or the node is a not_estimable record."""
    for k in keys:
        if not isinstance(d, dict) or d.get("status") == "not_estimable" or k not in d:
            return None
        d = d[k]
    return d if isinstance(d, dict) and d.get("status") != "not_estimable" else None


def cohort(per_spec: dict) -> dict:
    """Per-family cohort lines; BH q values are attached in place to the per-specimen records (the same dict objects)."""
    out = {}
    for ph in PHASES:
        for form in ("primary", "adjusted"):
            out[f"C_mdk_np|{ph}|{form}"] = _family(per_spec, lambda r, ph=ph, form=form: _dig(r, "partC", ph, "mdk_np", form))
            out[f"C_mdk_np|{ph}|{form}"]["control"] = _ctrl_family(per_spec, lambda r, ph=ph, form=form: _dig(r, "partC", ph, "matched_gene_control", form))
            out[f"C_sender|{ph}|{form}"] = _family(per_spec, lambda r, ph=ph, form=form: _dig(r, "partC", ph, "sender_mdk_ncl", form))
    for form in ("primary", "adjusted"):
        out[f"B_kpres_vs_mdk_np|{form}"] = _family(per_spec, lambda r, form=form: _dig(r, "partB", "spearman_keratin_presence_vs_mdk_np", form))
        out[f"B_kpres_vs_mdk_np|{form}"]["control"] = _ctrl_family(per_spec, lambda r, form=form: _dig(r, "partB", "matched_gene_control", form))
    return out


def flat_rows(per_spec: dict, genotype: dict) -> list:
    rows = []

    def add(spec, part, test, form, v, ctrl=None):
        if not isinstance(v, dict):
            return
        row = {"specimen": spec, "genotype": genotype.get(spec, ""), "n_eligible": per_spec[spec]["n_eligible"], "part": part, "test": test, "form": form,
               "status": v.get("status", ""), "rho": v.get("rho", ""), "p_two_sided": v.get("p_two_sided", ""), "q_bh": v.get("q_bh", "")}
        if isinstance(ctrl, dict) and ctrl.get("status") == "ok":
            row.update({"ctrl_percentile": ctrl["percentile"], "ctrl_p95": ctrl["control_p95"], "real_gt_p95": ctrl["real_gt_p95"]})
        rows.append(row)

    for s, r in per_spec.items():
        pb, pcx = r.get("partB"), r.get("partC")
        if isinstance(pb, dict) and pb.get("status") == "ok":
            for f in ("primary", "adjusted"):
                add(s, "B", "kpres_vs_mdk_np", f, pb["spearman_keratin_presence_vs_mdk_np"][f], pb["matched_gene_control"].get(f))
        if isinstance(pcx, dict) and pcx.get("status") == "ok":
            for ph in PHASES:
                x = pcx[ph]
                if x.get("status") == "not_estimable":
                    continue
                for f in ("primary", "adjusted"):
                    add(s, "C", f"mdk_np_vs_{ph}_e2", f, x["mdk_np"][f], x["matched_gene_control"].get(f))
                    if "status" not in x["sender_mdk_ncl"]:
                        add(s, "C", f"sender_vs_{ph}_e2", f, x["sender_mdk_ncl"][f])
    return rows


# --------------------------------------------------------------------------- A07 gate
def a07_reference(hurdle_run) -> dict:
    """A07 of the headline specimen from the stamped B8 run (``<hurdle_run>/<headline>/sweep/sweep_<headline>.json``; the story-level B8 json has no group sizes)."""
    p = Path(hurdle_run) / HEADLINE / "sweep" / f"sweep_{HEADLINE}.json"
    return json.loads(p.read_text())["variants"]["A07"]


def gate_a07(rec: dict, a07: dict) -> list:
    """P4-S2 keratin-present rule reproduces A07: group sizes and active-sender counts (fail-closed)."""
    b = rec["partB"]
    return [CM.gate_equal("A07 n_high", b["n_present"], a07["label"]["n_high"], rtol=0, atol=0),
            CM.gate_equal("A07 n_low", b["n_absent"], a07["label"]["n_low"], rtol=0, atol=0),
            CM.gate_equal("A07 n_active_high", b["signalling_mdk_ncl"]["n_active_present"], a07["magnitude"]["n_active_high"], rtol=0, atol=0),
            CM.gate_equal("A07 n_active_low", b["signalling_mdk_ncl"]["n_active_absent"], a07["magnitude"]["n_active_low"], rtol=0, atol=0)]
