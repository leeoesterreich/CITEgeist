"""B22 part 1 (check 9): is the depth adjustment right? Part 1 diagnosis (covariate correlations, leave-one-covariate-out) and Part 2 equal-depth test by
Poisson thinning to a common count target of 200 (no covariate regression). Descriptive only.

Promoted from the exploratory check 9 (``exploratory_v42/check9_depth.py``, git history). The analysis functions are that code, unchanged: ``spear``, ``thin``,
``log1p_sparse``, ``loo_subsets``, ``adjusted_rho``, ``agg_draws``, ``sign``, ``reading``, ``msr9``, ``ulm_sparse``, ``e2_sparse``, ``part1``, ``prop_table``, ``one_draw``,
``part2``, ``cohort9``, ``stored_originals``, ``three_way``, ``flat_part2`` (tests/midkine_remediation/test_moran_robustness.py checks the AST). Changed: ``load9`` takes the B5
compartments run as an argument (it was a module constant), the two stored-summary constants are gone (the stored check 7/8 values are the per-specimen dicts read from the
stamped B21 artifact), ``main`` became ``compute`` (inputs as arguments, no run directory of its own), and the table headers no longer print a dirty flag.
Rules: ``exploratory_v42/check9_DECISION_RULE.md`` (git history). The 500-count sensitivity run of the exploratory check is not reproduced (the rule fixes the target at 200).

Reused, not rewritten: keratin_e2 (the promoted check 7: E2 scores, odds_ratio, spearman_cols, bh), er_regulon (the promoted check 8: build_regulon, ulm_t, the committed regulon tables),
and what they reuse.
"""
from __future__ import annotations

import csv
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import a12_10_real as R  # noqa: E402
import compartments as C  # noqa: E402
import control_pairs as CP  # noqa: E402
import er_regulon as C8  # noqa: E402
import hurdle as H  # noqa: E402
import keratin_e2 as C7  # noqa: E402
import keratin_secretory_association as K  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

E2 = C7.E2
SEED9 = 20261109
TARGET = 200
N_DRAWS = 20
MIN_VALID = 15
B = 9999
N_CTRL = 200
MIN_EL, MIN_GROUP, MIN_NZ = C7.MIN_EL, C7.MIN_GROUP, C7.MIN_NZ
RNA_COVARIATE = "total"                       # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md)
COVS = tuple(H.PRIMARY_COVARIATES_TOTAL)      # log_spot_umi, log_rna_total, qp_prop
MDK = C.MDK
HEADLINE = "HCC22-088-P4-S2_1i_rep"
# test name -> (fixed, field, control vector, control transform). Conventions of checks 7/8: keratin or MDK is fixed, the second variable is the surrogated field
# (keratin tests: MDK is the field; MDK vs activity/score: the activity/score is the field).
TESTS = {
    "K1": ("kp", "mp", "kp", "present"), "K2": ("kp", "ml", "kp", "level"),
    "E1": ("ml", "act", "act", "level"), "E1p": ("mp", "act", "act", "present"),
    "P1": ("ml", "early", "early", "level"), "P1p": ("mp", "early", "early", "present"),
    "L1": ("ml", "late", "late", "level"), "L1p": ("mp", "late", "late", "present"),
}
# comparison -> (original-form source, equal-depth tests (level, present))
COMPARISONS = {"kpres_vs_mdk": ("c7_B", "K2", "K1"), "mdk_vs_esr1": ("c8", "E1", "E1p"),
               "mdk_vs_early_e2": ("c7_C_early", "P1", "P1p"), "mdk_vs_late_e2": ("c7_C_late", "L1", "L1p")}


# --------------------------------------------------------------------------- pure helpers (tested in test_check9.py)
def spear(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() < 1e-12 or b.std() < 1e-12 or not (np.isfinite(a).all() and np.isfinite(b).all()):
        return float("nan")
    return K.spearman(a, b)


def thin(E: sp.csr_matrix, tot, target: float, rng) -> sp.csr_matrix:
    """Poisson thinning of expected counts to ``target`` total: X' ~ Poisson(E * target / tot) per gene (approximation: E is non-integer)."""
    lam = (sp.diags(target / np.asarray(tot, float)) @ E).tocsr()
    out = lam.copy()
    out.data = rng.poisson(lam.data).astype(float)
    out.eliminate_zeros()
    return out


def log1p_sparse(X: sp.csr_matrix) -> sp.csr_matrix:
    out = X.copy()
    out.data = np.log1p(out.data)
    return out


def loo_subsets(names=COVS) -> dict:
    """Covariate subsets for the leave-one-out table: none (= primary), all (= adjusted), drop each, only each. Values = column indices of the covariate matrix."""
    n = len(names)
    d = {"none": (), "all": tuple(range(n))}
    for i, k in enumerate(names):
        d[f"drop_{k}"] = tuple(j for j in range(n) if j != i)
    for i, k in enumerate(names):
        d[f"only_{k}"] = (i,)
    return d


def adjusted_rho(a, b, Z, idx) -> float:
    """Spearman of a and b after both are residualised on the Z columns ``idx`` (none: raw Spearman)."""
    if len(idx) == 0:
        return spear(a, b)
    Zs = Z[:, list(idx)]
    return spear(K.residualize(a, Zs), K.residualize(b, Zs))


def agg_draws(draws: list) -> dict:
    """Aggregate per-draw test records (each ``{"valid": bool, "reason": str, ...}``) into the specimen-test record."""
    ok = [d for d in draws if d["valid"]]
    reasons: dict = {}
    for d in draws:
        if not d["valid"]:
            reasons[d["reason"]] = reasons.get(d["reason"], 0) + 1
    if len(ok) < MIN_VALID:
        return {"status": "not_estimable", "reason": f"valid draws {len(ok)} < {MIN_VALID} (too sparse)", "n_valid": len(ok), "invalid_reasons": reasons}
    rho = np.array([d["rho"] for d in ok])
    p = np.array([d["p"] for d in ok])
    out = {"status": "ok", "n_valid": len(ok), "invalid_reasons": reasons, "rho_mean": float(rho.mean()), "rho_min": float(rho.min()), "rho_max": float(rho.max()), "rho_sd": float(rho.std()),
           "p_combined": float(min(1.0, 2 * np.median(p))), "frac_draws_p_lt_0.05": float(np.mean(p < 0.05))}
    cp = [d["ctrl_pct"] for d in ok if d.get("ctrl_pct") is not None]
    if cp:
        out.update({"ctrl_pct_mean": float(np.mean(cp)), "n_draws_gt_p95": int(sum(d["gt95"] for d in ok if d.get("ctrl_pct") is not None)), "n_draws_with_ctrl": len(cp)})
    pr = [d["partial_qp"] for d in ok if d.get("partial_qp") is not None and np.isfinite(d["partial_qp"])]
    if pr:
        out["partial_qp_rho_mean"] = float(np.mean(pr))
    return out


def sign(x: float) -> int:
    return 0 if (x is None or not np.isfinite(x) or x == 0) else (1 if x > 0 else -1)


def reading(rows: list) -> dict:
    """Fixed reading rule (check9_DECISION_RULE.md). ``rows``: dicts with primary, adjusted, ed (equal-depth mean rho, level)."""
    disc = [r for r in rows if sign(r["primary"]) * sign(r["adjusted"]) < 0 and sign(r["ed"]) != 0]
    if len(disc) < 4:
        return {"n_discordant": len(disc), "verdict": "not decidable by sign (< 4 discordant specimens)", "n_ed_matches_primary": None, "n_ed_matches_adjusted": None}
    mp = sum(sign(r["ed"]) == sign(r["primary"]) for r in disc)
    ma = sum(sign(r["ed"]) == sign(r["adjusted"]) for r in disc)
    v = "adjusted over-corrects" if mp >= 2 * len(disc) / 3 else ("adjusted supported" if ma >= 2 * len(disc) / 3 else "mixed")
    return {"n_discordant": len(disc), "n_ed_matches_primary": int(mp), "n_ed_matches_adjusted": int(ma), "verdict": v}


# --------------------------------------------------------------------------- data
def load9(spec: str, target: float, comp_run) -> dict:
    import anndata as ad

    a = ad.read_h5ad(str(Path(comp_run) / f"{spec}.h5ad"), backed="r")
    obs_all = a.obs
    el = C.eligible_cancer(obs_all).to_numpy().astype(bool)
    rows = el & (obs_all["total_counts"].to_numpy(float) >= target)
    genes = list(a.var_names)
    X = R._dense(a.X[rows])
    Ed = R._dense(a.layers["E"][rows])
    xy = np.asarray(a.obsm["spatial"], float)[rows]
    obs = obs_all[rows]
    a.file.close()
    return {"obs": obs, "genes": genes, "X": X, "E": sp.csr_matrix(Ed), "xy": xy, "n": int(rows.sum()), "n_eligible": int(el.sum())}


def msr9(fixed, field, sn, label: str, B: int = B) -> dict:
    fixed, field = np.asarray(fixed, float), np.asarray(field, float)
    if not (np.isfinite(fixed).all() and np.isfinite(field).all()):
        return {"status": "not_estimable", "reason": "non-finite"}
    if fixed.std() < 1e-12 or field.std() < 1e-12:
        return {"status": "not_estimable", "reason": "constant variable"}
    seed = int(np.random.SeedSequence([SEED9, zlib.crc32(label.encode())]).generate_state(1)[0])
    r = K.msr_spearman(fixed, field, sn, B=B, seed=seed)
    return {"status": "ok", "rho": r["rho"], "p": r["p_two_sided"]}


def ulm_sparse(Xl: sp.csr_matrix, reg: dict, chunk: int = 1000) -> np.ndarray:
    out = np.empty(Xl.shape[0])
    for i in range(0, Xl.shape[0], chunk):
        out[i:i + chunk] = C8.ulm_t(Xl[i:i + chunk].toarray()[:, reg["feat_idx"]], reg["w"])
    return out


def e2_sparse(obs, genes, Xl) -> dict:
    import anndata as ad

    a = ad.AnnData(X=Xl.astype(np.float32), obs=pd.DataFrame({"spot_barcode": list(obs.index)}, index=list(obs.index)))
    a.var_names = pd.Index(genes)
    out = {}
    for ph in C7.PHASES:
        s = E2.compute_directional_score(a, ph)
        out[ph] = None if s is None else s.reindex(list(obs.index)).to_numpy(float)
    return out


# --------------------------------------------------------------------------- Part 1
def part1(spec, d, reg, Z, zerr, stored7, stored8) -> dict:
    obs, genes, X, E = d["obs"], d["genes"], d["X"], d["E"]
    n = d["n"]
    kc = obs["keratin_counts"].to_numpy(float)
    mdk_np = R._column(X, genes, MDK)
    mdk_e = np.asarray(E[:, genes.index(MDK)].todense()).ravel()
    act = C8.ulm_t(X[:, reg["feat_idx"]], reg["w"]) if reg["n_measured"] >= C8.MIN_TARGETS else np.full(n, np.nan)
    e2 = C7.e2_scores(obs, genes, X)
    sec = R._mean_z(X, genes, C8.SECRETORY)[0]
    var = {"mdk_np": mdk_np, "mdk_present": (mdk_e >= 1).astype(float), "keratin_presence": (kc >= 1).astype(float), "esr1_activity": act,
           "early_e2": e2["early"], "late_e2": e2["late"], "secretory": sec}
    var = {k: v for k, v in var.items() if v is not None}
    out = {"n": n, "covariate_error": zerr}
    if Z is None:
        return {**out, "status": "not_estimable", "reason": "covariate matrix unavailable: " + zerr}
    covs = {k: Z[:, i] for i, k in enumerate(COVS)}
    covs["log_total_counts"] = np.log(obs["total_counts"].to_numpy(float))
    out["corr_cov_outcome"] = {c: {o: spear(v, w) for o, w in var.items()} for c, v in covs.items()}
    out["corr_cov_cov"] = {c: {c2: spear(v, w) for c2, w in covs.items()} for c, v in covs.items()}
    pairs = {"a_kpres_vs_mdk": (var["keratin_presence"], mdk_np), "b_mdk_vs_esr1": (mdk_np, act)}
    if "early_e2" in var:
        pairs["c_mdk_vs_early"] = (mdk_np, var["early_e2"])
    if "late_e2" in var:
        pairs["d_mdk_vs_late"] = (mdk_np, var["late_e2"])
    out["loo"] = {pn: {sn_: adjusted_rho(a_, b_, Z, idx) for sn_, idx in loo_subsets().items()} for pn, (a_, b_) in pairs.items()}
    gates = []
    st = {"a_kpres_vs_mdk": C7._dig(stored7.get(spec, {}), "partB", "spearman_keratin_presence_vs_mdk_np"),
          "b_mdk_vs_esr1": C7._dig(stored8.get(spec, {}), "variants", "collectri", "main_mdk"),
          "c_mdk_vs_early": C7._dig(stored7.get(spec, {}), "partC", "early", "mdk_np"),
          "d_mdk_vs_late": C7._dig(stored7.get(spec, {}), "partC", "late", "mdk_np")}
    for pn, s in st.items():
        if not s or pn not in out["loo"]:
            continue
        for form, key in (("primary", "none"), ("adjusted", "all")):
            if s.get(form, {}).get("status") == "ok":
                gates.append(CM.gate_equal(f"{spec} {pn} {form} vs stored", out["loo"][pn][key], s[form]["rho"], rtol=1e-6, atol=1e-9))
    out["gates"] = gates
    out["status"] = "ok"
    return out


# --------------------------------------------------------------------------- Part 2
def prop_table(kp, mp) -> dict:
    n1, n0 = int(kp.sum()), int((~kp).sum())
    if n1 == 0 or n0 == 0:
        return {"n_present": n1, "n_absent": n0, "status": "one group empty"}
    m1, m0 = int(mp[kp].sum()), int(mp[~kp].sum())
    return {"n_present": n1, "n_absent": n0, "pct_mdk_present_keratin_present": 100.0 * m1 / n1, "pct_mdk_present_keratin_absent": 100.0 * m0 / n0, **C7.odds_ratio(m1, n1, m0, n0)}


def one_draw(spec, d, target, draw, reg, sn, ctrl_idx, B) -> dict:
    genes, E, obs = d["genes"], d["E"], d["obs"]
    n = d["n"]
    rng = np.random.default_rng(np.random.SeedSequence([SEED9, zlib.crc32(spec.encode()), int(target), draw]))
    tot = obs["total_counts"].to_numpy(float)
    Xc = thin(E, tot, target, rng)
    kidx = [genes.index(g) for g in C.KERATINS]
    kc = np.asarray(Xc[:, kidx].sum(axis=1)).ravel()
    mc = np.asarray(Xc[:, genes.index(MDK)].todense()).ravel()
    kp, mp = kc >= 1, mc >= 1
    Xl = log1p_sparse(Xc)
    act = ulm_sparse(Xl, reg) if reg["n_measured"] >= C8.MIN_TARGETS else np.full(n, np.nan)
    e2 = e2_sparse(obs, genes, Xl)
    V = {"kp": kp.astype(float), "mp": mp.astype(float), "ml": np.log1p(mc), "act": act, "early": e2["early"], "late": e2["late"]}
    cols = Xc[:, ctrl_idx].toarray() if ctrl_idx is not None else None
    qp = obs["qp_prop"].to_numpy(float)[:, None]
    mdk_nz, n1, n0 = int(mp.sum()), int(kp.sum()), int(n - kp.sum())
    rec = {"tests": {}, "sparsity": {"pct_mdk_ge1": 100 * float(mp.mean()), "pct_keratin_ge1": 100 * float(kp.mean()), "pct_both": 100 * float((mp & kp).mean()),
                                     "total_mean": float(np.asarray(Xc.sum(axis=1)).mean()), "total_sd": float(np.asarray(Xc.sum(axis=1)).std()), "median_genes_detected": float(np.median(np.diff(Xc.indptr)))},
           "table": prop_table(kp, mp)}
    for tn, (fx, fd, cv, tr) in TESTS.items():
        reason = ""
        if V[fd] is None or V[fx] is None or V[cv] is None:
            reason = "E2 phase not scorable"
        elif tn.startswith("K") and min(n1, n0) < MIN_GROUP:
            reason = f"keratin groups < {MIN_GROUP}"
        elif mdk_nz < MIN_NZ or mdk_nz == n:
            reason = f"MDK present in {mdk_nz} compartments (< {MIN_NZ} or all)"
        if not reason:
            r = msr9(V[fx], V[fd], sn, f"{spec}|{tn}|d{draw}|t{int(target)}", B)
            if r["status"] != "ok":
                reason = r["reason"]
        if reason:
            rec["tests"][tn] = {"valid": False, "reason": reason}
            continue
        t = {"valid": True, "rho": r["rho"], "p": r["p"], "ctrl_pct": None, "gt95": False}
        if cols is not None:
            cc = (cols >= 1).astype(float) if tr == "present" else np.log1p(cols)
            cr = C7.spearman_cols(V[cv], cc)
            cr = cr[np.isfinite(cr)]
            if len(cr):
                t.update({"ctrl_pct": C7.pct_below(r["rho"], cr), "gt95": bool(r["rho"] > np.percentile(cr, 95))})
        t["partial_qp"] = spear(K.residualize(V[fx], qp), K.residualize(V[fd], qp))
        rec["tests"][tn] = t
    return rec


def part2(spec, d, net, target, n_draws, B, nets) -> dict:
    genes, X, obs = d["genes"], d["X"], d["obs"]
    n = d["n"]
    res = {"specimen": spec, "n_rows": n, "n_eligible": d["n_eligible"], "target": target}
    if n < MIN_EL:
        return {**res, "status": "not_estimable", "reason": f"{n} compartments < {MIN_EL}"}
    reg = C8.build_regulon(net, genes, {MDK})
    sn = M.SpectralNull(knn_weights(d["xy"], K.K)) if n <= M.DENSE_EIGH_MAX_N else None
    if sn is None:
        return {**res, "status": "not_estimable", "reason": "n above dense-eigh ceiling"}
    ctrl_idx, ctrl_err = None, ""
    try:
        mdk_i = genes.index(MDK)
        dec, pool = CP.prepare_pool([X], genes, X, mdk_i, [])
        excl = set(E2.signature_union()) | set().union(*[set(n_[n_["source"] == C8.TF]["target"]) for n_ in nets.values()])
        pool = np.array([i for i in pool if genes[i] not in excl])
        ctrl_idx = CP.make_draws(genes, dec, pool, mdk_i, [], spec + "|check9", n_ctrl=N_CTRL)["a"]
    except ValueError as e:
        ctrl_err = str(e)
    draws = [one_draw(spec, d, target, k, reg, sn, ctrl_idx, B) for k in range(n_draws)]
    res["status"] = "ok"
    res["control_error"] = ctrl_err
    res["tests"] = {tn: agg_draws([dr["tests"][tn] for dr in draws]) for tn in TESTS}
    sp_keys = draws[0]["sparsity"].keys()
    res["sparsity_thinned"] = {k: float(np.mean([dr["sparsity"][k] for dr in draws])) for k in sp_keys}
    kc0 = obs["keratin_counts"].to_numpy(float)
    mc0 = np.asarray(d["E"][:, genes.index(MDK)].todense()).ravel()
    res["sparsity_original"] = {"pct_mdk_ge1": 100 * float((mc0 >= 1).mean()), "pct_keratin_ge1": 100 * float((kc0 >= 1).mean()), "pct_both": 100 * float(((mc0 >= 1) & (kc0 >= 1)).mean())}
    tabs = [dr["table"] for dr in draws if "pct_mdk_present_keratin_present" in dr["table"]]
    res["table_original"] = prop_table(kc0 >= 1, mc0 >= 1)
    res["table_thinned"] = {"n_draws_with_both_groups": len(tabs)}
    for k in ("n_present", "n_absent", "pct_mdk_present_keratin_present", "pct_mdk_present_keratin_absent", "odds_ratio"):
        v = [t[k] for t in tabs]
        if v:
            res["table_thinned"][k] = {"mean": float(np.mean(v)), "min": float(np.min(v)), "max": float(np.max(v))}
    res["table_thinned"]["n_draws_haldane"] = int(sum(t["haldane_anscombe"] for t in tabs))
    return res


# --------------------------------------------------------------------------- cohort, tables
def cohort9(per_spec: dict) -> dict:
    out = {}
    for tn in TESTS:
        ok = {s: r["tests"][tn] for s, r in per_spec.items() if r.get("status") == "ok" and r["tests"][tn]["status"] == "ok"}
        q = C7.bh([v["p_combined"] for v in ok.values()])
        for (s, v), qi in zip(ok.items(), q):
            v["q_bh"] = float(qi)
        rhos = [v["rho_mean"] for v in ok.values()]
        ranked = sorted(ok, key=lambda s: -ok[s]["rho_mean"])
        gt = [v for v in ok.values() if "n_draws_gt_p95" in v]
        out[tn] = {"n_estimable": len(ok), "n_positive": int(sum(r > 0 for r in rhos)), "n_q_lt_0.05": int(sum(v["q_bh"] < 0.05 for v in ok.values())),
                   "mean_rho": float(np.mean(rhos)) if rhos else None, "median_rho": float(np.median(rhos)) if rhos else None,
                   "n_specimens_gt_p95_in_majority_of_draws": int(sum(v["n_draws_gt_p95"] > v["n_valid"] / 2 for v in gt)),
                   "rank_by_rho_desc": {s: i + 1 for i, s in enumerate(ranked)}}
    return out


def stored_originals(stored7: dict, stored8: dict, spec: str, source: str):
    """(primary rho, adjusted rho, primary q, adjusted q) from the stored check 7/8 outputs."""
    if source == "c7_B":
        s = C7._dig(stored7.get(spec, {}), "partB", "spearman_keratin_presence_vs_mdk_np")
    elif source == "c8":
        s = C7._dig(stored8.get(spec, {}), "variants", "collectri", "main_mdk")
    else:
        s = C7._dig(stored7.get(spec, {}), "partC", source.split("_")[-1], "mdk_np")
    g = lambda f, k: (s[f][k] if s and s.get(f, {}).get("status") == "ok" and k in s[f] else None)  # noqa: E731
    return g("primary", "rho"), g("adjusted", "rho"), g("primary", "q_bh"), g("adjusted", "q_bh")


def three_way(per_spec: dict, stored7: dict, stored8: dict) -> dict:
    out = {}
    for cn, (src, tl, tp) in COMPARISONS.items():
        rows = []
        for s, r in per_spec.items():
            if r.get("status") != "ok":
                continue
            p, a, pq, aq = stored_originals(stored7, stored8, s, src)
            lv, pr = r["tests"][tl], r["tests"][tp]
            rows.append({"specimen": s, "n": r["n_rows"], "primary": p, "adjusted": a, "primary_q": pq, "adjusted_q": aq,
                         "ed": lv.get("rho_mean") if lv["status"] == "ok" else None, "ed_min": lv.get("rho_min"), "ed_max": lv.get("rho_max"), "ed_q": lv.get("q_bh"), "ed_n_valid": lv["n_valid"],
                         "ed_present": pr.get("rho_mean") if pr["status"] == "ok" else None, "ed_present_min": pr.get("rho_min"), "ed_present_max": pr.get("rho_max"), "ed_present_q": pr.get("q_bh"),
                         "ed_status": lv["status"] if lv["status"] != "ok" else "ok"})
        use = [x for x in rows if x["primary"] is not None and x["adjusted"] is not None and x["ed"] is not None]
        coh = {}
        for key in ("primary", "adjusted", "ed", "ed_present"):
            v = [x[key] for x in rows if x[key] is not None and x["primary"] is not None and x["adjusted"] is not None and x["ed"] is not None]
            coh[key] = {"n": len(v), "mean": float(np.mean(v)) if v else None, "median": float(np.median(v)) if v else None, "n_positive": int(sum(t > 0 for t in v))}
        coh["n_q_lt_0.05"] = {"primary": int(sum(x["primary_q"] is not None and x["primary_q"] < 0.05 for x in use)), "adjusted": int(sum(x["adjusted_q"] is not None and x["adjusted_q"] < 0.05 for x in use)),
                              "ed": int(sum(x["ed_q"] is not None and x["ed_q"] < 0.05 for x in use))}
        coh["ed_sign_matches_primary"] = int(sum(sign(x["ed"]) == sign(x["primary"]) for x in use))
        coh["ed_sign_matches_adjusted"] = int(sum(sign(x["ed"]) == sign(x["adjusted"]) for x in use))
        under = (coh["ed"]["mean"] is not None and abs(coh["ed"]["mean"]) < abs(coh["adjusted"]["mean"]) and coh["n_q_lt_0.05"]["ed"] < coh["n_q_lt_0.05"]["adjusted"])
        out[cn] = {"rows": rows, "cohort": coh, "reading": reading(use), "under_corrects_flag": bool(under)}
    return out


def f(x, nd=3):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def md_tables(meta, per_spec, coh, tw, p1) -> str:
    L = [f"# Check 9 generated tables (target {meta['target']}, draws {meta['n_draws']}, B {meta['n_surrogates']})", "",
         f"Output dir `{meta['out']}`; code commit {meta['git_commit'][:8]}.", ""]
    for cn, blk in tw.items():
        L += [f"## Three-way: {cn}", "", "| specimen | n | primary | adjusted | equal-depth level (mean [min,max]) | q | equal-depth present (mean [min,max]) | q | n valid |", "|---|---|---|---|---|---|---|---|---|"]
        for x in blk["rows"]:
            L.append(f"| {x['specimen']} | {x['n']} | {f(x['primary'])} | {f(x['adjusted'])} | {f(x['ed'])} [{f(x['ed_min'])}, {f(x['ed_max'])}] | {f(x['ed_q'], 4)} | {f(x['ed_present'])} [{f(x['ed_present_min'])}, {f(x['ed_present_max'])}] | {f(x['ed_present_q'], 4)} | {x['ed_n_valid']} ({x['ed_status']}) |")
        c = blk["cohort"]
        L += ["", f"Cohort (specimens with all three): primary mean {f(c['primary']['mean'])} median {f(c['primary']['median'])} n>0 {c['primary']['n_positive']}; adjusted mean {f(c['adjusted']['mean'])} median {f(c['adjusted']['median'])} n>0 {c['adjusted']['n_positive']}; "
                  f"equal-depth level mean {f(c['ed']['mean'])} median {f(c['ed']['median'])} n>0 {c['ed']['n_positive']}; equal-depth present mean {f(c['ed_present']['mean'])}. "
                  f"q<0.05: primary {c['n_q_lt_0.05']['primary']}, adjusted {c['n_q_lt_0.05']['adjusted']}, equal-depth {c['n_q_lt_0.05']['ed']}. Equal-depth sign matches primary in {c['ed_sign_matches_primary']}, adjusted in {c['ed_sign_matches_adjusted']} of {c['ed']['n']}.",
              f"Reading: {blk['reading']}; under-corrects flag {blk['under_corrects_flag']}.", ""]
    L += ["## Cohort lines per equal-depth test", ""]
    for tn, c in coh.items():
        L.append(f"- {tn}: n est {c['n_estimable']}, rho>0 {c['n_positive']}, q<0.05 {c['n_q_lt_0.05']}, mean {f(c['mean_rho'])}, median {f(c['median_rho'])}; real > control p95 in the majority of draws: {c['n_specimens_gt_p95_in_majority_of_draws']}")
    L += ["", "## Sparsity after thinning (mean over draws) vs original", "", "| specimen | MDK>=1 % thin (orig) | keratin>=1 % thin (orig) | both % thin (orig) | total mean (sd) | median genes detected | valid draws min over tests | not estimable tests |", "|---|---|---|---|---|---|---|---|"]
    for s, r in per_spec.items():
        if r.get("status") != "ok":
            L.append(f"| {s} | not estimable: {r.get('reason')} |||||||")
            continue
        t, o = r["sparsity_thinned"], r["sparsity_original"]
        ne = [tn for tn, v in r["tests"].items() if v["status"] != "ok"]
        L.append(f"| {s} | {f(t['pct_mdk_ge1'], 1)} ({f(o['pct_mdk_ge1'], 1)}) | {f(t['pct_keratin_ge1'], 1)} ({f(o['pct_keratin_ge1'], 1)}) | {f(t['pct_both'], 1)} ({f(o['pct_both'], 1)}) | {f(t['total_mean'], 1)} ({f(t['total_sd'], 1)}) | {f(t['median_genes_detected'], 0)} | {min(v['n_valid'] for v in r['tests'].values())} | {', '.join(ne) or '-'} |")
    h = HEADLINE
    if h in per_spec and per_spec[h].get("status") == "ok":
        r = per_spec[h]
        to, tt = r["table_original"], r["table_thinned"]
        L += ["", f"## {h} (P4-S2) exact", "", f"- original table: {to}", f"- thinned table (mean [min,max] over {tt['n_draws_with_both_groups']} draws): "
              + "; ".join(f"{k} {f(v['mean'], 3)} [{f(v['min'], 3)}, {f(v['max'], 3)}]" for k, v in tt.items() if isinstance(v, dict)) + f"; draws with a zero cell (Haldane-Anscombe): {tt['n_draws_haldane']}"]
        for tn, v in r["tests"].items():
            L.append(f"- {tn}: " + (f"rho mean {f(v['rho_mean'], 4)} [{f(v['rho_min'], 4)}, {f(v['rho_max'], 4)}], combined p {f(v['p_combined'], 4)}, q {f(v.get('q_bh'), 4)}, rank {coh[tn]['rank_by_rho_desc'].get(h)} of {coh[tn]['n_estimable']}, "
                                          f"ctrl pct mean {f(v.get('ctrl_pct_mean'), 1)}, draws > p95 {v.get('n_draws_gt_p95')}/{v['n_valid']}, partial(qp) {f(v.get('partial_qp_rho_mean'), 4)}" if v["status"] == "ok" else f"{v['status']}: {v['reason']}"))
    if p1:
        L += ["", "## Part 1: leave-one-covariate-out (Spearman after residualising both variables on the subset)", ""]
        for pn in ("a_kpres_vs_mdk", "b_mdk_vs_esr1", "c_mdk_vs_early", "d_mdk_vs_late"):
            names = list(loo_subsets())
            L += [f"### {pn}", "", "| specimen | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
            for s, r in p1.items():
                if r.get("status") == "ok" and pn in r["loo"]:
                    L.append(f"| {s} | " + " | ".join(f(r["loo"][pn][k]) for k in names) + " |")
            L.append("")
        L += ["### Spearman of each depth covariate with outcomes", ""]
        for s, r in p1.items():
            if r.get("status") != "ok":
                continue
            outs = list(next(iter(r["corr_cov_outcome"].values())))
            L += [f"{s} (n {r['n']})", "", "| covariate | " + " | ".join(outs) + " |", "|---|" + "---|" * len(outs)]
            for c, v in r["corr_cov_outcome"].items():
                L.append(f"| {c} | " + " | ".join(f(v[o]) for o in outs) + " |")
            cn_ = list(r["corr_cov_cov"])
            L += ["", "| cov x cov | " + " | ".join(cn_) + " |", "|---|" + "---|" * len(cn_)]
            for c, v in r["corr_cov_cov"].items():
                L.append(f"| {c} | " + " | ".join(f(v[o]) for o in cn_) + " |")
            L.append("")
    return "\n".join(L) + "\n"


def flat_part2(per_spec: dict) -> list:
    rows = []
    for s, r in per_spec.items():
        if r.get("status") != "ok":
            continue
        for tn, v in r["tests"].items():
            rows.append({"specimen": s, "n": r["n_rows"], "test": tn, "status": v["status"], "n_valid": v["n_valid"], **{k: v.get(k, "") for k in (
                "rho_mean", "rho_min", "rho_max", "rho_sd", "p_combined", "q_bh", "frac_draws_p_lt_0.05", "ctrl_pct_mean", "n_draws_gt_p95", "partial_qp_rho_mean")}})
    return rows


def compute(comp_run, stored7: dict, stored8: dict, specs, out, meta_common: dict, target: float = TARGET, n_draws: int = N_DRAWS, n_surrogates: int = B) -> tuple:
    """Check 9 (was ``main``): Part 1 + Part 2 on the B5 run ``comp_run``. ``stored7`` / ``stored8`` = per-specimen dicts of the stored check 7 / 8 results (from the stamped B21
    artifact); Part 1 is gated against them (fail closed: ``GateFailure``). Writes check9_part2_per_specimen.csv, check9_part1.csv and check9_tables.md into ``out``.
    Returns ``(summary, files)``."""
    nets = {v: pd.read_csv(C8.REG_DIR / fn, sep="\t") for v, (fn, _) in C8.VARIANTS.items() if v in ("collectri", "dorothea_ABC")}
    meta = {"check": "check9_depth", "seed": SEED9, "target": target, "n_draws": n_draws, "min_valid_draws": MIN_VALID, "n_surrogates": n_surrogates, "n_controls": N_CTRL, "k": K.K, **meta_common}
    p1, p2, gates = {}, {}, []
    for s in specs:
        t = time.time()
        d = load9(s, target, comp_run)
        if d["n"] >= MIN_EL:
            reg = C8.build_regulon(nets["collectri"], d["genes"], {MDK})
            try:
                Z, zerr = H.covariate_matrix(d["obs"], rna=RNA_COVARIATE), ""
            except H.CovariateDomainError as e:
                Z, zerr = None, str(e)
            p1[s] = part1(s, d, reg, Z, zerr, stored7, stored8)
            gates += p1[s].get("gates", [])
        else:
            p1[s] = {"status": "not_estimable", "reason": f"{d['n']} eligible < {MIN_EL}"}
        p2[s] = part2(s, d, nets["collectri"], target, n_draws, n_surrogates, nets)
        print(f"{s}: n={d['n']} part1={p1.get(s, {}).get('status')} part2={p2[s]['status']} {time.time()-t:.0f}s", flush=True)
        del d
    print(f"GATE PASSED ({len(gates)} comparisons with stored check 7/8 values)", flush=True)
    coh = cohort9(p2)
    tw = three_way(p2, stored7, stored8)
    summary = {**meta, "gates": gates, "cohort": coh, "three_way": tw, "part1": p1, "part2": p2}
    out = Path(out)
    rows = flat_part2(p2)
    with open(out / "check9_part2_per_specimen.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["specimen"])
        w.writeheader()
        w.writerows(rows)
    files = ["check9_part2_per_specimen.csv"]
    if p1:
        long = [{"specimen": s, "section": sec, "row": rr, "col": cc, "value": v} for s, r in p1.items() if r.get("status") == "ok"
                for sec in ("corr_cov_outcome", "corr_cov_cov", "loo") for rr, cd in r[sec].items() for cc, v in cd.items()]
        with open(out / "check9_part1.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["specimen", "section", "row", "col", "value"])
            w.writeheader()
            w.writerows(long)
        files.append("check9_part1.csv")
    (out / "check9_tables.md").write_text(md_tables(meta, p2, coh, tw, p1))
    files.append("check9_tables.md")
    return summary, files
