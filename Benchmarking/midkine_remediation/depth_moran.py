"""B22 part 2 (check 10): is the MDK-secretory bivariate Moran's I (B9), its matched controls (B17) and the keratin-secretory Spearman (B12) partly depth? Equal depth by
Poisson thinning to 200 counts, with a planted-association power check. Descriptive only.

Promoted from the exploratory check 10 (``exploratory_v42/check10_depth_moran.py``, git history). The analysis functions are that code, unchanged (``seed_for``, ``moran_pair``,
``b9_reference``, ``fidelity_gate``, ``mean_z_cols``, ``keratin_share``, ``pool_sparse``, ``control_I``, ``ctrl_block``, ``agg``, ``rec``, ``bh_attach``, ``noisefree_I``,
``calibrate_b``, ``plant_setup``, ``power_scenario``, ``make_ctx``, ``one_draw``, ``agg_controls``, ``cohort``, ``f``, ``agg_cell``, ``md_tables``, ``flat_rows``, ``write_csv``;
tests/midkine_remediation/test_moran_robustness.py checks the AST). Changed: ``run_specimen`` reads the B5 run from ``a.comp_run`` (it was a module constant), the B17 file
constant is gone (the stored B9 / B17 / B12 results are arguments), ``main`` became ``compute``, the table header no longer prints a dirty flag.
Rules: ``exploratory_v42/check10_DECISION_RULE.md`` (git history). The 500-count sensitivity run of the exploratory check is not reproduced (the rule fixes the target at 200).
"""
from __future__ import annotations

import csv
import sys
import time
import types
import zlib
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import a12_10_real as R  # noqa: E402
import compartments as C  # noqa: E402
import control_pairs as CP  # noqa: E402
import depth_thin as C9  # noqa: E402
import keratin_secretory_association as K  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

SEED10 = 20261110
TARGET = 200
N_DRAWS = 20
MIN_VALID = 15
B1 = 999               # Moran MSR surrogates, as B9
B3 = 9999              # Spearman MSR surrogates, as B12
N_CTRL = 200
MIN_EL, MIN_NZ = R.MIN_EL, R.MIN_NZ
MIN_KNZ = 30           # keratin share > 0 in >= 30 compartments (check 9 MIN_GROUP)
SEC = list(R.SECRETORY)
MDK = C.MDK
HEADLINE = "HCC22-088-P4-S2_1i_rep"
SCENARIOS = ("S0", "S1", "S2", "S3")
B_MAX, B_STEPS = 8.0, 40


# --------------------------------------------------------------------------- pure helpers (tested in test_check10.py)
def seed_for(label: str) -> int:
    return int(np.random.SeedSequence([SEED10, zlib.crc32(label.encode())]).generate_state(1)[0])


def moran_pair(x, y, W, sn, label: str, B: int = B1) -> dict:
    """B9 recipe: I = morans_bv(x, y, W); MSR p (two-sided) with x fixed and y surrogated. Typed not_estimable on non-finite or constant input."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        return {"status": "not_estimable", "reason": "non-finite"}
    if x.std() < 1e-12 or y.std() < 1e-12:
        return {"status": "not_estimable", "reason": "constant variable"}
    I = M.morans_bv(x, y, W)
    p = sn.pvalue(x, y, B=B, seed=seed_for(label))["p_two_sided"]
    return {"status": "ok", "I": float(I), "p": float(p)}


def b9_reference(X, genes, W, sn) -> dict:
    """Exactly what B9 (a12_10_real.run_analysis) computes for Definition C: x = MDK N_P, y = secretory score on N_P, B = 999, seed 42."""
    mdk = R._column(X, genes, MDK)
    y, nsec = R._mean_z(X, genes, SEC)
    return {"I": float(M.morans_bv(mdk, y, W)), "p_two_sided": float(sn.pvalue(mdk, y, B=R.B, seed=R.SEED)["p_two_sided"]), "n_secretory_genes": nsec}


def fidelity_gate(spec: str, got: dict, stored_row: dict) -> dict:
    """Fail closed: unthinned recompute must equal the stored B9 I to 1e-9 and the stored p exactly."""
    g = CM.gate_equal(f"{spec} B9 I (no thinning) vs stored", got["I"], stored_row["I"], rtol=1e-9, atol=1e-9)
    if got["p_two_sided"] != stored_row["p_two_sided"]:
        raise CM.GateFailure(f"GATE FAILED {spec} B9 p_two_sided (no thinning): got {got['p_two_sided']!r} expected {stored_row['p_two_sided']!r}")
    return {**g, "p_got": got["p_two_sided"], "p_expected": stored_row["p_two_sided"]}


def mean_z_cols(cols: np.ndarray) -> np.ndarray:
    """Per-column z (ddof 0, a constant column contributes 0), mean over columns: the B9 secretory score for an (n, g) matrix."""
    return R._mean_z_vecs([cols[:, j] for j in range(cols.shape[1])])


def keratin_share(Xc, kidx, mdk_i):
    """compartments.k_share on thinned counts: six keratins / (total - MDK), 0 when the denominator is 0. ``Xc`` is a CSC count matrix."""
    kc = np.asarray(Xc[:, kidx].sum(axis=1)).ravel()
    tot = np.asarray(Xc.sum(axis=1)).ravel()
    denom = tot - np.asarray(Xc[:, mdk_i].todense()).ravel()
    return np.where(denom > 0, kc / np.where(denom > 0, denom, 1.0), 0.0), kc, tot


def pool_sparse(Xl, genes) -> tuple:
    """(decile of each gene's mean over compartments, pool indices) on a thinned log1p matrix: pool = not in B17 EXCLUDE and not constant."""
    n = Xl.shape[0]
    mean = np.asarray(Xl.sum(axis=0)).ravel() / n
    var = np.asarray(Xl.multiply(Xl).sum(axis=0)).ravel() / n - mean ** 2
    dec = CP.deciles(np.round(mean, 12))        # means rounded to 12 decimals so that mathematically tied genes tie (sparse and dense summation differ in the last bit)
    pool = np.array([i for i, g in enumerate(genes) if g not in CP.EXCLUDE and var[i] > 1e-12])
    return dec, pool


def control_I(Xl_csc, mdk_l, sec_y, a_idx, b_idx, W) -> tuple:
    """Vectorised B17 controls: (I_a, I_b) = morans_bv of the control-a genes vs the real score, and of real MDK vs the control-b scores."""
    n = len(mdk_l)
    lag = W @ M.zscore(sec_y)
    ca = Xl_csc[:, np.asarray(a_idx)].toarray()
    Ia = ((ca - ca.mean(0)) / ca.std(0)).T @ lag / n
    nb, ns = b_idx.shape
    cb = Xl_csc[:, b_idx.ravel()].toarray()
    sd = cb.std(0)
    zb = np.where(sd > 0, (cb - cb.mean(0)) / np.where(sd > 0, sd, 1.0), 0.0)
    Y = zb.reshape(n, nb, ns).mean(axis=2)
    sdy = Y.std(0)
    Ib = np.full(nb, np.nan)
    ok = sdy > 0
    ZY = (Y[:, ok] - Y[:, ok].mean(0)) / sdy[ok]
    Ib[ok] = M.zscore(mdk_l) @ (W @ ZY) / n
    return Ia, Ib


def ctrl_block(real: float, I) -> dict:
    I = np.asarray(I, float)
    I = I[np.isfinite(I)]
    p95 = float(np.percentile(I, 95))
    return {"mean": float(I.mean()), "p5": float(np.percentile(I, 5)), "p95": p95, "pct": float(100.0 * np.mean(I < real)), "gt95": bool(real > p95), "n": int(len(I))}


def agg(draws: list) -> dict:
    """Aggregate per-draw records ``{"valid", "reason", "stat", "p"}`` into the specimen-test record (rule: >= MIN_VALID valid draws; combined p = min(1, 2 x median p))."""
    ok = [d for d in draws if d["valid"]]
    reasons = dict(Counter(d["reason"] for d in draws if not d["valid"]))
    if len(ok) < MIN_VALID:
        return {"status": "not_estimable", "reason": f"valid draws {len(ok)} < {MIN_VALID} (too sparse)", "n_valid": len(ok), "invalid_reasons": reasons}
    s = np.array([d["stat"] for d in ok])
    p = np.array([d["p"] for d in ok])
    return {"status": "ok", "n_valid": len(ok), "invalid_reasons": reasons, "mean": float(s.mean()), "min": float(s.min()), "max": float(s.max()), "sd": float(s.std()),
            "p_combined": float(min(1.0, 2 * np.median(p))), "frac_p_lt_0.05": float(np.mean(p < 0.05))}


def rec(r: dict, reason: str = "") -> dict:
    if reason:
        return {"valid": False, "reason": reason}
    if r["status"] != "ok":
        return {"valid": False, "reason": r["reason"]}
    return {"valid": True, "reason": "", "stat": r["I"] if "I" in r else r["rho"], "p": r["p"]}


def bh_attach(per_spec: dict, getter) -> None:
    """BH over the estimable specimens of one family; ``getter(res)`` returns the aggregate dict (q_bh is added in place)."""
    ok = [(s, getter(r)) for s, r in per_spec.items() if r.get("status") == "ok" and getter(r) and getter(r)["status"] == "ok"]
    q = R._bh([v["p_combined"] for _, v in ok]) if ok else []
    for (s, v), qi in zip(ok, q):
        v["q_bh"] = float(qi)


# --------------------------------------------------------------------------- planted association (power check)
def noisefree_I(R0, f, b, W) -> float:
    rate = R0 * np.exp(b * f)[:, None]
    return float(M.morans_bv(np.log1p(rate[:, 0]), mean_z_cols(np.log1p(rate[:, 1:])), W))


def calibrate_b(R0, f, target: float, W) -> dict:
    """Bisection of b in [0, B_MAX] so that the noise-free I equals ``target`` (target 0 -> b = 0). Monotonicity is assumed on the bracket; the achieved I is returned."""
    if target == 0:
        return {"b": 0.0, "I_noisefree": noisefree_I(R0, f, 0.0, W), "reachable": True}
    lo, hi = 0.0, B_MAX
    if noisefree_I(R0, f, lo, W) >= target:
        return {"b": 0.0, "I_noisefree": noisefree_I(R0, f, 0.0, W), "reachable": True, "note": "I(0) already >= target"}
    if noisefree_I(R0, f, hi, W) < target:
        return {"b": hi, "I_noisefree": noisefree_I(R0, f, hi, W), "reachable": False}
    for _ in range(B_STEPS):
        mid = 0.5 * (lo + hi)
        if noisefree_I(R0, f, mid, W) < target:
            lo = mid
        else:
            hi = mid
    b = 0.5 * (lo + hi)
    return {"b": float(b), "I_noisefree": noisefree_I(R0, f, b, W), "reachable": True}


def plant_setup(spec: str, E, tot, target: float, W, cols) -> tuple:
    """Permuted base rates (n, 12) and the latent field f (independent of depth by construction)."""
    n = E.shape[0]
    rng = np.random.default_rng(np.random.SeedSequence([SEED10, zlib.crc32((spec + "|plant").encode()), int(target)]))
    R0 = E[:, cols].toarray() * (target / np.asarray(tot, float))[:, None]
    for j in range(R0.shape[1]):
        R0[:, j] = R0[rng.permutation(n), j]
    e = rng.normal(size=n)
    f = W @ (W @ (W @ e))
    return R0, (f - f.mean()) / f.std()


def power_scenario(spec, R0, f, W, sn, target_I: float, scen: str, T: float, n_draws: int, B: int) -> dict:
    cal = calibrate_b(R0, f, target_I, W)
    if not cal["reachable"]:
        return {"status": "not_estimable", "reason": "target not reachable by b = %g" % B_MAX, "target": target_I, **cal}
    rate = R0 * np.exp(cal["b"] * f)[:, None]
    n = R0.shape[0]
    draws = []
    for k in range(n_draws):
        rng = np.random.default_rng(np.random.SeedSequence([SEED10, zlib.crc32((spec + "|plantdraw").encode()), int(T), SCENARIOS.index(scen), k]))
        cnt = rng.poisson(rate)
        nz = int((cnt[:, 0] >= 1).sum())
        if nz < MIN_NZ or nz == n:
            draws.append(rec({}, f"MDK present in {nz} compartments (< {MIN_NZ} or all)"))
            continue
        draws.append(rec(moran_pair(np.log1p(cnt[:, 0]), mean_z_cols(np.log1p(cnt[:, 1:])), W, sn, f"{spec}|plant{scen}|d{k}|t{int(T)}", B)))
    return {"target": target_I, **cal, **agg(draws)}


# --------------------------------------------------------------------------- one thinned draw
def make_ctx(spec: str, d: dict) -> dict:
    genes = list(d["genes"])
    W = knn_weights(np.asarray(d["xy"], float), R.K)
    return {"spec": spec, "genes": genes, "n": d["n"], "W": W, "sn": M.SpectralNull(W), "mdk_i": genes.index(MDK), "sec_idx": [genes.index(g) for g in SEC],
            "kidx": [genes.index(g) for g in C.KERATINS], "depth": np.log(d["obs"]["spot_umi"].to_numpy(float))}


def one_draw(d: dict, ctx: dict, target: float, k: int, B1_: int = B1, B3_: int = B3, n_ctrl: int = N_CTRL, controls: bool = True) -> dict:
    spec, n, W, sn, genes = ctx["spec"], ctx["n"], ctx["W"], ctx["sn"], ctx["genes"]
    rng = np.random.default_rng(np.random.SeedSequence([SEED10, zlib.crc32(spec.encode()), int(target), k]))
    Xc = C9.thin(d["E"], d["obs"]["total_counts"].to_numpy(float), target, rng).tocsc()
    Xl = Xc.copy()
    Xl.data = np.log1p(Xl.data)
    mc = np.asarray(Xc[:, ctx["mdk_i"]].todense()).ravel()
    mdk_l = np.log1p(mc)
    sec_cnt = Xc[:, ctx["sec_idx"]].toarray()
    sec = mean_z_cols(np.log1p(sec_cnt))
    ks, kc, tot = keratin_share(Xc, ctx["kidx"], ctx["mdk_i"])
    nz = int((mc >= 1).sum())
    t = f"t{int(target)}"
    out = {"detection": {"pct_mdk_ge1": 100.0 * float((mc >= 1).mean()), "pct_keratin_ge1": 100.0 * float((kc >= 1).mean()), "pct_kshare_gt0": 100.0 * float((ks > 0).mean()),
                         "pct_secretory_gene_ge1_mean": 100.0 * float((sec_cnt >= 1).mean()), "total_mean": float(tot.mean()), "total_sd": float(tot.std())}, "tests": {}}
    mdk_bad = f"MDK present in {nz} compartments (< {MIN_NZ} or all)" if (nz < MIN_NZ or nz == n) else ""
    out["tests"]["m1"] = rec(moran_pair(mdk_l, sec, W, sn, f"{spec}|m1|d{k}|{t}", B1_), mdk_bad) if not mdk_bad else rec({}, mdk_bad)
    out["tests"]["d_mdk"] = rec({}, mdk_bad) if mdk_bad else rec(moran_pair(ctx["depth"], mdk_l, W, sn, f"{spec}|dmdk|d{k}|{t}", B1_))
    out["tests"]["d_sec"] = rec(moran_pair(ctx["depth"], sec, W, sn, f"{spec}|dsec|d{k}|{t}", B1_))
    nk = int((ks > 0).sum())
    if nk < MIN_KNZ:
        out["tests"]["k3"] = rec({}, f"keratin share > 0 in {nk} compartments (< {MIN_KNZ})")
    elif sec.std() < 1e-12:
        out["tests"]["k3"] = rec({}, "constant secretory score")
    else:
        r = K.msr_spearman(ks, sec, sn, B=B3_, seed=seed_for(f"{spec}|k3|d{k}|{t}"))
        out["tests"]["k3"] = {"valid": True, "reason": "", "stat": r["rho"], "p": r["p_two_sided"]}
    out["control"] = None
    if controls and out["tests"]["m1"]["valid"]:
        real = out["tests"]["m1"]["stat"]
        try:
            dec, pool = pool_sparse(Xl, genes)
            dr = CP.make_draws(genes, dec, pool, ctx["mdk_i"], ctx["sec_idx"], f"{spec}|check10|{t}|d{k}", n_ctrl=n_ctrl, n_uni=1, seed=SEED10)
            Ia, Ib = control_I(Xl, mdk_l, sec, dr["a"], dr["b"], W)
            out["control"] = {"real": real, "a": ctrl_block(real, Ia), "b": ctrl_block(real, Ib), "pool_size": int(len(pool))}
        except ValueError as e:
            out["control"] = {"error": str(e)}
    return out


def agg_controls(draws: list) -> dict:
    ok = [d["control"] for d in draws if d["control"] and "error" not in d["control"]]
    errs = dict(Counter(d["control"]["error"] for d in draws if d["control"] and "error" in d["control"]))
    if len(ok) < MIN_VALID:
        return {"status": "not_estimable", "n_valid": len(ok), "reason": f"valid control draws {len(ok)} < {MIN_VALID}", "errors": errs}
    out = {"status": "ok", "n_valid": len(ok), "errors": errs, "real_mean": float(np.mean([c["real"] for c in ok]))}
    for tag in ("a", "b"):
        n_gt = int(sum(c[tag]["gt95"] for c in ok))
        out[tag] = {"mean": float(np.mean([c[tag]["mean"] for c in ok])), "p5": float(np.mean([c[tag]["p5"] for c in ok])), "p95": float(np.mean([c[tag]["p95"] for c in ok])),
                    "pct_mean": float(np.mean([c[tag]["pct"] for c in ok])), "n_draws_gt_p95": n_gt, "specimen_gt_p95": bool(n_gt > len(ok) / 2)}
    return out


# --------------------------------------------------------------------------- one specimen
def run_specimen(spec: str, a, b9: dict, gate_on: bool) -> dict:
    t0 = time.time()
    d = C9.load9(spec, a.target, a.comp_run)
    n = d["n"]
    res = {"specimen": spec, "n": n, "n_eligible_all": d["n_eligible"], "target": a.target}
    if n < MIN_EL:
        return {**res, "status": "not_estimable", "reason": f"{n} compartments < {MIN_EL}"}
    if n > M.DENSE_EIGH_MAX_N:
        return {**res, "status": "not_estimable", "reason": "n above dense-eigh ceiling"}
    ctx = make_ctx(spec, d)
    X, genes, obs = d["X"], ctx["genes"], d["obs"]
    ref = b9_reference(X, genes, ctx["W"], ctx["sn"])
    res["unthinned_b9_recipe"] = ref
    stored = b9["per_specimen"].get(spec, {})
    if gate_on and stored.get("status") == "ok":
        if stored["n_eligible"] != n:
            raise CM.GateFailure(f"GATE FAILED {spec}: n {n} != stored n_eligible {stored['n_eligible']}")
        res["gate"] = fidelity_gate(spec, ref, stored)
    mdk_np = R._column(X, genes, MDK)
    sec_np = R._mean_z(X, genes, SEC)[0]
    res["unthinned"] = {"d_mdk": moran_pair(ctx["depth"], mdk_np, ctx["W"], ctx["sn"], f"{spec}|dmdk_unthin|t{int(a.target)}", a.b1),
                        "d_sec": moran_pair(ctx["depth"], sec_np, ctx["W"], ctx["sn"], f"{spec}|dsec_unthin|t{int(a.target)}", a.b1)}
    # unthinned detection / original keratin share
    kc0 = obs["keratin_counts"].to_numpy(float)
    mc0 = np.asarray(d["E"][:, ctx["mdk_i"]].todense()).ravel()
    res["detection_original"] = {"pct_mdk_ge1": 100.0 * float((mc0 >= 1).mean()), "pct_keratin_ge1": 100.0 * float((kc0 >= 1).mean()), "pct_kshare_gt0": 100.0 * float((obs["k_share"].to_numpy(float) > 0).mean()),
                                 "pct_mdk_np_gt0": 100.0 * float((mdk_np > 0).mean())}
    draws = [one_draw(d, ctx, a.target, k, a.b1, a.b3, a.n_ctrl) for k in range(a.n_draws)]
    res["status"] = "ok"
    res["tests"] = {tn: agg([dr["tests"][tn] for dr in draws]) for tn in ("m1", "d_mdk", "d_sec", "k3")}
    res["controls"] = agg_controls(draws)
    res["detection_thinned"] = {k: float(np.mean([dr["detection"][k] for dr in draws])) for k in draws[0]["detection"]}
    res["detection_thinned_range"] = {k: [float(min(dr["detection"][k] for dr in draws)), float(max(dr["detection"][k] for dr in draws))] for k in ("pct_mdk_ge1", "pct_keratin_ge1", "pct_kshare_gt0")}
    # power check
    cols = [ctx["mdk_i"]] + ctx["sec_idx"]
    R0, f = plant_setup(spec, d["E"], obs["total_counts"].to_numpy(float), a.target, ctx["W"], cols)
    own = stored.get("I") if stored.get("status") == "ok" else None
    targets = {"S0": 0.0, "S1": b9["aggregate"]["mean"], "S2": 2 * b9["aggregate"]["mean"], "S3": own if (own is not None and own > 0) else None}
    res["power"] = {"depth_spearman_of_f": C9.spear(f, ctx["depth"])}
    for sc in SCENARIOS:
        res["power"][sc] = power_scenario(spec, R0, f, ctx["W"], ctx["sn"], targets[sc], sc, a.target, a.n_draws, a.b1) if targets[sc] is not None else {"status": "not_applicable", "reason": "stored I <= 0 or not estimable"}
    print(f"{spec}: n={n} m1={res['tests']['m1']['status']} k3={res['tests']['k3']['status']} ctrl={res['controls']['status']} {time.time()-t0:.0f}s", flush=True)
    del d
    return res


# --------------------------------------------------------------------------- cohort, reading rule, tables
def cohort(per_spec: dict, b9: dict, b17: dict, b12: dict, gate_on: bool) -> dict:
    for tn in ("m1", "d_mdk", "d_sec", "k3"):
        bh_attach(per_spec, lambda r, tn=tn: r["tests"][tn])
    for tn in ("d_mdk", "d_sec"):
        ok = {s: r["unthinned"][tn] for s, r in per_spec.items() if r.get("status") == "ok" and r["unthinned"][tn]["status"] == "ok"}
        q = R._bh([v["p"] for v in ok.values()]) if ok else []
        for (s, v), qi in zip(ok.items(), q):
            v["q_bh"] = float(qi)
    for sc in SCENARIOS:
        bh_attach(per_spec, lambda r, sc=sc: r["power"][sc])
    est = {s: r for s, r in per_spec.items() if r.get("status") == "ok"}
    out = {"stored": {"b9_mean_I": b9["aggregate"]["mean"], "b9_sig": b9["bh"]["significant_q_lt_0.05"], "b17_cohort": b17["cohort"]["citegeist_alone_definition_C"]}}
    # R1
    r1 = {}
    for s in b9["bh"]["significant_q_lt_0.05"]:
        v = est.get(s, {}).get("tests", {}).get("m1")
        if v is None or v["status"] != "ok":
            r1[s] = "not estimable"
        else:
            r1[s] = "remains q<0.05" if v["q_bh"] < 0.05 else "lost"
    out["R1"] = {"per_specimen": r1, "n_remaining": sum(x == "remains q<0.05" for x in r1.values()), "n_stored_sig": len(r1)}
    # R2/R3
    both = [s for s, r in est.items() if r["tests"]["m1"]["status"] == "ok" and r["controls"]["status"] == "ok"]
    if both:
        out["R2"] = {"n_specimens": len(both), "mean_equal_depth_I": float(np.mean([est[s]["tests"]["m1"]["mean"] for s in both])), "mean_ctrl_a_mean": float(np.mean([est[s]["controls"]["a"]["mean"] for s in both])),
                     "mean_ctrl_b_mean": float(np.mean([est[s]["controls"]["b"]["mean"] for s in both])), "mean_ctrl_a_p95": float(np.mean([est[s]["controls"]["a"]["p95"] for s in both])),
                     "mean_ctrl_b_p95": float(np.mean([est[s]["controls"]["b"]["p95"] for s in both])), "mean_real_in_control_draws": float(np.mean([est[s]["controls"]["real_mean"] for s in both]))}
        out["R2"]["above_controls"] = bool(out["R2"]["mean_equal_depth_I"] > out["R2"]["mean_ctrl_a_mean"] and out["R2"]["mean_equal_depth_I"] > out["R2"]["mean_ctrl_b_mean"])
        out["R3"] = {"n_specimens": len(both), "n_gt_p95_a": int(sum(est[s]["controls"]["a"]["specimen_gt_p95"] for s in both)), "n_gt_p95_b": int(sum(est[s]["controls"]["b"]["specimen_gt_p95"] for s in both))}
    # R4
    r4 = {}
    for s, st in out["R1"]["per_specimen"].items():
        if st == "lost":
            r4[s] = {sc: (est[s]["power"][sc].get("p_combined", 1.0) < 0.05 if est[s]["power"][sc].get("status") == "ok" else None) for sc in ("S1", "S2", "S3")}
            r4[s]["verdict_S1"] = "informative" if r4[s]["S1"] else "not informative (thinning could not see an association of that size)"
    out["R4"] = {"lost": r4, "n_informative_S1": int(sum(bool(v["S1"]) for v in r4.values())), "n_not_informative_S1": int(sum(not v["S1"] for v in r4.values()))}
    out["power_cohort"] = {sc: {"n_estimable": int(sum(r["power"][sc].get("status") == "ok" for r in est.values())), "n_combined_p_lt_0.05": int(sum(r["power"][sc].get("status") == "ok" and r["power"][sc]["p_combined"] < 0.05 for r in est.values())),
                                "n_q_lt_0.05": int(sum(r["power"][sc].get("status") == "ok" and r["power"][sc].get("q_bh", 1) < 0.05 for r in est.values()))} for sc in SCENARIOS}
    # R5
    r5 = {}
    for tn in ("d_mdk", "d_sec"):
        un = [r["unthinned"][tn] for r in est.values() if r["unthinned"][tn]["status"] == "ok"]
        th = [r["tests"][tn] for r in est.values() if r["tests"][tn]["status"] == "ok"]
        r5[tn] = {"unthinned_n_q_lt_0.05": int(sum(u.get("q_bh", 1) < 0.05 for u in un)), "unthinned_mean_I": float(np.mean([u["I"] for u in un])) if un else None, "unthinned_n_positive": int(sum(u["I"] > 0 for u in un)),
                  "thinned_n_q_lt_0.05": int(sum(t.get("q_bh", 1) < 0.05 for t in th)), "thinned_mean_I": float(np.mean([t["mean"] for t in th])) if th else None, "thinned_n_estimable": len(th)}
    r5["n_both_unthinned_positive"] = int(sum(r["unthinned"]["d_mdk"].get("I", -1) > 0 and r["unthinned"]["d_sec"].get("I", -1) > 0 for r in est.values()))
    out["R5"] = r5
    # R6
    r6 = {}
    ks = {s: r for s, r in est.items() if r["tests"]["k3"]["status"] == "ok"}
    for form in ("primary", "partial_primary"):
        rows = {s: b12["per_specimen"][form].get(s) for s in ks}
        rows = {s: v for s, v in rows.items() if v and v.get("status") == "ok"}
        sgn = lambda x: 0 if x == 0 else (1 if x > 0 else -1)  # noqa: E731
        r6[form] = {"n": len(rows), "stored_mean_rho": float(np.mean([v["rho"] for v in rows.values()])) if rows else None, "stored_n_q_lt_0.05": int(sum(v["q"] < 0.05 for v in rows.values())),
                    "equal_depth_mean_rho": float(np.mean([ks[s]["tests"]["k3"]["mean"] for s in rows])) if rows else None, "equal_depth_median_rho": float(np.median([ks[s]["tests"]["k3"]["mean"] for s in rows])) if rows else None,
                    "equal_depth_n_q_lt_0.05": int(sum(ks[s]["tests"]["k3"]["q_bh"] < 0.05 for s in rows)),
                    "sign_agree": int(sum(sgn(ks[s]["tests"]["k3"]["mean"]) == sgn(v["rho"]) for s, v in rows.items()))}
    r6["n_estimable_equal_depth"] = len(ks)
    r6["not_estimable"] = {s: r["tests"]["k3"].get("reason") for s, r in est.items() if r["tests"]["k3"]["status"] != "ok"}
    out["R6"] = r6
    return out


def f(x, nd=3):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def agg_cell(v, nd=3):
    return "not est. (%s valid)" % v.get("n_valid", 0) if v["status"] != "ok" else f"{f(v['mean'], nd)} [{f(v['min'], nd)}, {f(v['max'], nd)}]"


def md_tables(meta, per_spec, coh, b9, b17, b12) -> str:
    est = {s: r for s, r in per_spec.items()}
    L = [f"# Check 10 generated tables (target {meta['target']}, draws {meta['n_draws']}, B1 {meta['b1']}, B3 {meta['b3']}, controls {meta['n_ctrl']})", "",
         f"Output dir `{meta['out']}`; code commit {meta['git_commit'][:8]}; fidelity gate {'PASSED for ' + str(meta['n_gates']) + ' specimens' if meta['gate_on'] else 'not applicable (subset target)'}.", ""]
    L += ["## Part 1: Moran's I, MDK vs secretory (B9 recipe)", "", "| specimen | n | stored B9 I | stored q | unthinned recompute I | equal-depth I mean [min, max] | combined p | frac draws p<0.05 | BH q | valid draws |", "|---|---|---|---|---|---|---|---|---|---|"]
    for s, r in est.items():
        if r.get("status") != "ok":
            L.append(f"| {s} | {r.get('n')} | not estimable: {r.get('reason')} |||||||")
            continue
        st = b9["per_specimen"].get(s, {})
        v = r["tests"]["m1"]
        qs = b9["bh"].get(s, {}).get("q")
        L.append(f"| {s} | {r['n']} | {f(st.get('I'), 4)} | {f(qs, 4)} | {f(r['unthinned_b9_recipe']['I'], 4)} | {agg_cell(v, 4)} | {f(v.get('p_combined'), 4)} | {f(v.get('frac_p_lt_0.05'), 2)} | {f(v.get('q_bh'), 4)} | {v['n_valid']} |")
    r1 = coh["R1"]
    L += ["", f"R1: B9 q<0.05 specimens and their equal-depth status: {r1['per_specimen']}; remaining {r1['n_remaining']} of {r1['n_stored_sig']}.", ""]
    L += ["## Part 1b: Moran's I of log_spot_umi with MDK and with the secretory score", "", "| specimen | unthinned I(depth, MDK N_P) | p | q | thinned mean [min, max] | q | unthinned I(depth, secretory N_P) | p | q | thinned mean [min, max] | q |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for s, r in est.items():
        if r.get("status") != "ok":
            continue
        um, us, tm, ts = r["unthinned"]["d_mdk"], r["unthinned"]["d_sec"], r["tests"]["d_mdk"], r["tests"]["d_sec"]
        L.append(f"| {s} | {f(um.get('I'), 4)} | {f(um.get('p'), 3)} | {f(um.get('q_bh'), 4)} | {agg_cell(tm, 4)} | {f(tm.get('q_bh'), 4)} | {f(us.get('I'), 4)} | {f(us.get('p'), 3)} | {f(us.get('q_bh'), 4)} | {agg_cell(ts, 4)} | {f(ts.get('q_bh'), 4)} |")
    L += ["", f"R5: {coh['R5']}", ""]
    L += ["## Part 2: matched controls at equal depth (B17 shape); stored B17 beside", "", "| specimen | equal-depth real I | ctrl a mean / p95 / real pct / draws>p95 | ctrl b mean / p95 / real pct / draws>p95 | specimen >p95 (a, b) | stored real I | stored a mean / p95 / gt95 | stored b mean / p95 / gt95 | valid ctrl draws |", "|---|---|---|---|---|---|---|---|---|"]
    for s, r in est.items():
        if r.get("status") != "ok":
            continue
        c = r["controls"]
        sb = b17["per_specimen"].get(s, {})
        sa, sbb = sb.get("control_a_random_mdk", {}), sb.get("control_b_random_secretory", {})
        if c["status"] != "ok":
            L.append(f"| {s} | not est.: {c.get('reason')} |||||||| ")
            continue
        L.append(f"| {s} | {f(c['real_mean'], 4)} | {f(c['a']['mean'], 4)} / {f(c['a']['p95'], 4)} / {f(c['a']['pct_mean'], 1)} / {c['a']['n_draws_gt_p95']} | {f(c['b']['mean'], 4)} / {f(c['b']['p95'], 4)} / {f(c['b']['pct_mean'], 1)} / {c['b']['n_draws_gt_p95']} | "
                 f"{c['a']['specimen_gt_p95']}, {c['b']['specimen_gt_p95']} | {f(sb.get('real_I'), 4)} | {f(sa.get('mean'), 4)} / {f(sa.get('p95'), 4)} / {sa.get('real_gt_p95')} | {f(sbb.get('mean'), 4)} / {f(sbb.get('p95'), 4)} / {sbb.get('real_gt_p95')} | {c['n_valid']} |")
    L += ["", f"R2: {coh.get('R2')}", f"R3: {coh.get('R3')}", f"Stored B17 cohort: {coh['stored']['b17_cohort']}", ""]
    L += ["## Power check (planted association in the real lattice and sparsity)", ""]
    for sc in SCENARIOS:
        L += [f"### {sc}", "", "| specimen | target I | b | noise-free I | realised I mean [min, max] | frac draws p<0.05 | combined p | BH q | valid |", "|---|---|---|---|---|---|---|---|---|"]
        for s, r in est.items():
            if r.get("status") != "ok":
                continue
            v = r["power"][sc]
            if v.get("status") == "ok":
                L.append(f"| {s} | {f(v['target'], 4)} | {f(v['b'], 3)} | {f(v['I_noisefree'], 4)} | {agg_cell(v, 4)} | {f(v['frac_p_lt_0.05'], 2)} | {f(v['p_combined'], 4)} | {f(v.get('q_bh'), 4)} | {v['n_valid']} |")
            else:
                L.append(f"| {s} | {f(v.get('target'), 4)} | | | {v.get('status')}: {v.get('reason')} (valid {v.get('n_valid', '')}) |||||")
        L.append("")
    L += [f"Power cohort: {coh['power_cohort']}", f"R4: {coh['R4']}", "",
          "Correlation of the planted latent field with log_spot_umi (Spearman): " + "; ".join(f"{s} {f(r['power']['depth_spearman_of_f'])}" for s, r in est.items() if r.get("status") == "ok"), ""]
    L += ["## Part 3: B12 at equal depth (keratin share vs secretory score, Spearman + MSR, no regression)", "", "| specimen | n | stored primary rho (q) | stored partial rho (q) | equal-depth rho mean [min, max] | combined p | BH q | valid |", "|---|---|---|---|---|---|---|---|"]
    for s, r in est.items():
        if r.get("status") != "ok":
            continue
        pr, pp = b12["per_specimen"]["primary"].get(s, {}), b12["per_specimen"]["partial_primary"].get(s, {})
        v = r["tests"]["k3"]
        L.append(f"| {s} | {r['n']} | {f(pr.get('rho'), 4)} ({f(pr.get('q'), 4)}) | {f(pp.get('rho'), 4)} ({f(pp.get('q'), 4)}) | {agg_cell(v, 4)} | {f(v.get('p_combined'), 4)} | {f(v.get('q_bh'), 4)} | {v['n_valid']}{'' if v['status'] == 'ok' else ' (' + str(v.get('reason')) + ')'} |")
    L += ["", f"R6: {coh['R6']}", ""]
    h = per_spec.get(HEADLINE)
    if h and h.get("status") == "ok":
        L += [f"P4-S2 ({HEADLINE}) Part 3 exact: {h['tests']['k3']}", ""]
    L += ["## Detection and per-specimen n", "", "| specimen | n | MDK>=1 % thin [min,max] (orig) | keratin>=1 % thin [min,max] (orig) | keratin share>0 % thin [min,max] (orig) | secretory gene >=1 % (mean over 11) | total mean (sd) | MDK N_P>0 % (orig) | valid draws m1 / d_mdk / d_sec / k3 |", "|---|---|---|---|---|---|---|---|---|"]
    for s, r in est.items():
        if r.get("status") != "ok":
            continue
        t, o, rg = r["detection_thinned"], r["detection_original"], r["detection_thinned_range"]
        L.append(f"| {s} | {r['n']} | {f(t['pct_mdk_ge1'], 1)} [{f(rg['pct_mdk_ge1'][0], 1)},{f(rg['pct_mdk_ge1'][1], 1)}] ({f(o['pct_mdk_ge1'], 1)}) | {f(t['pct_keratin_ge1'], 1)} [{f(rg['pct_keratin_ge1'][0], 1)},{f(rg['pct_keratin_ge1'][1], 1)}] ({f(o['pct_keratin_ge1'], 1)}) | "
                 f"{f(t['pct_kshare_gt0'], 1)} [{f(rg['pct_kshare_gt0'][0], 1)},{f(rg['pct_kshare_gt0'][1], 1)}] ({f(o['pct_kshare_gt0'], 1)}) | {f(t['pct_secretory_gene_ge1_mean'], 1)} | {f(t['total_mean'], 1)} ({f(t['total_sd'], 1)}) | {f(o['pct_mdk_np_gt0'], 1)} | "
                 f"{r['tests']['m1']['n_valid']} / {r['tests']['d_mdk']['n_valid']} / {r['tests']['d_sec']['n_valid']} / {r['tests']['k3']['n_valid']} |")
    return "\n".join(L) + "\n"


def flat_rows(per_spec: dict) -> dict:
    p1, ct, p3, pw = [], [], [], []
    for s, r in per_spec.items():
        if r.get("status") != "ok":
            continue
        for tn in ("m1", "d_mdk", "d_sec"):
            v = r["tests"][tn]
            p1.append({"specimen": s, "n": r["n"], "test": tn, **{k: v.get(k, "") for k in ("status", "n_valid", "mean", "min", "max", "sd", "p_combined", "q_bh", "frac_p_lt_0.05")}})
        v = r["tests"]["k3"]
        p3.append({"specimen": s, "n": r["n"], "status": v["status"], **{k: v.get(k, "") for k in ("n_valid", "mean", "min", "max", "sd", "p_combined", "q_bh", "frac_p_lt_0.05")}})
        c = r["controls"]
        if c["status"] == "ok":
            for tag in ("a", "b"):
                ct.append({"specimen": s, "control": tag, "real_mean": c["real_mean"], "n_valid": c["n_valid"], **c[tag]})
        for sc in SCENARIOS:
            v = r["power"][sc]
            pw.append({"specimen": s, "scenario": sc, "status": v.get("status"), **{k: v.get(k, "") for k in ("target", "b", "I_noisefree", "n_valid", "mean", "min", "max", "p_combined", "q_bh", "frac_p_lt_0.05")}})
    return {"part1": p1, "controls": ct, "part3": p3, "power": pw}


def write_csv(path, rows) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["empty"])
        w.writeheader()
        w.writerows(rows)


def compute(comp_run, b9: dict, b17: dict, b12: dict, specs, out, meta_common: dict, target: float = TARGET, n_draws: int = N_DRAWS, b1: int = B1, b3: int = B3, n_ctrl: int = N_CTRL) -> tuple:
    """Check 10 (was ``main``) at the rule target: every specimen gated against the stored B9 (I to 1e-9, p exactly). ``b9`` / ``b17`` / ``b12`` = the stamped stage results.
    Writes check10_{part1,controls,part3,power}.csv and check10_tables.md into ``out``. Returns ``(summary, files)``; raises ``GateFailure``."""
    a = types.SimpleNamespace(target=target, n_draws=n_draws, b1=b1, b3=b3, n_ctrl=n_ctrl, comp_run=comp_run)
    gate_on = a.target == TARGET
    meta = {"check": "check10_depth_moran", "seed": SEED10, "target": target, "n_draws": n_draws, "min_valid_draws": MIN_VALID, "b1": b1, "b3": b3, "n_ctrl": n_ctrl, "k": R.K,
            "gate_on": gate_on, **meta_common}
    per_spec = {}
    for s in specs:
        per_spec[s] = run_specimen(s, a, b9, gate_on)
    gates = [r["gate"] for r in per_spec.values() if "gate" in r]
    meta["n_gates"] = len(gates)
    if gate_on:
        print(f"GATE PASSED ({len(gates)} specimens: B9 I within 1e-9, p equal)", flush=True)
        if not gates:
            raise CM.GateFailure("no gate ran: refusing")
    coh = cohort(per_spec, b9, b17, b12, gate_on)
    out = Path(out)
    files = []
    for name, rows in flat_rows(per_spec).items():
        write_csv(out / f"check10_{name}.csv", rows)
        files.append(f"check10_{name}.csv")
    (out / "check10_tables.md").write_text(md_tables(meta, per_spec, coh, b9, b17, b12))
    files.append("check10_tables.md")
    return {**meta, "gates": gates, "cohort": coh, "per_specimen": per_spec}, files
