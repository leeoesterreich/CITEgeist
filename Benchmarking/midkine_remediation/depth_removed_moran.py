"""B22 part 3 (check 11): is the MDK-secretory bivariate Moran's I (B9) and its matched controls (B17) depth? Full counts, depth removed by a linear projection (forms none, D1,
D2, D3), with matched controls and a planted-association power check. Descriptive only.

Promoted from the exploratory check 11 (``exploratory_v42/check11_depth_removed_moran.py``, git history). The analysis functions are that code, unchanged (``basis``,
``make_projector``, ``apply``, ``seed_for``, ``moran_pair``, ``pearson``, ``score``, ``control_I``, ``stored_ctrl_gate``, ``build``, ``raw_I``, ``bisect``, ``plant_realisation``,
``power_specimen``, ``cohort``, ``f``, ``md_tables``, ``write_csvs``; tests/midkine_remediation/test_moran_robustness.py checks the AST). Changed: ``run_specimen`` reads the B5 run
from ``a.comp_run`` (it was a module constant), the B17 file constant is gone (the stored B9 / B17 results are arguments), ``main`` became ``compute``, the table header no longer
prints a dirty flag. Rules: ``exploratory_v42/check11_DECISION_RULE.md`` (git history).
"""
from __future__ import annotations

import csv
import json
import sys
import time
import types
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import a12_10_real as R  # noqa: E402
import compartments as C  # noqa: E402
import control_pairs as CP  # noqa: E402
import depth_moran as C10  # noqa: E402
import hurdle as H  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

SEED11 = 20261111
N_REAL = 20
N_CTRL = 200
FORMS = ("none", "D1", "D2", "D3")
RNA_COVARIATE = "total"          # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md): D3 = depth + RNA content + tumour content
SCEN = ("S0", "S1", "S2", "S3", "A", "S1A")
B_MAX, B_STEPS = 8.0, 40
SEC = list(R.SECRETORY)
MDK = C.MDK
TOL = 1e-9


# --------------------------------------------------------------------------- projectors (tested in test_check11.py)
def basis(Z) -> np.ndarray:
    """Orthonormal basis of the column span of Z (rank-deficient directions dropped)."""
    U, s, _ = np.linalg.svd(np.asarray(Z, float), full_matrices=False)
    return U[:, s > 1e-10 * s[0]]


def make_projector(form: str, obs, depth) -> np.ndarray | None:
    """Orthonormal basis Q of the design span for ``form`` (None for 'none'); residual = v - Q Q^T v. May raise hurdle.CovariateDomainError for D3."""
    n = len(depth)
    if form == "none":
        return None
    if form == "D1":
        return basis(np.column_stack([np.ones(n), depth]))
    if form == "D2":
        g = CP.deciles(np.asarray(depth, float))
        return basis(np.column_stack([(g == k).astype(float) for k in np.unique(g)]))
    if form == "D3":
        return basis(np.column_stack([np.ones(n), H.covariate_matrix(obs, rna=RNA_COVARIATE)]))
    raise ValueError(form)


def apply(Q, V):
    """Residual of the columns of V (1-D or 2-D) after projection off span(Q)."""
    return V if Q is None else V - Q @ (Q.T @ V)


def seed_for(label: str) -> int:
    return int(np.random.SeedSequence([SEED11, zlib.crc32(label.encode())]).generate_state(1)[0])


def moran_pair(x, y, W, sn, seed: int, B: int) -> dict:
    """B9 recipe: I = morans_bv(x, y, W), MSR p two-sided with x fixed and y surrogated. Typed not_estimable for a constant variable."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        return {"status": "not_estimable", "reason": "non-finite"}
    if x.std() < 1e-12 or y.std() < 1e-12:
        return {"status": "not_estimable", "reason": "constant variable"}
    return {"status": "ok", "I": float(M.morans_bv(x, y, W)), "p": float(sn.pvalue(x, y, B=B, seed=seed)["p_two_sided"])}


def pearson(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float("nan") if a.std() < 1e-12 or b.std() < 1e-12 else float(np.corrcoef(a, b)[0, 1])


def score(cols) -> np.ndarray:
    """B9 secretory score of an (n, g) matrix."""
    return C10.mean_z_cols(np.asarray(cols, float))


# --------------------------------------------------------------------------- controls (B17 design, same removal)
def control_I(X, x, y, a_idx, b_idx, W, Q) -> tuple:
    """(I_a, I_b). x, y are the (already projected) real variables; control genes (a) and control scores (b) are projected with Q before I.
    I = z(x) . (W z(y)) / n as B17 (morans_bv). For Q = None this is exactly the B17 loop."""
    n = len(x)
    ca = apply(Q, X[:, np.asarray(a_idx)])
    sda = ca.std(0)
    oka = sda > 1e-12
    lag = W @ M.zscore(y)
    Ia = np.full(ca.shape[1], np.nan)
    Ia[oka] = ((ca[:, oka] - ca[:, oka].mean(0)) / sda[oka]).T @ lag / n
    nb, ns = b_idx.shape
    cb = X[:, b_idx.ravel()]
    sd = cb.std(0)
    zb = np.where(sd > 0, (cb - cb.mean(0)) / np.where(sd > 0, sd, 1.0), 0.0)
    Y = apply(Q, zb.reshape(n, nb, ns).mean(axis=2))
    sdy = Y.std(0)
    okb = sdy > 1e-12
    Ib = np.full(nb, np.nan)
    ZY = (Y[:, okb] - Y[:, okb].mean(0)) / sdy[okb]
    Ib[okb] = M.zscore(x) @ (W @ ZY) / n
    return Ia, Ib


def stored_ctrl_gate(spec: str, real: float, Ia, Ib, st: dict) -> dict:
    """Fail closed: the no-removal controls must reproduce the stored v2b171 block."""
    out = {}
    for tag, I, key in (("a", Ia, "control_a_random_mdk"), ("b", Ib, "control_b_random_secretory")):
        blk = C10.ctrl_block(real, I)
        s = st[key]
        for fld, got, exp in (("mean", blk["mean"], s["mean"]), ("p5", blk["p5"], s["p5"]), ("p95", blk["p95"], s["p95"]), ("pct", blk["pct"], s["real_percentile_in_controls"])):
            CM.gate_equal(f"{spec} control {tag} {fld} (no removal) vs stored B17", got, exp, rtol=TOL, atol=TOL)
        if blk["gt95"] != s["real_gt_p95"]:
            raise CM.GateFailure(f"GATE FAILED {spec} control {tag} real_gt_p95: got {blk['gt95']} expected {s['real_gt_p95']}")
        out[tag] = {"mean": blk["mean"], "p95": blk["p95"], "gt95": blk["gt95"]}
    return out


# --------------------------------------------------------------------------- planted association / artefact (power check)
def build(P, f, zd, b: float, c: float) -> np.ndarray:
    """Planted N_P columns: base x exp(b f) x exp(c z_depth); zeros stay zero."""
    return P * np.exp(b * f + c * zd)[:, None]


def raw_I(P, f, zd, b, c, W) -> float:
    out = build(P, f, zd, b, c)
    return float(M.morans_bv(out[:, 0], score(out[:, 1:]), W))


def bisect(fn, target: float) -> dict:
    """Smallest-effort bisection of param in [0, B_MAX] with fn(param) = noise-free raw I. fn(0) >= target -> param 0 (noted); fn(B_MAX) < target -> unreachable."""
    if target == 0:
        return {"param": 0.0, "I_noisefree": fn(0.0), "reachable": True}
    if fn(0.0) >= target:
        return {"param": 0.0, "I_noisefree": fn(0.0), "reachable": True, "note": "I(0) already >= target"}
    if fn(B_MAX) < target:
        return {"param": B_MAX, "I_noisefree": fn(B_MAX), "reachable": False}
    lo, hi = 0.0, B_MAX
    for _ in range(B_STEPS):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if fn(mid) < target else (lo, mid)
    p = 0.5 * (lo + hi)
    return {"param": float(p), "I_noisefree": fn(p), "reachable": True}


def plant_realisation(spec: str, N0, zd, Qf, W, k: int) -> tuple:
    """Permuted base columns (n, 12) and the depth-orthogonal latent field f for realisation k."""
    n = N0.shape[0]
    rng = np.random.default_rng(np.random.SeedSequence([SEED11, zlib.crc32((spec + "|plant").encode()), k]))
    P = N0.copy()
    for j in range(P.shape[1]):
        P[:, j] = P[rng.permutation(n), j]
    f = W @ (W @ (W @ rng.normal(size=n)))
    f = apply(Qf, f)
    return P, (f - f.mean()) / f.std()


def power_specimen(spec: str, X, genes, depth, W, sn, projs: dict, targets: dict, n_real: int, B: int) -> dict:
    cols = [genes.index(MDK)] + [genes.index(g) for g in SEC]
    N0 = np.asarray(X[:, cols], float)
    zd = (depth - depth.mean()) / depth.std()
    g = CP.deciles(depth)
    Qf = basis(np.column_stack([np.ones(len(depth)), zd] + [(g == k).astype(float) for k in np.unique(g)]))
    rec = {sc: {fm: [] for fm in FORMS} for sc in SCEN}
    calib = {sc: [] for sc in SCEN}
    sp_f = []
    for k in range(n_real):
        P, f = plant_realisation(spec, N0, zd, Qf, W, k)
        sp_f.append(float(abs(pearson(f, depth))))
        art = bisect(lambda c: raw_I(P, f, zd, 0.0, c, W), targets["S1"])
        pars = {"S0": (0.0, 0.0), "A": (0.0, art["param"])}
        cal = {"A": art}
        for sc in ("S1", "S2", "S3"):
            if targets.get(sc) is None:
                continue
            cal[sc] = bisect(lambda b: raw_I(P, f, zd, b, 0.0, W), targets[sc])
            pars[sc] = (cal[sc]["param"], 0.0)
        pars["S1A"] = (cal["S1"]["param"], art["param"])
        for sc in SCEN:
            if sc not in pars or (sc in cal and not cal[sc]["reachable"]) or (sc == "S1A" and not (cal["S1"]["reachable"] and art["reachable"])):
                for fm in FORMS:
                    rec[sc][fm].append(C10.rec({}, "target not reachable"))
                continue
            b, c = pars[sc]
            calib[sc].append({"param": c if sc == "A" else b, "I_noisefree": cal[sc]["I_noisefree"] if sc in cal else float("nan")})
            out = build(P, f, zd, b, c)
            x, y = out[:, 0], score(out[:, 1:])
            for fm in FORMS:
                Q = projs[fm]
                if fm != "none" and Q is None:
                    rec[sc][fm].append(C10.rec({}, "projector not estimable"))
                    continue
                rec[sc][fm].append(C10.rec(moran_pair(apply(Q, x), apply(Q, y), W, sn, seed_for(f"{spec}|{sc}|{fm.removesuffix('_percell')}|k{k}"), B)))   # a diagnostic per-cell form draws the SAME MSR seeds as the form it reproduces (D5_percell = D5)
    res = {"depth_abs_pearson_of_f_mean": float(np.mean(sp_f)), "depth_abs_pearson_of_f_max": float(np.max(sp_f)), "scenarios": {}}
    for sc in SCEN:
        if sc == "S3" and targets.get("S3") is None:
            res["scenarios"][sc] = {"status": "not_applicable", "reason": "stored I <= 0"}
            continue
        res["scenarios"][sc] = {"calibration_mean_param": float(np.mean([c_["param"] for c_ in calib[sc]])) if calib[sc] else None,
                                "calibration_mean_I_noisefree": float(np.mean([c_["I_noisefree"] for c_ in calib[sc]])) if calib[sc] else None,
                                "forms": {fm: C10.agg(rec[sc][fm]) for fm in FORMS}}
    return res


# --------------------------------------------------------------------------- one specimen
def run_specimen(spec: str, a, b9: dict, b17: dict) -> dict:
    t0 = time.time()
    d = R.load_cgeist_eligible(a.comp_run, spec)
    X, genes, obs, n = d["X"], list(d["genes"]), d["obs"], d["n"]
    print(f"  {spec}: loaded n={n} {time.time() - t0:.0f}s", flush=True)
    res = {"specimen": spec, "n": n}
    if n < R.MIN_EL or n > M.DENSE_EIGH_MAX_N:
        return {**res, "status": "not_estimable", "reason": f"n = {n}"}
    W = knn_weights(np.asarray(d["xy"], float), R.K)
    sn = M.SpectralNull(W)
    print(f"  {spec}: spectral null {time.time() - t0:.0f}s", flush=True)
    mdk_i = genes.index(MDK)
    sec_idx = [genes.index(g) for g in SEC]
    depth = np.log(obs["spot_umi"].to_numpy(float))
    x0 = X[:, mdk_i]
    y0 = score(X[:, sec_idx])
    stored = b9["per_specimen"][spec]
    stc = b17["per_specimen"][spec]
    pool_dec, pool = CP.prepare_pool([X], genes, X, mdk_i, sec_idx)
    dr = CP.make_draws(genes, pool_dec, pool, mdk_i, sec_idx, spec + "|cg_alone", n_ctrl=a.n_ctrl, n_uni=1)
    projs, form_err = {}, {}
    for fm in FORMS:
        try:
            projs[fm] = make_projector(fm, obs, depth)
        except H.CovariateDomainError as e:
            projs[fm], form_err[fm] = None, str(e)
    res["forms"] = {}
    for fm in FORMS:
        if fm in form_err:
            res["forms"][fm] = {"status": "not_estimable", "reason": form_err[fm]}
            continue
        Q = projs[fm]
        x, y = apply(Q, x0), apply(Q, y0)
        r = moran_pair(x, y, W, sn, R.SEED, R.B)
        out = {**r}
        if r["status"] == "ok":
            Ia, Ib = control_I(X, x, y, dr["a"], dr["b"], W, Q)
            out["controls"] = {"a": C10.ctrl_block(r["I"], Ia), "b": C10.ctrl_block(r["I"], Ib), "n_a_nan": int(np.isnan(Ia).sum()), "n_b_nan": int(np.isnan(Ib).sum())}
            if fm == "none":                                         # fidelity gates: stored B9 (I, p exactly) and stored B17 controls
                out["gate_b9"] = CM.gate_equal(f"{spec} B9 I (no removal) vs stored", r["I"], stored["I"], rtol=TOL, atol=TOL)
                if r["p"] != stored["p_two_sided"]:
                    raise CM.GateFailure(f"GATE FAILED {spec} B9 p_two_sided (no removal): got {r['p']!r} expected {stored['p_two_sided']!r}")
                out["gate_b17"] = stored_ctrl_gate(spec, r["I"], Ia, Ib, stc)
            out["depth"] = {}
            for nm, v in (("x", x), ("y", y)):
                out["depth"][nm] = {"pearson_r": pearson(depth, v), **moran_pair(depth, v, W, sn, R.SEED, R.B)}
        res["forms"][fm] = out
    print(f"  {spec}: forms and controls done {time.time() - t0:.0f}s", flush=True)
    targets = {"S1": b9["aggregate"]["mean"], "S2": 2 * b9["aggregate"]["mean"], "S3": stored["I"] if stored["I"] > 0 else None}
    res["power"] = power_specimen(spec, X, genes, depth, W, sn, projs, targets, a.n_real, a.b_power)
    res["status"] = "ok"
    print(f"{spec}: n={n} gates passed; I none/D1/D2/D3 = " + " ".join(f"{res['forms'][fm].get('I', float('nan')):.4f}" for fm in FORMS) + f" {time.time() - t0:.0f}s", flush=True)
    return res


# --------------------------------------------------------------------------- cohort and tables
def cohort(per: dict, b9: dict, b17: dict) -> dict:
    est = {s: r for s, r in per.items() if r.get("status") == "ok"}
    out = {"stored": {"b9_mean_I": b9["aggregate"]["mean"], "b9_sig": b9["bh"]["significant_q_lt_0.05"], "b17": b17["cohort"]["citegeist_alone_definition_C"]}}
    sig = b9["bh"]["significant_q_lt_0.05"]
    for fm in FORMS:
        ok = {s: r["forms"][fm] for s, r in est.items() if r["forms"][fm].get("status") == "ok"}
        q = R._bh([v["p"] for v in ok.values()]) if ok else []
        for (s, v), qi in zip(ok.items(), q):
            v["q_bh"] = float(qi)
        for nm in ("x", "y"):                                          # BH of the depth-check p across specimens, per column
            dd = {s: v["depth"][nm] for s, v in ok.items() if v["depth"][nm].get("status") == "ok"}
            qq = R._bh([v["p"] for v in dd.values()]) if dd else []
            for (s, v), qi in zip(dd.items(), qq):
                v["q_bh"] = float(qi)
    for sc in SCEN:
        for fm in FORMS:
            C10.bh_attach(est, lambda r, sc=sc, fm=fm: r["power"]["scenarios"][sc].get("forms", {}).get(fm) if r["power"]["scenarios"][sc].get("forms") else None)
    co = {}
    for fm in FORMS:
        ok = {s: est[s]["forms"][fm] for s in est if est[s]["forms"][fm].get("status") == "ok"}
        if not ok:
            co[fm] = {"status": "not_estimable"}
            continue
        r1 = {s: ("remains q<0.05" if (s in ok and ok[s]["q_bh"] < 0.05) else ("lost" if s in ok else "not estimable")) for s in sig}

        def pw(sc):
            r_ = {}
            for s in est:
                v = est[s]["power"]["scenarios"][sc].get("forms", {}).get(fm)
                r_[s] = v if v else {"status": "not_applicable"}
            return r_

        rec_ = {sc: pw(sc) for sc in SCEN}
        n_rec = {sc: int(sum(v.get("status") == "ok" and v["p_combined"] < 0.05 for v in rec_[sc].values())) for sc in SCEN}
        n_est = {sc: int(sum(v.get("status") == "ok" for v in rec_[sc].values())) for sc in SCEN}
        lost = {s: {sc: (bool(rec_[sc][s]["p_combined"] < 0.05) if rec_[sc][s].get("status") == "ok" else None) for sc in ("S1", "S2", "S3")} for s, v in r1.items() if v == "lost"}
        art = {s: {"raw_visible": bool(est[s]["power"]["scenarios"]["A"]["forms"]["none"].get("p_combined", 1) < 0.05),
                   "form_p_ge_0.05": bool(rec_["A"][s].get("p_combined", 0) >= 0.05) if rec_["A"][s].get("status") == "ok" else None} for s in est}
        co[fm] = {"n_estimable": len(ok), "R1": {"per_specimen": r1, "n_remaining": int(sum(v == "remains q<0.05" for v in r1.values())), "n_stored_sig": len(sig)},
                  "R2": {"mean_real_I": float(np.mean([v["I"] for v in ok.values()])), "mean_ctrl_a_mean": float(np.mean([v["controls"]["a"]["mean"] for v in ok.values()])),
                         "mean_ctrl_b_mean": float(np.mean([v["controls"]["b"]["mean"] for v in ok.values()])), "mean_ctrl_a_p95": float(np.mean([v["controls"]["a"]["p95"] for v in ok.values()])),
                         "mean_ctrl_b_p95": float(np.mean([v["controls"]["b"]["p95"] for v in ok.values()]))},
                  "R3": {"n_gt_p95_a": int(sum(v["controls"]["a"]["gt95"] for v in ok.values())), "n_gt_p95_b": int(sum(v["controls"]["b"]["gt95"] for v in ok.values()))},
                  "R4": {"lost": lost, "n_lost": len(lost), "n_informative_S1": int(sum(bool(v["S1"]) for v in lost.values())), "n_not_informative_S1": int(sum(not v["S1"] for v in lost.values()))},
                  "recovery_combined_p_lt_0.05": {sc: f"{n_rec[sc]} of {n_est[sc]}" for sc in SCEN},
                  "artefact": {"n_raw_visible": int(sum(v["raw_visible"] for v in art.values())), "n_removed_in_form_of_raw_visible": int(sum(v["raw_visible"] and bool(v["form_p_ge_0.05"]) for v in art.values())), "per_specimen": art},
                  "R5": {nm: {"n_q_lt_0.05": int(sum(v["depth"][nm].get("q_bh", 1) < 0.05 for v in ok.values())), "mean_I": float(np.mean([v["depth"][nm]["I"] for v in ok.values() if v["depth"][nm]["status"] == "ok"])),
                              "mean_pearson_r": float(np.mean([v["depth"][nm]["pearson_r"] for v in ok.values()]))} for nm in ("x", "y")}}
        co[fm]["R2"]["above_both_control_means"] = bool(co[fm]["R2"]["mean_real_I"] > co[fm]["R2"]["mean_ctrl_a_mean"] and co[fm]["R2"]["mean_real_I"] > co[fm]["R2"]["mean_ctrl_b_mean"])
    out["forms"] = co
    return out


def f(x, nd=4):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def md_tables(meta, per, coh, b9) -> str:
    L = [f"# Check 11 generated tables (n_real {meta['n_real']}, B power {meta['b_power']}, controls {meta['n_ctrl']})", "", f"Output dir `{meta['out']}`; commit {meta['git_commit'][:8]}; fidelity gates passed for {meta['n_gates']} specimens.", "",
         "## Part 1: Moran's I, MDK vs secretory, per form", "", "| specimen | n | stored I | stored q | none I (p) | D1 I (p, q) | D2 I (p, q) | D3 I (p, q) |", "|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            L.append(f"| {s} | {r.get('n')} | not estimable: {r.get('reason')} ||||||")
            continue
        st = b9["per_specimen"][s]
        cell = lambda fm: ("not est." if r["forms"][fm].get("status") != "ok" else f"{f(r['forms'][fm]['I'])} (p {f(r['forms'][fm]['p'], 3)}, q {f(r['forms'][fm].get('q_bh'))})")  # noqa: E731
        L.append(f"| {s} | {r['n']} | {f(st['I'])} | {f(b9['bh'][s]['q'])} | {cell('none')} | {cell('D1')} | {cell('D2')} | {cell('D3')} |")
    L += ["", "## Part 2: controls per form (real I; control a mean/p95/real pct; control b mean/p95/real pct; real>p95 a,b)", "", "| specimen | form | real I | a mean | a p95 | a pct | b mean | b p95 | b pct | >p95 a | >p95 b |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for fm in FORMS:
            v = r["forms"][fm]
            if v.get("status") != "ok":
                continue
            c = v["controls"]
            L.append(f"| {s} | {fm} | {f(v['I'])} | {f(c['a']['mean'])} | {f(c['a']['p95'])} | {f(c['a']['pct'], 1)} | {f(c['b']['mean'])} | {f(c['b']['p95'])} | {f(c['b']['pct'], 1)} | {c['a']['gt95']} | {c['b']['gt95']} |")
    L += ["", "## Part 1b: depth after removal (Pearson r and Moran's I of log_spot_umi with x, with y)", "", "| specimen | form | r(depth,x) | I(depth,x) | q | r(depth,y) | I(depth,y) | q |", "|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for fm in FORMS:
            v = r["forms"][fm]
            if v.get("status") != "ok":
                continue
            dx, dy = v["depth"]["x"], v["depth"]["y"]
            L.append(f"| {s} | {fm} | {f(dx['pearson_r'])} | {f(dx.get('I'))} | {f(dx.get('q_bh'))} | {f(dy['pearson_r'])} | {f(dy.get('I'))} | {f(dy.get('q_bh'))} |")
    L += ["", "## Power: per specimen, scenario and form: mean I [combined p] (fraction of realisations p<0.05)", "", "| specimen | scenario | none | D1 | D2 | D3 |", "|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for sc in SCEN:
            blk = r["power"]["scenarios"][sc]
            if "forms" not in blk:
                L.append(f"| {s} | {sc} | {blk.get('status')} ||||")
                continue
            cells = []
            for fm in FORMS:
                v = blk["forms"][fm]
                cells.append(f"not est. ({v.get('n_valid', 0)} valid)" if v["status"] != "ok" else f"{f(v['mean'])} [{f(v['p_combined'], 3)}] ({f(v['frac_p_lt_0.05'], 2)})")
            L.append(f"| {s} | {sc} | " + " | ".join(cells) + " |")
    L += ["", "## Cohort and reading rule", "", "```", json.dumps(coh, indent=1, default=float, sort_keys=True), "```", ""]
    return "\n".join(L)


def write_csvs(out, per):
    with open(out / "check11_part1.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["specimen", "n", "form", "status", "I", "p", "q_bh"])
        for s, r in per.items():
            for fm in FORMS:
                v = r.get("forms", {}).get(fm, {"status": r.get("status")})
                w.writerow([s, r.get("n"), fm, v.get("status"), v.get("I"), v.get("p"), v.get("q_bh")])
    with open(out / "check11_controls.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["specimen", "form", "real_I", "a_mean", "a_p5", "a_p95", "a_pct", "a_gt95", "b_mean", "b_p5", "b_p95", "b_pct", "b_gt95"])
        for s, r in per.items():
            for fm in FORMS:
                v = r.get("forms", {}).get(fm, {})
                if v.get("status") == "ok":
                    c = v["controls"]
                    w.writerow([s, fm, v["I"]] + [c[t][k] for t in ("a", "b") for k in ("mean", "p5", "p95", "pct", "gt95")])
    with open(out / "check11_depthcheck.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["specimen", "form", "var", "pearson_r", "I", "p", "q_bh"])
        for s, r in per.items():
            for fm in FORMS:
                v = r.get("forms", {}).get(fm, {})
                if v.get("status") == "ok":
                    for nm in ("x", "y"):
                        d = v["depth"][nm]
                        w.writerow([s, fm, nm, d["pearson_r"], d.get("I"), d.get("p"), d.get("q_bh")])
    with open(out / "check11_power.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["specimen", "scenario", "form", "status", "n_valid", "mean_I", "p_combined", "frac_p_lt_0.05", "q_bh", "calib_param_mean", "calib_I_noisefree_mean"])
        for s, r in per.items():
            if r.get("status") != "ok":
                continue
            for sc in SCEN:
                blk = r["power"]["scenarios"][sc]
                for fm in FORMS:
                    v = blk.get("forms", {}).get(fm)
                    if v:
                        w.writerow([s, sc, fm, v["status"], v.get("n_valid"), v.get("mean"), v.get("p_combined"), v.get("frac_p_lt_0.05"), v.get("q_bh"), blk["calibration_mean_param"], blk["calibration_mean_I_noisefree"]])


def compute(comp_run, b9: dict, b17: dict, specs, out, meta_common: dict, n_real: int = N_REAL, n_ctrl: int = N_CTRL, b_power: int = R.B) -> tuple:
    """Check 11 (was ``main``): every B9-estimable specimen, gated against the stored B9 (I to 1e-9, p exactly) and the stored B17 control blocks. Writes check11_{part1,controls,
    depthcheck,power}.csv and check11_tables.md into ``out``. Returns ``(summary, files)``; raises ``GateFailure``."""
    a = types.SimpleNamespace(n_real=n_real, n_ctrl=n_ctrl, b_power=b_power, comp_run=comp_run)
    per = {}
    for s in specs:
        per[s] = run_specimen(s, a, b9, b17)
    coh = cohort(per, b9, b17)
    meta = {"check": "check11_depth_removed_moran", "n_real": n_real, "n_ctrl": n_ctrl, "b_power": b_power, "seed": SEED11, "k": R.K,
            "n_gates": sum(r.get("status") == "ok" for r in per.values()), **meta_common}
    out = Path(out)
    write_csvs(out, per)
    (out / "check11_tables.md").write_text(md_tables(meta, per, coh, b9))
    return {**meta, "per_specimen": per, "cohort": coh}, ["check11_part1.csv", "check11_controls.csv", "check11_depthcheck.csv", "check11_power.csv", "check11_tables.md"]
