"""B22 part 5 (check 14): is the baseline Moran's I (about 0.03, also for random matched genes) a shared-denominator effect? Share forms none, S_A, S_B, S_C, S_D
(full counts, no thinning, no regression), matched controls, a planted-association power check and a shared-denominator artefact plant. Descriptive only.

Promoted from the exploratory check 14 (``exploratory_v42/check14_share_moran.py``, git history). The analysis functions are that code, unchanged (``seed_for``, ``share_pair``,
``form_matrix``, ``real_cols``, ``_I``, ``control_I_share_b``, ``_xy``, ``label14``, ``plant_realisation``, ``_planted_I``, ``power_form``, ``artefact_specimen``,
``reconstruction_gate``, ``denominator_block``, ``cohort``, ``f``, ``tables``, ``write_csvs``; tests/midkine_remediation/test_moran_robustness.py checks the AST). Changed: ``load``
takes the B5 run as an argument (it was a module constant), ``run_specimen`` passes ``a.comp_run``, ``main`` became ``compute``, the table header no longer prints a dirty flag.
The B17 file constant is gone (the stored B9 / B17 results are arguments). Rules: ``exploratory_v42/check14_DECISION_RULE.md`` (git history).
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
import depth_removed_moran as C11  # noqa: E402
import hurdle as H  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

SEED14 = 20261114
FORMS14 = ("none", "S_A", "S_B", "S_C", "S_D")
SCEN14 = ("S0", "S1", "S2")
S1_SIZE = 0.048217128908550164
CTRL_BASE = 0.034
EXCESS_MIN = 0.007
RECON_TOL = 1e-4
SEC, MDK = list(R.SECRETORY), C.MDK
DEN_UNI = ("I_log_nst", "I_log_dr", "I_log_den")
DEN_BI = ("nst~MDK", "nst~sec", "nst~lrna", "den~MDK", "den~sec", "den~lrna")


def seed_for(label: str) -> int:
    return int(np.random.SeedSequence([SEED14, zlib.crc32(label.encode())]).generate_state(1)[0])


# --------------------------------------------------------------------------- forms
def share_pair(E, tot, xi: int, yi) -> tuple:
    """S_B for ONE pair: x = value of gene xi, Y = (n, len(yi)) values of genes yi, all log1p(1e4 E_g / (tot - E_xi - sum E_yi)); 0 where that denominator <= 0."""
    yi = np.asarray(yi)
    den = tot - E[:, xi] - E[:, yi].sum(axis=1)
    ok = den > 0
    scale = np.where(ok, 1e4 / np.where(ok, den, 1.0), 0.0)
    return np.log1p(E[:, xi] * scale), np.log1p(E[:, yi] * scale[:, None])


def form_matrix(form: str, E, tot, ns, dr):
    """(n, G) values for the gene-independent forms. S_A is the stored N_S, S_C the stored N_A; S_B depends on the pair (``share_pair``)."""
    if form == "S_A":
        scale = np.where(tot > 0, 1e4 / np.where(tot > 0, tot, 1.0), 0.0)
        return np.log1p(E * scale[:, None])
    if form == "S_C":
        return np.log1p(E / dr[:, None])
    if form == "S_D":
        return np.log1p(E / ns[:, None])
    if form == "none":
        return np.log1p(E / ns[:, None] / dr[:, None])
    raise ValueError(form)


def real_cols(form: str, Xf, E, tot, mdk_i: int, sec_idx) -> np.ndarray:
    """(n, 12) columns [MDK, 11 secretory] of a form (S_B: the real pair's own denominator)."""
    if form == "S_B":
        x, Y = share_pair(E, tot, mdk_i, sec_idx)
        return np.column_stack([x, Y])
    return np.column_stack([Xf[:, mdk_i], Xf[:, sec_idx]])


def _I(x, y, W) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.std() < 1e-12 or y.std() < 1e-12:
        return float("nan")
    return float(M.morans_bv(x, y, W))


def control_I_share_b(E, tot, mdk_i, sec_idx, a_idx, b_idx, W) -> tuple:
    """S_B controls, pair by pair (each pair on its own denominator). a: control gene vs the 11 real secretory genes; b: real MDK vs 11 control genes."""
    Ia = np.array([_I(*_xy(share_pair(E, tot, int(c), sec_idx)), W) for c in a_idx])
    Ib = np.array([_I(*_xy(share_pair(E, tot, mdk_i, b_idx[j])), W) for j in range(b_idx.shape[0])])
    return Ia, Ib


def _xy(pair) -> tuple:
    x, Y = pair
    return x, C10.mean_z_cols(Y)


# --------------------------------------------------------------------------- reading rule
def label14(sa: dict, sb: dict) -> dict:
    """sa, sb = {'ctrl_a','ctrl_b','real','n_remaining'} cohort numbers (S_B needs real and n_remaining)."""
    den_thr = CTRL_BASE / 3
    baseline = bool(all(v <= den_thr for v in (sa["ctrl_a"], sa["ctrl_b"], sb["ctrl_a"], sb["ctrl_b"])))
    ex_a, ex_b = sb["real"] - sb["ctrl_a"], sb["real"] - sb["ctrl_b"]
    survives = bool(ex_a >= EXCESS_MIN and ex_b >= EXCESS_MIN and sb["n_remaining"] >= 3)
    lab = "+".join(n for n, f in (("DENOMINATOR_BASELINE", baseline), ("SIGNAL_SURVIVES", survives)) if f) or "NEITHER"
    return {"label": lab, "DENOMINATOR_BASELINE": baseline, "SIGNAL_SURVIVES": survives, "baseline_threshold": den_thr, "excess_min": EXCESS_MIN, "excess_S_B_vs_a": ex_a, "excess_S_B_vs_b": ex_b,
            "n_remaining_S_B": sb["n_remaining"]}


# --------------------------------------------------------------------------- power
def plant_realisation(spec: str, N0, W, k: int) -> tuple:
    n = N0.shape[0]
    rng = np.random.default_rng(np.random.SeedSequence([SEED14, zlib.crc32((spec + "|plant").encode()), k]))
    P = N0.copy()
    for j in range(P.shape[1]):
        P[:, j] = P[rng.permutation(n), j]
    f = W @ (W @ (W @ rng.normal(size=n)))
    return P, (f - f.mean()) / f.std()


def _planted_I(P, f, b, W) -> float:
    out = P * np.exp(b * f)[:, None]
    return float(M.morans_bv(out[:, 0], C10.mean_z_cols(out[:, 1:]), W))


def power_form(spec: str, N0, W, sn, n_real: int, B: int) -> dict:
    """S0/S1/S2 for one form's 12 columns (real sparsity), as check 11 without the depth projection."""
    targets = {"S1": S1_SIZE, "S2": 2 * S1_SIZE}
    rec = {sc: [] for sc in SCEN14}
    for k in range(n_real):
        P, f = plant_realisation(spec, N0, W, k)
        for sc in SCEN14:
            if sc == "S0":
                b, ok = 0.0, True
            else:
                cal = C11.bisect(lambda bb: _planted_I(P, f, bb, W), targets[sc])
                b, ok = cal["param"], cal["reachable"]
            if not ok:
                rec[sc].append(C10.rec({}, "target not reachable"))
                continue
            out = P * np.exp(b * f)[:, None]
            rec[sc].append(C10.rec(C11.moran_pair(out[:, 0], C10.mean_z_cols(out[:, 1:]), W, sn, seed_for(f"{spec}|{sc}|k{k}"), B)))
    return {sc: C10.agg(rec[sc]) for sc in SCEN14}


def artefact_specimen(spec: str, E12, T, den_med: float, W, sn, n_real: int, B: int) -> dict:
    """Independent genes (each column of E12 permuted) divided by a spatially structured common denominator D = den_med exp(c g): raw N_P-style vs S_A-style vs S_B-style."""
    n = E12.shape[0]
    rec = {k: [] for k in ("raw", "S_A", "S_B")}
    cs = []
    for k in range(n_real):
        rng = np.random.default_rng(np.random.SeedSequence([SEED14, zlib.crc32((spec + "|artefact").encode()), k]))
        Ep = E12.copy()
        for j in range(Ep.shape[1]):
            Ep[:, j] = Ep[rng.permutation(n), j]
        Tp = T[rng.permutation(n)]
        g = W @ (W @ (W @ rng.normal(size=n)))
        g = (g - g.mean()) / g.std()

        def raw(c, Ep=Ep, g=g):
            return np.log1p(Ep / (den_med * np.exp(c * g))[:, None])

        cal = C11.bisect(lambda c: float(M.morans_bv(raw(c)[:, 0], C10.mean_z_cols(raw(c)[:, 1:]), W)), S1_SIZE)
        if not cal["reachable"]:
            for kk in rec:
                rec[kk].append(C10.rec({}, "target not reachable"))
            continue
        cs.append(cal["param"])
        v_raw = raw(cal["param"])
        v_a = np.log1p(1e4 * Ep / Tp[:, None])
        den = Tp - Ep.sum(axis=1)
        okd = den > 0
        v_b = np.log1p(Ep * np.where(okd, 1e4 / np.where(okd, den, 1.0), 0.0)[:, None])
        for kk, v in (("raw", v_raw), ("S_A", v_a), ("S_B", v_b)):
            rec[kk].append(C10.rec(C11.moran_pair(v[:, 0], C10.mean_z_cols(v[:, 1:]), W, sn, seed_for(f"{spec}|art|{kk}|k{k}"), B)))
    return {"c_mean": float(np.mean(cs)) if cs else None, **{kk: C10.agg(rec[kk]) for kk in rec}}


# --------------------------------------------------------------------------- data
def load(spec: str, comp_run) -> dict:
    """Eligible cancer rows: stored N_P (X), E, N_S, N_A, obs, coordinates (a12_10_real.load_cgeist_eligible plus layers E, N_S, N_A)."""
    import anndata as ad

    a = ad.read_h5ad(str(Path(comp_run) / f"{spec}.h5ad"), backed="r")
    el = C.eligible_cancer(a.obs).to_numpy().astype(bool)
    out = {"obs": a.obs[el], "genes": list(a.var_names), "X": R._dense(a.X[el]), "E": R._dense(a.layers["E"][el]), "N_S": R._dense(a.layers["N_S"][el]), "N_A": R._dense(a.layers["N_A"][el]),
           "xy": np.asarray(a.obsm["spatial"], float)[el], "n": int(el.sum())}
    a.file.close()
    return out


def reconstruction_gate(spec: str, d: dict) -> dict:
    """The denominators are understood: N_P, N_S and N_A recomputed from E and obs equal the stored layers (float32 storage)."""
    E, o = d["E"], d["obs"]
    tot = E.sum(axis=1)
    ns, dr = o["n_st_eff"].to_numpy(float), o["depth_ratio"].to_numpy(float)
    res = {}
    for nm, got, stored in (("N_P", form_matrix("none", E, tot, ns, dr), d["X"]), ("N_S/S_A", form_matrix("S_A", E, tot, ns, dr), d["N_S"]), ("N_A/S_C", form_matrix("S_C", E, tot, ns, dr), d["N_A"])):
        m = float(np.abs(got - stored).max())
        if not m <= RECON_TOL:
            raise CM.GateFailure(f"GATE FAILED {spec} reconstruction {nm}: max |recomputed - stored| = {m!r} > {RECON_TOL}")
        res[nm] = m
    return res


def denominator_block(obs, x, y, W, sn) -> dict:
    ns, dr = obs["n_st_eff"].to_numpy(float), obs["depth_ratio"].to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lrna = np.log(H.rna_total(obs))      # protocol deviation 2026-11-06: "lrna" = log RNA content (rna_total), the covariate of the total mode; key names unchanged
    v = {"nst": np.log(ns), "dr": np.log(dr), "den": np.log(ns) + np.log(dr)}
    out = {}
    for key, nm in (("nst", "I_log_nst"), ("dr", "I_log_dr"), ("den", "I_log_den")):
        out[nm] = {**C11.moran_pair(v[key], v[key], W, sn, R.SEED, R.B)}
    for nm, (a_, b_) in {"nst~MDK": ("nst", x), "nst~sec": ("nst", y), "nst~lrna": ("nst", lrna), "den~MDK": ("den", x), "den~sec": ("den", y), "den~lrna": ("den", lrna)}.items():
        va = v[a_]
        out[nm] = {**C11.moran_pair(va, b_, W, sn, R.SEED, R.B), "pearson_r": C11.pearson(va, b_)}
    return out


def run_specimen(spec: str, a, b9: dict, b17: dict) -> dict:
    t0 = time.time()
    d = load(spec, a.comp_run)
    n, genes, obs = d["n"], list(d["genes"]), d["obs"]
    print(f"  {spec}: loaded n={n} {time.time() - t0:.0f}s", flush=True)
    res = {"specimen": spec, "n": n}
    if n < R.MIN_EL or n > M.DENSE_EIGH_MAX_N:
        return {**res, "status": "not_estimable", "reason": f"n = {n}"}
    res["reconstruction_max_abs_diff"] = reconstruction_gate(spec, d)
    X, E = d["X"], np.asfortranarray(d["E"])
    tot = E.sum(axis=1)
    ns, dr = obs["n_st_eff"].to_numpy(float), obs["depth_ratio"].to_numpy(float)
    W = knn_weights(np.asarray(d["xy"], float), R.K)
    sn = M.SpectralNull(W)
    mdk_i = genes.index(MDK)
    sec_idx = [genes.index(g) for g in SEC]
    stored, stc = b9["per_specimen"][spec], b17["per_specimen"][spec]
    pool_dec, pool = CP.prepare_pool([X], genes, X, mdk_i, sec_idx)
    dr_ = CP.make_draws(genes, pool_dec, pool, mdk_i, sec_idx, spec + "|cg_alone", n_ctrl=a.n_ctrl, n_uni=1)
    res["forms"] = {}
    cols = {}
    for fm in FORMS14:
        Xf = X if fm == "none" else (None if fm == "S_B" else form_matrix(fm, E, tot, ns, dr))
        cols[fm] = real_cols(fm, Xf, E, tot, mdk_i, sec_idx)
        x, y = cols[fm][:, 0], C10.mean_z_cols(cols[fm][:, 1:])
        r = C11.moran_pair(x, y, W, sn, R.SEED, R.B)
        out = {**r}
        if r["status"] == "ok":
            Ia, Ib = control_I_share_b(E, tot, mdk_i, sec_idx, dr_["a"], dr_["b"], W) if fm == "S_B" else C11.control_I(Xf, x, y, dr_["a"], dr_["b"], W, None)
            out["controls"] = {"a": C10.ctrl_block(r["I"], Ia), "b": C10.ctrl_block(r["I"], Ib), "n_a_nan": int(np.isnan(Ia).sum()), "n_b_nan": int(np.isnan(Ib).sum())}
            if fm == "none":                                   # fidelity gates: stored B9 (I to 1e-9, p exactly) and stored B17 controls (needs the 200 stored draws)
                out["gate_b9"] = CM.gate_equal(f"{spec} B9 I (no transform) vs stored", r["I"], stored["I"], rtol=C11.TOL, atol=C11.TOL)
                if r["p"] != stored["p_two_sided"]:
                    raise CM.GateFailure(f"GATE FAILED {spec} B9 p_two_sided (no transform): got {r['p']!r} expected {stored['p_two_sided']!r}")
                if a.n_ctrl == C11.N_CTRL:
                    out["gate_b17"] = C11.stored_ctrl_gate(spec, r["I"], Ia, Ib, stc)
        res["forms"][fm] = out
        del Xf
    print(f"  {spec}: forms and controls done {time.time() - t0:.0f}s", flush=True)
    res["denominator"] = denominator_block(obs, cols["none"][:, 0], C10.mean_z_cols(cols["none"][:, 1:]), W, sn)
    res["power"] = {fm: power_form(spec, cols[fm], W, sn, a.n_real, a.b_power) for fm in FORMS14}
    print(f"  {spec}: power done {time.time() - t0:.0f}s", flush=True)
    res["artefact"] = artefact_specimen(spec, np.asarray(E[:, [mdk_i] + sec_idx]), tot, float(np.median(ns * dr)), W, sn, a.n_real, a.b_power)
    res["status"] = "ok"
    print(f"{spec}: n={n} gates passed; I " + " ".join(f"{fm} {res['forms'][fm].get('I', float('nan')):.4f}" for fm in FORMS14) + f" {time.time() - t0:.0f}s", flush=True)
    return res


# --------------------------------------------------------------------------- cohort and tables
def cohort(per: dict, b9: dict, b17: dict) -> dict:
    est = {s: r for s, r in per.items() if r.get("status") == "ok"}
    sig = b9["bh"]["significant_q_lt_0.05"]
    out = {"stored": {"b9_mean_I": b9["aggregate"]["mean"], "b9_sig": sig, "b17": b17["cohort"]["citegeist_alone_definition_C"]}}
    for fm in FORMS14:
        ok = {s: r["forms"][fm] for s, r in est.items() if r["forms"][fm].get("status") == "ok"}
        q = R._bh([v["p"] for v in ok.values()]) if ok else []
        for (s, v), qi in zip(ok.items(), q):
            v["q_bh"] = float(qi)
        for sc in SCEN14:
            C10.bh_attach(est, lambda r, sc=sc, fm=fm: r["power"][fm][sc])
    for nm in DEN_UNI + DEN_BI:                                        # BH across specimens per denominator measure
        dd = {s: r["denominator"][nm] for s, r in est.items() if r["denominator"][nm].get("status") == "ok"}
        qq = R._bh([v["p"] for v in dd.values()]) if dd else []
        for (s, v), qi in zip(dd.items(), qq):
            v["q_bh"] = float(qi)
    co = {}
    for fm in FORMS14:
        ok = {s: est[s]["forms"][fm] for s in est if est[s]["forms"][fm].get("status") == "ok"}
        if not ok:
            co[fm] = {"status": "not_estimable"}
            continue
        r1 = {s: ("remains q<0.05" if (s in ok and ok[s]["q_bh"] < 0.05) else ("lost" if s in ok else "not estimable")) for s in sig}
        rec_ = {sc: {s: est[s]["power"][fm][sc] for s in est} for sc in SCEN14}
        n_rec = {sc: int(sum(v.get("status") == "ok" and v["p_combined"] < 0.05 for v in rec_[sc].values())) for sc in SCEN14}
        n_est = {sc: int(sum(v.get("status") == "ok" for v in rec_[sc].values())) for sc in SCEN14}
        lost = {s: {sc: (bool(rec_[sc][s]["p_combined"] < 0.05) if rec_[sc][s].get("status") == "ok" else None) for sc in ("S1", "S2")} for s, v in r1.items() if v == "lost"}
        mr = float(np.mean([v["I"] for v in ok.values()]))
        ma, mb = float(np.mean([v["controls"]["a"]["mean"] for v in ok.values()])), float(np.mean([v["controls"]["b"]["mean"] for v in ok.values()]))
        co[fm] = {"n_estimable": len(ok), "R1": {"per_specimen": r1, "n_remaining": int(sum(v == "remains q<0.05" for v in r1.values())), "n_stored_sig": len(sig)},
                  "R2": {"mean_real_I": mr, "mean_ctrl_a_mean": ma, "mean_ctrl_b_mean": mb, "excess_vs_a": mr - ma, "excess_vs_b": mr - mb,
                         "mean_ctrl_a_p95": float(np.mean([v["controls"]["a"]["p95"] for v in ok.values()])), "mean_ctrl_b_p95": float(np.mean([v["controls"]["b"]["p95"] for v in ok.values()]))},
                  "R3": {"n_gt_p95_a": int(sum(v["controls"]["a"]["gt95"] for v in ok.values())), "n_gt_p95_b": int(sum(v["controls"]["b"]["gt95"] for v in ok.values()))},
                  "R4": {"lost": lost, "n_lost": len(lost), "n_informative_S1": int(sum(bool(v["S1"]) for v in lost.values())), "n_not_informative_S1": int(sum(not v["S1"] for v in lost.values()))},
                  "recovery_combined_p_lt_0.05": {sc: f"{n_rec[sc]} of {n_est[sc]}" for sc in SCEN14}}
    out["forms"] = co
    out["reading"] = label14({"ctrl_a": co["S_A"]["R2"]["mean_ctrl_a_mean"], "ctrl_b": co["S_A"]["R2"]["mean_ctrl_b_mean"]},
                             {"ctrl_a": co["S_B"]["R2"]["mean_ctrl_a_mean"], "ctrl_b": co["S_B"]["R2"]["mean_ctrl_b_mean"], "real": co["S_B"]["R2"]["mean_real_I"], "n_remaining": co["S_B"]["R1"]["n_remaining"]})
    den = {}
    for nm in DEN_UNI + DEN_BI:
        v = [r["denominator"][nm] for r in est.values() if r["denominator"][nm].get("status") == "ok"]
        den[nm] = {"n": len(v), "mean_I": float(np.mean([x["I"] for x in v])) if v else None, "n_q_lt_0.05": int(sum(x.get("q_bh", 1) < 0.05 for x in v)),
                   "mean_pearson_r": float(np.mean([x["pearson_r"] for x in v])) if v and "pearson_r" in v[0] else None}
    out["denominator"] = den
    art = {}
    for kk in ("raw", "S_A", "S_B"):
        C10.bh_attach(est, lambda r, kk=kk: r["artefact"][kk])
        v = [r["artefact"][kk] for r in est.values() if r["artefact"][kk]["status"] == "ok"]
        art[kk] = {"n": len(v), "mean_I": float(np.mean([x["mean"] for x in v])) if v else None, "n_combined_p_lt_0.05": int(sum(x["p_combined"] < 0.05 for x in v))}
    out["artefact"] = art
    return out


def f(x, nd=4):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def tables(meta, per, coh, b9) -> str:
    L = [f"# Check 14 generated tables (n_real {meta['n_real']}, B power {meta['b_power']}, controls {meta['n_ctrl']})", "", f"Output dir `{meta['out']}`; commit {meta['git_commit'][:8]}; gates passed for {meta['n_gates']} specimens.", "",
         "## Reading rule", "", "```", json.dumps(coh["reading"], indent=1, default=float), "```", "",
         "## Part 1: I (MSR p, BH q within form)", "", "| specimen | n | stored I (q) | " + " | ".join(FORMS14) + " |", "|---|---|---|" + "---|" * len(FORMS14)]
    for s, r in per.items():
        if r.get("status") != "ok":
            L.append(f"| {s} | {r.get('n')} | not estimable: {r.get('reason')} |" + "|" * len(FORMS14))
            continue
        cell = lambda v: "not est." if v.get("status") != "ok" else f"{f(v['I'])} (p {f(v['p'], 3)}, q {f(v.get('q_bh'))})"  # noqa: E731
        L.append(f"| {s} | {r['n']} | {f(b9['per_specimen'][s]['I'])} ({f(b9['bh'][s]['q'])}) | " + " | ".join(cell(r["forms"][fm]) for fm in FORMS14) + " |")
    L += ["", "## Part 2: controls per form", "", "| specimen | form | real I | a mean | a p95 | a pct | b mean | b p95 | b pct | >p95 a | >p95 b |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for fm in FORMS14:
            v = r["forms"][fm]
            if v.get("status") == "ok":
                c = v["controls"]
                L.append(f"| {s} | {fm} | {f(v['I'])} | {f(c['a']['mean'])} | {f(c['a']['p95'])} | {f(c['a']['pct'], 1)} | {f(c['b']['mean'])} | {f(c['b']['p95'])} | {f(c['b']['pct'], 1)} | {c['a']['gt95']} | {c['b']['gt95']} |")
    L += ["", "## Denominator: Moran's I (MSR p, BH q across specimens) and Pearson r, stored N_P", "", "| specimen | " + " | ".join(DEN_UNI + DEN_BI) + " |", "|---|" + "---|" * len(DEN_UNI + DEN_BI)]
    for s, r in per.items():
        if r.get("status") == "ok":
            L.append(f"| {s} | " + " | ".join("not est." if r["denominator"][nm].get("status") != "ok" else f"{f(r['denominator'][nm]['I'])} (q {f(r['denominator'][nm].get('q_bh'), 3)}" + (f", r {f(r['denominator'][nm]['pearson_r'], 3)})" if "pearson_r" in r["denominator"][nm] else ")") for nm in DEN_UNI + DEN_BI) + " |")
    L += ["", "## Power: mean I [combined p] (frac p<0.05)", "", "| specimen | scenario | " + " | ".join(FORMS14) + " |", "|---|---|" + "---|" * len(FORMS14)]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for sc in SCEN14:
            cells = []
            for fm in FORMS14:
                v = r["power"][fm][sc]
                cells.append(f"not est. ({v.get('n_valid', 0)} valid)" if v["status"] != "ok" else f"{f(v['mean'])} [{f(v['p_combined'], 3)}] ({f(v['frac_p_lt_0.05'], 2)})")
            L.append(f"| {s} | {sc} | " + " | ".join(cells) + " |")
    L += ["", "## Shared-denominator artefact plant: mean I [combined p] (frac p<0.05); c = mean calibrated structure strength", "", "| specimen | c | raw N_P-style | S_A-style | S_B-style |", "|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") == "ok":
            ar = r["artefact"]
            L.append(f"| {s} | {f(ar['c_mean'], 3)} | " + " | ".join("not est." if ar[k]["status"] != "ok" else f"{f(ar[k]['mean'])} [{f(ar[k]['p_combined'], 3)}] ({f(ar[k]['frac_p_lt_0.05'], 2)})" for k in ("raw", "S_A", "S_B")) + " |")
    L += ["", "## Cohort", "", "```", json.dumps(coh, indent=1, default=float, sort_keys=True), "```", ""]
    return "\n".join(L)


def write_csvs(out, per):
    def wr(name, header, rows):
        with open(out / name, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)

    est = {s: r for s, r in per.items() if r.get("status") == "ok"}
    wr("check14_part1.csv", ["specimen", "n", "form", "status", "I", "p", "q_bh"], [[s, r.get("n"), fm, r["forms"][fm].get("status"), r["forms"][fm].get("I"), r["forms"][fm].get("p"), r["forms"][fm].get("q_bh")] for s, r in est.items() for fm in FORMS14])
    wr("check14_controls.csv", ["specimen", "form", "real_I"] + [f"{t}_{k}" for t in "ab" for k in ("mean", "p5", "p95", "pct", "gt95")],
       [[s, fm, v["I"]] + [v["controls"][t][k] for t in "ab" for k in ("mean", "p5", "p95", "pct", "gt95")] for s, r in est.items() for fm in FORMS14 for v in [r["forms"][fm]] if v.get("status") == "ok"])
    wr("check14_power.csv", ["specimen", "form", "scenario", "status", "n_valid", "mean_I", "p_combined", "frac_p_lt_0.05", "q_bh"],
       [[s, fm, sc, v["status"], v.get("n_valid"), v.get("mean"), v.get("p_combined"), v.get("frac_p_lt_0.05"), v.get("q_bh")] for s, r in est.items() for fm in FORMS14 for sc in SCEN14 for v in [r["power"][fm][sc]]])
    wr("check14_denominator.csv", ["specimen", "measure", "status", "I", "p", "q_bh", "pearson_r"], [[s, nm, v.get("status"), v.get("I"), v.get("p"), v.get("q_bh"), v.get("pearson_r")] for s, r in est.items() for nm, v in r["denominator"].items()])
    wr("check14_artefact.csv", ["specimen", "style", "c_mean", "status", "n_valid", "mean_I", "p_combined", "frac_p_lt_0.05", "q_bh"],
       [[s, k, r["artefact"]["c_mean"], v["status"], v.get("n_valid"), v.get("mean"), v.get("p_combined"), v.get("frac_p_lt_0.05"), v.get("q_bh")] for s, r in est.items() for k in ("raw", "S_A", "S_B") for v in [r["artefact"][k]]])


def compute(comp_run, b9: dict, b17: dict, specs, out, meta_common: dict, n_real: int = C11.N_REAL, n_ctrl: int = C11.N_CTRL, b_power: int = R.B) -> tuple:
    """Check 14 (was ``main``): gated against the stored B9 (I to 1e-9, p exactly), the stored B17 control blocks and the reconstruction of N_P / N_S / N_A. Writes
    check14_{part1,controls,power,denominator,artefact}.csv and check14_tables.md into ``out``. Returns ``(summary, files)``; raises ``GateFailure``."""
    a = types.SimpleNamespace(n_real=n_real, n_ctrl=n_ctrl, b_power=b_power, comp_run=comp_run)
    per = {}
    for s in specs:
        per[s] = run_specimen(s, a, b9, b17)
    coh = cohort(per, b9, b17)
    meta = {"check": "check14_share_moran", "n_real": n_real, "n_ctrl": n_ctrl, "b_power": b_power, "seed": SEED14, "k": R.K,
            "n_gates": sum(r.get("status") == "ok" for r in per.values()), **meta_common}
    out = Path(out)
    write_csvs(out, per)
    (out / "check14_tables.md").write_text(tables(meta, per, coh, b9))
    return {**meta, "per_specimen": per, "cohort": coh}, [f"check14_{nm}.csv" for nm in ("part1", "controls", "power", "denominator", "artefact")] + ["check14_tables.md"]
