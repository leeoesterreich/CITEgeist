"""B22 part 7 (check 17): is check 13's D5 loss of the MDK-secretory Moran's I an artefact of the covariate (self-inclusion, shared denominator)? Forms none, D5 (check 13
exactly), D5x (RNA-per-cell covariate recomputed per pair without the pair's genes) and D5xc (covariate = raw allocated-count total without the pair's genes, no division by the
shared n_st_eff x depth_ratio). Descriptive only.

Promoted from the exploratory check 17 (``exploratory_v42/check17_d5_fixed_moran.py``, git history). The analysis functions are that code, unchanged (``extra_genes``,
``pair_cov``, ``design_basis``, ``pair_basis``, ``control_I_pairwise``, ``label17``, ``gate_vs_check13``, ``covariate_gate``, ``CovCtx``; tests/midkine_remediation/test_moran_robustness.py
checks the AST). Changed: ``run_specimen`` reads the B5 run from ``a.comp_run`` (it was a module constant), the check 13 summary file the exploratory check gated against is the
check 13 result computed earlier IN THE SAME RUN (``c13``, an argument), ``main`` became ``compute`` (check 11's module global ``FORMS`` is patched for the duration of the call and
restored after, so nothing computed later sees check 11 changed), the table header no longer prints a dirty flag.
Rules: ``exploratory_v42/check17_DECISION_RULE.md`` (git history).
Protocol deviation 2026-11-06 (RNA-total covariate, ``PROTOCOL_DEVIATION_2026-11-06_rna_total.md``): check 13's D5 is now depth + RNA content (``hurdle.covariate_matrix(rna="total")``, ``rna_total`` =
rna_per_cell_equiv x n_st_eff x depth_ratio, the raw allocated-count total over genes not in keratins + MDK: no per-pair exclusion, no division by the shared denominator), so the form D5r of the
first version of this module is gone (D5 is that form). The old per-cell D5 is the check 13 arm ``D5_percell`` (``split_d3``). D5x and D5xc stay as cross-checks; ``run_specimen`` records
``d5xc_minus_d5_I`` (D5xc differs from D5 only by the pair's own genes). ``rna_total`` lives in ``hurdle`` and is asserted equal to the row sum of layer E over the non-keratin,
non-MDK genes (relative 1e-6) for every specimen (``rna_total_gate``).
"""
from __future__ import annotations

import contextlib
import json
import os
import sys
import tempfile
import time
import types
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

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
import share_moran as C14  # noqa: E402
import split_d3 as C13  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

FORMS17 = ("none", "D5", "D5x", "D5xc")
RNA_TOTAL_REL_TOL = 1e-6
EXCESS_MIN = 0.007
SEC, MDK = list(R.SECRETORY), C.MDK
REL_TOL_COV = 1e-4


# --------------------------------------------------------------------------- per-pair covariate
class CovCtx:
    """Everything needed to rebuild the RNA-per-cell covariate for any pair: the per-gene values E/D and E, and the always-excluded genes (keratins, MDK)."""

    def __init__(self, E, ns, dr, genes):
        genes = list(genes)
        self.E = np.asarray(E, float)
        self.pc = self.E / np.asarray(ns, float)[:, None] / np.asarray(dr, float)[:, None]      # E / n_st_eff / depth_ratio, as compartments.build_compartments
        self.base_excl = np.zeros(len(genes), bool)
        self.base_excl[[genes.index(g) for g in C.KERATINS + (C.MDK,)]] = True
        self.other = ~self.base_excl
        self.base_pc = self.pc[:, self.other].sum(axis=1)       # == obs rna_per_cell_equiv
        self.base_raw = self.E[:, self.other].sum(axis=1)       # raw allocated-count total over the same genes


def extra_genes(pair_idx, base_excl) -> np.ndarray:
    """The pair's genes that are not already excluded from the covariate (keratins, MDK)."""
    return np.array(sorted({int(g) for g in pair_idx if not base_excl[int(g)]}), dtype=int)


def pair_cov(ctx: CovCtx, kind: str, pair_idx) -> np.ndarray:
    """kind 'x': sum over S of E/D; kind 'xc': sum over S of E; S = genes not in keratins, MDK, or the pair. Rows where the subtraction leaves a rounding residue are summed directly."""
    ex = extra_genes(pair_idx, ctx.base_excl)
    M_, base = (ctx.pc, ctx.base_pc) if kind == "x" else (ctx.E, ctx.base_raw)
    out = base - M_[:, ex].sum(axis=1) if len(ex) else base.copy()
    bad = out <= 1e-9 * base
    if bad.any():
        keep = ctx.other.copy()
        keep[ex] = False
        out[bad] = M_[bad][:, keep].sum(axis=1)
    return out


def design_basis(depth, cov) -> np.ndarray | None:
    """Orthonormal basis of span[1, log spot_umi, log cov] (columns standardised as hurdle.covariate_matrix; the span is unchanged). None if cov is not finite and > 0."""
    cov = np.asarray(cov, float)
    if not (np.isfinite(cov).all() and (cov > 0).all()):
        return None
    return C11.basis(np.column_stack([np.ones(len(cov)), H.standardise(np.column_stack([depth, np.log(cov)]))]))


def pair_basis(ctx: CovCtx, depth, kind: str, pair_idx):
    return design_basis(depth, pair_cov(ctx, kind, pair_idx))


# --------------------------------------------------------------------------- controls
def _I_zz(x, y, W) -> float:
    sx, sy = x.std(), y.std()
    if not (sx > 1e-12 and sy > 1e-12):
        return float("nan")
    return float(M.zscore(x) @ (W @ M.zscore(y)) / len(x))


def control_I_pairwise(X, x0, y0, a_idx, b_idx, W, ctx: CovCtx, depth, kind: str, sec_idx, mdk_i: int) -> tuple:
    """(I_a, I_b) with the design rebuilt for every control pair. x0, y0 = UNPROJECTED real MDK N_P and real secretory score.
    a: control gene vs the real score, pair = {c} + 11 secretory. b: real MDK vs the control score, pair = the 11 control genes. Nan where the pair's covariate is not usable."""
    n = len(x0)
    Ia = np.full(len(a_idx), np.nan)
    for j, c in enumerate(np.asarray(a_idx)):
        Q = pair_basis(ctx, depth, kind, list(sec_idx) + [int(c)])
        if Q is not None:
            Ia[j] = _I_zz(C11.apply(Q, X[:, int(c)]), C11.apply(Q, y0), W)
    nb, ns_ = b_idx.shape
    cb = X[:, b_idx.ravel()]
    sd = cb.std(0)
    zb = np.where(sd > 0, (cb - cb.mean(0)) / np.where(sd > 0, sd, 1.0), 0.0)
    Y = zb.reshape(n, nb, ns_).mean(axis=2)
    Ib = np.full(nb, np.nan)
    for j in range(nb):
        Q = pair_basis(ctx, depth, kind, list(b_idx[j]))
        if Q is not None:
            Ib[j] = _I_zz(C11.apply(Q, x0), C11.apply(Q, Y[:, j]), W)
    return Ia, Ib


# --------------------------------------------------------------------------- reading rule
def label17(d5x: dict, d5xc: dict) -> dict:
    """d5x, d5xc = {'real','ctrl_a','ctrl_b','n_remaining'} cohort numbers."""
    out, flags = {}, {}
    for nm, d in (("SELF_INCLUSION", d5x), ("COUPLING", d5xc)):
        ea, eb = d["real"] - d["ctrl_a"], d["real"] - d["ctrl_b"]
        flags[nm] = bool(ea >= EXCESS_MIN and eb >= EXCESS_MIN and d["n_remaining"] >= 3)
        out[nm] = {"excess_vs_a": ea, "excess_vs_b": eb, "n_remaining": d["n_remaining"], "mean_real_I": d["real"]}
    return {"label": "+".join(k for k, v in flags.items() if v) or "NEITHER", **flags, "detail": out, "excess_min": EXCESS_MIN}


# --------------------------------------------------------------------------- gates
def gate_vs_check13(spec: str, res: dict, c13: dict) -> None:
    """D5 must reproduce check 13's stored status, I, p and control a/b blocks within 1e-9."""
    got, exp = res["forms"]["D5"], c13["per_specimen"][spec]["forms"]["D5"]
    if got.get("status") != exp.get("status"):
        raise CM.GateFailure(f"GATE FAILED {spec} D5 status {got.get('status')} vs check13 {exp.get('status')}")
    if exp.get("status") != "ok":
        return
    for k in ("I", "p"):
        CM.gate_equal(f"{spec} D5 {k} vs check13 stored", got[k], exp[k], rtol=1e-9, atol=1e-9)
    for t in "ab":
        for k in ("mean", "p5", "p95", "pct"):
            CM.gate_equal(f"{spec} D5 control {t} {k} vs check13 stored", got["controls"][t][k], exp["controls"][t][k], rtol=1e-9, atol=1e-9)
        if got["controls"][t]["gt95"] != exp["controls"][t]["gt95"]:
            raise CM.GateFailure(f"GATE FAILED {spec} D5 control {t} gt95 {got['controls'][t]['gt95']} vs check13 {exp['controls'][t]['gt95']}")


def covariate_gate(spec: str, ctx: CovCtx, obs) -> float:
    st = obs["rna_per_cell_equiv"].to_numpy(float)
    m = float(np.max(np.abs(ctx.base_pc - st) / np.maximum(np.abs(st), 1e-300)))
    if not m <= REL_TOL_COV:
        raise CM.GateFailure(f"GATE FAILED {spec} rna_per_cell_equiv reconstruction: max relative diff {m!r} > {REL_TOL_COV}")
    return m


rna_total = H.rna_total      # one definition for every stage (src/hurdle.py)


def rna_total_gate(spec: str, ctx: CovCtx, obs) -> float:
    """rna_total equals the row sum of layer E over the genes not in keratins + MDK (relative 1e-6), else fail closed. Returns the max relative difference."""
    got, ref = rna_total(obs), ctx.base_raw
    m = float(np.max(np.abs(got - ref) / np.maximum(np.abs(ref), 1e-300)))
    if not m <= RNA_TOTAL_REL_TOL:
        raise CM.GateFailure(f"GATE FAILED {spec} rna_total vs row sum of E over non-keratin non-MDK genes: max relative diff {m!r} > {RNA_TOTAL_REL_TOL}")
    return m


# --------------------------------------------------------------------------- one specimen
def run_specimen(spec: str, a, b9: dict, b17: dict, c13: dict) -> dict:
    t0 = time.time()
    d = C14.load(spec, a.comp_run)
    X, E, genes, obs, n = d["X"], d["E"], list(d["genes"]), d["obs"], d["n"]
    print(f"  {spec}: loaded n={n} {time.time() - t0:.0f}s", flush=True)
    res = {"specimen": spec, "n": n}
    if n < R.MIN_EL or n > M.DENSE_EIGH_MAX_N:
        return {**res, "status": "not_estimable", "reason": f"n = {n}"}
    W = knn_weights(np.asarray(d["xy"], float), R.K)
    sn = M.SpectralNull(W)
    mdk_i = genes.index(MDK)
    sec_idx = [genes.index(g) for g in SEC]
    depth = np.log(obs["spot_umi"].to_numpy(float))
    ctx = CovCtx(E, obs["n_st_eff"].to_numpy(float), obs["depth_ratio"].to_numpy(float), genes)
    res["covariate_reconstruction_max_rel_diff"] = covariate_gate(spec, ctx, obs)
    res["rna_total_max_rel_diff"] = rna_total_gate(spec, ctx, obs)
    x0, y0 = X[:, mdk_i], C11.score(X[:, sec_idx])
    stored, stc = b9["per_specimen"][spec], b17["per_specimen"][spec]
    pool_dec, pool = CP.prepare_pool([X], genes, X, mdk_i, sec_idx)
    dr = CP.make_draws(genes, pool_dec, pool, mdk_i, sec_idx, spec + "|cg_alone", n_ctrl=a.n_ctrl, n_uni=1)
    pair_real = [mdk_i] + sec_idx
    projs, form_err = {}, {}
    for fm in FORMS17:
        try:
            if fm in ("none", "D5"):
                projs[fm] = C13.make_projector13(fm, obs, depth)
            else:
                projs[fm] = pair_basis(ctx, depth, "x" if fm == "D5x" else "xc", pair_real)
                if projs[fm] is None:
                    form_err[fm] = "covariate not finite and > 0"
        except H.CovariateDomainError as e:
            projs[fm], form_err[fm] = None, str(e)
    res["forms"] = {}
    for fm in FORMS17:
        if fm in form_err:
            res["forms"][fm] = {"status": "not_estimable", "reason": form_err[fm]}
            continue
        Q = projs[fm]
        x, y = C11.apply(Q, x0), C11.apply(Q, y0)
        r = C11.moran_pair(x, y, W, sn, R.SEED, R.B)
        out = {**r}
        if r["status"] == "ok":
            if fm in ("none", "D5"):
                Ia, Ib = C11.control_I(X, x, y, dr["a"], dr["b"], W, Q)
            else:
                Ia, Ib = control_I_pairwise(X, x0, y0, dr["a"], dr["b"], W, ctx, depth, "x" if fm == "D5x" else "xc", sec_idx, mdk_i)
            out["controls"] = {"a": C10.ctrl_block(r["I"], Ia), "b": C10.ctrl_block(r["I"], Ib), "n_a_nan": int(np.isnan(Ia).sum()), "n_b_nan": int(np.isnan(Ib).sum())}
            if fm == "none":
                out["gate_b9"] = CM.gate_equal(f"{spec} B9 I (no removal) vs stored", r["I"], stored["I"], rtol=C11.TOL, atol=C11.TOL)
                if r["p"] != stored["p_two_sided"]:
                    raise CM.GateFailure(f"GATE FAILED {spec} B9 p_two_sided (no removal): got {r['p']!r} expected {stored['p_two_sided']!r}")
                if a.n_ctrl == C11.N_CTRL:
                    out["gate_b17"] = C11.stored_ctrl_gate(spec, r["I"], Ia, Ib, stc)
            out["depth"] = {nm: {"pearson_r": C11.pearson(depth, v), **C11.moran_pair(depth, v, W, sn, R.SEED, R.B)} for nm, v in (("x", x), ("y", y))}
        res["forms"][fm] = out
    gate_vs_check13(spec, res, c13)
    f5, fc = res["forms"]["D5"], res["forms"]["D5xc"]
    res["d5xc_minus_d5_I"] = float(fc["I"] - f5["I"]) if f5.get("status") == "ok" and fc.get("status") == "ok" else None
    print(f"  {spec}: forms and controls done {time.time() - t0:.0f}s", flush=True)
    old = np.log(obs["rna_per_cell_equiv"].to_numpy(float))
    cx, cxc = pair_cov(ctx, "x", pair_real), pair_cov(ctx, "xc", pair_real)
    share = ctx.pc[:, sec_idx].sum(axis=1) / ctx.base_pc
    res["covariate"] = {"spearman_old_vs_D5x": float(spearmanr(old, np.log(cx))[0]) if (cx > 0).all() else None,
                        "spearman_old_vs_D5xc": float(spearmanr(old, np.log(cxc))[0]) if (cxc > 0).all() else None,
                        "secretory_share_of_old_median": float(np.median(share)), "secretory_share_of_old_min": float(share.min()), "secretory_share_of_old_max": float(share.max())}
    targets = {"S1": b9["aggregate"]["mean"], "S2": 2 * b9["aggregate"]["mean"], "S3": stored["I"] if stored["I"] > 0 else None}
    res["power"] = C11.power_specimen(spec, X, genes, depth, W, sn, projs, targets, a.n_real, a.b_power)
    res["status"] = "ok"
    print(f"{spec}: n={n} gates passed; I " + " ".join(f"{fm} {res['forms'][fm].get('I', float('nan')):.4f}" for fm in FORMS17) + f" {time.time() - t0:.0f}s", flush=True)
    return res


# --------------------------------------------------------------------------- tables
def tables(meta, per, coh, b9) -> str:
    f = C11.f
    F = FORMS17
    L = [f"# Check 17 generated tables (n_real {meta['n_real']}, B power {meta['b_power']}, controls {meta['n_ctrl']})", "", f"Output dir `{meta['out']}`; commit {meta['git_commit'][:8]}; gates passed for {meta['n_gates']} specimens.", "",
         "## Reading rule", "", "```", json.dumps(coh["reading"], indent=1, default=float), "```", "",
         "## Part 1: I (MSR p, BH q within form)", "", "| specimen | n | stored I (q) | " + " | ".join(F) + " |", "|---|---|---|" + "---|" * len(F)]
    for s, r in per.items():
        if r.get("status") != "ok":
            L.append(f"| {s} | {r.get('n')} | not estimable: {r.get('reason')} |" + "|" * len(F))
            continue
        cell = lambda v: "not est." if v.get("status") != "ok" else f"{f(v['I'])} (p {f(v['p'], 3)}, q {f(v.get('q_bh'))})"  # noqa: E731
        L.append(f"| {s} | {r['n']} | {f(b9['per_specimen'][s]['I'])} ({f(b9['bh'][s]['q'])}) | " + " | ".join(cell(r["forms"][fm]) for fm in F) + " |")
    L += ["", "## Old vs new covariate", "", "| specimen | Spearman old~D5x | Spearman old~D5xc | secretory share of old (median) | min | max |", "|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") == "ok":
            c = r["covariate"]
            L.append(f"| {s} | {f(c['spearman_old_vs_D5x'], 3)} | {f(c['spearman_old_vs_D5xc'], 3)} | {f(c['secretory_share_of_old_median'], 4)} | {f(c['secretory_share_of_old_min'], 4)} | {f(c['secretory_share_of_old_max'], 4)} |")
    L += ["", "## Part 2: controls per form", "", "| specimen | form | real I | a mean | a p5 | a p95 | a pct | b mean | b p5 | b p95 | b pct | >p95 a | >p95 b |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for fm in F:
            v = r["forms"][fm]
            if v.get("status") == "ok":
                c = v["controls"]
                L.append(f"| {s} | {fm} | {f(v['I'])} | {f(c['a']['mean'])} | {f(c['a']['p5'])} | {f(c['a']['p95'])} | {f(c['a']['pct'], 1)} | {f(c['b']['mean'])} | {f(c['b']['p5'])} | {f(c['b']['p95'])} | {f(c['b']['pct'], 1)} | {c['a']['gt95']} | {c['b']['gt95']} |")
    L += ["", "## Cohort per form", "", "| form | n | mean real I | mean ctrl a | mean ctrl b | excess a | excess b | n >p95 a | n >p95 b | n q<0.05 | of 6 B9 sig remaining |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for fm in F:
        c = coh["forms"][fm]
        if "R2" not in c:
            L.append(f"| {fm} | not estimable |||||||||")
            continue
        r2 = c["R2"]
        L.append(f"| {fm} | {c['n_estimable']} | {f(r2['mean_real_I'])} | {f(r2['mean_ctrl_a_mean'])} | {f(r2['mean_ctrl_b_mean'])} | {f(r2['mean_real_I'] - r2['mean_ctrl_a_mean'])} | {f(r2['mean_real_I'] - r2['mean_ctrl_b_mean'])} | {c['R3']['n_gt_p95_a']} | {c['R3']['n_gt_p95_b']} | {coh['n_q_lt_0.05'][fm]} | {c['R1']['n_remaining']} |")
    L += ["", "## Part 1b: depth after removal", "", "| specimen | form | r(depth,x) | I(depth,x) | q | r(depth,y) | I(depth,y) | q |", "|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for fm in F:
            v = r["forms"][fm]
            if v.get("status") == "ok":
                dx, dy = v["depth"]["x"], v["depth"]["y"]
                L.append(f"| {s} | {fm} | {f(dx['pearson_r'])} | {f(dx.get('I'))} | {f(dx.get('q_bh'))} | {f(dy['pearson_r'])} | {f(dy.get('I'))} | {f(dy.get('q_bh'))} |")
    L += ["", "## Power: mean I [combined p] (frac p<0.05)", "", "| specimen | scenario | " + " | ".join(F) + " |", "|---|---|" + "---|" * len(F)]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for sc in C11.SCEN:
            blk = r["power"]["scenarios"][sc]
            if "forms" not in blk:
                L.append(f"| {s} | {sc} | {blk.get('status')} |" + "|" * (len(F) - 1))
                continue
            cells = [f"not est. ({blk['forms'][fm].get('n_valid', 0)} valid)" if blk["forms"][fm]["status"] != "ok" else
                     f"{f(blk['forms'][fm]['mean'])} [{f(blk['forms'][fm]['p_combined'], 3)}] ({f(blk['forms'][fm]['frac_p_lt_0.05'], 2)})" for fm in F]
            L.append(f"| {s} | {sc} | " + " | ".join(cells) + " |")
    L += ["", "## Cohort (full)", "", "```", json.dumps(coh, indent=1, default=float, sort_keys=True), "```", ""]
    return "\n".join(L)


@contextlib.contextmanager
def patched_forms():
    """check 11's power_specimen / cohort / csv writers read the module global ``FORMS``: swap in this check's forms and restore on exit."""
    keep = C11.FORMS
    C11.FORMS = FORMS17
    try:
        yield
    finally:
        C11.FORMS = keep


def compute(comp_run, b9: dict, b17: dict, c13: dict, specs, out, meta_common: dict, n_real: int = C11.N_REAL, n_ctrl: int = C11.N_CTRL, b_power: int = R.B) -> tuple:
    """Check 17 (was ``main``). ``c13`` = the check 13 summary computed earlier in the same run (``{"per_specimen": ...}``); D5 is gated against it (status, I, p and the control blocks, 1e-9).
    Writes check17_{part1,controls,depthcheck,power}.csv and check17_tables.md into ``out``. Returns ``(summary, files)``; raises ``GateFailure``."""
    a = types.SimpleNamespace(n_real=n_real, n_ctrl=n_ctrl, b_power=b_power, comp_run=comp_run)
    per = {}
    with patched_forms():
        for s in specs:
            per[s] = run_specimen(s, a, b9, b17, c13)
        coh = C11.cohort(per, b9, b17)
        est = [r for r in per.values() if r.get("status") == "ok"]
        coh["n_q_lt_0.05"] = {fm: int(sum(r["forms"][fm].get("q_bh", 1) < 0.05 for r in est if r["forms"][fm].get("status") == "ok")) for fm in FORMS17}

        def num(fm):
            c = coh["forms"][fm]
            return {"real": c["R2"]["mean_real_I"], "ctrl_a": c["R2"]["mean_ctrl_a_mean"], "ctrl_b": c["R2"]["mean_ctrl_b_mean"], "n_remaining": c["R1"]["n_remaining"]}

        coh["reading"] = label17(num("D5x"), num("D5xc"))
        sp = {k: [r["covariate"][k] for r in est if r["covariate"][k] is not None] for k in ("spearman_old_vs_D5x", "spearman_old_vs_D5xc", "secretory_share_of_old_median")}
        coh["covariate"] = {k: {"mean": float(np.mean(v)), "min": float(np.min(v)), "max": float(np.max(v))} for k, v in sp.items() if v}
        dd = [r["d5xc_minus_d5_I"] for r in est if r.get("d5xc_minus_d5_I") is not None]
        coh["d5xc_minus_d5"] = {"n": len(dd), "mean_I_diff": float(np.mean(dd)) if dd else None, "max_abs_I_diff": float(np.max(np.abs(dd))) if dd else None,
                                "per_specimen": {r["specimen"]: r.get("d5xc_minus_d5_I") for r in est},
                                "note": "D5xc = D5 without the pair's own genes in the covariate (about 0.2% of it); the difference is the self-inclusion effect on the new D5"}
        meta = {"check": "check17_d5_fixed_moran", "n_real": n_real, "n_ctrl": n_ctrl, "b_power": b_power, "seed": C11.SEED11, "k": R.K, "n_gates": len(est), **meta_common}
        out = Path(out)
        with tempfile.TemporaryDirectory(dir=out) as tmp:                                  # check 11's writer names its files check11_*
            C11.write_csvs(Path(tmp), per)
            for nm in ("part1", "controls", "depthcheck", "power"):
                os.replace(Path(tmp) / f"check11_{nm}.csv", out / f"check17_{nm}.csv")
        (out / "check17_tables.md").write_text(tables(meta, per, coh, b9))
    return {**meta, "per_specimen": per, "cohort": coh}, [f"check17_{nm}.csv" for nm in ("part1", "controls", "depthcheck", "power")] + ["check17_tables.md"]
