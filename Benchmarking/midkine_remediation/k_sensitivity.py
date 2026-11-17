"""B22 part 6 (check 15): neighbourhood size k in {6, 15, 18} for the MDK-secretory bivariate Moran's I and its matched controls (forms none, S_B, D1), with a
planted-association power check. k=6 stays the primary and is frozen in the protocol. Descriptive only.

Promoted from the exploratory check 15 (``exploratory_v42/check15_k_sensitivity.py``, git history). The analysis functions are that code, unchanged (``seed_for``, ``build_W``,
``neighbour_distance``, ``label15``, ``plant_realisation``, ``_xy_of``, ``planted_I``, ``power_form``, ``run_form``, ``gate_against``, ``cohort``, ``f``, ``tables``, ``write_csvs``;
tests/midkine_remediation/test_moran_robustness.py checks the AST). Changed: ``run_specimen`` reads the B5 run from ``a.comp_run`` (it was a module constant), the two reference
summary files the exploratory check gated against (check 14, check 11) are the check 14 and check 11 results computed earlier IN THE SAME RUN (arguments ``ref14`` / ``ref11``),
``main`` became ``compute``, the table header no longer prints a dirty flag. Rules: ``exploratory_v42/check15_DECISION_RULE.md`` (git history).
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
import control_pairs as CP  # noqa: E402
import depth_moran as C10  # noqa: E402
import depth_removed_moran as C11  # noqa: E402
import hurdle as H  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
import share_moran as C14  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

SEED15 = 20261115
KS = (6, 15, 18)
FORMS15 = ("none", "S_B", "D1")
LABEL_FORMS = ("none", "S_B")
SCEN15 = ("S0", "S1")
S1_SIZE = 0.048217128908550164
SEC, MDK = list(R.SECRETORY), C14.MDK
ARMS = ("a", "b")


def seed_for(label: str) -> int:
    return int(np.random.SeedSequence([SEED15, zlib.crc32(label.encode())]).generate_state(1)[0])


def build_W(xy, k: int):
    """The k-specific graph (row-standardised KNN, self excluded)."""
    return knn_weights(np.asarray(xy, float), int(k))


def neighbour_distance(xy, k: int, pitch: float) -> dict:
    """Distance to the k-th nearest other compartment, in pitch units: median and 95th percentile over compartments."""
    d, _ = cKDTree(np.asarray(xy, float)).query(np.asarray(xy, float), k=k + 1)
    v = d[:, k] / pitch
    return {"median": float(np.median(v)), "p95": float(np.percentile(v, 95)), "max": float(v.max())}


# --------------------------------------------------------------------------- reading rule
def label15(coh: dict) -> dict:
    """coh[form][k] = {'gt_a','gt_b','ratio_a','ratio_b'} (ratio None if undefined). WIDER_HELPS if some k in (15, 18) beats k=6 in BOTH none and S_B on all four counts."""
    detail = {}
    for k in (15, 18):
        ok_k = True
        for fm in LABEL_FORMS:
            c, b = coh[fm][k], coh[fm][6]
            conds = {"gt_a": c["gt_a"] >= b["gt_a"] + 2, "gt_b": c["gt_b"] >= b["gt_b"] + 2,
                     "ratio_a": c["ratio_a"] is not None and b["ratio_a"] is not None and c["ratio_a"] > b["ratio_a"],
                     "ratio_b": c["ratio_b"] is not None and b["ratio_b"] is not None and c["ratio_b"] > b["ratio_b"]}
            detail[f"k{k}_{fm}"] = {kk: bool(v) for kk, v in conds.items()}
            ok_k &= all(conds.values())
        detail[f"k{k}_all"] = bool(ok_k)
    helps = bool(detail["k15_all"] or detail["k18_all"])
    return {"label": "WIDER_HELPS" if helps else "WIDER_NO_GAIN", "detail": detail}


# --------------------------------------------------------------------------- power
def plant_realisation(spec: str, N0, W, j: int, Qf=None) -> tuple:
    n = N0.shape[0]
    rng = np.random.default_rng(np.random.SeedSequence([SEED15, zlib.crc32((spec + "|plant").encode()), j]))
    P = N0.copy()
    for c in range(P.shape[1]):
        P[:, c] = P[rng.permutation(n), c]
    f = W @ (W @ (W @ rng.normal(size=n)))
    f = C11.apply(Qf, f)
    return P, (f - f.mean()) / f.std()


def _xy_of(out, Q):
    return C11.apply(Q, out[:, 0]), C11.apply(Q, C10.mean_z_cols(out[:, 1:]))


def planted_I(P, f, b, W, Q) -> float:
    x, y = _xy_of(P * np.exp(b * f)[:, None], Q)
    return float(M.morans_bv(x, y, W))


def power_form(spec: str, N0, W, sn, Q, n_real: int, B: int, tag: str) -> dict:
    """S0 (nothing planted) and S1 (stored-size plant) for one form at one k. Q = None except D1 (projection applied after planting, f made orthogonal to the design)."""
    rec = {sc: [] for sc in SCEN15}
    bs = []
    for j in range(n_real):
        P, f = plant_realisation(spec, N0, W, j, Q)
        for sc in SCEN15:
            if sc == "S0":
                b, ok = 0.0, True
            else:
                cal = C11.bisect(lambda bb: planted_I(P, f, bb, W, Q), S1_SIZE)
                b, ok = cal["param"], cal["reachable"]
                bs.append(b)
            if not ok:
                rec[sc].append(C10.rec({}, "target not reachable"))
                continue
            x, y = _xy_of(P * np.exp(b * f)[:, None], Q)
            rec[sc].append(C10.rec(C11.moran_pair(x, y, W, sn, seed_for(f"{spec}|{tag}|{sc}|j{j}"), B)))
    return {"b_mean": float(np.mean(bs)) if bs else None, **{sc: C10.agg(rec[sc]) for sc in SCEN15}}


# --------------------------------------------------------------------------- one form at one k
def run_form(fm: str, k: int, d: dict, W, sn, dr_, stored_x) -> tuple:
    """Returns (record, N0 columns of the form for the power plant, projector or None)."""
    X, E, tot, genes, obs = stored_x["X"], stored_x["E"], stored_x["tot"], stored_x["genes"], d["obs"]
    mdk_i, sec_idx = stored_x["mdk_i"], stored_x["sec_idx"]
    Q = None
    if fm == "none":
        cols = np.column_stack([X[:, mdk_i], X[:, sec_idx]])
        x, y = cols[:, 0], C10.mean_z_cols(cols[:, 1:])
    elif fm == "S_B":
        cols = C14.real_cols("S_B", None, E, tot, mdk_i, sec_idx)
        x, y = cols[:, 0], C10.mean_z_cols(cols[:, 1:])
    else:
        cols = np.column_stack([X[:, mdk_i], X[:, sec_idx]])
        Q = C11.make_projector("D1", obs, np.log(obs["spot_umi"].to_numpy(float)))
        x, y = C11.apply(Q, cols[:, 0]), C11.apply(Q, C10.mean_z_cols(cols[:, 1:]))
    r = C11.moran_pair(x, y, W, sn, R.SEED, R.B)
    out = {**r}
    if r["status"] == "ok":
        if fm == "S_B":
            Ia, Ib = C14.control_I_share_b(E, tot, mdk_i, sec_idx, dr_["a"], dr_["b"], W)
        else:
            Ia, Ib = C11.control_I(X, x, y, dr_["a"], dr_["b"], W, Q)
        out["controls"] = {"a": C10.ctrl_block(r["I"], Ia), "b": C10.ctrl_block(r["I"], Ib), "n_a_nan": int(np.isnan(Ia).sum()), "n_b_nan": int(np.isnan(Ib).sum())}
        out["_Ia"], out["_Ib"] = Ia, Ib
    return out, cols, Q


def gate_against(spec: str, name: str, got: dict, ref: dict, check_ctrl: bool) -> dict:
    """Fail closed: I (1e-9), p exactly, and the control a/b blocks, against a reference output."""
    CM.gate_equal(f"{spec} {name} I", got["I"], ref["I"], rtol=C11.TOL, atol=C11.TOL)
    if got["p"] != ref["p"]:
        raise CM.GateFailure(f"GATE FAILED {spec} {name} p: got {got['p']!r} expected {ref['p']!r}")
    if check_ctrl:
        for t in ARMS:
            for fld in ("mean", "p5", "p95", "pct"):
                CM.gate_equal(f"{spec} {name} control {t} {fld}", got["controls"][t][fld], ref["controls"][t][fld], rtol=C11.TOL, atol=C11.TOL)
            if got["controls"][t]["gt95"] != ref["controls"][t]["gt95"]:
                raise CM.GateFailure(f"GATE FAILED {spec} {name} control {t} gt95: got {got['controls'][t]['gt95']} expected {ref['controls'][t]['gt95']}")
    return {"name": name, "I": got["I"], "p": got["p"], "ctrl_checked": bool(check_ctrl)}


def run_specimen(spec: str, a, b9: dict, b17: dict, ref14: dict, ref11: dict) -> dict:
    t0 = time.time()
    d = C14.load(spec, a.comp_run)
    n, genes, obs = d["n"], list(d["genes"]), d["obs"]
    print(f"  {spec}: loaded n={n} {time.time() - t0:.0f}s", flush=True)
    res = {"specimen": spec, "n": n}
    if n < R.MIN_EL or n > M.DENSE_EIGH_MAX_N:
        return {**res, "status": "not_estimable", "reason": f"n = {n}"}
    X, E = d["X"], np.asfortranarray(d["E"])
    sx = {"X": X, "E": E, "tot": E.sum(axis=1), "genes": genes, "mdk_i": genes.index(MDK), "sec_idx": [genes.index(g) for g in SEC]}
    stored, stc = b9["per_specimen"][spec], b17["per_specimen"][spec]
    pool_dec, pool = CP.prepare_pool([X], genes, X, sx["mdk_i"], sx["sec_idx"])
    dr_ = CP.make_draws(genes, pool_dec, pool, sx["mdk_i"], sx["sec_idx"], spec + "|cg_alone", n_ctrl=a.n_ctrl, n_uni=1)
    xy = np.asarray(d["xy"], float)
    pitch = H.median_pitch(xy)
    res.update({"pitch": pitch, "gates": {}, "k": {}})
    full = a.n_ctrl == C11.N_CTRL
    for k in KS:
        W = build_W(xy, k)
        sn = M.SpectralNull(W)
        rk = {"distance_pitch": neighbour_distance(xy, k, pitch), "forms": {}, "power": {}}
        for fm in FORMS15:
            out, cols, Q = run_form(fm, k, d, W, sn, dr_, sx)
            if k == 6:                                     # fidelity gates
                if out["status"] != "ok":
                    raise CM.GateFailure(f"GATE FAILED {spec} {fm} k=6 not estimable")
                if fm == "none":
                    res["gates"]["none_b9"] = gate_against(spec, "none vs stored B9", out, {"I": stored["I"], "p": stored["p_two_sided"]}, False)
                    if full:
                        res["gates"]["none_b17"] = C11.stored_ctrl_gate(spec, out["I"], out["_Ia"], out["_Ib"], stc)
                elif fm == "S_B":
                    res["gates"]["S_B_check14"] = gate_against(spec, "S_B vs check 14", out, ref14["per_specimen"][spec]["forms"]["S_B"], full)
                else:
                    res["gates"]["D1_check11"] = gate_against(spec, "D1 vs check 11", out, ref11["per_specimen"][spec]["forms"]["D1"], full)
            out.pop("_Ia", None)
            out.pop("_Ib", None)
            rk["forms"][fm] = out
            rk["power"][fm] = power_form(spec, cols, W, sn, Q, a.n_real, a.b_power, f"k{k}|{fm}")
        res["k"][str(k)] = rk
        print(f"  {spec}: k={k} done {time.time() - t0:.0f}s; I " + " ".join(f"{fm} {rk['forms'][fm].get('I', float('nan')):.4f}" for fm in FORMS15), flush=True)
    res["status"] = "ok"
    return res


# --------------------------------------------------------------------------- cohort and tables
def cohort(per: dict, b9: dict) -> dict:
    est = {s: r for s, r in per.items() if r.get("status") == "ok"}
    sig = b9["bh"]["significant_q_lt_0.05"]
    out = {"stored_q_lt_0.05": sig, "forms": {}}
    for fm in FORMS15:
        out["forms"][fm] = {}
        for k in KS:
            ok = {s: r["k"][str(k)]["forms"][fm] for s, r in est.items() if r["k"][str(k)]["forms"][fm].get("status") == "ok"}
            q = R._bh([v["p"] for v in ok.values()]) if ok else []
            for (s, v), qi in zip(ok.items(), q):
                v["q_bh"] = float(qi)
            for sc in SCEN15:
                C10.bh_attach(est, lambda r, k=k, fm=fm, sc=sc: r["k"][str(k)]["power"][fm][sc])
            mr = float(np.mean([v["I"] for v in ok.values()])) if ok else None
            cm = {t: float(np.mean([v["controls"][t]["mean"] for v in ok.values()])) for t in ARMS} if ok else {}
            rat = {t: (mr / cm[t] if ok and cm[t] > 0 else None) for t in ARMS}
            pw = {sc: [est[s]["k"][str(k)]["power"][fm][sc] for s in est] for sc in SCEN15}
            out["forms"][fm][k] = {
                "n_estimable": len(ok), "mean_real_I": mr, "mean_ctrl_a": cm.get("a"), "mean_ctrl_b": cm.get("b"),
                "ratio_a": rat["a"], "ratio_b": rat["b"], "diff_a": (mr - cm["a"]) if ok else None, "diff_b": (mr - cm["b"]) if ok else None,
                "gt_a": int(sum(v["controls"]["a"]["gt95"] for v in ok.values())), "gt_b": int(sum(v["controls"]["b"]["gt95"] for v in ok.values())),
                "mean_ctrl_a_p95": float(np.mean([v["controls"]["a"]["p95"] for v in ok.values()])) if ok else None,
                "mean_ctrl_b_p95": float(np.mean([v["controls"]["b"]["p95"] for v in ok.values()])) if ok else None,
                "n_q_lt_0.05": int(sum(v["q_bh"] < 0.05 for v in ok.values())), "q_lt_0.05": [s for s, v in ok.items() if v["q_bh"] < 0.05],
                "stored_sig_remaining": int(sum(s in ok and ok[s]["q_bh"] < 0.05 for s in sig)),
                "power": {sc: {"n_combined_p_lt_0.05": int(sum(v.get("status") == "ok" and v["p_combined"] < 0.05 for v in pw[sc])),
                               "n_estimable": int(sum(v.get("status") == "ok" for v in pw[sc]))} for sc in SCEN15}}
    out["reading"] = label15({fm: {k: out["forms"][fm][k] for k in KS} for fm in LABEL_FORMS})
    out["distance"] = {str(k): {"median_of_medians": float(np.median([r["k"][str(k)]["distance_pitch"]["median"] for r in est.values()])),
                                "max_of_p95": float(max(r["k"][str(k)]["distance_pitch"]["p95"] for r in est.values()))} for k in KS}
    return out


def f(x, nd=4):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def tables(meta, per, coh, b9) -> str:
    L = [f"# Check 15 generated tables (n_real {meta['n_real']}, B power {meta['b_power']}, controls {meta['n_ctrl']})", "", f"Output dir `{meta['out']}`; commit {meta['git_commit'][:8]}; gates passed for {meta['n_gates']} specimens.", "",
         "## Reading rule", "", "```", json.dumps(coh["reading"], indent=1, default=float), "```", "", "## Cohort per form and k", "",
         "| form | k | mean real I | ctrl a mean | ctrl b mean | ratio a | ratio b | diff a | diff b | >p95 a | >p95 b | q<0.05 | stored-sig remaining | S0 FP | S1 recovered |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for fm in FORMS15:
        for k in KS:
            c = coh["forms"][fm][k]
            L.append(f"| {fm} | {k} | {f(c['mean_real_I'])} | {f(c['mean_ctrl_a'])} | {f(c['mean_ctrl_b'])} | {f(c['ratio_a'], 3)} | {f(c['ratio_b'], 3)} | {f(c['diff_a'])} | {f(c['diff_b'])} | {c['gt_a']} | {c['gt_b']} | {c['n_q_lt_0.05']} | {c['stored_sig_remaining']} | "
                     f"{c['power']['S0']['n_combined_p_lt_0.05']} of {c['power']['S0']['n_estimable']} | {c['power']['S1']['n_combined_p_lt_0.05']} of {c['power']['S1']['n_estimable']} |")
    for fm in FORMS15:
        L += ["", f"## Per specimen, form {fm}: I (MSR p, BH q) and control a / b mean [p95] (real > p95)", "", "| specimen | n | k | I (p, q) | ctrl a mean [p95] | >p95 a | ctrl b mean [p95] | >p95 b |", "|---|---|---|---|---|---|---|---|"]
        for s, r in per.items():
            if r.get("status") != "ok":
                L.append(f"| {s} | {r.get('n')} | | not estimable: {r.get('reason')} | | | | |")
                continue
            for k in KS:
                v = r["k"][str(k)]["forms"][fm]
                if v.get("status") != "ok":
                    L.append(f"| {s} | {r['n']} | {k} | not est. | | | | |")
                    continue
                c = v["controls"]
                L.append(f"| {s} | {r['n']} | {k} | {f(v['I'])} ({f(v['p'], 3)}, {f(v.get('q_bh'))}) | {f(c['a']['mean'])} [{f(c['a']['p95'])}] | {c['a']['gt95']} | {f(c['b']['mean'])} [{f(c['b']['p95'])}] | {c['b']['gt95']} |")
    L += ["", "## Realised k-th neighbour distance in pitch units (median / 95th percentile / max over compartments)", "", "| specimen | n | pitch | " + " | ".join(f"k={k}" for k in KS) + " |", "|---|---|---|" + "---|" * len(KS)]
    for s, r in per.items():
        if r.get("status") == "ok":
            L.append(f"| {s} | {r['n']} | {f(r['pitch'], 3)} | " + " | ".join(f"{f(r['k'][str(k)]['distance_pitch']['median'], 3)} / {f(r['k'][str(k)]['distance_pitch']['p95'], 3)} / {f(r['k'][str(k)]['distance_pitch']['max'], 3)}" for k in KS) + " |")
    L += ["", "## Power: mean I [combined p] (frac p<0.05); b = mean planted strength for S1", "", "| specimen | k | scenario | " + " | ".join(FORMS15) + " |", "|---|---|---|" + "---|" * len(FORMS15)]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for k in KS:
            for sc in SCEN15:
                cells = []
                for fm in FORMS15:
                    v = r["k"][str(k)]["power"][fm][sc]
                    cells.append(f"not est. ({v.get('n_valid', 0)} valid)" if v["status"] != "ok" else f"{f(v['mean'])} [{f(v['p_combined'], 3)}] ({f(v['frac_p_lt_0.05'], 2)})")
                L.append(f"| {s} | {k} | {sc} | " + " | ".join(cells) + " |")
    L += ["", "## Cohort (json)", "", "```", json.dumps(coh, indent=1, default=float, sort_keys=True), "```", ""]
    return "\n".join(L)


def write_csvs(out, per):
    def wr(name, header, rows):
        with open(out / name, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)

    est = {s: r for s, r in per.items() if r.get("status") == "ok"}
    wr("check15_part1.csv", ["specimen", "n", "k", "form", "status", "I", "p", "q_bh"], [[s, r["n"], k, fm, v.get("status"), v.get("I"), v.get("p"), v.get("q_bh")] for s, r in est.items() for k in KS for fm in FORMS15 for v in [r["k"][str(k)]["forms"][fm]]])
    wr("check15_controls.csv", ["specimen", "k", "form", "real_I"] + [f"{t}_{x}" for t in ARMS for x in ("mean", "p5", "p95", "pct", "gt95")],
       [[s, k, fm, v["I"]] + [v["controls"][t][x] for t in ARMS for x in ("mean", "p5", "p95", "pct", "gt95")] for s, r in est.items() for k in KS for fm in FORMS15 for v in [r["k"][str(k)]["forms"][fm]] if v.get("status") == "ok"])
    wr("check15_distance.csv", ["specimen", "n", "pitch", "k", "median_pitch_units", "p95_pitch_units", "max_pitch_units"], [[s, r["n"], r["pitch"], k, *(r["k"][str(k)]["distance_pitch"][x] for x in ("median", "p95", "max"))] for s, r in est.items() for k in KS])
    wr("check15_power.csv", ["specimen", "k", "form", "scenario", "status", "n_valid", "mean_I", "p_combined", "frac_p_lt_0.05", "q_bh"],
       [[s, k, fm, sc, v["status"], v.get("n_valid"), v.get("mean"), v.get("p_combined"), v.get("frac_p_lt_0.05"), v.get("q_bh")] for s, r in est.items() for k in KS for fm in FORMS15 for sc in SCEN15 for v in [r["k"][str(k)]["power"][fm][sc]]])


def compute(comp_run, b9: dict, b17: dict, ref14: dict, ref11: dict, specs, out, meta_common: dict, n_real: int = C11.N_REAL, n_ctrl: int = C11.N_CTRL, b_power: int = R.B) -> tuple:
    """Check 15 (was ``main``). ``ref14`` / ``ref11`` = the check 14 / check 11 summaries computed earlier in the same run (``{"per_specimen": ...}``); the k=6 blocks are gated against
    them and against the stored B9 / B17. Writes check15_{part1,controls,distance,power}.csv and check15_tables.md into ``out``. Returns ``(summary, files)``; raises ``GateFailure``."""
    a = types.SimpleNamespace(n_real=n_real, n_ctrl=n_ctrl, b_power=b_power, comp_run=comp_run)
    if not np.isclose(b9["aggregate"]["mean"], S1_SIZE, rtol=0, atol=1e-15):
        raise CM.GateFailure(f"stored cohort mean I {b9['aggregate']['mean']!r} != fixed S1 size {S1_SIZE!r}")
    per = {}
    for s in specs:
        per[s] = run_specimen(s, a, b9, b17, ref14, ref11)
    coh = cohort(per, b9)
    meta = {"check": "check15_k_sensitivity", "n_real": n_real, "n_ctrl": n_ctrl, "b_power": b_power, "seed": SEED15, "ks": list(KS),
            "n_gates": sum(r.get("status") == "ok" for r in per.values()), **meta_common}
    out = Path(out)
    write_csvs(out, per)
    (out / "check15_tables.md").write_text(tables(meta, per, coh, b9))
    return {**meta, "per_specimen": per, "cohort": coh}, [f"check15_{nm}.csv" for nm in ("part1", "controls", "distance", "power")] + ["check15_tables.md"]
