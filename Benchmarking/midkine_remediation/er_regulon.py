"""B21 part 2 (check 8): ESR1 regulon activity (ULM) vs MDK on eligible cancer compartments. Descriptive only.

Promoted from the exploratory check 8 (``exploratory_v42/check8_er_regulon.py``, git history). The analysis functions are that code, unchanged
(``ulm_t``, ``build_regulon``, ``auc_mwu``, ``variant_block``, ``cohort8``, ``pair_table``, ``flat_rows8`` are AST-identical; tests/midkine_remediation/test_er_regulon.py checks it).
Changed: ``load8`` and ``analyze8`` take the B5 compartments run as the keyword ``comp_run`` (it was a module constant); ``gate_vs_check7`` compares with the check 7
per-specimen records computed earlier in the same stage run (it read the exploratory check 7 summary file). ``main`` is gone: the stage entry is ``e2_covariation.py``.
The saved regulon tables (CollecTRI primary, DoRothEA A-C) are the committed files in ``regulon_inputs/`` (fetched once by ``fetch_regulons.py``; the analysis never
touches the network). Rules: ``exploratory_v42/check8_DECISION_RULE.md`` (git history).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))

import a12_9_inputs as A  # noqa: E402
import a12_10_real as R  # noqa: E402
import compartments as C  # noqa: E402
import control_pairs as CP  # noqa: E402
import hurdle as H  # noqa: E402
import keratin_e2 as C7  # noqa: E402
import keratin_secretory_association as K  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

E2 = C7.E2
REG_DIR = HERE / "regulon_inputs"
B = C7.B
N_CTRL = C7.N_CTRL
MIN_EL = C7.MIN_EL
MIN_NZ = C7.MIN_NZ
RNA_COVARIATE = "total"   # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md)
MIN_TARGETS = 5           # decoupler min_n / tmin default
MIN_B_GROUP = 10
B_MIN_PROP, B_MIN_COUNTS = 0.05, 200
TF = "ESR1"
KERATINS = list(C.KERATINS)
SECRETORY = list(R.SECRETORY)
# name -> (network file, genes removed from the feature set). MDK is always removed from the features (rule file).
VARIANTS = {"collectri": ("collectri_human.tsv", {C.MDK}),
            "dorothea_ABC": ("dorothea_ABC_human.tsv", {C.MDK}),
            "collectri_ex_ks": ("collectri_human.tsv", {C.MDK, *KERATINS, *SECRETORY})}


# --------------------------------------------------------------------------- ULM (pure, tested against decoupler in tests/midkine_remediation/test_er_regulon.py)
def ulm_t(Y, w, chunk: int = 2000) -> np.ndarray:
    """decoupler-style univariate linear model: per row of Y (n, G) OLS of y on [1, w]; returns the slope's t statistic (df = G - 2)."""
    w = np.asarray(w, float)
    G = len(w)
    xc = w - w.mean()
    sxx = float(xc @ xc)
    if sxx <= 0 or G < 3:
        return np.full(Y.shape[0], np.nan)
    out = np.empty(Y.shape[0])
    for i in range(0, Y.shape[0], chunk):
        y = np.asarray(Y[i:i + chunk], float)
        yc = y - y.mean(axis=1, keepdims=True)
        slope = (yc @ xc) / sxx
        sse = np.maximum((yc * yc).sum(axis=1) - slope * slope * sxx, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[i:i + chunk] = slope / np.sqrt(sse / (G - 2) / sxx)
    return out


def build_regulon(net: pd.DataFrame, genes, exclude, min_n: int = MIN_TARGETS) -> dict:
    """As decoupler: the regression runs over EVERY gene of the data (``exclude`` removed), ESR1 weight 0 for non-targets; ESR1 is dropped
    (all-zero weights) when it has fewer than ``min_n`` targets on the panel (decoupler's ``tmin``)."""
    gset = set(genes)
    e = net[net["source"] == TF].drop_duplicates("target").set_index("target")["weight"]
    on_panel = [t for t in e.index if t in gset]
    kept = len(on_panel) >= min_n
    feat_idx = np.array([i for i, g in enumerate(genes) if g not in set(exclude)], dtype=int)
    fn = [genes[i] for i in feat_idx]
    w = np.array([float(e.get(g, 0.0)) if kept else 0.0 for g in fn])
    return {"feat_idx": feat_idx, "feat_names": fn, "w": w, "n_regulon": int(len(e)), "n_on_panel": int(len(on_panel)), "n_measured": int((w != 0).sum()), "targets": list(e.index), "esr1_kept_by_min_n": bool(kept)}


def auc_mwu(a, b) -> dict:
    """P(a > b) (+ half ties) = U / (n1 n0) and the two-sided Mann-Whitney p."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    u = mannwhitneyu(a, b, alternative="two-sided")
    return {"auc": float(u.statistic / (len(a) * len(b))), "p_mwu_two_sided": float(u.pvalue)}


# --------------------------------------------------------------------------- data
def load8(spec: str, *, comp_run) -> dict:
    import anndata as ad

    a = ad.read_h5ad(str(Path(comp_run) / f"{spec}.h5ad"), backed="r")
    obs_all = a.obs
    el = C.eligible_cancer(obs_all).to_numpy().astype(bool)
    mb = ((obs_all["qp_prop"] >= B_MIN_PROP) & (obs_all["total_counts"] >= B_MIN_COUNTS)).to_numpy().astype(bool)
    need = el | mb
    genes = list(a.var_names)
    X = R._dense(a.X[need])
    xy = np.asarray(a.obsm["spatial"], float)[el]
    obs = obs_all[need]
    a.file.close()
    el_in = el[need]
    return {"obs": obs, "obs_el": obs[el_in], "genes": genes, "X": X, "X_el": X[el_in], "el_in": el_in, "mb_in": mb[need], "xy": xy, "n": int(el.sum())}


# --------------------------------------------------------------------------- one specimen
def variant_block(spec, vname, reg, d, act_all, scores, sec, mdk_np, sender, sn, Z, ctrl_cols, ctrl_err, pc, mdk_ok, B) -> dict:
    el_in, mb_in, obs = d["el_in"], d["mb_in"], d["obs"]
    act = act_all[el_in]
    out = {"regulon": {"n_regulon_targets": reg["n_regulon"], "n_measured_on_panel": reg["n_measured"], "n_measured_nonzero_in_eligible": int((d["X_el"][:, reg["feat_idx"][reg["w"] != 0]] > 0).any(axis=0).sum()),
                       "n_regression_features": int(len(reg["feat_idx"])), "n_nan_activity": int((~np.isfinite(act)).sum())}}
    if reg["n_measured"] < MIN_TARGETS or not np.isfinite(act).all() or act.std() < 1e-12:
        why = f"{reg['n_measured']} measured targets < {MIN_TARGETS}" if reg["n_measured"] < MIN_TARGETS else "activity non-finite or constant"
        return {**out, "status": "not_estimable", "reason": why}
    out["status"] = "ok"
    out["activity_eligible"] = {"median": float(np.median(act)), "mean": float(act.mean()), "sd": float(act.std()), "median_total_counts": float(np.median(d["obs_el"]["total_counts"]))}
    L = lambda t: f"{spec}|{vname}|{t}"  # noqa: E731
    # (a) vs E2 scores
    a_blk = {}
    for ph in C7.PHASES:
        s = scores[ph]
        a_blk[ph] = C7._nonest("E2 phase not scorable") if s is None else C7.both_forms(s, act, sn, Z, L(f"a|{ph}"), B)
    out["a_e2"] = a_blk if not pc else C7._nonest("n above dense-eigh ceiling")
    # (b) cancer vs non-cancer
    ct = obs["cell_type"].astype(str).to_numpy()
    fin = np.isfinite(act_all)
    can, non = mb_in & (ct == "Cancer_Cells") & fin, mb_in & (ct != "Cancer_Cells") & fin
    if min(can.sum(), non.sum()) < MIN_B_GROUP:
        out["b_cancer_vs_noncancer"] = C7._nonest(f"groups {int(can.sum())}/{int(non.sum())} < {MIN_B_GROUP}")
    else:
        bb = {"status": "ok", "n_cancer": int(can.sum()), "n_noncancer": int(non.sum()), "median_cancer": float(np.median(act_all[can])), "median_noncancer": float(np.median(act_all[non])),
              **auc_mwu(act_all[can], act_all[non])}
        bb["median_diff"] = bb["median_cancer"] - bb["median_noncancer"]
        bb["by_noncancer_type"] = {t: {"n": int((non & (ct == t)).sum()), "median": float(np.median(act_all[non & (ct == t)]))} for t in sorted(set(ct[non])) if (non & (ct == t)).sum() >= MIN_B_GROUP}
        out["b_cancer_vs_noncancer"] = bb
    # main question and secondaries
    if pc:
        out["main_mdk"] = out["sender"] = out["secretory"] = C7._nonest("n above dense-eigh ceiling")
        return out
    out["secretory"] = C7.both_forms(sec, act, sn, Z, L("secretory"), B) if sec.std() > 0 else C7._nonest("constant secretory score")
    if not mdk_ok:
        out["main_mdk"] = out["sender"] = C7._nonest("MDK>0 in < 10 eligible compartments or constant")
        return out
    real = C7.both_forms(mdk_np, act, sn, Z, L("main_mdk"), B)
    out["main_mdk"] = {**real, "matched_gene_control": C7.control_block(real, act, ctrl_cols, Z) if ctrl_cols is not None else C7._nonest(ctrl_err)}
    out["sender"] = C7.both_forms(sender, act, sn, Z, L("sender"), B) if sender is not None else C7._nonest("sender score unavailable")
    return out


def analyze8(spec: str, nets: dict, commot_root, B: int = B, *, comp_run) -> dict:
    d = load8(spec, comp_run=comp_run)
    n, genes, X_el, obs_el = d["n"], d["genes"], d["X_el"], d["obs_el"]
    res = {"specimen": spec, "n_eligible": int(n), "n_b_candidates": int(d["mb_in"].sum())}
    if n < MIN_EL:
        res["status"] = "not_estimable"
        res["reason"] = f"{n} eligible < {MIN_EL}"
        return res
    res["status"] = "ok"
    pc = n > M.DENSE_EIGH_MAX_N
    mdk_np = R._column(X_el, genes, C.MDK)
    mdk_nz = int((mdk_np > 0).sum())
    res["mdk_nonzero"] = mdk_nz
    mdk_ok = mdk_nz >= MIN_NZ and mdk_np.std() > 0
    try:
        Z, zerr = H.covariate_matrix(obs_el, rna=RNA_COVARIATE), ""
    except H.CovariateDomainError as e:
        Z, zerr = None, str(e)
    res["covariate_error"] = zerr
    sc = A.load_score(commot_root / spec, "R00", list(obs_el.index))
    sender = None if isinstance(sc, dict) else np.asarray(sc, float)
    sn = None if pc else M.SpectralNull(knn_weights(d["xy"], K.K))
    sec = R._mean_z(X_el, genes, SECRETORY)[0]
    scores = C7.e2_scores(obs_el, genes, X_el) if not pc else {ph: None for ph in C7.PHASES}
    ctrl_cols, ctrl_err = None, ""
    try:
        mdk_i = genes.index(C.MDK)
        dec, pool = CP.prepare_pool([X_el], genes, X_el, mdk_i, [])
        excl = set(E2.signature_union()) | set().union(*[set(n_[n_["source"] == TF]["target"]) for n_ in nets.values()])
        pool = np.array([i for i in pool if genes[i] not in excl])
        dr = CP.make_draws(genes, dec, pool, mdk_i, [], spec + "|check8", n_ctrl=N_CTRL)
        ctrl_cols = X_el[:, dr["a"]]
        res["control_pool_size"] = int(len(pool))
    except ValueError as e:
        ctrl_err = str(e)
    res["variants"] = {}
    for vname, (_f, excl_v) in VARIANTS.items():
        net = nets[vname]
        reg = build_regulon(net, genes, excl_v)
        act_all = ulm_t(d["X"][:, reg["feat_idx"]], reg["w"]) if reg["n_measured"] >= MIN_TARGETS else np.full(len(d["X"]), np.nan)
        res["variants"][vname] = variant_block(spec, vname, reg, d, act_all, scores, sec, mdk_np, sender, sn, Z, ctrl_cols, ctrl_err, pc, mdk_ok, B)
        res["variants"][vname]["_panel_regulon"] = [{"target": t, "weight": float(net[(net.source == TF) & (net.target == t)]["weight"].iloc[0]), "on_panel_feature_set": bool(t in set(reg["feat_names"]))} for t in reg["targets"]]
    return res


# --------------------------------------------------------------------------- cohort
def cohort8(per_spec: dict, pairs) -> dict:
    out = {}
    for v in VARIANTS:
        for form in ("primary", "adjusted"):
            for ph in C7.PHASES:
                out[f"{v}|a_e2|{ph}|{form}"] = C7._family(per_spec, lambda r, ph=ph, form=form: C7._dig(r, "variants", v, "a_e2", ph, form))
            out[f"{v}|main_mdk|{form}"] = C7._family(per_spec, lambda r, form=form: C7._dig(r, "variants", v, "main_mdk", form))
            out[f"{v}|main_mdk|{form}"]["control"] = C7._ctrl_family(per_spec, lambda r, form=form: C7._dig(r, "variants", v, "main_mdk", "matched_gene_control", form))
            out[f"{v}|sender|{form}"] = C7._family(per_spec, lambda r, form=form: C7._dig(r, "variants", v, "sender", form))
            out[f"{v}|secretory|{form}"] = C7._family(per_spec, lambda r, form=form: C7._dig(r, "variants", v, "secretory", form))
        bs = {s: C7._dig(r, "variants", v, "b_cancer_vs_noncancer") for s, r in per_spec.items()}
        bs = {s: x for s, x in bs.items() if x and x.get("status") == "ok"}
        out[f"{v}|b_cancer_vs_noncancer"] = {"n_estimable": len(bs), "n_cancer_higher": int(sum(x["median_diff"] > 0 for x in bs.values())),
                                             "n_auc_gt_0.5": int(sum(x["auc"] > 0.5 for x in bs.values())),
                                             "median_auc": float(np.median([x["auc"] for x in bs.values()])) if bs else None, "median_diff_of_medians": float(np.median([x["median_diff"] for x in bs.values()])) if bs else None}
        out[f"{v}|c_biopsy_surgery"] = pair_table(per_spec, pairs, v)
    return out


def pair_table(per_spec: dict, pairs, v: str) -> dict:
    rows = []
    for bio, sur in pairs:
        gb, gs = C7._dig(per_spec.get(bio, {}), "variants", v, "activity_eligible"), C7._dig(per_spec.get(sur, {}), "variants", v, "activity_eligible")
        if not gb or not gs:
            rows.append({"biopsy": bio, "surgery": sur, "status": "not_estimable", "reason": "a member has < 30 eligible compartments or no activity"})
            continue
        rows.append({"biopsy": bio, "surgery": sur, "status": "ok", "n_eligible_biopsy": per_spec[bio]["n_eligible"], "n_eligible_surgery": per_spec[sur]["n_eligible"],
                     "median_activity_biopsy": gb["median"], "median_activity_surgery": gs["median"], "delta_surgery_minus_biopsy": gs["median"] - gb["median"],
                     "median_total_counts_biopsy": gb["median_total_counts"], "median_total_counts_surgery": gs["median_total_counts"]})
    ok = [r for r in rows if r["status"] == "ok"]
    return {"pairs": rows, "n_estimable_pairs": len(ok), "n_delta_lt0": int(sum(r["delta_surgery_minus_biopsy"] < 0 for r in ok)), "n_delta_gt0": int(sum(r["delta_surgery_minus_biopsy"] > 0 for r in ok)),
            "median_delta": float(np.median([r["delta_surgery_minus_biopsy"] for r in ok])) if ok else None}


def flat_rows8(per_spec: dict, genotype: dict) -> list:
    rows = []

    def add(spec, v, test, form, rec, ctrl=None):
        if not isinstance(rec, dict):
            return
        row = {"specimen": spec, "genotype": genotype.get(spec, ""), "n_eligible": per_spec[spec]["n_eligible"], "variant": v, "test": test, "form": form,
               "status": rec.get("status", ""), "rho": rec.get("rho", ""), "p_two_sided": rec.get("p_two_sided", ""), "q_bh": rec.get("q_bh", "")}
        if isinstance(ctrl, dict) and ctrl.get("status") == "ok":
            row.update({"ctrl_percentile": ctrl["percentile"], "ctrl_p95": ctrl["control_p95"], "real_gt_p95": ctrl["real_gt_p95"]})
        rows.append(row)

    for s, r in per_spec.items():
        for v, blk in r.get("variants", {}).items():
            if blk.get("status") != "ok":
                continue
            for form in ("primary", "adjusted"):
                for ph in C7.PHASES:
                    add(s, v, f"a_activity_vs_{ph}_e2", form, C7._dig(blk, "a_e2", ph, form))
                m = C7._dig(blk, "main_mdk")
                if m:
                    add(s, v, "main_mdk_np_vs_esr1_activity", form, m[form], m.get("matched_gene_control", {}).get(form))
                add(s, v, "sender_mdk_ncl_vs_esr1_activity", form, C7._dig(blk, "sender", form))
                add(s, v, "secretory_vs_esr1_activity", form, C7._dig(blk, "secretory", form))
    return rows


def gate_vs_check7(per_spec: dict, c7: dict) -> list:
    """Eligible counts and MDK>0 counts must equal check 7's (fail closed). ``c7`` = the check 7 per-specimen records computed in the same run."""
    g = []
    for s, r in per_spec.items():
        g.append(CM.gate_equal(f"{s} n_eligible vs check7", r["n_eligible"], c7[s]["n_eligible"], rtol=0, atol=0))
        pc = c7[s].get("partC", {})
        if r.get("status") == "ok" and pc.get("status") == "ok":
            g.append(CM.gate_equal(f"{s} mdk_nonzero vs check7", r["mdk_nonzero"], pc["mdk_nonzero"], rtol=0, atol=0))
    return g
