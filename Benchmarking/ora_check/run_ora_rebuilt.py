"""Corrected ORA Test B on rebuilt breast spots (rule: ora_rebuilt_DECISION_RULE.md). CPU only (htc).

  task N RUN    N = 0..11 breast (rebuilt inputs), 12..16 Xenium regions 0..4, 17 = legacy gate (breast index 0, OLD inputs, writes RUN/gate_legacy)
  gate RUN      compares the legacy gate with the stored sense calls and the audit's sens_calls.csv; writes RUN/gate_legacy/G1_PASS.json
  summary RUN   recomputes all mapping x call variants from the saved ORA tables, writes readout_rebuilt.{json,md}, variants.csv, repro.json

Breast is an adapted copy of run_ora_sense.breast_B (arms v41, v2 only; rebuilt input roots; the primary-library ORA table of every call kept).
Xenium is an adapted copy of run_ora_audit.xenium (oracle, v2 cross-fit, v41; real and null; no full-fit arm).
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE.parent))
from ora_check import ora_corrected_mapping as M  # noqa: E402
from ora_check import ora_lib as L  # noqa: E402
from ora_check import ora_rebuilt_lib as RL  # noqa: E402
from ora_check import ora_sense_lib as S  # noqa: E402
from ora_check import run_ora_sense as RS  # noqa: E402

V2_TS = "20261102T200000Z"
REBUILT_M3 = REPO / "output" / f"module3_rna_refined_7type_{V2_TS}" / "v2b1"
REBUILT_NUC = REPO / "output" / f"nuclei_rebuilt_7type_{V2_TS}" / "v2b2"
SENSE_RUN = HERE / "results_sense_20261101_1550"
AUDIT_RUN = HERE / "results_audit_20261101_1624"
N_BREAST, N_XEN, LEGACY_TASK = 12, 5, 17
G2_TOL = 1e-4
PRIMARY = RS.PRIMARY


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fail(run, tag, info):
    p = Path(run) / "FAIL" / f"{tag}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(info, indent=2, default=str))
    print("STOP", tag, info, flush=True)
    sys.exit(3)


def write_new(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "x") as f:  # refuses to overwrite
        f.write(text)


# ------------------------------------------------------------------ recording of the primary-library ORA table
class Rec:
    def __init__(self):
        self.dfs = []
        self.orig = L.ora
        L.ora = self._ora

    def _ora(self, genes, gmt, background):
        df = self.orig(genes, gmt, background)
        self.dfs.append(df)
        return df

    def row(self, genes, background, label, mapping):
        """RS.ora_row unchanged; returns (row, primary-library table)."""
        self.dfs = []
        row = RS.ora_row(genes, background, label, mapping)
        names = list(RS.libs().keys())
        tb = self.dfs[names.index(PRIMARY)][["Term", "Adjusted P-value"]].copy()
        self.dfs = []
        return row, tb


def write_tables(path: Path, tables: dict):
    parts = []
    for k, tb in tables.items():
        t = tb.copy()
        t.insert(0, "table_key", k)
        parts.append(t)
    big = pd.concat(parts, ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    big.to_parquet(path, index=False)


def read_tables(path: Path) -> dict:
    big = pd.read_parquet(path)
    return {k: g[["Term", "Adjusted P-value"]].reset_index(drop=True) for k, g in big.groupby("table_key", sort=False)}


# ------------------------------------------------------------------ breast
def breast_inputs(sample, mode):
    if mode == "legacy":
        return dict(
            p0=RS.PLAIN_DIR / f"{sample}_cell_prop_global_results.csv",
            v2=RS.V2_PROPS / sample / "cell_prop.csv",
            nuc=RS.NUCLEI_DIR / sample / f"{sample}_nuclei_counts_fullres_patchwise.csv",
        )
    return dict(
        p0=REBUILT_M3 / "p0" / sample / f"{sample}_cell_prop_global_results.csv",
        v2=REBUILT_M3 / f"{sample}_cell_prop_finetuned_results.csv",
        nuc=REBUILT_NUC / sample / f"{sample}_nuclei_counts_fullres_patchwise.csv",
    )


def breast(idx, run, mode):
    sys.path.insert(0, str(REPO / "Benchmarking" / "brightness_estimator"))
    import breast_predict as bp
    from CITEgeist.model.deconvolution.rna_refinement_adata import run_rna_refinement
    from ora_check.run_sace_layers import load_v41_profile

    run = Path(run)
    sample = RS.SPECIMENS[idx]
    out = run / ("gate_legacy" if mode == "legacy" else "")
    if mode == "rebuilt" and not (run / "gate_legacy" / "G1_PASS.json").exists():
        fail(run, f"breast_{idx:02d}_no_G1", {"error": "gate_legacy/G1_PASS.json missing: run the legacy gate first"})
    seed = RS.BASE_SEED + idx
    mapping = RS.load_mapping("term_mapping_7type.yaml")
    markers = L.marker_gene_set(load_v41_profile())
    inp = breast_inputs(sample, mode)
    for k, p in inp.items():
        assert p.is_file(), f"missing input {k}: {p}"
    rec = Rec()
    p0_full = pd.read_csv(inp["p0"], index_col=0)
    nuc = pd.read_csv(inp["nuc"], index_col=0).iloc[:, 0]
    stored_v2 = pd.read_csv(inp["v2"], index_col=0)
    keep = p0_full.index[p0_full.index.isin(nuc.index)]
    if mode == "rebuilt":  # G4
        ref = pd.read_csv(RS.PLAIN_DIR / f"{sample}_cell_prop_global_results.csv", index_col=0)
        common = p0_full.index.intersection(ref.index)
        g4 = {
            "same_spots": bool(set(stored_v2.index) == set(p0_full.index) == set(nuc.index) == set(ref.index)),
            "p0_vs_marker_share_max_abs": float((p0_full.loc[common] - ref.loc[common]).abs().max().max()),
            "n_common_with_marker_share": int(len(common)),
            "rows_sum_to_1_v2": bool(np.allclose(stored_v2.sum(axis=1), 1.0, atol=1e-6)),
            "rows_sum_to_1_p0": bool(np.allclose(p0_full.sum(axis=1), 1.0, atol=1e-6)),
            "columns_equal": bool(list(stored_v2.columns) == list(p0_full.columns)),
        }
        if not (g4["same_spots"] and g4["p0_vs_marker_share_max_abs"] == 0.0 and g4["rows_sum_to_1_v2"] and g4["rows_sum_to_1_p0"] and g4["columns_equal"]):
            fail(run, f"breast_{idx:02d}_G4", g4)
    else:
        g4 = {}
        assert set(stored_v2.index) == set(keep)
    spots = list(stored_v2.index)
    uni_file = RS.ARM_DIRS["v41"] / "sace_layers" / sample / "layers_real_marker.parquet"
    layer_genes = list(pd.read_parquet(uni_file).columns)
    gex_full = bp.load_raw_gex(sample, list(p0_full.index))
    gex_full.layers["raw_counts"] = gex_full.X.copy()
    gex = gex_full[list(keep)].copy()
    p0 = p0_full.loc[list(keep)]
    types = list(p0.columns)
    gn = list(gex.var_names)
    raw = gex.layers["raw_counts"]
    Y = np.asarray(raw.toarray() if hasattr(raw, "toarray") else raw, float)
    pos = {g: i for i, g in enumerate(gn)}
    obs_pos = {b: i for i, b in enumerate(gex.obs_names)}
    Yc = Y[np.array([obs_pos[b] for b in spots])]
    lib = Yc.sum(axis=1)
    ok = lib > 0
    uni_genes = [g for g in layer_genes if g in pos]
    mean = Yc[ok].mean(axis=0)
    expressed = [g for g in uni_genes if mean[pos[g]] >= S.MIN_MEAN]
    cand = [g for g in expressed if g not in markers]
    bg = [g.upper() for g in expressed]
    gi = np.array([pos[g] for g in uni_genes])
    Xl = np.log1p(1e4 * Yc[ok][:, gi] / lib[ok][:, None])
    ex, rl = S.prep_gene_ranks(Xl, lib[ok])
    folds = S.make_folds(len(gn), seed)
    gene_pos_in_uni = np.full(len(gn), -1)
    for i, g in enumerate(uni_genes):
        gene_pos_in_uni[pos[g]] = i
    meta = {
        "sample": sample, "mode": mode, "n_spots": len(spots), "n_spots_p0": int(len(p0_full)), "n_genes_universe": len(uni_genes),
        "n_expressed": len(expressed), "n_candidates": len(cand), "n_zero_nuclei_spots": int((nuc.loc[list(keep)] == 0).sum()), "G4": g4,
    }

    def refit(keepmask):
        g2 = gex[:, np.flatnonzero(keepmask)].copy()
        p, _ = run_rna_refinement(g2, p0, nuc.loc[list(keep)], n_hvg=2000, c=0.1, n_passes=3)
        return p.loc[spots].to_numpy(float)

    full = refit(np.ones(len(gn), dtype=bool))
    meta["repro_max_abs_v2"] = float(np.abs(full - stored_v2.loc[spots][types].to_numpy(float)).max())
    if mode == "rebuilt" and meta["repro_max_abs_v2"] > G2_TOL:
        fail(run, f"breast_{idx:02d}_G2", meta)
    shares = {"v41": [p0.loc[spots].to_numpy(float)] * len(folds), "v2": []}
    for f in folds:
        km = np.ones(len(gn), dtype=bool)
        km[f] = False
        shares["v2"].append(refit(km))
    print("refits done", meta, flush=True)

    def rho_for(arm, transform):
        o = np.full((len(uni_genes), len(types)), np.nan)
        for k, f in enumerate(folds):
            sh = transform(shares[arm][k])[ok]
            fu = gene_pos_in_uni[f]
            fu = fu[fu >= 0]
            o[fu] = S.partial_rho(ex[:, fu], rl, sh)
        return pd.DataFrame(o, index=uni_genes, columns=types)

    rows, tables = [], {}
    for arm in ("v41", "v2"):
        evaluable = {t: all(S.evaluable_share(sh[:, j]) for sh in shares[arm]) for j, t in enumerate(types)}
        conds = {"real": lambda x: x, "null": lambda x: S.permute_fraction(x, 100, seed, null_seed=seed)}
        for w in (25, 50, 75):
            conds[f"w{w}"] = (lambda ww: (lambda x: S.permute_fraction(x, ww, seed + 1000 * ww)))(w)
        for cond, tf in conds.items():
            rho = rho_for(arm, tf)
            genes = RS.gene_lists_from_scores(rho, cand)
            for t in types:
                row, tb = rec.row(genes[t], bg, t, mapping)
                tk = f"breast|{arm}|{cond}|{sample}|{t}"
                tables[tk] = tb
                rows.append(dict(test="B", dataset="breast", arm=arm, condition=cond, unit=sample, type=t, evaluable=evaluable[t], table_key=tk, **row))
                if cond == "null":
                    rows.append({**rows[-1], "condition": "w100"})  # w = 100 is the full permutation (rule); same table
        print("B done", sample, arm, flush=True)
    df = pd.DataFrame(rows)
    write_new(out / "calls" / f"B_breast_{idx:02d}.csv", df.to_csv(index=False))
    write_new(out / "calls" / f"B_breast_{idx:02d}_meta.json", json.dumps(meta, indent=2))
    write_tables(out / "tables" / f"B_breast_{idx:02d}.parquet", tables)
    inputs = {k: {"path": str(v), "sha256": sha(v)} for k, v in inp.items()}
    inputs["universe_parquet"] = {"path": str(uni_file), "sha256": sha(uni_file)}
    write_new(out / "inputs" / f"breast_{idx:02d}.json", json.dumps(inputs, indent=2))
    if mode == "legacy":
        chk = check_vs_stored(df, SENSE_RUN / "calls" / f"B_breast_{idx:02d}.csv", "legacy breast vs stored Test B")
        write_new(out / "G1_breast_vs_stored.json", json.dumps(chk, indent=2))
        if chk["n_mismatch"] or not chk["compared"]:
            fail(run, f"legacy_breast_{idx:02d}", chk)


def check_vs_stored(mine: pd.DataFrame, orig_csv, label):
    orig = pd.read_csv(orig_csv, keep_default_na=False, na_values=[""])
    mine = mine.copy()
    mine["condition"] = mine["condition"].astype(str)
    orig["condition"] = orig["condition"].astype(str)
    key = ["arm", "condition", "unit", "type"]
    m = mine.merge(orig, on=key, suffixes=("", "_orig"))
    bad = []
    for _, r in m.iterrows():
        for col in [c for c in orig.columns if c.startswith("pass_")] + ["passed"]:
            if bool(r[col]) != bool(r[col + "_orig"]):
                bad.append((r["arm"], r["condition"], r["unit"], r["type"], col))
        a, b = float(r["nlp"]), float(r["nlp_orig"])
        if abs(a - b) > 1e-6 * max(1.0, abs(b)):
            bad.append((r["arm"], r["condition"], r["unit"], r["type"], "nlp", a, b))
    return {"label": label, "compared": int(len(m)), "n_mismatch": len(bad), "mismatches": bad[:20]}


# ------------------------------------------------------------------ Xenium
def xenium(r, run):
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    run = Path(run)
    if not (run / "gate_legacy").exists():
        fail(run, f"xenium_{r}_no_run", {"error": "run dir not initialised by the legacy task"})
    rec = Rec()
    mapping = RS.load_mapping("term_mapping_xenium7.yaml")
    markers = RS._xen_markers()
    seed = RS.XEN_SEED + r
    _, data = data_io.load_xenium(7)
    d = data[r]
    gt = RS._gt(r)
    obs = list(d["obs"])
    genes = [g.upper() for g in d["genes"]]
    p0 = np.asarray(rs.ms_share(d["X"], d["names"], rs.PROF[7]))
    Y = d["Y"]
    lib_all = Y.sum(axis=1)
    sel = np.array([i for i, b in enumerate(obs) if b in gt.index and lib_all[i] > 0])
    mean = Y[sel].mean(axis=0)
    expressed = [g for g, m in zip(genes, mean) if m >= S.MIN_MEAN]
    cand = [g for g in expressed if g not in markers]
    tb_frozen = RS._testable(genes, mapping)
    ex, rl = S.prep_gene_ranks(np.log1p(1e4 * Y[sel] / lib_all[sel][:, None]), lib_all[sel])
    folds = S.make_folds(len(genes), seed)
    shares = {"v41": [p0] * 2, "oracle": [gt.reindex(obs).to_numpy(float)] * 2, "v2": []}
    for f in folds:
        mask = np.ones(len(genes), dtype=bool)
        mask[f] = False
        p, _ = refine_proportions(Y, d["nc"], p0, c=0.1, n_passes=3, model="multinomial", gene_mask=mask)
        shares["v2"].append(p)

    def rho_for(arm, tf):
        o = np.full((len(genes), len(RS.XTYPES)), np.nan)
        for k, f in enumerate(folds):
            o[f] = S.partial_rho(ex[:, f], rl, tf(shares[arm][k][sel]))
        return pd.DataFrame(o, index=genes, columns=RS.XTYPES)

    rows, tables = [], {}
    for arm in ("oracle", "v2", "v41"):
        for cond, tf in {"real": lambda x: x, "null": lambda x: S.permute_fraction(x, 100, seed, null_seed=seed)}.items():
            gl = RS.gene_lists_from_scores(rho_for(arm, tf), cand)
            for j, t in enumerate(RS.XTYPES):
                row, tb = rec.row(gl[t], expressed, t, mapping)  # ORA for every type; testability is applied per mapping later
                tk = f"xenium|{arm}|{cond}|region_{r}|{t}"
                tables[tk] = tb
                ev = bool(all(S.evaluable_share(sh[sel][:, j]) for sh in shares[arm]))
                rows.append(dict(test="B", dataset="xenium", arm=arm, condition=cond, unit=f"region_{r}", type=t, evaluable=ev,
                                 testable=bool(tb_frozen[t] >= S.MIN_PANEL_GENES), n_panel_genes=tb_frozen[t], table_key=tk, **row))
    df = pd.DataFrame(rows)
    # reproduction against the stored Test B rows (frozen mapping; every type is testable there, so the stored rows have ORA for all)
    stored = pd.read_csv(SENSE_RUN / "calls" / f"B_xenium_{r}.csv", keep_default_na=False, na_values=[""])
    chk = check_vs_stored(df, SENSE_RUN / "calls" / f"B_xenium_{r}.csv", f"xenium region {r} vs stored Test B")
    chk["stored_all_testable"] = bool(stored["testable"].astype(bool).all())
    write_new(run / "calls" / f"B_xenium_{r}.csv", df.to_csv(index=False))
    write_new(run / "calls" / f"B_xenium_{r}_meta.json", json.dumps({"region": r, "n_spots_selected": int(len(sel)), "n_expressed": len(expressed),
                                                                    "n_candidates": len(cand), "panel_genes": genes, "repro": chk}, indent=2))
    write_tables(run / "tables" / f"B_xenium_{r}.parquet", tables)
    if chk["n_mismatch"] or not chk["compared"] or not chk["stored_all_testable"]:
        fail(run, f"xenium_{r}_repro", chk)


def _b(x):
    return str(x).strip().lower() == "true"


# ------------------------------------------------------------------ legacy gate (G1) comparison with the audit's sens_calls.csv
def gate(run):
    run = Path(run)
    g = run / "gate_legacy"
    calls = pd.read_csv(g / "calls" / "B_breast_00.csv", keep_default_na=False, na_values=[""])
    calls["condition"] = calls["condition"].astype(str)
    tables = read_tables(g / "tables" / "B_breast_00.parquet")
    mapping = RS.load_mapping("term_mapping_7type.yaml")
    lin = RL.LineageCache(L.term_lineage, mapping)
    sens = pd.read_csv(AUDIT_RUN / "sens_calls.csv", keep_default_na=False, na_values=[""])
    sens["condition"] = sens["condition"].astype(str)
    sens = sens[(sens["dataset"] == "breast") & (sens["unit"] == RS.SPECIMENS[0])]
    bad, n = [], 0
    for r in calls.itertuples():
        if r.condition not in ("real", "null") or r.arm not in ("v41", "v2"):
            continue
        s = sens[(sens["arm"] == r.arm) & (sens["condition"] == r.condition) & (sens["type"] == r.type)]
        if len(s) != 1:
            bad.append(("no audit row", r.arm, r.condition, r.type))
            continue
        tb = tables[r.table_key]
        target = mapping["type_to_lineage"][r.type]
        cur = RL.current_call(tb, target, lin)
        col = RL.collapsed_pass(tb, target, lin)
        n += 1
        for name, mine, theirs in (("current", cur["passed"], _b(s["current"].iloc[0])), ("S_collapsed", col, _b(s["S_collapsed"].iloc[0])),
                                   ("nlp", cur["nlp"], float(r.nlp)), ("passed", cur["passed"], bool(r.passed))):
            ok = (abs(mine - theirs) <= 1e-9 * max(1.0, abs(theirs))) if name == "nlp" else (bool(mine) == bool(theirs))
            if not ok:
                bad.append((name, r.arm, r.condition, r.type, mine, theirs))
    stored = json.loads((g / "G1_breast_vs_stored.json").read_text())
    info = {"audit_rows_compared": n, "audit_mismatches": bad[:20], "n_audit_mismatch": len(bad), "breast_vs_stored": stored}
    if n == 0 or bad or stored["n_mismatch"] or not stored["compared"]:
        fail(run, "G1", info)
    write_new(g / "G1_PASS.json", json.dumps(info, indent=2))
    print("G1 PASS", info, flush=True)


# ------------------------------------------------------------------ summary
MAPPINGS = {"frozen": 0, "L1": 1, "L2": 2, "L3": 3}
VARIANTS = ("current", "collapsed")
OLD_B_KEY = "B"


def summary(run):
    run = Path(run)
    ov = M.load_overrides()
    mp_b, mp_x = RS.load_mapping("term_mapping_7type.yaml"), RS.load_mapping("term_mapping_xenium7.yaml")
    bcalls = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in sorted((run / "calls").glob("B_breast_[0-9][0-9].csv"))], ignore_index=True)
    xcalls = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in sorted((run / "calls").glob("B_xenium_[0-9].csv"))], ignore_index=True)
    for df in (bcalls, xcalls):
        df["condition"] = df["condition"].astype(str)
        for c in ("evaluable", "passed"):
            df[c] = df[c].astype(bool)
    if bcalls["unit"].nunique() != S.N_SPECIMENS or xcalls["unit"].nunique() != S.N_REGIONS:
        fail(run, "summary_units", {"breast_units": int(bcalls["unit"].nunique()), "xenium_units": int(xcalls["unit"].nunique())})
    tables = {}
    for f in sorted((run / "tables").glob("B_*.parquet")):
        tables.update(read_tables(f))
    panel = {f"region_{r}": json.loads((run / "calls" / f"B_xenium_{r}_meta.json").read_text())["panel_genes"] for r in range(S.N_REGIONS)}
    gmt = RS.libs()[PRIMARY]
    repro = {"G2_max_abs_v2": {}, "G3": {}, "G4": {}}
    for f in sorted((run / "calls").glob("B_breast_[0-9][0-9]_meta.json")):
        m = json.loads(f.read_text())
        repro["G2_max_abs_v2"][m["sample"]] = m["repro_max_abs_v2"]
        repro["G4"][m["sample"]] = m["G4"]
    repro["G1"] = json.loads((run / "gate_legacy" / "G1_PASS.json").read_text())
    repro["xenium_repro"] = {r: json.loads((run / "calls" / f"B_xenium_{r}_meta.json").read_text())["repro"] for r in range(S.N_REGIONS)}

    old = json.loads((SENSE_RUN / "readout_sense.json").read_text())["B"]
    res, allv = {}, []
    for mname, lv in MAPPINGS.items():
        fn = M.corrected_lineage_fn(ov, lv)
        lb, lx = RL.LineageCache(fn, mp_b), RL.LineageCache(fn, mp_x)
        tcount = {u: RL.testable_counts([g.upper() for g in panel[u]], gmt, mp_x, lx) for u in panel}
        for var in VARIANTS:
            b = RL.variant_table(bcalls, tables, mp_b, lb, var)
            x = RL.variant_table(xcalls, tables, mp_x, lx, var, testable=tcount)
            if mname == "frozen" and var == "current":  # G3: the table path reproduces the in-job calls exactly
                dp = int((b["passed"].to_numpy() != bcalls["passed"].to_numpy()).sum() + (x["passed"].to_numpy() != xcalls["passed"].to_numpy()).sum())
                dn = float(max(np.abs(b["nlp"].to_numpy(float) - bcalls["nlp"].to_numpy(float)).max(), np.abs(x["nlp"].to_numpy(float) - xcalls["nlp"].to_numpy(float)).max()))
                repro["G3"] = {"n_pass_differences": dp, "max_abs_nlp_difference": dn}
                if dp or dn > 1e-9:
                    fail(run, "G3", repro["G3"])
            r = S.readout(b, x)
            res[f"{mname}|{var}"] = {"readout": r, "xenium_testable_counts": tcount, "xenium_pass_totals": x[x["condition"] == "real"].groupby("arm")["passed"].sum().to_dict(),
                                     "xenium_null_oracle": int(x[(x["arm"] == "oracle") & (x["condition"] == "null")]["passed"].sum()),
                                     "breast_null_totals": {a: int(b[(b["arm"] == a) & (b["condition"] == "null")]["passed"].sum()) for a in ("v41", "v2")},
                                     "breast_pass_by_type": {a: b[(b["arm"] == a) & (b["condition"] == "real")].groupby("type")["passed"].sum().to_dict() for a in ("v41", "v2")}}
            keep = ["dataset", "arm", "condition", "unit", "type", "evaluable", "passed", "best_rank", "nlp"]
            allv.append(pd.concat([b[keep], x[keep]]).assign(mapping=mname, variant=var))
    pd.concat(allv, ignore_index=True).to_csv(run / "variants.csv", index=False)
    flags = changed_flags(res, old)
    comp_libs = {a: {c: int(bcalls[(bcalls["arm"] == a) & (bcalls["condition"] == "real")][c].sum()) for c in bcalls.columns if c.startswith("pass_")} for a in ("v41", "v2")}
    out = {"results": res, "changed_flags": flags, "old_frozen_B": {k: old[k] for k in ("verdict", "total", "per_type", "loss", "delta", "V1", "V2", "V3")},
           "comparison_library_passes_frozen_mapping_real": comp_libs}
    (run / "readout_rebuilt.json").write_text(json.dumps(out, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    (run / "repro.json").write_text(json.dumps(repro, indent=2, default=str))
    (run / "readout_rebuilt.md").write_text(render_md(out, repro))
    print(json.dumps({k: v["readout"]["verdict"] for k, v in res.items()}), flush=True)
    print(json.dumps(flags, indent=1), flush=True)


def changed_flags(res, old):
    """Rule section 'what counts as a changed conclusion': (A) rebuild, (B) mapping, (C) per-type sign, (D) calibration."""
    f = {}
    base = res["frozen|current"]["readout"]
    f["A_rebuild_verdict"] = {"baseline": old["verdict"], "rebuilt": base["verdict"], "flag": "SAME" if base["verdict"] == old["verdict"] else "CHANGED"}
    types = sorted(base["per_type"]["v2"])

    def signs(r):
        return {t: S_sign(r["per_type"]["v2"][t] - r["per_type"]["v41"][t]) for t in types}

    old_signs = {t: S_sign(old["per_type"]["v2"][t] - old["per_type"]["v41"][t]) for t in types}
    f["C_rebuild_vs_old_signs"] = {"old": old_signs, "rebuilt": signs(base), "flag": "SAME" if old_signs == signs(base) else "CHANGED"}
    for var in VARIANTS:
        ref = res[f"frozen|{var}"]["readout"]
        for m in ("L1", "L2", "L3"):
            r = res[f"{m}|{var}"]["readout"]
            v3r, v3c = ref["V3"], r["V3"]
            f[f"{m}|{var}"] = {
                "B_verdict": {"frozen": ref["verdict"], "corrected": r["verdict"], "flag": "SAME" if ref["verdict"] == r["verdict"] else "CHANGED"},
                "C_signs": {"frozen": signs(ref), "corrected": signs(r), "flag": "SAME" if signs(ref) == signs(r) else "CHANGED"},
                "D_calibration": {k: {"frozen": v3r.get(k), "corrected": v3c.get(k)} for k in ("status", "oracle_ge_methods", "gap_ok")},
            }
            f[f"{m}|{var}"]["D_calibration"]["flag"] = "SAME" if all(v["frozen"] == v["corrected"] for v in f[f"{m}|{var}"]["D_calibration"].values() if isinstance(v, dict)) else "CHANGED"
    return f


def S_sign(x):
    return int(np.sign(x))


def render_md(out, repro):
    res, old = out["results"], out["old_frozen_B"]
    L_ = ["# Corrected ORA Test B on rebuilt breast spots (primary library CellMarker_2024) -- post-hoc mapping sensitivity", "",
          "Rule: `ora_rebuilt_DECISION_RULE.md`. Plausibility check, not an accuracy measure. The author decides; nothing here chooses a mapping.", "",
          f"Old frozen Test B (stored, `results_sense_20261101_1550`): verdict {old['verdict']}, totals {json.dumps(old['total'])}, delta {old['delta']}, "
          f"loss v41-v2 {json.dumps(old['loss'])}.", "", "## Totals of 84, side by side", "",
          "| mapping | call | v41 | v2 | delta | verdict | V1 gap v41/v2 | V2 ok v41/v2 | V3 | null v41/v2 | Xenium oracle/v2/v41 (of 35) oracle null |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, v in res.items():
        r = v["readout"]
        m, c = k.split("|")
        xp = v["xenium_pass_totals"]
        L_.append(f"| {m} | {c} | {r['total']['v41']} | {r['total']['v2']} | {r['delta']} | {r['verdict']} | {r['V1']['v41']['gap']}/{r['V1']['v2']['gap']} | "
                  f"{r['V2']['v41']['ok']}/{r['V2']['v2']['ok']} | {r['V3'].get('status')} | {v['breast_null_totals']['v41']}/{v['breast_null_totals']['v2']} | "
                  f"{xp.get('oracle')}/{xp.get('v2')}/{xp.get('v41')} null {v['xenium_null_oracle']} |")
    L_ += ["", "## Passes of 12 by type (v41 / v2), real shares", ""]
    types = sorted(res["frozen|current"]["readout"]["per_type"]["v2"])
    L_ += ["| type | " + " | ".join(res.keys()) + " |", "|---|" + "---|" * len(res)]
    for t in types:
        L_.append(f"| {t} | " + " | ".join(f"{v['readout']['per_type']['v41'][t]} / {v['readout']['per_type']['v2'][t]}" for v in res.values()) + " |")
    L_ += ["", "## Changed-conclusion flags (rule)", "", "```", json.dumps(out["changed_flags"], indent=1, default=str), "```", "",
           "## Reproduction (G1-G4, Xenium)", "", "```", json.dumps({k: (v if k != "G1" else {"n_audit_mismatch": v["n_audit_mismatch"], "audit_rows_compared": v["audit_rows_compared"]}) for k, v in repro.items()}, indent=1, default=str)[:6000], "```", ""]
    return "\n".join(L_) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["task", "gate", "summary"])
    ap.add_argument("a")
    ap.add_argument("b", nargs="?")
    a = ap.parse_args()
    if a.cmd == "task":
        n, run = int(a.a), Path(a.b)
        if n == LEGACY_TASK:
            run.mkdir(parents=True, exist_ok=False)  # refuses an existing run dir
            (run / "gate_legacy").mkdir()
            breast(0, run, "legacy")
        elif n < N_BREAST:
            breast(n, run, "rebuilt")
        elif n < N_BREAST + N_XEN:
            xenium(n - N_BREAST, run)
        else:
            raise SystemExit("task id out of range")
    elif a.cmd == "gate":
        gate(Path(a.a))
    else:
        summary(Path(a.a))


if __name__ == "__main__":
    main()
