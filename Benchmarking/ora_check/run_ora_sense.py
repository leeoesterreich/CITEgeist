"""Proportion-sensitive ORA sense check (rule: ora_sense_DECISION_RULE.md). CPU only (htc).

Subcommands (all write rows to <run>/calls/<name>.csv):
  breast-A  IDX RUN   Test A (stored + blend SACE layers)      breast-B  IDX RUN   Test B (cross-fit partial Spearman)
  xenium-sace R RUN   fresh no-round SACE for V3 (Test A)      xenium-A  R RUN     Test A on those layers
  xenium-B  R RUN     Test B on Xenium
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE.parent))
from ora_check import ora_lib as L  # noqa: E402
from ora_check import ora_sense_lib as S  # noqa: E402
from ora_check import run_ora as R  # noqa: E402

V2_RUN = HERE / "results_v2_20261101_1501"
ARM_DIRS = {"v41": V2_RUN / "v41", "v2": V2_RUN / "v2", "v2_noparam": V2_RUN / "arm_c" / "v2_noparam"}
BASE_SEED = 20261014
XEN_SEED = 20261114
PRIMARY = "CellMarker_2024"
SPECIMENS = R.PROTO["specimens"]
PLAIN_DIR = REPO / "output" / "module3_marker_share_7type"
NUCLEI_DIR = REPO / "output" / "morphology_assignment_v3"
V2_PROPS = REPO / "Benchmarking/midkine_remediation/results_engine_v2/r1_20261031_1442/a_V3_multinomial"
NOPARAM_PROPS = REPO / "Benchmarking/midkine_remediation/results_engine_v2/noparam_20261101_1455/v2_noparam"
XB = REPO / "Benchmarking" / "xenium_benchmarking"
GATE_RUN = XB / "rna_refinement" / "runs" / "gexgate_20261031_1832"
P138_PROPS = XB / "CITEgeist" / "output_marker_share_p138_7type"
GT_PROP = XB / "ground_truth_xenium_rnaprot" / "ground_truth_7type" / "ground_truth"
CD4, CD8 = "CD4+ T cells", "CD8+ T cells"
XTYPES = ["B cells", CD4, CD8, "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]
XFILE = {"B cells": "B_cells", CD4: "CD4pos_T_cells", CD8: "CD8pos_T_cells", "Macrophages": "Macrophages",
         "Endothelial": "Endothelial", "Epithelial": "Epithelial", "Fibroblasts": "Fibroblasts"}
LIBS = None


def libs():
    global LIBS
    if LIBS is None:
        proto = R.deep_merge(R.PROTO, json.loads((HERE / "ora_protocol_v41_7type.json").read_text()))
        LIBS = {}
        for name, spec in proto["ora"]["libraries"].items():
            f = HERE / spec["file"]
            assert R.sha(f) == spec["sha256"], f"library hash mismatch: {name}"
            LIBS[name] = R.read_gmt(f)
    return LIBS


def load_mapping(name):
    return yaml.safe_load((HERE / name).read_text())


def ora_row(genes, background, label, mapping):
    """ORA in every library; primary-library values in passed/best_rank/nlp, the rest in suffixed columns."""
    row = {}
    for lib, gmt in libs().items():
        df = L.ora([g.upper() for g in genes], gmt, background)
        sc = S.score_call(df, mapping["type_to_lineage"][label], mapping, L.term_lineage)
        row[f"pass_{lib}"], row[f"rank_{lib}"], row[f"nlp_{lib}"] = sc["passed"], sc["best_rank"], sc["nlp"]
    row["passed"], row["best_rank"], row["nlp"] = row[f"pass_{PRIMARY}"], row[f"rank_{PRIMARY}"], row[f"nlp_{PRIMARY}"]
    return row


def emit(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    print("wrote", path, len(rows), "rows", flush=True)


def gene_lists_from_scores(score_df, candidates, n=S.TOP_N):
    """score_df: genes x types DataFrame -> {type: top-n genes}."""
    return {t: S.top_by_score(score_df[t], candidates, n) for t in score_df.columns}


# ============================================================ breast Test A
def breast_A(idx, run):
    from ora_check.run_sace_layers import load_v41_profile

    sample = SPECIMENS[idx]
    mapping = load_mapping("term_mapping_7type.yaml")
    markers = L.marker_gene_set(load_v41_profile())
    rows = []
    for arm, adir in ARM_DIRS.items():
        d = adir / "sace_layers" / sample
        real = pd.read_parquet(d / "layers_real_marker.parquet")
        uni = pd.read_parquet(d / "layers_uniform_marker.parquet")
        props = pd.read_csv(d / "props.csv", index_col=0)
        expressed = real.sum(axis=0) >= S.MIN_MEAN
        cand = [g for g in real.columns[expressed.to_numpy()] if g not in markers]
        bg = [g.upper() for g in real.columns[expressed.to_numpy()]]
        conds = {"real": d / "layers_real_marker.parquet", "null": d / "layers_spot_permuted_marker.parquet",
                 "w100": d / "layers_uniform_marker.parquet"}
        if arm in ("v41", "v2"):
            for w in (25, 50, 75):
                conds[f"w{w}"] = Path(run) / "blend" / arm / "sace_layers" / sample / f"layers_blend{w}_marker.parquet"
        for cond, f in conds.items():
            lay = pd.read_parquet(f).reindex(columns=real.columns)
            dl = S.delta_table(lay, uni)
            genes = gene_lists_from_scores(dl.T, cand)
            for t in real.index:
                rows.append(dict(test="A", dataset="breast", arm=arm, condition=cond, unit=sample, type=t,
                                 evaluable=S.evaluable_share(props[t]), **ora_row(genes[t], bg, t, mapping)))
        print("A done", sample, arm, flush=True)
    emit(Path(run) / "calls" / f"A_breast_{idx:02d}.csv", rows)


# ============================================================ breast Test B
def _p0_nuclei(sample):
    p0_full = pd.read_csv(PLAIN_DIR / f"{sample}_cell_prop_global_results.csv", index_col=0)
    nuc = pd.read_csv(NUCLEI_DIR / sample / f"{sample}_nuclei_counts_fullres_patchwise.csv", index_col=0).iloc[:, 0]
    keep = p0_full.index[p0_full.index.isin(nuc.index)]
    return p0_full, nuc, keep


def breast_B(idx, run):
    sys.path.insert(0, str(REPO / "Benchmarking" / "brightness_estimator"))
    import breast_predict as bp
    from CITEgeist.model.deconvolution.rna_refinement_adata import run_rna_refinement
    from ora_check.run_sace_layers import load_v41_profile

    sample = SPECIMENS[idx]
    seed = BASE_SEED + idx
    mapping = load_mapping("term_mapping_7type.yaml")
    markers = L.marker_gene_set(load_v41_profile())
    p0_full, nuc, keep = _p0_nuclei(sample)
    stored = {arm: pd.read_csv(p / sample / "cell_prop.csv", index_col=0) for arm, p in (("v2", V2_PROPS), ("v2_noparam", NOPARAM_PROPS))}
    for arm, st in stored.items():
        assert set(st.index) == set(keep), f"{arm}: spot set differs from p0 & nuclei"
    spots = list(stored["v2"].index)  # common spot set (the v2 set); same order for every arm
    layer_genes = list(pd.read_parquet(ARM_DIRS["v41"] / "sace_layers" / sample / "layers_real_marker.parquet").columns)
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
    meta = {"sample": sample, "n_spots": len(spots), "n_genes_universe": len(uni_genes), "n_expressed": len(expressed),
            "n_candidates": len(cand)}

    def refit(arm, keepmask):
        kw = dict(c=0.1, n_passes=3) if arm == "v2" else dict(c=0.0, n_passes=1)
        g2 = gex[:, np.flatnonzero(keepmask)].copy()
        p, _ = run_rna_refinement(g2, p0, nuc.loc[list(keep)], n_hvg=2000, **kw)
        return p.loc[spots].to_numpy(float)

    # full-gene reproduction check (recorded, not gated)
    for arm in ("v2", "v2_noparam"):
        full = refit(arm, np.ones(len(gn), dtype=bool))
        meta[f"repro_max_abs_{arm}"] = float(np.abs(full - stored[arm].loc[spots][types].to_numpy(float)).max())
    shares_by_fold = {"v41": [p0.loc[spots].to_numpy(float)] * len(folds)}
    for arm in ("v2", "v2_noparam"):
        shares_by_fold[arm] = []
        for f in folds:
            keepmask = np.ones(len(gn), dtype=bool)
            keepmask[f] = False
            shares_by_fold[arm].append(refit(arm, keepmask))
    print("refits done", sample, meta, flush=True)

    def rho_for(arm, transform):
        out = np.full((len(uni_genes), len(types)), np.nan)
        for k, f in enumerate(folds):
            sh = transform(shares_by_fold[arm][k])[ok]
            fu = gene_pos_in_uni[f]
            fu = fu[fu >= 0]
            out[fu] = S.partial_rho(ex[:, fu], rl, sh)
        return pd.DataFrame(out, index=uni_genes, columns=types)

    rows = []
    for arm in ("v41", "v2", "v2_noparam"):
        evaluable = {t: all(S.evaluable_share(sh[:, j]) for sh in shares_by_fold[arm]) for j, t in enumerate(types)}
        conds = {"real": lambda x: x, "null": lambda x: S.permute_fraction(x, 100, seed, null_seed=seed)}
        if arm in ("v41", "v2"):
            for w in (25, 50, 75):
                conds[f"w{w}"] = (lambda ww: (lambda x: S.permute_fraction(x, ww, seed + 1000 * ww)))(w)
        for cond, tf in conds.items():
            rho = rho_for(arm, tf)
            genes = gene_lists_from_scores(rho, cand)
            for t in types:
                rows.append(dict(test="B", dataset="breast", arm=arm, condition=cond, unit=sample, type=t,
                                 evaluable=evaluable[t], **ora_row(genes[t], bg, t, mapping)))
                if cond == "null" and arm in ("v41", "v2"):
                    rows.append({**rows[-1], "condition": "w100"})  # w = 100 is the full permutation (rule)
        print("B done", sample, arm, flush=True)
    emit(Path(run) / "calls" / f"B_breast_{idx:02d}.csv", rows)
    (Path(run) / "calls" / f"B_breast_{idx:02d}_meta.json").write_text(json.dumps(meta, indent=2))


# ============================================================ Xenium
def _xen_markers():
    sys.path.insert(0, str(XB))
    import benchmark_constants as C

    out = set()
    for spec in C.MS_V41_7_CELL_PROFILE_DICT.values():
        for tier in spec.values():
            for ab in tier:
                ab = str(ab)
                out.add(ab[:-2] if ab.endswith("-1") else ab)
    return {m.upper() for m in out}


def _testable(panel_genes, mapping):
    """Per type: number of panel genes in the union of primary-library terms mapping to its lineage."""
    gmt = libs()[PRIMARY]
    pg = {g.upper() for g in panel_genes}
    out = {}
    for t, lin in mapping["type_to_lineage"].items():
        gs = set()
        for term, genes in gmt.items():
            if L.term_lineage(term, mapping) == lin:
                gs |= {g.upper() for g in genes}
        out[t] = len(gs & pg)
    return out


def _gt(r):
    return pd.read_csv(GT_PROP / f"Xenium_region_{r}_prop.csv", index_col=0)[XTYPES]


def xenium_sace(r, run):
    sys.path.insert(0, str(XB / "spots_baseline"))
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import run_gex_gate as gg
    import sace_props_ablation as sp

    run = Path(run)
    _, data = data_io.load_xenium(7)
    gt = _gt(r)
    props_root = run / "props"
    gg._write_props(gt, props_root, "oracle", r)
    perm, _ = L.spot_permuted(gt, XEN_SEED + r)
    gg._write_props(perm, props_root, "oracle_perm", r)
    uni = pd.DataFrame(1.0 / len(XTYPES), index=data[r]["obs"], columns=XTYPES)
    gg._write_props(uni, props_root, "uniform", r)
    dirs = {"oracle": props_root / "oracle", "oracle_perm": props_root / "oracle_perm", "uniform": props_root / "uniform",
            "v2": GATE_RUN / "props" / "v2_full", "v41": P138_PROPS}
    for arm, pdir in dirs.items():
        sp.run_region("marker_share", r, run / "sace" / arm, str(pdir), "ms_v41", cell_mode="no_round")
        print("sace done", arm, r, flush=True)


def _xen_layer_means(run, arm, r):
    base = Path(run) / "sace" / arm / f"Xenium_region_{r}"
    props = pd.read_csv(base / "props.csv", index_col=0)
    n = len(props)
    out = {}
    genes = None
    for t in XTYPES:
        df = pd.read_csv(base / "layers" / f"{XFILE[t]}_layer.csv", index_col=0)
        df.columns = df.columns.str.upper()
        df = df.loc[:, ~df.columns.duplicated()]
        out[t] = df.sum(axis=0) / n
        genes = out[t].index if genes is None else genes.union(out[t].index)
    return pd.DataFrame({t: out[t].reindex(genes).fillna(0.0) for t in XTYPES}).T, props


def xenium_A(r, run):
    mapping = load_mapping("term_mapping_xenium7.yaml")
    markers = _xen_markers()
    uni, _ = _xen_layer_means(run, "uniform", r)
    oracle_lay, _ = _xen_layer_means(run, "oracle", r)
    rows = []
    for arm, src, cond in (("oracle", "oracle", "real"), ("oracle", "oracle_perm", "null"), ("v2", "v2", "real"), ("v41", "v41", "real")):
        lay, props = _xen_layer_means(run, src, r)
        real_lay = lay if cond == "real" else oracle_lay
        expressed = real_lay.sum(axis=0) >= S.MIN_MEAN
        cand = [g for g in real_lay.columns[expressed.to_numpy()] if g not in markers]
        bg = list(real_lay.columns[expressed.to_numpy()])
        tb = _testable(real_lay.columns, mapping)
        dl = S.delta_table(lay, uni)
        genes = gene_lists_from_scores(dl.T, cand)
        for t in XTYPES:
            testable = tb[t] >= S.MIN_PANEL_GENES
            row = ora_row(genes[t], bg, t, mapping) if testable else dict(passed=False, best_rank=np.nan, nlp=0.0)
            rows.append(dict(test="A", dataset="xenium", arm=arm, condition=cond, unit=f"region_{r}", type=t,
                             evaluable=bool(S.evaluable_share(props[t]) and testable), testable=bool(testable),
                             n_panel_genes=tb[t], **row))
    emit(Path(run) / "calls" / f"A_xenium_{r}.csv", rows)


def xenium_B(r, run):
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    mapping = load_mapping("term_mapping_xenium7.yaml")
    markers = _xen_markers()
    seed = XEN_SEED + r
    _, data = data_io.load_xenium(7)
    d = data[r]
    gt = _gt(r)
    obs = list(d["obs"])
    genes = [g.upper() for g in d["genes"]]
    p0 = np.asarray(rs.ms_share(d["X"], d["names"], rs.PROF[7]))
    Y = d["Y"]
    lib_all = Y.sum(axis=1)
    sel = np.array([i for i, b in enumerate(obs) if b in gt.index and lib_all[i] > 0])
    mean = Y[sel].mean(axis=0)
    expressed = [g for g, m in zip(genes, mean) if m >= S.MIN_MEAN]
    cand = [g for g in expressed if g not in markers]
    tb = _testable(genes, mapping)
    ex, rl = S.prep_gene_ranks(np.log1p(1e4 * Y[sel] / lib_all[sel][:, None]), lib_all[sel])
    folds = S.make_folds(len(genes), seed)
    shares = {"v41": [p0] * 2, "oracle": [gt.reindex(obs).to_numpy(float)] * 2, "v2": []}
    for f in folds:
        mask = np.ones(len(genes), dtype=bool)
        mask[f] = False
        p, _ = refine_proportions(Y, d["nc"], p0, c=0.1, n_passes=3, model="multinomial", gene_mask=mask)
        shares["v2"].append(p)

    def rho_for(arm, tf):
        out = np.full((len(genes), len(XTYPES)), np.nan)
        for k, f in enumerate(folds):
            out[f] = S.partial_rho(ex[:, f], rl, tf(shares[arm][k][sel]))
        return pd.DataFrame(out, index=genes, columns=XTYPES)

    rows = []
    for arm in ("oracle", "v2", "v41"):
        conds = {"real": lambda x: x}
        if arm == "oracle":
            conds["null"] = lambda x: S.permute_fraction(x, 100, seed, null_seed=seed)
        for cond, tf in conds.items():
            gl = gene_lists_from_scores(rho_for(arm, tf), cand)
            for j, t in enumerate(XTYPES):
                testable = tb[t] >= S.MIN_PANEL_GENES
                row = ora_row(gl[t], expressed, t, mapping) if testable else dict(passed=False, best_rank=np.nan, nlp=0.0)
                ev = bool(all(S.evaluable_share(sh[sel][:, j]) for sh in shares[arm]) and testable)
                rows.append(dict(test="B", dataset="xenium", arm=arm, condition=cond, unit=f"region_{r}", type=t,
                                 evaluable=ev, testable=bool(testable), n_panel_genes=tb[t], **row))
    emit(Path(run) / "calls" / f"B_xenium_{r}.csv", rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["breast-A", "breast-B", "xenium-sace", "xenium-A", "xenium-B"])
    ap.add_argument("index", type=int)
    ap.add_argument("run")
    a = ap.parse_args()
    {"breast-A": breast_A, "breast-B": breast_B, "xenium-sace": xenium_sace, "xenium-A": xenium_A, "xenium-B": xenium_B}[a.cmd](a.index, a.run)


if __name__ == "__main__":
    main()
