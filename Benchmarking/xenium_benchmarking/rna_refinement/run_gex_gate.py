"""Task 7: GEX gate (Gate 4), gene cross-fit and RNA-GT annotation sensitivity for engine v2.

Modes (run from the repo root with the CITEgeist_env python):
  region RUN_DIR R   refine proportions (full fit + 5 gene-withheld fits) for Xenium region R and run SACE on
                     v2-full, the 5 cross-fit arms and the v41 marker-share parity arm (props_dir = official p138).
  score  RUN_DIR     score every arm with the v41 GEX scorer, apply Gate 4, write RUN_DIR/gex_gate.json.
RUN_DIR must be NEW (created by `region 0`/first call; never an existing result folder).

Ground truth is never a model input: proportions come from protein (marker share) + raw RNA + nuclei only.
Gate 4 (spec gate): v2 SACE has the lowest per-type RMSE in >= 5 of 7 types vs {v41 CITEgeist (P138), Cell2Location,
Tangram} read from the v41 comparison CSV. No retuning.
"""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
FINAL_VARIANT = REPO / "Benchmarking/midkine_remediation/results_engine_v2/r1_20261031_1442/final_variant.json"
VARIANTS_RUN = HERE / "runs" / "variants_20261031_1404_r8"
GEX_EVAL = XB / "evaluation" / "rnaprot_v41"
P138_PROPS = XB / "CITEgeist" / "output_marker_share_p138_7type"
P138_SACE = XB / "CITEgeist" / "results_sace_marker_share_p138_7type"
K_FOLDS = 5
N_REGIONS = 5
GATE_MIN_TYPES = 5
COMPETITORS = ("Cell2Location", "Tangram")
CITEGEIST_ROWS = ("CITEgeist_SACE", "CITEgeist_v41p_SACE", "CITEgeist_p138_SACE")
PARITY_TOL = 0.005


def gene_folds(n_genes, k=K_FOLDS, seed=0):
    """Deterministic gene folds: one default_rng(seed) permutation split into k nearly equal parts."""
    return [np.sort(f) for f in np.array_split(np.random.default_rng(seed).permutation(n_genes), k)]


def crossfit_masks(n_genes, k=K_FOLDS, seed=0):
    """List of (fold_gene_indices, fit_mask) where fit_mask is bool[G] with the fold's genes set False."""
    out = []
    for fold in gene_folds(n_genes, k, seed):
        mask = np.ones(n_genes, dtype=bool)
        mask[fold] = False
        out.append((fold, mask))
    return out


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _variant():
    v = json.loads(FINAL_VARIANT.read_text())
    assert v["variant"] == "a_V3_multinomial" and v["model"] == "multinomial", v
    return v


def _write_props(props, props_root, arm, region):
    d = Path(props_root) / arm / f"Xenium_region_{region}" / "predictions"
    d.mkdir(parents=True, exist_ok=True)
    props.to_csv(d / "baseline.csv.gz")


# ------------------------------------------------------------------ region mode
def region_job(run_dir, region):
    sys.path.insert(0, str(REPO))
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    run_dir = Path(run_dir)
    v = _variant()
    types, data = data_io.load_xenium(7)
    genes = list(data[region]["genes"])
    for r in range(N_REGIONS):
        assert list(data[r]["genes"]) == genes, "gene panel differs across regions; folds would not be shared"
    d = data[region]
    p0 = rs.ms_share(d["X"], d["names"], rs.PROF[7])
    p0_df = pd.DataFrame(np.asarray(p0), index=d["obs"], columns=types)
    kw = dict(c=v["c"], n_passes=v["n_passes"], model=v["model"])

    checks = {}
    p138 = pd.read_csv(P138_PROPS / f"Xenium_region_{region}" / "predictions" / "baseline.csv.gz", index_col=0)
    checks["p0_vs_official_p138_max_abs"] = float(np.abs(p138.reindex(d["obs"])[types].values - p0_df.values).max())

    p_full, _ = refine_proportions(d["Y"], d["nc"], np.asarray(p0), **kw)
    saved = pd.read_csv(VARIANTS_RUN / "artifacts" / v["variant"] / "7type" / f"region_{region}" / "cell_prop.csv", index_col=0)
    checks["full_fit_vs_saved_artifact_max_abs"] = float(np.abs(saved.reindex(d["obs"])[types].values - p_full).max())
    assert checks["full_fit_vs_saved_artifact_max_abs"] < 1e-8, checks

    props_root = run_dir / "props"
    _write_props(pd.DataFrame(p_full, index=d["obs"], columns=types), props_root, "v2_full", region)
    folds = {}
    for k, (fold, mask) in enumerate(crossfit_masks(len(genes))):
        assert not mask[fold].any()
        p_k, _ = refine_proportions(d["Y"], d["nc"], np.asarray(p0), gene_mask=mask, **kw)
        _write_props(pd.DataFrame(p_k, index=d["obs"], columns=types), props_root, f"fold{k}", region)
        folds[k] = dict(n_withheld=int(len(fold)), n_fit_genes=int(mask.sum()),
                        max_abs_change_vs_full=float(np.abs(p_k - p_full).max()))
    (run_dir / "region_checks").mkdir(exist_ok=True)
    (run_dir / "region_checks" / f"region_{region}.json").write_text(json.dumps(dict(checks=checks, folds=folds), indent=2))
    print("region", region, checks, flush=True)

    sys.path.insert(0, str(XB / "spots_baseline"))
    import sace_props_ablation as sp

    arms = [("v2_full", props_root / "v2_full")] + [(f"fold{k}", props_root / f"fold{k}") for k in range(K_FOLDS)]
    arms.append(("p138_parity", P138_PROPS))
    for arm, pdir in arms:
        sp.run_region("marker_share", region, run_dir / "sace" / arm, str(pdir), "ms_v41")


# ------------------------------------------------------------------ scoring
def _load_scorer():
    os.environ["XEN_LABELSET"] = "rnaprot_v41"
    sys.path[:0] = [str(REPO), str(XB), str(XB / "evaluation" / "src")]
    import benchmark_paths as bp
    import compare_xenium_gex_methods as cx

    return bp, cx


def score_layers(layers_root, bp, cx, genes=None):
    """Mirror of cx.evaluate_method_on_region (same metrics, same filters) without the strict sha gate.

    ``genes`` (uppercase names) restricts scoring to those genes; the GT-expressed-gene filter is then applied
    within the subset (it is per-gene, so equivalent to filtering first). Returns (rows, excluded).
    """
    gt_root = bp.gt_gex_dir("rnaprot_v41")
    rows, excluded = [], []
    for rid in range(N_REGIONS):
        ldir = Path(layers_root) / f"Xenium_region_{rid}" / "layers"
        for stem in cx.CELL_TYPES:
            gt_df = cx.load_csv_spots_genes(gt_root / f"Xenium_region_{rid}" / f"{stem}_GT.csv")
            lp = cx.find_layer_file(ldir, stem)
            if lp is None:
                excluded.append(dict(region_id=rid, cell_type=stem, reason="layer_missing"))
                continue
            pred_df = cx.load_csv_spots_genes(lp)
            gt_df.columns, pred_df.columns = gt_df.columns.str.upper(), pred_df.columns.str.upper()
            gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
            pred_df = pred_df.loc[:, ~pred_df.columns.duplicated()]
            spots = sorted(set(gt_df.index) & set(pred_df.index))
            cg = sorted(set(gt_df.columns) & set(pred_df.columns))
            if genes is not None:
                cg = [g for g in cg if g in genes]
            if len(spots) < 5 or len(cg) < 5:
                excluded.append(dict(region_id=rid, cell_type=stem, reason="insufficient_overlap"))
                continue
            gv = gt_df.loc[spots, cg].values.astype(np.float64)
            pv = pred_df.loc[spots, cg].values.astype(np.float64)
            if gv.sum() == 0.0:
                excluded.append(dict(region_id=rid, cell_type=stem, reason="gt_all_zero"))
                continue
            keep = gv.sum(axis=0) > 0
            if keep.sum() < 5:
                excluded.append(dict(region_id=rid, cell_type=stem, reason="few_expressed_genes"))
                continue
            m = cx.compute_metrics(pv[:, keep], gv[:, keep])
            rows.append(dict(region_id=rid, cell_type=stem, n_gt_spots=len(gt_df.index), n_shared_spots=len(spots),
                             n_shared_genes=int(keep.sum()), **m))
    return rows, excluded


def per_type(rows, col="rmse"):
    df = pd.DataFrame(rows)
    return {t: float(df[df.cell_type == t][col].dropna().mean()) for t in df.cell_type.unique()}


def _baseline_tables():
    det = pd.read_csv(GEX_EVAL / "xenium_gex_comparison_7type.csv")
    summ = pd.read_csv(GEX_EVAL / "xenium_gex_comparison_7type_summary.csv")
    return det, summ


def _baseline_rmse(summ, method):
    s = summ[(summ.method == method) & (summ.cell_type != "ALL")]
    return {r.cell_type: float(r.rmse_mean) for r in s.itertuples()}


def gate4(v2, v41, comps, types):
    lowest = {}
    for t in types:
        others = {"v41": v41[t], **{m: c[t] for m, c in comps.items() if t in c}}
        lowest[t] = bool(v2[t] < min(others.values()))
    return lowest


def annotation_sensitivity(run_dir):
    """Proportions of v2 / marker share scored on the RNA-only GT (xenium.gt_rnagt) with the strict scorer."""
    sys.path[:0] = [str(REPO / "scripts"), str(XB / "evaluation" / "src"), str(XB)]
    import benchmark_paths as bp
    from evaluate_benchmark import common_spot_masks, evaluability_table, evaluate_strict
    from gt_rebuild.score_methods import N_REGIONS as NR, PREFIX, gt_types, pred_dirs_for

    out = {}
    for n in (6, 7):
        gtd = bp.gt_dir(n, "rnagt")
        types = gt_types(gtd, n)
        std = {m: str(p) for m, p in pred_dirs_for("rnagt", n).items()}
        ev = evaluability_table(gtd, types, NR, PREFIX)
        masks = common_spot_masks(gtd, std, ev, NR, PREFIX)
        row = {}
        for arm, fname in (("v2", "cell_prop.csv"), ("marker_share_p0", "p0.csv")):
            odir = Path(run_dir) / "annotation_sensitivity" / f"{arm}_{n}type"
            odir.mkdir(parents=True, exist_ok=True)
            for r in range(NR):
                df = pd.read_csv(VARIANTS_RUN / "artifacts" / "a_V3_multinomial" / f"{n}type" / f"region_{r}" / fname, index_col=0)
                df[types].to_csv(odir / f"Xenium_region_{r}_cell_prop_global_results.csv")
            res = evaluate_strict(gtd, str(odir), types, masks, ev, NR, PREFIX)
            row[arm] = float(res["summary"]["per_type_mean_r"])
        comp = {}
        for m, p in std.items():
            res = evaluate_strict(gtd, p, types, masks, ev, NR, PREFIX)
            comp[m] = float(res["summary"]["per_type_mean_r"])
        row["competitors_and_v41_same_masks"] = comp
        out[f"{n}type"] = row
    return out


def score_job(run_dir):
    run_dir = Path(run_dir)
    bp, cx = _load_scorer()
    det, summ = _baseline_tables()
    types = list(cx.CELL_TYPES)
    sys.path.insert(0, str(REPO))
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io

    genes = [g.upper() for g in data_io.load_xenium(7)[1][0]["genes"]]

    # 1. scorer validation: re-score the OFFICIAL p138 layers and compare to the shipped CSV
    off_rows, off_exc = score_layers(P138_SACE, bp, cx)
    off = per_type(off_rows)
    csv_p138 = _baseline_rmse(summ, "CITEgeist_p138_SACE")
    scorer_diff = {t: abs(off[t] - csv_p138[t]) for t in types}
    scorer_ok = max(scorer_diff.values()) < 1e-3  # CSV rounds to 4 d.p.

    # 2. standard scoring
    v2_rows, v2_exc = score_layers(run_dir / "sace" / "v2_full", bp, cx)
    par_rows, par_exc = score_layers(run_dir / "sace" / "p138_parity", bp, cx)
    v2, par = per_type(v2_rows), per_type(par_rows)
    v41 = csv_p138
    comps = {m: _baseline_rmse(summ, m) for m in COMPETITORS}
    others_cg = {m: _baseline_rmse(summ, m) for m in CITEGEIST_ROWS}
    lowest = gate4(v2, v41, comps, types)
    n_lowest = sum(lowest.values())
    parity = {t: par[t] - v41[t] for t in types}
    parity_ok = max(abs(x) for x in parity.values()) <= PARITY_TOL
    table = {t: dict(v2=v2[t], v41_p138=v41[t], **{m: comps[m][t] for m in COMPETITORS},
                     best_competitor=min(comps[m][t] for m in COMPETITORS),
                     v2_minus_v41=v2[t] - v41[t], v2_lowest=lowest[t]) for t in types}

    # 3. cross-fit
    gold = {}
    crossfit = {"folds": {}}
    folds = gene_folds(len(genes))
    ref_rows, _ = off_rows, off_exc
    for k, fold in enumerate(folds):
        fg = {genes[i] for i in fold}
        rows_k, exc_k = score_layers(run_dir / "sace" / f"fold{k}", bp, cx, genes=fg)
        rows_full, _ = score_layers(run_dir / "sace" / "v2_full", bp, cx, genes=fg)
        rows_v41, _ = score_layers(P138_SACE, bp, cx, genes=fg)
        crossfit["folds"][str(k)] = dict(
            genes_withheld=[genes[i] for i in fold],
            rmse_withheld_genes_crossfit=per_type(rows_k),
            rmse_withheld_genes_v2_full_fit=per_type(rows_full),
            rmse_withheld_genes_v41_p138=per_type(rows_v41),
            n_blocks=len(rows_k), excluded=exc_k,
        )
    cf_mean = {t: float(np.nanmean([crossfit["folds"][str(k)]["rmse_withheld_genes_crossfit"].get(t, np.nan) for k in range(K_FOLDS)])) for t in types}
    v41_mean = {t: float(np.nanmean([crossfit["folds"][str(k)]["rmse_withheld_genes_v41_p138"].get(t, np.nan) for k in range(K_FOLDS)])) for t in types}
    full_mean = {t: float(np.nanmean([crossfit["folds"][str(k)]["rmse_withheld_genes_v2_full_fit"].get(t, np.nan) for k in range(K_FOLDS)])) for t in types}
    crossfit["mean_over_folds"] = dict(crossfit=cf_mean, v2_full_fit=full_mean, v41_p138=v41_mean,
                                       types_crossfit_below_v41=sum(cf_mean[t] < v41_mean[t] for t in types))

    # 4. annotation sensitivity (RNA-only GT)
    annot = annotation_sensitivity(run_dir)

    blocks = lambda rows: [dict(region=r["region_id"], cell_type=r["cell_type"], n_gt_spots=r["n_gt_spots"],
                                n_shared_spots=r["n_shared_spots"], n_shared_genes=r["n_shared_genes"]) for r in rows]
    base_blocks = sorted((int(r.region_id), r.cell_type) for r in det[det.method == "CITEgeist_p138_SACE"].itertuples())
    result = dict(
        gate="Gate 4: v2 SACE lowest RMSE in >= 5 of 7 types vs {v41 P138, Cell2Location, Tangram} (standard scoring)",
        variant=_variant()["variant"],
        gate4_pass=bool(n_lowest >= GATE_MIN_TYPES),
        gate4_status="PASS" if n_lowest >= GATE_MIN_TYPES else "FAIL",
        n_types_v2_lowest=int(n_lowest),
        per_type_rmse=table,
        other_citegeist_rows_rmse=others_cg,
        parity=dict(note="same pipeline, official P138 marker-share props -> SACE, vs v41 CSV", per_type_diff=parity,
                    tolerance=PARITY_TOL, ok=bool(parity_ok)),
        scorer_validation=dict(note="official P138 layers re-scored with this script vs CSV (CSV rounded to 4 dp)",
                               per_type_abs_diff=scorer_diff, ok=bool(scorer_ok)),
        crossfit=crossfit,
        annotation_sensitivity=annot,
        producer=dict(script="Benchmarking/xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py",
                      script_sha256=_sha(XB / "evaluation/src/compare_xenium_gex_methods.py"),
                      run_gex_gate_sha256=_sha(__file__),
                      git_head=subprocess.check_output(["git", "-C", str(REPO), "--no-optional-locks", "rev-parse", "HEAD"], text=True).strip(),
                      complete_json=json.loads((GEX_EVAL / "_COMPLETE.json").read_text())),
        baseline_csv_sha256=dict(detail=_sha(GEX_EVAL / "xenium_gex_comparison_7type.csv"),
                                 summary=_sha(GEX_EVAL / "xenium_gex_comparison_7type_summary.csv")),
        block_manifest=dict(v2=blocks(v2_rows), v2_excluded=v2_exc, p138_parity=blocks(par_rows),
                            baseline_p138_blocks_match_v2=bool(base_blocks == sorted((r["region_id"], r["cell_type"]) for r in v2_rows))),
    )
    (run_dir / "gex_gate.json").write_text(json.dumps(result, indent=2, default=float))
    print(json.dumps({k: result[k] for k in ("gate4_status", "n_types_v2_lowest")}), flush=True)
    for t in types:
        print(f"{t:16s} " + " ".join(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}" for k, v in table[t].items()), flush=True)
    print("parity", parity_ok, "scorer", scorer_ok, flush=True)


if __name__ == "__main__":
    mode, run_dir = sys.argv[1], sys.argv[2]
    if mode == "region":
        Path(run_dir).mkdir(parents=True, exist_ok=True)
        region_job(run_dir, int(sys.argv[3]))
    elif mode == "score":
        score_job(run_dir)
    else:
        raise SystemExit(f"unknown mode {mode}")
