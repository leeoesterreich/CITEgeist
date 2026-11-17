"""Column (per-antibody) cosine scaling test for the marker-share estimator. Information + decision script.

Arms: P plain, B GMM brightness, C column cosine (each channel / its L2 norm across spots of the region) then plain
marker share, C2 = explicit row-cosine of the column-scaled matrix then marker share (identical to C: row scale cancels).
Reuses marker_share_eval.py (RCC, old SingleR GT, evaluate_all_regions) and marker_share_sim_eval.py (sims,
protein_sensitivity.baseline) logic verbatim. New RNA-cluster GT scoring is information only; writes nothing to results_rnagt/.

Usage: column_scaling_test.py OUT_DIR
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
PRODUCED_BY = "Benchmarking/xenium_benchmarking/spots_baseline/column_scaling_test.py"
ARMS = ("P", "B", "C", "C2")


def column_cosine_scale(X):
    """Divide each column (antibody) by its L2 norm across spots; all-zero columns stay zero."""
    X = np.asarray(X, dtype=float)
    norm = np.linalg.norm(X, axis=0, keepdims=True)
    return np.divide(X, norm, out=np.zeros_like(X), where=norm > 0)


def arm_props(arm, X, vn, profile, markers, ms, bright):
    """Return (DataFrame, n_zero) for one arm."""
    if arm == "P":
        return ms.marker_share_proportions(X, vn, profile)
    if arm == "B":
        return ms.marker_share_proportions(X, vn, profile, brightness=bright)
    Xc = column_cosine_scale(X)
    if arm == "C":  # no row normalisation at all: the row scale cancels in the renormalisation anyway
        types = list(profile)
        sc = np.column_stack([Xc[:, [vn.index(m) for m in ms._major_markers(profile[t])]].mean(axis=1) for t in types])
        den = sc.sum(axis=1)
        fr = np.divide(sc, den[:, None], out=np.zeros_like(sc), where=den[:, None] > 0)
        fr[den == 0] = 1.0 / len(types)
        return pd.DataFrame(fr, columns=types), int((den == 0).sum())
    return ms.marker_share_proportions(ms.cosine_normalize(Xc), vn, profile)  # C2


def main(out_dir):
    import scanpy as sc

    sys.path[:0] = [str(REPO), str(XB), str(XB / "evaluation" / "src"), str(REPO / "scripts"), str(REPO / "Benchmarking")]
    from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT, ACHIEVABLE_7_CELL_PROFILE_DICT
    import benchmark_paths as bp
    from evaluate_benchmark import (common_spot_masks, evaluability_table, evaluate_all_regions, evaluate_strict)
    from protein_sensitivity import baseline as B
    from provenance import write_stamped_json

    from CITEgeist.model.deconvolution import marker_share as ms

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cites = {}
    for i in range(5):
        cite = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
        X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
        cites[i] = (X, list(cite.var_names), cite.obs_names)
    summary = {"rcc_old_gt": {}, "sim": {}, "rnagt_info": {}, "composition": {}}
    stages = {}  # (n, arm) -> staged dir of predictions
    for n, profile in ((6, ACHIEVABLE_6_CELL_PROFILE_DICT), (7, ACHIEVABLE_7_CELL_PROFILE_DICT)):
        markers = sorted({m for v in profile.values() for m in v["Major"]})
        bright = {i: ms.estimate_marker_brightness(X, vn, markers) for i, (X, vn, _) in cites.items()}
        gt = XB / "ground_truth_singler" / f"ground_truth_{n}type"
        comp = {}
        for arm in ARMS:
            stage = Path(tempfile.mkdtemp(prefix=f"cs_{arm}_{n}_"))
            stages[(n, arm)] = stage
            allp = []
            for i, (X, vn, obs) in cites.items():
                df, _ = arm_props(arm, X, vn, profile, markers, ms, bright[i])
                df.index = obs
                df.to_csv(stage / f"Xenium_region_{i}_cell_prop_global_results.csv")
                allp.append(df)
            comp[arm] = pd.concat(allp).mean().round(5).to_dict()
            combine = n == 6
            res = evaluate_all_regions(gt_dir=str(gt), pred_dir=str(stage), n_regions=5, output_path=None,
                                       prefix="Xenium", use_achievable_7=not combine, combine_t_cells=combine)
            summary["rcc_old_gt"][f"{n}_{arm}"] = {"pooled_r": res["summary"]["overall_mean_pearson_r"]}
            print(f"RCC {n}-type {arm}: r={res['summary']['overall_mean_pearson_r']:.6f}", flush=True)
        summary["composition"][str(n)] = comp
        assert np.allclose(pd.Series(comp["C"]).values, pd.Series(comp["C2"]).values, atol=1e-4)
    # sims
    sprofile = B.sim_profile()
    smarkers = sorted({m for v in sprofile.values() for m in v["Major"]})
    units = {a: [] for a in ARMS}
    for cond in ("high_seg", "mixed"):
        for rep in range(5):
            _g, cite = B.load_sim(cond, rep)
            X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
            vn = list(cite.var_names)
            br = ms.estimate_marker_brightness(X, vn, smarkers)
            for arm in ARMS:
                df, _ = arm_props(arm, X, vn, sprofile, smarkers, ms, br)
                df.index = cite.obs_names
                units[arm].append({"condition": cond, "replicate": rep, **B.sim_prop_metrics(df, cond, rep)})
    for arm, rows in units.items():
        summary["sim"][arm] = {"pooled_r": float(np.mean([r["pearson_r"] for r in rows])),
                               "high_seg_r": float(np.mean([r["pearson_r"] for r in rows if r["condition"] == "high_seg"])),
                               "mixed_r": float(np.mean([r["pearson_r"] for r in rows if r["condition"] == "mixed"]))}
        print(f"SIM {arm}: {summary['sim'][arm]}", flush=True)
    # new RNA-cluster GT, strict evaluator, information only (reads only; writes only under out_dir)
    for n in (6, 7):
        gdir = bp.gt_dir(n, "rnagt")
        gcols = pd.read_csv(Path(gdir) / "ground_truth" / "Xenium_region_0_prop.csv", index_col=0, nrows=1).columns
        types = [c for c in gcols if c not in ("n_cells", "spot_x", "spot_y")]
        evaluable = evaluability_table(gdir, types, 5, "Xenium")
        methods = {m: (bp.citegeist_pred_dir(n) if m == "CITEgeist" else bp.pred_dir(m, n, "rnagt"))
                   for m in ["CITEgeist", "RCTD", "Seurat", "CARD", "Cell2Location", "Tangram"]}
        dirs = {**{m: str(d) for m, d in methods.items()}, **{f"arm_{a}": str(stages[(n, a)]) for a in ARMS}}
        masks = common_spot_masks(gdir, dirs, evaluable, 5, "Xenium")
        for arm in ARMS:
            r = evaluate_strict(gdir, str(stages[(n, arm)]), types, masks, evaluable, 5, "Xenium")["summary"]
            summary["rnagt_info"][f"{n}_{arm}"] = {"per_type_mean_r": r["per_type_mean_r"],
                                                  "flat_r": r["overall_mean_pearson_r"],
                                                  "per_type": {t: r.get(f"{t}_mean_r") for t in types}}
            print(f"RNAGT {n}-type {arm}: per_type_mean_r={r['per_type_mean_r']:.4f} flat_r={r['overall_mean_pearson_r']:.4f}", flush=True)
        # cross-check P vs the official CITEgeist rnagt result (read-only)
        off = bp.eval_results_dir("rnagt") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json"
        if off.is_file():
            summary["rnagt_info"][f"{n}_official_CITEgeist"] = json.loads(off.read_text())["summary"]["per_type_mean_r"]
    for s in stages.values():
        shutil.rmtree(s)
    write_stamped_json(out_dir / "column_scaling_results.json", summary, PRODUCED_BY)
    print("wrote", out_dir)


if __name__ == "__main__":
    main(sys.argv[1])
