"""Write CITEgeist v41 marker-share proportions for the Xenium RCC regions (6- and 7-type) plus stamped results.

Input : raw CITE .X of xenium_pseudovisium/data_protein_gt/h5ad_objects/Xenium_region_{i}_CITE.h5ad (same as
        marker_share_eval.py), ACHIEVABLE_{6,7}_CELL_PROFILE_DICT Major markers, plain estimator (brightness=None).
Output: <MAIN>/Benchmarking/xenium_benchmarking/CITEgeist/output_marker_share_{N}type/Xenium_region_{i}/predictions/baseline.csv.gz
        <MAIN>/Benchmarking/xenium_benchmarking/evaluation/results_v41/method_comparison_singler_{N}type/CITEgeist_{N}type_results.json
        (schema/stamp identical to protein_sensitivity/write_citegeist_results.py).
MAIN (env, default = this checkout) is the main repo root outputs are written under, so a clean worktree can run it.

Usage: write_marker_share_preds.py N_TYPES   (N_TYPES in 6, 7)   [legacy: ACHIEVABLE profile, no brightness, scores + stamps]
       write_marker_share_preds.py N_TYPES --out-root DIR [--profile {achievable,ms_v41}] [--brightness {none,gmm}]
           pure mode: only writes DIR/Xenium_region_{i}/predictions/baseline.csv.gz (no scoring, no results json).
           brightness=gmm estimates marker brightness PER REGION on that region's spots (column_scaling_test.py convention).
"""
import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
PRODUCED_BY = "Benchmarking/xenium_benchmarking/spots_baseline/write_marker_share_preds.py"
N_REGIONS = 5


def get_profile(name, n):
    """Profile dict for --profile {achievable,ms_v41} and n types (6 or 7)."""
    sys.path.insert(0, str(XB))
    import benchmark_constants as C

    prefix = {"achievable": "ACHIEVABLE", "ms_v41": "MS_V41"}[name]
    return getattr(C, f"{prefix}_{n}_CELL_PROFILE_DICT")


def write_region_preds(X, var_names, obs_names, profile, region_dir, brightness="none"):
    """Write region_dir/predictions/baseline.csv.gz (spot-indexed, one column per type); return the DataFrame.

    brightness="gmm": b_m estimated here, on THIS region's X, for the sorted Major markers of the profile.
    """
    from CITEgeist.model.deconvolution import marker_share as ms

    assert brightness in ("none", "gmm"), brightness
    bright = None
    if brightness == "gmm":
        markers = sorted({m for v in profile.values() for m in ms._major_markers(v)})
        bright = ms.estimate_marker_brightness(X, list(var_names), markers)
    df, _ = ms.marker_share_proportions(X, var_names, profile, brightness=bright)
    df.index = pd.Index(obs_names)
    pred = Path(region_dir) / "predictions"
    pred.mkdir(parents=True, exist_ok=True)
    df.to_csv(pred / "baseline.csv.gz")
    return df


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("n", type=int, choices=(6, 7))
    ap.add_argument("--profile", choices=("achievable", "ms_v41"), default="achievable")
    ap.add_argument("--brightness", choices=("none", "gmm"), default="none")
    ap.add_argument("--out-root", default=None, help="pure mode: write predictions only, under this dir")
    a = ap.parse_args(argv)
    if a.out_root is None and (a.profile != "achievable" or a.brightness != "none"):
        ap.error("non-default --profile/--brightness require --out-root (legacy mode overwrites the v41 results json)")
    return a


def write_preds_only(n, profile_name, brightness, out_root):
    """Pure step: per-region predictions for all regions; returns {region: DataFrame}. No scoring."""
    import scanpy as sc

    sys.path.insert(0, str(REPO))
    profile = get_profile(profile_name, n)
    out = {}
    for i in range(N_REGIONS):
        cite = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
        X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
        out[i] = write_region_preds(X, list(cite.var_names), cite.obs_names, profile,
                                    Path(out_root) / f"Xenium_region_{i}", brightness=brightness)
        print(f"region {i}: {out[i].shape} mean={out[i].mean().round(4).to_dict()}", flush=True)
    return out


def main(n):
    import scanpy as sc

    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(XB))
    sys.path.insert(0, str(XB / "evaluation" / "src"))
    sys.path.insert(0, str(REPO / "scripts"))
    from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT, ACHIEVABLE_7_CELL_PROFILE_DICT
    from benchmark_paths import citegeist_pred_dir, eval_results_dir
    from evaluate_benchmark import evaluate_all_regions
    from provenance import write_stamped_json

    assert n in (6, 7)
    main_root = Path(os.environ.get("MAIN", REPO)).resolve()
    rebase = lambda p: main_root / Path(p).relative_to(REPO)  # noqa: E731
    pred_root = rebase(citegeist_pred_dir(n))
    profile = {6: ACHIEVABLE_6_CELL_PROFILE_DICT, 7: ACHIEVABLE_7_CELL_PROFILE_DICT}[n]
    gt = XB / "ground_truth_singler" / f"ground_truth_{n}type"
    stage = Path(tempfile.mkdtemp(prefix=f"msp_{n}_"))
    run_dirs = {}
    for i in range(N_REGIONS):
        cite = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
        X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
        rd = pred_root / f"Xenium_region_{i}"
        write_region_preds(X, list(cite.var_names), cite.obs_names, profile, rd)
        # evaluate from the file actually written to disk
        pd.read_csv(rd / "predictions" / "baseline.csv.gz", index_col=0).to_csv(
            stage / f"Xenium_region_{i}_cell_prop_global_results.csv"
        )
        run_dirs[str(i)] = str(rd.relative_to(main_root))
    combine = n == 6
    res = evaluate_all_regions(
        gt_dir=str(gt), pred_dir=str(stage), n_regions=N_REGIONS, output_path=None, prefix="Xenium",
        use_achievable_7=not combine, combine_t_cells=combine,
    )
    shutil.rmtree(stage)
    assert res["summary"]["n_regions"] == N_REGIONS, res["summary"]["n_regions"]
    s = res["summary"]
    print(f"{n}-type marker-share: r={s['overall_mean_pearson_r']:.6f} rmse={s['overall_mean_rmse']:.4f} "
          f"mae={s['overall_mean_mae']:.4f}", flush=True)
    res["_inputs"] = {
        "run_dirs": run_dirs,
        "estimator": "marker_share_proportions (CITEgeist.model.deconvolution.marker_share), brightness=None",
        "marker_dict": f"benchmark_constants.ACHIEVABLE_{n}_CELL_PROFILE_DICT (Major markers)",
        "evaluator": "evaluate_benchmark.evaluate_all_regions (compare_all_methods.py settings)",
        "gt_dir": str(gt.relative_to(REPO)),
        "cite_dir": str(DATA.relative_to(REPO)),
    }
    out = rebase(eval_results_dir("march")) / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    write_stamped_json(out, res, PRODUCED_BY)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    if args.out_root is None:
        main(args.n)
    else:
        write_preds_only(args.n, args.profile, args.brightness, args.out_root)
