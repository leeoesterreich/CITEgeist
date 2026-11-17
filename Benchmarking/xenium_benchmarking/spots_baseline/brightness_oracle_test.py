"""Xenium brightness oracle test, run EXACTLY as fixed before the run.

Arms (MS_V41 dictionary, B Major = CD20 + CD138; strict per-type mean r on the results_rnaprot shared masks):
  plain  : marker share, b_m = 1 (= P138)
  gmm    : GMM brightness, max_components 3, per region (= v41p production)
  oracle : b*_m = OLS slope of the cosine-normalised marker value on the summed GT proportion of its Major type(s),
           all spots per region (rule AMENDMENT 1); slope <= 0.001 -> median b* of other markers (reported). DIAGNOSTIC ONLY.
Gates (abort before writing results): plain == P138 files and published SPOTS_P138 / brightness_isolation values;
gmm == published CITEgeist_v41p (1e-12) and v41p files.
Usage: brightness_oracle_test.py OUT_DIR  (must not exist; run from a clean worktree; env MAIN = main repo root)
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import settings_v41p as S  # noqa: E402  (sets sys.path, imports C, bp, provenance, ms, evaluator)

C, bp, ms, provenance = S.C, S.bp, S.ms, S.provenance
N_REGIONS, PREFIX = S.N_REGIONS, S.PREFIX
PRODUCER = "Benchmarking/xenium_benchmarking/spots_baseline/brightness_oracle_test.py"
MIN_SLOPE = 0.001
ARMS = ("plain", "gmm", "oracle")


def oracle_brightness(Xc, vn, profile, gt):
    """b*_m per marker (AMENDMENT 1): OLS slope (with intercept) of the marker's cosine-normalised value on the summed
    GT proportion of the type(s) listing it as Major, across all GT spots of the region. Slope <= 0.001 -> median b* of
    the other markers (reported). GT only; diagnostic."""
    slope, info = {}, {}
    for m in sorted({m for v in profile.values() for m in ms._major_markers(v)}):
        types = [t for t, v in profile.items() if m in ms._major_markers(v)]
        ab = gt[types].sum(axis=1).values
        col = Xc[:, vn.index(m)]
        slope[m] = float(np.polyfit(ab, col, 1)[0])
        info[m] = {"types": types, "n_gt_spots_used": int(len(gt)), "raw_slope": slope[m],
                   "mean_gt_abundance": float(ab.mean())}
    bstar = {}
    for m in slope:
        if slope[m] > MIN_SLOPE:
            bstar[m] = slope[m]
            info[m]["median_fallback"] = False
    for m in slope:
        if m not in bstar:
            others = [slope[k] for k in slope if k != m and slope[k] > MIN_SLOPE]
            if not others:
                raise RuntimeError(f"no positive slopes to take a median from (marker {m})")
            bstar[m] = float(np.median(others))
            info[m]["median_fallback"] = True
    return bstar, info


def write_pred(df, root, arm, n, i):
    d = root / "preds" / f"{arm}_{n}type" / f"Xenium_region_{i}" / "predictions"
    d.mkdir(parents=True, exist_ok=True)
    df.to_csv(d / "baseline.csv.gz")


def run(n, root):
    import scanpy as sc

    gtd = bp.gt_dir(n, S.LABELSET)
    types = S.gt_types(gtd, n)
    profile = getattr(C, f"MS_V41_{n}_CELL_PROFILE_DICT")
    assert sorted(profile) == sorted(types), (sorted(profile), sorted(types))
    diag = {}
    for i in range(N_REGIONS):
        cite = sc.read_h5ad(S.DATA / f"Xenium_region_{i}_CITE.h5ad")
        X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
        vn = list(cite.var_names)
        gt = pd.read_csv(Path(gtd) / "ground_truth" / f"{PREFIX}_region_{i}_prop.csv", index_col=0)
        gt = gt.reindex(cite.obs_names)
        has_gt = gt[types].notna().all(axis=1).values  # spots without a GT row cannot enter the oracle fit
        gt = gt[has_gt]
        markers = sorted({m for v in profile.values() for m in ms._major_markers(v)})
        bgmm, _ = ms.estimate_marker_brightness_with_diagnostics(X, vn, markers, max_components=3)
        Xc = ms.cosine_normalize(X)  # same normalisation marker_share uses internally
        bstar, info = oracle_brightness(Xc[has_gt], vn, profile, gt)
        if not all(np.isfinite(v) and v > 0 for v in bstar.values()):
            raise RuntimeError(f"region {i}: non-positive/non-finite b* {bstar}; oracle undefined as written")
        for arm, b in (("plain", None), ("gmm", bgmm), ("oracle", bstar)):
            df, _ = ms.marker_share_proportions(X, vn, profile, brightness=b)
            df.index = pd.Index(cite.obs_names)
            write_pred(df, root, arm, n, i)
        ms_ = list(bgmm)
        rho_ab = spearmanr([bgmm[m] for m in ms_], [info[m]["mean_gt_abundance"] for m in ms_])[0]
        rho_bs = spearmanr([bgmm[m] for m in ms_], [bstar[m] for m in ms_])[0]
        diag[str(i)] = {"spearman_b_vs_mean_gt_abundance": float(rho_ab), "spearman_b_vs_bstar": float(rho_bs),
                        "markers": {m: {"b_gmm": bgmm[m], "b_star": bstar[m], **info[m]} for m in ms_}}
    dirs = {a: str(root / "preds" / f"{a}_{n}type") for a in ARMS}
    std = {m: str(d) for m, d in S.pred_dirs_for(S.LABELSET, n).items()}
    ev = S.evaluability_table(gtd, types, N_REGIONS, PREFIX)
    masks = S.common_spot_masks(gtd, std, ev, N_REGIONS, PREFIX)
    masks2 = S.common_spot_masks(gtd, {**std, **dirs}, ev, N_REGIONS, PREFIX)
    assert all(masks[r].equals(masks2[r]) for r in range(N_REGIONS)), "adding the arms changed the shared masks"
    res = {}
    for a, d in dirs.items():
        r = S.evaluate_strict(gtd, d, types, masks, ev, N_REGIONS, PREFIX)
        r.pop("region_type_table")
        res[a] = r["summary"]
    # gates
    pub = json.loads((S.XB / "evaluation" / f"results_{S.LABELSET}" / f"method_comparison_singler_{n}type"
                      / f"CITEgeist_v41p_{n}type_results.json").read_text())["summary"]
    for k in ("per_type_mean_r", "overall_mean_pearson_r"):
        assert abs(res["gmm"][k] - pub[k]) < 1e-12, f"GATE FAIL {n}t gmm {k}: {res['gmm'][k]} vs {pub[k]}"
    iso = json.loads((S.XB / "evaluation" / "brightness_isolation" / f"{S.LABELSET}_{n}type.json").read_text())["summary"]["P138"]
    for k in ("per_type_mean_r", "overall_mean_pearson_r"):
        assert abs(res["plain"][k] - iso[k]) < 1e-12, f"GATE FAIL {n}t plain {k}: {res['plain'][k]} vs {iso[k]}"
    assert S.same_preds(dirs["gmm"], S.CG / f"output_marker_share_v41p_{n}type"), "gmm preds != v41p files"
    assert S.same_preds(dirs["plain"], S.CG / f"output_marker_share_p138_{n}type"), "plain preds != P138 files"
    print(f"GATES OK {n}t plain={res['plain']['per_type_mean_r']:.6f} gmm={res['gmm']['per_type_mean_r']:.6f} "
          f"flat gmm={res['gmm']['overall_mean_pearson_r']:.6f}", flush=True)
    g_gmm = res["gmm"]["per_type_mean_r"] - res["plain"]["per_type_mean_r"]
    g_or = res["oracle"]["per_type_mean_r"] - res["plain"]["per_type_mean_r"]
    return {
        "n_types": n, "types": types,
        "arms": {a: {"per_type_mean_r": s["per_type_mean_r"], "flat_r": s["overall_mean_pearson_r"],
                     "per_type_r": {t: s[f"{t}_mean_r"] for t in types}} for a, s in res.items()},
        "G_gmm": g_gmm, "G_or": g_or, "G_or_ge_half_G_gmm": bool(g_or >= 0.5 * g_gmm),
        "region_diagnostics": diag,
    }


if __name__ == "__main__":
    root = Path(sys.argv[1])
    if root.exists():
        sys.exit(f"{root} exists; refuse")
    (subprocess.run([str(S.MAIN / ".commands" / "_refuse_protected_outdir.sh"), str(root)], check=True) if (S.MAIN / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))
    root.mkdir(parents=True)
    payload = {"rule": "docs/superpowers/handoffs/2026-10-25-brightness-oracle-test-rule.md",
               "by_n_types": {str(n): run(n, root) for n in (6, 7)}}
    provenance.write_stamped_json(root / "brightness_oracle_test.json", payload, PRODUCER)
    print("wrote", root / "brightness_oracle_test.json", flush=True)
