"""CITEgeist v41p one-at-a-time settings check, strict per-type r on the rnaprot GT (replaces the QP-era Table S5 stress test).

Production = marker share, ms_v41 dictionary (B Major = CD20 + CD138), plain (brightness OFF) — the manuscript
headline estimator, equal to the P138 arm. The GMM brightness build is evaluated as a one-at-a-time sensitivity.
Each setting changes ONE thing from the build it is measured against. Brightness arms (brightness_off, brightness_gmm)
are measured against plain production (== the manuscript headline); the sub-arms (gmm_max_components_1/2/4,
dictionary_achievable) keep brightness ON and are measured against the brightness-ON ms_v41 3-comp build
(brightness_gmm), so each is a single-axis change within the Gaussian-mixture arm:
  brightness_gmm       on, max 3 components       (= GMM brightness build, the sensitivity; was the QP-era "production")
  brightness            off                       (= production, plain marker share)
  gmm_max_components    1, 2, 4                   (1 = one-component GMM; every marker falls back to its top-decile mean, since 7fc0f2fe)
  dictionary            achievable                (B Major = CD20 + CD45RA; brightness on)
Scored with evaluate_strict on the SAME shared spot masks as results_rnaprot (built from the standard roster).
Gates (abort before writing): production == published plain P138 (1e-12) and == output_marker_share_p138_* files;
brightness_gmm == published CITEgeist_v41p and == output_marker_share_v41p_* files.
Usage: settings_v41p.py OUT_DIR   (must not exist; run from a clean worktree; env MAIN = main repo root for outputs)
Output: OUT_DIR/preds/<setting>_<N>type/Xenium_region_<i>/predictions/baseline.csv.gz and OUT_DIR/settings_v41p_rnaprot.json (stamped)
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
for p in (REPO, REPO / "scripts", XB / "evaluation" / "src", XB):
    sys.path.insert(0, str(p))
import benchmark_constants as C  # noqa: E402
import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from evaluate_benchmark import common_spot_masks, evaluability_table, evaluate_strict  # noqa: E402
from gt_rebuild.score_methods import N_REGIONS, PREFIX, gt_types, pred_dirs_for  # noqa: E402

from CITEgeist.model.deconvolution import marker_share as ms  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/spots_baseline/settings_v41p.py"
LABELSET = "rnaprot"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
MAIN = Path(os.environ.get("MAIN", REPO))
CG = MAIN / "Benchmarking" / "xenium_benchmarking" / "CITEgeist"

# name -> (axis, level, profile, brightness, max_components)
SETTINGS = {
    "production": ("production", "marker share, ms_v41 dictionary, plain (no brightness correction)", "MS_V41", False, 3),
    "brightness_gmm": ("brightness", "on (GMM, max 3 components)", "MS_V41", True, 3),
    "brightness_off": ("brightness", "off", "MS_V41", False, 3),
    "gmm_max_components_1": ("gmm_max_components", "1", "MS_V41", True, 1),
    "gmm_max_components_2": ("gmm_max_components", "2", "MS_V41", True, 2),
    "gmm_max_components_4": ("gmm_max_components", "4", "MS_V41", True, 4),
    "dictionary_achievable": ("dictionary", "ACHIEVABLE (B Major = CD20 + CD45RA)", "ACHIEVABLE", True, 3),
}


def write_setting(n, name, root):
    _, _, prof, bright_on, mc = SETTINGS[name]
    profile = getattr(C, f"{prof}_{n}_CELL_PROFILE_DICT")
    import scanpy as sc

    diag_all = {}
    for i in range(N_REGIONS):
        cite = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
        X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
        vn = list(cite.var_names)
        bright = None
        if bright_on:
            markers = sorted({m for v in profile.values() for m in ms._major_markers(v)})
            bright, diag = ms.estimate_marker_brightness_with_diagnostics(X, vn, markers, max_components=mc)
            diag_all[str(i)] = diag
        df, _ = ms.marker_share_proportions(X, vn, profile, brightness=bright)
        df.index = pd.Index(cite.obs_names)
        d = root / "preds" / f"{name}_{n}type" / f"Xenium_region_{i}" / "predictions"
        d.mkdir(parents=True, exist_ok=True)
        df.to_csv(d / "baseline.csv.gz")
    return diag_all


def same_preds(a, b):
    for i in range(N_REGIONS):
        fa = pd.read_csv(Path(a) / f"Xenium_region_{i}" / "predictions" / "baseline.csv.gz", index_col=0)
        fb = pd.read_csv(Path(b) / f"Xenium_region_{i}" / "predictions" / "baseline.csv.gz", index_col=0)
        if not (fa.shape == fb.shape and (fa.index == fb.index).all() and np.allclose(fa.values, fb.values, atol=1e-12)):
            return False
    return True


def run(n, root):
    gtd = bp.gt_dir(n, LABELSET)
    types = gt_types(gtd, n)
    std = {m: str(d) for m, d in pred_dirs_for(LABELSET, n).items()}
    diag = {name: write_setting(n, name, root) for name in SETTINGS}
    dirs = {name: str(root / "preds" / f"{name}_{n}type") for name in SETTINGS}
    ev = evaluability_table(gtd, types, N_REGIONS, PREFIX)
    masks = common_spot_masks(gtd, std, ev, N_REGIONS, PREFIX)
    masks2 = common_spot_masks(gtd, {**std, **dirs}, ev, N_REGIONS, PREFIX)
    assert all(masks[r].equals(masks2[r]) for r in range(N_REGIONS)), "adding the settings changed the shared masks"
    res = {}
    for name, d in dirs.items():
        r = evaluate_strict(gtd, d, types, masks, ev, N_REGIONS, PREFIX)
        tab = r.pop("region_type_table")
        res[name] = (r["summary"], tab)
    # production == the plain P138 published headline (the SPOTS_P138 arm); the GMM build is a sensitivity.
    pub = json.loads((XB / "evaluation" / f"results_{LABELSET}_p138" / f"method_comparison_singler_{n}type"
                      / f"SPOTS_P138_{n}type_results.json").read_text())["summary"]
    prod = res["production"][0]
    for k in ("per_type_mean_r", "overall_mean_pearson_r"):
        assert abs(prod[k] - pub[k]) < 1e-12, f"GATE FAIL {n}t production {k}: {prod[k]} vs published plain P138 {pub[k]}"
    assert same_preds(dirs["production"], CG / f"output_marker_share_p138_{n}type"), "production preds != P138 files"
    # the GMM brightness build (the old production arm) stays one-at-a-time: gate it against its published v41p result.
    pub_gmm = json.loads((XB / "evaluation" / f"results_{LABELSET}" / f"method_comparison_singler_{n}type"
                          / f"CITEgeist_v41p_{n}type_results.json").read_text())["summary"]
    gmm = res["brightness_gmm"][0]
    for k in ("per_type_mean_r", "overall_mean_pearson_r"):
        assert abs(gmm[k] - pub_gmm[k]) < 1e-12, f"GATE FAIL {n}t brightness_gmm {k}: {gmm[k]} vs published v41p {pub_gmm[k]}"
    assert same_preds(dirs["brightness_gmm"], CG / f"output_marker_share_v41p_{n}type"), "brightness_gmm preds != v41p files"
    assert same_preds(dirs["brightness_off"], CG / f"output_marker_share_p138_{n}type"), "brightness_off preds != P138 files"
    print(f"GATES OK {n}t production(plain P138) ={prod['per_type_mean_r']:.6f} flat={prod['overall_mean_pearson_r']:.6f} "
          f"brightness_gmm={gmm['per_type_mean_r']:.6f}", flush=True)
    # Two reference axes so each delta is a clean one-at-a-time change. The brightness arms measure against plain
    # production (the headline); the sub-arms keep brightness ON and measure against the brightness-ON ms_v41 3-comp
    # build (`gmm`), so "instead of 3" / "the CD20+CD45RA dictionary for B cells" read as single-axis changes inside
    # the Gaussian-mixture arm. Delta to `prod` for a sub-arm would be a two-axis change (brightness AND the sub-axis).
    refmap = {
        "production": prod, "brightness_off": prod, "brightness_gmm": prod,
        "gmm_max_components_1": gmm, "gmm_max_components_2": gmm, "gmm_max_components_4": gmm,
        "dictionary_achievable": gmm,
    }
    out = {}
    for name, (s, tab) in res.items():
        axis, level = SETTINGS[name][:2]
        ref = refmap[name]
        per_type = {t: s[f"{t}_mean_r"] for t in types}
        out[name] = {
            "axis": axis, "level": level,
            "per_type_mean_r": s["per_type_mean_r"], "flat_r": s["overall_mean_pearson_r"],
            "per_type_r": per_type,
            "delta_per_type_mean_r": s["per_type_mean_r"] - ref["per_type_mean_r"],
            "delta_flat_r": s["overall_mean_pearson_r"] - ref["overall_mean_pearson_r"],
            "delta_per_type_r": {t: per_type[t] - ref[f"{t}_mean_r"] for t in types},
            "region_type_r": [{"region": int(a), "type": b, "r": float(c)} for a, b, c in
                              zip(tab["region"], tab["type"], tab["r"])],
            "brightness_diagnostics": diag[name],
        }
    return {"n_types": n, "types": types, "dirs": dirs, "settings": out}


if __name__ == "__main__":
    root = Path(sys.argv[1])
    if root.exists():
        sys.exit(f"{root} exists; refuse")
    import subprocess
    (subprocess.run([str(MAIN / ".commands" / "_refuse_protected_outdir.sh"), str(root)], check=True) if (MAIN / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))
    root.mkdir(parents=True)
    payload = {"labelset": LABELSET,
               "production": "marker share, ms_v41 dictionary, plain (no brightness correction)",
               "by_n_types": {str(n): run(n, root) for n in (6, 7)}}
    provenance.write_stamped_json(root / "settings_v41p_rnaprot.json", payload, PRODUCER)
    print("wrote", root / "settings_v41p_rnaprot.json", flush=True)
