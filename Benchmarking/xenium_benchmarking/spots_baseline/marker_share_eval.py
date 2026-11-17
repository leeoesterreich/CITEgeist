"""Does the GMM brightness correction beat the plain marker-share estimator? Xenium RCC (6- and 7-type).

Mirrors spots_baseline.py main() exactly (same CITEgeist ACHIEVABLE_{6,7} dictionaries, Major markers, SingleR GT dirs,
evaluate_all_regions call and settings); only the estimator is the package function
CITEgeist.model.deconvolution.marker_share.marker_share_proportions. Two arms:
  plain     brightness=None (must equal the SPOTS reading-a baseline: asserted below, exit non-zero otherwise)
  corrected brightness estimated per region from THAT region's CITE matrix only, over that dictionary's Major markers
            (estimate_marker_brightness, default seed); stored per region in _inputs.brightness_per_region.
No setting is tuned on the ground truth.

Usage: marker_share_eval.py N_TYPES OUT_DIR   (N_TYPES in 6, 7)
Writes OUT_DIR/{plain,corrected}_{N}type_results.json
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
SPOTS_REF = XB / "evaluation" / "results" / "spots_baseline"
ARMS = ("plain", "corrected")
PRODUCED_BY = "Benchmarking/xenium_benchmarking/spots_baseline/marker_share_eval.py"


def main(n, out_dir):
    import scanpy as sc

    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(XB))
    sys.path.insert(0, str(XB / "evaluation" / "src"))
    sys.path.insert(0, str(REPO / "scripts"))
    from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT, ACHIEVABLE_7_CELL_PROFILE_DICT
    from evaluate_benchmark import evaluate_all_regions
    from provenance import write_stamped_json

    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness, marker_share_proportions

    profile = {6: ACHIEVABLE_6_CELL_PROFILE_DICT, 7: ACHIEVABLE_7_CELL_PROFILE_DICT}[n]
    markers = sorted({m for v in profile.values() for m in v["Major"]})
    gt = XB / "ground_truth_singler" / f"ground_truth_{n}type"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cites = {}
    for i in range(5):
        cite = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
        X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
        cites[i] = (X, list(cite.var_names), cite.obs_names)
    bright = {str(i): estimate_marker_brightness(X, vn, markers) for i, (X, vn, _) in cites.items()}
    results = {}
    for arm in ARMS:
        stage = Path(tempfile.mkdtemp(prefix=f"ms_{arm}_{n}_"))
        zeros = {}
        for i, (X, vn, obs) in cites.items():
            b = bright[str(i)] if arm == "corrected" else None
            df, nz = marker_share_proportions(X, vn, profile, brightness=b)
            df.index = obs
            df.to_csv(stage / f"Xenium_region_{i}_cell_prop_global_results.csv")
            zeros[str(i)] = {"n_spots": int(len(df)), "n_all_zero_marker_spots": nz}
        combine = n == 6
        res = evaluate_all_regions(
            gt_dir=str(gt), pred_dir=str(stage), n_regions=5, output_path=None, prefix="Xenium",
            use_achievable_7=not combine, combine_t_cells=combine,
        )
        assert res["summary"]["n_regions"] == 5
        res["_inputs"] = {
            "estimator": "marker_share_proportions (CITEgeist.model.deconvolution.marker_share)",
            "arm": arm,
            "brightness": "none" if arm == "plain" else "per-region GMM positive-mode mean of cosine-normalized channel",
            "brightness_per_region": bright if arm == "corrected" else None,
            "marker_dict": f"benchmark_constants.ACHIEVABLE_{n}_CELL_PROFILE_DICT (Major markers)",
            "all_zero_marker_spots": zeros,
            "evaluator": "evaluate_benchmark.evaluate_all_regions (same settings as spots_baseline.py)",
            "gt_dir": str(gt.relative_to(REPO)),
            "cite_dir": str(DATA.relative_to(REPO)),
        }
        results[arm] = res
        shutil.rmtree(stage)
        s = res["summary"]
        print(f"{n}-type {arm}: r={s['overall_mean_pearson_r']:.6f} rmse={s['overall_mean_rmse']:.4f} "
              f"mae={s['overall_mean_mae']:.4f}", flush=True)
    ref = json.loads((SPOTS_REF / f"SPOTS_{n}type_results.json").read_text())["summary"]["overall_mean_pearson_r"]
    got = results["plain"]["summary"]["overall_mean_pearson_r"]
    if not abs(got - ref) < 1e-9:
        print(f"PARITY FAILED {n}-type: plain r={got!r} vs SPOTS reading-a r={ref!r}", flush=True)
        sys.exit(1)
    print(f"PARITY OK {n}-type: plain r={got:.12f} == SPOTS reading-a r={ref:.12f} (|diff|={abs(got - ref):.2e})", flush=True)
    for arm in ARMS:
        write_stamped_json(out_dir / f"{arm}_{n}type_results.json", results[arm], PRODUCED_BY)
    print("wrote", out_dir)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2])
