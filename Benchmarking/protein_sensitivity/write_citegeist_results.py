"""Write the stamped CITEgeist Xenium results files from the saved baseline predictions (A3 step 2).

Input : results/baseline_preds/xenium_{N}_{region}_race_*/ (predictions/baseline.csv.gz + metrics.csv),
        produced by `OUT_NAME=baseline_preds .commands/submit-protein-sensitivity.sh gate xenium {6|7}`.
Output: evaluation/results/method_comparison_singler_{N}type/CITEgeist_{N}type_results.json (stamped).

Settings mirror compare_all_methods.py. Hard asserts (never loosened): n_regions == 5, and each region's
Pearson r equals that run's own metrics.csv within 1e-4. Drift vs the old (pre-prediction-saving) gate
run is reported only.
Usage: write_citegeist_results.py N_TYPES
"""
import glob
import shutil
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
XEN = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XEN / "evaluation" / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from evaluate_benchmark import evaluate_all_regions  # noqa: E402
from provenance import write_stamped_json  # noqa: E402

TOL = 1e-4
PS = REPO / "Benchmarking" / "protein_sensitivity" / "results"


def _run_dir(root: Path, n: int, region: int) -> Path:
    """The single completed run dir (baseline.csv.gz + metrics.csv + run_meta.json) for a region."""
    cands = [
        d
        for d in glob.glob(str(root / f"xenium_{n}_{region}_race_*"))
        if (Path(d) / "predictions" / "baseline.csv.gz").exists() and (Path(d) / "run_meta.json").exists()
    ]
    assert len(cands) == 1, f"region {region}: expected exactly 1 completed run dir, got {cands}"
    return Path(cands[0])


def main(n: int) -> None:
    assert n in (6, 7)
    root = PS / "baseline_preds"
    stage = root / f"eval_{n}type"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir()
    runs = {}
    for i in range(5):
        d = _run_dir(root, n, i)
        shutil.copy(d / "predictions" / "baseline.csv.gz", stage / "tmp.csv.gz")
        pd.read_csv(stage / "tmp.csv.gz", index_col=0).to_csv(stage / f"Xenium_region_{i}_cell_prop_global_results.csv")
        (stage / "tmp.csv.gz").unlink()
        runs[i] = d
    gt = XEN / "ground_truth_singler" / f"ground_truth_{n}type"
    combine = n == 6
    res = evaluate_all_regions(
        gt_dir=str(gt),
        pred_dir=str(stage),
        n_regions=5,
        output_path=None,
        prefix="Xenium",
        use_achievable_7=not combine,
        combine_t_cells=combine,
    )
    assert res["summary"]["n_regions"] == 5, res["summary"]["n_regions"]
    print(f"{'region':>6} {'eval_r':>9} {'own_metrics_r':>14} {'abs_diff':>10} {'old_gate_r':>11} {'drift_vs_old':>13}")
    for i, reg in enumerate(res["regions"]):
        assert reg["region_id"] == i
        own = pd.read_csv(runs[i] / "metrics.csv")
        own = own[own["axis"] == "baseline"]
        assert len(own) == 1 and own["status"].iloc[0] == "ok", f"region {i}: bad metrics.csv"
        own_r = float(own["pearson_r"].iloc[0])
        diff = abs(reg["overall_pearson_r"] - own_r)
        old = glob.glob(str(PS / "gate" / f"xenium_{n}_{i}_race_*" / "metrics.csv"))
        old_r = float(pd.read_csv(old[0])["pearson_r"].iloc[0]) if old else float("nan")
        print(f"{i:>6} {reg['overall_pearson_r']:>9.5f} {own_r:>14.5f} {diff:>10.2e} {old_r:>11.5f} {own_r - old_r:>13.2e}")
        assert diff <= TOL, f"region {i}: evaluator r != run's own metrics.csv r (diff {diff}); STOP"
    s = res["summary"]
    print(f"{n}-type summary: r={s['overall_mean_pearson_r']:.4f} rmse={s['overall_mean_rmse']:.4f} mae={s['overall_mean_mae']:.4f}")
    res["_inputs"] = {
        "run_dirs": {str(i): str(runs[i].relative_to(REPO)) for i in runs},
        "evaluator": "evaluate_benchmark.evaluate_all_regions (compare_all_methods.py settings)",
        "gt_dir": str(gt.relative_to(REPO)),
    }
    out = XEN / "evaluation" / "results" / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json"
    write_stamped_json(out, res, "Benchmarking/protein_sensitivity/write_citegeist_results.py")
    print("wrote", out)


if __name__ == "__main__":
    main(int(sys.argv[1]))
