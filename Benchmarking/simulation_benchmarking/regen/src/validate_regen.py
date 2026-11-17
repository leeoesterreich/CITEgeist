"""A10 validation: score regenerated predictions with the repo's own evaluators and compare to the cached
2026-04-28-generation means (the cached simulated_benchmark/{prop,gex}_comparison.json).

Evaluators (unchanged): benchmarking_spot_deconv.main (proportions) and benchmarking_gex.calculate_rmse(normalize="range").
Hard asserts per prediction: 9/9 cell types, all GT spots present (no silent intersection), 9/9 GEX layers.
Tolerances (binding): prop |d corr| <= 0.01 ; GEX flat_r |d| <= 0.02. Never loosened here.

Usage: python validate_regen.py --methods RCTD CARD Seurat C2L Tangram [--results-root DIR] [--out DIR]
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import contextlib
import io
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

B = Path(str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking"))
sys.path.insert(0, str(B / "src"))
from benchmarking_gex import calculate_rmse  # noqa: E402
from benchmarking_spot_deconv import main as spot_deconv_main  # noqa: E402

CACHE = Path(str(_cg.root("DATA") / "simulated_benchmark_cache"))
PROP_TOL, GEX_TOL = 0.01, 0.02
CELL_TYPES = sorted(
    ["B-cells", "CAFs", "Cancer Epithelial", "Endothelial", "Myeloid", "Normal Epithelial", "PVL", "Plasmablasts", "T-cells"]
)
CACHE_KEY = {"RCTD": "RCTD", "CARD": "CARD", "Seurat": "Seurat", "C2L": "Cell2Location", "Tangram": "Tangram"}
PRED = {
    "RCTD": "{root}/rctd/{cond}/Wu_rep_{rep}/Wu_ST_{rep}_RCTD_deconv_predictions.csv",
    "CARD": "{root}/card/{cond}/Wu_rep_{rep}/Wu_ST_{rep}_CARD_deconv_predictions.csv",
    "Seurat": "{root}/seurat/{cond}/Wu_rep_{rep}/Wu_rep_{rep}_Seurat_deconv_predictions.csv",
    "C2L": "{root}/c2l/{cond}/Wu_rep_{rep}/cell2loc_deconv_predictions.csv",
    "Tangram": "{root}/tangram/{cond}/Wu_rep_{rep}/cell_type_proportions.csv",
}
LAYERS = {"C2L": "{root}/c2l/{cond}/Wu_rep_{rep}/layers", "Tangram": "{root}/tangram/{cond}/Wu_rep_{rep}/layers"}


def prop_one(pred_file, gt_file):
    gt = pd.read_csv(gt_file, index_col=0).drop(columns=["spot_x", "spot_y"], errors="ignore")
    pr = pd.read_csv(pred_file, index_col=0).drop(columns=["spot_x", "spot_y", "recon_error"], errors="ignore")
    assert sorted(gt.columns) == CELL_TYPES, f"GT types {sorted(gt.columns)}"
    assert sorted(pr.columns) == CELL_TYPES, f"{pred_file}: pred types {sorted(pr.columns)} (need 9/9)"
    missing = gt.index.difference(pr.index)
    assert len(missing) == 0, f"{pred_file}: {len(missing)} GT spots missing in prediction"
    pr = pr.loc[gt.index]
    with tempfile.TemporaryDirectory(dir="/tmp") as d:
        g, p = os.path.join(d, "gt.csv"), os.path.join(d, "pred.csv")
        gt.to_csv(g)
        pr.to_csv(p)
        with contextlib.redirect_stdout(io.StringIO()):
            r = spot_deconv_main(g, p)
    return {"n_spots": len(gt), "corr": float(r["corr"]), "jsd": float(r["JSD"]), "rmse": float(r["Sum_RMSE"]),
            "mae": float(r["Sum_MAE"])}


def gex_one(layers_dir, gt_dir, gene_subset=None):
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        r = calculate_rmse(str(gt_dir), str(layers_dir), normalize="range", gene_subset=gene_subset)
    n_pred = len([f for f in os.listdir(layers_dir) if f.endswith("_layer.csv")])
    assert n_pred == 9, f"{layers_dir}: {n_pred}/9 prediction layers present"
    # calculate_rmse silently skips types whose GT layer is all-zero; allow that ONLY when verified.
    scored = {k.replace(" ", "_") for k in r["metrics_per_cell_type"]}
    for gtf in sorted(Path(gt_dir).glob("*_GT.csv")):
        ct = gtf.name[: -len("_GT.csv")]
        if ct.replace(" ", "_") not in scored:
            tot = pd.read_csv(gtf, index_col=0).values.sum()
            assert tot == 0, f"{layers_dir}: {ct} unscored although GT total={tot}"
    return {"flat_r": float(r["mean_pearson_r"]), "spot_layer_r": float(r["mean_spot_layer_r"]),
            "avg_rmse": float(r["average_rmse"]), "avg_nrmse": float(r["average_nrmse"]), "avg_mae": float(r["average_mae"]),
            "n_types_scored": len(scored),
            "per_type": {k: {m: float(v) for m, v in d.items()} for k, d in r["metrics_per_cell_type"].items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods", nargs="+", default=list(PRED))
    ap.add_argument("--results-root", default=str(B / "regen/results"))
    ap.add_argument("--out", default=str(B / "regen/validation"))
    ap.add_argument("--gex-variant", default="layers", help="layers dir name (e.g. layers_legacy for the Tangram diagnostic)")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    cache_p = json.load(open(CACHE / "prop_comparison.json"))
    cache_g = json.load(open(CACHE / "gex_comparison.json"))
    rows, summ = [], []
    for m in a.methods:
        for cond in ["high_seg", "mixed"]:
            pr_rows, gx_rows = [], []
            for rep in range(5):
                gt_prop = B / f"replicates/{cond}/ST_sim/Wu_ST_{rep}_prop.csv"
                pf = Path(PRED[m].format(root=a.results_root, cond=cond, rep=rep))
                if pf.exists():
                    r = prop_one(pf, gt_prop)
                    pr_rows.append(r)
                    rows.append({"method": m, "cond": cond, "rep": rep, "kind": "prop", **r})
                if m in LAYERS:
                    ld = Path(LAYERS[m].format(root=a.results_root, cond=cond, rep=rep)).parent / a.gex_variant
                    if ld.is_dir():
                        r = gex_one(ld, B / f"replicates/{cond}/ST_GEX_sim/sample_{rep}/layers")
                        gx_rows.append(r)
                        rows.append({"method": m, "cond": cond, "rep": rep, "kind": "gex", **r})
            ck = f"{CACHE_KEY[m]}_{cond}"
            if pr_rows:
                new = {k: float(np.mean([x[k] for x in pr_rows])) for k in ["corr", "jsd", "rmse"]}
                old = cache_p[ck]
                d = new["corr"] - old["r"]
                summ.append({"method": m, "cond": cond, "kind": "prop", "n_reps": len(pr_rows),
                             "new_corr": new["corr"], "cached_r": old["r"], "delta_corr": d, "tol": PROP_TOL,
                             "pass": bool(abs(d) <= PROP_TOL), "new_jsd": new["jsd"], "cached_jsd": old["jsd"],
                             "new_rmse": new["rmse"], "cached_rmse": old["rmse"]})
            if gx_rows and ck in cache_g:
                new = {k: float(np.mean([x[k] for x in gx_rows])) for k in gx_rows[0]}
                old = cache_g[ck]
                d = new["flat_r"] - old["Mean_Pearson_r"]
                summ.append({"method": m, "cond": cond, "kind": "gex", "n_reps": len(gx_rows),
                             "new_flat_r": new["flat_r"], "cached_flat_r": old["Mean_Pearson_r"], "delta_flat_r": d,
                             "tol": GEX_TOL, "pass": bool(abs(d) <= GEX_TOL), "new_nrmse": new["avg_nrmse"],
                             "cached_nrmse": old["Avg_NRMSE"], "new_rmse": new["avg_rmse"], "cached_rmse": old["Avg_RMSE"],
                             "new_spot_layer_r": new["spot_layer_r"], "cached_spot_layer_r": old["Mean_Spot_Layer_r"]})
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = a.tag or "_".join(a.methods)
    pd.DataFrame(rows).to_csv(out / f"per_rep_{tag}.csv", index=False)
    (out / f"summary_{tag}.json").write_text(json.dumps(summ, indent=2))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(pd.DataFrame(summ).to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("ALL PASS" if summ and all(s["pass"] for s in summ) else "SOME FAIL / EMPTY")


if __name__ == "__main__":
    main()
