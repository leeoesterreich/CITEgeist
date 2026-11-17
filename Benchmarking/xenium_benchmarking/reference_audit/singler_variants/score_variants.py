"""Rebuild GT from variant SingleR labels (same code path as aggregate_to_spots.py) and re-score existing 6-type predictions
with the canonical evaluator (evaluate_benchmark.evaluate_all_regions, compare_all_methods.py settings). Read-only on canonical files."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import json, shutil, sys, glob
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(str(_cg.REPO_ROOT))
XB = REPO / "Benchmarking" / "xenium_benchmarking"
SV = XB / "reference_audit" / "singler_variants"
MAP = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_to_spot_mapping.csv"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
PS = REPO / "Benchmarking" / "protein_sensitivity" / "results" / "baseline_preds"
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(XB)); sys.path.insert(0, str(XB / "evaluation" / "src"))
from evaluate_benchmark import evaluate_all_regions
CT6 = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]

canon = pd.read_csv(XB / "ground_truth_singler/singler_labels/singler_all_cells.csv").set_index("cell_id")
mapping = pd.read_csv(MAP, index_col=0)
var = {n: pd.read_csv(SV / "singler_out" / n / "singler_all_cells.csv").set_index("cell_id") for n in ["control", "A", "B"]}

# ---- 1. confusion tables
conf = {}
for n in ["control", "A", "B"]:
    assert (var[n].index == canon.index).all()
    ct = pd.crosstab(canon.singler_label, var[n].singler_label)
    conf[n] = ct; ct.to_csv(SV / f"confusion_canonical_vs_{n}.csv")
    print(f"\n== confusion canonical (rows) vs {n} (cols) ==\n{ct.to_string()}")
    print("total cells", len(canon), " agreement (label identical):", int((canon.singler_label == var[n].singler_label).sum()))

# ---- GT builder (identical logic to aggregate_to_spots.py, cell types fixed)
def build_gt(labels: pd.Series, outdir: Path):
    merged = mapping.join(labels.rename("singler_label"), how="inner")
    mapped = merged[merged["spot_idx"] >= 0]
    mapped = mapped[mapped["singler_label"].isin(CT6)]
    (outdir / "ground_truth").mkdir(parents=True, exist_ok=True)
    cellfrac = {}
    for r in range(5):
        rc = mapped[mapped["region_id"] == r]
        c = rc.groupby(["spot_id", "singler_label"]).size().unstack(fill_value=0)
        for ct in CT6:
            if ct not in c.columns: c[ct] = 0
        c = c[CT6]
        props = c.div(c.sum(axis=1), axis=0).fillna(0)
        props = props.join(rc.groupby("spot_id").agg(n_cells=("singler_label", "count"), spot_x=("spot_x", "first"), spot_y=("spot_y", "first")))
        props.to_csv(outdir / "ground_truth" / f"Xenium_region_{r}_prop.csv")
        cellfrac[r] = dict(n_cells=len(rc), T_cell_frac=float((rc.singler_label == "T cells").mean()),
                           n_spots=len(props), T_spot_mean=float(props["T cells"].mean()))
    return cellfrac

# ---- GT variants
lab = {"recomputed_canonical": canon.singler_label,
       "control": var["control"].singler_label,
       "A_NKmerged": var["A"].singler_label.replace({"NK cells": "T cells"}),
       "A_NKremoved": var["A"].singler_label,       # NK not in CT6 -> dropped before proportions
       "B": var["B"].singler_label}
GT = {}; frac = {}
for k, s in lab.items():
    GT[k] = SV / "gt" / k
    frac[k] = build_gt(s, GT[k])
GT["canonical"] = XB / "ground_truth_singler" / "ground_truth_6type"
# GT recompute parity
for r in range(5):
    a = pd.read_csv(GT["canonical"] / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_col=0)
    b = pd.read_csv(GT["recomputed_canonical"] / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_col=0)
    pd.testing.assert_frame_equal(a.sort_index(), b.sort_index(), check_exact=False, rtol=1e-9, atol=1e-12)
print("\nGT rebuild reproduces canonical GT files: OK")
# NK-merged-in-T proportion region table (cells), plus NK fraction
rows = []
for r in range(5):
    row = {"region": r}
    for k in ["recomputed_canonical", "control", "A_NKremoved", "A_NKmerged", "B"]:
        row[f"{k}_Tcellfrac"] = frac[k][r]["T_cell_frac"]; row[f"{k}_Tspotmean"] = frac[k][r]["T_spot_mean"]
        row[f"{k}_ncells"] = frac[k][r]["n_cells"]; row[f"{k}_nspots"] = frac[k][r]["n_spots"]
    m = mapping[(mapping.spot_idx >= 0) & (mapping.region_id == r)].join(var["A"].singler_label)
    row["A_NK_cellfrac_of_6type+NK"] = float((m.singler_label == "NK cells").sum() / m.singler_label.isin(CT6 + ["NK cells"]).sum())
    rows.append(row)
pd.DataFrame(rows).to_csv(SV / "per_region_T_proportions.csv", index=False)
print(pd.DataFrame(rows).round(4).T.to_string())

# ---- predictions
stage = SV / "pred_stage"; shutil.rmtree(stage, ignore_errors=True)
PRED = {m: XB / m / "output_singler_6type" for m in ["Cell2Location", "RCTD", "Tangram", "Seurat", "CARD"]}
d = stage / "CITEgeist"; d.mkdir(parents=True)
for i in range(5):
    run = [x for x in glob.glob(str(PS / f"xenium_6_{i}_race_*")) if (Path(x) / "predictions" / "baseline.csv.gz").exists() and (Path(x) / "run_meta.json").exists()]
    assert len(run) == 1, run
    pd.read_csv(Path(run[0]) / "predictions" / "baseline.csv.gz", index_col=0).to_csv(d / f"Xenium_region_{i}_cell_prop_global_results.csv")
PRED["CITEgeist"] = d
# marker share (recomputed: no stored prediction files). plain = reported 0.789 arm; corrected also.
import scanpy as sc
from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT
from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness, marker_share_proportions
markers = sorted({m for v in ACHIEVABLE_6_CELL_PROFILE_DICT.values() for m in v["Major"]})
cites = {}
for i in range(5):
    c = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
    X = c.X.toarray() if hasattr(c.X, "toarray") else np.asarray(c.X); cites[i] = (X, list(c.var_names), c.obs_names)
bright = {i: estimate_marker_brightness(X, vn, markers) for i, (X, vn, _) in cites.items()}
for arm in ["plain", "corrected"]:
    d = stage / f"MarkerShare_{arm}"; d.mkdir()
    for i, (X, vn, obs) in cites.items():
        df, _ = marker_share_proportions(X, vn, ACHIEVABLE_6_CELL_PROFILE_DICT, brightness=bright[i] if arm == "corrected" else None)
        df.index = obs; df.to_csv(d / f"Xenium_region_{i}_cell_prop_global_results.csv")
    PRED[f"MarkerShare_{arm}"] = d

def ev(gt, pd_):
    return evaluate_all_regions(gt_dir=str(gt), pred_dir=str(pd_), n_regions=5, output_path=None, prefix="Xenium",
                                use_achievable_7=False, combine_t_cells=True)

# canonical reference numbers
ref = {}
for m in ["CITEgeist", "Cell2Location", "RCTD", "Tangram", "Seurat", "CARD"]:
    ref[m] = json.load(open(XB / "evaluation/results/method_comparison_singler_6type" / (f"{m}_6type_results.json" if m != "CARD" else "CARD_results.json")))["summary"]["overall_mean_pearson_r"]
ref["MarkerShare_plain"] = json.load(open(XB / "evaluation/results/marker_share/plain_6type_results.json"))["summary"]["overall_mean_pearson_r"]
ref["MarkerShare_corrected"] = json.load(open(XB / "evaluation/results/marker_share/corrected_6type_results.json"))["summary"]["overall_mean_pearson_r"]

out = []; spots = {}
for gtn in ["canonical", "recomputed_canonical", "control", "A_NKmerged", "A_NKremoved", "B"]:
    for m, pdir in PRED.items():
        res = ev(GT[gtn], pdir); s = res["summary"]
        out.append(dict(gt=gtn, method=m, pooled_r=s["overall_mean_pearson_r"], std_r=s["overall_std_pearson_r"], rmse=s["overall_mean_rmse"],
                        mae=s["overall_mean_mae"], jsd=s.get("overall_mean_jsd", np.nan), T_r=s["T cells_mean_r"], total_spots=s["total_spots"],
                        **{f"{c}_r": s[f"{c}_mean_r"] for c in CT6 if c != "T cells"}))
df = pd.DataFrame(out); df.to_csv(SV / "scores_all.csv", index=False)
# sanity
print("\n== SANITY: canonical GT file rescore vs canonical JSON ==")
bad = False
for m in PRED:
    got = df[(df["gt"] == "canonical") & (df["method"] == m)].pooled_r.iloc[0]
    ok = abs(got - ref[m]) < 1e-9; bad |= not ok
    print(f"{m:24s} rescored {got:.6f} canonical {ref[m]:.6f} {'OK' if ok else 'MISMATCH'}")
print("SANITY", "FAILED" if bad else "PASSED")
pd.set_option("display.width", 250)
print("\n== pooled r ==\n", df.pivot(index="method", columns="gt", values="pooled_r").round(4).to_string())
print("\n== T cell r ==\n", df.pivot(index="method", columns="gt", values="T_r").round(4).to_string())
print("\n== rmse ==\n", df.pivot(index="method", columns="gt", values="rmse").round(4).to_string())
print("\n== total spots ==\n", df.pivot(index="method", columns="gt", values="total_spots").to_string())
print("\n== ranking by pooled r ==")
for g in ["canonical", "control", "A_NKmerged", "A_NKremoved", "B"]:
    sub = df[df["gt"] == g].sort_values("pooled_r", ascending=False)
    print(g, " > ".join(f"{m}({r:.3f})" for m, r in zip(sub.method, sub.pooled_r)))
shutil.rmtree(stage)
if bad: sys.exit(1)
