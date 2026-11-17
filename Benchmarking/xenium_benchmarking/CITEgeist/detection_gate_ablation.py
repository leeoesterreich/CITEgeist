#!/usr/bin/env python -u
"""Detection gating ablation — compare B cell over-prediction with/without gate.

Runs QP deconvolution on all 5 Xenium RCC regions with detection gating
enabled (production) and disabled, then computes fold over-prediction
for B cells relative to ground truth.

Output: detection_gate_ablation_results.json
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

BENCHMARK_ROOT = REPO_ROOT / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

from CITEgeist.model.citegeist_model import CitegeistModel

DATA_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"
GT_DIR = BENCHMARK_ROOT / "ground_truth_singler" / "ground_truth_6type" / "ground_truth"
OUTPUT_DIR = Path(__file__).parent / "output_gate_ablation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CELL_TYPES_6 = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]

# Merge CD4+/CD8+ T cells into "T cells" for 6-type comparison
PROFILE_6TYPE = {}
for k, v in ACHIEVABLE_7_CELL_PROFILE_DICT.items():
    if k in ("CD4+ T cells", "CD8+ T cells"):
        if "T cells" not in PROFILE_6TYPE:
            PROFILE_6TYPE["T cells"] = {"Major": ["CD3E"], "Minor": ["CD4", "CD8A", "CD45RO", "GranzymeB"]}
    else:
        PROFILE_6TYPE[k] = v


def load_region(region_id: int):
    h5ad_dir = DATA_DIR / "h5ad_objects"
    gex = sc.read_h5ad(h5ad_dir / f"Xenium_region_{region_id}_GEX.h5ad")
    cite = sc.read_h5ad(h5ad_dir / f"Xenium_region_{region_id}_CITE.h5ad")
    return gex, cite


def load_gt_6type(region_id: int) -> pd.DataFrame:
    gt_path = GT_DIR / f"Xenium_region_{region_id}_prop.csv"
    return pd.read_csv(gt_path, index_col=0)


def run_qp(region_id: int, use_detection_gating: bool) -> pd.DataFrame:
    """Run QP deconvolution and return predicted proportions."""
    gex, cite = load_region(region_id)
    tag = "gate_on" if use_detection_gating else "gate_off"
    sample_name = f"ablation_{tag}_region_{region_id}"

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(OUTPUT_DIR / tag),
        simulation=True,
        gene_expression_adata=gex,
        antibody_capture_adata=cite,
    )
    model.load_cell_profile_dict(PROFILE_6TYPE)
    model.filter_gex(min_counts=25)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    model.run_cell_proportion_model(
        method="qp",
        use_detection_gating=use_detection_gating,
        use_gex_detection=use_detection_gating,
        detection_gate_ub=0.01,
        lambda_reg=1.0,
        validation_warn_only=True,
    )

    prop_df = model.results["cell_prop"]
    # Merge T cell subtypes if present
    if "CD4+ T cells" in prop_df.columns and "CD8+ T cells" in prop_df.columns:
        prop_df["T cells"] = prop_df["CD4+ T cells"] + prop_df["CD8+ T cells"]
        prop_df = prop_df.drop(columns=["CD4+ T cells", "CD8+ T cells"])
    return prop_df


def main():
    results = {"regions": {}, "summary": {}}

    all_bcell_gt = []
    all_bcell_gate_on = []
    all_bcell_gate_off = []

    for region_id in range(5):
        print(f"\n{'='*60}")
        print(f"Region {region_id}")
        print(f"{'='*60}")

        gt_df = load_gt_6type(region_id)

        # With detection gating (production)
        print(f"  Running WITH detection gating...")
        t0 = time.time()
        pred_on = run_qp(region_id, use_detection_gating=True)
        t_on = time.time() - t0
        print(f"  Done in {t_on:.1f}s")

        # Without detection gating
        print(f"  Running WITHOUT detection gating...")
        t0 = time.time()
        pred_off = run_qp(region_id, use_detection_gating=False)
        t_off = time.time() - t0
        print(f"  Done in {t_off:.1f}s")

        # Align indices
        common = sorted(set(gt_df.index) & set(pred_on.index) & set(pred_off.index))
        gt_aligned = gt_df.loc[common, CELL_TYPES_6]
        # Normalize GT to proportions if needed
        row_sums = gt_aligned.sum(axis=1)
        if row_sums.max() > 1.5:
            gt_aligned = gt_aligned.div(row_sums, axis=0)

        pred_on_aligned = pred_on.loc[common, [c for c in CELL_TYPES_6 if c in pred_on.columns]]
        pred_off_aligned = pred_off.loc[common, [c for c in CELL_TYPES_6 if c in pred_off.columns]]

        # B cell statistics
        bcell_gt = gt_aligned["B cells"].values
        bcell_on = pred_on_aligned["B cells"].values
        bcell_off = pred_off_aligned["B cells"].values

        all_bcell_gt.extend(bcell_gt)
        all_bcell_gate_on.extend(bcell_on)
        all_bcell_gate_off.extend(bcell_off)

        # Per-region metrics
        r_on, _ = stats.pearsonr(gt_aligned.values.flatten(), pred_on_aligned.values.flatten())
        r_off, _ = stats.pearsonr(gt_aligned.values.flatten(), pred_off_aligned.values.flatten())
        r_bcell_on, _ = stats.pearsonr(bcell_gt, bcell_on)
        r_bcell_off, _ = stats.pearsonr(bcell_gt, bcell_off)

        region_result = {
            "n_spots": len(common),
            "overall_r_with_gate": round(float(r_on), 4),
            "overall_r_without_gate": round(float(r_off), 4),
            "bcell_r_with_gate": round(float(r_bcell_on), 4),
            "bcell_r_without_gate": round(float(r_bcell_off), 4),
            "bcell_mean_gt": round(float(bcell_gt.mean()), 5),
            "bcell_mean_pred_with_gate": round(float(bcell_on.mean()), 5),
            "bcell_mean_pred_without_gate": round(float(bcell_off.mean()), 5),
        }
        results["regions"][f"region_{region_id}"] = region_result
        print(f"  Overall r: {r_on:.3f} (gate on) vs {r_off:.3f} (gate off)")
        print(f"  B cell r: {r_bcell_on:.3f} (gate on) vs {r_bcell_off:.3f} (gate off)")
        print(f"  B cell mean: GT={bcell_gt.mean():.4f}, pred_on={bcell_on.mean():.4f}, pred_off={bcell_off.mean():.4f}")

    # Summary across all regions
    all_bcell_gt = np.array(all_bcell_gt)
    all_bcell_gate_on = np.array(all_bcell_gate_on)
    all_bcell_gate_off = np.array(all_bcell_gate_off)

    mean_gt = float(all_bcell_gt.mean())
    mean_on = float(all_bcell_gate_on.mean())
    mean_off = float(all_bcell_gate_off.mean())

    fold_overpred_without_gate = mean_off / mean_gt if mean_gt > 0 else float("nan")
    fold_overpred_with_gate = mean_on / mean_gt if mean_gt > 0 else float("nan")
    fold_reduction = fold_overpred_without_gate / fold_overpred_with_gate if fold_overpred_with_gate > 0 else float("nan")

    r_all_on, _ = stats.pearsonr(all_bcell_gt, all_bcell_gate_on)
    r_all_off, _ = stats.pearsonr(all_bcell_gt, all_bcell_gate_off)

    results["summary"] = {
        "total_spots": len(all_bcell_gt),
        "bcell_mean_gt": round(mean_gt, 5),
        "bcell_mean_pred_with_gate": round(mean_on, 5),
        "bcell_mean_pred_without_gate": round(mean_off, 5),
        "fold_overprediction_without_gate": round(fold_overpred_without_gate, 2),
        "fold_overprediction_with_gate": round(fold_overpred_with_gate, 2),
        "fold_reduction_by_gating": round(fold_reduction, 2),
        "bcell_r_with_gate": round(float(r_all_on), 4),
        "bcell_r_without_gate": round(float(r_all_off), 4),
    }

    print(f"\n{'='*60}")
    print("SUMMARY ACROSS ALL REGIONS")
    print(f"{'='*60}")
    print(f"  B cell mean GT: {mean_gt:.4f}")
    print(f"  B cell mean pred WITH gate: {mean_on:.4f}")
    print(f"  B cell mean pred WITHOUT gate: {mean_off:.4f}")
    print(f"  Fold over-prediction without gate: {fold_overpred_without_gate:.2f}x")
    print(f"  Fold over-prediction with gate: {fold_overpred_with_gate:.2f}x")
    print(f"  Fold REDUCTION by gating: {fold_reduction:.2f}x")
    print(f"  B cell Pearson r WITH gate: {r_all_on:.4f}")
    print(f"  B cell Pearson r WITHOUT gate: {r_all_off:.4f}")

    out_path = OUTPUT_DIR / "detection_gate_ablation_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
