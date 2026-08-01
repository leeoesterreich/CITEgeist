"""
Full method comparison against ground truth.

Evaluates CITEgeist, Cell2Location, Tangram, RCTD, Seurat, and CARD
against either protein-gated or RNA-based ground truth.

Supports:
- --gt-dir: Specify ground truth directory (default: protein GT)
- --combine-t-cells: Merge CD4+ and CD8+ T cells in predictions (for RNA GT)

Note on CARD:
- Uses VIM-excluded results (output_protein_gt_novim) because VIM is a pan-mesenchymal
  marker expressed in RCC cells, not fibroblast-specific in kidney tissue.
- CARD ref-free is not compatible with Xenium's 405-gene panel (<20 markers/cell type).
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

# Add evaluation src to path
sys.path.insert(0, str(Path(__file__).parent))
from evaluate_benchmark import evaluate_all_regions, ACHIEVABLE_7_CELL_TYPES

# Add benchmarking root to path for benchmark_constants
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from benchmark_constants import GT_TO_ACHIEVABLE_6_MAPPING, GT_TO_ACHIEVABLE_7_MAPPING

# RNA GT uses 6 cell types (T cells combined)
RNA_6_CELL_TYPES = [
    "B cells",
    "T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path("/path/to/CITEgeist_analysis/benchmarks/xenium_benchmarking")
GT_DIR = Path("/path/to/CITEgeist_analysis/benchmarks/xenium_pseudovisium/data_protein_gt")
GT_DIR_SINGLER_6TYPE = BASE_DIR / "ground_truth_singler" / "ground_truth_6type"
GT_DIR_SINGLER_7TYPE = BASE_DIR / "ground_truth_singler" / "ground_truth_7type"
RESULTS_DIR = BASE_DIR / "evaluation" / "results" / "method_comparison"

# =============================================================================
# METHOD DOCUMENTATION
# =============================================================================

METHOD_NOTES = {
    "CITEgeist": "Full pipeline: QP + cellularity prior + multimodal EM + PC-MIL SC assignment",
}

CITEGEIST_COMMON_PARAMS = {
    "data_source": "Xenium pseudo-Visium (SingleR ground truth)",
    "min_counts": 25,
    "gex_preprocessing": "filter_gex + preprocess_gex(target_sum=10000)",
    "cell_profile_dict": "ACHIEVABLE_6/7_CELL_PROFILE_DICT (alphaSMA-only fibroblasts)",
    "n_regions": 5,
    "pipeline": "StarDist nuclei → QP + cellularity prior → multimodal EM → GMM detection → PC-MIL",
}

METHODS = {
    # CITEgeist (unified pipeline: QP + EM + PC-MIL)
    "CITEgeist": BASE_DIR / "CITEgeist" / "output",
    # Other methods
    "Cell2Location": BASE_DIR / "Cell2Location" / "output_protein_gt",
    "Tangram": BASE_DIR / "Tangram" / "output_protein_gt",
    "RCTD": BASE_DIR / "RCTD" / "output_protein_gt",
    "Seurat": BASE_DIR / "Seurat" / "output_protein_gt",
    # CARD: VIM excluded - VIM is a pan-mesenchymal marker expressed in RCC cells,
    # not fibroblast-specific in kidney tissue. Including it caused severe bias (93% fibroblasts).
    "CARD": BASE_DIR / "CARD" / "output_protein_gt_novim",
}

# SingleR GT method entries (6-type: primary comparison, 7-type: secondary)
METHODS_SINGLER_6TYPE = {
    "CITEgeist": BASE_DIR / "CITEgeist" / "output_singler_6type",
    "Cell2Location": BASE_DIR / "Cell2Location" / "output_singler_6type",
    "RCTD": BASE_DIR / "RCTD" / "output_singler_6type",
    "Tangram": BASE_DIR / "Tangram" / "output_singler_6type",
    "Seurat": BASE_DIR / "Seurat" / "output_singler_6type",
    "CARD": BASE_DIR / "CARD" / "output_singler_6type",
}

METHODS_SINGLER_7TYPE = {
    "CITEgeist": BASE_DIR / "CITEgeist" / "output_singler_7type",
    "Cell2Location": BASE_DIR / "Cell2Location" / "output_singler_7type",
    "RCTD": BASE_DIR / "RCTD" / "output_singler_7type",
    "Tangram": BASE_DIR / "Tangram" / "output_singler_7type",
    "Seurat": BASE_DIR / "Seurat" / "output_singler_7type",
    "CARD": BASE_DIR / "CARD" / "output_singler_7type",
}


def evaluate_single_cell_accuracy(
    assignments_dir: Path,
    gt_dir: Path,
    cell_type_mapping: Dict[str, str],
    region_names: List[str],
) -> Dict[str, Any]:
    """Evaluate single-cell assignment accuracy against Xenium GT labels.

    Args:
        assignments_dir: Directory with {region}_single_cell_assignments.csv
        gt_dir: Directory with Xenium GT cell labels
        cell_type_mapping: GT label -> evaluation label mapping
        region_names: List of region names to evaluate

    Returns:
        Dict with overall_accuracy, per_type_accuracy, per_region results
    """
    all_correct = 0
    all_total = 0
    per_type_correct = {}
    per_type_total = {}

    for region_name in region_names:
        assignments_path = assignments_dir / region_name / f"{region_name}_single_cell_assignments.csv"
        if not assignments_path.exists():
            logger.warning(f"No assignments found for {region_name}")
            continue

        assignments = pd.read_csv(assignments_path)
        gt_path = gt_dir / f"{region_name}_gt_labels.csv"
        if not gt_path.exists():
            logger.warning(f"No GT labels found for {region_name}")
            continue

        gt_labels = pd.read_csv(gt_path)

        # Match by nucleus_id
        merged = assignments.merge(gt_labels, on="nucleus_id", how="inner", suffixes=("_pred", "_gt"))

        for _, row in merged.iterrows():
            gt_type = cell_type_mapping.get(row["cell_type_gt"])
            if gt_type is None:
                continue  # Skip unknown GT types

            pred_type = row["cell_type"]
            correct = pred_type == gt_type

            all_total += 1
            all_correct += int(correct)
            per_type_total[gt_type] = per_type_total.get(gt_type, 0) + 1
            per_type_correct[gt_type] = per_type_correct.get(gt_type, 0) + int(correct)

    overall_accuracy = all_correct / max(all_total, 1)
    per_type_accuracy = {
        ct: per_type_correct.get(ct, 0) / max(per_type_total.get(ct, 0), 1) for ct in sorted(per_type_total.keys())
    }

    return {
        "overall_accuracy": overall_accuracy,
        "per_type_accuracy": per_type_accuracy,
        "n_cells_evaluated": all_total,
    }


def main():
    parser = argparse.ArgumentParser(description="Compare all deconvolution methods")
    parser.add_argument(
        "--gt-dir",
        type=str,
        default=str(GT_DIR),
        help="Ground truth directory (default: protein GT)",
    )
    parser.add_argument(
        "--combine-t-cells",
        action="store_true",
        help="Combine CD4+ and CD8+ T cells in predictions (for RNA GT)",
    )
    parser.add_argument(
        "--output-suffix",
        type=str,
        default="",
        help="Suffix for output files (e.g., '_rna_gt')",
    )
    parser.add_argument(
        "--gt-type",
        choices=["protein", "singler_6type", "singler_7type"],
        default=None,
        help="Convenience flag to select GT type (overrides --gt-dir and --combine-t-cells)",
    )
    args = parser.parse_args()

    # --gt-type convenience overrides
    if args.gt_type == "singler_6type":
        gt_dir = GT_DIR_SINGLER_6TYPE
        combine_t = True
        methods = METHODS_SINGLER_6TYPE
        suffix = "_singler_6type"
    elif args.gt_type == "singler_7type":
        gt_dir = GT_DIR_SINGLER_7TYPE
        combine_t = False
        methods = METHODS_SINGLER_7TYPE
        suffix = "_singler_7type"
    else:
        gt_dir = Path(args.gt_dir)
        combine_t = args.combine_t_cells
        methods = METHODS
        suffix = args.output_suffix

    # Determine cell types based on GT
    if combine_t:
        cell_types = RNA_6_CELL_TYPES
        print("Using RNA 6-cell-type mode (T cells combined)")
    else:
        cell_types = ACHIEVABLE_7_CELL_TYPES
        print("Using protein 7-cell-type mode")

    results_dir = RESULTS_DIR
    if suffix:
        results_dir = RESULTS_DIR.parent / f"method_comparison{suffix}"
    results_dir.mkdir(parents=True, exist_ok=True)

    all_results = {}

    for method_name, pred_dir in methods.items():
        print(f"\n{'='*60}")
        print(f"Evaluating: {method_name}")
        print(f"{'='*60}")

        if not pred_dir.exists():
            print(f"  SKIPPED: {pred_dir} does not exist")
            continue

        try:
            results = evaluate_all_regions(
                gt_dir=str(gt_dir),
                pred_dir=str(pred_dir),
                n_regions=5,
                output_path=str(results_dir / f"{method_name}_results.json"),
                prefix="Xenium",
                use_achievable_7=not combine_t,  # Use raw GT columns if combining T cells
                combine_t_cells=combine_t,  # Pass flag to evaluation
            )
            all_results[method_name] = results
            print(f"  Overall Pearson r: {results['summary']['overall_mean_pearson_r']:.4f}")
            print(f"  Overall RMSE:     {results['summary']['overall_mean_rmse']:.4f}")
            print(f"  Overall MAE:      {results['summary']['overall_mean_mae']:.4f}")
            if "overall_mean_jsd" in results["summary"]:
                print(f"  Overall JSD:      {results['summary']['overall_mean_jsd']:.4f}")
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback

            traceback.print_exc()

    if not all_results:
        print("\nNo methods evaluated successfully!")
        return

    # =========================================================================
    # Single-cell accuracy evaluation (CITEgeist only)
    # =========================================================================
    sc_results = None
    if "CITEgeist" in all_results:
        citegeist_dir = methods["CITEgeist"]
        region_names = [f"Xenium_region_{i}" for i in range(5)]

        # Select cell type mapping based on GT type
        if combine_t:
            sc_mapping = GT_TO_ACHIEVABLE_6_MAPPING
        else:
            sc_mapping = GT_TO_ACHIEVABLE_7_MAPPING

        sc_results = evaluate_single_cell_accuracy(
            assignments_dir=citegeist_dir,
            gt_dir=gt_dir,
            cell_type_mapping=sc_mapping,
            region_names=region_names,
        )

        if sc_results["n_cells_evaluated"] > 0:
            print(f"\n{'='*60}")
            print("CITEgeist Single-Cell Accuracy")
            print(f"{'='*60}")
            print(f"  Overall accuracy: {sc_results['overall_accuracy']:.4f}")
            print(f"  Cells evaluated:  {sc_results['n_cells_evaluated']}")
            print(f"  Per-type accuracy:")
            for ct, acc in sc_results["per_type_accuracy"].items():
                print(f"    {ct:<20}: {acc:.4f}")
            all_results["CITEgeist"]["single_cell"] = sc_results
        else:
            print("\n  No single-cell assignments found for CITEgeist")

    # =========================================================================
    # Summary comparison table
    # =========================================================================
    print("\n\n" + "=" * 80)
    print("FULL METHOD COMPARISON - Protein-Gated Ground Truth (Achievable-7)")
    print("=" * 80)

    # Print method notes
    print("\n--- Method Notes ---")
    for method, note in METHOD_NOTES.items():
        print(f"  {method}: {note}")

    print("\n--- Common CITEgeist Parameters ---")
    for key, val in CITEGEIST_COMMON_PARAMS.items():
        print(f"  {key}: {val}")

    # Overall metrics table
    print("\n--- Overall Metrics ---")
    print(f"{'Method':<24} {'Pearson r':>10} {'RMSE':>10} {'MAE':>10} {'JSD':>10}")
    print("-" * 64)

    # Sort by Pearson r descending
    sorted_methods = sorted(
        all_results.keys(),
        key=lambda m: all_results[m]["summary"]["overall_mean_pearson_r"],
        reverse=True,
    )

    for method in sorted_methods:
        s = all_results[method]["summary"]
        jsd_str = f"{s['overall_mean_jsd']:.4f}" if "overall_mean_jsd" in s else "N/A"
        print(
            f"{method:<24} "
            f"{s['overall_mean_pearson_r']:>10.4f} "
            f"{s['overall_mean_rmse']:>10.4f} "
            f"{s['overall_mean_mae']:>10.4f} "
            f"{jsd_str:>10}"
        )

    # Per-cell-type Pearson r comparison
    print("\n--- Per-Cell-Type Pearson r (mean ± std across regions) ---")
    header = f"{'Cell Type':<20}"
    for method in sorted_methods:
        header += f" {method:>16}"
    print(header)
    print("-" * (20 + 17 * len(sorted_methods)))

    for ct in cell_types:
        row = f"{ct:<20}"
        for method in sorted_methods:
            s = all_results[method]["summary"]
            mean_key = f"{ct}_mean_r"
            std_key = f"{ct}_std_r"
            if mean_key in s:
                mean_r = s[mean_key]
                std_r = s.get(std_key, 0)
                row += f" {mean_r:>7.3f}±{std_r:.3f}"
            else:
                row += f" {'N/A':>16}"
        print(row)

    # Mean of per-cell-type r
    print(f"{'MEAN':.<20}", end="")
    for method in sorted_methods:
        s = all_results[method]["summary"]
        ct_rs = []
        for ct in cell_types:
            mean_key = f"{ct}_mean_r"
            if mean_key in s:
                ct_rs.append(s[mean_key])
        if ct_rs:
            print(f" {np.mean(ct_rs):>7.3f}      ", end="")
        else:
            print(f" {'N/A':>16}", end="")
    print()

    # Per-region breakdown
    print(f"\n--- Per-Region Overall Pearson r ---")
    header = f"{'Region':<10}"
    for method in sorted_methods:
        header += f" {method:>14}"
    print(header)
    print("-" * (10 + 15 * len(sorted_methods)))

    for region_id in range(5):
        row = f"Region {region_id:<3}"
        for method in sorted_methods:
            regions = all_results[method]["regions"]
            region_data = [r for r in regions if r["region_id"] == region_id]
            if region_data:
                row += f" {region_data[0]['overall_pearson_r']:>14.4f}"
            else:
                row += f" {'N/A':>14}"
        print(row)

    # Wins table
    print(f"\n--- Category Wins (best per-cell-type r) ---")
    wins = {m: 0 for m in sorted_methods}
    for ct in cell_types:
        best_method = None
        best_r = -1
        for method in sorted_methods:
            s = all_results[method]["summary"]
            mean_key = f"{ct}_mean_r"
            if mean_key in s and s[mean_key] > best_r:
                best_r = s[mean_key]
                best_method = method
        if best_method:
            wins[best_method] += 1
            print(f"  {ct:<20}: {best_method} (r={best_r:.4f})")

    print(f"\nTotal wins: ", end="")
    for method in sorted_methods:
        print(f"{method}={wins[method]}  ", end="")
    print()

    # Save combined results
    combined_output = RESULTS_DIR / "full_comparison.json"

    def convert_numpy(obj):
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy(v) for v in obj]
        return obj

    output_data = {
        "method_notes": METHOD_NOTES,
        "citegeist_common_params": CITEGEIST_COMMON_PARAMS,
        "results": convert_numpy(all_results),
    }

    with open(combined_output, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"\nFull results saved to: {combined_output}")

    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
