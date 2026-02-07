"""
Sweep lambda_sparse values for cell-resolution Pass 1 on a single region.

Tests how sparsity penalty strength affects assignment sharpness and accuracy.
"""
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))

from model.gurobi_impl import (
    optimize_cell_proportions_per_marker,
    map_antibodies_to_profiles_v2,
)
from load_xenium_singlecell import load_xenium_singlecell
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

logger = logging.getLogger(__name__)
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_cell_resolution"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_type_assignments.csv"

LAMBDA_VALUES = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]


def run_sweep(region_id: int = 2, max_cells: int = 5000):
    """Sweep lambda_sparse on one region, Pass 1 only."""
    output_dir = OUTPUT_DIR / "lambda_sweep"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading Xenium region {region_id} (max {max_cells} cells)")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id, max_cells=max_cells, seed=42
    )
    logger.info(f"Loaded {adata_gex.shape[0]} cells")

    gt_df = pd.read_csv(GT_PATH, index_col=0)
    common_cells = sorted(set(adata_protein.obs_names) & set(gt_df.index))
    gt_labels = gt_df.loc[common_cells, "cell_type"].values
    logger.info(f"Ground truth for {len(common_cells)} cells")

    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)
    coords = adata_protein.obsm["spatial"]

    marker_data, mapped_marker_names, assignment_matrix, type_names = map_antibodies_to_profiles_v2(
        adata_protein, ACHIEVABLE_7_CELL_PROFILE_DICT
    )

    gt_types_in_profiles = set(type_names)
    evaluable_mask = np.array([gt in gt_types_in_profiles for gt in gt_labels])

    from sklearn.metrics import accuracy_score, f1_score

    all_results = []

    for lam in LAMBDA_VALUES:
        logger.info("=" * 60)
        logger.info(f"lambda_sparse = {lam}")

        Y, beta_values, beta_dict, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=marker_data,
            marker_names=mapped_marker_names,
            assignment_matrix=assignment_matrix,
            cell_type_names=type_names,
            lambda_sparse=lam,
            lambda_laplacian=0.01,
            coords=coords,
            laplacian_k=50,
            max_iterations=5,
        )

        max_Y = Y.max(axis=1)
        dominant_type = np.argmax(Y, axis=1)
        predicted_all = [type_names[dt] for dt in dominant_type]

        cell_to_pred = dict(zip(adata_protein.obs_names, predicted_all))
        predicted_gt = np.array([cell_to_pred[c] for c in common_cells])

        gt_eval = gt_labels[evaluable_mask]
        pred_eval = predicted_gt[evaluable_mask]

        acc = float(accuracy_score(gt_eval, pred_eval))
        f1_mac = float(f1_score(gt_eval, pred_eval, average="macro", zero_division=0))
        f1_wt = float(f1_score(gt_eval, pred_eval, average="weighted", zero_division=0))

        # Assignment sharpness stats
        entropy = -np.sum(Y * np.log(Y + 1e-10), axis=1)
        max_entropy = np.log(len(type_names))

        result = {
            "lambda_sparse": lam,
            "accuracy": acc,
            "f1_macro": f1_mac,
            "f1_weighted": f1_wt,
            "max_Y_mean": float(max_Y.mean()),
            "max_Y_median": float(np.median(max_Y)),
            "max_Y_p90": float(np.percentile(max_Y, 90)),
            "doublet_fraction": float(np.mean(max_Y < 0.6)),
            "one_hot_fraction": float(np.mean(max_Y > 0.9)),
            "entropy_mean": float(entropy.mean()),
            "entropy_ratio": float(entropy.mean() / max_entropy),
        }
        all_results.append(result)

        logger.info(
            f"  acc={acc:.3f}  f1_mac={f1_mac:.3f}  max_Y={max_Y.mean():.3f}  "
            f"one_hot={result['one_hot_fraction']:.3f}  entropy_ratio={result['entropy_ratio']:.3f}"
        )

    with open(output_dir / "lambda_sweep_results.json", "w") as f:
        json.dump(all_results, f, indent=2)

    # Print summary table
    logger.info("\n" + "=" * 90)
    logger.info("LAMBDA SWEEP SUMMARY")
    logger.info("=" * 90)
    logger.info(f"{'lambda':>8} {'acc':>6} {'f1_mac':>7} {'maxY':>6} {'medY':>6} {'p90Y':>6} {'1hot%':>6} {'ent_r':>6}")
    logger.info("-" * 90)
    for r in all_results:
        logger.info(
            f"{r['lambda_sparse']:8.1f} {r['accuracy']:6.3f} {r['f1_macro']:7.3f} "
            f"{r['max_Y_mean']:6.3f} {r['max_Y_median']:6.3f} {r['max_Y_p90']:6.3f} "
            f"{r['one_hot_fraction']:6.3f} {r['entropy_ratio']:6.3f}"
        )

    logger.info(f"\nResults saved to {output_dir / 'lambda_sweep_results.json'}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_sweep()
