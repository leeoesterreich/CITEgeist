#!/usr/bin/env python
"""Validate Module 3.5 per-cell functional calls against Xenium GT protein expression.

For each (cell_type, marker) pair in the projected calls, compares Module 3.5
binary functional calls against ground-truth derived from Xenium protein expression
using 2-component GMM thresholding.

Reads:  output_module3_5_benchmark_array/Xenium_region_{i}/projected_functional_calls.parquet
Reads:  ground_truth_singler/singler_labels/xenium_cell_proteins.csv
Reads:  xenium_pseudovisium/data_protein_gt/cell_type_assignments.csv
Reads:  xenium_pseudovisium/data_protein_gt/cell_to_spot_mapping.csv
Writes: output_module3_5_benchmark_array/Xenium_region_{i}/singlecell_validation.json

Usage:
    python validate_module3_5_singlecell.py --region-id 0
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.mixture import GaussianMixture

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
GT_PROTEINS_CSV = (
    REPO_ROOT
    / "Benchmarking"
    / "xenium_benchmarking"
    / "ground_truth_singler"
    / "singler_labels"
    / "xenium_cell_proteins.csv"
)
CELL_TYPES_CSV = (
    REPO_ROOT
    / "Benchmarking"
    / "xenium_benchmarking"
    / "xenium_pseudovisium"
    / "data_protein_gt"
    / "cell_type_assignments.csv"
)
CELL_MAPPING_CSV = (
    REPO_ROOT
    / "Benchmarking"
    / "xenium_benchmarking"
    / "xenium_pseudovisium"
    / "data_protein_gt"
    / "cell_to_spot_mapping.csv"
)
DEFAULT_BENCHMARK_DIR = (
    REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_module3_5_benchmark_array"
)

# Module 3.5 marker name -> Xenium protein CSV column name (where they differ)
FUNCTIONAL_TO_XENIUM_PROTEIN: dict[str, str] = {
    "CD274": "PD-L1",
    "PDCD1": "PD-1",
    "MKI67": "Ki-67",
    "VIM": "Vimentin",
    # Add further mappings as needed; unlisted markers use their name directly
}
MIN_CELLS = 20


def load_gt(region_id: int) -> pd.DataFrame:
    """Return per-cell DataFrame for region_id with columns: cell_type, <proteins...>."""
    proteins = pd.read_csv(GT_PROTEINS_CSV, index_col=0)
    cell_types = pd.read_csv(CELL_TYPES_CSV, index_col=0)
    mapping = pd.read_csv(CELL_MAPPING_CSV, index_col=0)
    region_cells = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)].index
    proteins = proteins.loc[proteins.index.intersection(region_cells)]
    cell_types = cell_types.loc[cell_types.index.intersection(region_cells)]
    gt = proteins.join(cell_types[["cell_type"]], how="inner")
    logger.info("GT: %d cells loaded for region %d", len(gt), region_id)
    return gt


def build_binary_gt(gt: pd.DataFrame, cell_type: str, xenium_col: str) -> tuple[pd.Series | None, str | None]:
    """GMM-threshold Xenium protein expression to binary GT calls.

    Returns (binary_series indexed by cell_id, None) on success,
    or (None, exclusion_reason) if not evaluable.
    """
    if xenium_col not in gt.columns:
        return None, f"column '{xenium_col}' not in GT protein data"
    cells = gt[gt["cell_type"] == cell_type]
    if len(cells) < MIN_CELLS:
        return None, f"only {len(cells)} cells of type '{cell_type}' (< {MIN_CELLS})"
    expr = np.log1p(cells[xenium_col].values.astype(float)).reshape(-1, 1)
    try:
        gmm = GaussianMixture(n_components=2, random_state=42)
        gmm.fit(expr)
        means = gmm.means_.ravel()
        stds = np.sqrt(gmm.covariances_.ravel())
        separation = abs(means[1] - means[0]) / max(stds.mean(), 1e-6)
        threshold = float(np.mean(sorted(means))) if separation >= 0.5 else float(np.percentile(expr, 90))
    except Exception as exc:
        logger.warning("GMM failed for (%s, %s): %s; using 90th pct", cell_type, xenium_col, exc)
        threshold = float(np.percentile(expr, 90))
    binary = (expr.ravel() >= threshold).astype(int)
    return pd.Series(binary, index=cells.index, name=xenium_col), None


def score_pair(
    pair_id: str,
    cell_type: str,
    marker: str,
    gt_binary: pd.Series,
    pred_df: pd.DataFrame,
    score_col: str,
    call_col: str,
) -> dict:
    """Compute sc_auroc, sc_auprc, sc_balanced_accuracy for one pair."""
    common = gt_binary.index.intersection(pred_df.index)
    gt_vals = gt_binary.loc[common].values
    base = {
        "pair_id": pair_id,
        "cell_type": cell_type,
        "marker": marker,
        "sc_n_positive": int(gt_vals.sum()) if len(gt_vals) else None,
        "sc_n_negative": int(len(gt_vals) - gt_vals.sum()) if len(gt_vals) else None,
    }

    if len(common) < MIN_CELLS:
        return {
            **base,
            "sc_evaluable": False,
            "sc_exclusion_reason": f"only {len(common)} common cells",
            "sc_auroc": None,
            "sc_auprc": None,
            "sc_balanced_accuracy": None,
        }
    if gt_vals.sum() < 5 or (len(gt_vals) - gt_vals.sum()) < 5:
        return {
            **base,
            "sc_evaluable": False,
            "sc_exclusion_reason": "too few positives or negatives",
            "sc_auroc": None,
            "sc_auprc": None,
            "sc_balanced_accuracy": None,
        }
    try:
        scores_raw = pred_df.loc[common, score_col].values.astype(float)
        calls_raw = pred_df.loc[common, call_col].values.astype(float)
        # Module 3.5 only populates scores for cells assigned to this type;
        # filter to cells with non-NaN scores to avoid spurious all-zero evaluations.
        valid_mask = np.isfinite(scores_raw) & np.isfinite(calls_raw)
        if valid_mask.sum() < MIN_CELLS:
            return {
                **base,
                "sc_evaluable": False,
                "sc_exclusion_reason": f"only {valid_mask.sum()} cells with non-NaN scores",
                "sc_auroc": None,
                "sc_auprc": None,
                "sc_balanced_accuracy": None,
            }
        gt_valid = gt_vals[valid_mask]
        scores = scores_raw[valid_mask]
        calls = calls_raw[valid_mask].astype(int)
        if gt_valid.sum() < 5 or (len(gt_valid) - gt_valid.sum()) < 5:
            return {
                **base,
                "sc_evaluable": False,
                "sc_exclusion_reason": "too few positives or negatives after NaN filter",
                "sc_auroc": None,
                "sc_auprc": None,
                "sc_balanced_accuracy": None,
            }
        return {
            **base,
            "sc_evaluable": True,
            "sc_exclusion_reason": None,
            "sc_auroc": float(roc_auc_score(gt_valid, scores)),
            "sc_auprc": float(average_precision_score(gt_valid, scores)),
            "sc_balanced_accuracy": float(balanced_accuracy_score(gt_valid, calls)),
        }
    except Exception as exc:
        return {
            **base,
            "sc_evaluable": False,
            "sc_exclusion_reason": str(exc),
            "sc_auroc": None,
            "sc_auprc": None,
            "sc_balanced_accuracy": None,
        }


def parse_pairs_from_columns(columns: list[str], known_cell_types: set[str]) -> list[tuple]:
    """Extract (cell_type, marker, call_col, score_col) from projected DataFrame columns.

    Column format: func_{marker}_{cell_type}_call / func_{marker}_{cell_type}_score
    """
    pairs = []
    for col in columns:
        if not (col.startswith("func_") and col.endswith("_call")):
            continue
        inner = col[len("func_") : -len("_call")]
        score_col = f"func_{inner}_score"
        if score_col not in columns:
            continue
        for ct in known_cell_types:
            if inner.endswith(f"_{ct}"):
                marker = inner[: -(len(ct) + 1)]
                pairs.append((ct, marker, col, score_col))
                break
    return pairs


def run_region(region_id: int, benchmark_dir: Path) -> None:
    projected_path = benchmark_dir / f"Xenium_region_{region_id}" / "projected_functional_calls.parquet"
    if not projected_path.exists():
        logger.error("Missing: %s — run the benchmark runner first", projected_path)
        sys.exit(1)

    pred_df = pd.read_parquet(projected_path).set_index("cell_id")
    gt = load_gt(region_id)
    known_cell_types = set(gt["cell_type"].unique())

    pairs = parse_pairs_from_columns(list(pred_df.columns), known_cell_types)
    logger.info("Found %d (cell_type, marker) pairs in projected calls", len(pairs))

    results = []
    for cell_type, marker, call_col, score_col in pairs:
        xenium_col = FUNCTIONAL_TO_XENIUM_PROTEIN.get(marker, marker)
        gt_binary, exclusion_reason = build_binary_gt(gt, cell_type, xenium_col)
        if gt_binary is None:
            logger.info("Skip (%s, %s): %s", cell_type, marker, exclusion_reason)
            results.append(
                {
                    "pair_id": f"{cell_type}__{marker}",
                    "cell_type": cell_type,
                    "marker": marker,
                    "sc_evaluable": False,
                    "sc_exclusion_reason": exclusion_reason,
                    "sc_auroc": None,
                    "sc_auprc": None,
                    "sc_balanced_accuracy": None,
                    "sc_n_positive": None,
                    "sc_n_negative": None,
                }
            )
            continue
        result = score_pair(f"{cell_type}__{marker}", cell_type, marker, gt_binary, pred_df, score_col, call_col)
        results.append(result)
        logger.info(
            "(%s, %s): auroc=%.3f auprc=%.3f bal_acc=%.3f",
            cell_type,
            marker,
            result.get("sc_auroc") or 0,
            result.get("sc_auprc") or 0,
            result.get("sc_balanced_accuracy") or 0,
        )

    n_eval = sum(1 for r in results if r["sc_evaluable"])
    out = {"region_id": region_id, "n_pairs": len(results), "n_evaluable": n_eval, "pairs": results}
    out_path = benchmark_dir / f"Xenium_region_{region_id}" / "singlecell_validation.json"
    out_path.write_text(json.dumps(out, indent=2))
    logger.info("Wrote %s (%d pairs, %d evaluable)", out_path, len(results), n_eval)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-id", type=int, required=True, choices=range(5))
    parser.add_argument("--benchmark-dir", type=Path, default=DEFAULT_BENCHMARK_DIR)
    args = parser.parse_args()
    run_region(args.region_id, args.benchmark_dir)


if __name__ == "__main__":
    main()
