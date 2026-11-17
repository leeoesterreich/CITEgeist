"""Re-score Module 3.5 benchmark using updated GT calling and pair list.

Reads existing projected_functional_calls.parquet (from the GPU benchmark run)
and re-evaluates with:
  - GMM-based GT calling (replaces quantile(0.8) which collapsed on sparse markers)
  - PD-1 and PD-L1 only (PCNA/Vimentin removed; both fail Module 1 kurtosis gate)

CPU-only.  Run on login node after the GPU benchmark has produced projected parquets.

Usage:
    python rescore_module3_5_benchmark.py
    python rescore_module3_5_benchmark.py --regions 0 1 2
"""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from CITEgeist.model.module3_5_benchmark import (
    aggregate_module3_5_results,
    build_gt_binary_calls,
    score_pair_predictions,
)
from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

BENCHMARK_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_module3_5_benchmark_array"
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"

THRESHOLDS = {"auroc": 0.75, "auprc": 0.50, "balanced_accuracy": 0.65}

GT_TO_6TYPE = {
    "B cells": "B cells",
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}

FUNCTIONAL_TABLE = {
    "PD-1": {"function": "T cell exhaustion", "active_types": ["T cells"]},
    "PD-L1": {"function": "Immune evasion", "active_types": ["Macrophages", "Epithelial"]},
}

HEADLINE_PAIRS = {"T cells__PD-1", "Macrophages__PD-L1", "Epithelial__PD-L1"}


def load_protein_gt(region_id: int) -> pd.DataFrame:
    """Load cell-level protein GT for the region, joined with 6-type labels."""
    _, adata_protein = load_xenium_singlecell(region_id=region_id)
    mapping = pd.read_csv(PSEUDOVISIUM_DIR / "cell_to_spot_mapping.csv", index_col=0)
    labels = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)
    mapping = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)].copy()
    mapping = mapping.join(labels, how="left")
    mapping = mapping[mapping["cell_type"].isin(GT_TO_6TYPE)].copy()
    mapping["cell_type_6"] = mapping["cell_type"].map(GT_TO_6TYPE)

    protein_df = adata_protein.to_df()
    protein_df["cell_type"] = mapping["cell_type_6"].reindex(protein_df.index)
    protein_df = protein_df.dropna(subset=["cell_type"]).copy()
    return protein_df


def rescore_region(region_id: int) -> dict[str, Any]:
    sample_name = f"Xenium_region_{region_id}"
    parquet_path = BENCHMARK_DIR / sample_name / "projected_functional_calls.parquet"

    if not parquet_path.exists():
        logger.error("Missing projected parquet for %s — run GPU benchmark first", sample_name)
        sys.exit(1)

    projected = pd.read_parquet(parquet_path)
    if "cell_id" in projected.columns:
        projected = projected.set_index("cell_id")

    logger.info("Region %d: loaded %d projected cells", region_id, len(projected))

    protein_df = load_protein_gt(region_id)
    logger.info("Region %d: %d cells with GT labels", region_id, len(protein_df))

    pair_results = []
    gt_summary = []

    for marker, config in FUNCTIONAL_TABLE.items():
        for cell_type in config["active_types"]:
            pair_id = f"{cell_type}__{marker}"
            call_col = f"func_{marker}_{cell_type}_call"
            score_col = f"func_{marker}_{cell_type}_score"

            if call_col not in projected.columns or score_col not in projected.columns:
                logger.warning("  %s: prediction columns missing, skipping", pair_id)
                continue

            gt_calls = build_gt_binary_calls(protein_df, cell_type=cell_type, marker=marker)

            # Report GT calling diagnostics
            n_pos = int(gt_calls.sum())
            n_total = len(gt_calls)
            gt_summary.append(
                {
                    "pair_id": pair_id,
                    "n_gt_cells": n_total,
                    "n_gt_positive": n_pos,
                    "gt_frac_positive": round(n_pos / n_total, 4) if n_total else 0,
                }
            )
            logger.info(
                "  %s GT: %d/%d positive (%.1f%%)", pair_id, n_pos, n_total, 100 * n_pos / n_total if n_total else 0
            )

            pred_scores = projected[score_col]
            pred_calls = projected[call_col]

            supporting_spots = int(
                protein_df.loc[protein_df["cell_type"] == cell_type].index.isin(projected.index).sum()
            )

            result = score_pair_predictions(
                pair_id=pair_id,
                cell_type=cell_type,
                marker=marker,
                gt_calls=gt_calls,
                pred_scores=pred_scores,
                pred_calls=pred_calls,
                n_supporting_spots=supporting_spots,
                headline=pair_id in HEADLINE_PAIRS,
            )
            pair_results.append(result)
            logger.info(
                "  %s: auroc=%.3f auprc=%.3f bal_acc=%.3f evaluable=%s",
                pair_id,
                result.auroc,
                result.auprc,
                result.balanced_accuracy,
                result.evaluable,
            )

    summary = aggregate_module3_5_results(pair_results, thresholds=THRESHOLDS)
    metrics_df = pd.DataFrame(summary["pairs"])

    out_dir = BENCHMARK_DIR / sample_name
    metrics_df.to_csv(out_dir / "module3_5_pair_metrics_gmm.csv", index=False)
    (out_dir / "module3_5_summary_gmm.json").write_text(
        json.dumps({**summary, "gt_summary": gt_summary}, indent=2, sort_keys=True) + "\n"
    )

    logger.info(
        "Region %d: %d pairs, %d evaluable, %d validated, passed=%s",
        region_id,
        summary["n_pairs"],
        summary["n_evaluable_pairs"],
        len(summary["validated_pairs"]),
        summary["benchmark_passed"],
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", type=int, nargs="+", default=list(range(5)), choices=range(5), metavar="REGION_ID")
    args = parser.parse_args()

    all_summaries = {}
    for rid in args.regions:
        all_summaries[rid] = rescore_region(rid)

    print("\n" + "=" * 65)
    print("RESCORE SUMMARY (GMM GT, PD-1/PD-L1 only)")
    print("=" * 65)
    print(f"{'Region':<8} {'Pairs':<8} {'Eval':<8} {'Valid':<8} {'Pass':<8} {'Headline failures'}")
    print("-" * 65)
    for rid, s in sorted(all_summaries.items()):
        print(
            f"{rid:<8} {s['n_pairs']:<8} {s['n_evaluable_pairs']:<8} "
            f"{len(s['validated_pairs']):<8} {str(s['benchmark_passed']):<8} "
            f"{s['headline_failures']}"
        )

    n_pass = sum(1 for s in all_summaries.values() if s["benchmark_passed"])
    print(f"\nOverall: {n_pass}/{len(all_summaries)} regions passed")

    # Write aggregate
    rows = []
    for rid, s in sorted(all_summaries.items()):
        rows.append(
            {
                "region_id": rid,
                "benchmark_passed": s["benchmark_passed"],
                "n_pairs": s["n_pairs"],
                "n_evaluable_pairs": s["n_evaluable_pairs"],
                "n_validated_pairs": len(s["validated_pairs"]),
                "n_headline_failures": len(s["headline_failures"]),
            }
        )
    pd.DataFrame(rows).to_csv(BENCHMARK_DIR / "benchmark_summary_gmm.csv", index=False)
    logger.info("Wrote %s", BENCHMARK_DIR / "benchmark_summary_gmm.csv")


if __name__ == "__main__":
    main()
