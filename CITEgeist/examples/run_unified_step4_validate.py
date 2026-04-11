#!/usr/bin/env python
"""Step 5: Marker gene validation of PC-MIL assignments.

Usage:
    python run_unified_step4_validate.py --sample HCC22-088-P1-S1
"""
import argparse
import json
import logging
import os
import sys

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from model.marker_validation import compute_marker_scores, summarize_validation
from model.unified_config import OUTPUT_BASE, RNA_MARKERS


def run_step4(sample_name):
    step3_marker = OUTPUT_BASE / sample_name / ".step3_complete"
    if not step3_marker.exists():
        logger.error(f"Step 3 not complete for {sample_name}")
        return

    val_dir = OUTPUT_BASE / sample_name / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)

    assignments = pd.read_csv(OUTPUT_BASE / sample_name / "pcmil" / "assignments.csv")

    m3_dir = OUTPUT_BASE / sample_name / "module3"
    parquet_files = list(m3_dir.glob("*gene_expression_pass1.parquet"))
    if not parquet_files:
        logger.error(f"No GEX parquet found in {m3_dir}")
        return
    gex_df = pd.read_parquet(parquet_files[0])

    scores = compute_marker_scores(assignments, gex_df, RNA_MARKERS)
    scores.to_csv(val_dir / "marker_gene_scores.csv", index=False)

    summary = summarize_validation(scores)
    with open(val_dir / "validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"=== Validation Summary for {sample_name} ===")
    logger.info(f"Total nuclei scored: {summary['overall']['n_nuclei']}")
    logger.info(f"Overall fraction correct: {summary['overall']['fraction_correct']:.3f}")
    logger.info(f"Overall median marker score: {summary['overall']['median_marker_score']:.3f}")
    for ct, metrics in summary["per_type"].items():
        logger.info(
            f"  {ct}: {metrics['fraction_correct']:.3f} correct "
            f"({metrics['n_nuclei']} nuclei, median score={metrics['median_marker_score']:.2f})"
        )

    step4_marker = OUTPUT_BASE / sample_name / ".step4_complete"
    step4_marker.touch()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 5: Validation")
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    run_step4(args.sample)
