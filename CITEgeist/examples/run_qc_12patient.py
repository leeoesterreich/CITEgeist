#!/usr/bin/env python
"""Run self-consistency QC on all 12 HCC22-088 patient samples.

For each sample, loads the SACE per-cell h5ad from morphology_assignment_v3
and the QP proportion CSV, then runs run_qc() in self_consistency mode.

Self-consistency mode checks intra-sample biology:
  - Single-cell QC: per-cell expression totals, genes detected, MT fraction
  - Marker enrichment: canonical RNA marker log2FC per cell type

Note: h5ad.X contains SACE-allocated continuous expression (not raw counts).
UMI thresholds are applied to SACE totals; semantics hold for empty cell
detection even without raw counts.

Usage:
    # Single sample (test):
    python run_qc_12patient.py --sample HCC22-088-P1-S1

    # All samples via SLURM array:
    sbatch CITEgeist/slurm/sbatch_qc_12patient.sh
"""

import argparse
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd

from model import run_qc

REPO = Path(__file__).resolve().parents[2]

SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

MORPH_DIR = REPO / "output" / "morphology_assignment_v3"
QP_DIR = REPO / "output" / "module3_cuopt_qp"
QC_OUTPUT_DIR = REPO / "output" / "qc_12patient"


def run_qc_sample(sample_name: str) -> None:
    """Run self-consistency QC for a single sample."""
    import anndata as ad

    out_dir = QC_OUTPUT_DIR / sample_name
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / "qc_summary.json"
    if summary_path.exists():
        logger.info("QC already complete for %s, skipping", sample_name)
        return

    # --- Load per-cell h5ad ---
    h5ad_path = MORPH_DIR / sample_name / f"{sample_name}_single_cell.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Per-cell h5ad not found: {h5ad_path}")

    logger.info("Loading per-cell h5ad for %s", sample_name)
    adata = ad.read_h5ad(str(h5ad_path))

    # QC framework expects obs["spot_id"]; h5ad uses "spot_barcode"
    if "spot_barcode" in adata.obs.columns and "spot_id" not in adata.obs.columns:
        adata.obs = adata.obs.rename(columns={"spot_barcode": "spot_id"})
        logger.info("Renamed obs['spot_barcode'] -> obs['spot_id']")

    logger.info(
        "%s: %d cells x %d genes, types: %s",
        sample_name,
        *adata.shape,
        dict(adata.obs["cell_type"].value_counts()),
    )

    # --- Load QP proportions ---
    props_path = QP_DIR / sample_name / f"{sample_name}_cell_prop_global_results.csv"
    if not props_path.exists():
        raise FileNotFoundError(f"Proportions CSV not found: {props_path}")

    props_df = pd.read_csv(props_path, index_col=0)
    # Drop solver diagnostic column if present
    props_df = props_df.drop(columns=["recon_error"], errors="ignore")
    logger.info(
        "Loaded proportions: %d spots x %d types", *props_df.shape
    )

    # --- Run self-consistency QC ---
    logger.info("Running self-consistency QC for %s", sample_name)
    run_qc(
        adata_per_cell=adata,
        proportions=props_df,
        mode="self_consistency",
        output_dir=str(out_dir),
    )
    logger.info("QC complete for %s -> %s", sample_name, out_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run self-consistency QC for HCC22-088 samples")
    parser.add_argument("--sample", type=str, default=None, help="Single sample to process")
    parser.add_argument("--array-idx", type=int, default=None, help="SLURM array index (0-11)")
    args = parser.parse_args()

    if args.sample:
        run_qc_sample(args.sample)
    elif args.array_idx is not None:
        sample = SAMPLES[args.array_idx]
        logger.info("Array task %d → %s", args.array_idx, sample)
        run_qc_sample(sample)
    else:
        for sample in SAMPLES:
            run_qc_sample(sample)


if __name__ == "__main__":
    main()
