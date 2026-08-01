#!/usr/bin/env python
"""Generate SingleR-based GEX ground truth for Xenium pseudo-Visium.

Joins singler_labels_7type.csv + cell_to_spot_mapping.csv + xenium_cell_gex.csv
on cell_id. Groups by region_id + spot_idx, sums GEX per gene per type.
Output: spots × genes CSV per cell type per region.

Usage:
    python generate_singler_gex_gt.py [--dry_run]
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path("/path/to/CITEgeist_analysis")
sys.path.insert(0, str(REPO_ROOT))

SINGLER_DIR = REPO_ROOT / "benchmarks" / "xenium_benchmarking" / "ground_truth_singler" / "singler_labels"
CELL_TO_SPOT = REPO_ROOT / "benchmarks" / "xenium_pseudovisium" / "data_protein_gt" / "cell_to_spot_mapping.csv"
OUT_ROOT = (
    REPO_ROOT
    / "benchmarks"
    / "xenium_benchmarking"
    / "ground_truth_singler"
    / "ground_truth_7type"
    / "ground_truth_gex"
)

CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
    "Macrophages",
]

# SingleR label → safe filename stem
LABEL_TO_FILENAME = {
    "B cells": "B_cells",
    "CD4+ T cells": "CD4pos_T_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
    "Macrophages": "Macrophages",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry_run", action="store_true", help="Check joins only, do not write files")
    args = parser.parse_args()

    # --- Load labels ---
    logger.info("Loading singler labels...")
    labels = pd.read_csv(SINGLER_DIR / "singler_labels_7type.csv")
    # Normalize: column is 'cell_id'
    labels = labels[["cell_id", "singler_label"]].copy()
    logger.info("  Labels shape: %s, unique types: %s", labels.shape, labels["singler_label"].unique().tolist())

    # --- Load cell-to-spot mapping ---
    logger.info("Loading cell-to-spot mapping...")
    mapping = pd.read_csv(CELL_TO_SPOT, index_col=0)
    mapping.index.name = "cell_id"
    mapping = mapping.reset_index()  # cell_id is now a column
    # Drop rows where region_id is NaN or spot_idx < 0 (183 known bad rows)
    n_before = len(mapping)
    mapping = mapping.dropna(subset=["region_id"])
    mapping = mapping[mapping["spot_idx"] >= 0].copy()
    mapping["region_id"] = mapping["region_id"].astype(int)
    mapping["spot_idx"] = mapping["spot_idx"].astype(int)
    logger.info(
        "  Mapping: %d rows before filter, %d after (dropped %d bad rows)",
        n_before,
        len(mapping),
        n_before - len(mapping),
    )

    # --- Inner join labels + mapping ---
    joined = labels.merge(mapping[["cell_id", "spot_idx", "spot_id", "region_id"]], on="cell_id", how="inner")
    logger.info("  After label-mapping join: %d cells", len(joined))

    # --- Load per-cell GEX (large file: 465K × 405 genes) ---
    logger.info("Loading per-cell GEX (may take ~60s)...")
    gex = pd.read_csv(SINGLER_DIR / "xenium_cell_gex.csv", index_col=0)
    gex.index.name = "cell_id"
    gene_names = list(gex.columns)
    logger.info("  GEX shape: %s", gex.shape)

    # --- Inner join with GEX ---
    # Only join cells that have both a label and a spot assignment
    cells_to_join = joined["cell_id"].values
    gex_subset = gex.loc[gex.index.isin(cells_to_join)].copy()
    gex_subset = gex_subset.reset_index()  # cell_id column
    full = joined.merge(gex_subset, on="cell_id", how="inner")
    logger.info("  After GEX join: %d cells", len(full))

    if args.dry_run:
        logger.info("Dry run complete. Exiting without writing files.")
        return

    # --- Generate GT files per region per type ---
    n_files_written = 0
    for region_id in range(5):
        region_data = full[full["region_id"] == region_id].copy()
        logger.info("Region %d: %d cells", region_id, len(region_data))

        out_dir = OUT_ROOT / f"Xenium_region_{region_id}"
        out_dir.mkdir(parents=True, exist_ok=True)

        for label in CELL_TYPES:
            fname_stem = LABEL_TO_FILENAME[label]
            type_data = region_data[region_data["singler_label"] == label].copy()

            if len(type_data) == 0:
                logger.warning("  Region %d / %s: 0 cells — skipping", region_id, label)
                continue

            # Sum GEX by spot_idx
            type_data_gex = type_data[["spot_idx"] + gene_names].copy()
            spot_sums = type_data_gex.groupby("spot_idx")[gene_names].sum()

            # Map spot_idx to canonical spot_{N} names
            spot_sums.index = [f"spot_{int(i)}" for i in spot_sums.index]
            spot_sums.index.name = "spot_id"

            # Save as spots × genes (spots as rows, genes as columns)
            out_path = out_dir / f"{fname_stem}_GT.csv"
            spot_sums.to_csv(out_path)
            n_files_written += 1
            logger.info("  Wrote %s: %s", out_path.name, spot_sums.shape)

    logger.info("Done. Wrote %d files.", n_files_written)

    # --- Validate ---
    logger.info("Validating output...")
    expected = 5 * len(CELL_TYPES)
    n_found = sum(
        1
        for r in range(5)
        for ct in LABEL_TO_FILENAME.values()
        if (OUT_ROOT / f"Xenium_region_{r}" / f"{ct}_GT.csv").exists()
    )
    logger.info("  Expected: %d files, found: %d files", expected, n_found)
    if n_found < expected:
        logger.warning("  MISSING %d files!", expected - n_found)
    else:
        logger.info("  All 35 files present.")

    # Check spot dimensions match pseudo-Visium (should be ~1000-1400 per region)
    for region_id in range(5):
        sample_file = OUT_ROOT / f"Xenium_region_{region_id}" / "Epithelial_GT.csv"
        if sample_file.exists():
            df = pd.read_csv(sample_file, index_col=0)
            logger.info("  Region %d Epithelial: %d spots × %d genes", region_id, df.shape[0], df.shape[1])


if __name__ == "__main__":
    main()
