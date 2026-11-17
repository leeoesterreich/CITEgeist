#!/usr/bin/env python
"""Generate SingleR-based GEX ground truth for Xenium pseudo-Visium.

Joins a labels CSV (cell_id, singler_label) + cell_to_spot_mapping.csv + xenium_cell_gex.csv
on cell_id. Groups by region_id + spot_idx, sums GEX per gene per type.
Output: spots × genes CSV per cell type per region.

Usage:
    python generate_singler_gex_gt.py --labels-csv L --out-root O --props-dir P [--dry_run]
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

REPO_ROOT = Path("/path/to/CITEgeist_analysis")
sys.path.insert(0, str(REPO_ROOT))

SINGLER_DIR = REPO_ROOT / "benchmarks" / "xenium_benchmarking" / "ground_truth_singler" / "singler_labels"
CELL_TO_SPOT = REPO_ROOT / "benchmarks" / "xenium_pseudovisium" / "data_protein_gt" / "cell_to_spot_mapping.csv"

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


def build_type_layer(
    region_cells: pd.DataFrame, label: str, canonical_spots: Sequence[str], gene_names: Sequence[str]
) -> pd.DataFrame:
    """Per-gene GEX sums of one cell type per spot, zero-filled to the canonical spot set.

    ``region_cells`` needs columns ``spot_idx``, ``singler_label`` and the genes. Returns exactly
    ``len(canonical_spots)`` rows (index ``spot_id``, canonical order); spots without the type are zero.
    Raises ValueError if a cell's ``spot_<idx>`` is not in ``canonical_spots``.
    """
    canonical_spots = list(canonical_spots)
    type_data = region_cells[region_cells["singler_label"] == label]
    sums = type_data.groupby("spot_idx")[list(gene_names)].sum()
    sums.index = [f"spot_{int(i)}" for i in sums.index]
    known = set(canonical_spots)
    unknown = [s for s in sums.index if s not in known]
    if unknown:
        raise ValueError(f"{len(unknown)} spot(s) of type {label!r} not in canonical spot set, e.g. {unknown[:5]}")
    layer = sums.reindex(canonical_spots, fill_value=0)
    layer.index.name = "spot_id"
    return layer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels-csv", required=True, type=Path, help="CSV with cell_id, singler_label columns")
    parser.add_argument("--out-root", required=True, type=Path, help="Output root (no default; never the v40 dir)")
    parser.add_argument(
        "--props-dir", required=True, type=Path, help="Dir with Xenium_region_{i}_prop.csv (canonical spot set)"
    )
    parser.add_argument("--gex-csv", type=Path, default=SINGLER_DIR / "xenium_cell_gex.csv")
    parser.add_argument("--mapping-csv", type=Path, default=CELL_TO_SPOT)
    parser.add_argument("--dry_run", action="store_true", help="Check joins only, do not write files")
    args = parser.parse_args()

    # --- Load labels ---
    logger.info("Loading singler labels...")
    labels = pd.read_csv(args.labels_csv)
    # Normalize: column is 'cell_id'
    labels = labels[["cell_id", "singler_label"]].copy()
    logger.info("  Labels shape: %s, unique types: %s", labels.shape, labels["singler_label"].unique().tolist())

    # --- Load cell-to-spot mapping ---
    logger.info("Loading cell-to-spot mapping...")
    mapping = pd.read_csv(args.mapping_csv, index_col=0)
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
    gex = pd.read_csv(args.gex_csv, index_col=0)
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

    # --- Build every layer first (validates all regions) so a failure writes nothing ---
    layers = {}
    for region_id in range(5):
        region_data = full[full["region_id"] == region_id]
        logger.info("Region %d: %d cells", region_id, len(region_data))
        canonical_spots = pd.read_csv(args.props_dir / f"Xenium_region_{region_id}_prop.csv", index_col=0).index
        for label in CELL_TYPES:
            if not (region_data["singler_label"] == label).any():
                logger.warning("  Region %d / %s: 0 cells — writing all-zero layer", region_id, label)
            layers[(region_id, label)] = build_type_layer(region_data, label, canonical_spots, gene_names)

    # --- Write ---
    n_files_written = 0
    for (region_id, label), spot_sums in layers.items():
        out_dir = args.out_root / f"Xenium_region_{region_id}"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{LABEL_TO_FILENAME[label]}_GT.csv"
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
        if (args.out_root / f"Xenium_region_{r}" / f"{ct}_GT.csv").exists()
    )
    logger.info("  Expected: %d files, found: %d files", expected, n_found)
    if n_found < expected:
        logger.warning("  MISSING %d files!", expected - n_found)
    else:
        logger.info("  All 35 files present.")

    # Check spot dimensions match pseudo-Visium (should be ~1000-1400 per region)
    for region_id in range(5):
        sample_file = args.out_root / f"Xenium_region_{region_id}" / "Epithelial_GT.csv"
        if sample_file.exists():
            df = pd.read_csv(sample_file, index_col=0)
            logger.info("  Region %d Epithelial: %d spots × %d genes", region_id, df.shape[0], df.shape[1])


if __name__ == "__main__":
    main()
