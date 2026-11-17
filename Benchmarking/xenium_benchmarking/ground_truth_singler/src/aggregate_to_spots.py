#!/usr/bin/env python
"""
Aggregate per-cell SingleR labels to pseudo-Visium spot proportions.

Produces ground truth CSVs matching evaluator convention:
  ground_truth/Xenium_region_{0-4}_prop.csv

Two modes:
  --n-types 6: Uses singler_all_cells.csv (T cells merged)
  --n-types 7: Uses singler_labels_7type.csv (CD4/CD8 split, ambiguous excluded)
"""

import argparse
import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

CELL_TYPES_6 = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]
CELL_TYPES_7 = ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels-csv", required=True, help="Path to singler labels CSV")
    parser.add_argument("--mapping-csv", required=True, help="Path to cell_to_spot_mapping.csv")
    parser.add_argument("--output-dir", required=True, help="Output directory for GT CSVs")
    parser.add_argument("--n-types", type=int, choices=[6, 7], required=True)
    parser.add_argument("--prefix", default="Xenium")
    args = parser.parse_args()

    output_dir = Path(args.output_dir) / "ground_truth"
    output_dir.mkdir(parents=True, exist_ok=True)

    cell_types = CELL_TYPES_6 if args.n_types == 6 else CELL_TYPES_7

    # Load labels and mapping
    labels_df = pd.read_csv(args.labels_csv)
    mapping_df = pd.read_csv(args.mapping_csv, index_col=0)
    logger.info(f"Labels: {len(labels_df)}, Mapping: {len(mapping_df)}")

    # Join labels to mapping (cell_id as index)
    labels_df = labels_df.set_index("cell_id")
    merged = mapping_df.join(labels_df[["singler_label"]], how="inner")

    # Filter to mapped cells (spot_idx >= 0)
    mapped = merged[merged["spot_idx"] >= 0].copy()
    logger.info(f"Mapped cells with labels: {len(mapped)}")

    # Filter to known cell types (exclude ambiguous "T cells" in 7-type mode)
    mapped = mapped[mapped["singler_label"].isin(cell_types)]
    logger.info(f"After filtering to {args.n_types} types: {len(mapped)}")

    # Aggregate per region
    for region_id in range(5):
        region_cells = mapped[mapped["region_id"] == region_id]

        # Count per spot per type
        spot_type_counts = region_cells.groupby(["spot_id", "singler_label"]).size().unstack(fill_value=0)

        # Ensure all cell types present
        for ct in cell_types:
            if ct not in spot_type_counts.columns:
                spot_type_counts[ct] = 0
        spot_type_counts = spot_type_counts[cell_types]

        # Compute proportions
        row_sums = spot_type_counts.sum(axis=1)
        props = spot_type_counts.div(row_sums, axis=0).fillna(0)

        # Add metadata
        spot_meta = region_cells.groupby("spot_id").agg(
            n_cells=("singler_label", "count"),
            spot_x=("spot_x", "first"),
            spot_y=("spot_y", "first"),
        )
        props = props.join(spot_meta)

        # Save
        out_path = output_dir / f"{args.prefix}_region_{region_id}_prop.csv"
        props.to_csv(out_path)
        logger.info(f"Region {region_id}: {len(props)} spots -> {out_path}")

    logger.info("Done.")


if __name__ == "__main__":
    main()
