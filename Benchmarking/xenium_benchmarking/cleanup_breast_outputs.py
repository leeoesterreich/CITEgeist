"""
Normalize breast benchmark output column names across all methods.

Fixes formatting artifacts: dotted names (Seurat), extra whitespace, casing.
Run after all methods complete, before evaluation.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

from breast_constants import BREAST_8_CELL_TYPES

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

DOT_REMAP = {
    "Cancer.Epithelial": "Cancer Epithelial",
    "Normal.Epithelial": "Normal Epithelial",
    "Plasma.cells": "Plasma cells",
    "T.cells": "T cells",
    "B.cells": "B-cells",
}

METHODS = ["RCTD", "Cell2Location", "Tangram", "CARD", "Seurat"]


def cleanup_method(method_dir: Path, n_regions: int = 5) -> None:
    """Normalize column names for one method's outputs."""
    for i in range(n_regions):
        # Try multiple output CSV patterns
        candidates = [
            method_dir / f"Xenium_region_{i}" / f"Xenium_region_{i}_cell_prop_finetuned_results.csv",
            method_dir / f"Xenium_region_{i}_cell_prop_finetuned_results.csv",
            method_dir / f"Xenium_region_{i}" / f"Xenium_region_{i}_deconv_predictions.csv",
        ]

        csv_path = None
        for c in candidates:
            if c.exists():
                csv_path = c
                break

        if csv_path is None:
            logger.warning(f"Missing output for region {i}")
            continue

        df = pd.read_csv(csv_path, index_col=0)

        renamed = {}
        for col in df.columns:
            new_col = col.strip()
            if new_col in DOT_REMAP:
                new_col = DOT_REMAP[new_col]
            elif "." in new_col:
                new_col = new_col.replace(".", " ")
            if new_col != col:
                renamed[col] = new_col

        if renamed:
            df = df.rename(columns=renamed)
            df.to_csv(csv_path)
            logger.info(f"  {csv_path.name}: renamed {renamed}")
        else:
            logger.info(f"  {csv_path.name}: no changes needed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--n-regions", type=int, default=5)
    args = parser.parse_args()

    base = Path(args.base_dir)
    for method in METHODS:
        method_dir = base / method / "output_breast"
        if method_dir.exists():
            logger.info(f"Cleaning {method}...")
            cleanup_method(method_dir, args.n_regions)
        else:
            logger.info(f"Skipping {method} (no output_breast/)")


if __name__ == "__main__":
    main()
