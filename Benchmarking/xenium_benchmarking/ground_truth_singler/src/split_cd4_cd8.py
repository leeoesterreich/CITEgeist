#!/usr/bin/env python
"""
Split SingleR 'T cells' into CD4+ and CD8+ using Xenium protein expression.

Loads:
- SingleR labels (singler_all_cells.csv)
- Xenium protein expression (xenium_cell_proteins.csv)

Outputs:
- singler_labels_7type.csv (all cells with CD4/CD8 split)
- cd4_cd8_split_report.json (ambiguity stats)
"""

import argparse
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.mixture import GaussianMixture

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)


def gmm_threshold(values: np.ndarray, n_components: int = 3) -> float:
    """
    N-component GMM to find background/signal threshold.
    Returns the midpoint between the lowest two component means
    (background vs first signal population).
    """
    nonzero = values[values > 0].reshape(-1, 1)
    if len(nonzero) < 50:
        return np.percentile(values[values > 0], 50) if (values > 0).sum() > 0 else 0
    gmm = GaussianMixture(n_components=n_components, random_state=42).fit(nonzero)
    means = sorted(gmm.means_.flatten())
    return (means[0] + means[1]) / 2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--singler-labels", required=True, help="Path to singler_all_cells.csv")
    parser.add_argument("--protein-csv", required=True, help="Path to Xenium protein expression CSV (cells x proteins)")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load SingleR labels
    labels_df = pd.read_csv(args.singler_labels)
    logger.info(f"Loaded {len(labels_df)} cell labels")

    # Load protein expression
    protein_df = pd.read_csv(args.protein_csv, index_col=0)
    logger.info(f"Loaded protein: {protein_df.shape}")

    # Find T cells
    t_mask = labels_df["singler_label"] == "T cells"
    n_tcells = t_mask.sum()
    logger.info(f"T cells to split: {n_tcells}")

    if n_tcells == 0:
        logger.warning("No T cells found, saving 6-type labels as 7-type")
        labels_df.to_csv(output_dir / "singler_labels_7type.csv", index=False)
        return

    # Get CD4 and CD8A protein for T cells
    t_cell_ids = labels_df.loc[t_mask, "cell_id"].values
    cd4_vals = protein_df.loc[t_cell_ids, "CD4"].values.astype(float)
    cd8a_vals = protein_df.loc[t_cell_ids, "CD8A"].values.astype(float)

    # 3-component GMM per marker to find background/signal threshold
    cd4_thresh = gmm_threshold(cd4_vals)
    cd8a_thresh = gmm_threshold(cd8a_vals)
    logger.info(f"3-GMM thresholds: CD4={cd4_thresh:.1f}, CD8A={cd8a_thresh:.1f}")

    # Classify confident calls
    cd4_above = cd4_vals > cd4_thresh
    cd8_above = cd8a_vals > cd8a_thresh

    cd4_only = cd4_above & ~cd8_above
    cd8_only = cd8_above & ~cd4_above
    both_above = cd4_above & cd8_above
    neither = ~cd4_above & ~cd8_above

    logger.info(
        f"Confident: CD4+ only={cd4_only.sum()}, CD8+ only={cd8_only.sum()}, "
        f"both={both_above.sum()}, neither={neither.sum()} ({100*neither.sum()/n_tcells:.1f}%)"
    )

    # Scale-corrected binarization for double-positive and neither cells
    # Normalize each marker by its threshold so they're on comparable scales
    cd4_scaled = cd4_vals / max(cd4_thresh, 1.0)
    cd8_scaled = cd8a_vals / max(cd8a_thresh, 1.0)

    # Final assignment: confident calls take priority, then scale-corrected binarize the rest
    cd4_pos = cd4_only | (both_above & (cd4_scaled > cd8_scaled)) | (neither & (cd4_scaled > cd8_scaled))
    cd8_pos = ~cd4_pos

    n_cd4 = cd4_pos.sum()
    n_cd8 = cd8_pos.sum()
    logger.info(f"Final: CD4+ {n_cd4} ({100*n_cd4/n_tcells:.1f}%), CD8+ {n_cd8} ({100*n_cd8/n_tcells:.1f}%)")

    # Apply split — all T cells assigned
    labels_7type = labels_df.copy()
    t_indices = labels_df.index[t_mask]
    labels_7type.loc[t_indices[cd4_pos], "singler_label"] = "CD4+ T cells"
    labels_7type.loc[t_indices[cd8_pos], "singler_label"] = "CD8+ T cells"

    labels_7type.to_csv(output_dir / "singler_labels_7type.csv", index=False)

    # Report
    report = {
        "total_t_cells": int(n_tcells),
        "cd4_positive": int(n_cd4),
        "cd8_positive": int(n_cd8),
        "cd4_confident_only": int(cd4_only.sum()),
        "cd8_confident_only": int(cd8_only.sum()),
        "double_positive": int(both_above.sum()),
        "neither_confident": int(neither.sum()),
        "method": "gmm3_threshold_then_scale_corrected_binarize",
        "cd4_threshold": float(cd4_thresh),
        "cd8a_threshold": float(cd8a_thresh),
        "cd4_pct": float(100 * n_cd4 / n_tcells),
        "cd8_pct": float(100 * n_cd8 / n_tcells),
    }
    with open(output_dir / "cd4_cd8_split_report.json", "w") as f:
        json.dump(report, f, indent=2)

    logger.info("Done.")


if __name__ == "__main__":
    main()
