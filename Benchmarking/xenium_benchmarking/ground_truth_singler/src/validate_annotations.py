#!/usr/bin/env python
"""
Validate SingleR annotations and compute concordance with existing GT approaches.

Outputs:
- validation_report.json (marker expression, pruned rates, confidence)
- concordance_report.json (ARI, NMI, F1 vs protein-gated and RNA k-means)
"""

import argparse
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, f1_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

# Marker genes expected to be highly expressed in each type
MARKER_GENES = {
    "B cells": ["CD79A", "MS4A1"],
    "T cells": ["CD3D", "CD3E"],
    "Macrophages": ["CD68", "CD163"],
    "Fibroblasts": ["DCN", "ACTA2"],
    "Epithelial": ["EPCAM", "KRT8"],
    "Endothelial": ["PECAM1", "VWF"],
}


def compute_marker_expression(labels_df, counts_df, cell_types):
    """For each cell type, compute mean expression of its markers vs others."""
    results = {}
    for ct, markers in MARKER_GENES.items():
        if ct not in cell_types:
            continue
        ct_mask = labels_df["singler_label"] == ct
        for marker in markers:
            if marker not in counts_df.columns:
                continue
            in_type = counts_df.loc[ct_mask, marker].mean()
            out_type = counts_df.loc[~ct_mask, marker].mean()
            ratio = in_type / out_type if out_type > 0 else float("inf")
            results[f"{ct}_{marker}"] = {
                "in_type_mean": float(in_type),
                "out_type_mean": float(out_type),
                "enrichment_ratio": float(ratio),
            }
    return results


def compute_concordance(singler_labels, other_labels):
    """Compute ARI, NMI, per-type F1 between two label sets."""
    # Align by cell ID
    common = singler_labels.index.intersection(other_labels.index)
    if len(common) == 0:
        return {"error": "no overlapping cells"}

    s = singler_labels.loc[common]
    o = other_labels.loc[common]

    # Drop NaN / Unknown
    valid = s.notna() & o.notna() & (o != "Unknown")
    s = s[valid]
    o = o[valid]

    if len(s) == 0:
        return {"error": "no valid overlapping cells after filtering"}

    ari = adjusted_rand_score(s, o)
    nmi = normalized_mutual_info_score(s, o)

    # Per-type F1
    all_types = sorted(set(s) | set(o))
    f1s = {}
    for ct in all_types:
        s_binary = (s == ct).astype(int)
        o_binary = (o == ct).astype(int)
        f1s[ct] = float(f1_score(s_binary, o_binary, zero_division=0))

    return {
        "n_cells": int(len(s)),
        "ARI": float(ari),
        "NMI": float(nmi),
        "per_type_F1": f1s,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--singler-labels", required=True)
    parser.add_argument("--xenium-counts", required=True, help="Cell-level counts CSV")
    parser.add_argument("--protein-gt-labels", help="Existing protein-gated labels CSV (optional)")
    parser.add_argument("--rna-gt-labels", help="Existing RNA k-means labels CSV (optional)")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    labels_df = pd.read_csv(args.singler_labels)
    counts_df = pd.read_csv(args.xenium_counts, index_col=0)

    # Align counts index with labels
    labels_df = labels_df.set_index("cell_id")
    common_cells = labels_df.index.intersection(counts_df.index)
    labels_df = labels_df.loc[common_cells]
    counts_df = counts_df.loc[common_cells]

    # Validation: marker expression
    logger.info("Checking marker expression...")
    marker_results = compute_marker_expression(labels_df, counts_df, set(labels_df["singler_label"]))

    # Confidence stats
    confidence = {
        "median_score": float(labels_df["singler_score"].median()),
        "median_delta_next": float(labels_df["delta_next"].median()),
        "pruned_pct": float(100 * labels_df["pruned"].mean()),
        "per_type_counts": labels_df["singler_label"].value_counts().to_dict(),
    }

    validation = {"marker_expression": marker_results, "confidence": confidence}
    with open(output_dir / "validation_report.json", "w") as f:
        json.dump(validation, f, indent=2)
    logger.info(f"Validation report saved. Pruned: {confidence['pruned_pct']:.1f}%")

    # Concordance
    concordance = {}
    singler_series = labels_df["singler_label"]

    if args.protein_gt_labels:
        logger.info("Computing concordance with protein-gated GT...")
        protein_labels = pd.read_csv(args.protein_gt_labels, index_col=0)["cell_type"]
        # Map protein 7-type -> 6-type for comparison
        protein_6type = protein_labels.replace({"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"})
        concordance["vs_protein_gated"] = compute_concordance(singler_series, protein_6type)

    if args.rna_gt_labels:
        logger.info("Computing concordance with RNA k-means GT...")
        rna_labels = pd.read_csv(args.rna_gt_labels, index_col=0)["cell_type"]
        concordance["vs_rna_kmeans"] = compute_concordance(singler_series, rna_labels)

    with open(output_dir / "concordance_report.json", "w") as f:
        json.dump(concordance, f, indent=2)
    logger.info("Concordance report saved.")


if __name__ == "__main__":
    main()
