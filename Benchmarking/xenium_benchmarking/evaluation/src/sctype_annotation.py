#!/usr/bin/env python3
"""
ScType-based per-cell RNA annotation for Xenium single-cell data.

Uses marker-based cell type scoring on 405 Xenium genes to produce
independent RNA-based cell annotations for validating Module 3 protein gating.

ScType assigns cell types using pre-defined marker databases weighted by
specificity. This provides an orthogonal annotation to protein-based gating.

Usage:
    python sctype_annotation.py --output_dir /path/to/output [--region REGION_ID]
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.preprocessing import scale

# Add paths
REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/src"))

from load_xenium import load_xenium_data, split_gex_protein

logger = logging.getLogger(__name__)

# Xenium data location
XENIUM_DATA_DIR = Path(
    "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/"
    "Xenium_RNA_Proteomic_RenalCellCarcinoma"
)

# ScType-style marker database for kidney/immune tissue
# Curated for the 405-gene Xenium panel.
# Format: {cell_type: {"positive": [markers], "negative": [markers]}}
# Only genes present in the 405-gene panel are included.
SCTYPE_MARKER_DB = {
    "T cells": {
        "positive": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
        "negative": ["CD19", "MS4A1", "CD14", "CD68"],
    },
    "CD4+ T cells": {
        "positive": ["CD3D", "CD3E", "CD4", "IL7R", "FOXP3", "IL2RA"],
        "negative": ["CD8A", "CD8B", "NKG7", "GNLY"],
    },
    "CD8+ T cells": {
        "positive": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "GZMA", "PRF1", "NKG7"],
        "negative": ["CD4", "FOXP3"],
    },
    "Exhausted T cells": {
        "positive": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "CD8A", "CD3E", "TOX", "CTLA4"],
        "negative": ["IL7R", "TCF7"],
    },
    "Regulatory T cells": {
        "positive": ["FOXP3", "IL2RA", "CTLA4", "CD3E", "CD4"],
        "negative": ["CD8A", "NKG7"],
    },
    "B cells": {
        "positive": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"],
        "negative": ["CD3D", "CD14", "CD68"],
    },
    "Plasma cells": {
        "positive": ["SDC1", "MZB1", "JCHAIN", "XBP1", "IGHG1"],
        "negative": ["MS4A1", "CD19", "PAX5"],
    },
    "Macrophages": {
        "positive": ["CD68", "CD14", "CD163", "CSF1R", "MSR1", "MARCO"],
        "negative": ["CD3D", "MS4A1", "EPCAM"],
    },
    "M2 Macrophages": {
        "positive": ["CD163", "MSR1", "MRC1", "CD68"],
        "negative": ["CD3D", "MS4A1"],
    },
    "Dendritic cells": {
        "positive": ["ITGAX", "HLA-DRA", "HLA-DRB1", "CD1C", "CLEC10A"],
        "negative": ["CD3D", "CD14", "MS4A1"],
    },
    "NK cells": {
        "positive": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "PRF1"],
        "negative": ["CD3D", "CD3E", "CD19"],
    },
    "Epithelial": {
        "positive": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
        "negative": ["VIM", "CD45", "PTPRC"],
    },
    "Clear cell RCC": {
        "positive": ["CA9", "PAX8", "EPCAM", "KRT8", "KRT18"],
        "negative": ["VIM", "CD3D", "CD68"],
    },
    "Papillary RCC": {
        "positive": ["PAX8", "EPCAM", "MET", "FH"],
        "negative": ["CA9"],
    },
    "Endothelial": {
        "positive": ["PECAM1", "VWF", "CDH5", "ERG", "FLT1", "KDR"],
        "negative": ["EPCAM", "CD3D", "CD68"],
    },
    "Fibroblasts": {
        "positive": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2", "VIM"],
        "negative": ["EPCAM", "CD3D", "CD68", "PECAM1"],
    },
    "Proliferating": {
        "positive": ["MKI67", "TOP2A", "PCNA", "STMN1"],
        "negative": [],
    },
}

# Broad lineage mapping for Tier 1 concordance
LINEAGE_MAP = {
    "T cells": "Immune",
    "CD4+ T cells": "Immune",
    "CD8+ T cells": "Immune",
    "Exhausted T cells": "Immune",
    "Regulatory T cells": "Immune",
    "B cells": "Immune",
    "Plasma cells": "Immune",
    "Macrophages": "Immune",
    "M2 Macrophages": "Immune",
    "Dendritic cells": "Immune",
    "NK cells": "Immune",
    "Epithelial": "Epithelial",
    "Clear cell RCC": "Epithelial",
    "Papillary RCC": "Epithelial",
    "Endothelial": "Stromal",
    "Fibroblasts": "Stromal",
    "Proliferating": "Other",
}


def filter_markers_to_panel(
    marker_db: Dict[str, Dict[str, List[str]]],
    available_genes: List[str],
) -> Dict[str, Dict[str, List[str]]]:
    """
    Filter marker database to only include genes present in the panel.

    Args:
        marker_db: ScType marker database.
        available_genes: List of gene names in the Xenium panel.

    Returns:
        Filtered marker database.
    """
    gene_set = set(available_genes)
    filtered = {}
    for cell_type, markers in marker_db.items():
        pos = [g for g in markers["positive"] if g in gene_set]
        neg = [g for g in markers["negative"] if g in gene_set]
        if len(pos) > 0:  # Need at least one positive marker
            filtered[cell_type] = {"positive": pos, "negative": neg}
            logger.info(
                f"  {cell_type}: {len(pos)}/{len(markers['positive'])} positive, "
                f"{len(neg)}/{len(markers['negative'])} negative markers available"
            )
        else:
            logger.warning(f"  {cell_type}: No positive markers in panel, skipping")
    return filtered


def sctype_score(
    adata: sc.AnnData,
    marker_db: Dict[str, Dict[str, List[str]]],
    scale_data: bool = True,
) -> pd.DataFrame:
    """
    Compute ScType-style cell type scores for each cell.

    For each cell, the score for a cell type is the mean z-scored expression
    of positive markers minus the mean z-scored expression of negative markers.

    Args:
        adata: AnnData with gene expression (preprocessed, log-normalized).
        marker_db: Filtered marker database.
        scale_data: Whether to z-score normalize genes before scoring.

    Returns:
        DataFrame (n_cells, n_cell_types) with ScType scores.
    """
    # Get dense matrix
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    gene_names = list(adata.var_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Z-score normalize per gene
    if scale_data:
        X_scaled = scale(X, axis=0)
    else:
        X_scaled = X.copy()

    n_cells = X_scaled.shape[0]
    cell_types = list(marker_db.keys())
    scores = np.zeros((n_cells, len(cell_types)), dtype=np.float64)

    for ct_idx, cell_type in enumerate(cell_types):
        pos_genes = marker_db[cell_type]["positive"]
        neg_genes = marker_db[cell_type]["negative"]

        # Positive marker score: mean of z-scored expression
        pos_indices = [gene_to_idx[g] for g in pos_genes if g in gene_to_idx]
        if pos_indices:
            pos_score = X_scaled[:, pos_indices].mean(axis=1)
        else:
            pos_score = np.zeros(n_cells)

        # Negative marker penalty: mean of z-scored expression (subtracted)
        neg_indices = [gene_to_idx[g] for g in neg_genes if g in gene_to_idx]
        if neg_indices:
            neg_score = X_scaled[:, neg_indices].mean(axis=1)
        else:
            neg_score = np.zeros(n_cells)

        # ScType score = positive - negative
        scores[:, ct_idx] = pos_score - neg_score

    return pd.DataFrame(scores, columns=cell_types, index=adata.obs_names)


def assign_cell_types(
    scores: pd.DataFrame,
    confidence_threshold: float = 0.0,
) -> pd.Series:
    """
    Assign each cell to the highest-scoring cell type.

    Args:
        scores: ScType scores DataFrame (n_cells, n_cell_types).
        confidence_threshold: Minimum score to assign (below = Unassigned).

    Returns:
        Series with cell type assignments.
    """
    best_type_idx = scores.values.argmax(axis=1)
    best_score = scores.values.max(axis=1)
    cell_types = scores.columns.tolist()

    assignments = pd.Series(
        [cell_types[i] if best_score[j] > confidence_threshold else "Unassigned"
         for j, i in enumerate(best_type_idx)],
        index=scores.index,
        name="sctype_annotation",
    )
    return assignments


def main(output_dir: str, verbose: bool = False):
    """Run ScType annotation on Xenium single-cell data."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading Xenium single-cell data...")
    adata = load_xenium_data(XENIUM_DATA_DIR)

    logger.info("Splitting into GEX and protein...")
    adata_gex, adata_protein = split_gex_protein(adata)

    logger.info(f"GEX matrix: {adata_gex.shape[0]} cells x {adata_gex.shape[1]} genes")
    logger.info(f"Protein matrix: {adata_protein.shape[0]} cells x {adata_protein.shape[1]} proteins")

    # Preprocess GEX
    logger.info("Preprocessing gene expression...")
    adata_gex_proc = adata_gex.copy()
    sc.pp.normalize_total(adata_gex_proc, target_sum=1e4)
    sc.pp.log1p(adata_gex_proc)

    # Filter marker database to available genes
    available_genes = list(adata_gex_proc.var_names)
    logger.info(f"Available genes in panel: {len(available_genes)}")
    logger.info("Filtering marker database to panel genes...")
    filtered_db = filter_markers_to_panel(SCTYPE_MARKER_DB, available_genes)
    logger.info(f"Cell types with markers: {len(filtered_db)}")

    # Save filtered database
    with open(output_dir / "filtered_marker_db.json", "w") as f:
        json.dump(filtered_db, f, indent=2)

    # Compute ScType scores
    logger.info("Computing ScType scores...")
    scores = sctype_score(adata_gex_proc, filtered_db)

    # Assign cell types
    logger.info("Assigning cell types...")
    assignments = assign_cell_types(scores, confidence_threshold=0.0)

    # Map to broad lineages
    lineage_assignments = assignments.map(
        lambda x: LINEAGE_MAP.get(x, "Other") if x != "Unassigned" else "Unassigned"
    )

    # Save results
    logger.info("Saving results...")
    results_df = pd.DataFrame({
        "sctype_annotation": assignments,
        "sctype_lineage": lineage_assignments,
    })
    # Add scores
    for col in scores.columns:
        results_df[f"score_{col}"] = scores[col].values

    results_df.to_csv(output_dir / "sctype_annotations.csv")
    scores.to_csv(output_dir / "sctype_scores.csv")

    # Summary statistics
    type_counts = assignments.value_counts()
    lineage_counts = lineage_assignments.value_counts()

    summary = {
        "n_cells": int(len(assignments)),
        "n_cell_types": int(len(filtered_db)),
        "n_genes_in_panel": len(available_genes),
        "type_counts": type_counts.to_dict(),
        "lineage_counts": lineage_counts.to_dict(),
        "unassigned_rate": float((assignments == "Unassigned").mean()),
    }

    with open(output_dir / "sctype_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Print summary
    print("\n" + "=" * 60)
    print("SCTYPE ANNOTATION SUMMARY")
    print("=" * 60)
    print(f"\nTotal cells: {len(assignments):,}")
    print(f"Cell types annotated: {len(filtered_db)}")
    print(f"Unassigned: {(assignments == 'Unassigned').sum():,} "
          f"({100 * (assignments == 'Unassigned').mean():.1f}%)")
    print("\nPer-type counts:")
    for ct, count in type_counts.items():
        pct = 100.0 * count / len(assignments)
        print(f"  {ct:25s}: {count:>7,} ({pct:5.1f}%)")
    print("\nLineage-level counts:")
    for lineage, count in lineage_counts.items():
        pct = 100.0 * count / len(assignments)
        print(f"  {lineage:15s}: {count:>7,} ({pct:5.1f}%)")
    print("=" * 60)

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ScType per-cell RNA annotation")
    parser.add_argument(
        "--output_dir",
        type=str,
        default=str(
            Path(__file__).resolve().parents[2]
            / "evaluation"
            / "results"
            / "sctype_annotation"
        ),
        help="Output directory",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logger.setLevel(logging.INFO)

    main(args.output_dir, args.verbose)
