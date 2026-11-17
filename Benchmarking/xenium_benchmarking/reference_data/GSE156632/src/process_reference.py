"""
Process GSE156632 scRNA-seq reference data for Xenium benchmarking.

This script:
1. Loads all 12 CSV.gz files from GSE156632
2. Merges into single AnnData object with batch information
3. Performs QC filtering (min_genes, max_genes, mt%)
4. Normalizes (log1p, target_sum=1e4)
5. Integrates batches using Harmony
6. Annotates cell types via CellTypist majority voting over Leiden clusters
7. Exports in formats for Cell2Location, Tangram, RCTD, Seurat

Target cell types (must match protein-gated ground truth):
- B cells
- T cells
- Macrophages
- Fibroblasts
- Epithelial
- Endothelial

Usage:
    python process_reference.py --input-dir /path/to/GSE156632 --output-dir /path/to/processed
"""

import argparse
import gc
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Try to import harmony - will be available in the cell2location_env
try:
    import harmonypy

    HARMONY_AVAILABLE = True
except ImportError:
    HARMONY_AVAILABLE = False
    print("Warning: harmonypy not available. Will skip batch integration.")

try:
    import celltypist

    CELLTYPIST_AVAILABLE = True
except ImportError:
    CELLTYPIST_AVAILABLE = False
    print("Warning: celltypist not available. Cannot annotate cell types.")

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


# CellTypist Immune_All_Low.pkl label map: all 98 fine-grained types → 6 coarse types
CELLTYPIST_LABEL_MAP: Dict[str, str] = {
    # B lineage → B cells
    "Age-associated B cells": "B cells",
    "B cells": "B cells",
    "Cycling B cells": "B cells",
    "Follicular B cells": "B cells",
    "Germinal center B cells": "B cells",
    "Large pre-B cells": "B cells",
    "Memory B cells": "B cells",
    "Naive B cells": "B cells",
    "Plasma cells": "B cells",
    "Plasmablasts": "B cells",
    "Pre-pro-B cells": "B cells",
    "Pro-B cells": "B cells",
    "Proliferative germinal center B cells": "B cells",
    "Small pre-B cells": "B cells",
    "Transitional B cells": "B cells",
    # T / NK / ILC lineage → T cells
    "CD16+ NK cells": "T cells",
    "CD16- NK cells": "T cells",
    "CD8a/a": "T cells",
    "CD8a/b(entry)": "T cells",
    "CRTAM+ gamma-delta T cells": "T cells",
    "Cycling NK cells": "T cells",
    "Cycling T cells": "T cells",
    "Cycling gamma-delta T cells": "T cells",
    "Follicular helper T cells": "T cells",
    "ILC": "T cells",
    "ILC precursor": "T cells",
    "ILC1": "T cells",
    "ILC2": "T cells",
    "ILC3": "T cells",
    "MAIT cells": "T cells",
    "Memory CD4+ cytotoxic T cells": "T cells",
    "NK cells": "T cells",
    "NKT cells": "T cells",
    "Regulatory T cells": "T cells",
    "T(agonist)": "T cells",
    "Tcm/Naive cytotoxic T cells": "T cells",
    "Tcm/Naive helper T cells": "T cells",
    "Tem/Effector helper T cells": "T cells",
    "Tem/Effector helper T cells PD1+": "T cells",
    "Tem/Temra cytotoxic T cells": "T cells",
    "Tem/Trm cytotoxic T cells": "T cells",
    "Transitional NK": "T cells",
    "Treg(diff)": "T cells",
    "Trm cytotoxic T cells": "T cells",
    "Type 1 helper T cells": "T cells",
    "Type 17 helper T cells": "T cells",
    "gamma-delta T cells": "T cells",
    # Myeloid → Macrophages
    "Alveolar macrophages": "Macrophages",
    "Classical monocytes": "Macrophages",
    "Cycling DCs": "Macrophages",
    "Cycling monocytes": "Macrophages",
    "DC": "Macrophages",
    "DC precursor": "Macrophages",
    "DC1": "Macrophages",
    "DC2": "Macrophages",
    "DC3": "Macrophages",
    "Erythrophagocytic macrophages": "Macrophages",
    "Granulocytes": "Macrophages",
    "Hofbauer cells": "Macrophages",
    "Intermediate macrophages": "Macrophages",
    "Intestinal macrophages": "Macrophages",
    "Kidney-resident macrophages": "Macrophages",
    "Kupffer cells": "Macrophages",
    "MNP": "Macrophages",
    "Macrophages": "Macrophages",
    "Mast cells": "Macrophages",
    "Migratory DCs": "Macrophages",
    "Mono-mac": "Macrophages",
    "Monocyte precursor": "Macrophages",
    "Monocytes": "Macrophages",
    "Myelocytes": "Macrophages",
    "Neutrophil-myeloid progenitor": "Macrophages",
    "Neutrophils": "Macrophages",
    "Non-classical monocytes": "Macrophages",
    "Promyelocytes": "Macrophages",
    "Transitional DC": "Macrophages",
    "pDC": "Macrophages",
    "pDC precursor": "Macrophages",
    # Stromal
    "Fibroblasts": "Fibroblasts",
    "Epithelial cells": "Epithelial",
    "Endothelial cells": "Endothelial",
    # Hematopoietic progenitors → Unknown (rare in adult kidney)
    "CMP": "Unknown",
    "Double-negative thymocytes": "Unknown",
    "Double-positive thymocytes": "Unknown",
    "ELP": "Unknown",
    "ETP": "Unknown",
    "Early MK": "Unknown",
    "Early erythroid": "Unknown",
    "Early lymphoid/T lymphoid": "Unknown",
    "GMP": "Unknown",
    "HSC/MPP": "Unknown",
    "Late erythroid": "Unknown",
    "MEMP": "Unknown",
    "Megakaryocyte precursor": "Unknown",
    "Megakaryocyte-erythroid-mast cell progenitor": "Unknown",
    "Megakaryocytes/platelets": "Unknown",
    "Mid erythroid": "Unknown",
    # Erythrocytes → Unknown (RBC contamination)
    "Erythrocytes": "Unknown",
}

# Canonical markers for DE-based QC validation
CANONICAL_MARKERS: Dict[str, List[str]] = {
    "B cells": ["CD79A", "CD79B", "MS4A1"],
    "T cells": ["CD3D", "CD3E", "CD3G"],
    "Macrophages": ["CD68", "CD163", "CSF1R"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN"],
    "Epithelial": ["EPCAM", "KRT18", "KRT8"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
}


def load_csv_files(input_dir: Path) -> ad.AnnData:
    """
    Load all GSE156632 CSV.gz files and merge into single AnnData.

    Args:
        input_dir: Directory containing CSV.gz files

    Returns:
        Merged AnnData object with batch information
    """
    csv_files = sorted(input_dir.glob("GSM*.csv.gz"))
    logger.info(f"Found {len(csv_files)} CSV files")

    adatas = []
    for csv_file in csv_files:
        logger.info(f"Loading {csv_file.name}...")

        # Extract sample name from filename
        # Format: GSM4735364_RCC1n.csv.gz -> RCC1n
        sample_name = csv_file.stem.replace(".csv", "").split("_")[1]

        # Determine if tumor or normal
        tissue_type = "normal" if sample_name.endswith("n") else "tumor"

        # Load CSV (genes as rows, cells as columns - need to transpose)
        # The CSV has Gene_ID, Symbol, then cell barcodes
        # We want to use Symbol (gene names) as var_names
        df = pd.read_csv(csv_file, index_col=0)
        logger.info(f"  Shape before processing: {df.shape}")

        # Extract gene symbols from the 'Symbol' column (first column after index)
        gene_symbols = df["Symbol"].values if "Symbol" in df.columns else df.index.values

        # Drop the Symbol column if it exists (keep only count data)
        if "Symbol" in df.columns:
            df = df.drop(columns=["Symbol"])

        logger.info(f"  Shape after dropping Symbol column: {df.shape}")

        # Transpose: cells should be rows, genes should be columns
        df = df.T
        logger.info(f"  Shape after transpose: {df.shape}")

        # Create AnnData
        adata = ad.AnnData(X=df.values.astype(np.float32))
        adata.obs_names = df.index.astype(str)
        adata.var_names = pd.Index(gene_symbols).astype(str)

        # Make var_names unique (some gene symbols may be duplicated)
        adata.var_names_make_unique()

        # Add metadata
        adata.obs["sample"] = sample_name
        adata.obs["tissue_type"] = tissue_type
        adata.obs["batch"] = sample_name  # For batch correction

        logger.info(f"  Created AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")
        adatas.append(adata)

        # Clean up
        del df
        gc.collect()

    # Merge all samples
    logger.info("Merging all samples...")
    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.obs_names_make_unique()

    logger.info(f"Merged AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")
    logger.info(f"Samples: {adata.obs['sample'].unique().tolist()}")

    # Clean up
    del adatas
    gc.collect()

    return adata


def qc_filtering(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 5000,
    max_mt_pct: float = 20.0,
    min_cells: int = 3,
) -> ad.AnnData:
    """
    Quality control filtering.

    Args:
        adata: Input AnnData
        min_genes: Minimum genes per cell
        max_genes: Maximum genes per cell
        max_mt_pct: Maximum mitochondrial percentage
        min_cells: Minimum cells per gene

    Returns:
        Filtered AnnData
    """
    logger.info("Performing QC filtering...")
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    logger.info(f"  After min_genes={min_genes}: {adata.n_obs} cells")

    # Filter by max genes
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes, :].copy()
    logger.info(f"  After max_genes={max_genes}: {adata.n_obs} cells")

    # Filter by MT percentage
    adata = adata[adata.obs["pct_counts_mt"] < max_mt_pct, :].copy()
    logger.info(f"  After max_mt_pct={max_mt_pct}: {adata.n_obs} cells")

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)
    logger.info(f"  After min_cells={min_cells}: {adata.n_vars} genes")

    logger.info(f"QC complete: {n_cells_before} -> {adata.n_obs} cells, " f"{n_genes_before} -> {adata.n_vars} genes")

    return adata


def normalize_and_scale(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """
    Normalize, log-transform, and identify highly variable genes.

    Args:
        adata: QC-filtered AnnData
        target_sum: Target sum for normalization

    Returns:
        Normalized AnnData with raw counts in .raw
    """
    logger.info("Normalizing and scaling...")

    # Store raw counts
    adata.raw = adata.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Find highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        batch_key="batch",
        flavor="seurat_v3",
        subset=False,
    )

    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"  Identified {n_hvg} highly variable genes")

    return adata


def run_harmony(adata: ad.AnnData, batch_key: str = "batch") -> ad.AnnData:
    """
    Run Harmony batch integration.

    Args:
        adata: Normalized AnnData
        batch_key: Column in obs for batch information

    Returns:
        AnnData with harmony-corrected embedding
    """
    if not HARMONY_AVAILABLE:
        logger.warning("Harmony not available, skipping batch integration")
        return adata

    logger.info("Running Harmony batch integration...")

    # PCA first
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

    # Run Harmony
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch_key,
        max_iter_harmony=20,
    )

    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    # Compute neighbors and UMAP on harmony-corrected embedding
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)

    logger.info(f"  Found {adata.obs['leiden'].nunique()} clusters")

    return adata


def annotate_with_celltypist(
    adata: ad.AnnData,
    model_name: str = "Immune_All_Low.pkl",
    label_map: Dict[str, str] = CELLTYPIST_LABEL_MAP,
) -> ad.AnnData:
    """
    Annotate cell types using CellTypist majority voting over Leiden clusters.

    CellTypist requires log-normalized data (normalize_total + log1p). Since
    adata.X is z-scaled after Harmony, we re-normalize from raw counts for
    CellTypist prediction, then store labels back in the original adata.

    Args:
        adata: Clustered AnnData (must have 'leiden' in .obs and .raw with counts)
        model_name: CellTypist model to use (default: Immune_All_Low.pkl)
        label_map: Mapping from CellTypist fine labels to coarse types

    Returns:
        AnnData with celltypist_label, celltypist_conf, and cell_type in .obs
    """
    if not CELLTYPIST_AVAILABLE:
        raise ImportError("celltypist is required for annotation. " "Install with: pip install celltypist")

    if "leiden" not in adata.obs.columns:
        raise ValueError("adata.obs must contain 'leiden' clusters before annotation")

    logger.info("Annotating cell types with CellTypist...")
    logger.info(f"  Model: {model_name}")

    # Re-normalize from raw for CellTypist (needs log-normalized, not z-scaled)
    logger.info("  Re-normalizing from raw counts for CellTypist...")
    if adata.raw is not None:
        adata_ct = adata.raw.to_adata()
    else:
        adata_ct = adata.copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    # Load model
    logger.info(f"  Loading CellTypist model: {model_name}")
    model = celltypist.models.Model.load(model=model_name)

    # Run CellTypist with majority voting over Leiden clusters
    logger.info("  Running CellTypist annotation with majority voting...")
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=True,
        over_clustering=adata.obs["leiden"],
    )

    # Extract results
    result_adata = predictions.to_adata()

    # Store fine-grained labels and confidence
    adata.obs["celltypist_label"] = result_adata.obs["majority_voting"].values
    adata.obs["celltypist_conf"] = result_adata.obs["conf_score"].values

    # Map to coarse types
    unmapped = set()
    coarse_labels = []
    for label in adata.obs["celltypist_label"]:
        if label in label_map:
            coarse_labels.append(label_map[label])
        else:
            unmapped.add(label)
            coarse_labels.append("Unknown")

    if unmapped:
        logger.warning(f"  {len(unmapped)} unmapped CellTypist labels → 'Unknown': " f"{sorted(unmapped)}")

    adata.obs["cell_type"] = pd.Categorical(coarse_labels)

    # Log fine-grained distribution
    fine_counts = adata.obs["celltypist_label"].value_counts()
    logger.info(f"  CellTypist fine-grained types found: {len(fine_counts)}")
    for ct, count in fine_counts.head(20).items():
        pct = 100 * count / adata.n_obs
        logger.info(f"    {ct}: {count} ({pct:.1f}%)")
    if len(fine_counts) > 20:
        logger.info(f"    ... and {len(fine_counts) - 20} more types")

    # Log coarse distribution
    ct_counts = adata.obs["cell_type"].value_counts()
    logger.info("  Coarse cell type distribution:")
    for ct, count in ct_counts.items():
        pct = 100 * count / adata.n_obs
        logger.info(f"    {ct}: {count} ({pct:.1f}%)")

    # Clean up
    del adata_ct, result_adata
    gc.collect()

    return adata


def run_annotation_qc(
    adata: ad.AnnData,
    output_dir: Path,
    canonical_markers: Dict[str, List[str]] = CANONICAL_MARKERS,
) -> None:
    """
    Run QC validation on CellTypist annotation and write qc_report.txt.

    Checks:
    1. Composition table (cell counts + percentages per type)
    2. Low-confidence clusters (Leiden clusters with <60% majority-vote agreement)
    3. DE marker sanity check (canonical markers in top 50 DE per type)

    Args:
        adata: Annotated AnnData with cell_type, celltypist_label, leiden in .obs
        output_dir: Directory to write qc_report.txt
        canonical_markers: Dict mapping coarse type to expected DE marker genes
    """
    logger.info("Running annotation QC...")
    lines: List[str] = []
    lines.append("=" * 70)
    lines.append("CellTypist Annotation QC Report")
    lines.append("=" * 70)
    lines.append("")

    # 1. Composition table
    lines.append("1. Cell Type Composition")
    lines.append("-" * 40)
    ct_counts = adata.obs["cell_type"].value_counts()
    for ct, count in ct_counts.items():
        pct = 100 * count / adata.n_obs
        lines.append(f"  {ct:20s}  {count:>7d}  ({pct:5.1f}%)")
    lines.append(f"  {'TOTAL':20s}  {adata.n_obs:>7d}")
    lines.append("")

    # 2. Low-confidence clusters
    lines.append("2. Low-Confidence Leiden Clusters (<60% majority-vote agreement)")
    lines.append("-" * 40)
    if "celltypist_label" in adata.obs.columns and "leiden" in adata.obs.columns:
        low_conf_found = False
        for cluster in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
            mask = adata.obs["leiden"] == cluster
            cluster_labels = adata.obs.loc[mask, "celltypist_label"]
            majority_label = cluster_labels.mode().iloc[0]
            agreement = (cluster_labels == majority_label).mean()
            n_cells = mask.sum()
            if agreement < 0.60:
                low_conf_found = True
                coarse = CELLTYPIST_LABEL_MAP.get(majority_label, "Unknown")
                lines.append(
                    f"  Cluster {cluster:>3s}: {agreement:.1%} agreement "
                    f"({n_cells} cells, majority={majority_label} → {coarse})"
                )
        if not low_conf_found:
            lines.append("  None — all clusters have >=60% agreement.")
    else:
        lines.append("  Skipped — missing celltypist_label or leiden columns.")
    lines.append("")

    # 3. DE marker sanity check
    lines.append("3. DE Marker Sanity Check (canonical markers in top 50 DE per type)")
    lines.append("-" * 40)

    # Filter out Unknown for DE analysis
    known_mask = adata.obs["cell_type"] != "Unknown"
    if known_mask.sum() < adata.n_obs:
        lines.append(f"  Note: {adata.n_obs - known_mask.sum()} 'Unknown' cells excluded from DE.")

    adata_known = adata[known_mask].copy()

    # Re-normalize from raw for DE (need log-normalized, not z-scaled)
    if adata_known.raw is not None:
        adata_de = adata_known.raw.to_adata()
    else:
        adata_de = adata_known.copy()
    adata_de.obs["cell_type"] = adata_known.obs["cell_type"]
    sc.pp.normalize_total(adata_de, target_sum=1e4)
    sc.pp.log1p(adata_de)

    # Run rank_genes_groups
    try:
        sc.tl.rank_genes_groups(adata_de, groupby="cell_type", method="wilcoxon")

        for ct, markers in canonical_markers.items():
            if ct not in adata_de.obs["cell_type"].cat.categories:
                lines.append(f"  {ct}: ABSENT from annotations — no cells assigned")
                continue

            # Get top 50 DE genes for this type
            de_result = sc.get.rank_genes_groups_df(adata_de, group=ct)
            top50 = set(de_result.head(50)["names"].tolist())

            found = [m for m in markers if m in top50]
            missing = [m for m in markers if m not in top50]

            status = "PASS" if len(found) >= 2 else "WARN"
            lines.append(f"  {ct}: {status} — {len(found)}/{len(markers)} canonical markers " f"in top 50 DE")
            if found:
                lines.append(f"    Found: {', '.join(found)}")
            if missing:
                lines.append(f"    Missing: {', '.join(missing)}")
    except Exception as e:
        lines.append(f"  DE analysis failed: {e}")

    lines.append("")
    lines.append("=" * 70)

    # Write report
    report_path = output_dir / "qc_report.txt"
    report_text = "\n".join(lines)
    report_path.write_text(report_text)
    logger.info(f"  QC report written to {report_path}")

    # Also log it
    for line in lines:
        logger.info(f"  {line}")

    # Clean up
    del adata_known, adata_de
    gc.collect()


def export_for_cell2location(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for Cell2Location reference training."""
    logger.info("Exporting for Cell2Location...")

    c2l_dir = output_dir / "cell2location"
    c2l_dir.mkdir(parents=True, exist_ok=True)

    # Cell2Location needs raw counts + cell type labels
    adata_export = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    adata_export.obs["cell_type"] = adata.obs["cell_type"]
    adata_export.obs["sample"] = adata.obs["sample"]
    adata_export.obs["batch"] = adata.obs["batch"]

    adata_export.write_h5ad(c2l_dir / "reference.h5ad")
    logger.info(f"  Saved to {c2l_dir / 'reference.h5ad'}")


def export_for_tangram(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for Tangram.

    Tangram expects log-normalized data (normalize_total + log1p), NOT z-scaled.
    The main adata.X is z-scaled after Harmony, so we need to re-normalize from raw.
    """
    logger.info("Exporting for Tangram...")

    tg_dir = output_dir / "tangram"
    tg_dir.mkdir(parents=True, exist_ok=True)

    # Tangram needs log-normalized counts (NOT z-scaled!)
    # Start from raw counts and apply normalize_total + log1p
    if adata.raw is not None:
        adata_export = adata.raw.to_adata()
    else:
        adata_export = adata.copy()

    # Copy over cell type annotations
    adata_export.obs["cell_type"] = adata.obs["cell_type"]
    adata_export.obs["sample"] = adata.obs["sample"]
    adata_export.obs["batch"] = adata.obs["batch"]

    # Apply log-normalization (matching simulated benchmark preprocessing)
    sc.pp.normalize_total(adata_export)
    sc.pp.log1p(adata_export)

    logger.info(f"  X range after log-norm: [{adata_export.X.min():.2f}, {adata_export.X.max():.2f}]")

    adata_export.write_h5ad(tg_dir / "reference.h5ad")
    logger.info(f"  Saved to {tg_dir / 'reference.h5ad'}")


def export_for_rctd(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for RCTD (R)."""
    logger.info("Exporting for RCTD...")

    rctd_dir = output_dir / "rctd"
    rctd_dir.mkdir(parents=True, exist_ok=True)

    # RCTD needs:
    # 1. counts matrix (genes x cells) as CSV
    # 2. cell types as CSV

    # Get raw counts
    if adata.raw is not None:
        counts = adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X
        gene_names = adata.raw.var_names
    else:
        counts = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        gene_names = adata.var_names

    # Save counts (genes x cells - transpose)
    counts_df = pd.DataFrame(
        counts.T,
        index=gene_names,
        columns=adata.obs_names,
    )
    counts_df.to_csv(rctd_dir / "counts.csv")
    logger.info(f"  Saved counts: {counts_df.shape}")

    # Save cell types
    cell_types = adata.obs[["cell_type"]].copy()
    cell_types.to_csv(rctd_dir / "cell_types.csv")
    logger.info(f"  Saved cell types: {cell_types.shape}")

    # Also save as h5ad for convenience
    adata_export = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    adata_export.obs["cell_type"] = adata.obs["cell_type"]
    adata_export.write_h5ad(rctd_dir / "reference.h5ad")


def export_for_seurat(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for Seurat."""
    logger.info("Exporting for Seurat...")

    seurat_dir = output_dir / "seurat"
    seurat_dir.mkdir(parents=True, exist_ok=True)

    # Seurat needs same format as RCTD
    # Get raw counts
    if adata.raw is not None:
        counts = adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X
        gene_names = adata.raw.var_names
    else:
        counts = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        gene_names = adata.var_names

    # Save counts (genes x cells - transpose)
    counts_df = pd.DataFrame(
        counts.T,
        index=gene_names,
        columns=adata.obs_names,
    )
    counts_df.to_csv(seurat_dir / "counts.csv")
    logger.info(f"  Saved counts: {counts_df.shape}")

    # Save cell types
    cell_types = adata.obs[["cell_type"]].copy()
    cell_types.to_csv(seurat_dir / "cell_types.csv")
    logger.info(f"  Saved cell types: {cell_types.shape}")

    # Also save as h5ad for convenience
    adata_export = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    adata_export.obs["cell_type"] = adata.obs["cell_type"]
    adata_export.write_h5ad(seurat_dir / "reference.h5ad")


def main():
    parser = argparse.ArgumentParser(description="Process GSE156632 reference for Xenium benchmarking")
    parser.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory containing GSE156632 CSV.gz files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for processed data",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Minimum genes per cell (default: 200)",
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        default=5000,
        help="Maximum genes per cell (default: 5000)",
    )
    parser.add_argument(
        "--max-mt-pct",
        type=float,
        default=20.0,
        help="Maximum mitochondrial percentage (default: 20.0)",
    )
    parser.add_argument(
        "--skip-harmony",
        action="store_true",
        help="Skip Harmony batch integration",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="Immune_All_Low.pkl",
        help="CellTypist model name (default: Immune_All_Low.pkl)",
    )
    parser.add_argument(
        "--model-cache-dir",
        type=str,
        default=None,
        help="Directory to cache CellTypist models (sets CELLTYPIST_FOLDER env var)",
    )

    args = parser.parse_args()

    # Set CellTypist cache directory if specified
    if args.model_cache_dir is not None:
        os.environ["CELLTYPIST_FOLDER"] = args.model_cache_dir

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("GSE156632 Reference Processing")
    logger.info("=" * 60)
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")

    # Step 1: Load data
    logger.info("\n" + "=" * 60)
    logger.info("Step 1: Loading CSV files")
    logger.info("=" * 60)
    adata = load_csv_files(input_dir)

    # Step 2: QC filtering
    logger.info("\n" + "=" * 60)
    logger.info("Step 2: QC Filtering")
    logger.info("=" * 60)
    adata = qc_filtering(
        adata,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_mt_pct=args.max_mt_pct,
    )

    # Step 3: Normalize
    logger.info("\n" + "=" * 60)
    logger.info("Step 3: Normalization")
    logger.info("=" * 60)
    adata = normalize_and_scale(adata)

    # Step 4: Harmony batch integration
    if not args.skip_harmony:
        logger.info("\n" + "=" * 60)
        logger.info("Step 4: Harmony Batch Integration")
        logger.info("=" * 60)
        adata = run_harmony(adata)
    else:
        logger.info("\n" + "=" * 60)
        logger.info("Step 4: Skipping Harmony (--skip-harmony)")
        logger.info("=" * 60)
        # Still need to do PCA, neighbors, clustering
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.8)

    # Step 5: Cell type annotation via CellTypist
    logger.info("\n" + "=" * 60)
    logger.info("Step 5: Cell Type Annotation (CellTypist)")
    logger.info("=" * 60)
    adata = annotate_with_celltypist(adata, model_name=args.model)

    # Step 5b: Annotation QC
    logger.info("\n" + "=" * 60)
    logger.info("Step 5b: Annotation QC Validation")
    logger.info("=" * 60)
    run_annotation_qc(adata, output_dir)

    # Step 6: Save combined dataset
    logger.info("\n" + "=" * 60)
    logger.info("Step 6: Saving Combined Dataset")
    logger.info("=" * 60)
    adata.write_h5ad(output_dir / "GSE156632_combined.h5ad")
    logger.info(f"Saved combined dataset to {output_dir / 'GSE156632_combined.h5ad'}")

    # Step 7: Export for each method
    logger.info("\n" + "=" * 60)
    logger.info("Step 7: Exporting for Deconvolution Methods")
    logger.info("=" * 60)

    export_for_cell2location(adata, output_dir)
    export_for_tangram(adata, output_dir)
    export_for_rctd(adata, output_dir)
    export_for_seurat(adata, output_dir)

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Processing Complete!")
    logger.info("=" * 60)
    logger.info(f"Final dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    logger.info(f"Cell types: {adata.obs['cell_type'].unique().tolist()}")
    logger.info(f"Output files:")
    logger.info(f"  - {output_dir / 'GSE156632_combined.h5ad'}")
    logger.info(f"  - {output_dir / 'cell2location' / 'reference.h5ad'}")
    logger.info(f"  - {output_dir / 'tangram' / 'reference.h5ad'}")
    logger.info(f"  - {output_dir / 'rctd' / 'counts.csv'}")
    logger.info(f"  - {output_dir / 'seurat' / 'counts.csv'}")


if __name__ == "__main__":
    main()
