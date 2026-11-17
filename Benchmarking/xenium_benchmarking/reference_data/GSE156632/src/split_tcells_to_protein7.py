#!/usr/bin/env python
"""
Derive the 7-type Xenium reference from the authoritative 6-type annotation.

Splits "T cells" -> "CD4+ T cells" / "CD8+ T cells" using discriminating
marker scoring on raw counts. All other cell types are copied verbatim,
ensuring the only difference between 6-type and 7-type is T cell labeling.

Also writes cell_types_with_sample.csv for both 6-type and 7-type rctd
directories, adding the sample column needed by CARD's sampleInfo parameter.

Usage:
    python split_tcells_to_protein7.py
    (Run via SLURM: sbatch run_split_tcells.sh)
"""

import argparse
import logging
import sys
from typing import Optional
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Import export functions from process_reference.py (same directory)
sys.path.insert(0, str(Path(__file__).parent))
from process_reference import (
    export_for_cell2location,
    export_for_rctd,
    export_for_seurat,
    export_for_tangram,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent
COMBINED_H5AD = BASE_DIR / "processed" / "GSE156632_combined.h5ad"
OUTPUT_6TYPE = BASE_DIR / "processed"
OUTPUT_7TYPE = BASE_DIR / "processed_protein7"

# Discriminating markers only — no shared CD3 genes to avoid score ties.
# CD4 markers: expressed in CD4+ helper T cells.
# CD8 markers: expressed in CD8+ cytotoxic T cells.
CD4_MARKERS = ["CD4", "IL7R", "TCF7"]
CD8_MARKERS = ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"]


def split_t_cells(adata: sc.AnnData) -> pd.Series:
    """
    Score T cells for CD4 vs CD8 marker expression and assign subtypes.

    Uses raw counts (adata.raw) for scoring. Non-T cell labels are unchanged.

    Returns:
        pd.Series with updated cell_type labels (same index as adata.obs)
    """
    logger.info("Splitting T cells into CD4+ and CD8+ subtypes...")

    # Use raw counts for scoring (pre-normalization, unbiased)
    if adata.raw is not None:
        expr = pd.DataFrame(
            adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X,
            index=adata.obs_names,
            columns=adata.raw.var_names,
        )
    else:
        expr = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            index=adata.obs_names,
            columns=adata.var_names,
        )

    t_mask = adata.obs["cell_type"] == "T cells"
    n_t = t_mask.sum()
    logger.info(f"  T cells to split: {n_t}")

    cd4_present = [m for m in CD4_MARKERS if m in expr.columns]
    cd8_present = [m for m in CD8_MARKERS if m in expr.columns]
    logger.info(f"  CD4 markers found: {cd4_present}")
    logger.info(f"  CD8 markers found: {cd8_present}")

    t_expr = expr.loc[t_mask]
    cd4_score = t_expr[cd4_present].mean(axis=1) if cd4_present else pd.Series(0.0, index=t_expr.index)
    cd8_score = t_expr[cd8_present].mean(axis=1) if cd8_present else pd.Series(0.0, index=t_expr.index)

    # Convert to plain str to allow adding new categories (CD4+/CD8+ T cells)
    new_labels = adata.obs["cell_type"].astype(str).copy()
    new_labels.loc[t_mask] = np.where(
        cd4_score.values >= cd8_score.values,
        "CD4+ T cells",
        "CD8+ T cells",
    )

    n_cd4 = (new_labels == "CD4+ T cells").sum()
    n_cd8 = (new_labels == "CD8+ T cells").sum()
    assert n_cd4 + n_cd8 == n_t, f"T cell count mismatch: {n_cd4} + {n_cd8} != {n_t}"
    logger.info(f"  -> CD4+ T cells: {n_cd4}, CD8+ T cells: {n_cd8}")

    return new_labels


def validate_7type_reference(labels_6: pd.Series, labels_7: pd.Series) -> None:
    """Assert non-T cell type counts are identical between 6-type and 7-type."""
    non_t_types = [ct for ct in labels_6.unique() if ct != "T cells"]
    for ct in non_t_types:
        count_6 = (labels_6 == ct).sum()
        count_7 = (labels_7 == ct).sum()
        assert count_6 == count_7, f"Non-T cell type mismatch for '{ct}': 6-type={count_6}, 7-type={count_7}"
    logger.info("Validation passed: all non-T cell type counts identical between 6-type and 7-type")


def write_cell_types_with_sample(adata: sc.AnnData, rctd_dir: Path) -> None:
    """
    Write cell_types_with_sample.csv to an rctd directory.

    Adds the 'sample' column to the cell type annotation file so CARD can
    use per-sample reference profiles via its sampleInfo parameter.
    """
    df = adata.obs[["cell_type", "sample"]].copy()
    path = rctd_dir / "cell_types_with_sample.csv"
    df.to_csv(path)
    samples = sorted(df["sample"].unique().tolist())
    logger.info(f"  Wrote {path.name}: {len(df)} cells, {len(samples)} samples: {samples}")


def export_labels_only(adata: sc.AnnData, out_dir: Path, tool: str, counts_src: Path) -> None:
    """Write cell_types.csv for rctd/seurat and symlink the (label-free) counts.csv.

    counts.csv is identical across label sets, so it is linked to the existing file
    instead of re-written (4.9 GB). reference.h5ad is a convenience copy and is skipped.
    """
    d = out_dir / tool
    d.mkdir(parents=True, exist_ok=True)
    link = d / "counts.csv"
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise FileExistsError(f"refusing to replace real file with symlink: {link}")
    link.symlink_to(Path(counts_src).resolve())
    adata.obs[["cell_type"]].to_csv(d / "cell_types.csv")
    logger.info(f"  {tool}: cell_types.csv written, counts.csv -> {link.resolve()}")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--combined-h5ad", type=Path, default=COMBINED_H5AD)
    ap.add_argument("--labels-csv", type=Path, default=None,
                    help="Overlay obs['cell_type'] from this CSV (index=barcode, column cell_type)")
    ap.add_argument("--out-6type", type=Path, default=OUTPUT_6TYPE)
    ap.add_argument("--out-7type", type=Path, default=OUTPUT_7TYPE)
    ap.add_argument("--link-counts-from", type=Path, default=None,
                    help="Existing counts.csv to symlink (skips re-writing rctd/seurat counts)")
    return ap


def main(argv: Optional[list] = None) -> None:
    args = build_parser().parse_args(argv)
    combined_h5ad, out6, out7 = args.combined_h5ad, args.out_6type, args.out_7type
    if not combined_h5ad.exists():
        raise FileNotFoundError(f"Authoritative combined h5ad not found: {combined_h5ad}")

    logger.info(f"Loading authoritative 6-type annotation from {combined_h5ad}...")
    adata = sc.read_h5ad(combined_h5ad)
    if args.labels_csv is not None:
        from build_march_reference import relabel_with_march

        adata = relabel_with_march(adata, args.labels_csv)
    logger.info(f"  Shape: {adata.shape}")
    logger.info(f"  6-type distribution:\n{adata.obs['cell_type'].value_counts().to_string()}")
    logger.info(f"  Samples in reference: {sorted(adata.obs['sample'].unique().tolist())}")
    run(adata, out6, out7, args.link_counts_from)


def run(adata: sc.AnnData, out6: Path, out7: Path, link_counts_from: Optional[Path] = None) -> None:
    """Write 6-type rctd sample info and all 7-type exports for an already-labelled adata."""
    # Write sample info for 6-type rctd (CARD uses this for sampleInfo)
    logger.info(f"\nWriting sample info for 6-type rctd...")
    (out6 / "rctd").mkdir(parents=True, exist_ok=True)
    write_cell_types_with_sample(adata, out6 / "rctd")

    # Derive 7-type labels by splitting T cells
    labels_7 = split_t_cells(adata)
    validate_7type_reference(adata.obs["cell_type"], labels_7)

    # Build 7-type adata (same cells, same counts, updated labels only)
    adata_7 = adata.copy()
    adata_7.obs["cell_type"] = labels_7
    logger.info(f"\n7-type distribution:\n{adata_7.obs['cell_type'].value_counts().to_string()}")

    # Export all method formats
    out7.mkdir(parents=True, exist_ok=True)
    logger.info(f"\nExporting 7-type reference formats to {out7}/...")
    export_for_cell2location(adata_7, out7)
    export_for_tangram(adata_7, out7)
    if link_counts_from is None:
        export_for_rctd(adata_7, out7)
        export_for_seurat(adata_7, out7)
    else:
        for tool in ("rctd", "seurat"):
            export_labels_only(adata_7, out7, tool, link_counts_from)

    # Write sample info for 7-type rctd
    logger.info(f"\nWriting sample info for 7-type rctd...")
    write_cell_types_with_sample(adata_7, out7 / "rctd")

    logger.info("\nDone. Files written:")
    logger.info(f"  {out6}/rctd/cell_types_with_sample.csv")
    logger.info(f"  {out7}/rctd/cell_types_with_sample.csv (+ all other formats)")


if __name__ == "__main__":
    main()
