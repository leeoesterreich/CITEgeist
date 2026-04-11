"""
Standard single-cell QC metrics and empty cell detection.

Computes UMI counts, genes detected, mitochondrial %, ribosomal %
per cell. Detects empty cells from protein-RNA sparsity mismatch
using configurable thresholds.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
# pylint: disable=wrong-import-position,wrong-import-order
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from ._types import QCResult

# pylint: enable=wrong-import-position,wrong-import-order

logger = logging.getLogger(__name__)

__all__ = ["run_single_cell_qc", "compute_cell_metrics", "detect_empty_cells", "check_compartment_emptiness"]

# ==================== Styling ====================

MIN_FONT_SIZE = 14

_STYLE_PARAMS = {
    "font.size": MIN_FONT_SIZE,
    "axes.titlesize": MIN_FONT_SIZE + 2,
    "axes.labelsize": MIN_FONT_SIZE,
    "xtick.labelsize": MIN_FONT_SIZE,
    "ytick.labelsize": MIN_FONT_SIZE,
    "legend.fontsize": MIN_FONT_SIZE,
}


def compute_cell_metrics(adata: AnnData) -> pd.DataFrame:
    """Compute per-cell QC metrics from raw count matrix.

    Args:
        adata: AnnData with raw counts in X, gene names in var_names.

    Returns:
        DataFrame indexed by cell with columns:
        total_umi, n_genes, pct_mt, pct_ribo.
    """
    X = adata.X
    total_umi = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()

    gene_names = adata.var_names
    mt_mask = gene_names.str.startswith("MT-")
    ribo_mask = gene_names.str.startswith("RPS") | gene_names.str.startswith("RPL")

    if sp.issparse(X):
        mt_sum = np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
        ribo_sum = np.asarray(X[:, ribo_mask].sum(axis=1)).ravel()
    else:
        mt_sum = X[:, mt_mask].sum(axis=1)
        ribo_sum = X[:, ribo_mask].sum(axis=1)

    # Avoid division by zero
    safe_total = np.where(total_umi > 0, total_umi, 1)
    pct_mt = (mt_sum / safe_total) * 100
    pct_ribo = (ribo_sum / safe_total) * 100

    return pd.DataFrame(
        {
            "total_umi": total_umi,
            "n_genes": n_genes,
            "pct_mt": pct_mt,
            "pct_ribo": pct_ribo,
        },
        index=adata.obs_names,
    )


def detect_empty_cells(
    metrics: pd.DataFrame,
    umi_threshold: int = 50,
    genes_threshold: int = 25,
) -> np.ndarray:
    """Flag cells below UMI or gene detection thresholds.

    Args:
        metrics: Output of compute_cell_metrics.
        umi_threshold: Minimum UMI count.
        genes_threshold: Minimum genes detected.

    Returns:
        Boolean array — True for cells flagged as empty.
    """
    is_empty = (metrics["total_umi"] < umi_threshold) | (metrics["n_genes"] < genes_threshold)
    return np.asarray(is_empty.values)


def check_compartment_emptiness(
    cell_types: pd.Series,
    is_empty: np.ndarray,
    threshold: float = 0.5,
) -> list[str]:
    """Flag cell types where >threshold fraction of cells are empty.

    Args:
        cell_types: Series of cell type labels per cell.
        is_empty: Boolean array from detect_empty_cells.
        threshold: Fraction above which to flag (default 0.5).

    Returns:
        List of warning flag strings.
    """
    flags = []
    for ct in cell_types.unique():
        mask = cell_types.values == ct
        n_total = mask.sum()
        if n_total == 0:
            continue
        n_empty = is_empty[mask].sum()
        frac = n_empty / n_total
        if frac > threshold:
            flags.append(f"{ct}: {n_empty}/{n_total} ({frac:.0%}) cells below QC thresholds")
    return flags


def _build_summary_table(
    metrics: pd.DataFrame,
    cell_types: pd.Series,
    is_empty: np.ndarray,
) -> pd.DataFrame:
    """Build per-cell-type summary table.

    Returns:
        DataFrame with columns: n_cells, median_UMI, median_genes,
        median_MT%, n_empty, pct_empty.
    """
    metrics = metrics.copy()
    metrics["cell_type"] = cell_types.values
    metrics["is_empty"] = is_empty

    rows = []
    for ct in sorted(cell_types.unique()):
        sub = metrics[metrics["cell_type"] == ct]
        rows.append(
            {
                "cell_type": ct,
                "n_cells": len(sub),
                "median_UMI": sub["total_umi"].median(),
                "median_genes": sub["n_genes"].median(),
                "median_MT%": sub["pct_mt"].median(),
                "n_empty": sub["is_empty"].sum(),
                "pct_empty": sub["is_empty"].mean() * 100,
            }
        )
    return pd.DataFrame(rows).set_index("cell_type")


def _plot_violin_qc(
    metrics: pd.DataFrame,
    cell_types: pd.Series,
) -> plt.Figure:
    """Violin plots of UMI, genes detected, MT% per cell type.

    Returns:
        matplotlib Figure with 3 subpanels.
    """
    with plt.rc_context(_STYLE_PARAMS):
        fig, axes = plt.subplots(1, 3, figsize=(24, 7))

        types_sorted = sorted(cell_types.unique())
        data_by_type = {ct: metrics[cell_types.values == ct] for ct in types_sorted}

        for ax, col, label in zip(
            axes,
            ["total_umi", "n_genes", "pct_mt"],
            ["Total UMI", "Genes Detected", "Mitochondrial %"],
        ):
            violin_data = [data_by_type[ct][col].values for ct in types_sorted]
            parts = ax.violinplot(violin_data, showmedians=True)
            for pc in parts["bodies"]:
                pc.set_alpha(0.7)
            ax.set_xticks(range(1, len(types_sorted) + 1))
            ax.set_xticklabels(types_sorted, rotation=45, ha="right")
            ax.set_ylabel(label)
            ax.set_title(label)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        fig.tight_layout()
    return fig


def _plot_empty_heatmap(
    summary_table: pd.DataFrame,
    patient_id: str | None = None,
) -> plt.Figure:
    """Heatmap of empty cell fraction per cell type.

    Args:
        summary_table: Output of _build_summary_table.
        patient_id: Optional label for the heatmap title.

    Returns:
        matplotlib Figure.
    """
    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

        pct = summary_table["pct_empty"].values.reshape(1, -1)
        im = ax.imshow(pct, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)
        ax.set_xticks(range(len(summary_table)))
        ax.set_xticklabels(summary_table.index, rotation=45, ha="right")
        ax.set_yticks([0])
        ax.set_yticklabels([patient_id or "All"])
        cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label("% Empty Cells")
        ax.set_title("Empty Cell Fraction by Cell Type")
        fig.tight_layout()
    return fig


def run_single_cell_qc(
    adata: AnnData,
    umi_threshold: int = 50,
    genes_threshold: int = 25,
    patient_id: str | None = None,
) -> QCResult:
    """Run single-cell QC: metrics, empty detection, summary, figures.

    Args:
        adata: AnnData with raw counts in X, obs['cell_type'], obs['spot_id'].
        umi_threshold: Min UMI count for empty cell detection.
        genes_threshold: Min genes detected for empty cell detection.
        patient_id: Optional label for plots.

    Returns:
        QCResult with metrics (summary_table, cell_metrics, is_empty),
        flags, and figures (violin_qc, empty_heatmap).
    """
    logger.info("Running single-cell QC...")
    cell_types = adata.obs["cell_type"]

    # Compute metrics
    cell_metrics = compute_cell_metrics(adata)
    is_empty = detect_empty_cells(cell_metrics, umi_threshold, genes_threshold)
    compartment_flags = check_compartment_emptiness(cell_types, is_empty)

    n_empty = is_empty.sum()
    logger.info(
        "Flagged %s/%s cells as empty (UMI<%s or genes<%s)", n_empty, len(adata), umi_threshold, genes_threshold
    )

    # Summary table
    summary = _build_summary_table(cell_metrics, cell_types, is_empty)

    # Flags
    flags = []
    if n_empty > 0:
        flags.append(
            f"{n_empty} cells with <{umi_threshold} UMIs or " f"<{genes_threshold} genes detected — flagged as empty"
        )
    flags.extend(compartment_flags)

    # Figures
    fig_violin = _plot_violin_qc(cell_metrics, cell_types)
    fig_heatmap = _plot_empty_heatmap(summary, patient_id)

    return QCResult(
        metrics={
            "summary_table": summary,
            "cell_metrics": cell_metrics,
            "is_empty": is_empty,
        },
        flags=flags,
        figures={
            "violin_qc": fig_violin,
            "empty_heatmap": fig_heatmap,
        },
    )
