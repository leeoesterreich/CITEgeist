"""
GEX error characterization against ground truth.

Computes per-type GEX correlations (Pearson/Spearman), per-gene
magnitude vs spatial fidelity analysis, and NMF program recovery
diagnostic. Also supports self-consistency mode with reference profiles.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
# pylint: disable=wrong-import-position,wrong-import-order
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity

from ._types import QCResult

# pylint: enable=wrong-import-position,wrong-import-order

logger = logging.getLogger(__name__)

__all__ = [
    "run_gex_qc",
    "compute_gex_correlations",
    "analyze_per_gene",
    "compare_nmf_programs",
]

MIN_FONT_SIZE = 14
_STYLE_PARAMS = {
    "font.size": MIN_FONT_SIZE,
    "axes.titlesize": MIN_FONT_SIZE + 2,
    "axes.labelsize": MIN_FONT_SIZE,
    "xtick.labelsize": MIN_FONT_SIZE,
    "ytick.labelsize": MIN_FONT_SIZE,
    "legend.fontsize": MIN_FONT_SIZE,
}

# Thresholds for gene classification
SPATIAL_R_THRESHOLD = 0.3  # Below this = spatially wrong


def compute_gex_correlations(
    pred_layers: dict[str, pd.DataFrame],
    gt_layers: dict[str, pd.DataFrame],
) -> dict:
    """Compute per-type GEX Pearson r and Spearman rho on log1p values.

    Args:
        pred_layers: type → DataFrame (spots × genes), raw counts.
        gt_layers: type → DataFrame (spots × genes), raw counts.

    Returns:
        Dict mapping type → {pearson_r, spearman_rho, pearson_spearman_gap}.
    """
    results = {}
    for ct in pred_layers:
        if ct not in gt_layers:
            continue
        pred = pred_layers[ct]
        gt = gt_layers[ct]

        common_genes = pred.columns.intersection(gt.columns)
        common_spots = pred.index.intersection(gt.index)
        if len(common_genes) < 5 or len(common_spots) < 5:
            results[ct] = {
                "pearson_r": np.nan,
                "spearman_rho": np.nan,
                "pearson_spearman_gap": np.nan,
                "n_genes": len(common_genes),
            }
            continue

        p = np.log1p(pred.loc[common_spots, common_genes].values).ravel()
        g = np.log1p(gt.loc[common_spots, common_genes].values).ravel()

        # Mask zeros for more meaningful correlation
        nonzero = (p > 0) | (g > 0)
        if nonzero.sum() < 10:
            results[ct] = {
                "pearson_r": np.nan,
                "spearman_rho": np.nan,
                "pearson_spearman_gap": np.nan,
                "n_genes": len(common_genes),
            }
            continue

        pr, _ = pearsonr(p[nonzero], g[nonzero])
        sr, _ = spearmanr(p[nonzero], g[nonzero])
        results[ct] = {
            "pearson_r": pr,
            "spearman_rho": sr,
            "pearson_spearman_gap": abs(pr - sr),
            "n_genes": len(common_genes),
        }

    return results


def analyze_per_gene(
    pred_layer: pd.DataFrame,
    gt_layer: pd.DataFrame,
) -> pd.DataFrame:
    """Per-gene spatial r and magnitude ratio for one cell type.

    Args:
        pred_layer: Predicted GEX (spots × genes).
        gt_layer: Ground truth GEX (spots × genes).

    Returns:
        DataFrame indexed by gene with spatial_r, log_mean_ratio, classification.
    """
    common_genes = pred_layer.columns.intersection(gt_layer.columns)
    common_spots = pred_layer.index.intersection(gt_layer.index)
    pred = pred_layer.loc[common_spots, common_genes]
    gt = gt_layer.loc[common_spots, common_genes]

    rows = []
    for gene in common_genes:
        p = np.log1p(pred[gene].values)
        g = np.log1p(gt[gene].values)

        gt_mean = gt[gene].mean()
        if gt_mean < 1e-10:
            rows.append(
                {
                    "gene": gene,
                    "spatial_r": np.nan,
                    "log_mean_ratio": np.nan,
                    "classification": "zero_gt",
                    "gt_variance": 0.0,
                }
            )
            continue

        if np.std(p) < 1e-10 or np.std(g) < 1e-10:
            spatial_r = 0.0
        else:
            spatial_r, _ = pearsonr(p, g)

        pred_mean = pred[gene].mean()
        log_ratio = np.log2((pred_mean + 1) / (gt_mean + 1))

        if spatial_r >= SPATIAL_R_THRESHOLD:
            classification = "magnitude_only" if abs(log_ratio) > 0.5 else "good"
        else:
            classification = "spatial_mismatch"

        rows.append(
            {
                "gene": gene,
                "spatial_r": spatial_r,
                "log_mean_ratio": log_ratio,
                "classification": classification,
                "gt_variance": gt[gene].var(),
            }
        )

    return pd.DataFrame(rows).set_index("gene")


def compare_nmf_programs(
    pred_layer: pd.DataFrame,
    gt_layer: pd.DataFrame,
    k: int = 5,
) -> np.ndarray:
    """Compare NMF programs between predicted and GT GEX.

    Args:
        pred_layer: Predicted GEX (spots × genes), raw counts.
        gt_layer: Ground truth GEX (spots × genes), raw counts.
        k: Number of NMF programs.

    Returns:
        (k, k) cosine similarity matrix between GT and predicted programs.
    """
    common_genes = pred_layer.columns.intersection(gt_layer.columns)
    common_spots = pred_layer.index.intersection(gt_layer.index)

    pred = np.clip(pred_layer.loc[common_spots, common_genes].values, 0, None)
    gt = np.clip(gt_layer.loc[common_spots, common_genes].values, 0, None)

    k = min(k, pred.shape[0], pred.shape[1])

    nmf_gt = NMF(n_components=k, random_state=42, max_iter=200)
    nmf_gt.fit(gt)
    W_gt = nmf_gt.components_  # (k, genes)

    nmf_pred = NMF(n_components=k, random_state=42, max_iter=200)
    nmf_pred.fit(pred)
    W_pred = nmf_pred.components_  # (k, genes)

    return np.asarray(cosine_similarity(W_gt, W_pred))


def _plot_bar_pearson_spearman(correlations: dict) -> plt.Figure:
    """Bar chart: Pearson r vs Spearman rho per cell type."""
    cell_types = [ct for ct in correlations if ct != "overall"]

    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        x = np.arange(len(cell_types))
        width = 0.35

        pearson_vals = [correlations[ct].get("pearson_r", 0) for ct in cell_types]
        spearman_vals = [correlations[ct].get("spearman_rho", 0) for ct in cell_types]

        ax.bar(x - width / 2, pearson_vals, width, label="Pearson r", color="#4878CF")
        ax.bar(x + width / 2, spearman_vals, width, label="Spearman ρ", color="#6ACC65")
        ax.set_xticks(x)
        ax.set_xticklabels(cell_types, rotation=45, ha="right")
        ax.set_ylabel("Correlation")
        ax.set_title("GEX: Pearson r vs Spearman ρ per Cell Type")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
    return fig


def _plot_gene_scatter(gene_analysis: pd.DataFrame, cell_type: str) -> plt.Figure:
    """Gene-level scatter: spatial r vs mean ratio, colored by variance quartile."""
    df = gene_analysis.dropna(subset=["spatial_r", "log_mean_ratio"])
    if len(df) == 0:
        fig, ax = plt.subplots(1, 1, figsize=(8, 7))
        ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=MIN_FONT_SIZE)
        return fig

    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(9, 7))

        df = df.copy()
        df["var_quartile"] = pd.qcut(df["gt_variance"].clip(lower=1e-10), q=4, labels=["Q1", "Q2", "Q3", "Q4"])
        colors = {"Q1": "#c0c0c0", "Q2": "#7eb0d5", "Q3": "#4878CF", "Q4": "#1f3d7a"}

        for q in ["Q1", "Q2", "Q3", "Q4"]:
            sub = df[df["var_quartile"] == q]
            ax.scatter(
                sub["spatial_r"],
                sub["log_mean_ratio"],
                label=f"Variance {q}",
                c=colors[q],
                alpha=0.6,
                s=20,
                edgecolors="none",
            )

        ax.axhline(0, color="k", linestyle="--", alpha=0.3)
        ax.axvline(SPATIAL_R_THRESHOLD, color="r", linestyle="--", alpha=0.3)
        ax.set_xlabel("Spatial Pearson r")
        ax.set_ylabel("log₂(Mean Ratio)")
        ax.set_title(f"{cell_type}: Per-Gene Spatial Fidelity vs Magnitude")
        ax.legend(loc="upper left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
    return fig


def _plot_nmf_cosine(cosine_matrix: np.ndarray) -> plt.Figure:
    """Heatmap of cosine similarity between GT and predicted NMF programs."""
    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        k = cosine_matrix.shape[0]
        im = ax.imshow(cosine_matrix, cmap="YlOrRd", vmin=0, vmax=1, aspect="equal")
        ax.set_xticks(range(k))
        ax.set_yticks(range(k))
        ax.set_xticklabels([f"Pred {i+1}" for i in range(k)])
        ax.set_yticklabels([f"GT {i+1}" for i in range(k)])
        ax.set_title("NMF Program Recovery\n(Cosine Similarity)")
        fig.colorbar(im, ax=ax, shrink=0.7)
        fig.tight_layout()
    return fig


def run_gex_qc(
    pred_gex_layers: dict[str, pd.DataFrame],
    gt_gex_layers: dict[str, pd.DataFrame] | None = None,
    reference_profiles: dict[str, pd.DataFrame] | None = None,
    nmf_k: int = 5,
    exemplar_type: str | None = None,
) -> QCResult:
    """Run GEX error characterization.

    Args:
        pred_gex_layers: type → DataFrame (spots × genes), predicted.
        gt_gex_layers: type → DataFrame (spots × genes), ground truth.
            If None, runs self-consistency mode.
        reference_profiles: type → DataFrame for self-consistency cosine
            similarity. Ignored in benchmark mode.
        nmf_k: Number of NMF programs for propagation diagnostic.
        exemplar_type: Cell type for gene-level scatter. If None, picks
            the type with the most genes in common with GT.

    Returns:
        QCResult with GEX metrics and figures.
    """
    logger.info("Running GEX QC...")
    metrics = {}
    flags = []
    figures = {}

    if gt_gex_layers is not None:
        # === Benchmark mode ===
        correlations = compute_gex_correlations(pred_gex_layers, gt_gex_layers)
        metrics["correlations"] = correlations

        for ct, vals in correlations.items():
            gap = vals.get("pearson_spearman_gap", 0)
            if gap and gap > 0.15:
                flags.append(f"GEX {ct}: Spearman-Pearson gap = {gap:.3f}")

        if exemplar_type is None:
            common_counts = {}
            for ct in pred_gex_layers:
                if ct in gt_gex_layers:
                    n = len(pred_gex_layers[ct].columns.intersection(gt_gex_layers[ct].columns))
                    common_counts[ct] = n
            exemplar_type = max(common_counts, key=common_counts.get) if common_counts else None

        if exemplar_type and exemplar_type in pred_gex_layers and exemplar_type in gt_gex_layers:
            gene_analysis = analyze_per_gene(pred_gex_layers[exemplar_type], gt_gex_layers[exemplar_type])
            metrics["gene_analysis"] = {exemplar_type: gene_analysis}

            if len(gene_analysis) > 0:
                vc = gene_analysis["classification"].value_counts()
                total = len(gene_analysis)
                metrics["gene_classification_summary"] = {cls: int(count) for cls, count in vc.items()}
                pct_spatial_ok = (vc.get("good", 0) + vc.get("magnitude_only", 0)) / total
                metrics["pct_spatially_correct"] = pct_spatial_ok

            figures["gene_scatter"] = _plot_gene_scatter(gene_analysis, exemplar_type)

        if exemplar_type and exemplar_type in pred_gex_layers and exemplar_type in gt_gex_layers:
            cosine_mat = compare_nmf_programs(
                pred_gex_layers[exemplar_type],
                gt_gex_layers[exemplar_type],
                k=nmf_k,
            )
            metrics["nmf_cosine"] = cosine_mat  # type: ignore[assignment]
            figures["nmf_recovery"] = _plot_nmf_cosine(cosine_mat)

        figures["bar_pearson_spearman"] = _plot_bar_pearson_spearman(correlations)

    elif reference_profiles is not None:
        # === Self-consistency mode ===
        cosine_results = {}
        for ct in pred_gex_layers:
            if ct not in reference_profiles:
                continue
            pred = pred_gex_layers[ct]
            ref = reference_profiles[ct]
            common_genes = pred.columns.intersection(ref.columns)
            if len(common_genes) < 5:
                cosine_results[ct] = np.nan
                continue

            pred_mean = pred[common_genes].mean(axis=0).values.reshape(1, -1)
            ref_mean = ref[common_genes].mean(axis=0).values.reshape(1, -1)
            cos = cosine_similarity(pred_mean, ref_mean)[0, 0]
            cosine_results[ct] = cos

        metrics["reference_cosine"] = cosine_results
    else:
        logger.warning("GEX QC: no GT or reference provided, skipping")

    return QCResult(metrics=metrics, flags=flags, figures=figures)
