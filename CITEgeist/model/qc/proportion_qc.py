"""
Proportion error characterization against ground truth.

Computes Pearson/Spearman correlations, error decomposition (scaling
vs spatial misallocation), discrete assignment round-trip, and
confusion analysis. Benchmark mode only.
"""

from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
# pylint: disable=wrong-import-position,wrong-import-order
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr, spearmanr

from ._types import QCResult

# pylint: enable=wrong-import-position,wrong-import-order

logger = logging.getLogger(__name__)

__all__ = [
    "run_proportion_qc",
    "compute_proportion_correlations",
    "decompose_scaling_error",
    "compute_discrete_round_trip",
    "compute_error_confusion",
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


def compute_proportion_correlations(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
) -> dict:
    """Compute Pearson r, Spearman rho, RMSE, MAE per cell type and overall.

    Args:
        pred: Predicted proportions (spots × types).
        gt: Ground truth proportions (spots × types), aligned index.

    Returns:
        Dict mapping cell_type → {pearson_r, spearman_rho, rmse, mae, pearson_spearman_gap}
        plus "overall" key with flattened metrics.
    """
    common_spots = pred.index.intersection(gt.index)
    common_types = [c for c in pred.columns if c in gt.columns]
    pred = pred.loc[common_spots, common_types]
    gt = gt.loc[common_spots, common_types]

    results = {}
    for ct in common_types:
        p, g = pred[ct].values, gt[ct].values
        if len(p) < 3 or np.std(g) < 1e-10:
            results[ct] = {
                "pearson_r": np.nan,
                "spearman_rho": np.nan,
                "rmse": np.nan,
                "mae": np.nan,
                "pearson_spearman_gap": np.nan,
            }
            continue
        pr, _ = pearsonr(p, g)
        sr, _ = spearmanr(p, g)
        rmse = np.sqrt(np.mean((p - g) ** 2))
        mae = np.mean(np.abs(p - g))
        results[ct] = {
            "pearson_r": pr,
            "spearman_rho": sr,
            "rmse": rmse,
            "mae": float(mae),
            "pearson_spearman_gap": abs(pr - sr),
        }

    # Overall (flattened)
    p_flat = pred.values.ravel()
    g_flat = gt.values.ravel()
    pr_all, _ = pearsonr(p_flat, g_flat)
    sr_all, _ = spearmanr(p_flat, g_flat)
    results["overall"] = {
        "pearson_r": pr_all,
        "spearman_rho": sr_all,
        "rmse": np.sqrt(np.mean((p_flat - g_flat) ** 2)),
        "mae": np.mean(np.abs(p_flat - g_flat)),
        "pearson_spearman_gap": abs(pr_all - sr_all),
    }

    # Per-spot JSD
    jsd_vals = []
    for i in range(len(pred)):
        p_row = pred.iloc[i].values
        g_row = gt.iloc[i].values
        p_row = np.clip(p_row, 1e-10, None)
        g_row = np.clip(g_row, 1e-10, None)
        p_row = p_row / p_row.sum()
        g_row = g_row / g_row.sum()
        jsd_vals.append(jensenshannon(p_row, g_row) ** 2)
    results["jsd_per_spot"] = np.array(jsd_vals)  # type: ignore[assignment]
    results["mean_jsd"] = np.mean(jsd_vals)

    return results


def decompose_scaling_error(pred: np.ndarray, gt: np.ndarray) -> dict:
    """Decompose prediction error into scaling and residual components.

    Fits linear regression pred ~ gt; scaling error = (fitted - gt),
    residual error = (pred - fitted).

    Args:
        pred: Predicted values.
        gt: Ground truth values.

    Returns:
        Dict with slope, scaling_fraction, residual_fraction.
    """
    if len(gt) < 3 or np.std(gt) < 1e-10:
        return {"slope": np.nan, "scaling_fraction": np.nan, "residual_fraction": np.nan}

    A = np.column_stack([gt, np.ones_like(gt)])
    params, _, _, _ = np.linalg.lstsq(A, pred, rcond=None)
    fitted = A @ params

    total_ss = np.sum((pred - gt) ** 2)
    if total_ss < 1e-10:
        return {"slope": params[0], "scaling_fraction": 0.0, "residual_fraction": 0.0}

    scaling_ss = np.sum((fitted - gt) ** 2)
    residual_ss = np.sum((pred - fitted) ** 2)

    return {
        "slope": params[0],
        "scaling_fraction": scaling_ss / total_ss,
        "residual_fraction": residual_ss / total_ss,
    }


def compute_discrete_round_trip(
    adata: AnnData,
    proportions: pd.DataFrame,
) -> dict:
    """Aggregate single-cell assignments to per-spot fractions.

    Args:
        adata: AnnData with obs['cell_type'] and obs['spot_id'].
        proportions: Continuous proportions for column alignment.

    Returns:
        Dict with discrete_proportions DataFrame and n_skipped_spots count.
    """
    cell_types = sorted(proportions.columns)

    counts = adata.obs.groupby(["spot_id", "cell_type"]).size().unstack(fill_value=0)
    for ct in cell_types:
        if ct not in counts.columns:
            counts[ct] = 0
    counts = counts[cell_types]

    totals = counts.sum(axis=1)
    valid = totals > 0
    n_skipped = (~valid).sum()

    if n_skipped > 0:
        logger.info("Discrete round-trip: skipping %s spots with zero nuclei", n_skipped)

    discrete_props = counts.loc[valid].div(totals[valid], axis=0)

    return {
        "discrete_proportions": discrete_props,
        "n_skipped_spots": int(n_skipped),
    }


def compute_error_confusion(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
) -> pd.DataFrame:
    """Build signed error confusion matrix.

    For each spot, find which type is most over-predicted and which is
    most under-predicted. Accumulate into a types × types matrix.

    Args:
        pred: Predicted proportions.
        gt: Ground truth proportions.

    Returns:
        DataFrame (types × types) — rows = over-predicted, cols = under-predicted.
    """
    common_spots = pred.index.intersection(gt.index)
    common_types = [c for c in pred.columns if c in gt.columns]
    diff = pred.loc[common_spots, common_types] - gt.loc[common_spots, common_types]

    confusion = pd.DataFrame(0.0, index=common_types, columns=common_types)
    for _, row in diff.iterrows():
        over_type = row.idxmax()
        under_type = row.idxmin()
        if over_type != under_type:
            confusion.loc[over_type, under_type] += row[over_type]

    confusion /= len(common_spots)
    return confusion


def _plot_scatter_grid(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
    correlations: dict,
) -> plt.Figure:
    """Scatter grid: predicted vs GT per cell type."""
    cell_types = [c for c in pred.columns if c in gt.columns]
    n_types = len(cell_types)
    ncols = min(n_types, 4)
    nrows = int(np.ceil(n_types / ncols))

    with plt.rc_context(_STYLE_PARAMS):
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows))
        axes = np.atleast_2d(axes)

        common = pred.index.intersection(gt.index)
        for idx, ct in enumerate(cell_types):
            ax = axes[idx // ncols, idx % ncols]
            x, y = gt.loc[common, ct].values, pred.loc[common, ct].values
            ax.scatter(x, y, alpha=0.3, s=10, edgecolors="none")
            ax.plot([0, 1], [0, 1], "k--", alpha=0.5, linewidth=1)

            metrics = correlations.get(ct, {})
            pr = metrics.get("pearson_r", np.nan)
            sr = metrics.get("spearman_rho", np.nan)
            ax.set_title(f"{ct}\nPearson={pr:.3f}  Spearman={sr:.3f}")
            ax.set_xlabel("Ground Truth")
            ax.set_ylabel("Predicted")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        for idx in range(n_types, nrows * ncols):
            axes[idx // ncols, idx % ncols].set_visible(False)

        fig.tight_layout()
    return fig


def _plot_spatial_errors(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
    spatial_coords: np.ndarray,
    spot_ids: pd.Index,
    key_types: list | None = None,
) -> plt.Figure:
    """Spatial error maps for key cell types."""
    cell_types = [c for c in pred.columns if c in gt.columns]
    if key_types is None:
        rmses = {ct: np.sqrt(np.mean((pred[ct] - gt[ct]) ** 2)) for ct in cell_types}
        key_types = sorted(rmses, key=rmses.get, reverse=True)[:3]

    common = pred.index.intersection(gt.index).intersection(spot_ids)
    spot_to_idx = {s: i for i, s in enumerate(spot_ids)}
    coord_idx = [spot_to_idx[s] for s in common]
    coords = spatial_coords[coord_idx]

    with plt.rc_context(_STYLE_PARAMS):
        fig, axes = plt.subplots(1, len(key_types), figsize=(8 * len(key_types), 7))
        if len(key_types) == 1:
            axes = [axes]

        for ax, ct in zip(axes, key_types):
            residual = pred.loc[common, ct].values - gt.loc[common, ct].values
            vmax = max(abs(residual.min()), abs(residual.max()))
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=residual,
                cmap="RdBu_r",
                vmin=-vmax,
                vmax=vmax,
                s=15,
                edgecolors="none",
            )
            fig.colorbar(scatter, ax=ax, shrink=0.7, label="Pred − GT")
            ax.set_title(f"{ct} Residual")
            ax.set_aspect("equal")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        fig.tight_layout()
    return fig


def _plot_round_trip(
    continuous_corr: dict,
    discrete_corr: dict,
    cell_types: list,
) -> plt.Figure:
    """Paired bar chart: continuous vs discrete proportions against GT."""
    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        x = np.arange(len(cell_types))
        width = 0.35

        cont_r = [continuous_corr.get(ct, {}).get("pearson_r", 0) for ct in cell_types]
        disc_r = [discrete_corr.get(ct, {}).get("pearson_r", 0) for ct in cell_types]

        ax.bar(x - width / 2, cont_r, width, label="Continuous", color="#4878CF")
        ax.bar(x + width / 2, disc_r, width, label="Discrete", color="#D65F5F")
        ax.set_xticks(x)
        ax.set_xticklabels(cell_types, rotation=45, ha="right")
        ax.set_ylabel("Pearson r vs GT")
        ax.set_title("Continuous vs Discrete Proportions")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
    return fig


def run_proportion_qc(
    adata: AnnData,
    proportions: pd.DataFrame,
    gt_proportions: pd.DataFrame,
    spatial_coords: np.ndarray | None = None,
) -> QCResult:
    """Run proportion error characterization against ground truth.

    Args:
        adata: AnnData with obs['cell_type'], obs['spot_id'].
        proportions: Predicted proportions (spots × types).
        gt_proportions: Ground truth proportions (spots × types).
        spatial_coords: Optional (n_spots, 2) array for spatial plots.
            If None, attempts to extract from adata unique spots.

    Returns:
        QCResult with correlation metrics, error decomposition,
        round-trip analysis, and 3 figures.
    """
    logger.info("Running proportion QC...")
    cell_types = sorted([c for c in proportions.columns if c in gt_proportions.columns])

    common_spots = proportions.index.intersection(gt_proportions.index)
    pred = proportions.loc[common_spots]
    gt = gt_proportions.loc[common_spots]

    correlations = compute_proportion_correlations(pred, gt)

    decomposition = {}
    for ct in cell_types:
        decomposition[ct] = decompose_scaling_error(pred[ct].values, gt[ct].values)

    error_confusion = compute_error_confusion(pred, gt)

    round_trip = compute_discrete_round_trip(adata, proportions)
    discrete_props = round_trip["discrete_proportions"]
    common_discrete = discrete_props.index.intersection(gt_proportions.index)
    if len(common_discrete) > 3:
        discrete_corr = compute_proportion_correlations(
            discrete_props.loc[common_discrete, cell_types],
            gt_proportions.loc[common_discrete, cell_types],
        )
    else:
        discrete_corr = {}

    flags = []
    for ct in cell_types:
        gap = correlations.get(ct, {}).get("pearson_spearman_gap", 0)
        if gap and gap > 0.15:
            flags.append(f"{ct} Spearman-Pearson gap = {gap:.3f} — check outliers")
    if round_trip["n_skipped_spots"] > 0:
        flags.append(f"Discrete round-trip: {round_trip['n_skipped_spots']} spots " f"with zero nuclei skipped")

    fig_scatter = _plot_scatter_grid(pred, gt, correlations)

    if spatial_coords is None:
        if "spatial" in adata.obsm:
            first_idx = adata.obs.groupby("spot_id").cumcount() == 0
            spatial_coords = adata[first_idx].obsm["spatial"]
            spot_ids = adata.obs.loc[first_idx, "spot_id"]
        else:
            spatial_coords = np.zeros((len(common_spots), 2))
            spot_ids = pd.Index(common_spots)
    else:
        spot_ids = pd.Index(common_spots)

    fig_spatial = _plot_spatial_errors(pred, gt, spatial_coords, spot_ids)
    fig_roundtrip = _plot_round_trip(correlations, discrete_corr, cell_types)

    return QCResult(
        metrics={
            "correlations": correlations,
            "decomposition": decomposition,
            "error_confusion": error_confusion,
            "round_trip": round_trip,
            "discrete_correlations": discrete_corr,
        },
        flags=flags,
        figures={
            "scatter_grid": fig_scatter,
            "spatial_errors": fig_spatial,
            "round_trip": fig_roundtrip,
        },
    )
