"""
Canonical marker enrichment, cross-patient consistency, and internal coherence.

Validates cell typing by checking that assigned cell types show expected
enrichment of canonical RNA markers. Works in both benchmark and
self-consistency modes.
"""

from __future__ import annotations

import logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

from . import QCResult
from .canonical_markers import CANONICAL_MARKERS, get_available_markers

logger = logging.getLogger(__name__)

__all__ = [
    "run_marker_enrichment",
    "compute_marker_enrichment",
    "check_cross_patient_consistency",
    "check_internal_coherence",
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


def compute_marker_enrichment(
    adata: AnnData,
    patient_id: str | None = None,
) -> pd.DataFrame:
    """Compute marker enrichment per cell type: log2FC, Wilcoxon p, AUC.

    Args:
        adata: AnnData with raw counts in X, obs['cell_type'].
        patient_id: Optional label for the patient column.

    Returns:
        DataFrame with columns: cell_type, marker, log2fc, pvalue, qvalue,
        auc, n_type, n_other, patient_id.
    """
    gene_names = list(adata.var_names)
    X = adata.X
    if sp.issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = np.asarray(X)

    rows = []
    for ct, markers in CANONICAL_MARKERS.items():
        available = get_available_markers(ct, gene_names)
        if not available:
            continue

        type_mask = adata.obs["cell_type"].values == ct
        other_mask = ~type_mask
        n_type = type_mask.sum()
        n_other = other_mask.sum()

        if n_type < 3 or n_other < 3:
            continue

        for marker in available:
            idx = gene_names.index(marker)
            vals_type = X_dense[type_mask, idx]
            vals_other = X_dense[other_mask, idx]

            mean_type = vals_type.mean() + 1
            mean_other = vals_other.mean() + 1
            log2fc = np.log2(mean_type / mean_other)

            try:
                stat, pval = mannwhitneyu(vals_type, vals_other, alternative="greater")
            except ValueError:
                pval = 1.0

            try:
                labels = np.concatenate([np.ones(n_type), np.zeros(n_other)])
                scores = np.concatenate([vals_type, vals_other])
                auc = roc_auc_score(labels, scores)
            except ValueError:
                auc = 0.5

            rows.append({
                "cell_type": ct,
                "marker": marker,
                "log2fc": log2fc,
                "pvalue": pval,
                "auc": auc,
                "n_type": n_type,
                "n_other": n_other,
                "patient_id": patient_id or "all",
            })

    df = pd.DataFrame(rows)
    if len(df) == 0:
        df = pd.DataFrame(columns=[
            "cell_type", "marker", "log2fc", "pvalue", "qvalue",
            "auc", "n_type", "n_other", "patient_id",
        ])
        return df

    # Benjamini-Hochberg FDR
    from scipy.stats import false_discovery_control
    try:
        df["qvalue"] = false_discovery_control(df["pvalue"].values, method="bh")
    except (ValueError, AttributeError):
        # Fallback for older scipy
        n_tests = len(df)
        ranked = df["pvalue"].rank()
        df["qvalue"] = df["pvalue"] * n_tests / ranked
        df["qvalue"] = df["qvalue"].clip(upper=1.0)
        df = df.sort_values("pvalue")
        df["qvalue"] = df["qvalue"].cummin()
        df = df.sort_index()

    return df


def check_cross_patient_consistency(
    enrichment_all_patients: pd.DataFrame,
    min_positive_fraction: float = 9 / 12,
    marker_fail_fraction: float = 0.3,
) -> list[str]:
    """Flag inconsistent markers and anomalous patients.

    Args:
        enrichment_all_patients: Concatenated enrichment tables with patient_id.
        min_positive_fraction: Minimum fraction of patients with positive log2FC.
        marker_fail_fraction: Flag patients where this fraction of markers fail.

    Returns:
        List of flag strings.
    """
    flags = []

    for (ct, marker), group in enrichment_all_patients.groupby(["cell_type", "marker"]):
        n_positive = (group["log2fc"] > 0).sum()
        n_total = len(group)
        if n_total == 0:
            continue
        frac = n_positive / n_total
        if frac < min_positive_fraction:
            flags.append(
                f"{ct}/{marker}: positive in {n_positive}/{n_total} patients "
                f"({frac:.0%}) — inconsistent enrichment"
            )

    for patient_id, group in enrichment_all_patients.groupby("patient_id"):
        n_markers = len(group)
        if n_markers == 0:
            continue
        n_failed = ((group["log2fc"] <= 0) | (group["qvalue"] >= 0.05)).sum()
        frac_failed = n_failed / n_markers
        if frac_failed > marker_fail_fraction:
            flags.append(
                f"Patient {patient_id}: {n_failed}/{n_markers} ({frac_failed:.0%}) "
                f"markers failed — anomalous sample"
            )

    return flags


def check_internal_coherence(
    proportions: pd.DataFrame,
    gex_layers: dict[str, pd.DataFrame],
    dominant_threshold: float = 0.3,
) -> dict:
    """Check that dominant-type spots have matching GEX marker expression.

    Args:
        proportions: Spot-level proportions (spots × types).
        gex_layers: type → DataFrame (spots × genes) from SACE.
        dominant_threshold: Minimum proportion to call a spot "dominant".

    Returns:
        Dict mapping type → {concordance_rate, n_dominant_spots, n_concordant}.
    """
    results = {}

    for ct in proportions.columns:
        markers = get_available_markers(
            ct, list(gex_layers.get(ct, pd.DataFrame()).columns)
        )
        if not markers or ct not in gex_layers:
            continue

        dominant = proportions[ct] > dominant_threshold
        dominant_spots = proportions.index[dominant]

        if len(dominant_spots) < 3:
            results[ct] = {"concordance_rate": np.nan, "n_dominant_spots": len(dominant_spots)}
            continue

        n_concordant = 0
        n_checked = 0

        for spot in dominant_spots:
            if spot not in gex_layers[ct].index:
                continue

            this_expr = 0.0
            other_expr = 0.0
            n_other = 0

            for marker in markers:
                if marker in gex_layers[ct].columns:
                    this_expr += gex_layers[ct].loc[spot, marker]

            for other_ct, layer in gex_layers.items():
                if other_ct == ct or spot not in layer.index:
                    continue
                for marker in markers:
                    if marker in layer.columns:
                        other_expr += layer.loc[spot, marker]
                        n_other += 1

            if n_other > 0:
                other_expr /= n_other
                n_checked += 1
                if this_expr > other_expr:
                    n_concordant += 1

        concordance = n_concordant / n_checked if n_checked > 0 else np.nan
        results[ct] = {
            "concordance_rate": concordance,
            "n_dominant_spots": len(dominant_spots),
            "n_concordant": n_concordant,
            "n_checked": n_checked,
        }

    return results


def _plot_enrichment_heatmap(
    enrichment_df: pd.DataFrame,
    cross_patient_df: pd.DataFrame | None = None,
) -> plt.Figure:
    """Heatmap of log2FC by cell type × marker, annotated with consistency."""
    avg = enrichment_df.groupby(["cell_type", "marker"])["log2fc"].mean().reset_index()
    pivot = avg.pivot(index="marker", columns="cell_type", values="log2fc")

    consistency = None
    if cross_patient_df is not None and "patient_id" in cross_patient_df.columns:
        cons = cross_patient_df.groupby(["cell_type", "marker"]).apply(
            lambda g: f"{(g['log2fc'] > 0).sum()}/{len(g)}"
        ).reset_index(name="consistency")
        consistency = cons.pivot(index="marker", columns="cell_type", values="consistency")

    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(14, max(8, len(pivot) * 0.5)))

        im = ax.imshow(pivot.values, cmap="RdBu_r", aspect="auto",
                       vmin=-3, vmax=3)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)

        if consistency is not None:
            for i, marker in enumerate(pivot.index):
                for j, ct in enumerate(pivot.columns):
                    if marker in consistency.index and ct in consistency.columns:
                        val = consistency.loc[marker, ct]
                        if pd.notna(val):
                            ax.text(j, i, str(val), ha="center", va="center",
                                    fontsize=MIN_FONT_SIZE - 2, color="black")

        fig.colorbar(im, ax=ax, shrink=0.6, label="log₂FC")
        ax.set_title("Canonical Marker Enrichment")
        fig.tight_layout()
    return fig


def _plot_cross_patient_strips(enrichment_df: pd.DataFrame) -> plt.Figure:
    """Strip plots of log2FC across patients per marker."""
    rng = np.random.default_rng(42)
    markers_by_type = enrichment_df.groupby("cell_type")["marker"].unique()
    all_pairs = []
    for ct, markers in markers_by_type.items():
        for m in markers:
            all_pairs.append(f"{ct}/{m}")

    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(max(14, len(all_pairs) * 0.6), 7))

        x_pos = 0
        x_ticks = []
        x_labels = []
        for ct, markers in markers_by_type.items():
            for marker in markers:
                sub = enrichment_df[
                    (enrichment_df["cell_type"] == ct) &
                    (enrichment_df["marker"] == marker)
                ]
                jitter = rng.random(len(sub)) * 0.3 - 0.15
                ax.scatter(
                    [x_pos] * len(sub) + jitter,
                    sub["log2fc"].values,
                    alpha=0.7, s=30, edgecolors="none",
                )
                ax.plot([x_pos - 0.3, x_pos + 0.3],
                        [sub["log2fc"].median()] * 2,
                        color="k", linewidth=2)
                x_ticks.append(x_pos)
                x_labels.append(f"{ct}\n{marker}")
                x_pos += 1

        ax.axhline(0, color="r", linestyle="--", alpha=0.5)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels, rotation=90, ha="center",
                           fontsize=max(8, MIN_FONT_SIZE - 4))
        ax.set_ylabel("log₂FC")
        ax.set_title("Cross-Patient Marker Consistency")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
    return fig


def _plot_coherence_bar(coherence: dict) -> plt.Figure:
    """Bar chart of internal coherence rate per cell type."""
    types = sorted(coherence.keys())
    rates = [coherence[ct].get("concordance_rate", 0) for ct in types]

    with plt.rc_context(_STYLE_PARAMS):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.bar(range(len(types)), rates, color="#4878CF")
        ax.set_xticks(range(len(types)))
        ax.set_xticklabels(types, rotation=45, ha="right")
        ax.set_ylabel("Concordance Rate")
        ax.set_ylim(0, 1.05)
        ax.set_title("Internal Coherence: Proportion ↔ GEX")
        ax.axhline(0.5, color="r", linestyle="--", alpha=0.5, label="Chance")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
    return fig


def run_marker_enrichment(
    adata: AnnData,
    proportions: pd.DataFrame,
    gex_layers: dict[str, pd.DataFrame] | None = None,
    multi_patient_enrichments: pd.DataFrame | None = None,
    patient_id: str | None = None,
) -> QCResult:
    """Run canonical marker validation, cross-patient, and internal coherence.

    Args:
        adata: AnnData with raw counts, obs['cell_type'].
        proportions: Spot-level proportions.
        gex_layers: Optional per-type GEX layers for internal coherence.
        multi_patient_enrichments: Pre-computed enrichment from multiple patients
            (concatenated DataFrames with patient_id column) for consistency check.
        patient_id: Label for this patient.

    Returns:
        QCResult with enrichment metrics and figures.
    """
    logger.info("Running marker enrichment QC...")
    metrics = {}
    flags = []
    figures = {}

    enrichment = compute_marker_enrichment(adata, patient_id)
    metrics["enrichment"] = enrichment

    if len(enrichment) > 0:
        if multi_patient_enrichments is not None and len(multi_patient_enrichments) > 0:
            consistency_flags = check_cross_patient_consistency(multi_patient_enrichments)
            flags.extend(consistency_flags)
            metrics["cross_patient_enrichments"] = multi_patient_enrichments
            figures["cross_patient_strips"] = _plot_cross_patient_strips(
                multi_patient_enrichments
            )
            figures["enrichment_heatmap"] = _plot_enrichment_heatmap(
                multi_patient_enrichments, multi_patient_enrichments
            )
        else:
            figures["enrichment_heatmap"] = _plot_enrichment_heatmap(enrichment)

    if gex_layers is not None:
        coherence = check_internal_coherence(proportions, gex_layers)
        metrics["internal_coherence"] = coherence
        figures["coherence_bar"] = _plot_coherence_bar(coherence)

        for ct, vals in coherence.items():
            rate = vals.get("concordance_rate", np.nan)
            if not np.isnan(rate) and rate < 0.5:
                flags.append(f"{ct}: internal coherence {rate:.0%} — below chance")

    return QCResult(metrics=metrics, flags=flags, figures=figures)
