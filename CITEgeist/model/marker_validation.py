"""Single-cell marker gene validation for PC-MIL assignments.

For each nucleus assigned a cell type, checks whether the deconvolved GEX
from its parent spot supports that assignment by comparing marker gene
expression for the assigned type vs other types.
"""
import logging
from typing import Dict, List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_marker_scores(
    assignments: pd.DataFrame,
    gex_df: pd.DataFrame,
    rna_markers: Dict[str, List[str]],
) -> pd.DataFrame:
    """Compute per-nucleus marker gene scores.

    Args:
        assignments: DataFrame with columns [nucleus_id, barcode, cell_type].
        gex_df: Deconvolved GEX DataFrame indexed as 'barcode:::cell_type',
                columns are gene names.
        rna_markers: Dict mapping cell type -> list of expected RNA marker genes.

    Returns:
        DataFrame with columns [nucleus_id, barcode, cell_type,
        assigned_marker_mean, other_marker_mean, marker_score,
        markers_above_others].
    """
    records = []
    all_types = list(rna_markers.keys())
    all_marker_genes = set()
    for genes in rna_markers.values():
        all_marker_genes.update(genes)
    available_genes = set(gex_df.columns) & all_marker_genes

    for _, row in assignments.iterrows():
        nid = row["nucleus_id"]
        barcode = row["barcode"]
        assigned_type = row["cell_type"]

        gex_key = f"{barcode}:::{assigned_type}"
        if gex_key not in gex_df.index:
            logger.debug(f"GEX key {gex_key} not found, skipping nucleus {nid}")
            continue

        gex_row = gex_df.loc[gex_key]

        assigned_markers = [g for g in rna_markers.get(assigned_type, []) if g in available_genes]
        if not assigned_markers:
            continue
        assigned_mean = gex_row[assigned_markers].mean()

        other_means = []
        for other_type in all_types:
            if other_type == assigned_type:
                continue
            other_markers = [g for g in rna_markers.get(other_type, []) if g in available_genes]
            if other_markers:
                other_means.append(gex_row[other_markers].mean())

        other_mean = np.mean(other_means) if other_means else 0.0
        marker_score = float(np.log2((assigned_mean + 1e-6) / (other_mean + 1e-6)))

        records.append({
            "nucleus_id": nid,
            "barcode": barcode,
            "cell_type": assigned_type,
            "assigned_marker_mean": float(assigned_mean),
            "other_marker_mean": float(other_mean),
            "marker_score": marker_score,
            "markers_above_others": assigned_mean > other_mean,
        })

    return pd.DataFrame(records)


def summarize_validation(scores_df: pd.DataFrame) -> dict:
    """Summarize marker validation results per type and overall.

    Args:
        scores_df: Output from compute_marker_scores().

    Returns:
        Dict with 'per_type' and 'overall' summary metrics.
    """
    per_type = {}
    if scores_df.empty or "cell_type" not in scores_df.columns:
        return {"per_type": {}, "overall": {"n_nuclei": 0, "fraction_correct": 0.0, "median_marker_score": 0.0}}
    for cell_type, group in scores_df.groupby("cell_type"):
        per_type[cell_type] = {
            "n_nuclei": len(group),
            "median_marker_score": float(group["marker_score"].median()),
            "mean_marker_score": float(group["marker_score"].mean()),
            "fraction_correct": float(group["markers_above_others"].mean()),
        }

    total_nuclei = len(scores_df)
    if total_nuclei > 0:
        overall_fraction = float(scores_df["markers_above_others"].mean())
        overall_median = float(scores_df["marker_score"].median())
    else:
        overall_fraction = 0.0
        overall_median = 0.0

    return {
        "per_type": per_type,
        "overall": {
            "n_nuclei": total_nuclei,
            "fraction_correct": overall_fraction,
            "median_marker_score": overall_median,
        },
    }
