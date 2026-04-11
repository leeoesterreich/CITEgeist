"""
QC orchestrator and figure composition.

Runs all QC modules in order, composes supplementary figures,
and writes summary artifacts.
"""

from __future__ import annotations

import json
import logging
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # pylint: disable=wrong-import-position
import numpy as np  # pylint: disable=wrong-import-position
import pandas as pd  # pylint: disable=wrong-import-position
from anndata import AnnData  # pylint: disable=wrong-import-position

from ._types import QCResult  # pylint: disable=wrong-import-position

logger = logging.getLogger(__name__)

__all__ = ["run_qc"]

MIN_FONT_SIZE = 14
_STYLE_PARAMS = {
    "font.size": MIN_FONT_SIZE,
    "axes.titlesize": MIN_FONT_SIZE + 2,
    "axes.labelsize": MIN_FONT_SIZE,
    "xtick.labelsize": MIN_FONT_SIZE,
    "ytick.labelsize": MIN_FONT_SIZE,
    "legend.fontsize": MIN_FONT_SIZE,
}


def _serialize_metrics(obj):
    """JSON-safe serialization of metrics dict."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, pd.DataFrame):
        return {"_type": "DataFrame", "shape": list(obj.shape), "columns": list(obj.columns)}
    if isinstance(obj, pd.Series):
        return {"_type": "Series", "values": obj.tolist()}
    if isinstance(obj, dict):
        return {k: _serialize_metrics(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_serialize_metrics(v) for v in obj]
    return obj


def run_qc(
    adata_per_cell: AnnData,
    proportions: pd.DataFrame,
    *,
    mode: str = "self_consistency",
    gt_proportions: pd.DataFrame | None = None,
    gt_gex_layers: dict[str, pd.DataFrame] | None = None,
    pred_gex_layers: dict[str, pd.DataFrame] | None = None,
    _reference_adata: AnnData | None = None,
    output_dir: str = "./qc_output",
    empty_umi_threshold: int = 50,
    empty_genes_threshold: int = 25,
) -> dict[str, QCResult]:
    """Orchestrate all QC checks.

    Args:
        adata_per_cell: AnnData with raw counts, obs['cell_type'], obs['spot_id'].
        proportions: Spot-level proportions (spots × types).
        mode: "benchmark" or "self_consistency".
        gt_proportions: Ground truth proportions (benchmark mode).
        gt_gex_layers: Ground truth GEX layers (benchmark mode).
        pred_gex_layers: Predicted GEX layers from SACE (benchmark mode).
            If None, GT layers are used as a scaffold (r=1 trivially).
        reference_adata: Reference scRNA-seq for self-consistency cosine sim.
        output_dir: Directory for output artifacts.
        empty_umi_threshold: Min UMI for empty cell detection.
        empty_genes_threshold: Min genes for empty cell detection.

    Returns:
        Dict mapping module_name → QCResult.

    Raises:
        ValueError: If mode="benchmark" and gt_proportions is None.
    """
    if mode == "benchmark" and gt_proportions is None:
        raise ValueError("Benchmark mode requires gt_proportions")

    os.makedirs(output_dir, exist_ok=True)
    results = {}
    all_flags = []

    # === Step 1: Single-cell QC ===
    from .single_cell_qc import run_single_cell_qc  # pylint: disable=import-outside-toplevel

    logger.info("=== Step 1: Single-cell QC ===")
    sc_result = run_single_cell_qc(
        adata_per_cell,
        umi_threshold=empty_umi_threshold,
        genes_threshold=empty_genes_threshold,
    )
    results["single_cell"] = sc_result
    all_flags.extend(sc_result.flags)

    # === Step 2: Filter empty cells ===
    is_empty = sc_result.metrics["is_empty"]
    n_empty = int(is_empty.sum())
    if n_empty > 0:
        logger.info("Filtering %s empty cells for downstream QC", n_empty)
        adata_filtered = adata_per_cell[~is_empty].copy()
    else:
        adata_filtered = adata_per_cell

    # === Step 3: Benchmark-only modules ===
    if mode == "benchmark":
        from .proportion_qc import run_proportion_qc  # pylint: disable=import-outside-toplevel

        logger.info("=== Step 3a: Proportion QC ===")
        prop_result = run_proportion_qc(adata_filtered, proportions, gt_proportions)
        results["proportion"] = prop_result
        all_flags.extend(prop_result.flags)

        if gt_gex_layers is not None:
            from .gex_qc import run_gex_qc  # pylint: disable=import-outside-toplevel

            logger.info("=== Step 3b: GEX QC ===")
            gex_result = run_gex_qc(
                pred_gex_layers=pred_gex_layers if pred_gex_layers is not None else gt_gex_layers,
                gt_gex_layers=gt_gex_layers,
            )
            results["gex"] = gex_result
            all_flags.extend(gex_result.flags)

    # === Step 4: Marker enrichment (both modes) ===
    from .marker_enrichment import run_marker_enrichment  # pylint: disable=import-outside-toplevel

    logger.info("=== Step 4: Marker enrichment ===")
    enrichment_result = run_marker_enrichment(adata_filtered, proportions)
    results["marker_enrichment"] = enrichment_result
    all_flags.extend(enrichment_result.flags)

    # === Step 5: Save artifacts ===
    logger.info("=== Step 5: Saving artifacts ===")

    for module_name, qc_result in results.items():
        for panel_id, fig in qc_result.figures.items():
            fig_path = os.path.join(output_dir, f"{module_name}_{panel_id}.pdf")
            try:
                fig.savefig(fig_path, dpi=150, bbox_inches="tight")
                logger.info("Saved: %s", fig_path)
            except (OSError, ValueError) as e:

                logger.warning("Could not save figure %s: %s", fig_path, e)
            finally:
                plt.close(fig)

    if n_empty > 0:
        empty_df = sc_result.metrics["cell_metrics"][is_empty].copy()
        empty_df["cell_type"] = adata_per_cell.obs.loc[empty_df.index, "cell_type"]
        empty_df.to_csv(os.path.join(output_dir, "empty_cells.csv"))

    summary = {
        "mode": mode,
        "n_cells_total": len(adata_per_cell),
        "n_cells_empty": n_empty,
        "n_cells_analyzed": len(adata_filtered),
        "flags": all_flags,
        "thresholds": {
            "empty_umi_threshold": empty_umi_threshold,
            "empty_genes_threshold": empty_genes_threshold,
        },
    }

    for module_name, qc_result in results.items():
        summary[module_name] = _serialize_metrics(qc_result.metrics)

    with open(os.path.join(output_dir, "qc_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    logger.info("QC complete. %s flags raised. Output: %s", len(all_flags), output_dir)
    return results
