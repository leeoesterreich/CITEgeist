"""
Unified evaluation for CITEgeist Xenium benchmark.

Evaluates both manual (achievable-7) and hierarchical (autodiscovery) modes
against the protein-gated achievable-7 ground truth (7 cell types from
single-cell protein classification).

Ground truth is protein-gated (hierarchical gating at single-cell level
before pseudo-Visium aggregation), so no collapse mapping is needed.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from breast_constants import TISSUE_CONFIGS

logger = logging.getLogger(__name__)

# =============================================================================
# ACHIEVABLE-7 GROUND TRUTH (Protein-Gated)
# =============================================================================
# Protein GT already outputs 7 cell types directly.
# No collapse mapping needed — types match achievable-7 exactly.

ACHIEVABLE_7_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def filter_gt_to_achievable_7(gt_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter protein GT to achievable-7 cell type columns.

    Protein GT already has the correct 7 types. This function simply
    selects those columns and drops metadata columns (n_cells, spot_x, etc.).
    """
    present_types = [ct for ct in ACHIEVABLE_7_CELL_TYPES if ct in gt_df.columns]
    return gt_df[present_types].copy()


def calculate_metrics(
    gt_df: pd.DataFrame,
    pred_df: pd.DataFrame,
    cell_types: List[str],
) -> Dict[str, float]:
    """
    Calculate evaluation metrics for a single region.

    Args:
        gt_df: Ground truth proportions (spots x cell_types)
        pred_df: Predicted proportions (spots x cell_types)
        cell_types: List of cell type column names

    Returns:
        Dict with metrics per cell type and overall
    """
    metrics = {}

    # Per-cell-type metrics
    for ct in cell_types:
        if ct in pred_df.columns and ct in gt_df.columns:
            gt_vals = gt_df[ct].values
            pred_vals = pred_df[ct].values

            # Pearson correlation
            r, p = stats.pearsonr(gt_vals, pred_vals)
            metrics[f"{ct}_pearson_r"] = r
            metrics[f"{ct}_pearson_p"] = p

            # RMSE
            rmse = np.sqrt(np.mean((gt_vals - pred_vals) ** 2))
            metrics[f"{ct}_rmse"] = rmse

            # MAE
            mae = np.mean(np.abs(gt_vals - pred_vals))
            metrics[f"{ct}_mae"] = mae

    # Overall metrics (across all cell types)
    gt_matrix = gt_df[cell_types].values
    pred_matrix = pred_df[[c for c in cell_types if c in pred_df.columns]].values

    # Ensure same shape
    common_types = [c for c in cell_types if c in pred_df.columns]
    gt_matrix = gt_df[common_types].values
    pred_matrix = pred_df[common_types].values

    # Per-spot JSD (Jensen-Shannon Divergence)
    jsd_values = []
    for i in range(len(gt_df)):
        gt_row = gt_matrix[i] + 1e-10  # Add epsilon to avoid log(0)
        pred_row = pred_matrix[i] + 1e-10
        # Normalize to sum to 1
        gt_row = gt_row / gt_row.sum()
        pred_row = pred_row / pred_row.sum()
        jsd = jensenshannon(gt_row, pred_row)
        if not np.isnan(jsd):
            jsd_values.append(jsd)

    metrics["mean_jsd"] = np.mean(jsd_values) if jsd_values else np.nan
    metrics["median_jsd"] = np.median(jsd_values) if jsd_values else np.nan

    # Overall Pearson (flattened)
    gt_flat = gt_matrix.flatten()
    pred_flat = pred_matrix.flatten()
    r_overall, _ = stats.pearsonr(gt_flat, pred_flat)
    metrics["overall_pearson_r"] = r_overall

    # Overall RMSE
    metrics["overall_rmse"] = np.sqrt(np.mean((gt_flat - pred_flat) ** 2))

    # Overall MAE
    metrics["overall_mae"] = np.mean(np.abs(gt_flat - pred_flat))

    return metrics


def combine_t_cells_in_df(df: pd.DataFrame) -> pd.DataFrame:
    """Combine CD4+ and CD8+ T cells into a single 'T cells' column."""
    df = df.copy()
    cd4_cols = [c for c in df.columns if "CD4" in c and "T" in c]
    cd8_cols = [c for c in df.columns if "CD8" in c and "T" in c]

    if cd4_cols and cd8_cols:
        df["T cells"] = df[cd4_cols[0]] + df[cd8_cols[0]]
        df = df.drop(columns=cd4_cols + cd8_cols)

    return df


def load_prediction_with_source(
    pred_dir,
    region_id: int,
    prefix: str = "Xenium",
    combine_t_cells: bool = False,
) -> Tuple[pd.DataFrame, Path]:
    """Load one region's predictions; return (frame, exact file path read).

    Read-only: never writes or stages a copy. The returned path lets callers
    record the sha256 of the file that was actually scored.
    """
    pred_dir = Path(pred_dir)
    pred_path = None
    sample_name = f"{prefix}_region_{region_id}"

    # Pattern 1: Finetuned format at top level (preferred for proportions)
    candidate = pred_dir / f"{sample_name}_cell_prop_finetuned_results.csv"
    if candidate.exists():
        pred_path = candidate

    # Pattern 2: Finetuned format in subdirectory
    if pred_path is None:
        candidate = pred_dir / sample_name / f"{sample_name}_cell_prop_finetuned_results.csv"
        if candidate.exists():
            pred_path = candidate

    # Pattern 3: Deconv predictions format (in subdirectory)
    if pred_path is None:
        candidate = pred_dir / sample_name / f"{sample_name}_deconv_predictions.csv"
        if candidate.exists():
            pred_path = candidate

    # Pattern 4: Global format at top level
    if pred_path is None:
        candidate = pred_dir / f"{sample_name}_cell_prop_global_results.csv"
        if candidate.exists():
            pred_path = candidate

    # Pattern 5: baseline.csv.gz in <sample>/predictions/ (read natively, read-only)
    if pred_path is None:
        candidate = pred_dir / sample_name / "predictions" / "baseline.csv.gz"
        if candidate.exists():
            pred_path = candidate

    if pred_path is None:
        raise FileNotFoundError(f"Predictions not found in {pred_dir}")

    pred_df = pd.read_csv(pred_path, index_col=0)

    # Normalize column names (R's make.names replaces +/spaces with dots)
    col_rename = {
        "CD4..T.cells": "CD4+ T cells",
        "CD8..T.cells": "CD8+ T cells",
        "B.cells": "B cells",
        "T.cells": "T cells",
    }
    pred_df = pred_df.rename(columns=col_rename)

    # Generic dot-to-space fallback for any remaining dotted column names
    dot_rename = {}
    for col in pred_df.columns:
        if "." in col and col not in col_rename.values():
            new_col = col.replace(".", " ")
            if new_col != col:
                dot_rename[col] = new_col
    if dot_rename:
        pred_df = pred_df.rename(columns=dot_rename)

    # Combine T cells if requested (for RNA GT which can't distinguish CD4/CD8)
    if combine_t_cells:
        pred_df = combine_t_cells_in_df(pred_df)
        logger.info("Combined CD4+ and CD8+ T cells into 'T cells'")

    # Normalize to proportions if predictions are cell counts (e.g. Tangram)
    row_sums = pred_df.select_dtypes(include=[np.number]).sum(axis=1)
    if row_sums.max() > 1.5:  # Not already proportions
        logger.info("Normalizing predictions to proportions (detected cell counts)")
        numeric_cols = pred_df.select_dtypes(include=[np.number]).columns
        pred_df[numeric_cols] = pred_df[numeric_cols].div(row_sums, axis=0)

    return pred_df, pred_path


def load_prediction(
    pred_dir,
    region_id: int,
    prefix: str = "Xenium",
    combine_t_cells: bool = False,
) -> pd.DataFrame:
    """Load one region's predictions (see load_prediction_with_source)."""
    return load_prediction_with_source(pred_dir, region_id, prefix, combine_t_cells)[0]


def evaluate_region(
    region_id: int,
    gt_dir: str,
    pred_dir: str,
    prefix: str = "Xenium",
    use_achievable_7: bool = True,
    combine_t_cells: bool = False,
    evaluated_types: List[str] = None,
) -> Dict[str, any]:
    """
    Evaluate a single region.

    Args:
        region_id: Region identifier
        gt_dir: Directory containing ground_truth/ folder
        pred_dir: Directory containing prediction CSVs
        prefix: Filename prefix
        use_achievable_7: Collapse GT to achievable-7 types (default True)
        evaluated_types: Optional list of cell types to evaluate; others become non-evaluated mass

    Returns:
        Dict with region metrics
    """
    gt_dir = Path(gt_dir)
    pred_dir = Path(pred_dir)

    # Load ground truth
    gt_path = gt_dir / "ground_truth" / f"{prefix}_region_{region_id}_prop.csv"
    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")
    gt_df = pd.read_csv(gt_path, index_col=0)

    # Collapse GT to achievable-7 if requested
    if use_achievable_7:
        gt_df = filter_gt_to_achievable_7(gt_df)
        logger.info(f"Filtered GT to achievable-7: {list(gt_df.columns)}")

    pred_df = load_prediction(pred_dir, region_id, prefix, combine_t_cells)

    # Align spots
    common_spots = pred_df.index.intersection(gt_df.index)
    if len(common_spots) == 0:
        raise ValueError(f"No common spots between GT and predictions for region {region_id}")

    pred_aligned = pred_df.loc[common_spots]
    gt_aligned = gt_df.loc[common_spots]

    # Get cell type columns (exclude metadata)
    metadata_cols = ["n_cells", "spot_x", "spot_y", "x", "y", "Unknown"]
    cell_types = [c for c in gt_aligned.columns if c not in metadata_cols]

    # Filter to evaluated types and compute non-evaluated mass
    non_eval_mass = 0.0
    if evaluated_types:
        all_pred_cols = [c for c in pred_aligned.columns if c not in metadata_cols]
        non_eval_cols = [c for c in all_pred_cols if c not in evaluated_types]
        if non_eval_cols:
            non_eval_mass = float(pred_aligned[non_eval_cols].sum(axis=1).mean())
            logger.info(f"Non-evaluated mass (mean): {non_eval_mass:.4f} from columns: {non_eval_cols}")
        eval_pred_cols = [c for c in evaluated_types if c in pred_aligned.columns]
        eval_gt_cols = [c for c in evaluated_types if c in gt_aligned.columns]
        pred_aligned = pred_aligned[eval_pred_cols]
        gt_aligned = gt_aligned[eval_gt_cols]
        cell_types = eval_gt_cols

    logger.info(f"Region {region_id}: {len(common_spots)} spots, {len(cell_types)} cell types")

    # Calculate metrics
    metrics = calculate_metrics(gt_aligned, pred_aligned, cell_types)
    metrics["region_id"] = region_id
    metrics["n_spots"] = len(common_spots)
    metrics["n_cell_types"] = len(cell_types)
    metrics["cell_types"] = cell_types
    metrics["non_eval_mass"] = non_eval_mass

    return metrics


def evaluate_all_regions(
    gt_dir: str,
    pred_dir: str,
    n_regions: int = 5,
    output_path: Optional[str] = None,
    prefix: str = "Xenium",
    use_achievable_7: bool = True,
    combine_t_cells: bool = False,
    evaluated_types: List[str] = None,
) -> Dict[str, any]:
    """
    Evaluate all regions and compute summary statistics.

    Args:
        gt_dir: Directory containing ground_truth/ folder
        pred_dir: Directory containing prediction CSVs
        n_regions: Number of regions
        output_path: Optional path to save results JSON
        prefix: Filename prefix
        use_achievable_7: Filter GT to achievable-7 cell types
        combine_t_cells: Combine CD4+ and CD8+ T cells in predictions
        evaluated_types: Optional list of cell types to evaluate; others become non-evaluated mass

    Returns:
        Dict with all metrics and summary
    """
    all_metrics = []

    for region_id in range(n_regions):
        try:
            metrics = evaluate_region(
                region_id,
                gt_dir,
                pred_dir,
                prefix,
                use_achievable_7,
                combine_t_cells,
                evaluated_types=evaluated_types,
            )
            all_metrics.append(metrics)
            logger.info(f"Region {region_id}: Pearson r = {metrics['overall_pearson_r']:.4f}")
        except Exception as e:
            logger.error(f"Failed to evaluate region {region_id}: {e}")

    if not all_metrics:
        raise ValueError("No regions evaluated successfully")

    # Get cell types from first region
    cell_types = all_metrics[0]["cell_types"]

    # Compute summary statistics
    summary = {
        "n_regions": len(all_metrics),
        "total_spots": sum(m["n_spots"] for m in all_metrics),
        "cell_types": cell_types,
    }

    # Per-cell-type summary
    for ct in cell_types:
        key = f"{ct}_pearson_r"
        values = [m[key] for m in all_metrics if key in m]
        if values:
            summary[f"{ct}_mean_r"] = np.mean(values)
            summary[f"{ct}_std_r"] = np.std(values)

    # Overall summary
    overall_r_values = [m["overall_pearson_r"] for m in all_metrics]
    summary["overall_mean_pearson_r"] = np.mean(overall_r_values)
    summary["overall_std_pearson_r"] = np.std(overall_r_values)

    jsd_values = [m["mean_jsd"] for m in all_metrics if not np.isnan(m["mean_jsd"])]
    if jsd_values:
        summary["overall_mean_jsd"] = np.mean(jsd_values)

    rmse_values = [m["overall_rmse"] for m in all_metrics]
    summary["overall_mean_rmse"] = np.mean(rmse_values)

    mae_values = [m["overall_mae"] for m in all_metrics]
    summary["overall_mean_mae"] = np.mean(mae_values)

    results = {
        "summary": summary,
        "regions": all_metrics,
    }

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert numpy types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_numpy(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(v) for v in obj]
            return obj

        with open(output_path, "w") as f:
            json.dump(convert_numpy(results), f, indent=2)
        logger.info(f"Results saved to {output_path}")

    return results


# =============================================================================
# STRICT CONTRACT (RNA-cluster GT scoring). The legacy path above is untouched.
# No silent skips: every anomaly either raises or is recorded in a report.
# =============================================================================


class MaskDegenerate(Exception):
    """The shared spot mask is unusable (too many GT spots lost, or an
    evaluable type lost all GT variance). Stop and ask the user; never
    silently exclude."""


_GT_CONST_TOL = 1e-12


def _read_gt_strict(gt_dir, region_id: int, prefix: str) -> pd.DataFrame:
    gt_path = Path(gt_dir) / "ground_truth" / f"{prefix}_region_{region_id}_prop.csv"
    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")
    gt_df = pd.read_csv(gt_path, index_col=0)
    if gt_df.index.has_duplicates:
        raise ValueError(f"Duplicate spot ids in GT: {gt_path}")
    return gt_df


def _is_constant(values) -> bool:
    v = np.asarray(values, dtype=float)
    return bool(np.ptp(v) <= _GT_CONST_TOL)


def evaluability_table(gt_dir, types: List[str], n_regions: int = 5, prefix: str = "Xenium") -> pd.DataFrame:
    """Region x type bool: False iff GT is constant on the FULL canonical GT spot set.

    Reads only GT files; call before any prediction is loaded. This is the
    fixed evaluability rule (it never depends on which methods are scored).
    """
    rows = {}
    for region_id in range(n_regions):
        gt_df = _read_gt_strict(gt_dir, region_id, prefix)
        missing = [t for t in types if t not in gt_df.columns]
        if missing:
            raise ValueError(f"GT region {region_id} lacks types {missing}")
        vals = gt_df[list(types)]
        if vals.isna().any().any():
            raise ValueError(f"NaN in GT region {region_id}")
        rows[region_id] = {t: not _is_constant(vals[t].values) for t in types}
    table = pd.DataFrame.from_dict(rows, orient="index")[list(types)]
    table.index.name = "region"
    return table.astype(bool)


def common_spot_masks(
    gt_dir,
    pred_dirs: Dict[str, str],
    evaluable: pd.DataFrame,
    n_regions: int = 5,
    prefix: str = "Xenium",
    max_drop_frac: float = 0.01,
) -> Dict[int, pd.Index]:
    """Per region, GT spots intersected with every method's spots.

    Raises FileNotFoundError/ValueError on a missing region file or an empty
    mask, and MaskDegenerate if the mask drops more than max_drop_frac of the
    GT spots, or if a type evaluable on the full GT is constant on the mask.
    """
    if not pred_dirs:
        raise ValueError("No prediction directories given")
    masks: Dict[int, pd.Index] = {}
    for region_id in range(n_regions):
        gt_df = _read_gt_strict(gt_dir, region_id, prefix)
        mask = gt_df.index
        for method, pdir in pred_dirs.items():
            pred_df = load_prediction(pdir, region_id, prefix)
            if pred_df.index.has_duplicates:
                raise ValueError(f"Duplicate spot ids in predictions of {method}, region {region_id}")
            mask = mask.intersection(pred_df.index)
            if len(mask) == 0:
                raise ValueError(
                    f"No common spots between GT and {method} for region {region_id} "
                    f"(spot-id format mismatch? GT e.g. {list(gt_df.index[:2])}, "
                    f"pred e.g. {list(pred_df.index[:2])})"
                )
        mask = gt_df.index[gt_df.index.isin(mask)]  # keep GT order
        n_gt = len(gt_df)
        dropped = n_gt - len(mask)
        if dropped / n_gt > max_drop_frac:
            raise MaskDegenerate(
                f"Region {region_id}: shared mask drops {dropped}/{n_gt} GT spots "
                f"(> {max_drop_frac:.1%}). Stop and ask the user."
            )
        for t in evaluable.columns:
            if bool(evaluable.loc[region_id, t]) and _is_constant(gt_df.loc[mask, t].values):
                raise MaskDegenerate(
                    f"Region {region_id}: type {t!r} is evaluable on the full GT but has "
                    f"zero GT variance on the shared mask. Stop and ask the user."
                )
        masks[region_id] = mask
    return masks


def _safe_pearson(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return float(np.corrcoef(x, y)[0, 1])


def evaluate_strict(
    gt_dir,
    pred_dir,
    types: List[str],
    masks: Dict[int, pd.Index],
    evaluable: pd.DataFrame,
    n_regions: int = 5,
    prefix: str = "Xenium",
    combine_t_cells: bool = False,
) -> Dict[str, any]:
    """Score one method on the shared masks. Per-type r is the headline metric.

    Per-type r: NaN (excluded) iff GT constant on the full GT; 0.0 if the
    prediction is constant or the type column is absent (zero-filled);
    otherwise Pearson r. Per-type mean over evaluable regions, then over types.
    """
    types = list(types)
    regions: List[Dict[str, any]] = []
    table_rows: List[Dict[str, any]] = []
    mask_report: Dict[int, Dict[str, int]] = {}

    for region_id in range(n_regions):
        if region_id not in masks:
            raise KeyError(f"No mask for region {region_id}")
        mask = masks[region_id]
        gt_df = _read_gt_strict(gt_dir, region_id, prefix)
        pred_df = load_prediction(pred_dir, region_id, prefix, combine_t_cells)
        if not mask.isin(pred_df.index).all() or not mask.isin(gt_df.index).all():
            raise ValueError(f"Mask for region {region_id} contains spots absent from GT or predictions")
        gt = gt_df.loc[mask, types]
        present = [t for t in types if t in pred_df.columns]
        pred = pred_df.loc[mask, present].reindex(columns=types, fill_value=0.0)
        if gt.isna().any().any() or pred.isna().any().any():
            raise ValueError(f"NaN in GT or predictions for region {region_id}")

        mask_report[region_id] = {
            "n_gt_spots": int(len(gt_df)),
            "n_mask_spots": int(len(mask)),
            "spots_dropped": int(len(gt_df) - len(mask)),
        }

        m: Dict[str, any] = {"region_id": region_id, "n_spots": int(len(mask)), "cell_types": types}
        for t in types:
            g = gt[t].values
            p = pred[t].values
            if not bool(evaluable.loc[region_id, t]):
                r, reason = np.nan, "gt_constant_excluded"
            elif t not in present:
                r, reason = 0.0, "missing_type_r0"
            elif _is_constant(p):
                r, reason = 0.0, "pred_constant_r0"
            else:
                r, reason = _safe_pearson(g, p), "ok"
            m[f"{t}_pearson_r"] = r
            m[f"{t}_rmse"] = float(np.sqrt(np.mean((g - p) ** 2)))
            m[f"{t}_mae"] = float(np.mean(np.abs(g - p)))
            table_rows.append(
                {
                    "region": region_id,
                    "type": t,
                    "evaluable": bool(evaluable.loc[region_id, t]),
                    "r": r,
                    "reason": reason,
                }
            )

        gm, pm = gt.values, pred.values
        jsd = []
        for i in range(len(gm)):
            g_row = gm[i] + 1e-10
            p_row = pm[i] + 1e-10
            v = jensenshannon(g_row / g_row.sum(), p_row / p_row.sum())
            if not np.isnan(v):
                jsd.append(v)
        m["mean_jsd"] = float(np.mean(jsd)) if jsd else np.nan
        m["median_jsd"] = float(np.median(jsd)) if jsd else np.nan

        g_flat, p_flat = gm.flatten(), pm.flatten()
        if _is_constant(g_flat):
            raise ValueError(f"Flattened GT is constant in region {region_id}")
        m["overall_pearson_r"] = 0.0 if _is_constant(p_flat) else _safe_pearson(g_flat, p_flat)
        m["overall_rmse"] = float(np.sqrt(np.mean((g_flat - p_flat) ** 2)))
        m["overall_mae"] = float(np.mean(np.abs(g_flat - p_flat)))
        regions.append(m)

    table = pd.DataFrame(table_rows, columns=["region", "type", "evaluable", "r", "reason"])

    summary: Dict[str, any] = {
        "n_regions": len(regions),
        "total_spots": int(sum(m["n_spots"] for m in regions)),
        "cell_types": types,
    }
    per_type_mean_r: Dict[str, float] = {}
    for t in types:
        vals = table.loc[table["type"] == t, "r"].dropna().values
        if len(vals):
            summary[f"{t}_mean_r"] = float(np.mean(vals))
            summary[f"{t}_std_r"] = float(np.std(vals))
            per_type_mean_r[t] = float(np.mean(vals))
        else:  # constant GT in every region: not scorable, reported not hidden
            summary[f"{t}_mean_r"] = np.nan
            summary[f"{t}_std_r"] = np.nan
    flat = [m["overall_pearson_r"] for m in regions]
    summary["overall_mean_pearson_r"] = float(np.mean(flat))
    summary["overall_std_pearson_r"] = float(np.std(flat))
    jsds = [m["mean_jsd"] for m in regions if not np.isnan(m["mean_jsd"])]
    if jsds:
        summary["overall_mean_jsd"] = float(np.mean(jsds))
    summary["overall_mean_rmse"] = float(np.mean([m["overall_rmse"] for m in regions]))
    summary["overall_mean_mae"] = float(np.mean([m["overall_mae"] for m in regions]))
    summary["per_type_mean_r"] = float(np.mean(list(per_type_mean_r.values()))) if per_type_mean_r else np.nan
    summary["n_types_scored"] = len(per_type_mean_r)
    summary["n_zero_filled_missing_type"] = int((table["reason"] == "missing_type_r0").sum())
    summary["n_zero_const_pred"] = int((table["reason"] == "pred_constant_r0").sum())

    return {"summary": summary, "regions": regions, "region_type_table": table, "mask_report": mask_report}


def print_summary(results: Dict[str, any]) -> None:
    """Print a formatted summary of results."""
    summary = results["summary"]

    print("\n" + "=" * 60)
    print("BENCHMARK EVALUATION RESULTS")
    print("=" * 60)

    print(f"\nRegions evaluated: {summary['n_regions']}")
    print(f"Total spots: {summary['total_spots']}")
    print(f"Cell types: {', '.join(summary['cell_types'])}")

    print("\n--- Per-Cell-Type Correlations ---")
    for ct in summary["cell_types"]:
        mean_key = f"{ct}_mean_r"
        std_key = f"{ct}_std_r"
        if mean_key in summary:
            mean_r = summary[mean_key]
            std_r = summary.get(std_key, 0)
            print(f"  {ct:20s}: r = {mean_r:.4f} ± {std_r:.4f}")

    print("\n--- Overall Metrics ---")
    print(f"  Mean Pearson r: {summary['overall_mean_pearson_r']:.4f} ± {summary['overall_std_pearson_r']:.4f}")
    if "overall_mean_jsd" in summary:
        print(f"  Mean JSD:       {summary['overall_mean_jsd']:.4f}")
    print(f"  Mean RMSE:      {summary['overall_mean_rmse']:.4f}")
    print(f"  Mean MAE:       {summary['overall_mean_mae']:.4f}")

    print("\n--- Per-Region Results ---")
    for region in results["regions"]:
        print(
            f"  Region {region['region_id']}: r = {region['overall_pearson_r']:.4f}, "
            f"JSD = {region['mean_jsd']:.4f}, spots = {region['n_spots']}"
        )

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Evaluate CITEgeist benchmark results")
    parser.add_argument(
        "--gt-dir",
        type=str,
        required=False,
        default=None,
        help="Directory containing ground_truth/ folder",
    )
    parser.add_argument(
        "--pred-dir",
        type=str,
        required=False,
        default=None,
        help="Directory containing prediction CSVs",
    )
    parser.add_argument(
        "--n-regions",
        type=int,
        default=5,
        help="Number of regions to evaluate",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON path for results",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="Xenium",
        help="Filename prefix",
    )
    parser.add_argument(
        "--no-achievable-7",
        action="store_true",
        help="Disable achievable-7 GT collapse (use raw 10-type GT)",
    )
    parser.add_argument(
        "--tissue",
        choices=["breast", "rcc"],
        default=None,
        help="Tissue shorthand — resolves --gt-dir, --pred-dir, --evaluated-types",
    )
    parser.add_argument(
        "--method",
        default=None,
        help="Method name (used with --tissue to resolve --pred-dir)",
    )
    parser.add_argument(
        "--evaluated-types",
        default=None,
        help="Comma-separated cell types to evaluate (others become non-evaluated mass)",
    )

    args = parser.parse_args()

    # Resolve tissue shorthand
    if args.tissue:
        cfg = TISSUE_CONFIGS[args.tissue]
        bench_dir = Path(__file__).parent.parent.parent
        if not args.gt_dir:
            args.gt_dir = cfg["data_dir"]
        if not args.pred_dir and args.method:
            args.pred_dir = str(bench_dir / args.method / f"output_{args.tissue}")
        if not args.evaluated_types:
            args.evaluated_types = ",".join(cfg["evaluated_types"])
        args.n_regions = cfg["n_regions"]
        args.prefix = cfg.get("prefix", "Xenium")

    if not args.gt_dir or not args.pred_dir:
        parser.error("--gt-dir and --pred-dir are required (or use --tissue to resolve them)")

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Parse evaluated types
    evaluated_types = None
    if args.evaluated_types:
        evaluated_types = [t.strip() for t in args.evaluated_types.split(",")]

    # Run evaluation
    results = evaluate_all_regions(
        gt_dir=args.gt_dir,
        pred_dir=args.pred_dir,
        n_regions=args.n_regions,
        output_path=args.output,
        prefix=args.prefix,
        use_achievable_7=not args.no_achievable_7,
        evaluated_types=evaluated_types,
    )

    # Print summary
    print_summary(results)


if __name__ == "__main__":
    main()
