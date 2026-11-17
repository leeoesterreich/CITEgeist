#!/usr/bin/env python
"""Spot×layer Module 3.5 validation benchmark on Xenium pseudo-Visium regions.

Evaluates M3.5 functional-gating predictions (per-spot, per-cell-type, per-marker)
against Xenium single-cell protein ground truth aggregated to the same spot×layer
resolution.  Produces per-region JSONs + spot CSVs and an aggregate summary.

Usage:
    # Run one region on a GPU node (sbatch):
    python run_module3_5_validation.py --region 0

    # Aggregate after all 5 regions finish:
    python run_module3_5_validation.py --aggregate
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr
from sklearn.metrics import balanced_accuracy_score, f1_score, precision_score, recall_score
from sklearn.mixture import GaussianMixture

# ---------------------------------------------------------------------------
# sys.path setup — must precede local imports
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent))  # xenium_benchmarking/
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent))  # repo root

from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT  # noqa: E402
from load_xenium_singlecell import load_xenium_singlecell  # noqa: E402

from CITEgeist.model.citegeist_model import CitegeistModel  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Path constants
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium"
EVALUATION_DIR = (
    REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "evaluation" / "module3_5_validation"
)
GT_PROTEIN_CSV = (
    REPO_ROOT
    / "Benchmarking"
    / "xenium_benchmarking"
    / "ground_truth_singler"
    / "singler_labels"
    / "xenium_cell_proteins.csv"
)
CELL_MAPPING_CSV = PSEUDOVISIUM_DIR / "data_protein_gt" / "cell_to_spot_mapping.csv"
CELL_TYPES_CSV = PSEUDOVISIUM_DIR / "data_protein_gt" / "cell_type_assignments.csv"

# ---------------------------------------------------------------------------
# Mapping constants
# ---------------------------------------------------------------------------
GT_TO_6TYPE: Dict[str, str] = {
    "B cells": "B cells",
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}

# CITEgeist M3.5 marker name → Xenium protein column name.
# In practice these are the same strings for this dataset; the mapping exists
# to decouple the two vocabularies if they ever diverge.
FUNCTIONAL_TO_XENIUM_PROTEIN: Dict[str, str] = {
    "PD-1": "PD-1",
    "PD-L1": "PD-L1",
}

# All markers that participate in QP profiles — these are *excluded* from
# M3.5 functional annotation (M3.5 only operates on unused panel markers).
QP_PROFILE_MARKERS = set()
for _ct, _tiers in ACHIEVABLE_6_CELL_PROFILE_DICT.items():
    for _tier_markers in _tiers.values():
        QP_PROFILE_MARKERS.update(_tier_markers)

# Functional marker table: all unused Xenium protein panel markers not in QP profiles.
XENIUM_VALIDATION_TABLE: Dict[str, Dict[str, Any]] = {
    "PD-L1": {"function": "Immune evasion", "active_types": ["Macrophages", "Epithelial"]},
    "PD-1": {"function": "T cell exhaustion", "active_types": ["T cells"]},
    "LAG-3": {"function": "T cell exhaustion", "active_types": ["T cells"]},
    "VISTA": {"function": "Immune checkpoint", "active_types": ["Macrophages", "T cells"]},
    "Vimentin": {"function": "Mesenchymal marker", "active_types": ["Fibroblasts", "Macrophages", "Epithelial"]},
    "Ki-67": {"function": "Proliferation", "active_types": ["Epithelial", "T cells", "Macrophages", "Fibroblasts"]},
    "HLA-DR": {"function": "Antigen presentation", "active_types": ["Macrophages", "B cells"]},
    "CD138": {"function": "Plasma cell marker", "active_types": ["B cells"]},
    "E-Cadherin": {"function": "Epithelial adhesion", "active_types": ["Epithelial"]},
    "Beta-catenin": {"function": "Wnt signaling", "active_types": ["Epithelial"]},
    "CD11c": {"function": "Dendritic/myeloid", "active_types": ["Macrophages"]},
}

MIN_CELLS_PER_LAYER = 3
MIN_SPOTS_FOR_CORRELATION = 10
MIN_REGIONS_FOR_SUMMARY = 3


# ===================================================================
# 1. load_gt_proteins_and_mapping
# ===================================================================
def load_gt_proteins_and_mapping(
    region_id: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load Xenium single-cell protein data and cell-to-spot mapping for a region.

    Returns:
        protein_df: DataFrame indexed by cell_id, columns = protein names,
                    with an extra ``cell_type`` column (6-type labels).
        mapping: DataFrame indexed by cell_id with columns
                 ``spot_id``, ``spot_idx``, ``region_id``, ``cell_type``,
                 ``cell_type_6``.
    """
    _, adata_protein = load_xenium_singlecell(region_id=region_id)
    protein_df = adata_protein.to_df()

    mapping = pd.read_csv(CELL_MAPPING_CSV, index_col=0)
    labels = pd.read_csv(CELL_TYPES_CSV, index_col=0)
    mapping = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)].copy()
    mapping = mapping.join(labels, how="left")
    mapping = mapping[mapping["cell_type"].isin(GT_TO_6TYPE)].copy()
    mapping["cell_type_6"] = mapping["cell_type"].map(GT_TO_6TYPE)

    # Align protein_df to mapped cells
    common = protein_df.index.intersection(mapping.index)
    protein_df = protein_df.loc[common].copy()
    mapping = mapping.loc[common].copy()
    protein_df["cell_type"] = mapping["cell_type_6"]

    logger.info(
        "Region %d: %d GT cells across %d spots",
        region_id,
        len(protein_df),
        mapping["spot_id"].nunique(),
    )
    return protein_df, mapping


# ===================================================================
# 2. run_m35_on_region
# ===================================================================
def run_m35_on_region(region_id: int, output_dir: Path) -> CitegeistModel:
    """Load pseudo-Visium data, run QP + M3.5, return the fitted model."""
    sample_name = f"Xenium_region_{region_id}"
    region_output = output_dir / sample_name
    region_output.mkdir(parents=True, exist_ok=True)

    data_h5ad_dir = PSEUDOVISIUM_DIR / "data_protein_gt" / "h5ad_objects"
    gex_adata = sc.read_h5ad(data_h5ad_dir / f"{sample_name}_GEX.h5ad")
    cite_adata = sc.read_h5ad(data_h5ad_dir / f"{sample_name}_CITE.h5ad")

    # Nuclei counts from cell-to-spot mapping
    mapping_full = pd.read_csv(CELL_MAPPING_CSV, index_col=0)
    labels_full = pd.read_csv(CELL_TYPES_CSV, index_col=0)
    region_map = mapping_full[(mapping_full["region_id"] == region_id) & (mapping_full["spot_idx"] != -1)].copy()
    region_map = region_map.join(labels_full, how="left")
    region_map = region_map[region_map["cell_type"].isin(GT_TO_6TYPE)]
    nuclei_counts = region_map.groupby("spot_id").size()

    logger.info(
        "Loaded pseudo-Visium region %d: %d spots, %d CITE markers",
        region_id,
        gex_adata.n_obs,
        cite_adata.n_vars,
    )

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(region_output),
        simulation=True,
        gene_expression_adata=gex_adata.copy(),
        antibody_capture_adata=cite_adata.copy(),
    )
    model.filter_gex(min_counts=25)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()
    model.load_cell_profile_dict(ACHIEVABLE_6_CELL_PROFILE_DICT)
    model.run_cell_proportion_model(method="qp", lambda_laplacian=0.0, nuclei_counts=nuclei_counts)
    model.run_module3_5_functional_annotation(functional_marker_table=XENIUM_VALIDATION_TABLE)

    logger.info("M3.5 completed for region %d", region_id)
    return model


# ===================================================================
# 3. build_gt_binary_calls_per_cell
# ===================================================================
def build_gt_binary_calls_per_cell(
    protein_df: pd.DataFrame,
    cell_type: str,
    xenium_col: str,
) -> pd.Series:
    """Build binary GT calls for cells of *cell_type* using GMM on log1p protein.

    A 2-component GMM is fitted to log1p(protein) values for cells of the
    given type.  If the two components are well-separated (difference of means
    >= 0.5 in log1p space), the threshold is the midpoint of the component
    means.  Otherwise, fall back to nonzero gating (threshold = 0.0 on raw).

    Returns:
        Series of 0/1 indexed by cell_id (only cells of *cell_type*).
    """
    mask = protein_df["cell_type"] == cell_type
    vals = protein_df.loc[mask, xenium_col].values.astype(float)
    log_vals = np.log1p(vals)

    if len(log_vals) < 20:
        # Too few cells — fall back to nonzero
        return pd.Series(
            (vals > 0.0).astype(int),
            index=protein_df.index[mask],
            name=f"{cell_type}_{xenium_col}_gt",
        )

    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(log_vals.reshape(-1, 1))
    means = gmm.means_.flatten()
    separation = abs(means[1] - means[0])

    if separation >= 0.5:
        threshold = np.mean(means)
        calls = (log_vals >= threshold).astype(int)
    else:
        # Unimodal — fall back to nonzero
        calls = (vals > 0.0).astype(int)

    return pd.Series(
        calls,
        index=protein_df.index[mask],
        name=f"{cell_type}_{xenium_col}_gt",
    )


# ===================================================================
# 4. build_spot_layer_gt
# ===================================================================
def build_spot_layer_gt(
    protein_df: pd.DataFrame,
    mapping: pd.DataFrame,
    evaluable_pairs: List[Tuple[str, str, str]],
) -> pd.DataFrame:
    """Aggregate single-cell GT to spot×layer resolution.

    For each (spot, cell_type, marker) triplet, compute:
      - gt_mean_expr: mean raw protein across cells of that type in that spot
      - gt_positive_frac: fraction of positive cells (by GMM threshold)
      - gt_binary: 1 if gt_positive_frac > 0.5, else 0

    Args:
        protein_df: Single-cell protein data with ``cell_type`` column.
        mapping: Cell-to-spot mapping with ``spot_id`` and ``cell_type_6``.
        evaluable_pairs: List of (cell_type, marker, xenium_col) tuples.

    Returns:
        DataFrame with columns: spot_id, cell_type, marker, gt_mean_expr,
        gt_positive_frac, gt_binary, gt_n_cells.
    """
    rows = []
    for cell_type, marker, xenium_col in evaluable_pairs:
        # Build per-cell binary calls
        gt_calls = build_gt_binary_calls_per_cell(protein_df, cell_type, xenium_col)

        # Get spot assignments for this cell type
        ct_mask = mapping["cell_type_6"] == cell_type
        ct_mapping = mapping.loc[ct_mask, ["spot_id"]].copy()

        # Align to cells that have both mapping and GT calls
        common_cells = ct_mapping.index.intersection(gt_calls.index)
        if len(common_cells) == 0:
            continue
        ct_mapping = ct_mapping.loc[common_cells]
        ct_calls = gt_calls.loc[common_cells]
        ct_expr = protein_df.loc[common_cells, xenium_col]

        # Group by spot
        grouped_expr = ct_expr.groupby(ct_mapping["spot_id"])
        grouped_calls = ct_calls.groupby(ct_mapping["spot_id"])

        for spot_id, expr_group in grouped_expr:
            n_cells = len(expr_group)
            if n_cells < MIN_CELLS_PER_LAYER:
                continue
            call_group = grouped_calls.get_group(spot_id)
            rows.append(
                {
                    "spot_id": spot_id,
                    "cell_type": cell_type,
                    "marker": marker,
                    "gt_mean_expr": float(expr_group.mean()),
                    "gt_positive_frac": float(call_group.mean()),
                    "gt_binary": int(call_group.mean() > 0.5),
                    "gt_n_cells": n_cells,
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=[
                "spot_id",
                "cell_type",
                "marker",
                "gt_mean_expr",
                "gt_positive_frac",
                "gt_binary",
                "gt_n_cells",
            ]
        )
    return pd.DataFrame(rows)


# ===================================================================
# 5. build_spot_layer_pred
# ===================================================================
def build_spot_layer_pred(
    model: CitegeistModel,
    evaluable_pairs: List[Tuple[str, str, str]],
) -> pd.DataFrame:
    """Extract M3.5 predictions at spot×layer resolution.

    Uses ``functional_gates`` (binary calls) and ``functional_lambda``
    (continuous allocation values) from the model results.

    Returns:
        DataFrame with columns: spot_id, cell_type, marker, pred_gate,
        pred_intensity.
    """
    gates_df = model.results.get("functional_gates")
    intensity_df = model.results.get("functional_intensity")

    if gates_df is None:
        logger.warning("No functional_gates in model results")
        return pd.DataFrame(columns=["spot_id", "cell_type", "marker", "pred_gate", "pred_intensity"])

    rows = []
    for cell_type, marker, _xenium_col in evaluable_pairs:
        gate_col = f"{cell_type}_{marker}_gate"
        ratio_col = f"{cell_type}_{marker}_ratio"

        if gate_col not in gates_df.columns:
            logger.debug("Gate column %s not found in model results", gate_col)
            continue

        for spot_id in gates_df.index:
            pred_gate = int(gates_df.loc[spot_id, gate_col])
            pred_intensity = float("nan")
            if intensity_df is not None and ratio_col in intensity_df.columns:
                pred_intensity = float(intensity_df.loc[spot_id, ratio_col])
            rows.append(
                {
                    "spot_id": spot_id,
                    "cell_type": cell_type,
                    "marker": marker,
                    "pred_gate": pred_gate,
                    "pred_intensity": pred_intensity,
                }
            )

    if not rows:
        return pd.DataFrame(columns=["spot_id", "cell_type", "marker", "pred_gate", "pred_intensity"])
    return pd.DataFrame(rows)


# ===================================================================
# 6. get_evaluable_pairs
# ===================================================================
def get_evaluable_pairs(
    model: CitegeistModel,
    protein_df: pd.DataFrame,
) -> List[Tuple[str, str, str]]:
    """Determine which (cell_type, marker, xenium_col) triples can be evaluated.

    A pair is evaluable if:
      1. The marker is in XENIUM_VALIDATION_TABLE
      2. The cell type is in the marker's active_types
      3. The corresponding gate column exists in model results
      4. The Xenium protein column exists in protein_df

    Returns:
        List of (cell_type, marker, xenium_col) tuples.
    """
    gates_df = model.results.get("functional_gates")
    if gates_df is None:
        logger.warning("No functional_gates — no evaluable pairs")
        return []

    pairs = []
    for marker, config in XENIUM_VALIDATION_TABLE.items():
        xenium_col = FUNCTIONAL_TO_XENIUM_PROTEIN.get(marker, marker)
        if xenium_col not in protein_df.columns:
            logger.warning("Xenium protein column %s not found in GT data", xenium_col)
            continue
        for cell_type in config["active_types"]:
            gate_col = f"{cell_type}_{marker}_gate"
            if gate_col in gates_df.columns:
                pairs.append((cell_type, marker, xenium_col))
            else:
                logger.warning(
                    "Gate column %s not in model results — pair not evaluable",
                    gate_col,
                )
    logger.info("Evaluable pairs: %s", pairs)
    return pairs


# ===================================================================
# 7. score_layer_gating
# ===================================================================
def score_layer_gating(
    gt_layers: pd.DataFrame,
    pred_layers: pd.DataFrame,
) -> Dict[str, Any]:
    """Score binary gating: balanced accuracy, precision, recall, F1.

    Merges GT and prediction DataFrames on (spot_id, cell_type, marker)
    and computes classification metrics per (cell_type, marker) pair
    as well as overall.

    Returns:
        Dict with ``per_pair`` (list of dicts) and ``overall`` (dict).
    """
    merged = gt_layers.merge(pred_layers, on=["spot_id", "cell_type", "marker"], how="inner")
    if merged.empty:
        return {"per_pair": [], "overall": {}}

    per_pair = []
    for (ct, mk), grp in merged.groupby(["cell_type", "marker"]):
        if len(grp) < MIN_SPOTS_FOR_CORRELATION:
            continue
        y_true = grp["gt_binary"].values
        y_pred = grp["pred_gate"].values
        # Need both classes present for balanced accuracy
        if len(np.unique(y_true)) < 2:
            ba = float(np.mean(y_true == y_pred))
        else:
            ba = float(balanced_accuracy_score(y_true, y_pred))
        per_pair.append(
            {
                "cell_type": ct,
                "marker": mk,
                "n_spots": len(grp),
                "balanced_accuracy": ba,
                "precision": float(precision_score(y_true, y_pred, zero_division=0)),
                "recall": float(recall_score(y_true, y_pred, zero_division=0)),
                "f1": float(f1_score(y_true, y_pred, zero_division=0)),
                "gt_positive_rate": float(y_true.mean()),
                "pred_positive_rate": float(y_pred.mean()),
            }
        )

    # Overall across all pairs
    overall: Dict[str, Any] = {}
    if len(merged) >= MIN_SPOTS_FOR_CORRELATION:
        y_true_all = merged["gt_binary"].values
        y_pred_all = merged["pred_gate"].values
        if len(np.unique(y_true_all)) < 2:
            ba_all = float(np.mean(y_true_all == y_pred_all))
        else:
            ba_all = float(balanced_accuracy_score(y_true_all, y_pred_all))
        overall = {
            "n_spots": len(merged),
            "balanced_accuracy": ba_all,
            "precision": float(precision_score(y_true_all, y_pred_all, zero_division=0)),
            "recall": float(recall_score(y_true_all, y_pred_all, zero_division=0)),
            "f1": float(f1_score(y_true_all, y_pred_all, zero_division=0)),
        }

    return {"per_pair": per_pair, "overall": overall}


# ===================================================================
# 8. score_layer_allocation
# ===================================================================
def score_layer_allocation(
    gt_layers: pd.DataFrame,
    pred_layers: pd.DataFrame,
) -> Dict[str, Any]:
    """Score continuous allocation: Pearson r and MAE.

    Compares gt_positive_frac (GT) against pred_intensity (M3.5 continuous
    allocation) at the spot×layer level.

    Returns:
        Dict with ``per_pair`` (list of dicts) and ``overall`` (dict).
    """
    merged = gt_layers.merge(pred_layers, on=["spot_id", "cell_type", "marker"], how="inner")
    if merged.empty:
        return {"per_pair": [], "overall": {}}

    # Drop rows with NaN lambda
    merged = merged.dropna(subset=["pred_intensity"]).copy()
    if merged.empty:
        return {"per_pair": [], "overall": {}}

    per_pair = []
    for (ct, mk), grp in merged.groupby(["cell_type", "marker"]):
        if len(grp) < MIN_SPOTS_FOR_CORRELATION:
            continue
        gt_vals = grp["gt_positive_frac"].values
        pred_vals = grp["pred_intensity"].values
        mae = float(np.mean(np.abs(gt_vals - pred_vals)))
        if np.std(gt_vals) > 0 and np.std(pred_vals) > 0:
            r, p = pearsonr(gt_vals, pred_vals)
            r_val, p_val = float(r), float(p)
        else:
            r_val, p_val = float("nan"), float("nan")
        per_pair.append(
            {
                "cell_type": ct,
                "marker": mk,
                "n_spots": len(grp),
                "pearson_r": r_val,
                "pearson_p": p_val,
                "mae": mae,
            }
        )

    overall: Dict[str, Any] = {}
    if len(merged) >= MIN_SPOTS_FOR_CORRELATION:
        gt_all = merged["gt_positive_frac"].values
        pred_all = merged["pred_intensity"].values
        mae_all = float(np.mean(np.abs(gt_all - pred_all)))
        if np.std(gt_all) > 0 and np.std(pred_all) > 0:
            r_all, p_all = pearsonr(gt_all, pred_all)
            r_overall, p_overall = float(r_all), float(p_all)
        else:
            r_overall, p_overall = float("nan"), float("nan")
        overall = {
            "n_spots": len(merged),
            "pearson_r": r_overall,
            "pearson_p": p_overall,
            "mae": mae_all,
        }

    return {"per_pair": per_pair, "overall": overall}


# ===================================================================
# 9. run_region — full pipeline for one region
# ===================================================================
def run_region(region_id: int) -> Dict[str, Any]:
    """Run the full M3.5 spot×layer validation for one region.

    Steps:
      1. Run M3.5 on pseudo-Visium data (QP + functional annotation).
      2. Load Xenium single-cell protein GT + mapping.
      3. Determine evaluable pairs.
      4. Build spot×layer GT and predictions.
      5. Score gating + allocation.
      6. Persist per-region results.

    Returns:
        Dict with gating_scores, allocation_scores, evaluable_pairs, region_id.
    """
    EVALUATION_DIR.mkdir(parents=True, exist_ok=True)
    output_dir = EVALUATION_DIR / "model_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: run M3.5
    model = run_m35_on_region(region_id, output_dir)

    # Step 2: load GT
    protein_df, mapping = load_gt_proteins_and_mapping(region_id)

    # Step 3: evaluable pairs
    evaluable_pairs = get_evaluable_pairs(model, protein_df)
    if not evaluable_pairs:
        logger.warning("No evaluable pairs for region %d", region_id)
        result = {
            "region_id": region_id,
            "evaluable_pairs": [],
            "gating_scores": {"per_pair": [], "overall": {}},
            "allocation_scores": {"per_pair": [], "overall": {}},
        }
        result_path = EVALUATION_DIR / f"region_{region_id}_results.json"
        result_path.write_text(json.dumps(_json_safe(result), indent=2) + "\n")
        return result

    # Step 4: build spot×layer data
    gt_layers = build_spot_layer_gt(protein_df, mapping, evaluable_pairs)
    pred_layers = build_spot_layer_pred(model, evaluable_pairs)

    logger.info(
        "Region %d: %d GT spot×layer rows, %d pred spot×layer rows",
        region_id,
        len(gt_layers),
        len(pred_layers),
    )

    # Step 5: score
    gating_scores = score_layer_gating(gt_layers, pred_layers)
    allocation_scores = score_layer_allocation(gt_layers, pred_layers)

    logger.info("Region %d gating overall: %s", region_id, gating_scores.get("overall"))
    logger.info(
        "Region %d allocation overall: %s",
        region_id,
        allocation_scores.get("overall"),
    )

    # Step 6: persist
    result = {
        "region_id": region_id,
        "evaluable_pairs": [{"cell_type": ct, "marker": mk, "xenium_col": xc} for ct, mk, xc in evaluable_pairs],
        "gating_scores": gating_scores,
        "allocation_scores": allocation_scores,
    }
    result_path = EVALUATION_DIR / f"region_{region_id}_results.json"
    result_path.write_text(json.dumps(_json_safe(result), indent=2) + "\n")
    logger.info("Saved region %d results to %s", region_id, result_path)

    # Save spot-level merged data for downstream analysis
    spot_merged = gt_layers.merge(pred_layers, on=["spot_id", "cell_type", "marker"], how="inner")
    spot_merged["region"] = region_id
    spot_csv = EVALUATION_DIR / f"region_{region_id}_spot_layers.csv"
    spot_merged.to_csv(spot_csv, index=False)
    logger.info("Saved %d spot×layer rows to %s", len(spot_merged), spot_csv)

    return result


# ===================================================================
# 10. _json_safe — numpy serialization helper
# ===================================================================
def _json_safe(obj: Any) -> Any:
    """Recursively convert numpy types to Python natives for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


# ===================================================================
# 11. aggregate_results
# ===================================================================
def aggregate_results() -> Dict[str, Any]:
    """Merge per-region JSONs and spot CSVs into a cross-region summary.

    Reads ``region_{i}_results.json`` and ``region_{i}_spot_layers.csv``
    for i in 0..4 from EVALUATION_DIR.

    Produces:
      - module3_5_validation_summary.json (aggregate metrics)
      - module3_5_spot_layers_all.csv (concatenated spot data)
    """
    all_results = []
    for region_id in range(5):
        result_path = EVALUATION_DIR / f"region_{region_id}_results.json"
        if result_path.exists():
            with open(result_path) as f:
                all_results.append(json.load(f))
            logger.info("Loaded region %d results", region_id)
        else:
            logger.warning("Region %d results not found at %s", region_id, result_path)

    # Merge spot-level CSVs
    spot_dfs = []
    for region_id in range(5):
        spot_csv = EVALUATION_DIR / f"region_{region_id}_spot_layers.csv"
        if spot_csv.exists():
            spot_dfs.append(pd.read_csv(spot_csv))
    if spot_dfs:
        all_spots = pd.concat(spot_dfs, ignore_index=True)
        all_spots_path = EVALUATION_DIR / "module3_5_spot_layers_all.csv"
        all_spots.to_csv(all_spots_path, index=False)
        logger.info(
            "Merged %d spot×layer rows across %d regions to %s",
            len(all_spots),
            len(spot_dfs),
            all_spots_path,
        )
    else:
        all_spots = pd.DataFrame()
        logger.warning("No spot-level CSVs found to merge")

    # Aggregate gating scores across regions
    gating_per_pair_all = []
    allocation_per_pair_all = []
    for res in all_results:
        rid = res["region_id"]
        for pp in res.get("gating_scores", {}).get("per_pair", []):
            pp_copy = dict(pp)
            pp_copy["region_id"] = rid
            gating_per_pair_all.append(pp_copy)
        for pp in res.get("allocation_scores", {}).get("per_pair", []):
            pp_copy = dict(pp)
            pp_copy["region_id"] = rid
            allocation_per_pair_all.append(pp_copy)

    # Compute cross-region means
    gating_summary = {}
    if gating_per_pair_all:
        gating_df = pd.DataFrame(gating_per_pair_all)
        for metric in ["balanced_accuracy", "precision", "recall", "f1"]:
            if metric in gating_df.columns:
                gating_summary[f"mean_{metric}"] = float(gating_df[metric].mean())
                gating_summary[f"std_{metric}"] = float(gating_df[metric].std())

    allocation_summary = {}
    if allocation_per_pair_all:
        alloc_df = pd.DataFrame(allocation_per_pair_all)
        for metric in ["pearson_r", "mae"]:
            if metric in alloc_df.columns:
                vals = alloc_df[metric].dropna()
                if len(vals) > 0:
                    allocation_summary[f"mean_{metric}"] = float(vals.mean())
                    allocation_summary[f"std_{metric}"] = float(vals.std())

    # Re-score on the merged spot data for a pooled overall metric
    pooled_gating = {}
    pooled_allocation = {}
    if not all_spots.empty and "gt_binary" in all_spots.columns and "pred_gate" in all_spots.columns:
        y_true = all_spots["gt_binary"].values
        y_pred = all_spots["pred_gate"].values
        if len(np.unique(y_true)) >= 2:
            pooled_gating = {
                "n_spots": len(all_spots),
                "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
                "precision": float(precision_score(y_true, y_pred, zero_division=0)),
                "recall": float(recall_score(y_true, y_pred, zero_division=0)),
                "f1": float(f1_score(y_true, y_pred, zero_division=0)),
            }
    if not all_spots.empty and "gt_positive_frac" in all_spots.columns and "pred_intensity" in all_spots.columns:
        valid = all_spots.dropna(subset=["pred_intensity"])
        if len(valid) >= MIN_SPOTS_FOR_CORRELATION:
            gt_vals = valid["gt_positive_frac"].values
            pred_vals = valid["pred_intensity"].values
            if np.std(gt_vals) > 0 and np.std(pred_vals) > 0:
                r, p = pearsonr(gt_vals, pred_vals)
                pooled_allocation = {
                    "n_spots": len(valid),
                    "pearson_r": float(r),
                    "pearson_p": float(p),
                    "mae": float(np.mean(np.abs(gt_vals - pred_vals))),
                }

    summary = {
        "n_regions": len(all_results),
        "n_total_spot_layers": len(all_spots),
        "gating_per_pair": gating_per_pair_all,
        "gating_cross_region_means": gating_summary,
        "gating_pooled": pooled_gating,
        "allocation_per_pair": allocation_per_pair_all,
        "allocation_cross_region_means": allocation_summary,
        "allocation_pooled": pooled_allocation,
    }

    summary_path = EVALUATION_DIR / "module3_5_validation_summary.json"
    summary_path.write_text(json.dumps(_json_safe(summary), indent=2) + "\n")
    logger.info("Saved aggregate summary to %s", summary_path)

    # Write per-region CSV for figure generation
    if gating_per_pair_all:
        gating_csv = pd.DataFrame(gating_per_pair_all)
        if allocation_per_pair_all:
            alloc_csv = pd.DataFrame(allocation_per_pair_all)
            per_region_csv = gating_csv.merge(alloc_csv, on=["cell_type", "marker", "region_id"], how="outer")
        else:
            per_region_csv = gating_csv
        per_region_path = EVALUATION_DIR / "module3_5_validation_per_region.csv"
        per_region_csv.to_csv(per_region_path, index=False)
        logger.info("Wrote per-region CSV: %d rows to %s", len(per_region_csv), per_region_path)

        # Pooled summary CSV: mean ± std across regions per (cell_type, marker)
        metric_cols = [
            "balanced_accuracy",
            "precision",
            "recall",
            "f1",
            "pearson_r",
            "mae",
        ]
        existing = [c for c in metric_cols if c in per_region_csv.columns]
        summary_rows = []
        for (ct, mk), grp in per_region_csv.groupby(["cell_type", "marker"]):
            if len(grp) < MIN_REGIONS_FOR_SUMMARY:
                continue
            row = {"cell_type": ct, "marker": mk, "n_regions": len(grp)}
            row["n_layers_total"] = int(grp["n_spots"].sum()) if "n_spots" in grp.columns else 0
            for col in existing:
                vals = grp[col].dropna()
                row[f"{col}_mean"] = float(vals.mean()) if len(vals) > 0 else None
                row[f"{col}_std"] = float(vals.std()) if len(vals) > 1 else None
            summary_rows.append(row)
        summary_csv = pd.DataFrame(summary_rows)
        summary_csv_path = EVALUATION_DIR / "module3_5_validation_summary.csv"
        summary_csv.to_csv(summary_csv_path, index=False)
        logger.info("Wrote summary CSV: %d pairs to %s", len(summary_csv), summary_csv_path)

    return summary


# ===================================================================
# 12. main — argparse entry point
# ===================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="M3.5 spot×layer validation benchmark on Xenium pseudo-Visium")
    parser.add_argument(
        "--region",
        type=int,
        choices=[0, 1, 2, 3, 4],
        default=None,
        help="Region ID to process (0-4). Requires GPU.",
    )
    parser.add_argument(
        "--aggregate",
        action="store_true",
        help="Aggregate per-region results into summary (CPU-only).",
    )
    args = parser.parse_args()

    if args.aggregate:
        aggregate_results()
    elif args.region is not None:
        run_region(args.region)
    else:
        parser.error("Specify --region N or --aggregate")


if __name__ == "__main__":
    main()
