#!/usr/bin/env python
"""Run benchmark-mode QC on CITEgeist Xenium results vs SingleR ground truth.

For each Xenium region (0–4):
  - Loads CITEgeist QP proportions vs SingleR GT proportions
  - Loads SACE per-type GEX layers vs GT GEX layers
  - Constructs a pseudo-cell AnnData (one cell per spot, dominant GT type)
  - Runs run_qc() in benchmark mode → PDFs + qc_summary.json

Then prints a cross-region aggregate table of proportion and GEX correlations.

Usage:
    python run_qc_xenium_benchmark.py [--region 0]  # single region test
    python run_qc_xenium_benchmark.py               # all 5 regions
"""

import argparse
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Repo root → model package
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "CITEgeist"))

import json
import numpy as np
import pandas as pd
import anndata as ad

from model import run_qc

# ── Paths ──────────────────────────────────────────────────────────────────
BENCH = Path(__file__).resolve().parents[2]
GT_PROP_DIR = BENCH / "ground_truth_singler/ground_truth_7type/ground_truth"
GT_GEX_DIR = BENCH / "ground_truth_singler/ground_truth_7type/ground_truth_gex"
QP_PRED_DIR = BENCH / "CITEgeist/output_qp_singler_7type"
SACE_DIR = BENCH / "CITEgeist/results_sace_singler_7type"
QC_OUT_DIR = BENCH / "CITEgeist/qc_benchmark"

REGIONS = [0, 1, 2, 3, 4]

# ── Cell-type name normalisation ───────────────────────────────────────────
# Map from file/column names → CANONICAL_MARKERS keys used by QC framework
_RENAME = {
    "B cells": "B_Cells",
    "B_cells": "B_Cells",
    "CD4+ T cells": "CD4_T_Cells",
    "CD4pos_T_cells": "CD4_T_Cells",
    "CD8+ T cells": "CD8_T_Cells",
    "CD8pos_T_cells": "CD8_T_Cells",
    "Macrophages": "Macrophages",
    "Fibroblasts": "Fibroblasts",
    "Epithelial": "Epithelial",
    "Endothelial": "Endothelial",
}
_PROP_EXTRA = {"n_cells", "spot_x", "spot_y"}


def _load_proportions(csv_path: Path) -> pd.DataFrame:
    """Load a proportion CSV, drop non-type columns, rename, normalise rows."""
    df = pd.read_csv(csv_path, index_col=0)
    drop_cols = [c for c in df.columns if c in _PROP_EXTRA]
    df = df.drop(columns=drop_cols)
    df = df.rename(columns=_RENAME)
    # Row-normalise to sum=1 (guard against floating point drift)
    row_sums = df.sum(axis=1).replace(0, np.nan)
    df = df.div(row_sums, axis=0).fillna(0)
    return df


def _load_gex_layers(region_dir: Path, suffix: str) -> dict[str, pd.DataFrame]:
    """Load per-type GEX layer CSVs from a directory.

    Args:
        region_dir: directory containing *_<suffix>.csv files.
        suffix: file suffix after the type name (e.g. "GT" or "layer").

    Returns:
        Dict mapping CANONICAL_MARKERS key → DataFrame (spots × genes).
    """
    layers = {}
    for f in sorted(region_dir.glob(f"*_{suffix}.csv")):
        raw_name = f.stem[: -(len(suffix) + 1)]  # strip "_GT" or "_layer"
        canon = _RENAME.get(raw_name)
        if canon is None:
            logger.warning("No canonical mapping for '%s', skipping", raw_name)
            continue
        df = pd.read_csv(f, index_col=0)
        layers[canon] = df
    return layers


def _build_adata(
    gt_props: pd.DataFrame, gt_gex_layers: dict[str, pd.DataFrame], region_label: str, gt_raw_csv: Path
) -> ad.AnnData:
    """Construct a pseudo-cell AnnData (one spot = one pseudo-cell).

    cell_type = dominant GT cell type per spot (argmax of GT proportions).
    X         = sum of all GT GEX layers (reconstructed spot-level GEX).
    spatial   = spot_x, spot_y from original GT CSV.
    """
    # Dominant GT type per spot
    dominant = gt_props.idxmax(axis=1)  # Series: spot_id → canonical type

    # X: sum GT layers → spots × genes reconstructed GEX
    all_dfs = list(gt_gex_layers.values())
    common_spots = gt_props.index
    common_genes = all_dfs[0].columns if all_dfs else pd.Index([])
    for df in all_dfs[1:]:
        common_genes = common_genes.intersection(df.columns)

    spot_gex = pd.DataFrame(0.0, index=common_spots, columns=common_genes)
    for df in all_dfs:
        shared_spots = common_spots.intersection(df.index)
        spot_gex.loc[shared_spots, :] += df.loc[shared_spots, common_genes].values

    adata = ad.AnnData(X=spot_gex.values.astype(np.float32), obs=pd.DataFrame(index=common_spots))
    adata.obs["cell_type"] = dominant.reindex(common_spots).values
    adata.obs["spot_id"] = common_spots.values
    adata.var_names = common_genes

    # Spatial coords from original GT prop CSV
    gt_raw = pd.read_csv(gt_raw_csv, index_col=0)
    if "spot_x" in gt_raw.columns and "spot_y" in gt_raw.columns:
        coords = gt_raw.loc[common_spots, ["spot_x", "spot_y"]].values
        adata.obsm["spatial"] = coords.astype(np.float32)

    adata.obs["region"] = region_label
    return adata


def run_qc_region(region: int) -> dict:
    """Run benchmark QC for one Xenium region. Returns parsed summary dict."""
    tag = f"Xenium_region_{region}"
    out_dir = QC_OUT_DIR / tag
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / "qc_summary.json"
    if summary_path.exists():
        logger.info("QC already done for %s, loading cached summary", tag)
        with open(summary_path) as f:
            return json.load(f)

    logger.info("=== Running benchmark QC: %s ===", tag)

    # --- Proportions ---
    gt_prop_csv = GT_PROP_DIR / f"{tag}_prop.csv"
    qp_pred_csv = QP_PRED_DIR / tag / f"{tag}_deconv_predictions.csv"
    gt_props = _load_proportions(gt_prop_csv)
    qp_props = _load_proportions(qp_pred_csv)

    # Align to common spots
    common_spots = gt_props.index.intersection(qp_props.index)
    gt_props = gt_props.loc[common_spots]
    qp_props = qp_props.loc[common_spots]

    # Align to common types
    common_types = gt_props.columns.intersection(qp_props.columns)
    gt_props = gt_props[common_types]
    qp_props = qp_props[common_types]
    logger.info("%s: %d spots, %d types", tag, len(common_spots), len(common_types))

    # --- GEX layers ---
    sace_layer_dir = SACE_DIR / tag / "layers"
    gt_gex_dir = GT_GEX_DIR / tag
    pred_gex = _load_gex_layers(sace_layer_dir, "layer")
    gt_gex = _load_gex_layers(gt_gex_dir, "GT")

    # Keep only types present in both pred and GT GEX
    gex_types = set(pred_gex) & set(gt_gex)
    pred_gex = {k: v for k, v in pred_gex.items() if k in gex_types}
    gt_gex = {k: v for k, v in gt_gex.items() if k in gex_types}

    # Align spot indices across GEX layers (intersect all three: props, pred, GT)
    for ct in list(pred_gex):
        shared = common_spots.intersection(pred_gex[ct].index).intersection(gt_gex[ct].index)
        pred_gex[ct] = pred_gex[ct].loc[shared]
        gt_gex[ct] = gt_gex[ct].loc[shared]

    # --- Pseudo-cell AnnData ---
    adata = _build_adata(gt_props, gt_gex, tag, gt_prop_csv)
    # Align adata to common_spots
    shared = common_spots.intersection(adata.obs_names)
    adata = adata[shared].copy()
    gt_props = gt_props.loc[shared]
    qp_props = qp_props.loc[shared]

    # --- Run QC ---
    run_qc(
        adata_per_cell=adata,
        proportions=qp_props,
        mode="benchmark",
        gt_proportions=gt_props,
        gt_gex_layers=gt_gex,
        pred_gex_layers=pred_gex,
        output_dir=str(out_dir),
    )

    with open(summary_path) as f:
        return json.load(f)


def print_aggregate_table(summaries: list[dict], regions: list[int]) -> None:
    """Print per-region + mean proportion and GEX correlation table."""
    rows = []
    for region, s in zip(regions, summaries):
        prop_corr = s.get("proportion", {}).get("correlations", {})
        gex_corr = s.get("gex", {}).get("correlations", {})
        overall_r = prop_corr.get("overall", {})
        prop_r = overall_r.get("pearson_r", float("nan")) if isinstance(overall_r, dict) else float("nan")
        gex_vals = [v["pearson_r"] for v in gex_corr.values() if isinstance(v, dict) and "pearson_r" in v]
        mean_gex_r = float(np.mean(gex_vals)) if gex_vals else float("nan")
        rows.append({"region": region, "prop_r_overall": prop_r, "mean_gex_r": mean_gex_r})

    df = pd.DataFrame(rows)
    print("\n=== XENIUM BENCHMARK — AGGREGATE SUMMARY ===")
    print(df.to_string(index=False, float_format="{:.3f}".format))
    print(f"\nMean across regions:  prop_r={df['prop_r_overall'].mean():.3f}  gex_r={df['mean_gex_r'].mean():.3f}")

    # Per-type proportion breakdown (region 0 as representative)
    if summaries:
        s0 = summaries[0]
        corr = s0.get("proportion", {}).get("correlations", {})
        print("\nPer-type proportion r (region 0):")
        for ct, v in sorted(corr.items()):
            if isinstance(v, dict) and "pearson_r" in v:
                print(f"  {ct:20s}: {v['pearson_r']:.3f}")
        gex_corr = s0.get("gex", {}).get("correlations", {})
        if gex_corr:
            print("\nPer-type GEX r (region 0):")
            for ct, v in sorted(gex_corr.items()):
                if isinstance(v, dict) and "pearson_r" in v:
                    print(f"  {ct:20s}: {v['pearson_r']:.3f}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--region", type=int, default=None, help="Single region index to process (0-4). Omit for all.")
    args = parser.parse_args()

    QC_OUT_DIR.mkdir(parents=True, exist_ok=True)

    regions = [args.region] if args.region is not None else REGIONS
    summaries = []
    for r in regions:
        s = run_qc_region(r)
        summaries.append(s)

    if len(summaries) > 1:
        print_aggregate_table(summaries, regions)
    else:
        # Single region: just print flags
        s = summaries[0]
        print(f"Flags: {s.get('flags', [])}")


if __name__ == "__main__":
    main()
