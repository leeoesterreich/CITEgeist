#!/usr/bin/env python
"""GEX evaluation of Tangram and Cell2Location on the breast Xenium benchmark.

Produces results_breast/gex_comparison.json consumed by the breast Fig 2D (σ panel).

Ground-truth strategy:
  The per-cell-type GT CSVs (data_breast/ground_truth_gex/) were saved at the
  INDIVIDUAL CELL level (rows = cells, no spot IDs).  Because cell→spot mappings
  are not persisted, we reconstruct spot-level GT using the standard deconvolution
  identity:

    GT_ct[s, g]  =  prop_ct[s]  ×  total_gex[s, g]

  where prop_ct is the ground-truth proportion of cell type *ct* in spot *s*
  (data_breast/ground_truth/Xenium_region_N_prop.csv) and total_gex is the
  aggregated spot GEX (data_breast/h5ad_objects/Xenium_region_N_GEX.h5ad).
  This is the standard target that deconvolution methods aim to recover.

Alignment: common spots × common genes, log1p-transform both, Pearson r + RMSE.

Output schema (matches _gex_data_for_metric in breast_benchmarking/generate.py):
  {
    "Tangram": {"regions": [{"per_celltype": {"<CellType>": {"rmse": float, "pearson_r": float}}}, ...]},
    "Cell2Location": {"regions": [...5 entries...]}
  }

Usage:
    python compare_breast_gex_methods.py [--output_dir /path/to/results_breast]
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sparse
import anndata as ad
from scipy.stats import pearsonr

REPO_ROOT = Path("/path/to/CITEgeist_analysis")
BENCH_ROOT = REPO_ROOT / "benchmarks" / "xenium_benchmarking"
DATA_ROOT = BENCH_ROOT / "data_breast"

PROP_DIR = DATA_ROOT / "ground_truth"
H5AD_DIR = DATA_ROOT / "h5ad_objects"

METHODS = {
    "Tangram": BENCH_ROOT / "Tangram" / "output_breast",
    "Cell2Location": BENCH_ROOT / "Cell2Location" / "output_breast",
}

BREAST_CELL_TYPES = [
    "Cancer Epithelial",
    "Macrophages",
    "Fibroblasts",
    "Perivascular",
    "Endothelial",
    "T cells",
    "Normal Epithelial",
    "Plasma cells",
]

N_REGIONS = 5
MIN_SPOTS = 2
MIN_EXPRESSED_GENES = 5

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def load_gt_spot(region_id: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load total spot GEX and GT proportions for a region.

    Returns:
        total_gex: DataFrame (spots × genes)
        prop: DataFrame (spots × cell_types)
    """
    h5ad_path = H5AD_DIR / f"Xenium_region_{region_id}_GEX.h5ad"
    prop_path = PROP_DIR / f"Xenium_region_{region_id}_prop.csv"

    adata = ad.read_h5ad(str(h5ad_path))
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    total_gex = pd.DataFrame(X, index=adata.obs.index, columns=adata.var_names)

    prop = pd.read_csv(str(prop_path), index_col=0)
    return total_gex, prop


def compute_gt_ct_spot(total_gex: pd.DataFrame, prop: pd.DataFrame, cell_type: str) -> pd.DataFrame:
    """Construct spot-level GT for one cell type via prop × total GEX.

    Args:
        total_gex: spots × genes raw counts
        prop: spots × cell_types proportions
        cell_type: cell type name

    Returns:
        DataFrame (spots × genes), uppercase gene names, deduped
    """
    if cell_type not in prop.columns:
        return pd.DataFrame()

    common_spots = sorted(set(total_gex.index) & set(prop.index))
    if len(common_spots) < MIN_SPOTS:
        return pd.DataFrame()

    ct_prop = prop.loc[common_spots, cell_type]
    gt = total_gex.loc[common_spots].multiply(ct_prop, axis=0)
    gt.columns = gt.columns.str.upper()
    gt = gt.loc[:, ~gt.columns.duplicated()]
    return gt


def load_pred(method_base_dir: Path, region_id: int, cell_type: str) -> pd.DataFrame:
    """Load prediction layer CSV for one (method, region, cell_type).

    Returns empty DataFrame if file not found.
    """
    layer_path = method_base_dir / f"Xenium_region_{region_id}" / "layers" / f"{cell_type}_layer.csv"
    if not layer_path.exists():
        return pd.DataFrame()

    pred = pd.read_csv(str(layer_path), index_col=0)
    pred.columns = pred.columns.str.upper()
    pred = pred.loc[:, ~pred.columns.duplicated()]
    return pred


def evaluate(pred: pd.DataFrame, gt: pd.DataFrame) -> dict | None:
    """Align pred and GT, apply log1p, compute pearson_r and rmse.

    Returns None if fewer than MIN_SPOTS shared spots or MIN_EXPRESSED_GENES genes.
    """
    common_spots = sorted(set(gt.index) & set(pred.index))
    common_genes = sorted(set(gt.columns) & set(pred.columns))

    if len(common_spots) < MIN_SPOTS or len(common_genes) < MIN_EXPRESSED_GENES:
        return None

    gt_vals = gt.loc[common_spots, common_genes].values.astype(np.float64)
    pred_vals = pred.loc[common_spots, common_genes].values.astype(np.float64)

    # Filter to GT-expressed genes (nonzero total)
    expressed = gt_vals.sum(axis=0) > 0
    n_expressed = int(expressed.sum())
    if n_expressed < MIN_EXPRESSED_GENES:
        return None

    gt_vals = gt_vals[:, expressed]
    pred_vals = pred_vals[:, expressed]

    # Rigor fix: uniform per-block GT-scale normalization for scale-fair RMSE.
    # Methods differ in their reference scale (e.g. Tangram uses raw counts 10-40x
    # higher than GT because its reference was not library-size normalised).
    # Pearson r is scale-invariant and unaffected; RMSE would otherwise be
    # confounded by scale differences across methods.
    # We rescale each (region, method, cell_type) prediction block to match
    # the GT block's total mass with a single global scalar, leaving GT unchanged.
    s_pred = float(pred_vals.sum())
    s_gt = float(gt_vals.sum())
    if s_pred > 0:
        pred_vals = pred_vals * (s_gt / s_pred)

    gt_log = np.log1p(gt_vals)
    pred_log = np.log1p(pred_vals)

    rmse = float(np.sqrt(np.mean((pred_log - gt_log) ** 2)))
    flat_r, _ = pearsonr(pred_log.ravel(), gt_log.ravel())

    return {
        "pearson_r": float(flat_r),
        "rmse": rmse,
        "n_spots": len(common_spots),
        "n_genes": n_expressed,
    }


def main(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "gex_comparison.json"

    results = {}

    for method_name, method_base_dir in METHODS.items():
        if not method_base_dir.exists():
            logger.warning("Method dir missing, skipping: %s (%s)", method_name, method_base_dir)
            continue

        logger.info("=" * 60)
        logger.info("Evaluating: %s", method_name)
        logger.info("=" * 60)

        method_regions = []

        for region_id in range(N_REGIONS):
            logger.info("  Region %d", region_id)
            total_gex, prop = load_gt_spot(region_id)
            region_entry = {"per_celltype": {}}

            for ct in BREAST_CELL_TYPES:
                gt = compute_gt_ct_spot(total_gex, prop, ct)
                if gt.empty:
                    logger.warning("    %s: no GT for region %d", ct, region_id)
                    continue

                pred = load_pred(method_base_dir, region_id, ct)
                if pred.empty:
                    logger.warning("    %s: no pred layer for region %d", ct, region_id)
                    continue

                metrics = evaluate(pred, gt)
                if metrics is None:
                    logger.warning("    %s region %d: insufficient overlap", ct, region_id)
                    continue

                region_entry["per_celltype"][ct] = {
                    "rmse": round(metrics["rmse"], 6),
                    "pearson_r": round(metrics["pearson_r"], 6),
                }
                logger.info(
                    "    %-22s r=%+.4f  rmse=%.4f  (%d spots, %d genes)",
                    ct,
                    metrics["pearson_r"],
                    metrics["rmse"],
                    metrics["n_spots"],
                    metrics["n_genes"],
                )

            method_regions.append(region_entry)
            logger.info("  Region %d: %d cell types evaluated", region_id, len(region_entry["per_celltype"]))

        results[method_name] = {"regions": method_regions}

    with open(str(output_path), "w") as f:
        json.dump(results, f, indent=2)
    logger.info("\nWrote %s", output_path)

    # Sanity summary: mean RMSE across all types and regions per method
    logger.info("\n%s", "=" * 50)
    logger.info("SANITY CHECK — mean RMSE per method (all cell types × regions)")
    logger.info("%s", "=" * 50)
    for method_name, data in results.items():
        all_rmse = []
        all_r = []
        for reg in data["regions"]:
            for ct_metrics in reg["per_celltype"].values():
                if ct_metrics.get("rmse") is not None:
                    all_rmse.append(ct_metrics["rmse"])
                if ct_metrics.get("pearson_r") is not None:
                    all_r.append(ct_metrics["pearson_r"])
        if all_rmse:
            logger.info(
                "  %-16s  mean_rmse=%.4f  mean_r=%+.4f  (n=%d)",
                method_name,
                float(np.mean(all_rmse)),
                float(np.mean(all_r)),
                len(all_rmse),
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Breast GEX evaluation: Tangram and Cell2Location vs prop×total GT")
    parser.add_argument(
        "--output_dir",
        default=str(BENCH_ROOT / "evaluation" / "results_breast"),
        help="Output directory (default: evaluation/results_breast)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(Path(args.output_dir))
