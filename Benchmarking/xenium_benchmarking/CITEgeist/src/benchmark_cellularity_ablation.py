#!/usr/bin/env python3
"""
Cellularity ablation: run cuOPT QP for one (region, condition) pair.

Conditions:
  0 = baseline (no nuclei_counts)
  1 = stardist_dapi (StarDist DAPI counts)
  2 = oracle_gt (Xenium GT cell counts)

Caches M1/M2 artifacts per region so only M3 varies.

Usage:
    python -u benchmark_cellularity_ablation.py \
        --region 0 --condition 0 --output-dir output/cellularity_ablation
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "CITEgeist"))

from model.citegeist_model import CitegeistModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

PSEUDOVISIUM = REPO / "Benchmarking/xenium_pseudovisium/data"
STARDIST_DIR = (
    REPO / "Benchmarking/xenium_benchmarking/CITEgeist/output/stardist_counts"
)
CONDITIONS = {
    0: "baseline",
    1: "stardist_dapi",
    2: "oracle_gt",
    3: "oracle_abundance_prior",
    4: "oracle_rowsum_bounds",
}

PROFILE_DICT = {
    "CD4+ T cells": {"Major": ["CD3E"]},
    "CD8+ T cells": {"Major": ["CD3E"]},
    "Macrophages": {"Major": ["CD68"]},
    "Epithelial": {"Major": ["EPCAM"]},
    "Fibroblasts": {"Major": ["VIM"]},
    "Endothelial": {"Major": ["PECAM1"]},
}

EVAL_TYPES = [
    "CD4+ T cells",
    "CD8+ T cells",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
    "Macrophages",
]


def load_nuclei_counts(region: int, condition: int):
    if condition == 0:
        return None
    elif condition == 1:
        csv = STARDIST_DIR / f"stardist_counts_region_{region}.csv"
        if not csv.exists():
            raise FileNotFoundError(f"StarDist counts not found: {csv}")
        s = pd.read_csv(csv, index_col=0).squeeze()
        logger.info("Loaded StarDist counts: %d spots, mean=%.1f", len(s), s.mean())
        return s
    elif condition in (2, 3, 4):
        csv = STARDIST_DIR / f"oracle_counts_region_{region}.csv"
        if not csv.exists():
            raise FileNotFoundError(f"Oracle counts not found: {csv}")
        s = pd.read_csv(csv, index_col=0).squeeze()
        logger.info("Loaded oracle GT counts: %d spots, mean=%.1f", len(s), s.mean())
        return s
    else:
        raise ValueError(f"Unknown condition: {condition}")


def load_profiles(model: CitegeistModel):
    logger.info("Loading hardcoded cell profiles (skipping M1/M2 discovery)...")
    model.load_cell_profile_dict(PROFILE_DICT)


def evaluate(
    pred_df: pd.DataFrame, gt_df: pd.DataFrame, region: int, condition_name: str
) -> dict:
    common_spots = pred_df.index.intersection(gt_df.index)
    common_types = [t for t in EVAL_TYPES if t in pred_df.columns and t in gt_df.columns]

    pred = pred_df.loc[common_spots, common_types].values
    gt = gt_df.loc[common_spots, common_types].values

    # Renormalize both to sum-to-1 over common types (drops Unassigned etc.)
    gt_sums = gt.sum(axis=1, keepdims=True)
    gt_sums = np.maximum(gt_sums, 1e-10)
    gt = gt / gt_sums
    pred_sums = pred.sum(axis=1, keepdims=True)
    pred_sums = np.maximum(pred_sums, 1e-10)
    pred = pred / pred_sums

    flat_r, _ = pearsonr(gt.ravel(), pred.ravel())
    flat_rmse = float(np.sqrt(np.mean((gt - pred) ** 2)))

    per_type = {}
    for j, t in enumerate(common_types):
        gt_col = gt[:, j]
        pred_col = pred[:, j]
        if gt_col.std() < 1e-6:
            per_type[t] = {
                "r": float("nan"),
                "rmse": float(np.sqrt(np.mean((gt_col - pred_col) ** 2))),
            }
        else:
            r, _ = pearsonr(gt_col, pred_col)
            per_type[t] = {
                "r": round(float(r), 4),
                "rmse": round(float(np.sqrt(np.mean((gt_col - pred_col) ** 2))), 4),
            }

    return {
        "region": region,
        "condition": condition_name,
        "n_spots": len(common_spots),
        "n_types": len(common_types),
        "flat_r": round(float(flat_r), 4),
        "flat_rmse": round(flat_rmse, 4),
        "per_type": per_type,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--region", type=int, required=True, choices=range(5))
    parser.add_argument("--condition", type=int, required=True, choices=range(5))
    parser.add_argument("--output-dir", type=str, required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    condition_name = CONDITIONS[args.condition]

    result_file = output_dir / f"result_region{args.region}_{condition_name}.json"
    if result_file.exists():
        logger.info("Skipping (already exists): %s", result_file)
        return

    gex_path = PSEUDOVISIUM / "h5ad_objects" / f"Xenium_region_{args.region}_GEX.h5ad"
    cite_path = PSEUDOVISIUM / "h5ad_objects" / f"Xenium_region_{args.region}_CITE.h5ad"
    gt_path = PSEUDOVISIUM / "ground_truth" / f"Xenium_region_{args.region}_prop.csv"

    gex_adata = sc.read_h5ad(str(gex_path))
    cite_adata = sc.read_h5ad(str(cite_path))
    gt_df = pd.read_csv(gt_path, index_col=0)
    gt_df = gt_df.drop(columns=["n_cells", "spot_x", "spot_y"], errors="ignore")

    logger.info(
        "Region %d: %d spots, %d genes, %d proteins",
        args.region,
        gex_adata.n_obs,
        gex_adata.n_vars,
        cite_adata.n_vars,
    )

    model = CitegeistModel(
        sample_name=f"xenium_region_{args.region}",
        output_folder=str(output_dir / f"region{args.region}_{condition_name}"),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=cite_adata,
    )
    model.preprocess_gex(target_sum=1000)
    model.preprocess_antibody()

    load_profiles(model)

    nuclei_counts = load_nuclei_counts(args.region, args.condition)

    if nuclei_counts is not None:
        spot_names = cite_adata.obs_names
        matched = nuclei_counts.index.intersection(spot_names)
        match_rate = len(matched) / len(spot_names)
        logger.info(
            "Coverage QC: %.1f%% spots matched (%d/%d)",
            match_rate * 100,
            len(matched),
            len(spot_names),
        )
        if match_rate < 0.95:
            raise ValueError(
                f"Coverage QC FAILED: only {match_rate:.1%} spots matched "
                f"(need >=95%). Check barcode alignment."
            )

    logger.info("Running M3 (condition=%s)...", condition_name)

    if args.condition == 3:
        # Abundance prior: inject normalized counts into obs, activate soft penalty
        median_n = float(np.median(nuclei_counts.values))
        normalized = nuclei_counts / median_n
        model.antibody_capture_adata.obs["nuclei_count_target"] = normalized.reindex(
            model.antibody_capture_adata.obs_names, fill_value=1.0
        ).values
        logger.info(
            "Abundance prior: median_N=%.1f, normalized target range=[%.2f, %.2f]",
            median_n, normalized.min(), normalized.max(),
        )
        model.run_cell_proportion_model(
            method="qp",
            nuclei_counts=nuclei_counts,
            max_iterations=50,
            use_nuclei_prior=True,
            nuclei_prior_lambda=1.0,
        )
    elif args.condition == 4:
        # Row-sum bounds: inject normalized counts, use tight per-spot bounds
        median_n = float(np.median(nuclei_counts.values))
        normalized = nuclei_counts / median_n
        model.antibody_capture_adata.obs["nuclei_count_target"] = normalized.reindex(
            model.antibody_capture_adata.obs_names, fill_value=1.0
        ).values
        logger.info(
            "Row-sum bounds: median_N=%.1f, normalized target range=[%.2f, %.2f]",
            median_n, normalized.min(), normalized.max(),
        )
        model.run_cell_proportion_model(
            method="qp",
            nuclei_counts=nuclei_counts,
            max_iterations=50,
            use_nuclei_prior=True,
            nuclei_prior_lambda=5.0,
            nuclei_prior_bounds=(0.3, 3.0),
        )
    else:
        model.run_cell_proportion_model(
            method="qp",
            nuclei_counts=nuclei_counts,
            max_iterations=50,
        )

    prop_df = model.results["cell_prop"]
    logger.info("Proportions: %s", prop_df.shape)

    result = evaluate(prop_df, gt_df, args.region, condition_name)

    if nuclei_counts is not None:
        result["coverage_rate"] = round(match_rate, 4)
        result["nuclei_mean"] = round(
            float(nuclei_counts.reindex(cite_adata.obs_names, fill_value=0).mean()), 2
        )

    with open(result_file, "w") as f:
        json.dump(result, f, indent=2)
    logger.info(
        "Result: flat_r=%.4f, flat_rmse=%.4f -> %s",
        result["flat_r"],
        result["flat_rmse"],
        result_file.name,
    )

    prop_csv = output_dir / f"props_region{args.region}_{condition_name}.csv"
    prop_df.to_csv(prop_csv)
    logger.info("Saved proportions: %s", prop_csv)


if __name__ == "__main__":
    main()
