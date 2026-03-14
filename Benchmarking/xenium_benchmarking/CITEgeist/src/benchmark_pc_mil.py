#!/usr/bin/env python
"""
Benchmark PC-MIL on Xenium pseudo-Visium with 5-fold cross-validation.

Pipeline:
1. Load pre-extracted ViT-S features per nucleus
2. Load Module 3 hybrid proportions (conditioning input)
3. Load CLR-normalized protein signals (reconstruction target)
4. Load ground truth proportions + single-cell labels
5. 5-fold CV: train PC-MIL, evaluate proportion r + single-cell accuracy

Usage:
    python benchmark_pc_mil.py \
        --output-dir output/pc_mil \
        --n-epochs 200 \
        --device cuda
"""
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch

# Setup paths
REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.pc_mil import PCMILModel, build_profile_matrix
from CITEgeist.model.pc_mil_training import (
    SpotDataset, compute_inverse_frequency_weights, train_pc_mil,
)
from CITEgeist.model.pc_mil_inference import pc_mil_infer_spot
from CITEgeist.model.detection import detect_cell_types
from CITEgeist.model.constrained_assignment import random_assign

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

# Data paths
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
CELLPOSE_XENIUM_MATCHES = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/vae_sinkhorn_2ch/evaluation_singlecell/cellpose_xenium_matches.csv"
HYBRID_OUTPUT_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_detection_filter"

# Cell types (achievable-7 for Xenium)
CELL_TYPES = [
    "B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages",
    "Endothelial", "Epithelial", "Fibroblasts",
]

# 5 regions, leave-one-out 5-fold CV
FOLDS = [
    {"train": [1, 2, 3, 4], "val": [0]},
    {"train": [0, 2, 3, 4], "val": [1]},
    {"train": [0, 1, 3, 4], "val": [2]},
    {"train": [0, 1, 2, 4], "val": [3]},
    {"train": [0, 1, 2, 3], "val": [4]},
]


def _load_cellpose_xenium_matches() -> pd.DataFrame:
    """Load and cache Cellpose-to-Xenium spatial matching (86% match rate <10px)."""
    if not hasattr(_load_cellpose_xenium_matches, "_cache"):
        _load_cellpose_xenium_matches._cache = pd.read_csv(CELLPOSE_XENIUM_MATCHES)
    return _load_cellpose_xenium_matches._cache


def load_region_data(region_id: int, features_dir: Path) -> Dict:
    """Load all data for one pseudo-Visium region.

    Uses Cellpose nucleus IDs from patches_v2 to align ViT features with
    spatially-matched Xenium GT labels (cellpose_xenium_matches.csv).

    Returns dict with keys:
        features_per_spot: List of (N_i, 384) arrays per spot
        protein_props: (n_spots, K) Module 3 hybrid proportions
        protein_signals: (n_spots, M) CLR-normalized protein signals
        true_props: (n_spots, K) ground truth proportions
        detected: (n_spots, K) boolean detection mask
        cell_type_labels: List of (N_i,) string arrays per spot (GT labels)
        spot_barcodes: List of spot ID strings
    """
    import scanpy as sc

    region_name = f"Xenium_region_{region_id}"

    # 1. Ground truth proportions
    gt_path = PSEUDOVISIUM_DIR / "ground_truth" / f"{region_name}_prop.csv"
    gt_df = pd.read_csv(gt_path, index_col=0)
    spot_ids = list(gt_df.index)
    true_props = gt_df[CELL_TYPES].values  # (n_spots, K)

    # 2. Module 3 hybrid proportions (conditioning input)
    hybrid_dir = HYBRID_OUTPUT_DIR / region_name
    prop_path = hybrid_dir / f"{region_name}_cell_prop_finetuned_results.csv"
    prop_df = pd.read_csv(prop_path, index_col=0)
    prop_df = prop_df.reindex(spot_ids)
    protein_props = prop_df[CELL_TYPES].values  # (n_spots, K)

    # 3. Detection mask
    det_path = hybrid_dir / f"{region_name}_detection_mask.csv"
    det_df = pd.read_csv(det_path, index_col=0)
    det_df = det_df.reindex(spot_ids)
    detected = det_df[CELL_TYPES].values.astype(bool)  # (n_spots, K)

    # 4. CLR-normalized protein signals from CITE h5ad
    cite_path = PSEUDOVISIUM_DIR / "h5ad_objects" / f"{region_name}_CITE.h5ad"
    cite_adata = sc.read_h5ad(cite_path)
    cite_adata = cite_adata[spot_ids, :]
    protein_signals = cite_adata.X.copy()
    if hasattr(protein_signals, 'toarray'):
        protein_signals = protein_signals.toarray()
    protein_signals = np.array(protein_signals, dtype=np.float32)  # (n_spots, M)

    # 5. Cellpose-to-Xenium spatial matches for per-nucleus GT labels
    matches_df = _load_cellpose_xenium_matches()
    region_matches = matches_df[matches_df["region"] == region_id].copy()
    # Build cellpose_id → gt_type lookup (only close matches with known GT)
    close_matches = region_matches[region_matches["distance"] < 15.0]
    cellpose_to_gt = dict(zip(close_matches["cellpose_id"], close_matches["gt_type"]))

    # 6. Cellpose nucleus IDs per spot (from patches_v2)
    nuc_id_dir = REPO_ROOT / f"Benchmarking/xenium_benchmarking/CITEgeist/output/patches_v2/region_{region_id}"

    # 7. ViT features + morphology sizes + aligned GT labels per spot
    features_per_spot = []
    morph_per_spot = []
    cell_type_labels = []

    for spot_id in spot_ids:
        # Load pre-extracted ViT features
        feat_file = features_dir / f"r{region_id}_spot_{spot_id}.npy"
        if feat_file.exists():
            feats = np.load(feat_file)  # (N_i, 384)
        else:
            feats = np.zeros((0, 384), dtype=np.float32)

        features_per_spot.append(feats)

        # Load morphology size features (3-dim: width, height, area proxy)
        size_file = nuc_id_dir / f"spot_{spot_id}_sizes.npy"
        if size_file.exists() and feats.shape[0] > 0:
            sizes = np.load(size_file).astype(np.float32)  # (N_i, 3)
        else:
            sizes = np.zeros((feats.shape[0], 3), dtype=np.float32)
        morph_per_spot.append(sizes)

        # Load Cellpose nucleus IDs (same order as patches/features)
        nuc_id_file = nuc_id_dir / f"spot_{spot_id}_nucleus_ids.npy"
        if nuc_id_file.exists() and feats.shape[0] > 0:
            nuc_ids = np.load(nuc_id_file)
            labels = np.array([
                cellpose_to_gt.get(int(nid), "Unknown") for nid in nuc_ids
            ])
        else:
            labels = np.array(["Unknown"] * feats.shape[0])

        cell_type_labels.append(labels)

    return {
        "features_per_spot": features_per_spot,
        "morph_per_spot": morph_per_spot,
        "protein_props": protein_props,
        "protein_signals": protein_signals,
        "true_props": true_props,
        "detected": detected,
        "cell_type_labels": cell_type_labels,
        "spot_barcodes": spot_ids,
    }


def evaluate_fold(
    model: PCMILModel,
    val_data: Dict,
    cell_type_names: List[str],
    device: str,
) -> Dict:
    """Evaluate model on validation fold.

    Evaluates both:
    - Proportion-level: Pearson r, per-type r, RMSE
    - Single-cell: accuracy using spatially-matched Cellpose-to-Xenium GT labels
      (only evaluates nuclei with known GT, skips "Unknown")
    """
    model.eval()

    all_pred_props = []
    all_true_props = []
    K = len(cell_type_names)

    # Single-cell accuracy tracking
    correct = 0
    total = 0
    per_type = {ct: {"correct": 0, "total": 0} for ct in cell_type_names}

    has_morph = "morph_per_spot" in val_data

    for i in range(len(val_data["features_per_spot"])):
        feats = val_data["features_per_spot"][i]
        if len(feats) == 0:
            continue

        img_feats = torch.tensor(feats, dtype=torch.float32)
        prot_props = torch.tensor(val_data["protein_props"][i], dtype=torch.float32)

        morph_feats = None
        if has_morph:
            morph_feats = torch.tensor(val_data["morph_per_spot"][i], dtype=torch.float32)

        # Detection mask
        detected = val_data.get("detected", np.ones(K, dtype=bool))
        if detected.ndim == 2:
            detected = detected[i]

        result_df = pc_mil_infer_spot(
            model=model,
            image_features=img_feats,
            protein_proportions=prot_props,
            detected_types=detected,
            cell_type_names=cell_type_names,
            morph_features=morph_feats,
        )

        # Proportion evaluation from discrete assignments
        type_counts = result_df["cell_type"].value_counts()
        pred_p = np.zeros(K)
        for k, ct in enumerate(cell_type_names):
            pred_p[k] = type_counts.get(ct, 0) / len(result_df)
        all_pred_props.append(pred_p)
        all_true_props.append(val_data["true_props"][i])

        # Single-cell accuracy (using spatially-matched GT labels)
        if "cell_type_labels" in val_data:
            gt_labels = val_data["cell_type_labels"][i]
            pred_labels = result_df["cell_type"].values
            for j in range(min(len(gt_labels), len(pred_labels))):
                gt = gt_labels[j]
                if gt == "Unknown" or gt not in cell_type_names:
                    continue  # Skip unmatched/unknown nuclei
                total += 1
                pred = pred_labels[j]
                if gt == pred:
                    correct += 1
                per_type[gt]["total"] += 1
                if gt == pred:
                    per_type[gt]["correct"] += 1

    # Stack for vectorized metrics
    pred_arr = np.array(all_pred_props)  # (n_spots, K)
    true_arr = np.array(all_true_props)  # (n_spots, K)

    # Overall Pearson r (flattened)
    pred_flat = pred_arr.flatten()
    true_flat = true_arr.flatten()
    if pred_flat.std() > 0 and true_flat.std() > 0:
        prop_r = float(np.corrcoef(pred_flat, true_flat)[0, 1])
    else:
        prop_r = 0.0

    # Per-type Pearson r
    per_type_r = {}
    for k, ct in enumerate(cell_type_names):
        p = pred_arr[:, k]
        t = true_arr[:, k]
        if p.std() > 0 and t.std() > 0:
            per_type_r[ct] = float(np.corrcoef(p, t)[0, 1])
        else:
            per_type_r[ct] = 0.0

    # RMSE
    rmse = float(np.sqrt(np.mean((pred_arr - true_arr) ** 2)))

    # Single-cell accuracy
    sc_accuracy = correct / total if total > 0 else 0.0
    per_type_acc = {}
    for ct in cell_type_names:
        t = per_type[ct]["total"]
        c = per_type[ct]["correct"]
        per_type_acc[ct] = c / t if t > 0 else 0.0

    return {
        "proportion_r": prop_r,
        "per_type_r": per_type_r,
        "rmse": rmse,
        "single_cell_accuracy": sc_accuracy,
        "per_type_accuracy": per_type_acc,
        "n_spots": len(all_pred_props),
        "n_evaluated_cells": total,
    }


def main():
    parser = argparse.ArgumentParser(description="PC-MIL Xenium Benchmark")
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--features-dir", type=str, required=True,
                        help="Directory with pre-extracted ViT features per region")
    parser.add_argument("--n-epochs", type=int, default=200)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--lambda-recon", type=float, default=1.0)
    parser.add_argument("--lambda-entropy", type=float, default=0.1)
    parser.add_argument("--lambda-diversity", type=float, default=0.5)
    parser.add_argument("--lambda-hungarian", type=float, default=1.0)
    parser.add_argument("--patience", type=int, default=30)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--fold", type=int, default=None, help="Run single fold (0-4)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    features_dir = Path(args.features_dir)

    folds_to_run = [args.fold] if args.fold is not None else range(5)
    all_results = []

    for fold_idx in folds_to_run:
        fold = FOLDS[fold_idx]
        logger.info(f"=== Fold {fold_idx}: train={fold['train']}, val={fold['val']} ===")

        # Load data for train/val regions
        train_features, train_morph, train_props, train_signals, train_true, train_detected = [], [], [], [], [], []
        for rid in fold["train"]:
            data = load_region_data(rid, features_dir)
            train_features.extend(data["features_per_spot"])
            train_morph.extend(data["morph_per_spot"])
            train_props.append(data["protein_props"])
            train_signals.append(data["protein_signals"])
            train_true.append(data["true_props"])
            train_detected.append(data["detected"])

        train_props = np.concatenate(train_props)
        train_signals = np.concatenate(train_signals)
        train_true = np.concatenate(train_true)
        train_detected_arr = np.concatenate(train_detected)

        # Inverse frequency weights
        weights = compute_inverse_frequency_weights(train_true)

        train_ds = SpotDataset(
            train_features, train_props, train_signals, train_true, weights,
            detected_types=train_detected_arr,
            morph_per_spot=train_morph,
        )

        # Val — collect per-spot lists
        val_features, val_morph = [], []
        val_props_list, val_signals_list, val_true_list, val_detected_list = [], [], [], []
        val_labels_list = []
        for rid in fold["val"]:
            data = load_region_data(rid, features_dir)
            val_features.extend(data["features_per_spot"])
            val_morph.extend(data["morph_per_spot"])
            val_props_list.append(data["protein_props"])
            val_signals_list.append(data["protein_signals"])
            val_true_list.append(data["true_props"])
            val_detected_list.append(data["detected"])
            val_labels_list.extend(data["cell_type_labels"])

        val_props = np.concatenate(val_props_list)
        val_true = np.concatenate(val_true_list)
        val_signals = np.concatenate(val_signals_list)
        val_detected = np.concatenate(val_detected_list)

        val_ds = SpotDataset(
            val_features, val_props, val_signals, val_true,
            detected_types=val_detected,
            morph_per_spot=val_morph,
        )

        # Build val_data dict for evaluate_fold (per-spot indexing)
        val_data_combined = {
            "features_per_spot": val_features,
            "morph_per_spot": val_morph,
            "protein_props": val_props,
            "true_props": val_true,
            "detected": val_detected,
            "cell_type_labels": val_labels_list,
        }

        # Model (with morphology branch: 3-dim size features)
        K = len(CELL_TYPES)
        M = train_signals.shape[1]
        model = PCMILModel(
            image_dim=384, n_types=K, n_markers=M, morph_dim=3,
        )

        save_path = str(output_dir / f"pc_mil_fold{fold_idx}.pt")

        history = train_pc_mil(
            model=model,
            train_dataset=train_ds,
            val_dataset=val_ds,
            n_epochs=args.n_epochs,
            lr=args.lr,
            lambda_recon=args.lambda_recon,
            lambda_entropy=args.lambda_entropy,
            lambda_diversity=args.lambda_diversity,
            lambda_hungarian=args.lambda_hungarian,
            patience=args.patience,
            device=args.device,
            save_path=save_path,
        )

        # Evaluate best model
        model.load_state_dict(torch.load(save_path, weights_only=True))
        model.to(args.device)

        fold_results = evaluate_fold(
            model, val_data_combined, CELL_TYPES, args.device,
        )
        fold_results["fold"] = fold_idx
        fold_results["history"] = history
        all_results.append(fold_results)

        logger.info(f"Fold {fold_idx}: prop_r={fold_results['proportion_r']:.4f}, "
                     f"rmse={fold_results['rmse']:.4f}, "
                     f"sc_acc={fold_results['single_cell_accuracy']:.4f} "
                     f"({fold_results['n_evaluated_cells']} cells)")

        # Save fold results
        with open(output_dir / f"fold{fold_idx}_results.json", "w") as f:
            json.dump({k: v for k, v in fold_results.items() if k != "history"},
                      f, indent=2, default=str)

    # Summary
    if all_results:
        mean_r = np.mean([r["proportion_r"] for r in all_results])
        mean_rmse = np.mean([r["rmse"] for r in all_results])
        mean_sc_acc = np.mean([r["single_cell_accuracy"] for r in all_results])
        total_cells = sum(r["n_evaluated_cells"] for r in all_results)
        logger.info(f"=== SUMMARY: mean_prop_r={mean_r:.4f}, mean_rmse={mean_rmse:.4f}, "
                     f"mean_sc_acc={mean_sc_acc:.4f} ({total_cells} cells) ===")

        # Aggregate per-type metrics across folds
        agg_per_type_r = {}
        agg_per_type_acc = {}
        for ct in CELL_TYPES:
            agg_per_type_r[ct] = float(np.mean([r["per_type_r"].get(ct, 0.0) for r in all_results]))
            agg_per_type_acc[ct] = float(np.mean([r["per_type_accuracy"].get(ct, 0.0) for r in all_results]))

        summary = {
            "mean_proportion_r": float(mean_r),
            "mean_rmse": float(mean_rmse),
            "mean_single_cell_accuracy": float(mean_sc_acc),
            "n_evaluated_cells": total_cells,
            "per_type_r": agg_per_type_r,
            "per_type_accuracy": agg_per_type_acc,
            "folds": [{k: v for k, v in r.items() if k != "history"} for r in all_results],
        }
        with open(output_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)


if __name__ == "__main__":
    main()
