#!/usr/bin/env python
"""
Benchmark Two-Stage VAE-Guided Assignment.

Runs Stage 1 (hybrid post-filter) + Stage 2 (VAE assignment) and evaluates
single-cell classification accuracy against Xenium ground truth.

Usage:
    python benchmark_two_stage.py --region 0 --output-dir ./output/two_stage
"""
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score, classification_report

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.vae import VAE
from CITEgeist.model.two_stage_pipeline import TwoStagePipeline

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
DATA_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
HYBRID_OUTPUT_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_detection_filter"
VAE_CHECKPOINT = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/vae_sinkhorn_v2/vae/vae_final.pt"
PATCHES_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/patches_v2"

CELL_TYPES = ['B cells', 'CD4+ T cells', 'CD8+ T cells', 'Macrophages',
              'Endothelial', 'Epithelial', 'Fibroblasts']


def load_vae(checkpoint_path: Path, device: str) -> VAE:
    """Load trained VAE from checkpoint."""
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    in_channels = checkpoint.get('in_channels', 2)
    latent_dim = checkpoint.get('latent_dim', 128)

    vae = VAE(in_channels=in_channels, latent_dim=latent_dim)
    vae.load_state_dict(checkpoint['model_state_dict'])
    vae = vae.to(device)
    vae.eval()

    logger.info(f"Loaded VAE: in_channels={in_channels}, latent_dim={latent_dim}")
    return vae


def load_stage1_results(region_id: int) -> tuple:
    """Load Stage 1 (hybrid) results."""
    region_dir = HYBRID_OUTPUT_DIR / f"Xenium_region_{region_id}"

    # Load proportions
    props_path = region_dir / f"Xenium_region_{region_id}_deconv_predictions.csv"
    proportions = pd.read_csv(props_path, index_col=0)

    # Load cell counts
    counts_path = region_dir / f"Xenium_region_{region_id}_cell_counts.csv"
    cell_counts = pd.read_csv(counts_path, index_col=0)

    return proportions, cell_counts


def load_patches_for_region(region_id: int):
    """Create patch loader function for a region."""
    region_patches_dir = PATCHES_DIR / f"Xenium_region_{region_id}" / "spot_patches"

    def load_patches(spot_id: str) -> torch.Tensor:
        patch_file = region_patches_dir / f"{spot_id}_patches.npy"
        if not patch_file.exists():
            return None
        patches = np.load(patch_file).astype(np.float32)
        return torch.from_numpy(patches)

    return load_patches


def load_ground_truth(region_id: int) -> pd.DataFrame:
    """Load Xenium single-cell ground truth."""
    # Ground truth from pseudo-Visium data
    gt_path = DATA_DIR / "single_cell_gt" / f"Xenium_region_{region_id}_cell_types.csv"
    if gt_path.exists():
        return pd.read_csv(gt_path, index_col=0)

    # Alternative: derive from Xenium cell annotations
    logger.warning(f"No ground truth file found at {gt_path}")
    return None


def evaluate_assignments(
    assignments: Dict[str, np.ndarray],
    ground_truth: pd.DataFrame,
    cell_types: list,
) -> dict:
    """Evaluate single-cell assignments against ground truth."""
    y_true = []
    y_pred = []

    for spot_id, pred_indices in assignments.items():
        # Get ground truth for this spot's nuclei
        spot_gt = ground_truth[ground_truth['spot_id'] == spot_id]
        if len(spot_gt) == 0:
            continue

        gt_labels = spot_gt['cell_type'].values

        # Convert predictions to labels
        pred_labels = [cell_types[i] for i in pred_indices[:len(gt_labels)]]

        y_true.extend(gt_labels)
        y_pred.extend(pred_labels)

    if len(y_true) == 0:
        return {'error': 'No ground truth matches'}

    accuracy = accuracy_score(y_true, y_pred)
    report = classification_report(y_true, y_pred, output_dict=True, zero_division=0)

    return {
        'accuracy': accuracy,
        'n_predictions': len(y_true),
        'classification_report': report,
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark Two-Stage Assignment")
    parser.add_argument('--region', type=int, required=True, help='Region ID (0-4)')
    parser.add_argument('--output-dir', type=str, required=True, help='Output directory')
    parser.add_argument('--device', type=str, default='cuda', help='Device')
    parser.add_argument('--skip-training', action='store_true', help='Skip Stage 2 training')
    parser.add_argument('--checkpoint', type=str, default=None, help='Stage 2 checkpoint to load')
    parser.add_argument('--n-epochs', type=int, default=50, help='Training epochs')
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = args.device if torch.cuda.is_available() else 'cpu'
    logger.info(f"Using device: {device}")

    # Load VAE
    vae = load_vae(VAE_CHECKPOINT, device)

    # Load Stage 1 results
    logger.info(f"Loading Stage 1 results for region {args.region}")
    stage1_props, stage1_counts = load_stage1_results(args.region)

    # Create pipeline
    pipeline = TwoStagePipeline(
        vae=vae,
        cell_types=CELL_TYPES,
        device=device,
    )

    # Load patches
    load_patches = load_patches_for_region(args.region)

    # Run Stage 2
    checkpoint_path = Path(args.checkpoint) if args.checkpoint else None

    if not args.skip_training and checkpoint_path is None:
        # Train Stage 2
        pipeline.train_stage2(
            stage1_proportions=stage1_props,
            stage1_cell_counts=stage1_counts,
            load_patches_fn=load_patches,
            n_epochs=args.n_epochs,
            checkpoint_dir=output_dir / f"region_{args.region}",
        )

    assignments = pipeline.run_stage2(
        stage1_proportions=stage1_props,
        stage1_cell_counts=stage1_counts,
        load_patches_fn=load_patches,
        skip_training=args.skip_training,
        checkpoint_path=checkpoint_path,
    )

    # Save assignments
    assignments_path = output_dir / f"region_{args.region}_assignments.json"
    with open(assignments_path, 'w') as f:
        json.dump({k: v.tolist() for k, v in assignments.items()}, f)
    logger.info(f"Saved assignments to {assignments_path}")

    # Evaluate
    ground_truth = load_ground_truth(args.region)
    if ground_truth is not None:
        results = evaluate_assignments(assignments, ground_truth, CELL_TYPES)
        logger.info(f"Accuracy: {results.get('accuracy', 'N/A'):.4f}")
    else:
        results = {'warning': 'No ground truth available for evaluation'}

    # Save results
    results_path = output_dir / f"region_{args.region}_results.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Results saved to {results_path}")


if __name__ == '__main__':
    main()
