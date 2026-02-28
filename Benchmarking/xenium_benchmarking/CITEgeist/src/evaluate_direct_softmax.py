#!/usr/bin/env python3
"""
Evaluate DirectSoftmaxModel for single-cell type assignment.

Maps predictions to Xenium ground truth cells and computes accuracy.
"""
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import sys

import numpy as np
import pandas as pd
import torch
from scipy.spatial import cKDTree
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from tqdm import tqdm

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.vae import VAEEncoder
from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Xenium data paths
XENIUM_DATA_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium"


def load_vae_encoder(checkpoint_path: str, device: torch.device) -> Tuple[VAEEncoder, int]:
    """Load frozen VAE encoder from checkpoint."""
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    in_channels = checkpoint.get("in_channels", 2)
    latent_dim = checkpoint.get("latent_dim", 128)

    encoder = VAEEncoder(in_channels=in_channels, latent_dim=latent_dim)
    state_dict = checkpoint["model_state_dict"]
    encoder_state = {
        k.replace("encoder.", ""): v
        for k, v in state_dict.items()
        if k.startswith("encoder.")
    }
    encoder.load_state_dict(encoder_state)
    encoder = encoder.to(device)
    encoder.eval()

    return encoder, latent_dim


def load_direct_softmax_model(
    encoder: VAEEncoder,
    checkpoint_path: str,
    device: torch.device,
) -> DirectSoftmaxModel:
    """Load trained DirectSoftmaxModel."""
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)

    n_types = checkpoint["n_types"]
    latent_dim = checkpoint["latent_dim"]
    projection_dim = checkpoint["projection_dim"]
    temperature = checkpoint.get("temperature", 0.1)

    # Load optional new parameters from checkpoint
    enable_unknown = checkpoint.get("enable_unknown", False)
    unknown_threshold = checkpoint.get("unknown_threshold", 0.15)
    use_attention = checkpoint.get("use_attention", False)
    use_per_class_attention = checkpoint.get("use_per_class_attention", False)
    attention_confidence_bias = checkpoint.get("attention_confidence_bias", True)
    use_size_features = checkpoint.get("use_size_features", False)

    model = DirectSoftmaxModel(
        encoder=encoder,
        n_types=n_types,
        latent_dim=latent_dim,
        projection_dim=projection_dim,
        temperature=temperature,
        enable_unknown=enable_unknown,
        unknown_threshold=unknown_threshold,
        use_attention=use_attention,
        use_per_class_attention=use_per_class_attention,
        attention_confidence_bias=attention_confidence_bias,
        use_size_features=use_size_features,
    )

    # Load projection head
    model.projection.load_state_dict(checkpoint["projection_state_dict"])

    # Load prototypes
    model.prototypes.data = checkpoint["prototypes"].to(device)

    # Load attention aggregator if present
    if use_attention and "aggregator_state_dict" in checkpoint:
        model.aggregator.load_state_dict(checkpoint["aggregator_state_dict"])

    model = model.to(device)
    model.eval()

    logger.info(f"Loaded DirectSoftmaxModel: n_types={n_types}, projection_dim={projection_dim}")
    logger.info(f"  enable_unknown={enable_unknown}, use_attention={use_attention}, use_size_features={use_size_features}")
    logger.info(f"Cell types: {checkpoint.get('celltype_cols', 'unknown')}")

    return model, checkpoint.get("celltype_cols", [f"Type_{i}" for i in range(n_types)]), use_size_features


def load_xenium_ground_truth() -> pd.DataFrame:
    """Load Xenium cells with ground truth cell types."""
    # Load cell coordinates
    cells_file = XENIUM_DATA_DIR / "cells.parquet"
    if cells_file.exists():
        cells_df = pd.read_parquet(cells_file)
    else:
        cells_file = XENIUM_DATA_DIR / "cells.csv.gz"
        cells_df = pd.read_csv(cells_file, compression='gzip')

    # Load ground truth assignments
    gt_file = PSEUDOVISIUM_DIR / "data_protein_gt" / "cell_type_assignments.csv"
    gt_df = pd.read_csv(gt_file, index_col=0)

    # Merge
    cells_df = cells_df.set_index('cell_id')
    cells_df = cells_df.join(gt_df)

    return cells_df.reset_index()


def run_inference_on_spot(
    model: DirectSoftmaxModel,
    patches: torch.Tensor,
    proportions: torch.Tensor,
    use_hungarian: bool = True,
    size_features: Optional[torch.Tensor] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run inference on patches for one spot."""
    with torch.no_grad():
        assignments, confidence = model.assign(
            patches,
            proportions if use_hungarian else None,
            use_hungarian=use_hungarian,
            size_features=size_features,
        )
    return assignments.cpu().numpy(), confidence.cpu().numpy()


def evaluate_single_cell_accuracy(
    model: DirectSoftmaxModel,
    patches_dir: Path,
    proportions_file: Path,
    celltype_cols: List[str],
    xenium_df: pd.DataFrame,
    regions: List[int],
    device: torch.device,
    use_hungarian: bool = True,
    nucleus_features_dir: Optional[Path] = None,
    use_size_features: bool = False,
) -> Dict:
    """
    Evaluate single-cell assignment accuracy.

    Maps Cellpose nuclei to nearest Xenium cells and compares types.
    """
    # Load proportions
    props_df = pd.read_csv(proportions_file)
    props_df["spot_id"] = props_df["spot_id"].astype(str)

    # Load nucleus features (has coordinates)
    all_predictions = []
    all_ground_truth = []
    all_confidence = []

    for region in regions:
        logger.info(f"Processing region {region}...")

        # Load nucleus features for coordinate mapping
        if nucleus_features_dir:
            features_file = nucleus_features_dir / f"region_{region}" / "nucleus_features.csv"
        else:
            features_file = patches_dir.parent / "patches" / f"region_{region}" / "nucleus_features.csv"

        if not features_file.exists():
            logger.warning(f"No nucleus features at {features_file}, skipping region {region}")
            continue

        nucleus_df = pd.read_csv(features_file)
        logger.info(f"Loaded {len(nucleus_df)} nuclei for region {region}")

        # Convert spot_id format: "spot_1032" -> "r0_spot_1032"
        nucleus_df["spot_id_full"] = f"r{region}_" + nucleus_df["spot_id"].astype(str)

        # Filter Xenium cells to region (approximate bounds)
        # Regions are ~2300 microns wide, starting at x=0
        region_width = 2300
        x_min = region * region_width
        x_max = (region + 1) * region_width

        region_xenium = xenium_df[
            (xenium_df['x_centroid'] >= x_min - 100) &
            (xenium_df['x_centroid'] <= x_max + 100) &
            xenium_df['cell_type'].notna()
        ].copy()

        if len(region_xenium) == 0:
            logger.warning(f"No Xenium cells with GT in region {region}")
            continue

        logger.info(f"Found {len(region_xenium)} Xenium cells in region {region}")

        # Build KDTree for Xenium cells
        xen_coords = region_xenium[['x_centroid', 'y_centroid']].values
        xen_tree = cKDTree(xen_coords)

        # Get unique spots in this region
        region_spots = [s for s in props_df["spot_id"].values if s.startswith(f"r{region}_")]
        logger.info(f"Found {len(region_spots)} spots in region {region}")

        for spot_id in tqdm(region_spots, desc=f"Region {region}"):
            # Load patches
            patch_file = patches_dir / f"{spot_id}_patches.npy"
            if not patch_file.exists():
                continue

            patches = np.load(patch_file).astype(np.float32)
            if len(patches) < 3:
                continue

            patches_tensor = torch.from_numpy(patches).to(device)

            # Load size features if needed
            size_features_tensor = None
            if use_size_features:
                size_file = patches_dir / f"{spot_id}_sizes.npy"
                if size_file.exists():
                    sizes = np.load(size_file).astype(np.float32)
                    size_features_tensor = torch.from_numpy(sizes).to(device)
                else:
                    logger.warning(f"Size features not found for {spot_id}, skipping")
                    continue

            # Get proportions
            spot_props = props_df[props_df["spot_id"] == spot_id]
            if len(spot_props) == 0:
                continue

            proportions = torch.tensor(
                spot_props[celltype_cols].values[0].astype(np.float32),
                device=device
            )

            # Run inference
            assignments, confidence = run_inference_on_spot(
                model, patches_tensor, proportions, use_hungarian,
                size_features=size_features_tensor,
            )

            # Get nucleus coordinates for this spot (use full spot_id)
            spot_nuclei = nucleus_df[nucleus_df["spot_id_full"] == spot_id].copy()

            if len(spot_nuclei) == 0:
                continue

            # Ensure we have same number of nuclei as patches
            if len(spot_nuclei) != len(assignments):
                # Take first N nuclei to match patches
                spot_nuclei = spot_nuclei.iloc[:len(assignments)]

            if len(spot_nuclei) == 0:
                continue

            # Transform coordinates (Cellpose pixels -> Xenium microns)
            # Cellpose coordinates are in image pixels, need to convert to tissue microns
            pixel_size = 0.2125  # um/pixel for Xenium morphology images
            cp_x = spot_nuclei["centroid_x"].values * pixel_size + x_min
            cp_y = spot_nuclei["centroid_y"].values * pixel_size

            # Find nearest Xenium cell
            cp_coords = np.column_stack([cp_x, cp_y])
            distances, indices = xen_tree.query(cp_coords, k=1)

            # Filter to close matches
            max_dist = 15.0  # microns
            valid_mask = distances <= max_dist

            for i, (valid, idx, conf) in enumerate(zip(valid_mask, indices, confidence)):
                if valid and i < len(assignments):
                    assignment_idx = assignments[i]
                    # Handle Unknown class (index = n_types when enable_unknown=True)
                    if assignment_idx >= len(celltype_cols):
                        pred_type = "Unknown"
                    else:
                        pred_type = celltype_cols[assignment_idx]
                    gt_type = region_xenium.iloc[idx]["cell_type"]

                    all_predictions.append(pred_type)
                    all_ground_truth.append(gt_type)
                    all_confidence.append(conf)

    if len(all_predictions) == 0:
        logger.error("No valid predictions!")
        return {"error": "No valid predictions"}

    # Compute metrics
    # Map predictions to same label space as ground truth
    unique_gt = sorted(set(all_ground_truth))
    unique_pred = sorted(set(all_predictions))

    logger.info(f"Ground truth types: {unique_gt}")
    logger.info(f"Predicted types: {unique_pred}")

    # Accuracy
    accuracy = accuracy_score(all_ground_truth, all_predictions)

    # Per-class metrics
    report = classification_report(
        all_ground_truth, all_predictions,
        output_dict=True, zero_division=0
    )

    # Confusion matrix
    conf_matrix = confusion_matrix(all_ground_truth, all_predictions, labels=unique_gt)

    results = {
        "n_predictions": len(all_predictions),
        "accuracy": accuracy,
        "mean_confidence": float(np.mean(all_confidence)),
        "classification_report": report,
        "confusion_matrix": conf_matrix.tolist(),
        "labels": unique_gt,
    }

    logger.info(f"\n=== Results ===")
    logger.info(f"Total predictions: {len(all_predictions)}")
    logger.info(f"Accuracy: {accuracy:.4f}")
    logger.info(f"Mean confidence: {np.mean(all_confidence):.4f}")
    logger.info(f"\nPer-class accuracy:")
    for cls in unique_gt:
        if cls in report:
            logger.info(f"  {cls}: precision={report[cls]['precision']:.3f}, recall={report[cls]['recall']:.3f}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Evaluate DirectSoftmaxModel")
    parser.add_argument("--vae-checkpoint", type=str, required=True)
    parser.add_argument("--model-checkpoint", type=str, required=True)
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--proportions-file", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--nucleus-features-dir", type=str, default=None,
                        help="Directory containing region_X/nucleus_features.csv")
    parser.add_argument("--regions", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--no-hungarian", action="store_true")
    args = parser.parse_args()

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    # Load models
    encoder, latent_dim = load_vae_encoder(args.vae_checkpoint, device)
    model, celltype_cols, use_size_features = load_direct_softmax_model(encoder, args.model_checkpoint, device)

    # Load ground truth
    xenium_df = load_xenium_ground_truth()
    logger.info(f"Loaded {len(xenium_df)} Xenium cells")

    # Run evaluation
    results = evaluate_single_cell_accuracy(
        model=model,
        patches_dir=Path(args.patches_dir),
        proportions_file=Path(args.proportions_file),
        celltype_cols=celltype_cols,
        xenium_df=xenium_df,
        regions=args.regions,
        device=device,
        use_hungarian=not args.no_hungarian,
        nucleus_features_dir=Path(args.nucleus_features_dir) if args.nucleus_features_dir else None,
        use_size_features=use_size_features,
    )

    # Save results
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results_file = output_dir / "evaluation_results.json"
    with open(results_file, "w") as f:
        # Convert numpy arrays to lists for JSON
        json_results = {k: v for k, v in results.items()}
        json.dump(json_results, f, indent=2, default=str)

    logger.info(f"Saved results to {results_file}")


if __name__ == "__main__":
    main()
