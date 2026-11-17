#!/usr/bin/env python
"""
Unified preprocessing for Xenium benchmark: StarDist segmentation,
DAPI+boundary patch extraction, SimCLR training, and feature extraction.

Run once before the benchmark array job.

The Xenium morphology data consists of full-tissue OME-TIFF files (one per channel).
This script crops per-region images using REGION_BOUNDS_UM from prepare_patches.py,
then runs StarDist segmentation and patch extraction per region.

Usage:
    python run_preprocess.py \
        --dapi /path/to/ch0000_dapi.ome.tif \
        --boundary /path/to/ch0001_atp1a1_cd45_e-cadherin.ome.tif \
        --input-dir ./data_protein_gt \
        --output-dir ./output/preprocessing
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.spatial import cKDTree

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.morphology.segmentation import StarDistSegmenter
from CITEgeist.model.morphology.simclr import SimCLR

# Import region loading utilities from sibling script
BENCHMARK_SRC = Path(__file__).parent
sys.path.insert(0, str(BENCHMARK_SRC))
from prepare_patches import (
    PIXEL_SIZE_UM,
    REGION_BOUNDS_UM,
    load_dual_channel_image,
    load_ome_tiff_region,
)
from train_simclr import SimCLRAugmentation

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

SPOT_RADIUS_UM = 27.5  # Pseudo-Visium spot radius
N_REGIONS = 5
PATCH_SIZE = 96


def run_stardist_segmentation(dapi_path: Path, input_dir: Path, output_dir: Path) -> None:
    """Stage 1: Run StarDist on DAPI for all regions.

    Loads full-tissue DAPI OME-TIFF, crops per region, runs StarDist,
    assigns nuclei to pseudo-Visium spots via KDTree.
    """
    nuclei_dir = output_dir / "nuclei_counts"
    nuclei_dir.mkdir(parents=True, exist_ok=True)

    segmenter = StarDistSegmenter(modality="dapi")

    for region_id in range(N_REGIONS):
        stage_start = time.time()
        logger.info(f"Segmenting region {region_id}...")

        # Load DAPI region from full-tissue OME-TIFF
        dapi = load_ome_tiff_region(dapi_path, region_id)
        logger.info(f"  DAPI region shape: {dapi.shape}")

        # StarDist expects uint8 or uint16
        if dapi.max() > 0:
            dapi_uint8 = ((dapi / dapi.max()) * 255).astype(np.uint8)
        else:
            dapi_uint8 = dapi.astype(np.uint8)

        # Compute n_tiles to avoid GPU OOM on large DAPI images
        h, w = dapi_uint8.shape[:2]
        n_tiles_y = max(1, h // 512)
        n_tiles_x = max(1, w // 512)
        n_tiles = (n_tiles_y, n_tiles_x)  # 2-tuple for 2D DAPI
        logger.info(f"  n_tiles: {n_tiles} for image {h}x{w}")
        masks, centroids_df = segmenter.segment(dapi_uint8, prob_thresh=0.5, n_tiles=n_tiles)
        logger.info(f"  {len(centroids_df)} nuclei detected")

        # Centroids are in region-local pixel coordinates
        # Convert to tissue-level micron coordinates for spot matching
        x_min_um, _, y_min_um, _ = REGION_BOUNDS_UM[region_id]
        # Account for padding (100 µm)
        padding_um = 100.0
        x_offset_um = x_min_um - padding_um
        y_offset_um = y_min_um - padding_um

        centroids_um = np.column_stack(
            [
                centroids_df["x_pixel"].values * PIXEL_SIZE_UM + x_offset_um,
                centroids_df["y_pixel"].values * PIXEL_SIZE_UM + y_offset_um,
            ]
        )

        # Load spot coordinates
        spots_path = input_dir / "h5ad_objects" / f"Xenium_region_{region_id}_CITE.h5ad"
        protein_adata = sc.read_h5ad(spots_path)
        spot_coords = protein_adata.obsm["spatial"]  # (n_spots, 2) in microns
        spot_barcodes = protein_adata.obs_names.tolist()

        # KDTree assignment: each nucleus → nearest spot within radius
        spot_tree = cKDTree(spot_coords)
        distances, spot_indices = spot_tree.query(centroids_um, k=1)

        in_spot = distances <= SPOT_RADIUS_UM
        centroids_df["spot_index"] = -1
        centroids_df.loc[in_spot, "spot_index"] = spot_indices[in_spot]
        centroids_df["barcode"] = ""
        centroids_df.loc[in_spot, "barcode"] = [spot_barcodes[i] for i in spot_indices[in_spot]]
        centroids_df["distance_um"] = distances

        # Per-spot nuclei counts
        assigned = centroids_df[centroids_df["spot_index"] >= 0]
        counts = assigned.groupby("barcode").size().reindex(spot_barcodes, fill_value=0)
        counts_df = pd.DataFrame({"barcode": spot_barcodes, "nuclei_count": counts.values})
        counts_df.to_csv(nuclei_dir / f"region_{region_id}_nuclei_counts.csv", index=False)

        # Save centroid info for patch extraction
        centroids_df.to_csv(nuclei_dir / f"region_{region_id}_centroids.csv", index=False)

        elapsed = time.time() - stage_start
        logger.info(f"  {in_spot.sum()}/{len(centroids_df)} nuclei assigned to spots ({elapsed:.1f}s)")


def extract_patches(dapi_path: Path, boundary_path: Path, output_dir: Path) -> None:
    """Stage 2: Extract DAPI+boundary patches around each centroid."""
    patches_dir = output_dir / "patches"
    patches_dir.mkdir(parents=True, exist_ok=True)
    features_dir = output_dir / "features"
    features_dir.mkdir(parents=True, exist_ok=True)
    nuclei_dir = output_dir / "nuclei_counts"

    half = PATCH_SIZE // 2

    for region_id in range(N_REGIONS):
        stage_start = time.time()
        logger.info(f"Extracting patches for region {region_id}...")

        # Load 2-channel image (DAPI + boundary), normalized to [0, 1]
        img = load_dual_channel_image(dapi_path, boundary_path, region_id)
        # img shape: (2, H, W)

        # Load centroids (region-local pixel coordinates)
        centroids_df = pd.read_csv(nuclei_dir / f"region_{region_id}_centroids.csv")
        assigned = centroids_df[centroids_df["spot_index"] >= 0].copy()

        H, W = img.shape[1], img.shape[2]
        patches = []
        valid_indices = []

        for idx, row in assigned.iterrows():
            cy = int(row["y_pixel"])
            cx = int(row["x_pixel"])

            if cy - half < 0 or cy + half > H or cx - half < 0 or cx + half > W:
                continue

            patch = img[:, cy - half : cy + half, cx - half : cx + half]
            patches.append(patch)
            valid_indices.append(idx)

        if not patches:
            logger.warning(f"  No valid patches for region {region_id}")
            continue

        patches_arr = np.stack(patches, axis=0).astype(np.float32)
        np.save(patches_dir / f"region_{region_id}_patches.npy", patches_arr)

        # Save nucleus-spot mapping to features dir (consumed by run_benchmark.py)
        valid_df = assigned.loc[valid_indices, ["nucleus_id", "barcode"]].reset_index(drop=True)
        valid_df.to_csv(features_dir / f"region_{region_id}_nucleus_spot_map.csv", index=False)

        elapsed = time.time() - stage_start
        logger.info(f"  Region {region_id}: {len(patches)} patches extracted ({elapsed:.1f}s)")


def train_simclr_backbone(output_dir: Path, epochs: int, batch_size: int, lr: float, temperature: float) -> Path:
    """Stage 3: Train SimCLR on pooled patches from all regions."""
    patches_dir = output_dir / "patches"
    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info(f"Training SimCLR on device: {device}")

    # Pool all patches
    all_patches = []
    for region_id in range(N_REGIONS):
        path = patches_dir / f"region_{region_id}_patches.npy"
        if not path.exists():
            logger.warning(f"  No patches for region {region_id}, skipping")
            continue
        p = np.load(path)
        all_patches.append(p)
        logger.info(f"  Loaded {len(p)} patches from region {region_id}")

    all_patches = np.concatenate(all_patches, axis=0)
    logger.info(f"Total patches for SimCLR training: {len(all_patches)}")

    model = SimCLR(
        in_channels=2,
        img_size=PATCH_SIZE,
        embed_dim=384,
        temperature=temperature,
    )
    model = model.to(device)

    augmentation = SimCLRAugmentation(img_size=PATCH_SIZE)

    # Patches from load_dual_channel_image are already [0, 1]
    patches_tensor = torch.tensor(all_patches, dtype=torch.float32)
    if patches_tensor.max() > 1.0:
        patches_tensor = patches_tensor / 255.0

    dataset = torch.utils.data.TensorDataset(patches_tensor)
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=4,
        drop_last=True,
    )

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    best_loss = float("inf")
    checkpoint_path = output_dir / "simclr_backbone.pt"

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        for (batch,) in loader:
            batch = batch.to(device)
            view1 = augmentation(batch)
            view2 = augmentation(batch)

            loss = model(view1, view2)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        scheduler.step()
        avg_loss = epoch_loss / max(n_batches, 1)

        if (epoch + 1) % 50 == 0:
            logger.info(f"  Epoch {epoch+1}/{epochs}: loss={avg_loss:.4f}")

        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save(model.encoder.state_dict(), checkpoint_path)

    logger.info(f"SimCLR training complete. Best loss: {best_loss:.4f}")
    logger.info(f"Backbone saved to {checkpoint_path}")
    return checkpoint_path


def extract_features(output_dir: Path) -> None:
    """Stage 4: Extract 384-dim features per nucleus using frozen SimCLR encoder."""
    patches_dir = output_dir / "patches"
    features_dir = output_dir / "features"
    features_dir.mkdir(parents=True, exist_ok=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    model = SimCLR(in_channels=2, img_size=PATCH_SIZE, embed_dim=384)
    model.encoder.load_state_dict(torch.load(output_dir / "simclr_backbone.pt", weights_only=True))
    model = model.to(device)
    model.eval()

    for region_id in range(N_REGIONS):
        stage_start = time.time()
        path = patches_dir / f"region_{region_id}_patches.npy"
        if not path.exists():
            logger.warning(f"  No patches for region {region_id}, skipping")
            continue

        logger.info(f"Extracting features for region {region_id}...")
        patches = np.load(path)
        patches_tensor = torch.tensor(patches, dtype=torch.float32)
        if patches_tensor.max() > 1.0:
            patches_tensor = patches_tensor / 255.0

        features_list = []
        batch_size = 256
        with torch.no_grad():
            for i in range(0, len(patches_tensor), batch_size):
                batch = patches_tensor[i : i + batch_size].to(device)
                feats = model.encoder(batch)
                features_list.append(feats.cpu().numpy())

        features = np.concatenate(features_list, axis=0)
        np.save(features_dir / f"region_{region_id}_features.npy", features)

        elapsed = time.time() - stage_start
        logger.info(f"  Region {region_id}: {len(features)} features, shape={features.shape} ({elapsed:.1f}s)")


def main():
    parser = argparse.ArgumentParser(description="Unified preprocessing for Xenium CITEgeist benchmark")
    parser.add_argument(
        "--dapi", type=str, required=True, help="Path to full-tissue DAPI OME-TIFF (ch0000_dapi.ome.tif)"
    )
    parser.add_argument(
        "--boundary", type=str, required=True, help="Path to full-tissue boundary OME-TIFF (ch0001_*.ome.tif)"
    )
    parser.add_argument(
        "--input-dir", type=str, required=True, help="Input directory with h5ad_objects/ (pseudo-Visium data)"
    )
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory for preprocessing results")
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.5)
    parser.add_argument("--skip-segmentation", action="store_true", help="Skip StarDist segmentation (reuse existing)")
    parser.add_argument("--skip-patches", action="store_true", help="Skip patch extraction (reuse existing)")
    parser.add_argument("--skip-training", action="store_true", help="Skip SimCLR training (reuse existing backbone)")
    parser.add_argument("--skip-features", action="store_true", help="Skip feature extraction (reuse existing)")

    args = parser.parse_args()
    dapi_path = Path(args.dapi)
    boundary_path = Path(args.boundary)
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    start = time.time()

    if not args.skip_segmentation:
        logger.info("=" * 60)
        logger.info("STAGE 1: StarDist Segmentation")
        logger.info("=" * 60)
        run_stardist_segmentation(dapi_path, input_dir, output_dir)

    if not args.skip_patches:
        logger.info("=" * 60)
        logger.info("STAGE 2: Patch Extraction")
        logger.info("=" * 60)
        extract_patches(dapi_path, boundary_path, output_dir)

    if not args.skip_training:
        logger.info("=" * 60)
        logger.info("STAGE 3: SimCLR Training")
        logger.info("=" * 60)
        train_simclr_backbone(output_dir, args.epochs, args.batch_size, args.lr, args.temperature)

    if not args.skip_features:
        logger.info("=" * 60)
        logger.info("STAGE 4: Feature Extraction")
        logger.info("=" * 60)
        extract_features(output_dir)

    elapsed = time.time() - start
    logger.info(f"Preprocessing complete in {elapsed / 60:.1f} minutes")


if __name__ == "__main__":
    main()
