#!/usr/bin/env python
"""Step 2+3: Cellpose patch extraction + ViT-S feature extraction.

Usage:
    python run_unified_step2_features.py --sample HCC22-088-P1-S1 --modality he
"""
import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.ndimage import center_of_mass
from scipy.spatial import cKDTree

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from model.unified_config import OUTPUT_BASE, VIT_MODEL, PATCH_SIZE, DATA_DIR
from model.vit_extractor import ViTFeatureExtractor


def _get_cellpose_model(model_type, use_gpu):
    """Cellpose v3/v4 compatibility wrapper."""
    try:
        from cellpose import models
        return models.Cellpose(model_type=model_type, gpu=use_gpu)
    except AttributeError:
        from cellpose.models import CellposeModel
        return CellposeModel(model_type=model_type, gpu=use_gpu)


def load_image_and_segment(sample_name, modality):
    """Load image + run/reuse Cellpose segmentation."""
    cellpose_dir = OUTPUT_BASE / sample_name / "cellpose"
    cellpose_dir.mkdir(parents=True, exist_ok=True)

    if modality == "he":
        import squidpy as sq
        import json as _json
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None  # Fullres images can be large

        sample_path = DATA_DIR / sample_name / "outs"
        spatial_dir = sample_path / "spatial"

        # Load fullres image for patch extraction
        fullres_path = spatial_dir / "tissue_fullres_image.tif"
        if not fullres_path.exists():
            fullres_path = spatial_dir / "tissue_fullres_image.png"
        logger.info(f"Loading fullres image from {fullres_path}")
        fullres_img = np.array(Image.open(fullres_path))
        if fullres_img.ndim == 2:
            fullres_img = np.stack([fullres_img] * 3, axis=-1)
        elif fullres_img.shape[2] == 4:
            fullres_img = fullres_img[:, :, :3]  # Drop alpha
        logger.info(f"Fullres image shape: {fullres_img.shape}")

        # Load scale factors
        with open(spatial_dir / "scalefactors_json.json") as f:
            scalefactors = _json.load(f)
        hires_scale = scalefactors["tissue_hires_scalef"]

        # Load spot coordinates — need load_images=True for obsm['spatial']
        adata = sq.read.visium(
            str(sample_path), counts_file="filtered_feature_bc_matrix.h5",
            load_images=True, gex_only=True,
        )
        # obsm['spatial'] is in fullres pixel coordinates
        raw_coords = adata.obsm["spatial"]
        finite_mask = np.isfinite(raw_coords).all(axis=1)
        if not finite_mask.all():
            n_nan = (~finite_mask).sum()
            logger.warning(f"Filtering {n_nan} spots with NaN spatial coords")
        spatial_coords = raw_coords[finite_mask]  # Fullres pixel coords
        barcodes = [b for b, m in zip(adata.obs_names, finite_mask) if m]

        # Cellpose on fullres (or reuse cached fullres masks)
        module3_dir = OUTPUT_BASE / sample_name / "module3"
        cached_fullres_masks = list(module3_dir.glob("*_cellpose_masks_fullres.npy"))
        cached_hires_masks = list(module3_dir.glob("*_cellpose_masks_hires.npy"))

        if cached_fullres_masks:
            logger.info(f"Reusing cached fullres Cellpose masks from {cached_fullres_masks[0]}")
            masks = np.load(cached_fullres_masks[0])
        else:
            logger.info("Running Cellpose on fullres H&E image")
            cp_model = _get_cellpose_model("cyto2", torch.cuda.is_available())
            masks, _, _, _ = cp_model.eval(fullres_img, channels=[0, 0], diameter=None)
            # Cache for future use
            np.save(module3_dir / f"{sample_name}_cellpose_masks_fullres.npy", masks)
            logger.info(f"Cached fullres masks: {len(np.unique(masks)) - 1} nuclei")

        image = fullres_img
    else:
        raise NotImplementedError("DAPI modality deferred to Xenium follow-up")

    # Extract centroids
    nucleus_ids = np.unique(masks)
    nucleus_ids = nucleus_ids[nucleus_ids > 0]
    centroids = center_of_mass(masks > 0, masks, nucleus_ids)
    centroids_df = pd.DataFrame(centroids, columns=["y_pixel", "x_pixel"])
    centroids_df["nucleus_id"] = nucleus_ids.astype(int)

    np.save(cellpose_dir / "nuclei_masks.npy", masks)
    centroids_df.to_csv(cellpose_dir / "nuclei_centroids.csv", index=False)

    # Spot radius in fullres pixels (half the spot diameter)
    spot_radius_px = scalefactors.get("spot_diameter_fullres", 140) / 2.0

    return masks, centroids_df, image, spatial_coords, barcodes, spot_radius_px


def assign_nuclei_to_spots(centroids_df, spatial_coords, barcodes, spot_radius=None):
    """Assign each nucleus to nearest Visium spot.

    Note: obsm['spatial'] from squidpy is (x, y) in fullres pixel coords.
    center_of_mass returns (row, col) = (y, x). We match by swapping nucleus
    coords to (x, y) before building the KDTree.
    """
    tree = cKDTree(spatial_coords)  # (x, y) from obsm['spatial']
    # Swap nucleus (y, x) → (x, y) to match spatial_coords convention
    nucleus_coords = centroids_df[["x_pixel", "y_pixel"]].values
    distances, indices = tree.query(nucleus_coords)

    centroids_df = centroids_df.copy()
    centroids_df["spot_barcode"] = [barcodes[i] for i in indices]
    centroids_df["distance_to_spot"] = distances

    if spot_radius is not None:
        n_before = len(centroids_df)
        centroids_df = centroids_df[centroids_df["distance_to_spot"] <= spot_radius]
        logger.info(f"Kept {len(centroids_df)}/{n_before} nuclei within radius {spot_radius}")

    return centroids_df


def extract_patches(image, centroids_df, modality, patch_size=224):
    """Extract patches centered on each nucleus centroid."""
    h, w = image.shape[:2]
    half = patch_size // 2
    patches = []
    valid_ids = []

    for _, row in centroids_df.iterrows():
        cy, cx = int(row["y_pixel"]), int(row["x_pixel"])
        y1, y2 = cy - half, cy + half
        x1, x2 = cx - half, cx + half

        if y1 < 0 or x1 < 0 or y2 > h or x2 > w:
            continue

        patch = image[y1:y2, x1:x2]

        if modality == "dapi" and patch.shape[2] == 2:
            pad = np.zeros((patch_size, patch_size, 1), dtype=patch.dtype)
            patch = np.concatenate([patch, pad], axis=2)

        patches.append(patch)
        valid_ids.append(int(row["nucleus_id"]))

    patches_arr = np.array(patches)
    logger.info(f"Extracted {len(patches_arr)} patches "
                f"({len(centroids_df) - len(patches_arr)} skipped, out of bounds)")
    return patches_arr, valid_ids


def extract_vit_features(patches, device="cuda"):
    """Extract ViT-S features from patches."""
    extractor = ViTFeatureExtractor(model_name=VIT_MODEL, device=device)

    patches_chw = np.transpose(patches, (0, 3, 1, 2)).astype(np.float32)
    if patches_chw.max() > 1.0:
        patches_chw = patches_chw / 255.0

    features = extractor.extract_numpy(patches_chw, batch_size=64)
    logger.info(f"Extracted features: {features.shape}")
    return features


def run_step2(sample_name, modality="he"):
    step1_marker = OUTPUT_BASE / sample_name / ".step1_complete"
    if not step1_marker.exists():
        logger.error(f"Step 1 not complete for {sample_name}")
        return

    step2_marker = OUTPUT_BASE / sample_name / ".step2_complete"
    if step2_marker.exists():
        logger.info(f"Step 2 already complete for {sample_name}, skipping")
        return

    masks, centroids_df, image, spatial_coords, barcodes, spot_radius_px = load_image_and_segment(
        sample_name, modality,
    )
    centroids_df = assign_nuclei_to_spots(
        centroids_df, spatial_coords, barcodes, spot_radius=spot_radius_px,
    )
    patches, valid_ids = extract_patches(image, centroids_df, modality, PATCH_SIZE)

    if len(patches) == 0:
        logger.error(f"No valid patches for {sample_name}")
        return

    device = "cuda" if torch.cuda.is_available() else "cpu"
    features = extract_vit_features(patches, device)

    feat_dir = OUTPUT_BASE / sample_name / "features"
    feat_dir.mkdir(parents=True, exist_ok=True)
    np.save(feat_dir / "vit_features.npy", features)
    np.save(feat_dir / "nucleus_ids.npy", np.array(valid_ids))

    cellpose_dir = OUTPUT_BASE / sample_name / "cellpose"
    valid_centroids = centroids_df[centroids_df["nucleus_id"].isin(valid_ids)]
    valid_centroids.to_csv(cellpose_dir / "nuclei_centroids.csv", index=False)

    if "spot_barcode" in valid_centroids.columns:
        nps = valid_centroids.groupby("spot_barcode").size().reset_index(name="n_nuclei")
        nps.to_csv(cellpose_dir / "nuclei_per_spot.csv", index=False)

    step2_marker.touch()
    logger.info(f"Step 2 complete for {sample_name}: {len(features)} nuclei")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 2+3: Features")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--modality", default="he", choices=["he", "dapi"])
    args = parser.parse_args()
    run_step2(args.sample, args.modality)
