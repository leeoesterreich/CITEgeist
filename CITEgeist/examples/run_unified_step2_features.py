#!/usr/bin/env python
"""Step 2+3: StarDist segmentation + patch extraction + ViT-S feature extraction.

Usage:
    python run_unified_step2_features.py --sample HCC22-088-P1-S1 --modality he
"""
import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("TF_FORCE_GPU_ALLOW_GROWTH", "true")

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
from model.segmentation import StarDistSegmenter


def load_image_and_segment(sample_name, modality):
    """Load image + run/reuse StarDist segmentation."""
    seg_dir = OUTPUT_BASE / sample_name / "segmentation"
    seg_dir.mkdir(parents=True, exist_ok=True)

    if modality == "he":
        import squidpy as sq
        import json as _json
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None  # Fullres images can be large

        sample_path = DATA_DIR / sample_name / "outs"
        spatial_dir = sample_path / "spatial"

        # Load scale factors
        with open(spatial_dir / "scalefactors_json.json") as f:
            scalefactors = _json.load(f)
        spot_diameter = scalefactors.get("spot_diameter_fullres", 140)

        # Load tissue positions to get tissue bounding box (fullres coords)
        tp = pd.read_csv(spatial_dir / "tissue_positions.csv")
        in_tissue = tp[tp["in_tissue"] == 1]
        tissue_rows = in_tissue["pxl_row_in_fullres"].values
        tissue_cols = in_tissue["pxl_col_in_fullres"].values

        # Crop bounding box with padding (1 spot diameter on each side)
        pad = int(spot_diameter * 2)
        fullres_raw = Image.open(
            spatial_dir / "tissue_fullres_image.tif"
            if (spatial_dir / "tissue_fullres_image.tif").exists()
            else spatial_dir / "tissue_fullres_image.png"
        )
        img_w, img_h = fullres_raw.size
        crop_y1 = max(0, int(tissue_rows.min()) - pad)
        crop_y2 = min(img_h, int(tissue_rows.max()) + pad)
        crop_x1 = max(0, int(tissue_cols.min()) - pad)
        crop_x2 = min(img_w, int(tissue_cols.max()) + pad)
        logger.info(f"Tissue crop: [{crop_y1}:{crop_y2}, {crop_x1}:{crop_x2}] "
                     f"(full image: {img_h}x{img_w})")

        # Crop fullres to tissue region only
        fullres_crop = np.array(fullres_raw.crop((crop_x1, crop_y1, crop_x2, crop_y2)))
        del fullres_raw  # Free memory
        if fullres_crop.ndim == 2:
            fullres_crop = np.stack([fullres_crop] * 3, axis=-1)
        elif fullres_crop.shape[2] == 4:
            fullres_crop = fullres_crop[:, :, :3]
        logger.info(f"Cropped tissue image: {fullres_crop.shape}")

        # Load spot coordinates — need load_images=True for obsm['spatial']
        adata = sq.read.visium(
            str(sample_path), counts_file="filtered_feature_bc_matrix.h5",
            load_images=True, gex_only=True,
        )
        # obsm['spatial'] is (col, row) = (x, y) in fullres pixel coords
        raw_coords = adata.obsm["spatial"]
        finite_mask = np.isfinite(raw_coords).all(axis=1)
        if not finite_mask.all():
            n_nan = (~finite_mask).sum()
            logger.warning(f"Filtering {n_nan} spots with NaN spatial coords")
        # Shift coordinates to crop space
        spatial_coords_fullres = raw_coords[finite_mask]
        spatial_coords = spatial_coords_fullres.copy()
        spatial_coords[:, 0] -= crop_x1  # x offset
        spatial_coords[:, 1] -= crop_y1  # y offset
        barcodes = [b for b, m in zip(adata.obs_names, finite_mask) if m]

        # StarDist on cropped tissue image
        module3_dir = OUTPUT_BASE / sample_name / "module3"
        cached_crop_masks = list(module3_dir.glob("*_stardist_masks_tissue.npy"))

        if cached_crop_masks:
            logger.info(f"Reusing cached tissue StarDist masks from {cached_crop_masks[0]}")
            masks = np.load(cached_crop_masks[0])
        else:
            logger.info(f"Running StarDist on tissue crop ({fullres_crop.shape})")
            segmenter = StarDistSegmenter(modality="he")
            masks = segmenter.segment(fullres_crop)
            np.save(module3_dir / f"{sample_name}_stardist_masks_tissue.npy", masks)
            n_nuclei = len(np.unique(masks)) - 1
            logger.info(f"StarDist found {n_nuclei} nuclei in tissue crop")
            # Free TF memory
            import tensorflow as tf
            tf.keras.backend.clear_session()
            del segmenter

        image = fullres_crop
    else:
        raise NotImplementedError("DAPI modality deferred to Xenium follow-up")

    # Extract centroids
    nucleus_ids = np.unique(masks)
    nucleus_ids = nucleus_ids[nucleus_ids > 0]
    centroids = center_of_mass(masks > 0, masks, nucleus_ids)
    centroids_df = pd.DataFrame(centroids, columns=["y_pixel", "x_pixel"])
    centroids_df["nucleus_id"] = nucleus_ids.astype(int)

    np.save(seg_dir / "nuclei_masks.npy", masks)
    centroids_df.to_csv(seg_dir / "nuclei_centroids.csv", index=False)

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

    seg_dir = OUTPUT_BASE / sample_name / "segmentation"
    valid_centroids = centroids_df[centroids_df["nucleus_id"].isin(valid_ids)]
    valid_centroids.to_csv(seg_dir / "nuclei_centroids.csv", index=False)

    if "spot_barcode" in valid_centroids.columns:
        nps = valid_centroids.groupby("spot_barcode").size().reset_index(name="n_nuclei")
        nps.to_csv(seg_dir / "nuclei_per_spot.csv", index=False)

    step2_marker.touch()
    logger.info(f"Step 2 complete for {sample_name}: {len(features)} nuclei")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 2+3: Features")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--modality", default="he", choices=["he", "dapi"])
    args = parser.parse_args()
    run_step2(args.sample, args.modality)
