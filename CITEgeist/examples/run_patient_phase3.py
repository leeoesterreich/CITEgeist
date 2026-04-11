#!/usr/bin/env python
"""Phase 3: ViT-Small feature extraction from H&E nucleus patches.

Loads the fullres H&E image and nucleus centroids from Phase 1, extracts
224x224 px patches centred on each nucleus, encodes them with an ImageNet
ViT-Small (vit_small_patch16_224), and saves 384-dim embeddings.

Artifacts saved to --output-dir/<sample_name>/:
    vit_features.npy   — (N, 384) float32 embeddings
    nucleus_ids.npy    — (N,) int nucleus IDs aligned with vit_features
    .phase3_complete

Usage:
    python run_patient_phase3.py --sample HCC22-088-P1-S1
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("TF_FORCE_GPU_ALLOW_GROWTH", "true")

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import pandas as pd
import torch
from model.unified_config import PATCH_SIZE, VIT_MODEL
from model.vit_extractor import ViTFeatureExtractor


def run_phase3(sample_name, output_dir, seg_dir, data_dir):
    """Extract ViT-Small features from H&E nucleus patches.

    Args:
        sample_name (str): Sample identifier, e.g. 'HCC22-088-P1-S1'.
        output_dir (str or Path): Root directory for phase3 feature outputs.
        seg_dir (str or Path): Root directory of phase1 segmentation outputs.
        data_dir (str or Path): Root of SpaceRanger processed_files directory.
    """
    output_dir = Path(output_dir)
    seg_dir = Path(seg_dir)
    data_dir = Path(data_dir)

    sample_out = output_dir / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    marker_file = sample_out / ".phase3_complete"
    if marker_file.exists():
        logger.info("Phase 3 already complete for %s, skipping", sample_name)
        return

    # ----------------------------------------------------------------
    # Validate Phase 1 completion
    # ----------------------------------------------------------------
    phase1_marker = seg_dir / sample_name / ".phase1_complete"
    if not phase1_marker.exists():
        raise RuntimeError(f"Phase 1 not complete for {sample_name}. " f"Expected marker: {phase1_marker}")

    seg_sample_dir = seg_dir / sample_name / "segmentation"

    # ----------------------------------------------------------------
    # Load nucleus centroids (in crop coordinates)
    # ----------------------------------------------------------------
    centroids_df = pd.read_csv(seg_sample_dir / "nuclei_centroids.csv")
    logger.info("Loaded %d nucleus centroids for %s", len(centroids_df), sample_name)

    # ----------------------------------------------------------------
    # Load crop offsets to map centroid coords back to fullres image
    # ----------------------------------------------------------------
    offsets_path = seg_sample_dir / "crop_offsets.json"
    if offsets_path.exists():
        with open(offsets_path) as f:
            offsets = json.load(f)
        crop_x1 = offsets["crop_x1"]
        crop_y1 = offsets["crop_y1"]
    else:
        logger.warning("crop_offsets.json not found; assuming centroids are in fullres coords")
        crop_x1 = 0
        crop_y1 = 0

    # ----------------------------------------------------------------
    # Load fullres H&E image (NEVER use hires — too small for 224px patches)
    # ----------------------------------------------------------------
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None

    sample_path = data_dir / sample_name / "outs" / "spatial"
    tif_path = sample_path / "tissue_fullres_image.tif"
    png_path = sample_path / "tissue_fullres_image.png"
    if tif_path.exists():
        fullres_raw = Image.open(tif_path)
    elif png_path.exists():
        fullres_raw = Image.open(png_path)
    else:
        raise FileNotFoundError(
            f"Fullres H&E image not found at {sample_path}. " "Expected tissue_fullres_image.tif or .png"
        )
    fullres_img = np.array(fullres_raw)
    del fullres_raw
    if fullres_img.ndim == 2:
        fullres_img = np.stack([fullres_img] * 3, axis=-1)
    elif fullres_img.shape[2] == 4:
        fullres_img = fullres_img[:, :, :3]
    logger.info("Fullres image shape: %s", fullres_img.shape)

    # ----------------------------------------------------------------
    # Extract 224px patches centred on each nucleus (fullres coords)
    # ----------------------------------------------------------------
    patches, valid_ids = _extract_patches(
        fullres_img,
        centroids_df,
        crop_x1,
        crop_y1,
        PATCH_SIZE,
    )

    if len(patches) == 0:
        logger.error("No valid patches extracted for %s; aborting Phase 3", sample_name)
        return

    del fullres_img  # Free memory before GPU encoding

    # ----------------------------------------------------------------
    # Encode with ImageNet ViT-Small → 384-dim embeddings
    # ----------------------------------------------------------------
    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info("Encoding %d patches with %s on %s", len(patches), VIT_MODEL, device)

    extractor = ViTFeatureExtractor(model_name=VIT_MODEL, device=device)
    patches_chw = np.transpose(patches, (0, 3, 1, 2)).astype(np.float32)
    if patches_chw.max() > 1.0:
        patches_chw = patches_chw / 255.0

    features = extractor.extract_numpy(patches_chw, batch_size=64)
    logger.info("Extracted features: %s", features.shape)

    # ----------------------------------------------------------------
    # Save outputs
    # ----------------------------------------------------------------
    np.save(sample_out / "vit_features.npy", features)
    np.save(sample_out / "nucleus_ids.npy", np.array(valid_ids, dtype=np.int64))

    marker_file.touch()
    logger.info("Phase 3 complete for %s: %d nuclei encoded", sample_name, len(features))


def _extract_patches(image, centroids_df, crop_x1, crop_y1, patch_size=224):
    """Extract square patches centered on nucleus centroids.

    Centroids are in crop-space coords; we add the crop offsets to obtain
    fullres pixel positions, then crop directly from the fullres image.

    Args:
        image (np.ndarray): Fullres H&E image (H, W, 3).
        centroids_df (pd.DataFrame): Has y_pixel, x_pixel, nucleus_id columns
            in crop-coordinate space.
        crop_x1 (int): X offset of crop origin in fullres image.
        crop_y1 (int): Y offset of crop origin in fullres image.
        patch_size (int): Side length of patch in pixels.

    Returns:
        Tuple[np.ndarray, list]: patches array (N, patch_size, patch_size, 3)
            and list of corresponding nucleus IDs.
    """
    img_h, img_w = image.shape[:2]
    half = patch_size // 2
    patches = []
    valid_ids = []

    for _, row in centroids_df.iterrows():
        # Convert crop coords → fullres coords
        cy = int(row["y_pixel"]) + crop_y1
        cx = int(row["x_pixel"]) + crop_x1

        y1, y2 = cy - half, cy + half
        x1, x2 = cx - half, cx + half

        if y1 < 0 or x1 < 0 or y2 > img_h or x2 > img_w:
            continue

        patch = image[y1:y2, x1:x2]
        patches.append(patch)
        valid_ids.append(int(row["nucleus_id"]))

    patches_arr = np.array(patches) if patches else np.empty((0, patch_size, patch_size, 3))
    logger.info(
        "Extracted %d patches (%d skipped, out of bounds)",
        len(patches_arr),
        len(centroids_df) - len(patches_arr),
    )
    return patches_arr, valid_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Patient pipeline Phase 3: ViT feature extraction")
    parser.add_argument("--sample", required=True, help="Sample name, e.g. HCC22-088-P1-S1")
    parser.add_argument(
        "--output-dir",
        default="output/patient_pipeline/phase3",
        help="Root directory for phase3 feature outputs",
    )
    parser.add_argument(
        "--seg-dir",
        default="output/patient_pipeline/phase1",
        help="Root directory of phase1 segmentation outputs",
    )
    parser.add_argument(
        "--data-dir",
        default="/ix1/alee/LO_LAB/General/Lab_Data/" "20250210_CITEGeistPublicData_GEO_Alex/processed_files",
        help="Root of SpaceRanger processed_files directory",
    )
    args = parser.parse_args()
    run_phase3(args.sample, args.output_dir, args.seg_dir, args.data_dir)
