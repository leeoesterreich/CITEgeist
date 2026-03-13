#!/usr/bin/env python
"""
Extract ViT-S features from Xenium nucleus patches.

Patches are flat files: patches_dir/r{region}_spot_{spot_id}_patches.npy
Each is (N_nuclei, 2, 96, 96) with DAPI + boundary channels.

Output: per-spot feature files: output_dir/r{region}_spot_{spot_id}.npy
Each is (N_nuclei, 384) ViT-S CLS token features.

Usage:
    python extract_vit_features.py \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/vit_features \
        --device cuda
"""
import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.vit_extractor import ViTFeatureExtractor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)


def convert_2ch_to_3ch(patches: np.ndarray) -> np.ndarray:
    """Convert 2-channel (DAPI+boundary) to 3-channel [DAPI, boundary, boundary].

    Args:
        patches: (N, 2, H, W) array

    Returns:
        (N, 3, H, W) array
    """
    dapi = patches[:, 0:1, :, :]
    boundary = patches[:, 1:2, :, :]
    return np.concatenate([dapi, boundary, boundary], axis=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--vit-model", type=str, default="vit_small_patch16_224")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", type=str, default="cuda")
    args = parser.parse_args()

    patches_dir = Path(args.patches_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize ViT extractor
    logger.info(f"Loading ViT model: {args.vit_model}")
    vit = ViTFeatureExtractor(
        model_name=args.vit_model,
        pretrained=True,
        device=args.device,
    )
    logger.info(f"ViT embed_dim: {vit.embed_dim}")

    # Find all patch files (flat dir: r{region}_spot_{spot_id}_patches.npy)
    patch_files = sorted(patches_dir.glob("r*_spot_*_patches.npy"))
    logger.info(f"Found {len(patch_files)} patch files")

    n_processed = 0
    n_skipped = 0

    for pf in patch_files:
        # Extract spot identifier: r{region}_spot_{spot_id}
        stem = pf.stem.replace("_patches", "")  # e.g., "r0_spot_spot_1032"
        out_file = output_dir / f"{stem}.npy"

        # Skip if already extracted
        if out_file.exists():
            n_skipped += 1
            continue

        patches = np.load(pf)
        if len(patches) == 0:
            np.save(out_file, np.zeros((0, 384), dtype=np.float32))
            n_processed += 1
            continue

        # Convert 2ch -> 3ch if needed
        if patches.ndim == 4 and patches.shape[1] == 2:
            patches = convert_2ch_to_3ch(patches)

        # Resize patches to 224x224 (ViT expects 224x224 input, patches are 96x96)
        patches_t = torch.from_numpy(patches.astype(np.float32))
        patches_t = F.interpolate(patches_t, size=(224, 224), mode="bilinear", align_corners=False)
        patches_resized = patches_t.numpy()

        # Extract features
        features = vit.extract_numpy(patches_resized, batch_size=args.batch_size)
        np.save(out_file, features)
        n_processed += 1

        if n_processed % 200 == 0:
            logger.info(f"  Processed {n_processed} spots ({n_skipped} skipped)")

    logger.info(f"Done. Processed {n_processed}, skipped {n_skipped}")


if __name__ == "__main__":
    main()
