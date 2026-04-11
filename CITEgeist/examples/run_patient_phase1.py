#!/usr/bin/env python
"""Phase 1: StarDist nuclei segmentation for patient H&E samples.

Loads fullres H&E image from SpaceRanger outputs, runs StarDist segmentation,
extracts nucleus centroids, and assigns each nucleus to its nearest Visium spot.

Artifacts saved to --output-dir/<sample_name>/:
    segmentation/nuclei_masks.npy         — integer label mask (crop coords)
    segmentation/nuclei_centroids.csv     — y_pixel, x_pixel, nucleus_id, spot_barcode
    segmentation/nuclei_per_spot.csv      — barcode, n_nuclei
    segmentation/crop_offsets.json        — crop_x1, crop_y1 for back-projection
    .phase1_complete                      — completion marker

Usage:
    python run_patient_phase1.py --sample HCC22-088-P1-S1
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
from scipy.ndimage import center_of_mass
from scipy.spatial import cKDTree


def run_phase1(sample_name, output_dir, data_dir):
    """Run StarDist segmentation for a single patient sample.

    Args:
        sample_name (str): Sample identifier, e.g. 'HCC22-088-P1-S1'.
        output_dir (str or Path): Root directory for phase1 outputs.
        data_dir (str or Path): Root of SpaceRanger processed_files directory.
    """
    output_dir = Path(output_dir)
    data_dir = Path(data_dir)

    sample_out = output_dir / sample_name
    seg_dir = sample_out / "segmentation"
    seg_dir.mkdir(parents=True, exist_ok=True)

    marker_file = sample_out / ".phase1_complete"
    if marker_file.exists():
        logger.info("Phase 1 already complete for %s, skipping", sample_name)
        return

    sample_path = data_dir / sample_name / "outs"
    spatial_dir = sample_path / "spatial"

    # ----------------------------------------------------------------
    # Load scale factors
    # ----------------------------------------------------------------
    with open(spatial_dir / "scalefactors_json.json") as f:
        scalefactors = json.load(f)
    spot_diameter_fullres = scalefactors.get("spot_diameter_fullres", 140.0)
    pixel_size_um = 55.0 / spot_diameter_fullres  # Visium spots are 55 um
    stardist_scale = pixel_size_um / 0.25  # StarDist trained at ~0.25 um/px
    logger.info(
        "Scale factors: spot_diameter_fullres=%.1f, pixel_size=%.3f um/px, " "StarDist scale=%.2f",
        spot_diameter_fullres,
        pixel_size_um,
        stardist_scale,
    )

    # ----------------------------------------------------------------
    # Load tissue positions for bounding-box crop
    # ----------------------------------------------------------------
    tp_path = spatial_dir / "tissue_positions.csv"
    if not tp_path.exists():
        tp_path = spatial_dir / "tissue_positions_list.csv"
    tp = pd.read_csv(tp_path)
    if "in_tissue" in tp.columns:
        in_tissue = tp[tp["in_tissue"] == 1]
    else:
        in_tissue = tp  # older SpaceRanger omits the column for filtered barcodes
    tissue_rows = in_tissue["pxl_row_in_fullres"].values
    tissue_cols = in_tissue["pxl_col_in_fullres"].values

    # ----------------------------------------------------------------
    # Load fullres H&E image — NEVER use hires (too small for 224px patches)
    # ----------------------------------------------------------------
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None  # Fullres images can exceed PIL default limit

    tif_path = spatial_dir / "tissue_fullres_image.tif"
    png_path = spatial_dir / "tissue_fullres_image.png"
    if tif_path.exists():
        fullres_raw = Image.open(tif_path)
    elif png_path.exists():
        fullres_raw = Image.open(png_path)
    else:
        raise FileNotFoundError(
            f"Fullres H&E image not found at {spatial_dir}. " "Expected tissue_fullres_image.tif or .png"
        )
    img_w, img_h = fullres_raw.size
    logger.info("Fullres image loaded: %d x %d px", img_h, img_w)

    # Crop to tissue bounding box with 2-spot padding
    pad = int(spot_diameter_fullres * 2)
    crop_y1 = max(0, int(tissue_rows.min()) - pad)
    crop_y2 = min(img_h, int(tissue_rows.max()) + pad)
    crop_x1 = max(0, int(tissue_cols.min()) - pad)
    crop_x2 = min(img_w, int(tissue_cols.max()) + pad)
    logger.info(
        "Tissue crop: rows [%d:%d], cols [%d:%d]",
        crop_y1,
        crop_y2,
        crop_x1,
        crop_x2,
    )

    fullres_crop = np.array(fullres_raw.crop((crop_x1, crop_y1, crop_x2, crop_y2)))
    del fullres_raw  # Release memory before StarDist
    if fullres_crop.ndim == 2:
        fullres_crop = np.stack([fullres_crop] * 3, axis=-1)
    elif fullres_crop.shape[2] == 4:
        fullres_crop = fullres_crop[:, :, :3]
    logger.info("Cropped tissue image: %s", fullres_crop.shape)

    # Save crop offsets for downstream back-projection
    with open(seg_dir / "crop_offsets.json", "w") as f:
        json.dump({"crop_x1": int(crop_x1), "crop_y1": int(crop_y1)}, f)

    # ----------------------------------------------------------------
    # Run StarDist segmentation
    # ----------------------------------------------------------------
    from model.segmentation import StarDistSegmenter

    h, w = fullres_crop.shape[:2]
    scaled_h = int(h * stardist_scale)
    scaled_w = int(w * stardist_scale)
    n_tiles_y = max(1, scaled_h // 1024)
    n_tiles_x = max(1, scaled_w // 1024)
    logger.info(
        "StarDist n_tiles: (%d, %d, 1) for scaled image %dx%d",
        n_tiles_y,
        n_tiles_x,
        scaled_h,
        scaled_w,
    )

    segmenter = StarDistSegmenter(modality="he")
    masks, _ = segmenter.segment(
        fullres_crop,
        scale=stardist_scale,
        n_tiles=(n_tiles_y, n_tiles_x, 1),
    )

    import tensorflow as tf

    tf.keras.backend.clear_session()
    del segmenter

    n_nuclei = int((np.unique(masks) > 0).sum())
    logger.info("StarDist found %d nuclei", n_nuclei)

    np.save(seg_dir / "nuclei_masks.npy", masks)

    # ----------------------------------------------------------------
    # Extract centroids (in crop coordinates)
    # ----------------------------------------------------------------
    nucleus_ids = np.unique(masks)
    nucleus_ids = nucleus_ids[nucleus_ids > 0]
    centroids = center_of_mass(masks > 0, masks, nucleus_ids)
    centroids_df = pd.DataFrame(centroids, columns=["y_pixel", "x_pixel"])
    centroids_df["nucleus_id"] = nucleus_ids.astype(int)

    # ----------------------------------------------------------------
    # Load spot coordinates and assign nuclei to spots
    # ----------------------------------------------------------------
    import squidpy as sq

    adata = sq.read.visium(
        str(sample_path),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=True,
        gex_only=True,
    )
    raw_coords = adata.obsm["spatial"]  # (x, y) in fullres pixel coords
    finite_mask = np.isfinite(raw_coords).all(axis=1)
    if not finite_mask.all():
        logger.warning("Filtering %d spots with NaN spatial coords", (~finite_mask).sum())
    spatial_coords_fullres = raw_coords[finite_mask]
    barcodes = [b for b, m in zip(adata.obs_names, finite_mask) if m]

    # Shift spot coords to crop space
    spatial_coords = spatial_coords_fullres.copy()
    spatial_coords[:, 0] -= crop_x1  # x
    spatial_coords[:, 1] -= crop_y1  # y

    spot_radius_px = spot_diameter_fullres / 2.0
    centroids_df = _assign_nuclei_to_spots(
        centroids_df,
        spatial_coords,
        barcodes,
        spot_radius=spot_radius_px,
    )

    centroids_df.to_csv(seg_dir / "nuclei_centroids.csv", index=False)

    # Nuclei count per spot
    nps = centroids_df.groupby("spot_barcode").size().reset_index(name="n_nuclei")
    nps.to_csv(seg_dir / "nuclei_per_spot.csv", index=False)
    logger.info(
        "Assigned %d nuclei across %d spots (median %.1f per spot)",
        len(centroids_df),
        len(nps),
        nps["n_nuclei"].median(),
    )

    marker_file.touch()
    logger.info("Phase 1 complete for %s", sample_name)


def _assign_nuclei_to_spots(centroids_df, spatial_coords, barcodes, spot_radius=None):
    """Assign each nucleus to its nearest Visium spot.

    obsm['spatial'] from squidpy is (x, y). center_of_mass returns (y, x),
    so we swap nucleus coords before querying the KDTree.

    Args:
        centroids_df (pd.DataFrame): Has columns y_pixel, x_pixel, nucleus_id.
        spatial_coords (np.ndarray): (n_spots, 2) in (x, y) fullres crop coords.
        barcodes (list): Spot barcodes aligned with spatial_coords.
        spot_radius (float, optional): Discard nuclei farther than this.

    Returns:
        pd.DataFrame: centroids_df with spot_barcode and distance_to_spot columns.
    """
    tree = cKDTree(spatial_coords)
    nucleus_coords = centroids_df[["x_pixel", "y_pixel"]].values  # swap to (x, y)
    distances, indices = tree.query(nucleus_coords)

    centroids_df = centroids_df.copy()
    centroids_df["spot_barcode"] = [barcodes[i] for i in indices]
    centroids_df["distance_to_spot"] = distances

    if spot_radius is not None:
        n_before = len(centroids_df)
        centroids_df = centroids_df[centroids_df["distance_to_spot"] <= spot_radius]
        logger.info(
            "Kept %d/%d nuclei within radius %.1f px",
            len(centroids_df),
            n_before,
            spot_radius,
        )

    return centroids_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Patient pipeline Phase 1: StarDist segmentation")
    parser.add_argument("--sample", required=True, help="Sample name, e.g. HCC22-088-P1-S1")
    parser.add_argument(
        "--output-dir",
        default="output/patient_pipeline/phase1",
        help="Root directory for phase1 outputs",
    )
    parser.add_argument(
        "--data-dir",
        default="/ix1/alee/LO_LAB/General/Lab_Data/" "20250210_CITEGeistPublicData_GEO_Alex/processed_files",
        help="Root of SpaceRanger processed_files directory",
    )
    args = parser.parse_args()
    run_phase1(args.sample, args.output_dir, args.data_dir)
