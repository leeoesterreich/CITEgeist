#!/usr/bin/env python
"""Diagnostic: compare StarDist nuclei counts with/without scale+prob_thresh tuning.

Runs StarDist on Region 1 DAPI image under 4 conditions:
  1. baseline (no scale, no prob_thresh) — current pipeline default
  2. scale=0.85 only
  3. prob_thresh=0.6 only
  4. scale=0.85 + prob_thresh=0.6

For each condition, reports:
  - Total nuclei detected
  - Nuclei assigned to spots (within 55um radius)
  - Mean/median nuclei per spot
  - Ratio vs Xenium ground truth (22,910 cells / 16.2 per spot)

Also computes per-nucleus area distribution to assess debris/fragment prevalence.

Usage:
    python diagnose_stardist_overcounting.py [--region 1]
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None  # Large Xenium images exceed default limit
from scipy.spatial import cKDTree

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Paths
IMAGE_DIR = Path(
    str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "scResolve" / "images" / "morphology_hires")
)
DATA_DIR = Path(
    str(_cg.REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt")
)
GT_SUMMARY = DATA_DIR / "dataset_summary.json"

# Xenium pixel size and StarDist training resolution
XENIUM_PIXEL_SIZE_UM = 0.2125
STARDIST_TRAINING_PIXEL_SIZE = 0.25
SPOT_DIAMETER_UM = 55.0

# Conditions to test
CONDITIONS = {
    "baseline": {},
    "scale_only": {"scale": XENIUM_PIXEL_SIZE_UM / STARDIST_TRAINING_PIXEL_SIZE},
    "prob_only": {"prob_thresh": 0.6},
    "scale_and_prob": {
        "scale": XENIUM_PIXEL_SIZE_UM / STARDIST_TRAINING_PIXEL_SIZE,
        "prob_thresh": 0.6,
    },
    "scale_and_prob_strict": {
        "scale": XENIUM_PIXEL_SIZE_UM / STARDIST_TRAINING_PIXEL_SIZE,
        "prob_thresh": 0.7,
    },
}


def load_xenium_gt(region_id: int) -> dict:
    """Load Xenium ground truth cell counts for a region."""
    with open(GT_SUMMARY) as f:
        summary = json.load(f)
    for r in summary["regions"]:
        if r["region_id"] == region_id:
            return r
    raise ValueError(f"Region {region_id} not found in GT summary")


def load_image_and_spots(region_id: int):
    """Load DAPI image and pseudo-Visium spot coordinates."""
    import scanpy as sc

    region_name = f"Xenium_region_{region_id}"

    # Load morphology image
    img_path = IMAGE_DIR / region_name / "morphology.png"
    rgb = np.array(Image.open(img_path).convert("RGB"))
    logger.info("Image: %s, shape=%s", img_path.name, rgb.shape)

    # Load coordinate info
    coord_path = IMAGE_DIR / region_name / "coord_info.json"
    with open(coord_path) as f:
        coord_info = json.load(f)
    pixel_size = float(coord_info["pixel_size"])

    # Load spot coordinates
    cite_path = DATA_DIR / "h5ad_objects" / f"{region_name}_CITE.h5ad"
    cite_adata = sc.read_h5ad(cite_path)
    spot_coords_micron = np.asarray(cite_adata.obsm["spatial"], dtype=np.float64)

    # Convert to pixel coords
    offset_x = float(coord_info["micron_bounds"]["x_min"])
    offset_y = float(coord_info["micron_bounds"]["y_min"])
    spot_coords_pixel = np.zeros_like(spot_coords_micron)
    spot_coords_pixel[:, 0] = (spot_coords_micron[:, 0] - offset_x) / pixel_size
    spot_coords_pixel[:, 1] = (spot_coords_micron[:, 1] - offset_y) / pixel_size

    spot_radius_px = (SPOT_DIAMETER_UM / pixel_size) / 2.0
    spot_names = cite_adata.obs_names

    return rgb, spot_coords_pixel, spot_radius_px, spot_names, pixel_size


def run_stardist_condition(image_rgb, condition_name, stardist_kwargs, model, img_norm, n_tiles):
    """Run StarDist with specific parameters and return centroids + areas.

    Uses StarDist's details['points'] directly instead of recomputing
    center_of_mass from the mask (which OOMs on 298M-pixel images with 166K labels).
    """
    logger.info("--- Condition: %s (kwargs=%s) ---", condition_name, stardist_kwargs)

    t0 = time.time()
    masks, details = model.predict_instances(img_norm, n_tiles=n_tiles, **stardist_kwargs)
    elapsed = time.time() - t0

    # Use StarDist's native centroid output (row, col) → (x, y)
    points = details["points"]  # (N, 2) in (row, col) order
    n_nuclei = len(points)
    logger.info("  Detected %d nuclei in %.1fs", n_nuclei, elapsed)

    if n_nuclei > 0:
        centroids_xy = np.column_stack([points[:, 1], points[:, 0]]).astype(np.float64)
        # Efficient area computation via bincount (no float32 copy of mask)
        masks_int = np.asarray(masks, dtype=np.int32)
        areas = np.bincount(masks_int.ravel())
        areas = areas[1:]  # skip background (label 0)
        # Probabilities for reporting
        probs = details.get("prob", np.ones(n_nuclei))
    else:
        centroids_xy = np.zeros((0, 2), dtype=np.float64)
        areas = np.array([], dtype=np.int64)
        probs = np.array([], dtype=np.float64)

    del masks  # Free ~1.2GB immediately
    return centroids_xy, areas, probs, elapsed


def assign_to_spots(centroids_xy, spot_coords_pixel, spot_radius_px):
    """Count nuclei assigned to each spot."""
    if centroids_xy.shape[0] == 0:
        return 0, np.zeros(spot_coords_pixel.shape[0])

    tree = cKDTree(spot_coords_pixel)
    dists, _ = tree.query(centroids_xy, k=1, workers=-1)
    assigned = dists <= spot_radius_px
    n_assigned = int(assigned.sum())

    # Per-spot counts
    _, idxs = tree.query(centroids_xy[assigned], k=1, workers=-1)
    counts = np.zeros(spot_coords_pixel.shape[0], dtype=np.int64)
    np.add.at(counts, idxs, 1)

    return n_assigned, counts


def compute_area_stats(areas, pixel_size):
    """Compute area distribution statistics in um^2."""
    if len(areas) == 0:
        return {}
    areas_um2 = areas * (pixel_size**2)
    return {
        "mean_um2": float(np.mean(areas_um2)),
        "median_um2": float(np.median(areas_um2)),
        "std_um2": float(np.std(areas_um2)),
        "p5_um2": float(np.percentile(areas_um2, 5)),
        "p25_um2": float(np.percentile(areas_um2, 25)),
        "p75_um2": float(np.percentile(areas_um2, 75)),
        "p95_um2": float(np.percentile(areas_um2, 95)),
        "n_below_20um2": int(np.sum(areas_um2 < 20)),
        "n_below_30um2": int(np.sum(areas_um2 < 30)),
        "pct_below_30um2": float(100 * np.mean(areas_um2 < 30)),
    }


def main():
    parser = argparse.ArgumentParser(description="Diagnose StarDist overcounting")
    parser.add_argument("--region", type=int, default=1, help="Region ID (default: 1)")
    args = parser.parse_args()

    region_id = args.region
    logger.info("=" * 70)
    logger.info("STARDIST OVERCOUNTING DIAGNOSTIC — Region %d", region_id)
    logger.info("=" * 70)

    # Load ground truth
    gt = load_xenium_gt(region_id)
    gt_n_cells = gt["n_cells"]
    gt_cells_per_spot = gt["mean_cells_per_spot"]
    gt_n_spots = gt["n_spots"]
    logger.info("Xenium GT: %d cells, %d spots, %.1f cells/spot", gt_n_cells, gt_n_spots, gt_cells_per_spot)

    # Load image and spots
    image_rgb, spot_coords_pixel, spot_radius_px, spot_names, pixel_size = load_image_and_spots(region_id)
    n_spots = len(spot_names)
    logger.info("Spots: %d, radius=%.1f px (%.1f um)", n_spots, spot_radius_px, spot_radius_px * pixel_size)

    # Load StarDist model once (reused across conditions)
    from stardist.models import StarDist2D
    from csbdeep.utils import normalize as csb_normalize

    logger.info("Loading StarDist model...")
    model = StarDist2D.from_pretrained("2D_versatile_fluo")

    # Pre-normalize image once (shared across conditions)
    img = image_rgb[..., 0] if image_rgb.ndim == 3 else image_rgb
    del image_rgb  # Free RGB copy (~900MB)
    img_norm = csb_normalize(img, 1, 99.8)
    del img

    h, w = img_norm.shape[:2]
    n_tiles = (max(1, h // 512), max(1, w // 512))
    logger.info("Image normalized. n_tiles=%s for %dx%d", n_tiles, h, w)

    # Run each condition
    results = {}
    for cond_name, kwargs in CONDITIONS.items():
        centroids_xy, areas, probs, elapsed = run_stardist_condition(
            image_rgb=None,  # unused now
            condition_name=cond_name,
            stardist_kwargs=kwargs,
            model=model,
            img_norm=img_norm,
            n_tiles=n_tiles,
        )
        n_total = len(centroids_xy)
        n_assigned, spot_counts = assign_to_spots(centroids_xy, spot_coords_pixel, spot_radius_px)
        area_stats = compute_area_stats(areas, pixel_size)

        mean_per_spot = float(np.mean(spot_counts)) if n_spots > 0 else 0
        median_per_spot = float(np.median(spot_counts)) if n_spots > 0 else 0

        # Probability distribution stats
        prob_stats = {}
        if len(probs) > 0:
            prob_stats = {
                "mean": float(np.mean(probs)),
                "median": float(np.median(probs)),
                "p10": float(np.percentile(probs, 10)),
                "p25": float(np.percentile(probs, 25)),
                "n_below_0.5": int(np.sum(probs < 0.5)),
                "n_below_0.6": int(np.sum(probs < 0.6)),
                "n_below_0.7": int(np.sum(probs < 0.7)),
            }

        result = {
            "condition": cond_name,
            "stardist_kwargs": {k: float(v) for k, v in kwargs.items()},
            "n_total_nuclei": n_total,
            "n_assigned_to_spots": n_assigned,
            "pct_assigned": float(100 * n_assigned / n_total) if n_total > 0 else 0,
            "mean_nuclei_per_spot": mean_per_spot,
            "median_nuclei_per_spot": median_per_spot,
            "ratio_vs_xenium_total": float(n_assigned / gt_n_cells) if gt_n_cells > 0 else 0,
            "ratio_vs_xenium_per_spot": float(mean_per_spot / gt_cells_per_spot) if gt_cells_per_spot > 0 else 0,
            "area_stats": area_stats,
            "prob_stats": prob_stats,
            "elapsed_sec": elapsed,
        }
        results[cond_name] = result

        logger.info(
            "  RESULT: %d total, %d assigned (%.0f%%), %.1f/spot, " "ratio=%.2fx vs Xenium",
            n_total,
            n_assigned,
            result["pct_assigned"],
            mean_per_spot,
            result["ratio_vs_xenium_per_spot"],
        )

    # Print comparison table
    logger.info("")
    logger.info("=" * 90)
    logger.info(
        "COMPARISON TABLE — Region %d (Xenium GT: %d cells, %.1f/spot)", region_id, gt_n_cells, gt_cells_per_spot
    )
    logger.info("-" * 90)
    logger.info(
        "%-22s %8s %10s %8s %8s %10s %8s", "Condition", "Total", "Assigned", "/spot", "Ratio", "Area<30", "Time(s)"
    )
    logger.info("-" * 90)
    for cond_name, r in results.items():
        area_below = r["area_stats"].get("pct_below_30um2", 0)
        logger.info(
            "%-22s %8d %10d %8.1f %8.2fx %9.1f%% %8.1f",
            cond_name,
            r["n_total_nuclei"],
            r["n_assigned_to_spots"],
            r["mean_nuclei_per_spot"],
            r["ratio_vs_xenium_per_spot"],
            area_below,
            r["elapsed_sec"],
        )
    logger.info("-" * 90)
    logger.info("Xenium GT: %8s %10d %8.1f %8s", "—", gt_n_cells, gt_cells_per_spot, "1.00x")
    logger.info("=" * 90)

    # Save results
    output_dir = Path(__file__).resolve().parent.parent / "output" / "stardist_diagnostic"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"region_{region_id}_diagnostic.json"

    output = {
        "region_id": region_id,
        "xenium_gt": {
            "n_cells": gt_n_cells,
            "n_spots": gt_n_spots,
            "mean_cells_per_spot": gt_cells_per_spot,
        },
        "geometry": {
            "pixel_size_um": pixel_size,
            "spot_diameter_um": SPOT_DIAMETER_UM,
            "spot_radius_px": float(spot_radius_px),
            "scale_factor": float(XENIUM_PIXEL_SIZE_UM / STARDIST_TRAINING_PIXEL_SIZE),
        },
        "conditions": results,
    }
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    logger.info("Results saved to: %s", output_path)


if __name__ == "__main__":
    main()
