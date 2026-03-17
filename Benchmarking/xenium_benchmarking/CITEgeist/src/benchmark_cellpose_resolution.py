#!/usr/bin/env python
"""
Benchmark nuclei-count fidelity across image resolutions on Xenium pseudo-Visium.

This version reads morphology images from explicit filesystem paths (not h5ad uns['spatial']).
Uses StarDist for nuclei segmentation via ``compute_spot_nuclei_counts()``.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import cv2
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.citegeist_model import CitegeistModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def _prepare_rgb_uint8(image: np.ndarray) -> np.ndarray:
    if image.ndim == 2:
        image = np.stack([image, image, image], axis=-1)
    if image.ndim != 3:
        raise ValueError(f"Expected 2D/3D image, got shape {image.shape}")
    if image.shape[-1] == 4:
        image = image[..., :3]
    if image.dtype == np.uint8:
        return image
    image = np.asarray(image, dtype=np.float32)
    vmin = float(np.nanmin(image))
    vmax = float(np.nanmax(image))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        return np.zeros(image.shape, dtype=np.uint8)
    return (255.0 * np.clip((image - vmin) / (vmax - vmin), 0.0, 1.0)).astype(np.uint8)


def _load_rgb_image(image_path: Path) -> np.ndarray:
    if not image_path.exists():
        raise FileNotFoundError(f"Missing image: {image_path}")
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    if bgr is None:
        raise ValueError(f"Failed to load image: {image_path}")
    return cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)


def _load_fullres_crop_if_available(fullres_path: Optional[Path], coord: Dict) -> Optional[np.ndarray]:
    if fullres_path is None or not fullres_path.exists():
        return None
    try:
        import tifffile
        import zarr
    except Exception:
        logger.warning("tifffile/zarr not available; skipping fullres crop loading.")
        return None

    x0 = int(coord["pixel_bounds"]["x_min"])
    x1 = int(coord["pixel_bounds"]["x_max"])
    y0 = int(coord["pixel_bounds"]["y_min"])
    y1 = int(coord["pixel_bounds"]["y_max"])

    try:
        with tifffile.TiffFile(str(fullres_path)) as tf:
            z = zarr.open(tf.series[0].aszarr(), mode="r")
            # tifffile.aszarr can expose a Group with level keys (e.g., "0","1",...).
            if hasattr(z, "keys"):
                keys = sorted(list(z.keys()), key=str)
                if len(keys) == 0:
                    raise ValueError("Empty zarr group from fullres OME store.")
                z = z[keys[0]]

            if z.ndim == 2:
                crop = np.asarray(z[y0:y1, x0:x1])
            elif z.ndim == 3:
                # Expected Xenium morphology stack is typically ZYX. Use middle Z plane.
                if z.shape[0] <= 64 and z.shape[1] > 256 and z.shape[2] > 256:
                    z_idx = int(z.shape[0] // 2)
                    crop = np.asarray(z[z_idx, y0:y1, x0:x1])
                elif z.shape[-1] <= 8:
                    crop = np.asarray(z[y0:y1, x0:x1, 0])
                else:
                    crop = np.asarray(z[y0:y1, x0:x1])
            elif z.ndim == 4:
                # Possible TZYX/ZYXC layout: prefer first T and middle Z.
                if z.shape[0] <= 8 and z.shape[1] <= 64:
                    z_idx = int(z.shape[1] // 2)
                    crop = np.asarray(z[0, z_idx, y0:y1, x0:x1])
                elif z.shape[-1] <= 8:
                    crop = np.asarray(z[0, y0:y1, x0:x1, 0])
                else:
                    raise ValueError(f"Unsupported 4D layout for fullres crop: {z.shape}")
            else:
                logger.warning("Unsupported fullres image ndim=%s; skipping fullres crop.", z.ndim)
                return None
        return _prepare_rgb_uint8(crop)
    except Exception as exc:
        logger.warning("Failed to read fullres crop from %s: %r", fullres_path, exc)
        return None


def _build_external_spatial_adata(
    adata: sc.AnnData,
    region_id: int,
    lowres_image_root: Path,
    hires_image_root: Path,
    coord_info_root: Path,
    fullres_image_path: Optional[Path],
    spot_diameter_um: float,
) -> sc.AnnData:
    """
    Build AnnData with spatial metadata from external image files.

    Args:
        adata: Source AnnData with spatial coordinates in microns
        region_id: Xenium region ID
        lowres_image_root: Root folder for lowres morphology images
        hires_image_root: Root folder for hires morphology images
        coord_info_root: Root folder for coord_info.json files
        fullres_image_path: Optional path to full-resolution morphology OME-TIFF
        spot_diameter_um: Spot diameter in microns (e.g., 55.0 for Visium geometry)

    Returns:
        AnnData with proper spatial metadata for segmentation
    """
    region_name = f"Xenium_region_{region_id}"
    coord_info_path = coord_info_root / region_name / "coord_info.json"
    lowres_path = lowres_image_root / region_name / "morphology.png"
    hires_path = hires_image_root / region_name / "morphology.png"

    if not coord_info_path.exists():
        raise FileNotFoundError(f"Missing coord_info: {coord_info_path}")

    with open(coord_info_path) as f:
        coord = json.load(f)

    lowres_img = _load_rgb_image(lowres_path)
    hires_img = _load_rgb_image(hires_path)
    fullres_img = _load_fullres_crop_if_available(fullres_image_path, coord)

    pixel_size = float(coord["pixel_size"])  # microns/pixel in hires crop frame
    x_min = float(coord["micron_bounds"]["x_min"])
    y_min = float(coord["micron_bounds"]["y_min"])
    full_w = float(coord["image_size"][0])

    # Convert spot coordinates from microns to pixels (in hires image frame)
    spot_microns = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    spot_fullres = np.empty_like(spot_microns, dtype=np.float64)
    spot_fullres[:, 0] = (spot_microns[:, 0] - x_min) / pixel_size
    spot_fullres[:, 1] = (spot_microns[:, 1] - y_min) / pixel_size

    # Convert spot diameter from microns to pixels
    spot_diameter_fullres = spot_diameter_um / pixel_size
    logger.info(
        "Spot diameter: %.1f µm = %.1f pixels (pixel_size=%.4f µm/px)",
        spot_diameter_um, spot_diameter_fullres, pixel_size
    )

    scalefactors = {
        "spot_diameter_fullres": spot_diameter_fullres,
        "tissue_hires_scalef": float(hires_img.shape[1]) / full_w,
        "tissue_lowres_scalef": float(lowres_img.shape[1]) / full_w,
    }

    images = {"lowres": _prepare_rgb_uint8(lowres_img), "hires": _prepare_rgb_uint8(hires_img)}
    if fullres_img is not None:
        images["fullres"] = fullres_img

    out = adata.copy()
    out.obsm["spatial"] = spot_fullres
    out.uns["spatial"] = {region_name: {"images": images, "scalefactors": scalefactors}}
    out.uns["platform"] = "xenium_pseudovisium"
    return out


def _align_true_counts_to_obs(obs_names: pd.Index, true_counts: pd.Series) -> pd.Series:
    """Align GT counts to AnnData obs names with simple spot_ prefix fallback."""
    if obs_names.isin(true_counts.index).sum() > 0:
        return true_counts.reindex(obs_names).fillna(0.0)

    if obs_names.str.startswith("spot_").all() and not true_counts.index.astype(str).str.startswith("spot_").all():
        tc = true_counts.copy()
        tc.index = "spot_" + tc.index.astype(str)
        return tc.reindex(obs_names).fillna(0.0)

    if true_counts.index.astype(str).str.startswith("spot_").all() and not obs_names.str.startswith("spot_").all():
        tc = true_counts.copy()
        tc.index = tc.index.astype(str).str.replace("spot_", "", regex=False)
        return tc.reindex(obs_names).fillna(0.0)

    return true_counts.reindex(obs_names).fillna(0.0)


def compute_true_cells_per_spot(mapping_csv: Path, region_id: int) -> pd.Series:
    df = pd.read_csv(mapping_csv, index_col=0)
    required = {"spot_idx", "spot_id", "region_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required mapping columns: {missing}")

    spot_idx_num = pd.to_numeric(df["spot_idx"], errors="coerce")
    region_num = pd.to_numeric(df["region_id"], errors="coerce")

    valid = spot_idx_num.notna() & (spot_idx_num >= 0)
    valid &= df["spot_id"].notna()
    valid &= region_num.notna()
    valid &= region_num == float(region_id)
    region_df = df.loc[valid].copy()
    counts = region_df.groupby("spot_id").size().astype(float)
    counts.name = "true_cells_per_spot"
    return counts


def _calc_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    eps = 1e-8
    mae = float(np.mean(np.abs(y_pred - y_true)))
    rmse = float(np.sqrt(np.mean((y_pred - y_true) ** 2)))
    mape = float(np.mean(np.abs((y_pred - y_true) / np.maximum(np.abs(y_true), eps))))
    pr = pearsonr(y_true, y_pred)
    sr = spearmanr(y_true, y_pred)
    return {
        "pearson_r": float(pr.statistic),
        "pearson_p": float(pr.pvalue),
        "spearman_rho": float(sr.statistic),
        "spearman_p": float(sr.pvalue),
        "mae": mae,
        "rmse": rmse,
        "mape": mape,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark Cellpose resolution fidelity on Xenium pseudo-Visium.")
    parser.add_argument("--region-id", type=int, default=0, help="Xenium pseudo-Visium region ID.")
    parser.add_argument(
        "--h5ad-dir",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects",
        help="Directory containing Xenium_region_<id>_CITE.h5ad/GEX.h5ad.",
    )
    parser.add_argument(
        "--mapping-csv",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_to_spot_mapping.csv",
        help="Path to cell_to_spot_mapping.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_cellpose_resolution",
        help="Output directory for benchmark metrics/artifacts.",
    )
    parser.add_argument(
        "--lowres-image-root",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "scResolve" / "images" / "morphology",
        help="Root folder containing lowres region images.",
    )
    parser.add_argument(
        "--hires-image-root",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "scResolve" / "images" / "morphology_hires",
        help="Root folder containing hires region images.",
    )
    parser.add_argument(
        "--coord-info-root",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "scResolve" / "images" / "morphology_hires",
        help="Root folder containing coord_info.json files.",
    )
    parser.add_argument(
        "--fullres-image-path",
        type=Path,
        default=Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/morphology.ome.tif"),
        help="Path to slide-level fullres morphology OME-TIFF; optional.",
    )
    parser.add_argument(
        "--spot-diameter-um",
        type=float,
        default=55.0,
        help="Spot diameter in microns. Standard Visium = 55.0 µm, Visium HD = 8.0 µm.",
    )
    parser.add_argument("--resolutions", nargs="+", default=["lowres", "hires", "fullres"], help="Modes to benchmark.")
    parser.add_argument("--max-fullres-side", type=int, default=9000, help="Max side for fallback upsampling.")
    parser.add_argument(
        "--fullres-patch-radius-multiplier",
        type=float,
        default=1.5,
        help="Patch radius multiplier (spot radius units) for per-spot fullres segmentation.",
    )
    parser.add_argument(
        "--fullres-patch-workers",
        type=int,
        default=4,
        help="Number of parallel workers for fullres patch segmentation.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cite_path = args.h5ad_dir / f"Xenium_region_{args.region_id}_CITE.h5ad"
    gex_path = args.h5ad_dir / f"Xenium_region_{args.region_id}_GEX.h5ad"
    if not cite_path.exists() or not gex_path.exists():
        raise FileNotFoundError(f"Missing CITE/GEX h5ad for region {args.region_id} in {args.h5ad_dir}")

    logger.info("Loading %s", cite_path)
    adata_cite = sc.read_h5ad(cite_path)
    adata_gex = sc.read_h5ad(gex_path)

    adata_for_seg = _build_external_spatial_adata(
        adata=adata_gex,
        region_id=args.region_id,
        lowres_image_root=args.lowres_image_root,
        hires_image_root=args.hires_image_root,
        coord_info_root=args.coord_info_root,
        fullres_image_path=args.fullres_image_path,
        spot_diameter_um=args.spot_diameter_um,
    )

    true_counts = compute_true_cells_per_spot(args.mapping_csv, args.region_id)
    true_aligned = _align_true_counts_to_obs(adata_cite.obs_names, true_counts)

    rows: List[Dict[str, float]] = []
    for mode in args.resolutions:
        logger.info("Running Cellpose resolution mode: %s", mode)
        mode_out = args.output_dir / f"region_{args.region_id}" / mode
        mode_out.mkdir(parents=True, exist_ok=True)

        model = CitegeistModel(
            sample_name=f"Xenium_region_{args.region_id}",
            output_folder=str(mode_out),
            simulation=True,
            gene_expression_adata=adata_for_seg.copy(),
            antibody_capture_adata=adata_for_seg.copy(),
        )
        t0 = time.time()
        pred_counts = model.compute_spot_nuclei_counts(
            resolution_mode=mode,
            max_fullres_side=args.max_fullres_side,
            fullres_patch_mode=(mode == "fullres"),
            fullres_patch_radius_multiplier=args.fullres_patch_radius_multiplier,
            fullres_patch_workers=args.fullres_patch_workers,
            save_masks=False,
        )
        runtime_sec = float(time.time() - t0)

        pred_aligned = pred_counts.reindex(adata_cite.obs_names).fillna(0.0).astype(float)
        y_true = true_aligned.to_numpy(dtype=float)
        y_pred = pred_aligned.to_numpy(dtype=float)
        metrics = _calc_metrics(y_true=y_true, y_pred=y_pred)
        metrics.update(
            {
                "region_id": int(args.region_id),
                "resolution_mode": mode,
                "runtime_sec": runtime_sec,
                "pred_total_nuclei": float(np.sum(y_pred)),
                "true_total_cells": float(np.sum(y_true)),
            }
        )
        rows.append(metrics)

        compare_df = pd.DataFrame(
            {"spot": adata_cite.obs_names, "true_cells_per_spot": y_true, "pred_nuclei_count": y_pred}
        ).set_index("spot")
        compare_df.to_csv(mode_out / f"region_{args.region_id}_counts_comparison.csv")

    metrics_df = pd.DataFrame(rows)
    metrics_csv = args.output_dir / f"region_{args.region_id}_resolution_fidelity_metrics.csv"
    metrics_json = args.output_dir / f"region_{args.region_id}_resolution_fidelity_summary.json"
    metrics_df.to_csv(metrics_csv, index=False)
    with open(metrics_json, "w") as f:
        json.dump(metrics_df.to_dict(orient="records"), f, indent=2)
    logger.info("Saved metrics: %s", metrics_csv)
    logger.info("Saved summary: %s", metrics_json)


if __name__ == "__main__":
    main()
