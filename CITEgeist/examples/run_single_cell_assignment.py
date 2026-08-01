"""
Single-Cell Assignment Pipeline for Patient Data.

Runs per-sample:
  Stage 1: StarDist nuclei segmentation on H&E
  Stage 5: Cell assignment (Hungarian) using pre-computed cuOPT proportions, writing
    cell_assignments.csv; then SACE per-cell GEX allocation
    (run_sace_allocation(output_mode="single_cell", use_morphology=False)), which also
    assigns by count-constrained Hungarian so {sample}_single_cell.h5ad agrees with the
    CSV on assignment method. Note: SACE's single_cell mode always re-runs its own
    internal StarDist segmentation, so cell IDs in the h5ad are not the same nuclei as
    Stage 1's nucleus_spot_mapping.csv / cell_assignments.csv.
  Stage 6: Visualization (spatial overlays)

Usage:
  python run_single_cell_assignment.py --sample sample-P1-S1 [--stages 1,5,6]
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
MODEL_DIR = REPO_ROOT / "CITEgeist" / "model"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(MODEL_DIR))
sys.path.insert(0, str(MODEL_DIR / "morphology"))

# CITEgeist public data — download from NCBI GEO accession GSE289326 and set this path
PATIENT_DATA_ROOT = Path("/path/to/CITEgeist_public_data/processed_files")
MODULE3_ROOT = REPO_ROOT / "output" / "module3_7type"
OUTPUT_ROOT = REPO_ROOT / "output" / "morphology_assignment"

# 12 deduplicated patient samples
SAMPLES = [
    "sample-P1-S1",
    "sample-P1-S2",
    "sample-P2-S1",
    "sample-P2-S2",
    "sample-P3-S1",
    "sample-P3-S2",
    "sample-P4-S1",
    "sample-P4-S2",
    "sample-P5-S1",
    "sample-P5-S2",
    "sample-P6-S1",
    "sample-P6-S2",
]

# Unified cell type profile (must match Module 3 output columns)
CELL_TYPES = [
    "Cancer_Cells",
    "Macrophages",
    "CD4_T_Cells",
    "CD8_T_Cells",
    "B_Cells",
    "Endothelial",
    "Fibroblasts",
]

# Column rename map: Profile_N -> cell type name (matches Module 1-2 discovery order)
PROFILE_TO_CELLTYPE = {f"Profile_{i+1}": ct for i, ct in enumerate(CELL_TYPES)}

# Unified marker profile dict (marker name -> list of cell types)
# Marker → cell type mapping (for Stage 4 binary reconstruction matrix)
# Matches Module 1-2 discovery output across all 12 patient samples
PROFILE_DICT = {
    "EPCAM-1": ["Cancer_Cells"],
    "CD68-1": ["Macrophages"],
    "CD14-1": ["Macrophages"],
    "CD3E-1": ["CD8_T_Cells", "CD4_T_Cells"],
    "CD8A-1": ["CD8_T_Cells"],
    "CD4-1": ["CD4_T_Cells"],
    "MS4A1-1": ["B_Cells"],
    "CD19-1": ["B_Cells"],
    "PECAM1-1": ["Endothelial"],
    "ACTA2-1": ["Fibroblasts"],
}

# CitegeistModel format: cell type → {"Major": [markers], "Minor": [markers]}
# EPCAM-1 is the sole Major cancer discriminator; SDC1-1 and KRT5-1 are Minor
# (lower weight in QP, available for M3.5 functional gating downstream)
MODEL_PROFILE_DICT = {
    "Cancer_Cells": {"Major": ["EPCAM-1"], "Minor": ["SDC1-1", "KRT5-1"]},
    "Macrophages": {"Major": ["CD68-1"], "Minor": ["CD14-1"]},
    "CD4_T_Cells": {"Major": ["CD3E-1", "CD4-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "B_Cells": {"Major": ["MS4A1-1", "CD19-1"]},
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
}

# Sample metadata for Module 3 preprocessing
SAMPLE_METADATA = {s: {"min_counts": 100 if "S1" in s else 25} for s in SAMPLES}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def _resolve_he_image(sample_name: str):
    """Resolve the best available H&E image and its coordinate scale factor.

    Prefers high-res tissue_fullres_image.tif (original microscopy, ~10-18k px)
    over the CytAssist image (3000x3000, often downsampled).

    Returns:
        (he_path, scale_factor) where scale_factor maps fullres coords -> image coords
    """
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    sf_path = patient_dir / "spatial" / "scalefactors_json.json"
    import json as _json

    with open(sf_path) as f:
        scalefactors = _json.load(f)

    # Prefer full-res original microscopy image
    fullres_path = patient_dir / "spatial" / "tissue_fullres_image.tif"
    if fullres_path.exists():
        # Spots are in full-res coordinates — no scaling needed
        return fullres_path, 1.0, scalefactors

    # Fallback to CytAssist (downsampled)
    cytassist_path = patient_dir / "spatial" / "cytassist_image.tiff"
    if cytassist_path.exists():
        scale = scalefactors.get("regist_target_img_scalef", 1.0)
        return cytassist_path, scale, scalefactors

    # Last resort: hires PNG
    hires_path = patient_dir / "spatial" / "tissue_hires_image.png"
    assert hires_path.exists(), f"No H&E image found for {sample_name}"
    scale = scalefactors["tissue_hires_scalef"]
    return hires_path, scale, scalefactors


# ===== STAGE 1: StarDist Nuclei Segmentation =============================


def stage1_segment(
    sample_name: str,
    output_dir: Path,
    gpu: bool = True,
    resolution_mode: str = "hires",
    **stardist_kwargs,
) -> Dict:
    """Segment nuclei from H&E image using StarDist with 3-layer defense.

    Uses the production StarDist segmenter from CITEgeist.model.segmentation
    with automatic pixel size estimation from Visium scalefactors:
      - Layer 1 (Physics): scale = pixel_size_um / 0.25 (auto from scalefactors)
      - Layer 2 (Model): prob_thresh = 0.6 (conservative, removes debris)
      - Layer 3 (Biology): area filter 20-500 µm² (removes sub-nuclear + merged)
    """
    import json

    import squidpy as sq
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None

    logger.info("Stage 1: StarDist segmentation for %s", sample_name)

    # Load spatial coordinates via squidpy (with images for segmentation)
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    adata = sq.read.visium(
        str(patient_dir), counts_file="filtered_feature_bc_matrix.h5", load_images=True, gex_only=False
    )
    adata.var_names_make_unique()

    # Filter spots with NaN spatial coordinates (common in surgical samples)
    if "spatial" in adata.obsm:
        coords = adata.obsm["spatial"]
        valid_mask = np.all(np.isfinite(coords), axis=1)
        n_invalid = (~valid_mask).sum()
        if n_invalid > 0:
            logger.warning("  Filtering %d spots with NaN/inf spatial coordinates", n_invalid)
            adata = adata[valid_mask].copy()

    # Run production StarDist segmentation pipeline.
    # pixel_size_um is auto-estimated from scalefactors (55µm / spot_diameter_fullres).
    # This enables the full 3-layer defense without manual pixel size input.
    from CITEgeist.model.morphology.segmentation import (
        compute_spot_nuclei_counts_patchwise,
        estimate_pixel_size_um,
        save_segmentation_artifacts,
    )

    # Estimate pixel size from scalefactors
    px_size = estimate_pixel_size_um(adata)
    if px_size is None:
        raise RuntimeError("Could not estimate pixel_size_um from scalefactors")
    logger.info("  Auto pixel_size: %.4f µm/px (fullres frame)", px_size)

    # Load fullres TIF for patchwise segmentation at native resolution.
    # Each 55µm spot patch is ~70px wide at 0.79 µm/px — StarDist processes
    # each patch individually, keeping memory bounded.
    from PIL import Image

    fullres_path = patient_dir / "spatial" / "tissue_fullres_image.tif"
    if not fullres_path.exists():
        raise FileNotFoundError(f"No fullres TIF found at {fullres_path}")

    logger.info("  Loading fullres image: %s", fullres_path)
    fullres_img = np.array(Image.open(fullres_path).convert("RGB"))
    logger.info("  Fullres shape: %s, StarDist scale: %.2f", fullres_img.shape, px_size / 0.25)

    # Get spot coordinates (fullres pixel space)
    spot_centers = np.asarray(adata.obsm["spatial"], dtype=np.float64)

    result = compute_spot_nuclei_counts_patchwise(
        fullres_image=fullres_img,
        spot_centers_fullres=spot_centers,
        spot_names=adata.obs_names,
        pixel_size_um=px_size,
        modality="he",
        **stardist_kwargs,
    )

    # Free fullres image
    del fullres_img

    # Save artifacts (mask, centroids, QC)
    artifacts = save_segmentation_artifacts(
        output_folder=str(output_dir),
        sample_name=sample_name,
        result=result,
        save_masks=True,
    )
    logger.info("  Segmentation artifacts saved: %s", list(artifacts.keys()))

    # Save nucleus-spot mapping in the format expected by downstream stages
    if result.nucleus_spot_map is not None:
        nsm = result.nucleus_spot_map.rename(
            columns={
                "spot_id": "spot_barcode",
                "x_pixel": "centroid_x",
                "y_pixel": "centroid_y",
            }
        )
        nsm.to_csv(output_dir / "nucleus_spot_mapping.csv", index=False)
        n_assigned = len(nsm)
        n_spots = nsm["spot_barcode"].nunique()
    else:
        n_assigned = 0
        n_spots = 0

    n_nuclei = int(result.nuclei_count_raw.sum())
    logger.info(
        "  Detected %d nuclei, %d assigned to %d spots (avg %.1f/spot)",
        n_nuclei,
        n_assigned,
        n_spots,
        n_assigned / max(1, n_spots),
    )

    # Log QC summary
    qc_dict = {}
    if result.qc is not None:
        qc_dict = result.qc.to_dict()
        logger.info(
            "  QC: %d raw → %d after filter (-%d small, -%d large), " "density=%.0f/mm² [%s], median area=%.1f µm²",
            result.qc.n_raw,
            result.qc.n_after_area_filter,
            result.qc.n_removed_small,
            result.qc.n_removed_large,
            result.qc.density_per_mm2,
            result.qc.density_flag,
            result.qc.area_median_um2,
        )
        # Save QC JSON
        qc_path = output_dir / f"{sample_name}_segmentation_qc.json"
        with open(qc_path, "w") as f:
            json.dump(qc_dict, f, indent=2)

    # Also save mask as .npy for Stage 2 ViT patch extraction compatibility
    mask_path = output_dir / "nuclei_mask.npy"
    np.save(mask_path, result.masks)

    return {
        "mask_path": str(mask_path),
        "n_nuclei": n_nuclei,
        "n_assigned": n_assigned,
        "adata_path": str(patient_dir),
        "image_shape": list(result.image_shape),
        "resolution_mode": resolution_mode,
        "qc": qc_dict,
    }


# ===== STAGE 5: Cell Assignment + GEX ===================================


def stage5_assign_and_gex(
    sample_name: str,
    output_dir: Path,
) -> Dict:
    """Cell assignment (Hungarian, default) + SACE GEX deconvolution.

    Loads pre-computed cuOPT proportions from module3_unified instead of
    re-running the GPU QP solve.
    """
    import sys as _sys

    import squidpy as sq

    _sys.path.insert(0, str(REPO_ROOT))
    from CITEgeist.model.citegeist_model import CitegeistModel

    logger.info("Stage 5: Cell assignment + GEX for %s", sample_name)

    # --- 1. Load pre-computed cuOPT proportions ---
    prop_path = MODULE3_ROOT / f"{sample_name}_cell_prop_finetuned_results.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Module 3 proportions not found: {prop_path}")
    prop_df = pd.read_csv(prop_path, index_col=0)
    if "Profile_1" in prop_df.columns:
        prop_df = prop_df.rename(columns=PROFILE_TO_CELLTYPE)
        if "Unknown" in prop_df.columns:
            prop_df = prop_df.drop(columns=["Unknown"])
    cell_types = list(prop_df.columns)
    n_spots = len(prop_df)
    logger.info("  Loaded proportions: %d spots x %d types", n_spots, len(cell_types))

    # --- 2. Load nucleus-spot mapping ---
    nuc_map = pd.read_csv(output_dir / "nucleus_spot_mapping.csv")
    # Build positional cell_to_spot: integer index into prop_df row order
    spot_order = {bc: i for i, bc in enumerate(prop_df.index)}
    nuc_map = nuc_map[nuc_map["spot_barcode"].isin(spot_order)].copy()
    nuc_map["spot_idx"] = nuc_map["spot_barcode"].map(spot_order)
    cell_ids = nuc_map["nucleus_id"].astype(str).values
    cell_to_spot = nuc_map["spot_idx"].values.astype(int)

    # nuclei_counts: per-spot count Series aligned to prop_df index
    counts_series = nuc_map.groupby("spot_barcode")["nucleus_id"].count()
    nuclei_counts = counts_series.reindex(prop_df.index, fill_value=0)

    logger.info("  %d nuclei across %d spots", len(cell_ids), (nuclei_counts > 0).sum())

    # --- 3. Build detection mask (full-ones fallback — no GPU QP re-run needed) ---
    detection_mask = np.ones((n_spots, len(cell_types)), dtype=bool)
    logger.info("  Using full-ones detection mask (fallback — no re-deconvolution)")

    # --- 4. Instantiate CitegeistModel and inject results ---
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    adata_raw = sq.read.visium(
        str(patient_dir),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=True,
        gex_only=False,
    )

    meta = SAMPLE_METADATA[sample_name]
    model = CitegeistModel(
        sample_name=sample_name,
        adata=adata_raw,
        output_folder=str(output_dir),
    )
    model.split_adata()

    # run_sace_allocation(output_mode="single_cell") re-runs its own StarDist
    # segmentation and hard-requires a true fullres H&E image — hires/lowres PNGs
    # are too coarse to pass its >=1 µm/px resolution gate. Register the same
    # fullres TIF Stage 1 uses so the single-cell GEX step can actually run.
    #
    # Registered on model.antibody_capture_adata only (the object
    # run_sace_allocation's `source_adata = self.antibody_capture_adata or
    # self.gene_expression_adata or self.adata` selection actually reads), and only
    # AFTER split_adata() — split_adata()'s var-axis .copy() deep-copies .uns, so
    # writing to adata_raw.uns before the split would put the (potentially multi-GB)
    # fullres array in three places (adata_raw, gene_expression_adata,
    # antibody_capture_adata) simultaneously, all persistent. Registering here keeps
    # it to one persistent copy (plus a second, transient one at the `[common].copy()`
    # subset below, which also deep-copies .uns — freed once the old
    # antibody_capture_adata object is garbage collected). Local `fullres_img`
    # binding freed with `del` once assigned, same discipline Stage 1 uses.
    fullres_path = patient_dir / "spatial" / "tissue_fullres_image.tif"
    if fullres_path.exists():
        from PIL import Image as _PILImage

        _PILImage.MAX_IMAGE_PIXELS = None
        lib_id = next(iter(model.antibody_capture_adata.uns["spatial"]))
        fullres_img = np.array(_PILImage.open(fullres_path).convert("RGB"))
        model.antibody_capture_adata.uns["spatial"][lib_id]["images"]["fullres"] = fullres_img
        del fullres_img
    else:
        logger.warning(
            "  No fullres TIF for %s — SACE single-cell GEX allocation will fail " "its H&E resolution gate",
            sample_name,
        )

    # Safety net: load_images=True already populates obsm["spatial"], so this
    # should be a no-op, but keep it in case a library ever omits tissue_positions.
    if "spatial" not in model.gene_expression_adata.obsm:
        pos_path = patient_dir / "spatial" / "tissue_positions.csv"
        pos_df = pd.read_csv(pos_path, index_col=0)
        common_bcs = model.gene_expression_adata.obs_names
        coords = (
            pos_df.loc[pos_df.index.intersection(common_bcs), ["pxl_col_in_fullres", "pxl_row_in_fullres"]]
            .reindex(common_bcs)
            .values
        )
        model.gene_expression_adata.obsm["spatial"] = coords
        logger.info("  Injected spatial coords from tissue_positions.csv")

    model.filter_gex(min_counts=meta["min_counts"])
    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(MODEL_PROFILE_DICT)

    # Subset model to spots that have proportions
    common = prop_df.index.intersection(model.gene_expression_adata.obs_names)
    common_pos = prop_df.index.get_indexer(common)  # positional rows in original prop_df
    model.gene_expression_adata = model.gene_expression_adata[common].copy()
    model.antibody_capture_adata = model.antibody_capture_adata[common].copy()
    prop_df = prop_df.loc[common]
    detection_mask = detection_mask[common_pos]
    # Rebuild cell_to_spot against the subsetted prop_df index
    new_spot_order = {bc: i for i, bc in enumerate(prop_df.index)}
    valid_cells = nuc_map["spot_barcode"].isin(new_spot_order)
    nuc_map = nuc_map[valid_cells].copy()
    cell_ids = nuc_map["nucleus_id"].astype(str).values
    cell_to_spot = nuc_map["spot_barcode"].map(new_spot_order).values.astype(int)

    # Inject pre-computed results
    model.results["cell_prop"] = prop_df
    model.results["detection_mask_bool"] = detection_mask

    # --- 5. Cell assignment ---
    # Bayesian single-cell assignment is available by passing morph_scores_precomputed
    # (a (C, n_types) array); a non-DL score producer is deferred.
    assignments_df = model.assign_cells(
        nuclei_counts=nuclei_counts.reindex(prop_df.index, fill_value=0),
        cell_to_spot=cell_to_spot,
        cell_ids=cell_ids,
        assignment_method="hungarian",
        detection_mask=detection_mask,
    )

    assignments_df.to_csv(output_dir / "cell_assignments.csv", index=False)
    logger.info("  Assigned %d nuclei across %d spots", len(assignments_df), assignments_df["spot_id"].nunique())

    # --- 6. SACE per-cell GEX allocation ---
    # resolution_mode="fullres" — the fullres TIF was registered on
    # model.antibody_capture_adata above (see the memory-discipline comment there);
    # hires/lowres are too coarse to pass the H&E resolution gate.
    # use_morphology=False — keeps this on the same count-constrained Hungarian
    # assignment as Stage 5's assign_cells() call above, so {sample}_single_cell.h5ad
    # agrees with cell_assignments.csv on assignment method (the default
    # use_morphology=True would instead assign by morphology-Bayesian posterior,
    # silently disagreeing with the CSV). Note this does NOT reconcile segmentation:
    # single_cell mode always re-runs its own internal StarDist pass, so cell IDs in
    # the h5ad are still a different nucleus set than nucleus_spot_mapping.csv.
    _, cell_adata, _ = model.run_sace_allocation(
        output_mode="single_cell",
        resolution_mode="fullres",
        use_morphology=False,
    )
    sc_path = output_dir / f"{sample_name}_single_cell.h5ad"
    cell_adata.write_h5ad(str(sc_path))
    logger.info(
        "  Saved single-cell AnnData: %d cells x %d genes → %s", cell_adata.shape[0], cell_adata.shape[1], sc_path
    )

    n_assigned = len(assignments_df)
    type_dist = assignments_df["assigned_type"].value_counts().to_dict()
    json_path = output_dir / "pipeline_results.json"
    import json as _json

    existing = {}
    if json_path.exists():
        with open(json_path) as f:
            existing = _json.load(f)
    existing["stage5"] = {"n_assigned": n_assigned, "type_distribution": type_dist}
    with open(json_path, "w") as f:
        _json.dump(existing, f, indent=2)

    return {"n_assigned": n_assigned, "type_distribution": type_dist}


# ===== STAGE 6: Visualization ===========================================


def stage6_visualize(
    sample_name: str,
    output_dir: Path,
) -> Dict:
    """Generate spatial overlays of cell type assignments."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    logger.info("Stage 6: Visualization for %s", sample_name)
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None

    # Load the image matching the segmentation resolution
    import json as _json

    pipeline_meta_path = output_dir / "pipeline_results.json"
    seg_resolution = "fullres"
    if pipeline_meta_path.exists():
        with open(pipeline_meta_path) as f:
            pm = _json.load(f)
        seg_resolution = pm.get("stage1", {}).get("resolution_mode", seg_resolution)

    if seg_resolution in ("hires", "lowres"):
        import squidpy as sq

        patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
        adata_img = sq.read.visium(str(patient_dir), load_images=True)
        lib = next(iter(adata_img.uns["spatial"].values()))
        he_image = lib["images"][seg_resolution]
        if he_image.dtype != np.uint8:
            he_image = (he_image * 255).astype(np.uint8)
    else:
        he_path, _, _ = _resolve_he_image(sample_name)
        he_image = np.array(Image.open(he_path))

    assignments_df = pd.read_csv(output_dir / "cell_assignments.csv")

    # Load nucleus mapping — try new patchwise name first, fall back to old
    nmap_path = output_dir / f"{sample_name}_nucleus_spot_map.csv"
    if not nmap_path.exists():
        nmap_path = output_dir / "nucleus_spot_mapping.csv"
    nucleus_spot_map = pd.read_csv(nmap_path)

    # Normalize column names for merge
    if "x_pixel" in nucleus_spot_map.columns:
        nucleus_spot_map = nucleus_spot_map.rename(columns={"x_pixel": "centroid_x", "y_pixel": "centroid_y"})
    if "spot_id" in nucleus_spot_map.columns and "spot_barcode" not in nucleus_spot_map.columns:
        nucleus_spot_map = nucleus_spot_map.rename(columns={"spot_id": "spot_barcode"})

    # Merge — assignments uses cell_id, mapping uses nucleus_id (same values)
    merge_key = "nucleus_id" if "nucleus_id" in assignments_df.columns else "cell_id"
    if merge_key == "cell_id":
        assignments_df = assignments_df.rename(columns={"cell_id": "nucleus_id"})
    merged = assignments_df.merge(nucleus_spot_map, on="nucleus_id", how="inner")

    # Color palette for cell types
    type_colors = {
        "Cancer_Cells": "#d62728",
        "Macrophages": "#ff7f0e",
        "CD4_T_Cells": "#e377c2",
        "CD8_T_Cells": "#8c564b",
        "B_Cells": "#2ca02c",
        "Endothelial": "#1f77b4",
        "Fibroblasts": "#9467bd",
    }

    # --- Spatial Overlay ---
    fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150)
    ax.imshow(he_image)
    for _, row in merged.iterrows():
        color = type_colors.get(row["assigned_type"], "#000000")
        ax.scatter(row["centroid_x"], row["centroid_y"], c=color, s=2, alpha=0.6, linewidths=0)
    # Legend
    handles = [mpatches.Patch(color=c, label=t) for t, c in type_colors.items() if t in merged["assigned_type"].values]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.8)
    ax.set_title(f"{sample_name} — Cell Type Assignments", fontsize=12)
    ax.axis("off")
    fig.savefig(output_dir / "spatial_overlay.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("  Saved spatial_overlay.png")

    return {"visualizations": ["spatial_overlay.png"]}


# ===== MAIN =============================================================


def main():
    parser = argparse.ArgumentParser(description="H&E Morphology Assignment Pipeline")
    parser.add_argument("--sample", required=True, help="Sample name (e.g., sample-P1-S1)")
    parser.add_argument("--stages", default="1,5,6", help="Comma-separated stages to run (default: 1,5,6)")
    parser.add_argument("--gpu", action="store_true", default=True, help="Use GPU for StarDist")
    parser.add_argument("--no-gpu", action="store_true", help="Disable GPU")
    parser.add_argument(
        "--resolution-mode",
        default="hires",
        choices=["lowres", "hires", "fullres"],
        help="Image resolution for StarDist segmentation (default: hires)",
    )
    parser.add_argument(
        "--module3-dir",
        type=str,
        default=None,
        help="Override Module 3 proportion directory " "(default: output/module3_unified)",
    )
    parser.add_argument(
        "--output-dir-override",
        type=str,
        default=None,
        help="Override output directory (default: output/morphology_assignment/{sample})",
    )
    args = parser.parse_args()

    if args.no_gpu:
        args.gpu = False

    stages = [int(s) for s in args.stages.split(",")]
    sample = args.sample

    assert sample in SAMPLES, f"Unknown sample: {sample}. Must be one of: {SAMPLES}"

    output_dir = OUTPUT_ROOT / sample
    output_dir.mkdir(parents=True, exist_ok=True)

    # Override paths if specified
    if args.module3_dir:
        global MODULE3_ROOT
        MODULE3_ROOT = Path(args.module3_dir)
        logger.info("Using Module 3 dir: %s", MODULE3_ROOT)

    if args.output_dir_override:
        output_dir = Path(args.output_dir_override)
        output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    if 1 in stages:
        results["stage1"] = stage1_segment(sample, output_dir, gpu=args.gpu, resolution_mode=args.resolution_mode)
    if 5 in stages:
        results["stage5"] = stage5_assign_and_gex(sample, output_dir)
    if 6 in stages:
        if not (output_dir / "cell_assignments.csv").exists():
            logger.error("No cell assignments — skipping stage 6.")
        else:
            results["stage6"] = stage6_visualize(sample, output_dir)

    # Save run metadata
    with open(output_dir / "pipeline_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info("Pipeline complete for %s. Results in %s", sample, output_dir)


if __name__ == "__main__":
    main()
