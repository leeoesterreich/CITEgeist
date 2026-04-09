"""
H&E Morphology Single-Cell Assignment Pipeline for Patient Data.

Runs per-sample:
  Stage 1: StarDist nuclei segmentation on H&E
  Stage 2: ViT-Small patch extraction + feature embedding
  Stage 5: Bayesian cell assignment using pre-computed cuOPT proportions + SACE GEX
  Stage 6: Visualization (spatial overlays)

Usage:
  python run_morphology_assignment.py --sample HCC22-088-P1-S1 [--stages 1,2,5,6]
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
MODEL_DIR = REPO_ROOT / "CITEgeist" / "model"
BENCHMARK_SRC = REPO_ROOT / "Benchmarking" / "visiumhd_benchmarking" / "src"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(MODEL_DIR))
sys.path.insert(0, str(BENCHMARK_SRC))

PATIENT_DATA_ROOT = Path(
    "/ix1/alee/LO_LAB/General/Lab_Data/"
    "20250210_CITEGeistPublicData_GEO_Alex/processed_files"
)
MODULE3_ROOT = REPO_ROOT / "output" / "module3_unified"
OUTPUT_ROOT = REPO_ROOT / "output" / "morphology_assignment"

# 12 deduplicated patient samples
SAMPLES = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]

# Unified cell type profile (must match Module 3 output columns)
CELL_TYPES = [
    "Endothelial", "Fibroblasts", "B_Cells", "Macrophages", "Monocytes",
    "CD8_T_Cells", "CD4_T_Cells", "Cancer_Luminal", "Cancer_Basal",
    "Dendritic_Cells",
]

# Unified marker profile dict (marker name -> list of cell types)
# Marker → cell type mapping (for Stage 4 binary reconstruction matrix)
PROFILE_DICT = {
    "PECAM1-1": ["Endothelial"],
    "ACTA2-1": ["Fibroblasts"],
    "CD19-1": ["B_Cells"],
    "CD68-1": ["Macrophages"],
    "CD163-1": ["Macrophages"],
    "CD14-1": ["Monocytes"],
    "CD8A-1": ["CD8_T_Cells"],
    "CD4-1": ["CD4_T_Cells"],
    "CD3E-1": ["CD8_T_Cells", "CD4_T_Cells"],
    "EPCAM-1": ["Cancer_Luminal", "Cancer_Basal"],
    "KRT5-1": ["Cancer_Basal"],
    "SDC1-1": ["Cancer_Basal"],
    "ITGAX-1": ["Dendritic_Cells"],
    "HLA-DRA-1": ["Dendritic_Cells"],
}

# CitegeistModel format: cell type → {"Major": [markers]} (for Stage 5 SACE)
MODEL_PROFILE_DICT = {
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "CD8_T_Cells": {"Major": ["CD8A-1", "CD3E-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1", "CD3E-1"]},
    "Cancer_Luminal": {"Major": ["EPCAM-1"]},
    "Cancer_Basal": {"Major": ["KRT5-1", "SDC1-1", "EPCAM-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1", "HLA-DRA-1"]},
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
    import squidpy as sq
    import json
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None

    logger.info("Stage 1: StarDist segmentation for %s", sample_name)

    # Load spatial coordinates via squidpy (with images for segmentation)
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    adata = sq.read.visium(str(patient_dir), counts_file="filtered_feature_bc_matrix.h5",
                           load_images=True, gex_only=False)
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
        save_segmentation_artifacts,
        estimate_pixel_size_um,
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
    logger.info("  Fullres shape: %s, StarDist scale: %.2f",
                fullres_img.shape, px_size / 0.25)

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
        nsm = result.nucleus_spot_map.rename(columns={
            "spot_id": "spot_barcode",
            "x_pixel": "centroid_x",
            "y_pixel": "centroid_y",
        })
        nsm.to_csv(output_dir / "nucleus_spot_mapping.csv", index=False)
        n_assigned = len(nsm)
        n_spots = nsm["spot_barcode"].nunique()
    else:
        n_assigned = 0
        n_spots = 0

    n_nuclei = int(result.nuclei_count_raw.sum())
    logger.info("  Detected %d nuclei, %d assigned to %d spots (avg %.1f/spot)",
                n_nuclei, n_assigned, n_spots,
                n_assigned / max(1, n_spots))

    # Log QC summary
    qc_dict = {}
    if result.qc is not None:
        qc_dict = result.qc.to_dict()
        logger.info("  QC: %d raw → %d after filter (-%d small, -%d large), "
                     "density=%.0f/mm² [%s], median area=%.1f µm²",
                     result.qc.n_raw, result.qc.n_after_area_filter,
                     result.qc.n_removed_small, result.qc.n_removed_large,
                     result.qc.density_per_mm2, result.qc.density_flag,
                     result.qc.area_median_um2)
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


# ===== STAGE 2: Patch Extraction + ViT Features =========================

def stage2_extract_features(
    sample_name: str,
    output_dir: Path,
    vit_model: str = "vit_small_patch16_224",
    batch_size: int = 256,
    device: str = "cuda",
    patch_size: int = 224,
    expansion: float = 0.75,
) -> Dict:
    """Extract ViT embeddings for each nucleus patch."""
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    from vit_extractor import ViTFeatureExtractor
    from extract_patches_he import extract_nucleus_patch, normalize_he_patch

    logger.info("Stage 2: ViT feature extraction for %s", sample_name)

    # Load fullres image for patch extraction at native resolution
    he_path, _, scalefactors = _resolve_he_image(sample_name)
    logger.info("  Using image: %s", he_path.name)
    image = np.array(Image.open(he_path).convert("RGB"))

    # Determine pixel size for centroid-based patch extraction
    spot_diam_fullres = scalefactors.get("spot_diameter_fullres", 69.2)
    # If using fullres image (scale=1.0), pixel_size = 55/spot_diam_fullres
    _, img_scale, _ = _resolve_he_image(sample_name)
    pixel_size_um = 55.0 / (spot_diam_fullres * img_scale)
    logger.info("  Image resolution: %.4f µm/px", pixel_size_um)

    # Check if mask-based or centroid-based extraction
    mask_path = output_dir / "nuclei_mask.npy"
    nucleus_spot_map = pd.read_csv(output_dir / "nucleus_spot_mapping.csv")
    has_centroids = "centroid_x" in nucleus_spot_map.columns

    use_centroid_mode = False
    if mask_path.exists():
        mask = np.load(mask_path)
        # Patchwise segmentation produces a dummy (1,1) mask
        if mask.shape[0] <= 1 or mask.shape[1] <= 1:
            use_centroid_mode = True
            logger.info("  Patchwise segmentation detected — using centroid-based patch extraction")
        else:
            logger.info("  Using mask-based patch extraction (mask shape: %s)", mask.shape)
    else:
        use_centroid_mode = True
        logger.info("  No mask found — using centroid-based patch extraction")

    if use_centroid_mode:
        from extract_patches_he import extract_centroid_patch, normalize_he_patch
    else:
        from extract_patches_he import extract_nucleus_patch, normalize_he_patch

    # Load Module 3 proportions for this sample
    # Support both flat (module3_unified/) and nested (module3_cuopt_qp/{sample}/) layouts
    prop_path = MODULE3_ROOT / f"{sample_name}_cell_prop_global_results.csv"
    if not prop_path.exists():
        prop_path = MODULE3_ROOT / sample_name / f"{sample_name}_cell_prop_global_results.csv"
    prop_df = pd.read_csv(prop_path, index_col=0)

    # Initialize ViT
    extractor = ViTFeatureExtractor(model_name=vit_model, pretrained=True, device=device)
    embed_dim = extractor.embed_dim
    logger.info("  ViT model: %s, embed_dim: %d", vit_model, embed_dim)

    # Process spot by spot
    embeddings_dir = output_dir / "embeddings"
    embeddings_dir.mkdir(exist_ok=True)
    n_spots_processed = 0

    for barcode, group in nucleus_spot_map.groupby("spot_barcode"):
        if barcode not in prop_df.index:
            continue

        spot_dir = embeddings_dir / barcode
        spot_dir.mkdir(exist_ok=True)

        # Skip if already extracted
        if (spot_dir / "embeddings.npy").exists():
            n_spots_processed += 1
            continue

        patches = []
        valid_ids = []

        if use_centroid_mode:
            for _, row in group.iterrows():
                try:
                    patch = extract_centroid_patch(
                        image,
                        centroid_x=row["centroid_x"],
                        centroid_y=row["centroid_y"],
                        output_size=patch_size,
                        pixel_size_um=pixel_size_um,
                    )
                    patches.append(normalize_he_patch(patch))
                    valid_ids.append(row["nucleus_id"])
                except (ValueError, IndexError):
                    continue
        else:
            nucleus_ids = group["nucleus_id"].values
            for nid in nucleus_ids:
                try:
                    patch = extract_nucleus_patch(
                        image, mask, nid,
                        output_size=patch_size, expansion=expansion,
                    )
                    patches.append(normalize_he_patch(patch))
                    valid_ids.append(nid)
                except (ValueError, IndexError):
                    continue

        if len(patches) == 0:
            continue

        patches_arr = np.stack(patches)  # (N, 3, H, W)
        embeds = extractor.extract_numpy(patches_arr, batch_size=batch_size)

        # Save per-spot
        np.save(spot_dir / "embeddings.npy", embeds)
        np.save(spot_dir / "nucleus_ids.npy", np.array(valid_ids))

        # Save Module 3 proportions for this spot (filter to cell type columns only)
        available_types = [ct for ct in CELL_TYPES if ct in prop_df.columns]
        props = prop_df.loc[barcode, available_types].values.astype(np.float32)
        np.save(spot_dir / "proportions.npy", props)

        n_spots_processed += 1

    logger.info("  Processed %d spots", n_spots_processed)
    return {"n_spots": n_spots_processed, "embed_dim": embed_dim}


# ===== STAGE 5: Bayesian Cell Assignment + GEX ==========================

def stage5_assign_and_gex(
    sample_name: str,
    output_dir: Path,
    device: str = "cpu",
    embeddings_dir_override: Path = None,
) -> Dict:
    """Bayesian cell assignment + SACE GEX deconvolution.

    Loads pre-computed cuOPT proportions from module3_unified instead of
    re-running the GPU QP solve. Uses morphology embeddings from prior Stage 2
    as the soft prior for bayesian assignment.
    """
    import squidpy as sq
    import sys as _sys
    _sys.path.insert(0, str(REPO_ROOT))
    from CITEgeist.model.citegeist_model import CitegeistModel

    logger.info("Stage 5: Bayesian cell assignment + GEX for %s", sample_name)

    # --- 1. Load pre-computed cuOPT proportions ---
    prop_path = MODULE3_ROOT / f"{sample_name}_cell_prop_finetuned_results.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Module 3 proportions not found: {prop_path}")
    prop_df = pd.read_csv(prop_path, index_col=0)
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

    # --- 3. Assemble morphology embeddings from per-spot dirs ---
    embeddings_dir = embeddings_dir_override if embeddings_dir_override else output_dir / "embeddings"
    emb_list = []
    id_list = []
    emb_id_to_row = {}  # populated below if embeddings found
    for barcode in prop_df.index:
        spot_dir = embeddings_dir / barcode
        emb_file = spot_dir / "embeddings.npy"
        nid_file = spot_dir / "nucleus_ids.npy"
        if not emb_file.exists() or not nid_file.exists():
            continue
        emb = np.load(emb_file)      # (N_cells_in_spot, D)
        nids = np.load(nid_file, allow_pickle=True).astype(str)
        emb_list.append(emb)
        id_list.extend(nids.tolist())

    if emb_list:
        all_embeddings = np.concatenate(emb_list, axis=0)
        # Reorder to match cell_ids order
        emb_id_to_row = {nid: i for i, nid in enumerate(id_list)}
        emb_aligned = np.zeros((len(cell_ids), all_embeddings.shape[1]), dtype=np.float32)
        has_emb = np.zeros(len(cell_ids), dtype=bool)
        for i, cid in enumerate(cell_ids):
            if cid in emb_id_to_row:
                emb_aligned[i] = all_embeddings[emb_id_to_row[cid]]
                has_emb[i] = True
        morphology_embeddings = emb_aligned
        logger.info("  Assembled embeddings: %d/%d cells have embeddings",
                    has_emb.sum(), len(cell_ids))
    else:
        logger.warning("  No embeddings found in %s — morphology prior disabled", embeddings_dir)
        morphology_embeddings = None

    # --- 4. Build detection mask (full-ones fallback — no GPU QP re-run needed) ---
    detection_mask = np.ones((n_spots, len(cell_types)), dtype=bool)
    logger.info("  Using full-ones detection mask (fallback — no re-deconvolution)")

    # --- 5. Instantiate CitegeistModel and inject results ---
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    adata_raw = sq.read.visium(
        str(patient_dir),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=False,
        gex_only=False,
    )
    meta = SAMPLE_METADATA[sample_name]
    model = CitegeistModel(
        sample_name=sample_name,
        adata=adata_raw,
        output_folder=str(output_dir),
    )
    model.split_adata()

    # Inject spatial coordinates — load_images=False omits obsm["spatial"]
    if "spatial" not in model.gene_expression_adata.obsm:
        pos_path = patient_dir / "spatial" / "tissue_positions.csv"
        pos_df = pd.read_csv(pos_path, index_col=0)
        common_bcs = model.gene_expression_adata.obs_names
        coords = pos_df.loc[
            pos_df.index.intersection(common_bcs),
            ["pxl_col_in_fullres", "pxl_row_in_fullres"]
        ].reindex(common_bcs).values
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

    # Realign morphology embeddings to rebuilt cell_ids
    if morphology_embeddings is not None:
        new_emb = np.zeros((len(cell_ids), morphology_embeddings.shape[1]), dtype=np.float32)
        for i, cid in enumerate(cell_ids):
            if cid in emb_id_to_row:
                new_emb[i] = morphology_embeddings[emb_id_to_row[cid]]
        morphology_embeddings = new_emb
        logger.info("  Realigned embeddings after spot subset: %d cells", len(cell_ids))

    # Inject pre-computed results
    model.results["cell_prop"] = prop_df
    model.results["detection_mask_bool"] = detection_mask

    # --- 6. Bayesian assignment ---
    _assignment_method = "bayesian" if morphology_embeddings is not None else "hungarian"
    if _assignment_method == "hungarian":
        logger.warning("  No embeddings — falling back to hungarian assignment")
    assignments_df = model.assign_cells(
        nuclei_counts=nuclei_counts.reindex(prop_df.index, fill_value=0),
        cell_to_spot=cell_to_spot,
        cell_ids=cell_ids,
        morphology_embeddings=morphology_embeddings,
        morphology_weight=0.5 if morphology_embeddings is not None else 0.0,
        assignment_method=_assignment_method,
        detection_mask=detection_mask,
        device=device,
    )

    assignments_df.to_csv(output_dir / "cell_assignments.csv", index=False)
    logger.info("  Assigned %d nuclei across %d spots",
                len(assignments_df), assignments_df["spot_id"].nunique())

    # --- 7. SACE per-cell GEX allocation ---
    cell_assignments_dict = dict(zip(
        assignments_df["cell_id"].astype(str),
        assignments_df["assigned_type"],
    ))

    cell_spot_map = nuc_map[nuc_map["nucleus_id"].astype(str).isin(cell_assignments_dict)].copy()
    cell_spot_map = cell_spot_map.rename(columns={
        "nucleus_id": "cell_id",
        "centroid_x": "x",
        "centroid_y": "y",
    })
    cell_spot_map["cell_id"] = cell_spot_map["cell_id"].astype(str)

    try:
        _, cell_adata, _ = model.run_sace_allocation(
            cell_assignments=cell_assignments_dict,
            cell_spot_map=cell_spot_map,
            max_iter=1,
        )
        sc_path = output_dir / f"{sample_name}_single_cell.h5ad"
        cell_adata.write_h5ad(str(sc_path))
        logger.info("  Saved single-cell AnnData: %d cells x %d genes → %s",
                    cell_adata.shape[0], cell_adata.shape[1], sc_path)
    except Exception as e:
        logger.error("  SACE GEX allocation failed: %s", e, exc_info=True)

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
    """Generate spatial overlays and attention heatmaps."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

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
    nucleus_spot_map = pd.read_csv(output_dir / "nucleus_spot_mapping.csv")

    # Merge assignments with coordinates
    merged = assignments_df.merge(nucleus_spot_map, on="nucleus_id", how="inner")

    # Color palette for cell types
    type_colors = {
        "Endothelial": "#1f77b4",
        "Fibroblasts": "#ff7f0e",
        "B_Cells": "#2ca02c",
        "Macrophages": "#d62728",
        "Monocytes": "#9467bd",
        "CD8_T_Cells": "#8c564b",
        "CD4_T_Cells": "#e377c2",
        "Cancer_Luminal": "#7f7f7f",
        "Cancer_Basal": "#bcbd22",
        "Dendritic_Cells": "#17becf",
    }

    # --- Spatial Overlay ---
    fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150)
    ax.imshow(he_image)
    for _, row in merged.iterrows():
        color = type_colors.get(row["assigned_type"], "#000000")
        ax.scatter(row["centroid_x"], row["centroid_y"],
                   c=color, s=2, alpha=0.6, linewidths=0)
    # Legend
    handles = [mpatches.Patch(color=c, label=t) for t, c in type_colors.items()
               if t in merged["assigned_type"].values]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.8)
    ax.set_title(f"{sample_name} — Cell Type Assignments", fontsize=12)
    ax.axis("off")
    fig.savefig(output_dir / "spatial_overlay.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("  Saved spatial_overlay.png")

    # --- Attention Heatmaps (per cell type) ---
    heatmap_dir = output_dir / "attention_heatmaps"
    heatmap_dir.mkdir(exist_ok=True)

    attentions = np.load(output_dir / "attention_weights.npz", allow_pickle=True)
    embeddings_dir = output_dir / "embeddings"

    for type_idx, ct_name in enumerate(CELL_TYPES):
        fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150)
        ax.imshow(he_image)

        xs, ys, ws = [], [], []
        for barcode in attentions.files:
            att = attentions[barcode]  # (N, K)
            spot_dir = embeddings_dir / barcode
            if not (spot_dir / "nucleus_ids.npy").exists():
                continue
            nids = np.load(spot_dir / "nucleus_ids.npy")
            spot_nuclei = nucleus_spot_map[nucleus_spot_map["nucleus_id"].isin(nids)]

            for j, nid in enumerate(nids):
                row = spot_nuclei[spot_nuclei["nucleus_id"] == nid]
                if len(row) == 0 or j >= att.shape[0]:
                    continue
                if type_idx >= att.shape[1]:
                    continue
                xs.append(row.iloc[0]["centroid_x"])
                ys.append(row.iloc[0]["centroid_y"])
                ws.append(att[j, type_idx])

        if len(xs) > 0:
            sc_plot = ax.scatter(xs, ys, c=ws, cmap="hot", s=2, alpha=0.7,
                                 vmin=0, vmax=np.percentile(ws, 95), linewidths=0)
            plt.colorbar(sc_plot, ax=ax, shrink=0.6, label="Attention weight")

        ax.set_title(f"{sample_name} — {ct_name} attention", fontsize=12)
        ax.axis("off")
        fig.savefig(heatmap_dir / f"{ct_name}_attention.png", bbox_inches="tight", dpi=150)
        plt.close(fig)

    logger.info("  Saved %d attention heatmaps", len(CELL_TYPES))

    return {"visualizations": ["spatial_overlay.png"] + [f"{ct}_attention.png" for ct in CELL_TYPES]}


# ===== MAIN =============================================================

def main():
    parser = argparse.ArgumentParser(description="H&E Morphology Assignment Pipeline")
    parser.add_argument("--sample", required=True, help="Sample name (e.g., HCC22-088-P1-S1)")
    parser.add_argument("--stages", default="1,2,5,6",
                        help="Comma-separated stages to run (default: 1,2,5,6)")
    parser.add_argument("--epochs", type=int, default=100, help="MIL training epochs")
    parser.add_argument("--lr", type=float, default=1e-3, help="MIL learning rate")
    parser.add_argument("--vit-model", default="vit_small_patch16_224", help="ViT model name")
    parser.add_argument("--device", default="cuda", help="Device for GPU stages")
    parser.add_argument("--gpu", action="store_true", default=True, help="Use GPU for StarDist/ViT")
    parser.add_argument("--no-gpu", action="store_true", help="Disable GPU")
    parser.add_argument("--batch-size", type=int, default=256, help="ViT batch size")
    parser.add_argument("--resolution-mode", default="hires",
                        choices=["lowres", "hires", "fullres"],
                        help="Image resolution for StarDist segmentation (default: hires)")
    parser.add_argument("--pooled-checkpoint", type=str, default=None,
                        help="Path to pooled MIL checkpoint (skips per-sample Stage 3)")
    parser.add_argument("--module3-dir", type=str, default=None,
                        help="Override Module 3 proportion directory "
                             "(default: output/module3_unified)")
    parser.add_argument("--embeddings-dir", type=str, default=None,
                        help="Read embeddings from this dir instead of output_dir "
                             "(for reusing Stage 1-2 outputs)")
    parser.add_argument("--output-dir-override", type=str, default=None,
                        help="Override output directory (default: output/morphology_assignment/{sample})")
    args = parser.parse_args()

    if args.no_gpu:
        args.device = "cpu"
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

    embeddings_source = Path(args.embeddings_dir) if args.embeddings_dir else output_dir

    results = {}

    if 1 in stages:
        results["stage1"] = stage1_segment(sample, output_dir, gpu=args.gpu,
                                              resolution_mode=args.resolution_mode)
    if 2 in stages:
        results["stage2"] = stage2_extract_features(
            sample, output_dir, vit_model=args.vit_model,
            batch_size=args.batch_size, device=args.device,
        )
    if 5 in stages:
        results["stage5"] = stage5_assign_and_gex(
            sample, output_dir, device=args.device,
            embeddings_dir_override=embeddings_source / "embeddings" if args.embeddings_dir else None,
        )
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
