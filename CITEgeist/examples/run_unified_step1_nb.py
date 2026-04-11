#!/usr/bin/env python
"""Step 1 NB: Run NB emission deconvolution on patient samples.

Replaces QP-only Module 3 with NB emission model for improved proportions
and downstream GEX deconvolution.

Pipeline:
1. Load patient Visium data (sq.read.visium, gex_only=False)
2. CitegeistModel: split_adata, preprocess_gex, preprocess_antibody
3. QP warm start (run_cell_proportion_model)
4. Extract raw protein counts, build NB marker config
5. Run NB optimization with QP warm-start proportions
6. Inject NB proportions into model, run GEX deconvolution
7. Save NB proportions CSV + GEX parquet + marker file

Usage:
    python run_unified_step1_nb.py --sample HCC22-088-P1-S1
    python run_unified_step1_nb.py --sample HCC22-088-P1-S1 --output-dir /path/to/output
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Absolute paths — no dirname tricks (compute nodes resolve to /var/spool/)
REPO_ROOT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

import squidpy as sq
from model import CitegeistModel
from model.unified_config import CELL_PROFILES_NESTED

from CITEgeist.model.nb_deconvolution import run_nb_optimization
from CITEgeist.model.nb_initialization import (
    initialize_background,
    initialize_dispersion,
    initialize_lambda,
)

# ---------------------------------------------------------------------------
# Patient data path
# ---------------------------------------------------------------------------
DATA_DIR = Path("/ix1/alee/LO_LAB/General/Lab_Data/" "20250210_CITEGeistPublicData_GEO_Alex/processed_files")

# ---------------------------------------------------------------------------
# Patient marker-to-type config for NB model
# ---------------------------------------------------------------------------
# 9 cell types from unified_config, canonical order.
# Marker names are AFTER _strip_suffix (no "-1").
# Each entry: marker_canonical -> list of (type_name, strength)
PATIENT_CELL_TYPES = [
    "Cancer",
    "Macrophages",
    "CD8_T_Cells",
    "CD4_T_Cells",
    "B_Cells",
    "Endothelial",
    "Fibroblasts",
    "Monocytes",
    "Dendritic_Cells",
]

PATIENT_MARKER_TYPE_TABLE = {
    "EPCAM": [("Cancer", "strong")],
    "CD68": [("Macrophages", "strong")],
    "CD163": [("Macrophages", "strong")],
    "CD3E": [("CD8_T_Cells", "strong"), ("CD4_T_Cells", "strong")],
    "CD8A": [("CD8_T_Cells", "strong")],
    "CD4": [("CD4_T_Cells", "strong")],
    "CD19": [("B_Cells", "strong")],
    "PECAM1": [("Endothelial", "strong")],
    "ACTA2": [("Fibroblasts", "strong")],
    "CD14": [("Monocytes", "strong")],
    "ITGAX": [("Dendritic_Cells", "strong"), ("Macrophages", "soft")],
    "HLA-DRA": [("Dendritic_Cells", "strong"), ("B_Cells", "soft"), ("Macrophages", "soft")],
}


def _strip_suffix(name):
    """Strip -1 suffix from patient marker names (e.g., CD68-1 -> CD68).

    Only strips when the base name is in our patient marker table.
    """
    if name.endswith("-1"):
        base = name[:-2]
        if base in PATIENT_MARKER_TYPE_TABLE:
            return base
    return name


def build_patient_marker_config(available_markers):
    """Build NB marker config for patient antibody panel.

    Returns:
        markers: list of canonical marker names found in data.
        active_mask: (T, M) bool array.
        type_names: list of cell type names.
        raw_indices: list of int indices into available_markers for each marker.
    """
    type_names = list(PATIENT_CELL_TYPES)
    T = len(type_names)

    markers = []
    raw_indices = []
    seen = set()

    for idx, raw_name in enumerate(available_markers):
        canonical = _strip_suffix(raw_name)
        if canonical in PATIENT_MARKER_TYPE_TABLE and canonical not in seen:
            markers.append(canonical)
            raw_indices.append(idx)
            seen.add(canonical)

    M = len(markers)
    active_mask = np.zeros((T, M), dtype=bool)

    for m_idx, marker in enumerate(markers):
        for type_name, _strength in PATIENT_MARKER_TYPE_TABLE[marker]:
            if type_name in type_names:
                t_idx = type_names.index(type_name)
                active_mask[t_idx, m_idx] = True

    return markers, active_mask, type_names, raw_indices


def run_step1_nb(sample_name, output_dir=None):
    """Run NB Module 3 for a single patient sample."""

    if output_dir is None:
        output_dir = REPO_ROOT / "output" / "module3_nb"
    else:
        output_dir = Path(output_dir)

    sample_output = output_dir / sample_name
    sample_output.mkdir(parents=True, exist_ok=True)

    marker_file = sample_output / ".step1_nb_complete"
    if marker_file.exists():
        logger.info("Step 1 NB already complete for %s, skipping", sample_name)
        return

    # ------------------------------------------------------------------
    # 1. Load patient Visium data
    # ------------------------------------------------------------------
    sample_path = DATA_DIR / sample_name / "outs"
    if not sample_path.exists():
        raise FileNotFoundError(f"Sample path not found: {sample_path}")

    logger.info("Loading patient data from %s", sample_path)
    adata = sq.read.visium(
        str(sample_path),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=True,
        gex_only=False,
    )

    model = CitegeistModel(
        sample_name=sample_name,
        adata=adata,
        output_folder=str(sample_output),
    )
    model.split_adata()

    # Filter spots with NaN spatial coordinates
    for adata_attr in ("gene_expression_adata", "antibody_capture_adata"):
        ad = getattr(model, adata_attr, None)
        if ad is not None and "spatial" in ad.obsm:
            finite_mask = np.isfinite(ad.obsm["spatial"]).all(axis=1)
            n_nan = (~finite_mask).sum()
            if n_nan > 0:
                logger.warning(
                    "Filtering %d spots with NaN spatial coords from %s",
                    n_nan,
                    adata_attr,
                )
                setattr(model, adata_attr, ad[finite_mask].copy())

    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(CELL_PROFILES_NESTED)

    # ------------------------------------------------------------------
    # 2. QP warm start
    # ------------------------------------------------------------------
    logger.info("--- Stage 1: QP warm start ---")
    t0 = time.time()
    model.run_cell_proportion_model(
        validation_warn_only=True,
        skip_finetuning=True,
        lambda_laplacian=0.0,  # No spatial smoothing — let NB handle its own sharpness
    )
    qp_time = time.time() - t0
    logger.info("QP warm start done in %.1fs", qp_time)

    qp_props = model.results["cell_prop"]
    logger.info(
        "QP proportions: %d spots x %d types",
        *qp_props.shape,
    )

    # ------------------------------------------------------------------
    # 3. Build NB marker config from patient antibody panel
    # ------------------------------------------------------------------
    available_markers = list(model.antibody_capture_adata.var_names)
    markers, active_mask, type_names, raw_indices = build_patient_marker_config(
        available_markers,
    )
    T = len(type_names)
    M = len(markers)
    logger.info(
        "NB config: %d types x %d markers, %d active pairs",
        T,
        M,
        int(active_mask.sum()),
    )
    logger.info("Markers found: %s", markers)

    if M == 0:
        raise RuntimeError("No NB markers found in antibody panel. " f"Available: {available_markers}")

    # ------------------------------------------------------------------
    # 4. Extract raw counts for NB markers
    # ------------------------------------------------------------------
    raw = model.antibody_capture_adata.layers.get("raw_counts", None)
    if raw is None:
        raise RuntimeError("No raw_counts layer — preprocess_antibody() must store it")
    raw = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)

    S = raw[:, raw_indices].astype(np.float64)
    S = np.maximum(S, 0.0)
    I = S.shape[0]  # noqa: E741 — number of spots (NB model notation)

    # ------------------------------------------------------------------
    # 5. Prepare cellularity (N=1, no nuclei counts for patient data)
    #    and QP warm-start proportions
    # ------------------------------------------------------------------
    N = np.ones(I, dtype=np.float64)

    # Map QP 9-type proportions to NB type order
    p_init = np.zeros((I, T), dtype=np.float64)
    for t_idx, ct in enumerate(type_names):
        if ct in qp_props.columns:
            p_init[:, t_idx] = qp_props[ct].values

    # Normalize to sum to 1
    row_sums = p_init.sum(axis=1, keepdims=True)
    p_init = np.where(row_sums > 0, p_init / row_sums, 1.0 / T)

    # ------------------------------------------------------------------
    # 6. Initialize NB parameters and run optimization
    # ------------------------------------------------------------------
    logger.info("--- Stage 2: NB emission optimization ---")
    t1 = time.time()

    median_N = float(np.median(N))  # 1.0 for patient data
    lam_init = initialize_lambda(S, active_mask, median_N)
    b_init = initialize_background(S)
    r_init = initialize_dispersion(S)

    logger.info(
        "NB init: %d spots, median_N=%.1f, lam_nz_mean=%.2f, b_mean=%.1f",
        I,
        median_N,
        lam_init[active_mask].mean(),
        b_init.mean(),
    )

    result = run_nb_optimization(
        S=S,
        N=N,
        active_mask=active_mask,
        n_iter=100,
        lam_init=lam_init,
        r_init=r_init,
        b_init=b_init,
        p_init=p_init,
        # Sprint hyperparams
        lam_prior_sigma=2.0,
        lr_global=0.02,
        lr_local=0.05,
        patience=3,
        n_starts=1,  # QP warm start only — multi-start selects smoother solutions
        device="cpu",
    )

    nb_time = time.time() - t1
    n_iters_run = len(result["held_out_nll"])
    logger.info(
        "NB optimization done in %.1fs (%d iterations)",
        nb_time,
        n_iters_run,
    )

    # ------------------------------------------------------------------
    # 7. Build NB proportions DataFrame (9-type, same columns as QP)
    # ------------------------------------------------------------------
    spot_names = model.antibody_capture_adata.obs_names
    nb_props = pd.DataFrame(
        result["proportions"],
        index=spot_names,
        columns=type_names,
    )

    # Ensure all QP columns present (should already match)
    for ct in qp_props.columns:
        if ct not in nb_props.columns:
            nb_props[ct] = 0.0
    nb_props = nb_props[list(qp_props.columns)]

    # Renormalize
    row_sums = nb_props.sum(axis=1)
    nb_props = nb_props.div(row_sums.replace(0, 1), axis=0)

    # Save NB proportions
    props_path = sample_output / f"{sample_name}_cell_prop_nb_results.csv"
    nb_props.to_csv(props_path)
    logger.info("Saved NB proportions: %s (%d x %d)", props_path, *nb_props.shape)

    # ------------------------------------------------------------------
    # 8. Inject NB proportions and run GEX deconvolution
    # ------------------------------------------------------------------
    logger.info("--- Stage 3: GEX deconvolution with NB proportions ---")
    t2 = time.time()

    model.results["cell_prop"] = nb_props
    model.run_cell_expression_pass1()

    gex_time = time.time() - t2
    logger.info("GEX deconvolution done in %.1fs", gex_time)

    # Save GEX results
    gex_path = sample_output / f"{sample_name}_gene_expression_nb.parquet"
    if "gene_expression" in model.results and model.results["gene_expression"] is not None:
        model.results["gene_expression"].to_parquet(gex_path)
        logger.info("Saved GEX: %s", gex_path)
    else:
        logger.warning("No gene_expression result to save")

    # Save results h5ad for downstream steps
    try:
        result_adata = model.get_results_adata()
        h5ad_path = sample_output / f"{sample_name}_module3_nb_results.h5ad"
        result_adata.write_h5ad(str(h5ad_path))
        logger.info("Saved results h5ad: %s", h5ad_path)
    except Exception as e:
        logger.warning("Could not save h5ad: %s", e)

    # ------------------------------------------------------------------
    # 9. Write completion marker
    # ------------------------------------------------------------------
    marker_file.touch()
    logger.info(
        "Step 1 NB complete for %s (QP=%.1fs, NB=%.1fs, GEX=%.1fs)",
        sample_name,
        qp_time,
        nb_time,
        gex_time,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Unified pipeline Step 1 NB: Module 3 with NB emission model",
    )
    parser.add_argument("--sample", required=True, help="Patient sample name")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: output/module3_nb/)",
    )
    args = parser.parse_args()
    run_step1_nb(args.sample, args.output_dir)
