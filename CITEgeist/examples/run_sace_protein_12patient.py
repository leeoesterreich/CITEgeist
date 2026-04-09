#!/usr/bin/env python
"""Run SACE per-cell protein annotation on all 12 HCC22-088 patient samples.

For each sample:
  1. Reconstruct CitegeistModel from saved outputs (no re-running QP)
  2. Call run_sace_protein() using saved cell_assignments and nucleus_spot_map
  3. Update {sample}_single_cell.h5ad with obs columns:
       func_{marker}_{type}_gate   (binary 0/1)
       func_{marker}_{type}_score  (continuous SACE-allocated level)
  4. Write per-sample protein_gates.csv for downstream analysis

Usage:
    # Single sample (test before submitting array):
    python run_sace_protein_12patient.py --sample HCC22-088-P1-S1

    # All samples via SLURM array:
    sbatch CITEgeist/slurm/sbatch_sace_protein_12patient.sh
"""

import argparse
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import json

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

from model import CitegeistModel
from model.unified_config import CELL_PROFILES_NESTED
from model.annotation.functional_annotation import DEFAULT_FUNCTIONAL_TABLE


REPO = Path(__file__).resolve().parents[2]

SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

MORPH_DIR = REPO / "output" / "morphology_assignment_v3"
QP_DIR = REPO / "output" / "module3_cuopt_qp"
M1_DIR = REPO / "output" / "module12_discovery"
RAW_DATA = Path(
    "/ix1/alee/LO_LAB/General/Lab_Data/"
    "20250210_CITEGeistPublicData_GEO_Alex/processed_files"
)
OUTPUT_DIR = REPO / "output" / "sace_protein_12patient"

# Cancer type aliases used in old QP outputs — merged to Epithelial in v2
_CANCER_TYPE_RENAME = {"Cancer_Basal": "Epithelial", "Cancer_Luminal": "Epithelial"}


def _load_m1_functional_table(sample_name: str) -> dict:
    """Return DEFAULT_FUNCTIONAL_TABLE filtered to M1-interesting markers for this sample.

    Loads the module12 discovery JSON for sample_name and keeps only entries
    whose marker key appears in M1's interesting_markers list.  Falls back to
    the full DEFAULT_FUNCTIONAL_TABLE if the JSON is missing.
    """
    discovery_path = M1_DIR / f"{sample_name}_module12_discovery.json"
    if not discovery_path.exists():
        logger.warning(
            "M1 discovery JSON not found for %s; using full DEFAULT_FUNCTIONAL_TABLE",
            sample_name,
        )
        return DEFAULT_FUNCTIONAL_TABLE

    with open(discovery_path) as f:
        disc = json.load(f)

    interesting = set(disc.get("module1", {}).get("interesting_markers", []))
    filtered = {k: v for k, v in DEFAULT_FUNCTIONAL_TABLE.items() if k in interesting}

    dropped = set(DEFAULT_FUNCTIONAL_TABLE) - set(filtered)
    logger.info(
        "%s: M1 filter kept %d/%d functional markers (dropped: %s)",
        sample_name, len(filtered), len(DEFAULT_FUNCTIONAL_TABLE),
        ", ".join(sorted(dropped)) if dropped else "none",
    )
    return filtered

# min_counts mirrors Phase 2 biopsy/surgical split
_BIOPSY_KEYWORDS = ("-S1",)


def _is_biopsy(sample_name: str) -> bool:
    return any(kw in sample_name for kw in _BIOPSY_KEYWORDS)


def run_sace_protein_sample(sample_name: str) -> None:
    """Run SACE protein annotation for a single sample.

    Args:
        sample_name: Sample identifier, e.g. 'HCC22-088-P1-S1'.
    """
    sample_out = OUTPUT_DIR / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    marker_file = sample_out / ".sace_v2_complete"
    if marker_file.exists():
        logger.info("SACE protein already complete for %s, skipping", sample_name)
        return

    # ------------------------------------------------------------------
    # 1. Load Visium data
    # ------------------------------------------------------------------
    visium_dir = RAW_DATA / sample_name / "outs"
    if not visium_dir.exists():
        raise FileNotFoundError(f"SpaceRanger outs not found: {visium_dir}")

    logger.info("Loading SpaceRanger data for %s", sample_name)
    adata = sq.read.visium(
        str(visium_dir),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=False,
        gex_only=False,
    )

    # sq.read.visium with load_images=False does not populate obsm["spatial"].
    # Manually inject spatial coordinates from tissue_positions.csv.
    pos_file = visium_dir / "spatial" / "tissue_positions.csv"
    if not pos_file.exists():
        pos_file = visium_dir / "spatial" / "tissue_positions_list.csv"
    if pos_file.exists():
        pos_df = pd.read_csv(pos_file, index_col=0)
        # SpaceRanger >= 2.0 uses pxl_row_in_fullres / pxl_col_in_fullres
        if "pxl_col_in_fullres" in pos_df.columns and "pxl_row_in_fullres" in pos_df.columns:
            coords = pos_df.loc[adata.obs_names, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
        else:
            # SpaceRanger < 2.0: columns are unnamed (0=barcode,1=in_tissue,2=row,3=col,4=y,5=x)
            coords = pos_df.loc[adata.obs_names].iloc[:, [4, 3]].values
        import numpy as _np
        adata.obsm["spatial"] = _np.array(coords, dtype=float)
        logger.info("Injected spatial coords from %s", pos_file.name)
    else:
        logger.warning("tissue_positions file not found for %s; spatial coords unavailable", sample_name)

    model = CitegeistModel(
        sample_name=sample_name,
        adata=adata,
        output_folder=str(sample_out),
    )
    model.split_adata()

    # Filter NaN spatial coordinates
    for attr in ("gene_expression_adata", "antibody_capture_adata"):
        ad = getattr(model, attr, None)
        if ad is not None and "spatial" in ad.obsm:
            finite_mask = np.isfinite(ad.obsm["spatial"]).all(axis=1)
            n_nan = int((~finite_mask).sum())
            if n_nan > 0:
                logger.warning(
                    "Filtering %d spots with NaN spatial coords from %s", n_nan, attr
                )
                setattr(model, attr, ad[finite_mask].copy())

    min_counts = 100 if _is_biopsy(sample_name) else 25
    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(CELL_PROFILES_NESTED)

    # ------------------------------------------------------------------
    # 2. Inject saved proportions from cuOPT QP re-run (Epithelial profile)
    # ------------------------------------------------------------------
    props_path = QP_DIR / sample_name / f"{sample_name}_cell_prop_global_results.csv"
    if not props_path.exists():
        raise FileNotFoundError(f"QP proportions not found: {props_path}")

    props_df = pd.read_csv(props_path, index_col=0)
    # Drop solver diagnostic column if present
    props_df = props_df.drop(columns=["recon_error"], errors="ignore")
    # Align to GEX spot order (filter_gex may have dropped some spots)
    gex_barcodes = model.gene_expression_adata.obs_names
    common = gex_barcodes.intersection(props_df.index)
    if len(common) < len(gex_barcodes):
        logger.warning(
            "%s: %d GEX spots not in saved proportions; subsetting",
            sample_name, len(gex_barcodes) - len(common),
        )
        model.gene_expression_adata = model.gene_expression_adata[common].copy()
        model.antibody_capture_adata = model.antibody_capture_adata[common].copy()
    model.results["cell_prop"] = props_df.loc[model.gene_expression_adata.obs_names]
    logger.info(
        "Injected proportions: %d spots x %d types",
        *model.results["cell_prop"].shape,
    )

    # ------------------------------------------------------------------
    # 3. Build cell_assignments dict from saved file (rename Cancer → Epithelial)
    # ------------------------------------------------------------------
    assignments_path = MORPH_DIR / sample_name / "cell_assignments.csv"
    if not assignments_path.exists():
        raise FileNotFoundError(f"cell_assignments.csv not found: {assignments_path}")

    assignments_df = pd.read_csv(assignments_path)
    assignments_df["assigned_type"] = assignments_df["assigned_type"].replace(_CANCER_TYPE_RENAME)
    cell_assignments = dict(
        zip(
            assignments_df["nucleus_id"].astype(str),
            assignments_df["assigned_type"],
        )
    )
    n_epithelial = sum(1 for v in cell_assignments.values() if v == "Epithelial")
    logger.info("Loaded %d cell assignments (%d Epithelial)", len(cell_assignments), n_epithelial)

    # ------------------------------------------------------------------
    # 4. Build cell_spot_map DataFrame from nucleus_spot_map
    # ------------------------------------------------------------------
    spot_map_path = MORPH_DIR / sample_name / f"{sample_name}_nucleus_spot_map.csv"
    if not spot_map_path.exists():
        raise FileNotFoundError(f"nucleus_spot_map not found: {spot_map_path}")

    spot_map_df = pd.read_csv(spot_map_path)
    # run_sace_protein() accepts nucleus_id / spot_id / x_pixel / y_pixel and
    # handles the rename internally — pass as-is.
    cell_spot_map = spot_map_df.copy()
    cell_spot_map["nucleus_id"] = cell_spot_map["nucleus_id"].astype(str)

    # ------------------------------------------------------------------
    # 5a. Run SACE GEX (allocates spot GEX to per-cell expression)
    # ------------------------------------------------------------------
    logger.info("Running SACE GEX for %s", sample_name)
    model.run_cell_expression_pass1(
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
    )
    gex_cell_adata = model.results.get("sace_cell_adata")
    if gex_cell_adata is not None:
        logger.info(
            "SACE GEX complete: %d cells x %d genes",
            *gex_cell_adata.shape,
        )
    else:
        logger.warning("SACE GEX returned no cell_adata for %s", sample_name)

    # ------------------------------------------------------------------
    # 5b. Run SACE protein (M1-filtered functional table)
    # ------------------------------------------------------------------
    m1_functional_table = _load_m1_functional_table(sample_name)
    logger.info("Running run_sace_protein() for %s", sample_name)
    result = model.run_sace_protein(
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        functional_table=m1_functional_table,
        max_iter=1,
        min_high_component_log_mean=1.0,
    )

    if not result:
        logger.warning("run_sace_protein() returned empty result for %s", sample_name)
        marker_file.touch()
        return

    # ------------------------------------------------------------------
    # 6. Rebuild single_cell.h5ad: updated GEX (X), cell_type, protein gates
    # ------------------------------------------------------------------
    sc_path = MORPH_DIR / sample_name / f"{sample_name}_single_cell.h5ad"
    if not sc_path.exists():
        raise FileNotFoundError(f"single_cell.h5ad not found: {sc_path}")

    adata_sc = sc.read_h5ad(sc_path)
    gates_df = result["protein_gates_df"]
    cell_protein = result["cell_protein"]       # (N_cells, M) float32 array
    protein_names = result["protein_names"]     # list of marker names (canonical)

    # 6a. Update cell_type: rename Cancer_Basal/Cancer_Luminal → Epithelial
    adata_sc.obs["cell_type"] = adata_sc.obs["cell_type"].replace(_CANCER_TYPE_RENAME)

    # 6b. Replace X with new SACE GEX allocation (Epithelial-aware proportions)
    if gex_cell_adata is not None:
        common_cells = adata_sc.obs_names.intersection(gex_cell_adata.obs_names)
        if len(common_cells) == adata_sc.n_obs:
            # Align genes: intersect to handle any minor var_names differences
            common_genes = adata_sc.var_names.intersection(gex_cell_adata.var_names)
            if len(common_genes) == adata_sc.n_vars:
                import scipy.sparse as _sp
                new_X = gex_cell_adata[common_cells].X
                if _sp.issparse(new_X):
                    new_X = new_X.toarray()
                adata_sc.X = new_X.astype(np.float32)
                logger.info("%s: replaced X with SACE GEX (%d cells × %d genes)",
                            sample_name, len(common_cells), len(common_genes))
            else:
                logger.warning("%s: gene mismatch (old %d, new %d) — X not updated",
                               sample_name, adata_sc.n_vars, len(common_genes))
        else:
            logger.warning("%s: cell count mismatch (h5ad %d, sace %d) — X not updated",
                           sample_name, adata_sc.n_obs, len(common_cells))

    # 6c. Drop stale func_* columns from previous M3.5 run
    stale_func = [c for c in adata_sc.obs.columns if c.startswith("func_")]
    if stale_func:
        adata_sc.obs = adata_sc.obs.drop(columns=stale_func)
        logger.info("%s: dropped %d stale func_* columns", sample_name, len(stale_func))

    # 6d. Add gate columns (func_{marker}_{type}_gate)
    n_added_gate = 0
    for col in gates_df.columns:
        adata_sc.obs[col] = gates_df[col].reindex(adata_sc.obs.index).values
        n_added_gate += 1

    # 6e. Add score columns (func_{marker}_{type}_score)
    protein_df = pd.DataFrame(
        cell_protein,
        index=gates_df.index,
        columns=protein_names,
    )
    n_added_score = 0
    for col in gates_df.columns:
        if not col.endswith("_gate"):
            continue
        inner = col[len("func_"):-len("_gate")]  # "{marker}_{type}"
        matched_marker = None
        for pname in sorted(protein_names, key=len, reverse=True):
            safe_pname = pname.replace("-", "_")
            if inner.startswith(safe_pname + "_") or inner == safe_pname:
                matched_marker = pname
                break
        if matched_marker is None:
            continue
        score_col = col[:-len("_gate")] + "_score"
        adata_sc.obs[score_col] = protein_df[matched_marker].reindex(adata_sc.obs.index).values
        n_added_score += 1

    # 6f. Update protein obsm
    protein_aligned = protein_df.reindex(adata_sc.obs.index)
    adata_sc.obsm["protein"] = protein_aligned.values.astype(np.float32)
    adata_sc.uns["protein_names"] = protein_names

    logger.info(
        "%s: updated h5ad — cell_type Epithelial, %d gate + %d score cols (%d cells)",
        sample_name, n_added_gate, n_added_score, adata_sc.n_obs,
    )

    adata_sc.write_h5ad(sc_path)
    logger.info("Updated single_cell.h5ad written to %s", sc_path)

    # Log GMM summary so we can check PD-L1 gating_method across samples
    gmm_summary = result.get("gmm_summary", [])
    if gmm_summary:
        gmm_df = pd.DataFrame(gmm_summary)
        gmm_path = sample_out / f"{sample_name}_gmm_summary.csv"
        gmm_df.to_csv(gmm_path, index=False)
        logger.info("GMM summary written to %s", gmm_path)
        # Flag PD-L1 Macrophages gating decision for downstream review
        pdl1_mac = gmm_df[
            (gmm_df.get("marker", gmm_df.columns[0]) == "PD-L1")
            & (gmm_df.get("cell_type", "") == "Macrophages")
        ] if {"marker", "cell_type"}.issubset(gmm_df.columns) else pd.DataFrame()
        if not pdl1_mac.empty:
            method = pdl1_mac.iloc[0].get("gating_method", "unknown")
            logger.info("PD-L1 Macrophages gating_method for %s: %s", sample_name, method)
        else:
            logger.info("PD-L1 Macrophages not found in GMM summary for %s", sample_name)

    marker_file.touch()
    logger.info("SACE protein complete for %s", sample_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SACE per-cell protein annotation for 12-patient HCC22-088 cohort"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample", help="Sample name, e.g. HCC22-088-P1-S1")
    group.add_argument(
        "--array-index",
        type=int,
        help="SLURM_ARRAY_TASK_ID (0-11) to select sample from SAMPLES list",
    )
    args = parser.parse_args()

    if args.array_index is not None:
        sample = SAMPLES[args.array_index]
    else:
        sample = args.sample

    logger.info("Processing sample: %s", sample)
    run_sace_protein_sample(sample)
