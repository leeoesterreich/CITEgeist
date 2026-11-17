"""
run_pseudo_geomx.py

Run the CITEgeist pipeline (M1 → M2 → M3 QP → SACE GEX) on pseudo-GeoMx
AOI data produced by create_pseudo_geomx.py.

Input
-----
  aoi_gex.h5ad     : AnnData  (n_aois × n_genes),  obsm["spatial"] required
  aoi_protein.h5ad : AnnData  (n_aois × n_proteins), obsm["spatial"] required

Output (written to --output-dir)
---------------------------------
  predicted_proportions.csv  : AOI × cell-type proportions
  discovered_profiles.json   : M2-selected protein profiles
  sace_cell_gex.h5ad         : per-AOI GEX AnnData from SACE
  pipeline_results.json      : timing, shapes, diagnostics
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd

# ── sys.path ─────────────────────────────────────────────────────────────────
_REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(_REPO_ROOT / "CITEgeist"))

from CITEgeist.model import analyze_marker_colocalization, identify_interesting_markers, select_profiles  # noqa: E402
from CITEgeist.model.citegeist_model import CitegeistModel  # noqa: E402
from CITEgeist.model.discovery.spatial_colocalization import discover_profiles_continuous  # noqa: E402

# ── logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Helper: build 1-to-1 AOI cell-spot map
# ─────────────────────────────────────────────────────────────────────────────


def build_1to1_cell_spot_map(adata: ad.AnnData) -> pd.DataFrame:
    """Build a 1-to-1 cell-spot map where each AOI is both a cell and a spot.

    Parameters
    ----------
    adata : AnnData
        AOI-level AnnData with ``obsm["spatial"]`` (n_aois × 2).

    Returns
    -------
    pd.DataFrame
        Columns: cell_id, spot_barcode, spot_idx, x, y.
    """
    aoi_names = list(adata.obs_names)
    coords = adata.obsm["spatial"]  # (n_aois, 2)
    return pd.DataFrame(
        {
            "cell_id": aoi_names,
            "spot_barcode": aoi_names,
            "spot_idx": list(range(len(aoi_names))),
            "x": coords[:, 0].astype(float),
            "y": coords[:, 1].astype(float),
        }
    )


# ─────────────────────────────────────────────────────────────────────────────
# Helper: M1 + M2
# ─────────────────────────────────────────────────────────────────────────────


def run_m1_m2(
    model: CitegeistModel,
    top_k: int = 3,
    n_permutations: int = 999,
    seed: int = 42,
) -> Dict:
    """Run M1 (marker interest) and M2 (colocalization + profile discovery).

    Parameters
    ----------
    model : CitegeistModel
        Initialised and preprocessed model.
    top_k : int
        Number of top profile candidates for M2b.
    n_permutations : int
        Permutation count for Moran's I and colocalization tests.
    seed : int
        Random seed.

    Returns
    -------
    dict
        ``cell_profile_dict`` ready for ``model.load_cell_profile_dict()``.
    """
    protein_adata = model.antibody_capture_adata
    n_aois = protein_adata.n_obs

    # Raw antibody matrix (dense)
    X_raw = protein_adata.X.toarray() if hasattr(protein_adata.X, "toarray") else np.array(protein_adata.X)
    coords = protein_adata.obsm["spatial"]
    marker_names = list(protein_adata.var_names)

    # Clamp k to avoid "not enough neighbours" errors with sparse AOI sets
    k_safe = min(8, n_aois - 1)
    if k_safe < 2:
        raise ValueError(f"Need at least 3 AOIs for M1/M2; got {n_aois}.")
    if k_safe < 8:
        logger.warning(
            "Only %d AOIs — clamping morans_k and neighbor_k to %d.",
            n_aois,
            k_safe,
        )

    # ── M1: marker interest ──────────────────────────────────────────────────
    logger.info("Running M1 — marker interest detection (%d markers).", len(marker_names))
    m1_result = identify_interesting_markers(
        X=X_raw,
        coords=coords,
        marker_names=marker_names,
        morans_k=k_safe,
        morans_n_perm=min(n_permutations, 99),  # M1 uses a lighter perm count
        seed=seed,
        verbose=True,
    )
    interesting = m1_result.interesting_markers
    logger.info("M1 found %d interesting markers: %s", len(interesting), interesting)

    if len(interesting) == 0:
        logger.warning("M1 found no interesting markers; falling back to all markers.")
        interesting = marker_names

    # ── M2a: colocalization ──────────────────────────────────────────────────
    logger.info("Running M2a — spatial colocalization.")
    coloc_result = analyze_marker_colocalization(
        X=X_raw,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting,
        neighbor_k=k_safe,
        n_permutations=n_permutations,
        seed=seed,
        verbose=True,
    )

    # ── M2b: profile discovery ───────────────────────────────────────────────
    logger.info("Running M2b — profile discovery (top_k=%d).", top_k)
    discovery_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        top_k=top_k,
        distance_metric="colocalization_score",
        seed=seed,
        verbose=True,
    )

    # ── M2c: profile selection ───────────────────────────────────────────────
    logger.info("Running M2c — profile selection.")
    candidate_profiles = discovery_result.profiles  # List[List[str]]
    selection_result = select_profiles(
        X=X_raw,
        coords=coords,
        marker_names=marker_names,
        profiles=candidate_profiles,
        _interesting_markers=interesting,
        _colocalization_result=coloc_result,
        verbose=True,
    )

    # ── Convert to cell_profile_dict ─────────────────────────────────────────
    cell_profile_dict: Dict = {}
    for i, marker_list in enumerate(selection_result.selected_profiles):
        cell_profile_dict[f"Type_{i}"] = {
            "Major": list(marker_list),
            "Soft": [],
        }
    logger.info(
        "M2 selected %d profiles: %s",
        len(cell_profile_dict),
        {k: v["Major"] for k, v in cell_profile_dict.items()},
    )
    return cell_profile_dict


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────


def run_pipeline(
    aoi_gex_path: str,
    aoi_protein_path: str,
    output_dir: str,
    skip_m1m2: bool = False,
    profile_dict_path: Optional[str] = None,
    seed: int = 42,
) -> Dict:
    """Run the full CITEgeist pipeline on pseudo-GeoMx AOI data.

    Parameters
    ----------
    aoi_gex_path : str
        Path to ``aoi_gex.h5ad`` (AOI × genes).
    aoi_protein_path : str
        Path to ``aoi_protein.h5ad`` (AOI × proteins).
    output_dir : str
        Directory where outputs will be written.
    skip_m1m2 : bool
        If True, skip M1/M2 and load profile dict from *profile_dict_path*.
    profile_dict_path : str, optional
        Path to a JSON file containing a pre-computed cell_profile_dict.
        Required when *skip_m1m2* is True.
    seed : int
        Global random seed.

    Returns
    -------
    dict
        Summary of results (shapes, timings, diagnostics).
    """
    t_total_start = time.time()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    results: Dict = {
        "aoi_gex_path": str(aoi_gex_path),
        "aoi_protein_path": str(aoi_protein_path),
        "output_dir": str(output_dir),
        "seed": seed,
        "timing": {},
        "shapes": {},
    }

    # ── Load inputs ───────────────────────────────────────────────────────────
    logger.info("Loading GEX AnnData from %s", aoi_gex_path)
    gex_adata = ad.read_h5ad(aoi_gex_path)
    logger.info("Loading protein AnnData from %s", aoi_protein_path)
    protein_adata = ad.read_h5ad(aoi_protein_path)

    logger.info(
        "Inputs: %d AOIs, %d genes, %d proteins.",
        gex_adata.n_obs,
        gex_adata.n_vars,
        protein_adata.n_vars,
    )
    results["shapes"]["n_aois"] = gex_adata.n_obs
    results["shapes"]["n_genes_raw"] = gex_adata.n_vars
    results["shapes"]["n_proteins"] = protein_adata.n_vars

    # ── Initialise model ──────────────────────────────────────────────────────
    logger.info("Initialising CitegeistModel.")
    model = CitegeistModel(
        sample_name="pseudo_geomx",
        output_folder=str(out),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    # ── Preprocessing ─────────────────────────────────────────────────────────
    logger.info("Preprocessing — filter_gex / preprocess_gex / preprocess_antibody.")
    model.filter_gex(min_counts=10)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    results["shapes"]["n_genes_postfilter"] = model.gene_expression_adata.n_vars

    # ── M1 + M2 (or load pre-computed profiles) ───────────────────────────────
    t_m1m2_start = time.time()
    if skip_m1m2:
        if profile_dict_path is None:
            raise ValueError("--profile-dict is required when --skip-m1m2 is set.")
        logger.info("Skipping M1/M2 — loading profiles from %s", profile_dict_path)
        with open(profile_dict_path) as fh:
            cell_profile_dict = json.load(fh)
    else:
        cell_profile_dict = run_m1_m2(model, seed=seed)

    results["timing"]["m1_m2_s"] = round(time.time() - t_m1m2_start, 2)
    results["discovered_profiles"] = {k: v["Major"] for k, v in cell_profile_dict.items()}

    # Save profiles
    profiles_path = out / "discovered_profiles.json"
    with open(profiles_path, "w") as fh:
        json.dump(cell_profile_dict, fh, indent=2)
    logger.info("Profiles saved to %s", profiles_path)

    # ── Load profiles into model ──────────────────────────────────────────────
    model.load_cell_profile_dict(cell_profile_dict)

    # ── M3 QP deconvolution ───────────────────────────────────────────────────
    logger.info("Running M3 — QP cell-proportion deconvolution.")
    t_m3_start = time.time()
    model.run_cell_proportion_model(
        method="qp",
        lambda_laplacian=0.0,
        validation_warn_only=True,
    )
    results["timing"]["m3_qp_s"] = round(time.time() - t_m3_start, 2)

    prop_df = model.results["cell_prop"]
    logger.info("M3 done — proportions shape: %s", prop_df.shape)
    results["shapes"]["prop_df"] = list(prop_df.shape)

    # Save proportions
    props_path = out / "predicted_proportions.csv"
    prop_df.to_csv(props_path)
    logger.info("Proportions saved to %s", props_path)

    # ── SACE GEX allocation ───────────────────────────────────────────────────
    logger.info("Running SACE GEX allocation.")
    t_sace_start = time.time()

    # Build 1:1 AOI cell-spot map using the (potentially filtered) GEX adata
    cell_spot_map = build_1to1_cell_spot_map(model.gene_expression_adata)

    # Assign each AOI to its dominant cell type
    # Only include AOIs present in prop_df (some may be filtered out)
    common_aois = [name for name in model.gene_expression_adata.obs_names if name in prop_df.index]
    cell_assignments: Dict[str, str] = {name: prop_df.loc[name].idxmax() for name in common_aois}
    # Restrict cell_spot_map to common AOIs
    cell_spot_map = cell_spot_map[cell_spot_map["cell_id"].isin(common_aois)].reset_index(drop=True)

    logger.info(
        "SACE: %d AOIs → %d cell types.",
        len(cell_assignments),
        len(set(cell_assignments.values())),
    )

    spot_type_gex, diagnostics = model.run_sace_allocation(
        output_mode="layers",
    )
    cell_adata = None
    results["timing"]["sace_s"] = round(time.time() - t_sace_start, 2)

    if cell_adata is not None:
        results["shapes"]["sace_cell_adata"] = list(cell_adata.shape)
        sace_path = out / "sace_cell_gex.h5ad"
        cell_adata.write_h5ad(sace_path)
        logger.info("SACE cell AnnData saved to %s", sace_path)
    else:
        logger.warning("SACE returned None for cell_adata — skipping h5ad save.")
        results["shapes"]["sace_cell_adata"] = None

    if diagnostics:
        results["sace_diagnostics"] = diagnostics

    # ── Finalise ──────────────────────────────────────────────────────────────
    results["timing"]["total_s"] = round(time.time() - t_total_start, 2)
    logger.info("Pipeline complete in %.1fs.", results["timing"]["total_s"])

    pipeline_results_path = out / "pipeline_results.json"
    with open(pipeline_results_path, "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    logger.info("Pipeline results saved to %s", pipeline_results_path)

    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run CITEgeist M1→M2→M3→SACE on pseudo-GeoMx AOI data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--gex",
        required=True,
        metavar="PATH",
        help="Path to aoi_gex.h5ad (AOI × genes).",
    )
    parser.add_argument(
        "--protein",
        required=True,
        metavar="PATH",
        help="Path to aoi_protein.h5ad (AOI × proteins).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        metavar="DIR",
        help="Directory where outputs will be written.",
    )
    parser.add_argument(
        "--skip-m1m2",
        action="store_true",
        default=False,
        help=("Skip M1/M2 discovery and load a pre-computed profile dict from " "--profile-dict instead."),
    )
    parser.add_argument(
        "--profile-dict",
        default=None,
        metavar="PATH",
        help="Path to a JSON cell_profile_dict (required with --skip-m1m2).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Global random seed.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_pipeline(
        aoi_gex_path=args.gex,
        aoi_protein_path=args.protein,
        output_dir=args.output_dir,
        skip_m1m2=args.skip_m1m2,
        profile_dict_path=args.profile_dict,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
