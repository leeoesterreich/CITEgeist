#!/usr/bin/env python
"""GEX enrichment diagnostic: magnitude-aware profile initialization.

Compares three SACE initialization strategies using oracle proportions:
  1. pearson_r  — current marker-guided init (Pearson r with protein markers)
  2. ols_beta   — rate-space multivariate OLS regression on oracle proportions
  3. confounded — proportion-weighted mean expression (magnitude-aware but confounded)

Experiment A: Profile-level comparison (rankings, distinctness, GT cosine)
Experiment B: SACE init-sensitivity benchmark (spot_layer_r)

Usage:
    python benchmark_gex_enrichment_diagnostic.py --condition high_seg
    python benchmark_gex_enrichment_diagnostic.py --condition mixed
    python benchmark_gex_enrichment_diagnostic.py --condition high_seg --experiment A
    python benchmark_gex_enrichment_diagnostic.py --condition high_seg --experiment B
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
import scanpy as sc
from scipy import stats

# ---------------------------------------------------------------------------
# Repo path setup
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "src"))

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_SCCUBE_DIR = Path(
    str(_cg.root("DATA") / "Wu_Visium" / "Simulations" / "scCube_12k" / "replicates")
)
DEFAULT_OUTPUT_BASE = (
    REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "results_enrichment_diagnostic"
)

SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
}
SIM_TYPES = list(SIMULATION_CELL_PROFILE_DICT.keys())

RARE_TYPE_MASS_THRESHOLD = 1.0  # sum of proportions across spots

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data loading (mirrors benchmark_oracle_ceiling.py)
# ---------------------------------------------------------------------------


def load_simulation_data(sccube_dir: Path, condition: str, replicate_id: int):
    """Load simulation GEX, CITE, and GT proportions for one replicate.

    Returns:
        Tuple of (Y, antibody_data, antibody_names, gt_props, spot_names,
                  gene_names, spot_coords, gt_dir)
        where Y is (N, G) raw counts, gt_props is (N, T) oracle proportions.
    """
    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    gex_path = h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad"
    cite_path = h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad"

    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)

    Y = np.asarray(adata_gex.X, dtype=float)
    if hasattr(Y, "toarray"):
        Y = Y.toarray()

    antibody_data = np.asarray(adata_cite.X, dtype=float)
    if hasattr(antibody_data, "toarray"):
        antibody_data = antibody_data.toarray()
    antibody_names = list(adata_cite.var_names)

    # GT proportions from obs columns
    gt_props = np.zeros((adata_gex.n_obs, len(SIM_TYPES)))
    for t_idx, ct in enumerate(SIM_TYPES):
        if ct in adata_gex.obs.columns:
            gt_props[:, t_idx] = adata_gex.obs[ct].values.astype(float)
    row_sums = gt_props.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    gt_props = gt_props / row_sums

    spot_names = list(adata_gex.obs_names)
    gene_names = list(adata_gex.var_names)

    # Spatial coordinates
    if "spot_x" in adata_gex.obs.columns:
        spot_coords = np.column_stack(
            [
                adata_gex.obs["spot_x"].values.astype(float),
                adata_gex.obs["spot_y"].values.astype(float),
            ]
        )
    elif "spatial" in adata_gex.obsm:
        spot_coords = adata_gex.obsm["spatial"]
    else:
        spot_coords = np.zeros((adata_gex.n_obs, 2))

    gt_dir = sccube_dir / condition / "ST_GEX_sim" / f"sample_{replicate_id}" / "layers"

    return (Y, antibody_data, antibody_names, gt_props, spot_names, gene_names, spot_coords, gt_dir)


# ---------------------------------------------------------------------------
# Profile initialization strategies
# ---------------------------------------------------------------------------


def _postprocess_profile(H_raw: np.ndarray, gene_mask: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """Apply identical postprocessing to all init strategies.

    Matches _marker_guided_init behavior:
    1. Floor expressed genes to 0.01
    2. Keep zero-expression genes at eps
    3. Normalize rows to simplex

    Args:
        H_raw: (T, G) raw profile matrix (non-negative).
        gene_mask: (G,) bool — True for genes with nonzero total expression.
        eps: Floor for zero-expression genes.

    Returns:
        (T, G) normalized profile matrix.
    """
    H = H_raw.copy()
    H = np.maximum(H, eps)
    H[:, gene_mask] = np.maximum(H[:, gene_mask], 0.01)
    H[:, ~gene_mask] = eps
    H = H / H.sum(axis=1, keepdims=True)
    return H


def compute_pearson_r_profile(
    Y: np.ndarray,
    antibody_data: np.ndarray,
    antibody_names: list,
    type_names: list,
    cell_profile_dict: dict,
) -> np.ndarray:
    """Pearson-r init: correlate genes with mean protein signal per type.

    Mirrors _marker_guided_init from sace_gex.py with identical postprocessing.

    Returns:
        (T, G) normalized profile matrix.
    """
    from CITEgeist.model.deconvolution.detection_refinement import compute_gene_type_correlations

    H = compute_gene_type_correlations(
        Y,
        antibody_data,
        antibody_names,
        cell_profile_dict,
        type_names,
    )
    gene_mask = Y.sum(axis=0) > 0
    return _postprocess_profile(H, gene_mask)


def compute_ols_beta_profile(
    Y: np.ndarray,
    proportions: np.ndarray,
    type_names: list,
    eps: float = 1e-10,
) -> tuple:
    """Rate-space OLS regression init.

    Regresses Y/lib on proportions to get deconfounded, magnitude-aware
    type profiles. Uses lstsq (not closed-form inverse) for numerical safety
    with compositional predictors.

    Args:
        Y: (N, G) raw spot counts.
        proportions: (N, T) oracle proportions.
        type_names: Ordered type names.
        eps: Numerical floor.

    Returns:
        Tuple of (mu_global, diagnostics) where:
        - mu_global: (T, G) normalized profile matrix
        - diagnostics: dict with condition_number, neg_frac_per_type, per_type_mass
    """
    N, G = Y.shape
    T = len(type_names)

    lib = Y.sum(axis=1, keepdims=True)  # (N, 1)
    lib = np.where(lib > 0, lib, 1.0)
    Y_rate = Y / lib  # (N, G)

    beta, residuals, rank, sv = np.linalg.lstsq(proportions, Y_rate, rcond=None)
    cond = float(sv[0] / sv[-1]) if len(sv) > 0 and sv[-1] > 0 else float("inf")

    neg_frac = (beta < 0).mean(axis=1).tolist()  # per-type
    per_type_mass = proportions.sum(axis=0).tolist()

    # Clip negatives
    beta = np.maximum(beta, eps)

    gene_mask = Y.sum(axis=0) > 0
    mu_global = _postprocess_profile(beta, gene_mask, eps=eps)

    diagnostics = {
        "condition_number": cond,
        "neg_frac_per_type": dict(zip(type_names, neg_frac)),
        "per_type_mass": dict(zip(type_names, per_type_mass)),
        "lib_size_cv": float(np.std(lib.ravel()) / np.mean(lib.ravel())),
        "rank": int(rank) if rank is not None else T,
    }
    return mu_global, diagnostics


def compute_confounded_profile(
    Y: np.ndarray,
    proportions: np.ndarray,
    type_names: list,
    eps: float = 1e-10,
) -> np.ndarray:
    """Confounded (proportion-weighted mean) init in rate space.

    mu[t, g] = sum_s(prop[s,t] * Y[s,g] / lib[s]) / sum_s(prop[s,t])

    This is magnitude-aware but confounded when proportions are correlated.

    Args:
        Y: (N, G) raw spot counts.
        proportions: (N, T) oracle proportions.
        type_names: Ordered type names.
        eps: Numerical floor.

    Returns:
        (T, G) normalized profile matrix.
    """
    T = len(type_names)
    G = Y.shape[1]

    lib = Y.sum(axis=1, keepdims=True)
    lib = np.where(lib > 0, lib, 1.0)
    Y_rate = Y / lib

    H = np.zeros((T, G))
    for t in range(T):
        weights = proportions[:, t]
        total_w = weights.sum()
        if total_w > 0:
            H[t] = (weights[:, None] * Y_rate).sum(axis=0) / total_w
        else:
            H[t] = 1.0 / G

    gene_mask = Y.sum(axis=0) > 0
    return _postprocess_profile(H, gene_mask, eps=eps)


# ---------------------------------------------------------------------------
# Experiment A: Profile-level diagnostic
# ---------------------------------------------------------------------------


def _gini_coefficient(x: np.ndarray) -> float:
    """Gini coefficient of a 1D array."""
    x = np.sort(np.abs(x))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return float((2 * np.sum(index * x) - (n + 1) * np.sum(x)) / (n * np.sum(x)))


def _cosine_sim(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine similarity between two vectors."""
    dot = np.dot(a, b)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-15 or nb < 1e-15:
        return 0.0
    return float(dot / (na * nb))


def compare_profiles(
    profiles: dict,
    type_names: list,
    gt_profile: np.ndarray = None,
    type_masses: np.ndarray = None,
    k: int = 50,
) -> dict:
    """Compare multiple profile matrices on ranking and distinctness metrics.

    Args:
        profiles: Dict[init_name -> (T, G) profile matrix].
        type_names: Ordered type names.
        gt_profile: Optional (T, G) ground truth profile (simulation only).
        type_masses: Optional (T,) total proportion mass per type.
        k: Number of top genes to compare.

    Returns:
        Dict with per-type and aggregate metrics.
    """
    init_names = list(profiles.keys())
    T = len(type_names)

    results = {"per_type": {}, "aggregate": {}}

    for t_idx, ct in enumerate(type_names):
        is_rare = type_masses is not None and type_masses[t_idx] < RARE_TYPE_MASS_THRESHOLD
        type_results = {"is_rare": is_rare, "mass": float(type_masses[t_idx]) if type_masses is not None else None}

        # Per-init metrics
        for init_name in init_names:
            H = profiles[init_name]
            row = H[t_idx]
            top_k_idx = np.argsort(-row)[:k]
            type_results[f"{init_name}_gini_top{k}"] = _gini_coefficient(row[top_k_idx])

            if gt_profile is not None:
                type_results[f"{init_name}_gt_cosine"] = _cosine_sim(row, gt_profile[t_idx])

        # Pairwise comparisons between inits
        for i, name_a in enumerate(init_names):
            for name_b in init_names[i + 1 :]:
                Ha = profiles[name_a][t_idx]
                Hb = profiles[name_b][t_idx]

                # Top-k overlap
                top_a = set(np.argsort(-Ha)[:k])
                top_b = set(np.argsort(-Hb)[:k])
                overlap = len(top_a & top_b) / k

                # Spearman rank correlation
                rho, _ = stats.spearmanr(Ha, Hb)

                pair_key = f"{name_a}_vs_{name_b}"
                type_results[f"{pair_key}_top{k}_overlap"] = overlap
                type_results[f"{pair_key}_spearman_rho"] = float(rho)

        results["per_type"][ct] = type_results

    # Profile distinctness: mean pairwise cosine distance of rows
    for init_name in init_names:
        H = profiles[init_name]
        H_norm = H / (np.linalg.norm(H, axis=1, keepdims=True) + 1e-15)
        cos_matrix = H_norm @ H_norm.T
        mask = ~np.eye(T, dtype=bool)
        mean_cos_dist = float(1.0 - cos_matrix[mask].mean())
        results["aggregate"][f"{init_name}_mean_cosine_distance"] = mean_cos_dist

    # Aggregate averages (unweighted and mass-weighted)
    metric_keys = []
    for key in list(results["per_type"].values())[0]:
        if key not in ("is_rare", "mass") and isinstance(list(results["per_type"].values())[0][key], float):
            metric_keys.append(key)

    for mk in metric_keys:
        vals = [results["per_type"][ct][mk] for ct in type_names if results["per_type"][ct][mk] is not None]
        if vals:
            results["aggregate"][f"mean_{mk}"] = float(np.mean(vals))
            if type_masses is not None:
                masses = np.array(
                    [type_masses[i] for i, ct in enumerate(type_names) if results["per_type"][ct][mk] is not None]
                )
                if masses.sum() > 0:
                    results["aggregate"][f"wmean_{mk}"] = float(np.average(vals, weights=masses / masses.sum()))

    return results


def run_experiment_a(
    Y: np.ndarray,
    antibody_data: np.ndarray,
    antibody_names: list,
    gt_props: np.ndarray,
    gene_names: list,
    gt_dir: Path,
    type_names: list,
    cell_profile_dict: dict,
) -> dict:
    """Run Experiment A: profile-level diagnostic.

    Args:
        Y: (N, G) raw spot counts.
        antibody_data: (N, M) antibody matrix.
        antibody_names: Marker names.
        gt_props: (N, T) oracle proportions.
        gene_names: Gene names.
        gt_dir: Path to GT layer CSVs for GT profile extraction.
        type_names: Ordered type names.
        cell_profile_dict: Marker-type mapping.

    Returns:
        Dict with profile comparison results and diagnostics.
    """
    logger.info("=== Experiment A: Profile-Level Diagnostic ===")

    # Compute three profiles
    H_pearson = compute_pearson_r_profile(
        Y,
        antibody_data,
        antibody_names,
        type_names,
        cell_profile_dict,
    )
    H_beta, ols_diag = compute_ols_beta_profile(Y, gt_props, type_names)
    H_conf = compute_confounded_profile(Y, gt_props, type_names)

    profiles = {"pearson_r": H_pearson, "ols_beta": H_beta, "confounded": H_conf}

    # GT profiles from simulation layers (if available)
    gt_profile = None
    if gt_dir.exists():
        from benchmark_oracle_ceiling import load_gt_layers, extract_global_profiles

        try:
            gt_layers = load_gt_layers(gt_dir, type_names)
            gt_profile = extract_global_profiles(gt_layers, type_names, gene_names)
            logger.info("Loaded GT profiles for cosine comparison")
        except Exception as e:
            logger.warning("Could not load GT profiles: %s", e)

    type_masses = gt_props.sum(axis=0)

    comparison = compare_profiles(
        profiles,
        type_names,
        gt_profile=gt_profile,
        type_masses=type_masses,
    )

    return {
        "comparison": comparison,
        "ols_diagnostics": ols_diag,
    }


# ---------------------------------------------------------------------------
# Experiment B: SACE init-sensitivity benchmark
# ---------------------------------------------------------------------------


def _mu_global_to_spotwise(mu_global: np.ndarray, n_spots: int) -> dict:
    """Convert a global (T, G) profile to spotwise_profiles_init format.

    run_sace uses spotwise_profiles_init as first-priority init. By setting
    every spot to the same mu_global, we override marker-guided and confounded
    init paths.

    Args:
        mu_global: (T, G) normalized type profiles.
        n_spots: Number of spots.

    Returns:
        Dict[spot_idx -> (T, G) array].
    """
    return {s: mu_global for s in range(n_spots)}


def _build_cell_assignments_from_props(
    gt_props: np.ndarray,
    spot_names: list,
    type_names: list,
    n_cells_per_spot: int = 5,
) -> tuple:
    """Build synthetic cell assignments from oracle proportions.

    Creates n_cells_per_spot cells per spot, assigned by proportion-weighted
    sampling (deterministic: round to nearest integer counts).

    Args:
        gt_props: (N, T) oracle proportions.
        spot_names: Spot names.
        type_names: Type names.
        n_cells_per_spot: Cells per spot.

    Returns:
        Tuple of (cell_assignments, cell_spot_map) where:
        - cell_assignments: Dict[cell_id -> type_name]
        - cell_spot_map: DataFrame with cell_id, spot_barcode, spot_idx, x, y
    """
    N = len(spot_names)
    T = len(type_names)
    cell_assignments = {}
    rows = []

    for s in range(N):
        # Deterministic rounding of proportions to cell counts
        raw_counts = gt_props[s] * n_cells_per_spot
        counts = np.round(raw_counts).astype(int)
        # Adjust to hit exact total
        diff = n_cells_per_spot - counts.sum()
        if diff > 0:
            # Add to largest fractional remainders
            remainders = raw_counts - np.floor(raw_counts)
            for _ in range(diff):
                idx = np.argmax(remainders)
                counts[idx] += 1
                remainders[idx] = -1
        elif diff < 0:
            # Remove from smallest counts > 0
            for _ in range(-diff):
                nonzero = np.where(counts > 0)[0]
                idx = nonzero[np.argmin(counts[nonzero])]
                counts[idx] -= 1

        cell_idx = 0
        for t in range(T):
            for _ in range(counts[t]):
                cid = f"cell_{s}_{cell_idx}"
                cell_assignments[cid] = type_names[t]
                rows.append(
                    {
                        "cell_id": cid,
                        "spot_barcode": spot_names[s],
                        "spot_idx": s,
                        "x": 0.0,  # Placeholder — SACE only uses spot_idx
                        "y": 0.0,
                    }
                )
                cell_idx += 1

    cell_spot_map = pd.DataFrame(rows)
    return cell_assignments, cell_spot_map


def run_experiment_b(
    Y: np.ndarray,
    antibody_data: np.ndarray,
    antibody_names: list,
    gt_props: np.ndarray,
    spot_names: list,
    gene_names: list,
    spot_coords: np.ndarray,
    gt_dir: Path,
    type_names: list,
    cell_profile_dict: dict,
) -> dict:
    """Run Experiment B: SACE init-sensitivity benchmark.

    Runs SACE with max_iter=1 using three different mu_global inits,
    then evaluates spot_layer_r against ground truth layers.

    Args:
        Y: (N, G) raw spot counts.
        antibody_data: (N, M) antibody matrix.
        antibody_names: Marker names.
        gt_props: (N, T) oracle proportions.
        spot_names: Spot names.
        gene_names: Gene names.
        spot_coords: (N, 2) spatial coordinates.
        gt_dir: Path to GT layer CSVs.
        type_names: Ordered type names.
        cell_profile_dict: Marker-type mapping.

    Returns:
        Dict with per-init metrics and diagnostics.
    """
    from CITEgeist.model.gex.sace_gex import run_sace
    from benchmarking_gex import calculate_rmse

    logger.info("=== Experiment B: SACE Init-Sensitivity Benchmark ===")

    N = Y.shape[0]
    props_df = pd.DataFrame(gt_props, index=spot_names, columns=type_names)

    # Build cell assignments from oracle proportions
    cell_assignments, cell_spot_map = _build_cell_assignments_from_props(
        gt_props,
        spot_names,
        type_names,
    )
    logger.info("Built %d cell assignments across %d spots", len(cell_assignments), N)

    # Compute the three mu_global profiles
    H_pearson = compute_pearson_r_profile(
        Y,
        antibody_data,
        antibody_names,
        type_names,
        cell_profile_dict,
    )
    H_beta, ols_diag = compute_ols_beta_profile(Y, gt_props, type_names)
    # H_confounded is triggered by passing antibody_data=None

    init_configs = {
        "pearson_r": {
            "spotwise_profiles_init": _mu_global_to_spotwise(H_pearson, N),
            "antibody_data": None,
            "antibody_names": None,
            "cell_profile_dict": None,
        },
        "ols_beta": {
            "spotwise_profiles_init": _mu_global_to_spotwise(H_beta, N),
            "antibody_data": None,
            "antibody_names": None,
            "cell_profile_dict": None,
        },
        "confounded": {
            "spotwise_profiles_init": None,
            "antibody_data": None,
            "antibody_names": None,
            "cell_profile_dict": None,
        },
    }

    results = {}
    for init_name, config in init_configs.items():
        logger.info("--- Running SACE with init: %s ---", init_name)
        t0 = time.time()

        spot_type_gex, cell_adata, diagnostics = run_sace(
            spot_counts=Y,
            proportions=props_df,
            cell_assignments=cell_assignments,
            cell_spot_map=cell_spot_map,
            spot_coords=spot_coords,
            gene_names=gene_names,
            max_iter=1,
            **config,
        )

        elapsed = time.time() - t0

        # Convert spot_type_gex Dict[spot_idx -> (T, G)] to per-type layer DataFrames
        # This is much faster than iterating over per-cell AnnData.
        layers = {}
        type_to_idx = {ct: i for i, ct in enumerate(type_names)}
        for ct in type_names:
            t_idx = type_to_idx[ct]
            layer_arr = np.array([spot_type_gex[s][t_idx] for s in range(len(spot_names))])
            layers[ct] = pd.DataFrame(layer_arr, index=spot_names, columns=gene_names)

        # Export and evaluate
        layers_dir = Path(f"/tmp/enrichment_diag_{init_name}")
        layers_dir.mkdir(parents=True, exist_ok=True)
        for ct, df in layers.items():
            df.to_csv(layers_dir / f"{ct}_layer.csv")

        metrics = calculate_rmse(
            ground_truth_dir=str(gt_dir),
            predictions_dir=str(layers_dir),
            normalize="range",
        )

        results[init_name] = {
            "spot_layer_r": metrics.get("mean_spot_layer_r", float("nan")),
            "pearson_r": metrics.get("mean_pearson_r", float("nan")),
            "avg_nrmse": metrics.get("average_nrmse", float("nan")),
            "elapsed_s": elapsed,
            "per_type": {
                ct: {k: float(v) for k, v in cm.items()} for ct, cm in metrics.get("metrics_per_cell_type", {}).items()
            },
        }

        logger.info(
            "%s: spot_layer_r=%.4f  (%.1fs)",
            init_name,
            results[init_name]["spot_layer_r"],
            elapsed,
        )

    results["ols_diagnostics"] = ols_diag
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def print_summary_table(all_results: dict) -> None:
    """Print formatted summary tables for Experiments A and B."""
    # --- Experiment A summary ---
    if "experiment_a" in all_results:
        print("\n" + "=" * 70)
        print("EXPERIMENT A: Profile-Level Diagnostic")
        print("=" * 70)

        for cond, cond_results in all_results["experiment_a"].items():
            print(f"\n--- Condition: {cond} ---")
            for rep_key, rep_result in cond_results.items():
                agg = rep_result["comparison"]["aggregate"]
                print(f"\n  {rep_key}:")

                # GT cosine (most important)
                for init in ["pearson_r", "ols_beta", "confounded"]:
                    gt_key = f"mean_{init}_gt_cosine"
                    if gt_key in agg:
                        print(f"    {init:12s} GT cosine: {agg[gt_key]:.4f}")

                # Profile distinctness
                for init in ["pearson_r", "ols_beta", "confounded"]:
                    dist_key = f"{init}_mean_cosine_distance"
                    if dist_key in agg:
                        print(f"    {init:12s} distinctness: {agg[dist_key]:.4f}")

                # Top-50 overlaps
                for pair in ["pearson_r_vs_ols_beta", "pearson_r_vs_confounded", "ols_beta_vs_confounded"]:
                    key = f"mean_{pair}_top50_overlap"
                    if key in agg:
                        print(f"    {pair} top50 overlap: {agg[key]:.3f}")

    # --- Experiment B summary ---
    if "experiment_b" in all_results:
        print("\n" + "=" * 70)
        print("EXPERIMENT B: SACE Init-Sensitivity Benchmark")
        print("=" * 70)
        print(f"\n{'Condition':<12} {'Rep':<8} {'pearson_r':>10} {'ols_beta':>10} {'confounded':>10}")
        print("-" * 55)

        for cond, cond_results in all_results["experiment_b"].items():
            for rep_key, rep_result in cond_results.items():
                vals = []
                for init in ["pearson_r", "ols_beta", "confounded"]:
                    vals.append(rep_result.get(init, {}).get("spot_layer_r", float("nan")))
                print(f"{cond:<12} {rep_key:<8} {vals[0]:>10.4f} {vals[1]:>10.4f} {vals[2]:>10.4f}")

        # Aggregate means
        print("-" * 55)
        for cond in all_results["experiment_b"]:
            if not all_results["experiment_b"][cond]:
                continue
            for init in ["pearson_r", "ols_beta", "confounded"]:
                vals = [
                    r.get(init, {}).get("spot_layer_r", float("nan"))
                    for r in all_results["experiment_b"][cond].values()
                ]
                finite_vals = [v for v in vals if not np.isnan(v)]
                mean_val = float(np.mean(finite_vals)) if finite_vals else float("nan")
                print(f"  {cond} mean {init}: {mean_val:.4f}")


def main():
    parser = argparse.ArgumentParser(description="GEX enrichment diagnostic: magnitude-aware init comparison")
    parser.add_argument(
        "--condition", type=str, default=None, help="Condition to run (high_seg, mixed). Default: both."
    )
    parser.add_argument("--replicate", type=int, default=None, help="Single replicate to run (0-4). Default: all 5.")
    parser.add_argument(
        "--experiment", type=str, default=None, choices=["A", "B"], help="Run only one experiment. Default: both."
    )
    parser.add_argument("--sccube_dir", type=str, default=str(DEFAULT_SCCUBE_DIR))
    parser.add_argument("--output_dir", type=str, default=str(DEFAULT_OUTPUT_BASE))
    args = parser.parse_args()

    sccube_dir = Path(args.sccube_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    conditions = [args.condition] if args.condition else ["high_seg", "mixed"]
    replicates = [args.replicate] if args.replicate is not None else list(range(5))

    all_results = {"experiment_a": {}, "experiment_b": {}}

    for condition in conditions:
        all_results["experiment_a"][condition] = {}
        all_results["experiment_b"][condition] = {}

        for rep_id in replicates:
            rep_key = f"rep_{rep_id}"
            logger.info("\n>>> Condition=%s, Replicate=%d <<<", condition, rep_id)

            Y, ab, ab_names, props, spots, genes, coords, gt_dir = load_simulation_data(
                sccube_dir,
                condition,
                rep_id,
            )

            if args.experiment is None or args.experiment == "A":
                result_a = run_experiment_a(
                    Y,
                    ab,
                    ab_names,
                    props,
                    genes,
                    gt_dir,
                    SIM_TYPES,
                    SIMULATION_CELL_PROFILE_DICT,
                )
                all_results["experiment_a"][condition][rep_key] = result_a

            if args.experiment is None or args.experiment == "B":
                result_b = run_experiment_b(
                    Y,
                    ab,
                    ab_names,
                    props,
                    spots,
                    genes,
                    coords,
                    gt_dir,
                    SIM_TYPES,
                    SIMULATION_CELL_PROFILE_DICT,
                )
                all_results["experiment_b"][condition][rep_key] = result_b

    # Print summary
    print_summary_table(all_results)

    # Save results JSON (condition-specific to avoid overwrites)
    cond_tag = conditions[0] if len(conditions) == 1 else "all"
    results_path = output_dir / f"enrichment_diagnostic_{cond_tag}.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    logger.info("Results saved to %s", results_path)


if __name__ == "__main__":
    main()
