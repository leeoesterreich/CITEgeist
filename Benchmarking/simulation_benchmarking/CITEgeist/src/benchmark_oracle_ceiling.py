#!/usr/bin/env python
"""Oracle ceiling experiments for SACE GEX gap decomposition.

Runs 4 experiments per replicate to decompose the SACE↔C2L gap into
profile quality vs. proportion quality vs. allocation mechanism.

Experiments:
  1. oracle_W    — GT type profiles + QP proportions
  2. oracle_H    — marker-guided W + GT proportions
  3. oracle_both — GT type profiles + GT proportions
  4. cross_method — C2L-derived profiles + QP proportions

Each experiment: construct W, H → single E-step → export layers → evaluate.

Usage:
    python benchmark_oracle_ceiling.py --replicate 0 --condition high_seg
    python benchmark_oracle_ceiling.py --replicate 0 --condition high_seg --experiments oracle_W oracle_H
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

# ---------------------------------------------------------------------------
# Repo path setup
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "src"))

# ---------------------------------------------------------------------------
# Default paths
# ---------------------------------------------------------------------------
DEFAULT_SCCUBE_DIR = Path(
    str(_cg.root("DATA") / "Wu_Visium" / "Simulations" / "scCube_12k" / "replicates")
)
DEFAULT_OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "results_oracle"
DEFAULT_C2L_BASE = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "Cell2Location"
DEFAULT_SACE_BASE = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "results_sace"

# Cell type profile dict — 9 types, 2 Major markers each (simulation)
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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------


def load_gt_layers(gt_dir: Path, type_names: list) -> dict:
    """Load ground truth per-type GEX layers as DataFrames.

    Args:
        gt_dir: Directory containing {CellType}_GT.csv files.
        type_names: Ordered list of cell type names.

    Returns:
        Dict[type_name -> DataFrame (spots x genes)].
    """
    layers = {}
    for ct in type_names:
        gt_path = gt_dir / f"{ct}_GT.csv"
        if not gt_path.exists():
            raise FileNotFoundError(f"GT file not found: {gt_path}")
        layers[ct] = pd.read_csv(gt_path, index_col=0)
    return layers


def extract_global_profiles(layers: dict, type_names: list, gene_names: list, eps: float = 1e-10) -> np.ndarray:
    """Extract global type profiles W from per-type GEX layers.

    W[t, g] = sum_s(layer[t][s, g]) / sum_s(sum_g(layer[t][s, g]))
    This is the MLE type profile given per-spot per-type counts.

    Args:
        layers: Dict[type_name -> DataFrame (spots x genes)].
        type_names: Ordered type names.
        gene_names: Ordered gene names.
        eps: Floor to avoid division by zero.

    Returns:
        (T, G) normalized type profile matrix.
    """
    T = len(type_names)
    G = len(gene_names)
    W = np.zeros((T, G))
    for t_idx, ct in enumerate(type_names):
        df = layers[ct].reindex(columns=gene_names, fill_value=0.0)
        row_sums = df.values.sum()
        if row_sums > 0:
            W[t_idx] = df.values.sum(axis=0) / row_sums
        else:
            W[t_idx] = 1.0 / G
    W = np.maximum(W, eps)
    W = W / W.sum(axis=1, keepdims=True)
    return W


def extract_gt_proportions(adata_gex: sc.AnnData, type_names: list) -> np.ndarray:
    """Extract ground truth proportions from GEX h5ad obs columns.

    Args:
        adata_gex: GEX AnnData with type columns in .obs.
        type_names: Ordered type names matching obs columns.

    Returns:
        (N, T) proportions matrix.
    """
    props = np.zeros((adata_gex.n_obs, len(type_names)))
    for t_idx, ct in enumerate(type_names):
        if ct in adata_gex.obs.columns:
            props[:, t_idx] = adata_gex.obs[ct].values.astype(float)
        else:
            logger.warning("Type '%s' not in obs columns", ct)
    # Normalize rows to sum to 1
    row_sums = props.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    props = props / row_sums
    return props


def load_qp_proportions(
    sace_dir: Path,
    condition: str,
    replicate_id: int,
    type_names: list,
    adata_gex: sc.AnnData,
) -> np.ndarray:
    """Load QP proportions from existing SACE benchmark results.

    Falls back to running the CitegeistModel pipeline if no cached results.

    Args:
        sace_dir: Base SACE results directory.
        condition: "high_seg" or "mixed".
        replicate_id: 0-4.
        type_names: Ordered type names.
        adata_gex: GEX AnnData for spot ordering.

    Returns:
        (N, T) proportions matrix aligned to adata_gex.obs_names.
    """
    # Try loading from SACE results (the model saves proportions during benchmark)
    model_dir = sace_dir / condition / f"Wu_rep_{replicate_id}" / "citegeist"
    prop_path = model_dir / "cell_prop_global_results.csv"
    if prop_path.exists():
        props_df = pd.read_csv(prop_path, index_col=0)
        # Keep only type columns, reindex to adata spot order
        type_cols = [c for c in type_names if c in props_df.columns]
        props_df = props_df[type_cols].reindex(adata_gex.obs_names)
        n_types = len(type_names)
        props_df = props_df.fillna(1.0 / n_types)
        return props_df.values

    # Fallback: run CitegeistModel to get proportions
    logger.info("No cached proportions found, running CitegeistModel...")
    return _run_proportion_pipeline(adata_gex, condition, replicate_id, type_names)


def _run_proportion_pipeline(
    adata_gex: sc.AnnData,
    condition: str,
    replicate_id: int,
    type_names: list,
) -> np.ndarray:
    """Run full CitegeistModel proportion pipeline as fallback.

    Returns:
        (N, T) proportions aligned to adata_gex.obs_names.
    """
    from CITEgeist.model.citegeist_model import CitegeistModel

    sccube_dir = DEFAULT_SCCUBE_DIR
    cite_path = sccube_dir / condition / "h5ad_objects" / f"Wu_rep_{replicate_id}_CITE.h5ad"
    adata_cite = sc.read_h5ad(cite_path)

    model = CitegeistModel(
        sample_name=f"Wu_rep_{replicate_id}",
        output_folder=f"/tmp/oracle_props_{condition}_{replicate_id}",
        simulation=True,
        gene_expression_adata=adata_gex.copy(),
        antibody_capture_adata=adata_cite,
    )
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()
    model.cell_profile_dict = SIMULATION_CELL_PROFILE_DICT

    global_props_df, _ = model.run_cell_proportion_model(method="qp", max_workers=8)
    props_df = global_props_df.drop(columns=["recon_error"], errors="ignore")
    props_df = props_df.reindex(adata_gex.obs_names)
    props_df = props_df.fillna(1.0 / len(type_names))
    return props_df[type_names].values


def load_c2l_layers(c2l_dir: Path, condition: str, replicate_id: int, type_names: list) -> dict:
    """Load Cell2Location per-type GEX layers.

    Args:
        c2l_dir: Base Cell2Location results directory.
        condition: "high_seg" or "mixed".
        replicate_id: 0-4.
        type_names: Ordered type names.

    Returns:
        Dict[type_name -> DataFrame (spots x genes)].
    """
    layers_dir = c2l_dir / condition / f"cell2location_map_{replicate_id}" / "layers"
    layers = {}
    for ct in type_names:
        path = layers_dir / f"{ct}_layer.csv"
        if not path.exists():
            raise FileNotFoundError(f"C2L layer not found: {path}")
        layers[ct] = pd.read_csv(path, index_col=0)
    return layers


# ---------------------------------------------------------------------------
# Core E-step allocation
# ---------------------------------------------------------------------------


def e_step_allocate(Y: np.ndarray, W: np.ndarray, H: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """Single E-step: allocate spot counts to types.

    E_x[s, t, g] = Y[s, g] * (B[s, t] * W[t, g]) / sum_t'(B[s, t'] * W[t', g])

    where B[s, t] = H[s, t] * lib_size[s] (library-size-scaled proportions).

    Args:
        Y: (N, G) raw spot counts.
        W: (T, G) normalized type profiles (rows sum to 1).
        H: (N, T) proportions (rows sum to 1).
        eps: Numerical floor.

    Returns:
        (N, T, G) allocated counts. sum over t equals Y for each (s, g).
    """
    N, G = Y.shape
    T = W.shape[0]
    lib_sizes = Y.sum(axis=1)  # (N,)
    B = H * lib_sizes[:, None]  # (N, T)

    # weights[s, t, g] = B[s, t] * W[t, g]
    weights = B[:, :, None] * W[None, :, :]  # (N, T, G)
    lam = weights.sum(axis=1)  # (N, G)
    lam_safe = np.where(lam > 0, lam, 1.0)

    E_x = Y[:, None, :] * weights / lam_safe[:, None, :]
    E_x = np.where(lam[:, None, :] > 0, E_x, 0.0)
    return E_x


def e_step_to_layers(E_x: np.ndarray, spot_names: list, gene_names: list, type_names: list) -> dict:
    """Convert E_x (N, T, G) to per-type DataFrames for evaluation.

    Args:
        E_x: (N, T, G) allocated counts.
        spot_names: Ordered spot names.
        gene_names: Ordered gene names.
        type_names: Ordered type names.

    Returns:
        Dict[type_name -> DataFrame (spots x genes)].
    """
    layers = {}
    for t_idx, ct in enumerate(type_names):
        layers[ct] = pd.DataFrame(
            E_x[:, t_idx, :],
            index=spot_names,
            columns=gene_names,
        )
    return layers


def export_layers(layers: dict, output_dir: Path) -> None:
    """Write per-cell-type layer CSVs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for ct, df in layers.items():
        df.to_csv(output_dir / f"{ct}_layer.csv")


# ---------------------------------------------------------------------------
# Experiment runners
# ---------------------------------------------------------------------------


def run_oracle_experiment(
    experiment_name: str,
    Y: np.ndarray,
    W: np.ndarray,
    H: np.ndarray,
    spot_names: list,
    gene_names: list,
    type_names: list,
    gt_dir: Path,
    output_dir: Path,
) -> dict:
    """Run one oracle experiment: E-step → export → evaluate.

    Args:
        experiment_name: Name for logging and output subdirectory.
        Y: (N, G) raw spot counts.
        W: (T, G) type profiles.
        H: (N, T) proportions.
        spot_names: Spot names for output DataFrames.
        gene_names: Gene names for output DataFrames.
        type_names: Cell type names.
        gt_dir: Ground truth directory for evaluation.
        output_dir: Output directory for this experiment's layers.

    Returns:
        Dict with experiment name, metrics, and timing.
    """
    from benchmarking_gex import calculate_rmse

    logger.info("=== Experiment: %s ===", experiment_name)
    t0 = time.time()

    # E-step allocation
    E_x = e_step_allocate(Y, W, H)

    # Convert to layers and export
    layers = e_step_to_layers(E_x, spot_names, gene_names, type_names)
    layers_dir = output_dir / experiment_name / "layers"
    export_layers(layers, layers_dir)

    # Evaluate
    metrics = calculate_rmse(
        ground_truth_dir=str(gt_dir),
        predictions_dir=str(layers_dir),
        normalize="range",
    )

    elapsed = time.time() - t0
    logger.info(
        "%s | spot_layer_r=%.4f  pearson_r=%.4f  nrmse=%.4f  (%.1fs)",
        experiment_name,
        metrics.get("mean_spot_layer_r", float("nan")),
        metrics.get("mean_pearson_r", float("nan")),
        metrics.get("average_nrmse", float("nan")),
        elapsed,
    )

    return {
        "experiment": experiment_name,
        "elapsed_s": elapsed,
        **{k: v for k, v in metrics.items() if k != "metrics_per_cell_type"},
        "metrics_per_cell_type": {
            ct: {mk: float(mv) for mk, mv in cm.items()} for ct, cm in metrics.get("metrics_per_cell_type", {}).items()
        },
    }


def run_all_experiments(
    replicate_id: int,
    condition: str,
    experiments: list,
    sccube_dir: Path,
    output_base: Path,
    c2l_base: Path,
    sace_base: Path,
) -> list:
    """Run all requested oracle experiments for one replicate.

    Args:
        replicate_id: 0-4.
        condition: "high_seg" or "mixed".
        experiments: List of experiment names to run.
        sccube_dir: Root scCube replicates directory.
        output_base: Root output directory.
        c2l_base: Cell2Location results base directory.
        sace_base: SACE results base directory.

    Returns:
        List of result dicts, one per experiment.
    """
    # Load data
    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    gex_path = h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad"
    cite_path = h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad"
    logger.info("Loading GEX: %s", gex_path)
    adata_gex = sc.read_h5ad(gex_path)
    logger.info("Loading CITE: %s", cite_path)
    adata_cite = sc.read_h5ad(cite_path)

    Y = np.asarray(adata_gex.X, dtype=float)  # (N, G) raw counts
    if hasattr(Y, "toarray"):
        Y = Y.toarray()
    N, G = Y.shape
    spot_names = list(adata_gex.obs_names)
    gene_names = list(adata_gex.var_names)
    type_names = SIM_TYPES
    T = len(type_names)

    gt_dir = sccube_dir / condition / "ST_GEX_sim" / f"sample_{replicate_id}" / "layers"

    output_dir = output_base / condition / f"Wu_rep_{replicate_id}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Pre-compute shared components ---

    # GT layers (needed for oracle_W, oracle_both, and as evaluation GT)
    gt_layers = load_gt_layers(gt_dir, type_names)
    W_oracle = extract_global_profiles(gt_layers, type_names, gene_names)
    logger.info("Oracle W: shape=%s, row sums=%s", W_oracle.shape, W_oracle.sum(axis=1)[:3].round(4))

    # GT proportions (needed for oracle_H, oracle_both)
    H_oracle = extract_gt_proportions(adata_gex, type_names)
    logger.info(
        "Oracle H: shape=%s, row sums range=[%.3f, %.3f]",
        H_oracle.shape,
        H_oracle.sum(axis=1).min(),
        H_oracle.sum(axis=1).max(),
    )

    # QP proportions (needed for oracle_W, cross_method)
    H_qp = load_qp_proportions(sace_base, condition, replicate_id, type_names, adata_gex)
    logger.info("QP H: shape=%s", H_qp.shape)

    # Marker-guided W (needed for oracle_H)
    from CITEgeist.model.gex.sace_gex import _marker_guided_init

    antibody_data = np.asarray(adata_cite.X, dtype=float)
    if hasattr(antibody_data, "toarray"):
        antibody_data = antibody_data.toarray()
    antibody_names = list(adata_cite.var_names)
    W_marker, mg_diag = _marker_guided_init(
        Y,
        antibody_data,
        antibody_names,
        SIMULATION_CELL_PROFILE_DICT,
        type_names,
    )
    if W_marker is None:
        logger.warning("Marker-guided init failed, using confounded init for oracle_H")
        # Confounded fallback
        W_marker = np.zeros((T, G))
        for t in range(T):
            weighted = H_qp[:, t : t + 1] * Y
            total = weighted.sum()
            if total > 0:
                W_marker[t] = weighted.sum(axis=0) / total
            else:
                W_marker[t] = 1.0 / G
        W_marker = W_marker / W_marker.sum(axis=1, keepdims=True)
    logger.info("Marker-guided W: mean_row_corr=%.4f", mg_diag.get("mean_row_corr", -1))

    # --- Run experiments ---
    results = []

    experiment_configs = {
        "oracle_W": (W_oracle, H_qp),
        "oracle_H": (W_marker, H_oracle),
        "oracle_both": (W_oracle, H_oracle),
    }

    # cross_method needs C2L layers
    if "cross_method" in experiments:
        try:
            c2l_layers = load_c2l_layers(c2l_base, condition, replicate_id, type_names)
            W_c2l = extract_global_profiles(c2l_layers, type_names, gene_names)
            experiment_configs["cross_method"] = (W_c2l, H_qp)
            logger.info("C2L profiles loaded for cross_method experiment")
        except FileNotFoundError as e:
            logger.warning("Skipping cross_method: %s", e)

    for exp_name in experiments:
        if exp_name not in experiment_configs:
            logger.warning("Skipping unknown experiment: %s", exp_name)
            continue
        W_exp, H_exp = experiment_configs[exp_name]
        result = run_oracle_experiment(
            experiment_name=exp_name,
            Y=Y,
            W=W_exp,
            H=H_exp,
            spot_names=spot_names,
            gene_names=gene_names,
            type_names=type_names,
            gt_dir=gt_dir,
            output_dir=output_dir,
        )
        results.append(result)

    # Save combined results
    summary_path = output_dir / "oracle_results.json"
    with open(summary_path, "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    logger.info("Results saved: %s", summary_path)

    # Print comparison table
    logger.info("\n=== ORACLE CEILING RESULTS (rep=%d, cond=%s) ===", replicate_id, condition)
    logger.info("%-15s  %12s  %12s  %12s", "Experiment", "spot_layer_r", "pearson_r", "avg_nrmse")
    for r in results:
        logger.info(
            "%-15s  %12.4f  %12.4f  %12.4f",
            r["experiment"],
            r.get("mean_spot_layer_r", float("nan")),
            r.get("mean_pearson_r", float("nan")),
            r.get("average_nrmse", float("nan")),
        )

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Oracle ceiling experiments for SACE GEX gap decomposition.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--replicate", type=int, default=0, help="Replicate ID (0-4)")
    parser.add_argument(
        "--condition",
        type=str,
        default="high_seg",
        choices=["high_seg", "mixed"],
        help="Simulation condition",
    )
    parser.add_argument(
        "--experiments",
        nargs="+",
        default=["oracle_W", "oracle_H", "oracle_both", "cross_method"],
        help="Experiments to run",
    )
    parser.add_argument("--sccube_dir", type=Path, default=DEFAULT_SCCUBE_DIR)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument("--c2l_dir", type=Path, default=DEFAULT_C2L_BASE)
    parser.add_argument("--sace_dir", type=Path, default=DEFAULT_SACE_BASE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger.info(
        "Oracle ceiling | rep=%d cond=%s experiments=%s",
        args.replicate,
        args.condition,
        args.experiments,
    )
    run_all_experiments(
        replicate_id=args.replicate,
        condition=args.condition,
        experiments=args.experiments,
        sccube_dir=args.sccube_dir,
        output_base=args.output_dir,
        c2l_base=args.c2l_dir,
        sace_base=args.sace_dir,
    )


if __name__ == "__main__":
    main()
