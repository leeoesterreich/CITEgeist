"""
Leiden clustering baseline for Module 1-2 discovery comparison.

Runs standard Scanpy Leiden clustering at multiple resolutions on Xenium
protein expression data (both single-cell and pseudo-Visium) and extracts
cluster-enriched marker signatures for head-to-head comparison with
Module 1-2's spatial co-expression discovery.

Usage:
    python leiden_baseline_comparison.py \
        --region 0 \
        --resolution-level spot \
        --output-dir results/discovery_comparison
"""
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[4]
PSEUDOVISIUM_DIR = (
    REPO_ROOT
    / "Benchmarking/xenium_pseudovisium/data_protein_gt/h5ad_objects"
)

sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))
from load_xenium_singlecell import load_xenium_singlecell

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

LEIDEN_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0, 1.5]


def load_data(
    region_id: int, resolution: str
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Load protein expression and spatial coordinates.

    Args:
        region_id: Region index (0-4).
        resolution: 'spot' for pseudo-Visium or 'cell' for single-cell.

    Returns:
        (X, coords, marker_names)
    """
    if resolution == "spot":
        cite_path = PSEUDOVISIUM_DIR / f"Xenium_region_{region_id}_CITE.h5ad"
        adata = sc.read_h5ad(str(cite_path))
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float64)
        coords = adata.obsm["spatial"]
        marker_names = list(adata.var_names)
    else:
        _, adata_protein = load_xenium_singlecell(region_id=region_id)
        X = adata_protein.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float64)
        coords = adata_protein.obsm["spatial"]
        marker_names = list(adata_protein.var_names)

    logger.info(
        f"Loaded region {region_id} ({resolution}): "
        f"{X.shape[0]} observations, {X.shape[1]} markers"
    )
    return X, coords, marker_names


def run_leiden_at_resolution(
    X: np.ndarray,
    marker_names: List[str],
    resolution: float,
    log2fc_threshold: float = 0.5,
    pval_threshold: float = 0.05,
) -> Dict[str, List[str]]:
    """
    Run Leiden clustering and extract cluster marker signatures.

    Args:
        X: Expression matrix (n_obs, n_markers).
        marker_names: Marker names.
        resolution: Leiden resolution parameter.
        log2fc_threshold: Minimum log2 fold change for marker inclusion.
        pval_threshold: Maximum adjusted p-value for marker inclusion.

    Returns:
        Dict mapping cluster_id -> list of enriched marker names.
    """
    adata = sc.AnnData(X)
    adata.var_names = marker_names

    # Standard Scanpy workflow
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=min(X.shape[1] - 1, 20))
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")

    n_clusters = len(adata.obs["leiden"].unique())
    logger.info(f"  Resolution {resolution}: {n_clusters} clusters")

    # Extract enriched markers per cluster
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    cluster_signatures: Dict[str, List[str]] = {}
    for cluster_id in sorted(adata.obs["leiden"].unique()):
        result = sc.get.rank_genes_groups_df(adata, group=cluster_id)
        sig_markers = result[
            (result["logfoldchanges"] > log2fc_threshold)
            & (result["pvals_adj"] < pval_threshold)
        ]["names"].tolist()
        cluster_signatures[f"cluster_{cluster_id}"] = sig_markers
        logger.info(
            f"    Cluster {cluster_id} "
            f"({(adata.obs['leiden'] == cluster_id).sum()} obs): "
            f"{sig_markers}"
        )

    return cluster_signatures


def run_leiden_all_resolutions(
    X: np.ndarray,
    marker_names: List[str],
) -> Dict[str, Dict[str, List[str]]]:
    """
    Run Leiden at all resolutions, return nested dict.

    Returns:
        {resolution_str: {cluster_id: [markers]}}
    """
    results = {}
    for res in LEIDEN_RESOLUTIONS:
        logger.info(f"Running Leiden at resolution {res}...")
        sigs = run_leiden_at_resolution(X, marker_names, res)
        results[str(res)] = sigs
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run Leiden baseline for discovery comparison"
    )
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument(
        "--resolution-level",
        choices=["spot", "cell"],
        required=True,
        help="Data resolution: 'spot' (pseudo-Visium) or 'cell' (single-cell)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(
            REPO_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results/discovery_comparison"
        ),
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    X, coords, marker_names = load_data(args.region, args.resolution_level)

    leiden_results = run_leiden_all_resolutions(X, marker_names)

    # Save results
    output_file = (
        output_dir / f"leiden_region_{args.region}_{args.resolution_level}.json"
    )
    with open(output_file, "w") as f:
        json.dump(
            {
                "region_id": args.region,
                "resolution_level": args.resolution_level,
                "n_observations": X.shape[0],
                "n_markers": X.shape[1],
                "marker_names": marker_names,
                "leiden_resolutions": LEIDEN_RESOLUTIONS,
                "results": leiden_results,
            },
            f,
            indent=2,
        )

    logger.info(f"Saved Leiden results to {output_file}")


if __name__ == "__main__":
    main()
