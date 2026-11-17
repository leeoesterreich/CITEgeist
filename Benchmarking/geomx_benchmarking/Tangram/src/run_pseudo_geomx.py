"""
Run Tangram on pseudo-GeoMx AOI data.

Adapted from xenium_benchmarking/Tangram/src/run_benchmark.py.
Single batch, density_prior='rna_count_based'.
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def find_marker_genes(adata_ref, labels_key="cell_type", n_markers=100):
    """Find top marker genes per cell type via Wilcoxon ranking."""
    logger.info(f"Finding top {n_markers} markers per cell type...")
    adata_work = adata_ref.copy()
    sc.pp.normalize_total(adata_work, target_sum=1e4)
    sc.pp.log1p(adata_work)
    sc.tl.rank_genes_groups(adata_work, groupby=labels_key, method="wilcoxon")

    markers = set()
    for ct in adata_work.obs[labels_key].unique():
        group_markers = adata_work.uns["rank_genes_groups"]["names"][ct][:n_markers]
        markers.update(group_markers.tolist())

    logger.info(f"Found {len(markers)} unique marker genes")
    return list(markers)


def extract_gex_layers(ad_map, adata_ref, adata_sp, labels_key="cell_type"):
    """Extract per-cell-type GEX layers from Tangram mapping."""
    logger.info("Extracting per-cell-type GEX layers...")

    mapping_matrix = ad_map.X
    cell_types = adata_ref.obs[labels_key].values
    unique_cts = np.unique(cell_types)

    ref_expr = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else adata_ref.X

    layer_dfs = {}
    for ct in unique_cts:
        ct_mask = cell_types == ct
        ct_mapping = mapping_matrix[ct_mask, :]
        ct_weights = ct_mapping / (ct_mapping.sum(axis=0, keepdims=True) + 1e-10)
        ct_expr = ref_expr[ct_mask, :]
        layer_expr = ct_weights.T @ ct_expr

        layer_dfs[ct] = pd.DataFrame(
            layer_expr,
            index=adata_sp.obs_names,
            columns=adata_ref.var_names,
        )
        logger.info(f"  {ct}: {layer_dfs[ct].shape}")

    return layer_dfs


def main():
    parser = argparse.ArgumentParser(description="Run Tangram on pseudo-GeoMx")
    parser.add_argument("--input-h5ad", type=str, required=True, help="Path to aoi_gex.h5ad")
    parser.add_argument("--ref-path", type=str, required=True, help="Path to reference.h5ad")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--n-markers", type=int, default=100)
    parser.add_argument("--num-epochs", type=int, default=1000)
    parser.add_argument("--density-prior", type=str, default="rna_count_based")
    parser.add_argument("--device", type=str, default="cuda:0")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info(f"Loading spatial data from {args.input_h5ad}")
    adata_sp = sc.read_h5ad(args.input_h5ad)
    logger.info(f"Spatial shape: {adata_sp.shape}")

    logger.info(f"Loading reference from {args.ref_path}")
    adata_ref = sc.read_h5ad(args.ref_path)
    logger.info(f"Reference shape: {adata_ref.shape}")

    # Find markers and preprocess
    markers = find_marker_genes(adata_ref, n_markers=args.n_markers)

    # Normalize spatial to match reference (normalize_total + log1p)
    sc.pp.normalize_total(adata_sp, target_sum=1e4)
    sc.pp.log1p(adata_sp)

    # Tangram preprocessing
    tg.pp_adatas(adata_ref, adata_sp, genes=markers)
    logger.info(f"After pp_adatas: ref={adata_ref.shape}, spatial={adata_sp.shape}")

    # Map cells to spots
    logger.info(f"Mapping (epochs={args.num_epochs}, device={args.device})...")
    ad_map = tg.map_cells_to_space(
        adata_sc=adata_ref,
        adata_sp=adata_sp,
        mode="cells",
        density_prior=args.density_prior,
        num_epochs=args.num_epochs,
        device=args.device,
    )

    # Extract proportions from mapping
    logger.info("Extracting proportions...")
    tg.project_cell_annotations(ad_map, adata_sp, annotation="cell_type")
    prop_cols = [c for c in adata_sp.obsm["tangram_ct_pred"].columns]
    proportions = adata_sp.obsm["tangram_ct_pred"][prop_cols].copy()
    # Normalize to sum to 1
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    prop_path = output_dir / "predicted_proportions.csv"
    proportions.to_csv(prop_path)
    logger.info(f"Wrote proportions: {prop_path} ({proportions.shape})")

    # Legacy naming
    proportions.to_csv(output_dir / "pseudo_geomx_cell_prop_finetuned_results.csv")

    # Extract GEX layers
    layer_dfs = extract_gex_layers(ad_map, adata_ref, adata_sp)
    layers_dir = output_dir / "layers"
    layers_dir.mkdir(parents=True, exist_ok=True)
    for ct, df in layer_dfs.items():
        layer_path = layers_dir / f"{ct}_layer.csv"
        df.to_csv(layer_path)

    # Save mapping
    ad_map.write(output_dir / "pseudo_geomx_mapping.h5ad")

    logger.info("Done.")


if __name__ == "__main__":
    main()
