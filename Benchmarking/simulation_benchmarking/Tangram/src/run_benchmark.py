"""
Run Tangram deconvolution on Wu BRCA simulation data.

This script runs Tangram spatial mapping using the Wu 2021 BRCA scRNA-seq
reference to predict cell type proportions and per-cell-type gene expression
on simulated Visium spots.

Tangram maps single cells to spatial locations by learning a probabilistic
cell-to-spot assignment that maximises gene expression correlation.

Usage:
    python run_benchmark.py \
        --replicate-id 0 \
        --condition mixed \
        --input-path /path/to/Wu_rep_0_GEX.h5ad \
        --wu-reference-path /path/to/Wu_etal_2021_BRCA_scRNASeq \
        --output-dir /path/to/Tangram/mixed/tangram_map_0
"""

import argparse
import json
import logging
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

try:
    from fast_matrix_market import mmread
except ImportError:
    from scipy.io import mmread

try:
    import tangram as tg
except ImportError:
    tg = None  # type: ignore  # only needed at runtime for GPU mapping

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def load_wu_reference(wu_reference_path: str) -> AnnData:
    """
    Load and preprocess the Wu 2021 BRCA scRNA-seq reference dataset.

    Reads the raw count matrix in Matrix Market format, attaches cell metadata,
    and applies standard preprocessing: filter low-quality cells/genes, then
    normalize library size and log-transform.

    Args:
        wu_reference_path: Path to directory containing matrix.mtx, genes.tsv,
            barcodes.tsv, and metadata.csv.

    Returns:
        Preprocessed AnnData with obs['cell_type'] set from celltype_major.
    """
    logger.info(f"Loading Wu reference from {wu_reference_path}...")

    matrix_path = os.path.join(wu_reference_path, "matrix.mtx")
    genes_path = os.path.join(wu_reference_path, "genes.tsv")
    barcodes_path = os.path.join(wu_reference_path, "barcodes.tsv")
    metadata_path = os.path.join(wu_reference_path, "metadata.csv")

    genes_df = pd.read_csv(genes_path, header=None, sep="\t")
    barcodes_df = pd.read_csv(barcodes_path, header=None)
    barcodes = barcodes_df[0]

    matrix = mmread(matrix_path).T.tocsc().tocsr().astype("float32")
    var_df = pd.DataFrame(index=genes_df[0])

    adata = AnnData(matrix, obs=pd.DataFrame(index=barcodes), var=var_df)
    adata.obs.index.name = None

    metadata_df = pd.read_csv(metadata_path, index_col=0)
    adata.obs["cell_type"] = metadata_df["celltype_major"].values
    adata.obs["subtype"] = metadata_df["subtype"].values

    logger.info(f"  Raw reference shape: {adata.shape}")
    logger.info(f"  Subtypes: {adata.obs['subtype'].value_counts().to_dict()}")

    er_mask = adata.obs["subtype"] == "ER+"
    adata = adata[er_mask].copy()
    logger.info(f"  After ER+ filter: {adata.shape}")
    logger.info(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

    sc.pp.filter_cells(adata, min_genes=10)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    logger.info(f"  Preprocessed reference shape: {adata.shape}")
    return adata


def find_marker_genes(
    adata_ref: AnnData,
    labels_key: str = "cell_type",
    n_markers: int = 100,
) -> list:
    """
    Find marker genes for each cell type using the Wilcoxon rank-sum test.

    Args:
        adata_ref: Reference AnnData with cell type labels.
        labels_key: obs column containing cell type labels.
        n_markers: Number of top markers per cell type.

    Returns:
        List of unique marker gene names (union across all cell types).
    """
    logger.info(f"Finding {n_markers} marker genes per cell type...")

    adata = adata_ref.copy()
    sc.tl.rank_genes_groups(adata, groupby=labels_key, method="wilcoxon")

    markers = set()
    for ct in adata.obs[labels_key].unique():
        ct_markers = sc.get.rank_genes_groups_df(adata, group=ct)
        top = ct_markers.head(n_markers)["names"].tolist()
        markers.update(top)
        logger.info(f"  {ct}: {len(top)} markers")

    marker_list = list(markers)
    logger.info(f"Total unique markers: {len(marker_list)}")
    return marker_list


def preprocess_for_tangram(
    adata_ref: AnnData,
    adata_sp: AnnData,
    markers: list,
    labels_key: str = "cell_type",
):
    """
    Preprocess reference and spatial data for Tangram.

    Restricts both objects to the intersection of marker genes present in
    both datasets, then calls tg.pp_adatas to set up training genes.

    Args:
        adata_ref: Reference AnnData.
        adata_sp: Spatial AnnData.
        markers: List of marker gene names.
        labels_key: obs column containing cell type labels.

    Returns:
        Tuple of (preprocessed adata_ref, preprocessed adata_sp).
    """
    logger.info("Preprocessing data for Tangram...")

    shared_genes = np.intersect1d(
        np.intersect1d(adata_ref.var_names, adata_sp.var_names),
        markers,
    )
    logger.info(f"  Shared marker genes: {len(shared_genes)}")

    if tg is not None:
        tg.pp_adatas(adata_ref, adata_sp, genes=list(shared_genes))
    else:
        # Manual fallback: subset both objects to shared genes (mirrors tg.pp_adatas)
        adata_ref._inplace_subset_var(shared_genes)
        adata_sp._inplace_subset_var(shared_genes)

    logger.info(f"  Reference after preprocessing: {adata_ref.shape}")
    logger.info(f"  Spatial after preprocessing:   {adata_sp.shape}")
    return adata_ref, adata_sp


def run_tangram_mapping(
    adata_ref: AnnData,
    adata_sp: AnnData,
    labels_key: str = "cell_type",
    num_epochs: int = 500,
    device: str = "cuda:0",
):
    """
    Run Tangram cell-to-space mapping.

    Args:
        adata_ref: Preprocessed reference AnnData.
        adata_sp: Preprocessed spatial AnnData.
        labels_key: obs column containing cell type labels.
        num_epochs: Number of training epochs.
        device: PyTorch device string (e.g. 'cuda:0' or 'cpu').

    Returns:
        Tuple of (mapping AnnData, training time in seconds).
    """
    logger.info(f"Running Tangram mapping ({num_epochs} epochs on {device})...")
    start_time = time.time()

    ad_map = tg.map_cells_to_space(
        adata_sc=adata_ref,
        adata_sp=adata_sp,
        mode="cells",
        density_prior="rna_count_based",
        num_epochs=num_epochs,
        device=device,
    )

    train_time = time.time() - start_time
    logger.info(f"Mapping completed in {train_time:.1f}s ({train_time / 60:.1f} min)")
    return ad_map, train_time


def extract_proportions(
    ad_map: AnnData,
    adata_ref: AnnData,
    adata_sp: AnnData,
    labels_key: str = "cell_type",
) -> pd.DataFrame:
    """
    Extract cell type proportions from Tangram mapping.

    Uses tg.project_cell_annotations to populate adata_sp.obsm['tangram_ct_pred'],
    then row-normalises so each spot sums to 1.

    Args:
        ad_map: Mapping AnnData returned by tg.map_cells_to_space.
        adata_ref: Reference AnnData.
        adata_sp: Spatial AnnData (modified in-place).
        labels_key: obs column containing cell type labels.

    Returns:
        DataFrame of shape (n_spots, n_cell_types), rows sum to 1.
    """
    logger.info("Extracting cell type proportions...")

    tg.project_cell_annotations(ad_map, adata_sp, annotation=labels_key)

    if "tangram_ct_pred" in adata_sp.obsm:
        proportions = pd.DataFrame(
            adata_sp.obsm["tangram_ct_pred"],
            index=adata_sp.obs_names,
            columns=adata_sp.uns.get("tangram_ct_pred_names"),
        )
    else:
        logger.info("  Falling back to mapping matrix for proportions...")
        cell_types = adata_ref.obs[labels_key]
        mapping_matrix = ad_map.X
        proportions = pd.DataFrame(index=adata_sp.obs_names, columns=cell_types.unique())
        for ct in cell_types.unique():
            ct_mask = (cell_types == ct).values
            proportions[ct] = mapping_matrix[ct_mask, :].sum(axis=0)

    proportions = proportions.div(proportions.sum(axis=1), axis=0).fillna(0)
    logger.info(f"  Proportions shape: {proportions.shape} " f"(row mean: {proportions.sum(axis=1).mean():.4f})")
    return proportions


def extract_gex_layers(
    ad_map: AnnData,
    adata_ref: AnnData,
    adata_sp: AnnData,
    labels_key: str = "cell_type",
) -> dict:
    """
    Extract per-cell-type gene expression layers using the Tangram mapping matrix.

    For each cell type, computes a weighted average of reference cell expressions
    where weights are the cells' mapping probabilities to each spot, normalised
    per spot so that within-type weights sum to 1.

    ad_map.X has shape (n_cells, n_spots).

    Args:
        ad_map: Mapping AnnData returned by tg.map_cells_to_space.
        adata_ref: Reference AnnData (cells × genes after pp_adatas).
        adata_sp: Spatial AnnData (spots).
        labels_key: obs column containing cell type labels.

    Returns:
        Dict mapping cell type name → DataFrame of shape (n_spots, n_genes).
    """
    logger.info("Extracting per-cell-type gene expression layers...")

    mapping_matrix = ad_map.X  # (n_cells, n_spots)
    cell_types = adata_ref.obs[labels_key].values
    unique_cts = np.unique(cell_types)

    ref_expr = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else adata_ref.X

    layer_dfs = {}
    for ct in unique_cts:
        logger.info(f"  Processing {ct}...")
        ct_mask = cell_types == ct
        ct_mapping = mapping_matrix[ct_mask, :]  # (n_cells_ct, n_spots)
        ct_weights = ct_mapping / (ct_mapping.sum(axis=0, keepdims=True) + 1e-10)
        ct_expr = ref_expr[ct_mask, :]  # (n_cells_ct, n_genes)
        layer_expr = ct_weights.T @ ct_expr  # (n_spots, n_genes)
        layer_dfs[ct] = pd.DataFrame(
            layer_expr,
            index=adata_sp.obs_names,
            columns=adata_ref.var_names,
        )

    return layer_dfs


def main():
    parser = argparse.ArgumentParser(description="Run Tangram deconvolution on Wu BRCA simulation data")
    parser.add_argument("--replicate-id", type=int, required=True, help="Replicate index (0–4)")
    parser.add_argument(
        "--condition", type=str, required=True, choices=["mixed", "high_seg"], help="Simulation condition"
    )
    parser.add_argument("--input-path", type=str, required=True, help="Path to spatial h5ad file")
    parser.add_argument("--wu-reference-path", type=str, required=True, help="Path to Wu 2021 MTX directory")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory for this replicate")
    parser.add_argument("--n-markers", type=int, default=100, help="Top marker genes per cell type (default: 100)")
    parser.add_argument("--num-epochs", type=int, default=500, help="Tangram training epochs (default: 500)")
    parser.add_argument("--no-gpu", action="store_true", help="Disable GPU; fall back to CPU")
    args = parser.parse_args()

    device = "cpu" if args.no_gpu else "cuda:0"
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"Tangram Simulation Benchmark — {args.condition} rep {args.replicate_id}")
    logger.info("=" * 60)
    logger.info(f"Input:     {args.input_path}")
    logger.info(f"Reference: {args.wu_reference_path}")
    logger.info(f"Output:    {output_dir}")
    logger.info(f"Device:    {device}")

    # Spatial data: raw counts → normalize library size → log-transform
    logger.info("\nLoading spatial data...")
    adata_sp = sc.read_h5ad(args.input_path)
    logger.info(f"  Spatial shape: {adata_sp.shape}")
    sc.pp.normalize_total(adata_sp)
    sc.pp.log1p(adata_sp)
    logger.info("  Applied normalize_total + log1p to spatial data")

    # Reference
    adata_ref = load_wu_reference(args.wu_reference_path)

    # Marker genes + preprocessing
    markers = find_marker_genes(adata_ref, n_markers=args.n_markers)
    adata_ref, adata_sp = preprocess_for_tangram(adata_ref, adata_sp, markers)

    # Tangram mapping
    ad_map, train_time = run_tangram_mapping(adata_ref, adata_sp, num_epochs=args.num_epochs, device=device)

    # Proportions
    proportions = extract_proportions(ad_map, adata_ref, adata_sp)
    prop_path = output_dir / "cell_type_proportions.csv"
    proportions.to_csv(prop_path)
    logger.info(f"Saved proportions to {prop_path}")

    # GEX layers — write to local /tmp first to avoid NFS partial-write corruption
    import shutil
    import tempfile

    layer_dfs = extract_gex_layers(ad_map, adata_ref, adata_sp)
    layers_dir = output_dir / "layers"
    layers_dir.mkdir(exist_ok=True)
    for ct, layer_df in layer_dfs.items():
        layer_path = layers_dir / f"{ct.replace(' ', '_')}_layer.csv"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, dir="/tmp") as tmp:
            tmp_path = Path(tmp.name)
        layer_df.to_csv(tmp_path)
        shutil.move(str(tmp_path), str(layer_path))
        logger.info(f"  Saved {ct}: {layer_df.shape}")

    # Run stats
    stats = {
        "replicate_id": args.replicate_id,
        "condition": args.condition,
        "n_spots": adata_sp.n_obs,
        "n_genes": adata_sp.n_vars,
        "n_markers": len(markers),
        "n_cell_types": len(proportions.columns),
        "cell_types": proportions.columns.tolist(),
        "train_time_sec": train_time,
        "num_epochs": args.num_epochs,
    }
    with open(output_dir / "run_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    # Mapping matrix
    sample_name = f"Tangram_{args.condition}_rep_{args.replicate_id}"
    ad_map.write(output_dir / f"{sample_name}_mapping.h5ad")

    logger.info("\n" + "=" * 60)
    logger.info("Deconvolution Complete!")
    logger.info("=" * 60)
    logger.info(
        f"  Spots: {adata_sp.n_obs} | " f"Cell types: {len(proportions.columns)} | " f"Train time: {train_time:.0f}s"
    )


if __name__ == "__main__":
    main()
