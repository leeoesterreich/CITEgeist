"""
Run Cell2Location on pseudo-GeoMx AOI data.

Adapted from xenium_benchmarking/Cell2Location/src/run_benchmark.py.
Key change: n_cells_per_location=50 (AOIs have ~50-500 cells vs ~5-10 for Visium).

GOTCHAS:
  - varm.copy() returns zeros → use .values
  - compute_expected_per_cell_type returns (spots,1) → use proportion × inf_aver
"""

import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run Cell2Location on pseudo-GeoMx")
    parser.add_argument("--input-h5ad", type=str, required=True, help="Path to aoi_gex.h5ad")
    parser.add_argument("--ref-dir", type=str, required=True, help="Path to C2L trained reference model dir")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--n-cells-per-location", type=int, default=50)
    parser.add_argument("--detection-alpha", type=int, default=200)
    parser.add_argument("--max-epochs", type=int, default=30000)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Import cell2location after args parsing to get faster --help
    from cell2location.models import Cell2location

    # Load spatial data
    logger.info(f"Loading spatial data from {args.input_h5ad}")
    adata_vis = sc.read_h5ad(args.input_h5ad)
    logger.info(f"Spatial shape: {adata_vis.shape}")

    # Load trained reference model
    ref_dir = Path(args.ref_dir)
    logger.info(f"Loading trained reference from {ref_dir}")

    # Load cell type signatures from trained model
    inf_aver = pd.read_csv(ref_dir / "cell_type_signatures.csv", index_col=0)
    logger.info(f"Reference signatures: {inf_aver.shape}")

    # Preprocess spatial: filter to shared genes, remove MT genes
    shared_genes = [g for g in adata_vis.var_names if g in inf_aver.index]
    logger.info(f"Shared genes: {len(shared_genes)} / {adata_vis.n_vars}")
    adata_vis = adata_vis[:, shared_genes].copy()
    inf_aver = inf_aver.loc[shared_genes]

    # Remove MT genes
    mt_genes = [g for g in adata_vis.var_names if g.startswith("MT-")]
    if mt_genes:
        non_mt = [g for g in adata_vis.var_names if g not in mt_genes]
        adata_vis = adata_vis[:, non_mt].copy()
        inf_aver = inf_aver.loc[non_mt]
        logger.info(f"Removed {len(mt_genes)} MT genes, {adata_vis.n_vars} remaining")

    # Setup and train Cell2Location model
    logger.info("Setting up Cell2Location model...")
    Cell2location.setup_anndata(adata=adata_vis)

    mod = Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=args.n_cells_per_location,
        detection_alpha=args.detection_alpha,
    )

    logger.info(f"Training (max_epochs={args.max_epochs})...")
    mod.train(
        max_epochs=args.max_epochs,
        batch_size=None,
        train_size=1,
        accelerator="gpu",
    )

    # Save training history
    plt.figure()
    mod.plot_history(iter_start=1000)
    plt.savefig(output_dir / "training_history.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Export posterior
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs},
    )

    # Save model
    mod.save(str(output_dir / "model"), overwrite=True)

    # --- Extract proportions ---
    logger.info("Extracting proportions...")
    q05_key = "q05_cell_abundance_w_sf"
    if q05_key in adata_vis.obsm:
        abundance = adata_vis.obsm[q05_key]
    else:
        available = list(adata_vis.obsm.keys())
        raise KeyError(f"'{q05_key}' not found in obsm. Available: {available}")

    if isinstance(abundance, pd.DataFrame):
        proportions = abundance.div(abundance.sum(axis=1), axis=0)
    else:
        proportions = pd.DataFrame(
            abundance / abundance.sum(axis=1, keepdims=True),
            index=adata_vis.obs_names,
            columns=inf_aver.columns,
        )

    prop_path = output_dir / "predicted_proportions.csv"
    proportions.to_csv(prop_path)
    logger.info(f"Wrote proportions: {prop_path} ({proportions.shape})")

    # Legacy naming
    proportions.to_csv(output_dir / "pseudo_geomx_cell_prop_finetuned_results.csv")

    # --- Extract per-cell-type GEX layers ---
    logger.info("Extracting per-cell-type GEX layers...")
    layers_dir = output_dir / "layers"
    layers_dir.mkdir(parents=True, exist_ok=True)

    # Use proportion × inf_aver outer product (NOT compute_expected_per_cell_type)
    # inf_aver: (genes, cell_types) — extract via .values (NOT .copy() which returns zeros)
    inf_aver_vals = adata_vis.varm["q05_per_cluster_mu_fg"].values  # (genes, n_types)
    cell_types = inf_aver.columns.tolist()

    for i, ct in enumerate(cell_types):
        # Per-spot expected expression for this cell type:
        # abundance_ct * gene_signature_ct (outer product per spot)
        ct_abundance = abundance.iloc[:, i].values if isinstance(abundance, pd.DataFrame) else abundance[:, i]
        ct_signature = inf_aver_vals[:, i]  # (genes,)
        layer_data = np.outer(ct_abundance, ct_signature)  # (spots, genes)

        df = pd.DataFrame(layer_data, index=adata_vis.obs_names, columns=adata_vis.var_names)
        layer_path = layers_dir / f"{ct}_layer.csv"
        df.to_csv(layer_path)
        logger.info(f"  {ct}: {df.shape}")

    # Save annotated AnnData
    adata_vis.write(output_dir / "pseudo_geomx_annotated.h5ad")

    logger.info("Done.")


if __name__ == "__main__":
    main()
