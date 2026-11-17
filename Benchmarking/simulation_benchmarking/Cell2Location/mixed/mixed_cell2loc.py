import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

print("Loading packages...")

import gc
import os
import sys

import cell2location
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# Set folders and paths
_BENCH_BASE = str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking")
ref_run_name = f"{_BENCH_BASE}/Cell2Location/reference_signatures"
results_folder = f"{_BENCH_BASE}/Cell2Location/mixed"
input_folder = f"{_BENCH_BASE}/replicates/mixed/h5ad_objects"

print("Loading reference signatures...")
adata_file = f"{ref_run_name}/reference_major.h5ad"
adata_ref = sc.read_h5ad(adata_file)

# Extract pre-computed signatures directly from h5ad (avoids Pyro version mismatch on RegressionModel.load)
inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
    [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
].copy()
inf_aver.columns = adata_ref.uns["mod"]["factor_names"]


# Process a single replicate
def process_replicate(replicate_name):

    # Setting replicate paths and name
    replicate_path = os.path.join(input_folder, replicate_name)
    run_name = os.path.join(results_folder, f"cell2location_map_{replicate_name.split('_')[2]}")

    # Create the output directory
    os.makedirs(run_name, exist_ok=True)

    # Load Visium query dataset
    adata_vis_0 = sc.read_h5ad(replicate_path)

    # Set "SYMBOL" variable
    adata_vis_0.var["SYMBOL"] = adata_vis_0.var.index
    adata_ref.var["SYMBOL"] = adata_ref.var.index

    # Find mitochondria-encoded (MT) genes
    adata_vis_0.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_vis_0.var["SYMBOL"]]

    # Remove MT genes for spatial mapping (keeping their counts in the object)
    adata_vis_0.obsm["MT"] = adata_vis_0[:, adata_vis_0.var["MT_gene"].values].X.toarray()
    adata_vis_0 = adata_vis_0[:, ~adata_vis_0.var["MT_gene"].values]

    # Prepare `adata_vis_0` for cell2location
    # adata_vis_0.X_norm = adata_vis_0.X
    # adata_vis_0.X      = np.expm1(adata_vis_0.X_norm).round()

    # Find shared genes and subset
    intersect = np.intersect1d(adata_vis_0.var_names, inf_aver.index)
    adata_vis_0 = adata_vis_0[:, adata_vis_0.var_names.isin(intersect)].copy()
    inf_aver_shared = inf_aver.loc[inf_aver.index.isin(intersect), :].copy()

    # Setup cell2location
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis_0)

    # Train model
    mod = cell2location.models.Cell2location(
        adata_vis_0,
        cell_state_df=inf_aver_shared,
        N_cells_per_location=5,
        detection_alpha=200,
    )
    mod.train(max_epochs=30000, batch_size=None, train_size=1, accelerator="gpu")

    # Save the loss curve
    loss_curve_path = os.path.join(run_name, "loss_curve.png")
    mod.plot_history(1000)  # Exclude the first 1000 epochs from the plot
    plt.legend(labels=["full data training"])
    plt.savefig(loss_curve_path)
    plt.close()

    # Export posterior and save results
    adata_vis_0 = mod.export_posterior(
        adata_vis_0, sample_kwargs={"num_samples": 3000, "batch_size": mod.adata.n_obs, "accelerator": "gpu"}
    )

    # Compute expected expression per cell type
    expected_dict = mod.module.model.compute_expected_per_cell_type(mod.samples["post_sample_q05"], mod.adata_manager)

    # Add to anndata layers
    for i, n in enumerate(mod.factor_names_):
        adata_vis_0.layers[n] = expected_dict["mu"][i]

    # Add cell abundance quantile to `obs`
    adata_vis_0.obs[adata_vis_0.uns["mod"]["factor_names"]] = adata_vis_0.obsm["q05_cell_abundance_w_sf"]

    # Save proportions
    df = adata_vis_0.obsm["q05_cell_abundance_w_sf"]
    total_abundance = df.sum(axis=1)
    proportions = df.div(total_abundance, axis=0)
    proportions.columns = [col.split("_")[-1] for col in proportions.columns]
    proportions.reset_index(inplace=True)
    proportions.rename(columns={"index": "spot"}, inplace=True)
    proportions.to_csv(f"{run_name}/cell2loc_deconv_predictions.csv", index=False)

    # Export layers to CSV
    layers_output_dir = f"{run_name}/layers"
    os.makedirs(layers_output_dir, exist_ok=True)

    # Saving GEX layer for each cell type
    for layer_name in adata_vis_0.layers.keys():
        layer_data = adata_vis_0.layers[layer_name].toarray()
        df = pd.DataFrame(layer_data, index=adata_vis_0.obs.index, columns=adata_vis_0.var.index)
        df.to_csv(f"{layers_output_dir}/{layer_name}_layer.csv")

    # Save cell2location results
    mod.save(f"{run_name}", overwrite=True)

    print(f"Processing complete for {replicate_name}")


if __name__ == "__main__":
    replicate_name = sys.argv[1]  # Get replicate name from command-line arguments
    print(f"Starting processing for {replicate_name}")
    process_replicate(replicate_name)
