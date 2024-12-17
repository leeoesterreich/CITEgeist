import os
import sys
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import cell2location

# Set folders and paths
ref_run_name = '/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/cell2location/reference_signatures'
results_folder = '/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/cell2location/high_seg'
input_folder = '/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/replicates/high_seg/h5ad_objects'

# Load reference signatures
adata_file = f"{ref_run_name}/reference_major.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[
        f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']
    ]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

# Process a single replicate
def process_replicate(replicate_name):
    replicate_path = os.path.join(input_folder, replicate_name)
    run_name = os.path.join(results_folder, f"cell2location_map_{replicate_name.split('_')[2]}")
    
    # Create output directory
    os.makedirs(run_name, exist_ok=True)

    # Load Visium query dataset
    adata_vis_0 = sc.read_h5ad(replicate_path)

    # Prepare `adata_vis_0`
    adata_vis_0.X_norm = adata_vis_0.X
    adata_vis_0.X = np.expm1(adata_vis_0.X_norm).round()

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
    mod.train(max_epochs=30000, batch_size=None, train_size=1, use_gpu=True)

    # Save the loss curve
    loss_curve_path = os.path.join(run_name, "loss_curve.png")
    mod.plot_history(1000)  # Exclude the first 1000 epochs from the plot
    plt.legend(labels=['full data training'])
    plt.savefig(loss_curve_path)
    plt.close()

    # Save results
    mod.save(f"{run_name}", overwrite=True)

    print(f"Processing complete for {replicate_name}")


if __name__ == "__main__":
    replicate_name = sys.argv[1]  # Get replicate name from command-line arguments
    process_replicate(replicate_name)
