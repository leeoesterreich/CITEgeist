import sys
import os
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
from fast_matrix_market import mmread
import scipy
import gc
import torch
import tangram as tg
import matplotlib.pyplot as plt

def load_wu_data(path, test_mode=False):
    matrix_path = os.path.join(path, "matrix.mtx")
    genes_path = os.path.join(path, "genes.tsv")
    barcodes_path = os.path.join(path, "barcodes.tsv")
    metadata_path = os.path.join(path, "metadata.csv")
    try:
        genes_df = pd.read_csv(genes_path, header=None, sep="\t")
        barcodes_df = pd.read_csv(barcodes_path, header=None)
        barcodes = barcodes_df[0]

        matrix = mmread(matrix_path).T.tocsc().tocsr().astype("float32")
        var_df = pd.DataFrame(index=genes_df[0])
        
        adata = AnnData(
            matrix, obs=pd.DataFrame(index=barcodes), var=var_df
        )

        metadata_df = pd.read_csv(metadata_path, index_col=0)
        adata.obs['patient_ID'] = metadata_df['orig.ident']
        adata.obs['celltype'] = metadata_df['celltype_major']
        adata.obs["tissue"] = "Cancer"
        adata.obs.index.name = None
        if test_mode:
            adata = adata[np.random.choice(adata.shape[0], 10000, replace=False)]

        return adata

    except Exception as e:
        print("Error in load_wu_data function:", e)
        raise

def process_reference_data(wu_adata):
    """Process the Wu reference dataset."""
    sc.pp.filter_cells(wu_adata, min_genes=100)
    sc.pp.filter_genes(wu_adata, min_cells=3)
    sc.pp.normalize_total(wu_adata)
    sc.pp.log1p(wu_adata)
    return wu_adata

def run_tangram_analysis(spatial_path, output_dir, wu_reference_path, device="cuda:0"):
    """Main function to run Tangram analysis."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load spatial data
    print("Loading spatial data...")
    ad_sp = sc.read_h5ad(spatial_path)
    sc.pp.log1p(ad_sp)
    
    # Load and process reference data
    print("Loading and processing reference data...")
    wu_adata = load_wu_data(wu_reference_path)
    ad_sc = process_reference_data(wu_adata)
    
    # Calculate markers
    print("Calculating marker genes...")
    sc.tl.rank_genes_groups(ad_sc, groupby="celltype", use_raw=False)
    markers_df = pd.DataFrame(ad_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
    markers = list(np.unique(markers_df.melt().value.values))
    
    # Pre-process AnnData objects
    print("Pre-processing data for Tangram...")
    tg.pp_adatas(ad_sc, ad_sp, genes=markers)
    
    # Run Tangram mapping
    print("Running Tangram mapping...")
    ad_map = tg.map_cells_to_space(
        ad_sc, 
        ad_sp,
        mode="cells",
        density_prior='rna_count_based',
        num_epochs=500,
        device=device
    )
    
    # Project cell annotations
    print("Projecting cell annotations...")
    tg.project_cell_annotations(ad_map, ad_sp, annotation="celltype")
    
    # Generate and save visualization
    print("Generating visualizations...")
    annotation_list = list(pd.unique(ad_sc.obs['celltype']))
    fig = tg.plot_cell_annotation_sc(ad_sp, annotation_list, perc=0.02, scale_factor=1, spot_size=1)
    plt.savefig(os.path.join(output_dir, 'cell_type_mapping.pdf'))
    plt.close()
    
    # Process cell type predictions
    print("Processing cell type predictions...")
    ad_sp.obsm["node_types"] = ad_sp.obsm["tangram_ct_pred"]
    tangram_ct_pred = ad_sp.obsm["tangram_ct_pred"]
    proportions = tangram_ct_pred.div(tangram_ct_pred.sum(axis=1), axis=0)
    proportions = proportions.fillna(0)
    ad_sp.obsm["proportions"] = proportions
    
    # Save proportions
    proportions.to_csv(os.path.join(output_dir, "cell_type_proportions.csv"))
    
    # Calculate weighted expression matrices
    print("Calculating weighted expression matrices...")
    weighted_expression_per_cell_type = {}
    for cell_type in ad_sp.obsm['tangram_ct_pred'].columns:
        cell_type_probs = ad_sp.obsm['tangram_ct_pred'][cell_type].values[:, None]
        weighted_expression_per_cell_type[cell_type] = ad_sp.X * cell_type_probs
    
    # Save weighted expression matrices
    print("Saving weighted expression matrices...")
    layers_dir = os.path.join(output_dir, "layers")
    os.makedirs(layers_dir, exist_ok=True)
    
    for cell_type, matrix in weighted_expression_per_cell_type.items():
        file_name = f"{cell_type.replace(' ', '_')}_layer.csv"
        file_path = os.path.join(layers_dir, file_name)
        matrix_df = pd.DataFrame(matrix, index=ad_sp.obs.index, columns=ad_sp.var.index.str.upper())
        matrix_df.to_csv(file_path)
    
    print("Analysis complete! Results saved to:", output_dir)

def main():
    parser = argparse.ArgumentParser(description='Run Tangram analysis on spatial transcriptomics data')
    parser.add_argument('--spatial_path', required=True, help='Path to spatial data h5ad file')
    parser.add_argument('--output_dir', required=True, help='Directory to save outputs')
    parser.add_argument('--wu_reference_path', required=True, help='Path to Wu reference data directory')
    parser.add_argument('--device', default='cuda:0', help='Device to run Tangram on (default: cuda:0)')
    
    args = parser.parse_args()
    
    run_tangram_analysis(
        args.spatial_path,
        args.output_dir,
        args.wu_reference_path,
        args.device
    )

if __name__ == "__main__":
    main() 
