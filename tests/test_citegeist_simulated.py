# tests/test_citegeist_simulated.py

import os
import sys
import argparse
import logging
import gc
from datetime import datetime
from typing import Dict, Any, List, Optional

import numpy as np
import scanpy as sc
import pandas as pd

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Now import using the full package path
from model.citegeist_model import CitegeistModel
from model.utils import benchmark_cell_proportions, calculate_expression_metrics, export_anndata_layers

def calculate_gex_metrics(ground_truth_dir, layer_dir, pass_number=None):
    """
    Calculate gene expression metrics and format them into a DataFrame.
    
    Args:
        ground_truth_dir (str): Directory containing ground truth files
        layer_dir (str): Directory containing prediction layers
        pass_number (int, optional): Pass number for logging
        
    Returns:
        pd.DataFrame: DataFrame containing metrics
    """
    metrics = calculate_expression_metrics(ground_truth_dir, layer_dir, normalize="range", pass_number=pass_number)
    
    # Create DataFrame with metrics
    metrics_values = {
        'Pass': [f"Pass {pass_number}" if pass_number else "Unknown"] * 6,
        'Metric': [
            'Average RMSE', 'Median RMSE', 'Average NRMSE',
            'Median NRMSE', 'Average MAE', 'Median MAE'
        ],
        'Value': [
            np.mean([m['RMSE'] for m in metrics.values()]),
            np.median([m['RMSE'] for m in metrics.values()]),
            np.mean([m['NRMSE'] for m in metrics.values()]),
            np.median([m['NRMSE'] for m in metrics.values()]),
            np.mean([m['MAE'] for m in metrics.values()]),
            np.median([m['MAE'] for m in metrics.values()])
        ]
    }
    return pd.DataFrame(metrics_values)

def calculate_improvements(pass1_metrics, pass2_metrics):
    """
    Calculate improvement percentages between passes.
    
    Args:
        pass1_metrics (pd.DataFrame): Metrics from pass 1
        pass2_metrics (pd.DataFrame): Metrics from pass 2
        
    Returns:
        pd.DataFrame: DataFrame containing improvement percentages
    """
    improvements = {}
    for metric in pass1_metrics['Metric'].unique():
        pass1_value = pass1_metrics[pass1_metrics['Metric'] == metric]['Value'].values[0]
        pass2_value = pass2_metrics[pass2_metrics['Metric'] == metric]['Value'].values[0]
        improvement = ((pass1_value - pass2_value) / pass1_value) * 100
        improvements[metric] = improvement
    
    return pd.DataFrame({
        'Metric': list(improvements.keys()),
        'Improvement_Percentage': list(improvements.values())
    })

##############################################################################
# Example cell-type profile dictionary for demonstration (adjust as needed).
##############################################################################
cell_type_profiles = {
    "B-cells": {
        "Major": ["B-cells_Protein_1", "B-cells_Protein_2"]
    },
    "CAFs": {
        "Major": ["CAFs_Protein_1", "CAFs_Protein_2"]
    },
    "Cancer Epithelial": {
        "Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]
    },
    "Endothelial": {
        "Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]
    },
    "Myeloid": {
        "Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]
    },
    "Normal Epithelial": {
        "Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]
    },
    "PVL": {
        "Major": ["PVL_Protein_1", "PVL_Protein_2"]
    },
    "Plasmablasts": {
        "Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]
    },
    "T-cells": {
        "Major": ["T-cells_Protein_1", "T-cells_Protein_2"]
    }
}

def main():
    """
    Example script to run CITEgeist on simulated data using the newer model-based
    implementation. This script demonstrates how to:
      1) Read in parameter arguments
      2) Load data
      3) Subset and prepare antibody capture
      4) Map antibodies to cell-type profiles
      5) Run the CITEgeist model to obtain cell-type proportions and gene-expression deconvolution
      6) Save the outputs
    """

    parser = argparse.ArgumentParser(description='Run CITEgeist on simulated data.')
    parser.add_argument('--radius', type=float, required=True, help='Radius for neighbor detection')
    parser.add_argument('--lambda_prior_weight', type=float, default=0.5, help='Weight for prior guidance in pass 2')
    parser.add_argument('--input_folder', type=str, default='.', help='Folder all requisite samples and ground truth')
    parser.add_argument('--output_folder', type=str, default='citegeist_output', help='Output folder')
    parser.add_argument('--sample_prefix', type=str, default='Wu_rep', help='Prefix to filter sample files')
    parser.add_argument('--profiling_only', action='store_true', default=False, 
                        help='If set, only compute cell-type proportions (no gene expression deconvolution).')
    args = parser.parse_args()

    radius = args.radius
    lambda_prior_weight = args.lambda_prior_weight
    
    variables = f"radius_{radius}_lambda_prior_weight_{lambda_prior_weight}"

    input_folder = args.input_folder
    output_folder = args.output_folder

    suffix = "FilteredRadiiArrayWinsorCLRDiscreteErrorMinimizing"

    output_folder = os.path.join(output_folder, f'test_results/{variables}.', suffix + "CITEgeistOutput")

    # Create an output directory
    os.makedirs(output_folder, exist_ok=True)

    # Initialize logging
    log_file = f"Simulated_CITEgeist_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{radius}_{lambda_prior_weight}.log"
    log_path = os.path.join(args.output_folder, log_file)
    logging.basicConfig(
        filename=log_path,
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG
    )
    logging.info(f"Starting CITEgeist run with parameters: radius={radius}, lambda_prior_weight={lambda_prior_weight}")

    # Find all unique sample numbers for Wu_rep_{number} pairs
    h5ad_dir = os.path.join(args.input_folder, "h5ad_objects")
    sample_numbers = []
    for f in os.listdir(h5ad_dir):
        if f.startswith(args.sample_prefix):
            number = f.split('_')[2].split('.')[0]
            if number not in sample_numbers:
                sample_numbers.append(number)
    
    if len(sample_numbers) == 0:
        logging.error(f"No samples found in {h5ad_dir} matching prefix: {args.sample_prefix}")
        print(f"No samples found in {h5ad_dir} matching prefix: {args.sample_prefix}")
        sys.exit(1)

    for number in sample_numbers:
        sample_name = f"{args.sample_prefix}_{number}"
        logging.info(f"Processing sample: {sample_name}")

        adata_cite_path = os.path.join(args.input_folder, "h5ad_objects", f"{sample_name}_CITE.h5ad")
        adata_gex_path = os.path.join(args.input_folder, "h5ad_objects", f"{sample_name}_GEX.h5ad")
        
        # Verify files exist before loading
        if not (os.path.exists(adata_cite_path) and os.path.exists(adata_gex_path)):
            logging.error(f"Missing required files for {sample_name}")
            continue
            
        adata_cite = sc.read_h5ad(adata_cite_path)
        logging.info(f"Loaded {adata_cite_path} with shape {adata_cite.shape}")
        
        adata_gex = sc.read_h5ad(adata_gex_path)
        logging.info(f"Loaded {adata_gex_path} with shape {adata_gex.shape}")

        if not np.issubdtype(adata_gex.X.dtype, np.number):
            logging.warning("adata_gex.X contains non-numeric values. Please check the data.")
        
        if np.any(np.isnan(adata_gex.X)):
            logging.warning("adata_gex.X contains NaN values. Please handle them before proceeding.")
        
        if np.any(adata_gex.X < 0):
            logging.warning("adata_gex.X contains negative values, which is not expected in count data.")
        

        ##############################################################################
        # Initialize the model
        ##############################################################################
        
        model = CitegeistModel(sample_name=sample_name, output_folder=output_folder, 
                               simulation=True, 
                               gene_expression_adata=adata_gex, antibody_capture_adata=adata_cite)
        
        # Load cell profile dictionary
        model.load_cell_profile_dict(cell_type_profiles)

        model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)

        # Preprocess datasets
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()

        # Register Gurobi license
        model.register_gurobi("/ihome/crc/install/gurobi/gurobi1102/linux64/lic/gurobi.lic")

        ##############################################################################
        # 1) Cell Proportion Inference
        ##############################################################################
        logging.info(f"Running cell proportion model for {sample_name} ...")
        model.run_cell_proportion_model()
        logging.info(f"Completed cell proportion inference for {sample_name}.")

        # Plot cell proportions (Append cell proportions) 
        model.append_proportions_to_adata()

        # Benchmarking Cell Proportions
        st_folder = os.path.join(input_folder, "ST_sim")
        proportions_path = os.path.join(output_folder, f"{sample_name}_cell_prop_results.csv")
        test_spots_df = pd.read_csv(proportions_path, index_col=0).sort_index().sort_index(axis=1)
        spot_composition_df = pd.read_csv(os.path.join(st_folder, f"Wu_ST_{number}_prop.csv"), index_col=0).sort_index().sort_index(axis=1)
        spot_composition_df = spot_composition_df.iloc[:, :-2]

        # Check if columns match
        if not np.array_equal(test_spots_df.columns, spot_composition_df.columns):
            logging.warning(f"test_spots_df columns: {test_spots_df.columns}, spot_composition_df columns: {spot_composition_df.columns}")
            raise ValueError("ERROR: The column names in the input CSV files do not match or are not in the same order!")

        # Convert DataFrames to numpy arrays
        test_spots_metadata_mtrx = test_spots_df.values
        spot_composition_mtrx = spot_composition_df.values
        column_names = test_spots_df.columns.tolist()

        results = benchmark_cell_proportions(test_spots_metadata_mtrx, spot_composition_mtrx, column_names)
        logging.info(f"Cell proportion benchmarking results: {results}")

        # Save cell proportion results
        prop_results_df = pd.DataFrame([results])
        prop_results_name = os.path.join(output_folder, f'{sample_name}_cellprop_results_summary_{suffix}_{radius}_.csv')
        
        
        prop_results_df.to_csv(prop_results_name, index=False)
        logging.info(f"Cell proportion results summary: \n{prop_results_df}")
        print(f"Cell proportion results summary: \n{prop_results_df}")


        if args.profiling_only:
            logging.info("Skipping gene-expression deconvolution (profiling_only=True).")
            continue

        ##############################################################################
        # 2) Gene Expression Deconvolution - Pass 1
        ##############################################################################
        logging.info(f"Running pass 1 gene expression optimization for {sample_name} ...")

        # Run first pass
        pass1_results = model.run_cell_expression_pass1(
            radius=radius, 
            alpha=0.5, 
            lambda_reg_gex=0.001,
            max_workers=None, 
            checkpoint_interval=100, 
            output_dir="checkpoints", 
            rerun=True
        )

        # Calculate pass 1 metrics
        ground_truth_folder = os.path.join(input_folder, "ST_GEX_sim")
        ground_truth_dir = os.path.join(ground_truth_folder, f"sample_{number}", "layers")
        layer_dir_pass1 = os.path.join(output_folder, f"{sample_name}_pass1/layers")

        if os.path.exists(ground_truth_dir):
            logging.info("Calculating metrics for pass 1...")
            pass1_metrics = calculate_gex_metrics(ground_truth_dir, layer_dir_pass1, pass_number=1)
            logging.info(f"Pass 1 metrics:\n{pass1_metrics}")
            
            # Save pass 1 metrics
            metrics_path_pass1 = os.path.join(output_folder, f"{sample_name}_gex_metrics_pass1.csv")
            pass1_metrics.to_csv(metrics_path_pass1, index=False)

            ##############################################################################
            # 3) Compute Prior and Run Pass 2
            ##############################################################################
            logging.info("Computing prior from pass 1 results...")
            prior_info = model.compute_expression_prior(
                spotwise_profiles_pass1=pass1_results['spotwise_profiles'],
                cell_type_numbers_array=model.results['cell_prop'].values
            )

            # Run pass 2 with prior guidance
            pass2_results = model.run_cell_expression_pass2(
                global_prior=prior_info['global_prior'],
                dimensions=pass1_results['dimensions'],
                radius=radius,
                alpha=0.5,
                lambda_reg_gex=0.001,
                lambda_prior_weight=lambda_prior_weight,
                max_workers=None,
                checkpoint_interval=100,
                output_dir="checkpoints",
                rerun=True
            )

            # Calculate pass 2 metrics
            layer_dir_pass2 = os.path.join(output_folder, f"{sample_name}_pass2/layers")

            pass2_metrics = calculate_gex_metrics(ground_truth_dir, layer_dir_pass2, pass_number=2)
            logging.info(f"Pass 2 metrics:\n{pass2_metrics}")
            
            # Save pass 2 metrics
            metrics_path_pass2 = os.path.join(output_folder, f"{sample_name}_gex_metrics_pass2.csv")
            pass2_metrics.to_csv(metrics_path_pass2, index=False)

            # Calculate and save improvements
            improvements_df = calculate_improvements(pass1_metrics, pass2_metrics)
            logging.info("Improvements from pass 1 to pass 2:")
            for _, row in improvements_df.iterrows():
                logging.info(f"{row['Metric']}: {row['Improvement_Percentage']:.2f}% improvement")

            # Save improvements
            improvements_df.to_csv(
                os.path.join(output_folder, f"{sample_name}_gex_improvements.csv"),
                index=False
            )

            # Save combined metrics
            pd.concat([pass1_metrics, pass2_metrics]).to_csv(
                os.path.join(output_folder, f"{sample_name}_gex_metrics_combined.csv"),
                index=False
            )

            ##############################################################################
            # 4) Run Phase 3 WGCNA-based Optimization
            ##############################################################################
            logging.info("Running Phase 3 WGCNA-based optimization...")

            # Compute spatial adjacency matrix if needed
            spatial_coords = model.gene_expression_adata.obsm.get('spatial', None)
            spatial_adjacency = None
            if spatial_coords is not None:
                from scipy.spatial.distance import cdist
                distances = cdist(spatial_coords, spatial_coords)
                spatial_adjacency = np.exp(-distances / radius)  # Gaussian kernel
                spatial_adjacency[distances > radius] = 0  # Threshold by radius

            # Run Phase 3
            model.run_multimodal_phase_3_wgcna(
                max_clusters=30,  # Number of gene modules
                alpha_gex=1.0,    # Weight for gene expression term
                alpha_antibody=1.0,  # Weight for antibody term
                alpha_spatial=0.5 if spatial_adjacency is not None else 0.0,  # Weight for spatial term
                lambda_reg_module=0.1,  # Regularization for module usage
                spatial_adjacency=spatial_adjacency
            )

            # Calculate Phase 3 metrics for cell proportions
            phase3_proportions_path = os.path.join(output_folder, f"{sample_name}_cell_prop_phase3_results.csv")
            phase3_spots_df = pd.read_csv(phase3_proportions_path, index_col=0).sort_index().sort_index(axis=1)
            
            # Convert to numpy arrays for benchmarking
            phase3_spots_mtrx = phase3_spots_df.values
            phase3_prop_results = benchmark_cell_proportions(phase3_spots_mtrx, spot_composition_mtrx, column_names)
            logging.info(f"Phase 3 cell proportion benchmarking results: {phase3_prop_results}")

            # Save Phase 3 cell proportion results
            phase3_prop_results_df = pd.DataFrame([phase3_prop_results])
            phase3_prop_results_name = os.path.join(output_folder, f'{sample_name}_cellprop_results_phase3_summary.csv')
            phase3_prop_results_df.to_csv(phase3_prop_results_name, index=False)

            # Calculate Phase 3 metrics for gene expression
            layer_dir_phase3 = os.path.join(output_folder, f"{sample_name}_phase3/layers")
            phase3_metrics = calculate_gex_metrics(ground_truth_dir, layer_dir_phase3, pass_number=3)
            logging.info(f"Phase 3 metrics:\n{phase3_metrics}")
            
            # Save Phase 3 metrics
            metrics_path_phase3 = os.path.join(output_folder, f"{sample_name}_gex_metrics_phase3.csv")
            phase3_metrics.to_csv(metrics_path_phase3, index=False)

            # Calculate improvements from Pass 2 to Phase 3
            phase3_improvements_df = calculate_improvements(pass2_metrics, phase3_metrics)
            logging.info("Improvements from Pass 2 to Phase 3:")
            for _, row in phase3_improvements_df.iterrows():
                logging.info(f"{row['Metric']}: {row['Improvement_Percentage']:.2f}% improvement")

            # Save Phase 3 improvements
            phase3_improvements_df.to_csv(
                os.path.join(output_folder, f"{sample_name}_gex_improvements_phase3.csv"),
                index=False
            )

            # Save all metrics combined
            pd.concat([pass1_metrics, pass2_metrics, phase3_metrics]).to_csv(
                os.path.join(output_folder, f"{sample_name}_gex_metrics_all_phases.csv"),
                index=False
            )

        else:
            logging.warning(f"Ground truth directory not found: {ground_truth_dir}")
            logging.warning("Skipping metric calculations and subsequent phases.")

        # Clean up
        gc.collect()

    logging.info("All samples processed successfully.")

if __name__ == "__main__":
    main()