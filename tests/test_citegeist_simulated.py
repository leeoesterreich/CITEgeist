# tests/test_citegeist_simulated.py

import os
import sys
import argparse
import logging
import gc
import json
from datetime import datetime
from typing import Dict, Any, List, Optional
import time

import numpy as np
import scanpy as sc
import pandas as pd
import scipy.sparse

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
    
    # Check if all metrics are not None or Nan, if they are, print the celltype
    for celltype, metric in metrics.items():
        if metric['RMSE'] is None or metric['NRMSE'] is None or metric['MAE'] is None or np.isnan(metric['RMSE']) or np.isnan(metric['NRMSE']) or np.isnan(metric['MAE']):
            print(f"Cell type {celltype} has None or NaN metrics")

    # Create DataFrame with metrics while excluding None or NaN values
    
    metrics_values = {
        'Pass': [f"Pass {pass_number}" if pass_number else "Unknown"] * 6,
        'Metric': [
            'Average RMSE', 'Median RMSE', 'Average NRMSE',
            'Median NRMSE', 'Average MAE', 'Median MAE'
        ],
        'Value': [
            np.nanmean([m['RMSE'] for m in metrics.values() if m['RMSE'] is not None]),
            np.nanmedian([m['RMSE'] for m in metrics.values() if m['RMSE'] is not None]),
            np.nanmean([m['NRMSE'] for m in metrics.values() if m['NRMSE'] is not None]),
            np.nanmedian([m['NRMSE'] for m in metrics.values() if m['NRMSE'] is not None]),
            np.nanmean([m['MAE'] for m in metrics.values() if m['MAE'] is not None]),
            np.nanmedian([m['MAE'] for m in metrics.values() if m['MAE'] is not None])
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
    parser.add_argument('--lambda_reg', type=float, required=True, 
                       help='Regularization strength for elastic net')
    parser.add_argument('--alpha_elastic', type=float, required=True, 
                       help='Elastic net mixing parameter (0=L2, 1=L1)')
    parser.add_argument('--max_y_change', type=float, required=True, 
                       help='Maximum allowed change in Y values between iterations (0,1)')
    parser.add_argument('--input_folder', type=str, default='.', help='Folder all requisite samples and ground truth')
    parser.add_argument('--output_folder', type=str, default='citegeist_output', help='Output folder')
    parser.add_argument('--sample_prefix', type=str, default='Wu_rep', help='Prefix to filter sample files')
    parser.add_argument('--profiling_only', action='store_true', default=False, 
                        help='If set, only compute cell-type proportions (no gene expression deconvolution).')
    parser.add_argument('--skip_pass2', action='store_true', default=False, 
                        help='If set, skip pass 2 and only run pass 1.')
    args = parser.parse_args()

    radius = args.radius
    lambda_reg = args.lambda_reg
    alpha_elastic = args.alpha_elastic
    max_y_change = args.max_y_change
    variables = f"radius_{radius}_lambda_{lambda_reg}_alpha_{alpha_elastic}_max_y_change_{max_y_change}"

    input_folder = args.input_folder
    output_folder = args.output_folder

    suffix = "FilteredRadiiArrayWinsorCLRDiscreteErrorMinimizing"

    output_folder = os.path.join(output_folder, f'test_results/{variables}.', suffix + "CITEgeistOutput")

    # Create an output directory
    os.makedirs(output_folder, exist_ok=True)

    # Initialize logging
    log_file = f"Simulated_CITEgeist_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{radius}_{lambda_reg}_{alpha_elastic}.log"
    log_path = os.path.join(args.output_folder, log_file)
    logging.basicConfig(
        filename=log_path,
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG
    )
    logging.info(f"Starting CITEgeist run with parameters: radius={radius}, lambda_reg={lambda_reg}, alpha_elastic={alpha_elastic}")

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
        start_time = time.time()  # Start timing for this sample
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

        X = adata_gex.X
        if X is None:
            logging.warning("adata_gex.X is None")
            return
        if scipy.sparse.issparse(X):
            X = X.toarray()  # type: ignore
        if isinstance(X, np.ndarray) and np.any(X < 0):
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

        global_cell_type_proportions_df, finetuned_cell_type_proportions_df = model.run_cell_proportion_model(
            radius=radius,
            tolerance=1e-4,
            max_iterations=20,
            lambda_reg=lambda_reg,
            alpha=alpha_elastic,
            max_workers=None,
            checkpoint_interval=100,
            max_y_change=max_y_change
        )

        
        logging.info(f"Completed cell proportion inference for {sample_name}.")

        # # Plot cell proportions (Append cell proportions) 
        # model.append_proportions_to_adata()

        # # Benchmarking Cell Proportions
        st_folder = os.path.join(input_folder, "ST_sim")

        # proportions_path = os.path.join(output_folder, f"{sample_name}_cell_prop_results.csv")
        # test_spots_df = pd.read_csv(proportions_path, index_col=0).sort_index().sort_index(axis=1)


        spot_composition_df = pd.read_csv(os.path.join(st_folder, f"Wu_ST_{number}_prop.csv"), index_col=0).sort_index().sort_index(axis=1)
        spot_composition_df = spot_composition_df.iloc[:, :-2]

        # Sort indices numerically by spot number
        def sort_spot_indices(df):
            # Extract numbers from spot names and sort
            df.index = pd.Index(df.index, name='spot')
            return df.reindex(sorted(df.index, key=lambda x: int(x.split('_')[1]) if '_' in x else float('inf')))

        spot_composition_df = sort_spot_indices(spot_composition_df)
        
        results_dict = {
            'global': global_cell_type_proportions_df,
            'finetune': finetuned_cell_type_proportions_df
        }

        for key, test_spots_df in results_dict.items():
            # Sort test spots numerically
            test_spots_df = sort_spot_indices(test_spots_df)

            # Sort both DataFrames by index and ensure indices match
            test_spots_df = test_spots_df.sort_index()
            spot_composition_df = spot_composition_df.sort_index()

            print(test_spots_df.index)
            print(spot_composition_df.index)

            # Verify that indices match
            if not np.array_equal(test_spots_df.index, spot_composition_df.index):
                logging.warning(f"test_spots_df indices: {test_spots_df.index}, spot_composition_df indices: {spot_composition_df.index}")
                raise ValueError("ERROR: The row indices in the input CSV files do not match or are not in the same order!")

            # Check if columns match
            if not np.array_equal(test_spots_df.columns, spot_composition_df.columns):
                logging.warning(f"test_spots_df columns: {test_spots_df.columns}, spot_composition_df columns: {spot_composition_df.columns}")
                raise ValueError("ERROR: The column names in the input CSV files do not match or are not in the same order!")

            # Convert DataFrames to numpy arrays
            test_spots_metadata_mtrx = test_spots_df.values
            spot_composition_mtrx = spot_composition_df.values


            column_names = test_spots_df.columns.tolist()

            results = benchmark_cell_proportions(test_spots_metadata_mtrx, spot_composition_mtrx, column_names)
            logging.info(f"{key.capitalize()} Cell proportion benchmarking results: {results}")

            # Save cell proportion results
            prop_results_df = pd.DataFrame([results])
            prop_results_name = os.path.join(output_folder, f'{sample_name}_cellprop_results_summary_{key}_{suffix}_{radius}_.csv')
            
            prop_results_df.to_csv(prop_results_name, index=False)
            logging.info(f"{key.capitalize()} Cell proportion results summary: \n{prop_results_df}")
            print(f"{key.capitalize()} Cell proportion results summary: \n{prop_results_df}")


        model.append_proportions_to_adata(key='finetuned')


        if args.profiling_only:
            logging.info("Skipping gene-expression deconvolution (profiling_only=True).")
            continue

        ##############################################################################
        # 2) Gene Expression Deconvolution - Pass 1
        ##############################################################################
        logging.info(f"Running pass 1 gene expression optimization for {sample_name} ...")

        # Run first pass with weight parameters
        pass1_results = model.run_cell_expression_pass1(
            radius=radius,
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

            
            assert pass1_metrics is not None, "Pass 1 metrics are None"
            
            print(f"Pass 1 metrics:\n{pass1_metrics}")
            print(f"Pass 1 metrics:\n{pass1_metrics}")
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

            
        end_time = time.time()
        runtime = end_time - start_time
        runtime_message = f"Runtime for sample {sample_name}: {runtime:.2f} seconds ({runtime/60:.2f} minutes)"
        print(runtime_message)
        logging.info(runtime_message)

if __name__ == "__main__":
    main()