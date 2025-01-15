# tests/test_citegeist_simulated.py

import os
import sys
import argparse
import logging
import gc
from datetime import datetime

import numpy as np
import scanpy as sc
import pandas as pd

# Here we import the new CITEgeist code (instead of re-defining large chunks of it).
# Adjust these imports according to your project structure.
import sys
import os

# Add the parent directory (not just model) to the system path
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Now import using the full package path
from model.citegeist_model import CitegeistModel
from model.utils import benchmark_cell_proportions, calculate_expression_metrics, export_anndata_layers

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
    # parser.add_argument('--alpha', type=float, required=True, help='Elastic-net alpha value')
    # parser.add_argument('--gex_lambda', type=float, required=True, help='Regularization lambda for gene expression')
    parser.add_argument('--input_folder', type=str, default='.', help='Folder all requisite samples and ground truth')
    parser.add_argument('--output_folder', type=str, default='citegeist_output', help='Output folder')
    parser.add_argument('--sample_prefix', type=str, default='Wu_rep', help='Prefix to filter sample files')
    parser.add_argument('--profiling_only', action='store_true', default=False, 
                        help='If set, only compute cell-type proportions (no gene expression deconvolution).')
    args = parser.parse_args()

    radius = args.radius
    input_folder = args.input_folder
    output_folder = args.output_folder

    suffix = "FilteredRadiiArrayWinsorCLRDiscrete"

    output_folder = os.path.join(output_folder, f'test_results/radius_{radius}.', suffix + "CITEgeistOutput")

    # Create an output directory
    os.makedirs(output_folder, exist_ok=True)

    # Initialize logging
    log_file = f"Simulated_CITEgeist_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{radius}.log"
    log_path = os.path.join(args.output_folder, log_file)
    logging.basicConfig(
        filename=log_path,
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG
    )
    logging.info(f"Starting CITEgeist run with parameters: radius={radius}")

    # -------------------------------------------------------------------------
    # 1) Load H5AD files (Gene Expression and Antibody Capture)
    #    Example: any .h5ad files that match your naming pattern
    # -------------------------------------------------------------------------
    # Find all unique sample numbers for Wu_rep_{number} pairs
    h5ad_dir = os.path.join(args.input_folder, "h5ad_objects")
    sample_numbers = []
    for f in os.listdir(h5ad_dir):
        if f.startswith(args.sample_prefix):
            # Extract number from Wu_rep_{number}
            number = f.split('_')[2].split('.')[0]  # Handle potential file extensions
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
        
        # Check for log1p transformation
        if np.any(adata_gex.X > 0) and np.all(adata_gex.X == np.round(np.expm1(adata_gex.X))):
            logging.warning("adata_gex.X appears to be in log1p space. Consider applying inverse transformation.")

        
        # Since the GEX data is log1p transformed, we will inverse that transformation.
        adata_gex.X = np.round(np.expm1(adata_gex.X))  # Round after inverse transformation
        # Inverse log1p transformation
        logging.info("Inverse log1p transformation applied to the GEX data.")


        # Our simulated data is pre-split into CITE and GEX since it doesn't have the idiosyncrasies of Visium data.

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
        model.preprocess_gex()
        model.preprocess_antibody()

        # Register Gurobi license
        model.register_gurobi("/ihome/crc/install/gurobi/gurobi1102/linux64/lic/gurobi.lic")

        logging.info(f"Running cell proportion model for {sample_name} ...")



        ##############################################################################
        # 1) Cell Proportion Inference
        ##############################################################################

        # Run cell proportion model
        model.run_cell_proportion_model()

        logging.info(f"Completed cell proportion inference for {sample_name}.")

        # Plot cell proportions (Append cell proportions) 
        model.append_proportions_to_adata()


        

        # Benchmarking Cell Proportions

        st_folder = os.path.join(input_folder, "ST_sim")
    
        number = sample_name.split('_')[2]


        proportions_path = os.path.join(output_folder, f"{sample_name}_cell_prop_results.csv")
        
        test_spots_df = pd.read_csv(proportions_path, index_col = 0).sort_index().sort_index(axis=1)

        spot_composition_df = pd.read_csv(os.path.join(st_folder, f"Wu_ST_{number}_prop.csv"), index_col = 0).sort_index().sort_index(axis=1)
        spot_composition_df = spot_composition_df.iloc[:, :-2]


        # Check if columns match
        if not np.array_equal(test_spots_df.columns, spot_composition_df.columns):
            logging.warning(f"test_spots_df columns: {test_spots_df.columns}, spot_composition_df columns: {spot_composition_df.columns}")
            raise ValueError("ERROR: The column names in the input CSV files do not match or are not in the same order!")

        # Check if rows match
        if not np.array_equal(test_spots_df.index, spot_composition_df.index):
            logging.warning(f"test_spots_df index: {test_spots_df.index}, spot_composition_df index: {spot_composition_df.index}")
            mismatches = [(i, test_spots_df.index[i], spot_composition_df.index[i]) 
                        for i in range(len(test_spots_df.index)) 
                        if test_spots_df.index[i] != spot_composition_df.index[i]]
            logging.warning("Mismatched indices:")
            for idx, test_index, comp_index in mismatches:
                logging.warning(f"Position {idx}: test_spots_df index = {test_index}, spot_composition_df index = {comp_index}")
            

        # Convert DataFrames to numpy arrays
        test_spots_metadata_mtrx = test_spots_df.values
        spot_composition_mtrx = spot_composition_df.values
        column_names = test_spots_df.columns.tolist()  # Get column names for RMSE and MAE

        
        results = benchmark_cell_proportions(test_spots_metadata_mtrx, spot_composition_mtrx, column_names)

        logging.info(f"Cell proportion benchmarking results: {results}")

        # Convert to a DataFrame
        prop_results_df = pd.DataFrame([results])  # Wrap in a list to create a single-row DataFrame
        
        # Save to a CSV file
        prop_results_name = os.path.join(output_folder, f'{sample_name}_cellprop_results_summary_{suffix}_{radius}_.csv')
        
        prop_results_df.to_csv(prop_results_name, index=False)

        print(prop_results_df)
        logging.info(f"Cell proportion results summary: \n{prop_results_df}")


        ##############################################################################
        # 2) Gene Expression Deconvolution
        ##############################################################################

        if args.profiling_only:
            logging.info("Skipping gene-expression deconvolution (profiling_only=True).")
            continue


        logging.info(f"Running gene expression optimization for {sample_name} ...")

        model.run_cell_expression_model(radius=radius, alpha=0.5, lambda_reg_gex=0.001, lambda_prior_weight=1,
                                max_workers=None, checkpoint_interval=100, 
                                output_dir="checkpoints", rerun=True)

        # Plot cell proportions (Append cell proportions) 
        model.append_gex_to_adata()

        prop_gex_adata = model.get_adata()

        print(prop_gex_adata)


        # Benchmarking Cell Gene Expression

        layer_dir = os.path.join(output_folder, f"{sample_name}_{suffix}_{radius}/layers")
    
        export_anndata_layers(prop_gex_adata, layer_dir)

        ground_truth_folder = os.path.join(input_folder,"ST_GEX_sim")

        # Set directories for ground truth and prediction files
        ground_truth_dir = os.path.join(ground_truth_folder,f"sample_{number}","layers")

        # Add these debug statements before calling calculate_expression_metrics
        logging.info(f"Ground truth directory: {ground_truth_dir}")
        logging.info(f"Ground truth files: {os.listdir(ground_truth_dir)}")
        logging.info(f"Layer directory: {layer_dir}")
        logging.info(f"Layer files: {os.listdir(layer_dir)}")

        # Add this debug code before calculate_expression_metrics
        for gt_file in os.listdir(ground_truth_dir):
            if gt_file.endswith("_GT.csv"):
                cell_type = gt_file.replace("_GT.csv", "")
                gt_path = os.path.join(ground_truth_dir, gt_file)
                pred_path = os.path.join(layer_dir, f"{cell_type}_layer.csv")
                
                if os.path.exists(pred_path):
                    gt_df = pd.read_csv(gt_path, index_col=0)
                    pred_df = pd.read_csv(pred_path, index_col=0)
                    logging.info(f"Cell type: {cell_type}")
                    logging.info(f"Ground truth shape: {gt_df.shape}")
                    logging.info(f"Prediction shape: {pred_df.shape}")
                    logging.info(f"Common genes: {len(gt_df.index.intersection(pred_df.index))}")
                    logging.info(f"Common spots: {len(gt_df.columns.intersection(pred_df.columns))}")

        # Calculate RMSE and Normalized RMSE
        metrics = calculate_expression_metrics(ground_truth_dir, layer_dir, normalize = "range")


        average_rmse  = metrics.get('average_rmse')
        median_rmse   = metrics.get('median_rmse')
        average_nrmse = metrics.get('average_nrmse')
        median_nrmse  = metrics.get('median_nrmse')
        average_mae   = metrics.get('average_mae')
        median_mae    = metrics.get('median_mae')
    
        # Create the metrics dictionary
        metrics_values = {
            'Average RMSE': average_rmse,
            'Median RMSE': median_rmse,
            'Average NRMSE': average_nrmse,
            'Median NRMSE': median_nrmse,
            'Average MAE': average_mae,
            'Median MAE': median_mae
        }
    
        # Convert the metrics dictionary to a DataFrame
        gex_results = pd.DataFrame(list(metrics_values.items()), columns=['Metric', 'Value'])

        print(gex_results)

        logging.info(f"Gene expression metrics: \n{gex_results}")

        # Save the DataFrame to a CSV file
        gex_results.to_csv(os.path.join(output_folder,f'{sample_name}_gex_metrics_summary_{suffix}_{radius}.csv'), index=False)

        ##############################################################################
        # 3) Phase 3: WGCNA-based Multimodal Integration
        ##############################################################################
        
        logging.info(f"Running Phase 3 WGCNA-based integration for {sample_name} ...")

        # Run Phase 3 optimization
        model.run_multimodal_phase_3_wgcna(
            max_clusters=30,  # Adjust based on expected number of gene modules
            alpha_gex=1.0,
            alpha_antibody=0.5,
            alpha_spatial=0.0,  # No spatial smoothing for now
            lambda_reg_module=0.1
        )

        # Get updated AnnData with Phase 3 layers
        prop_gex_adata_phase3 = model.get_adata()

        # Export Phase 3 layers
        layer_dir_phase3 = os.path.join(output_folder, f"{sample_name}_{suffix}_{radius}_phase3/layers")
        export_anndata_layers(prop_gex_adata_phase3, layer_dir_phase3)

        # Calculate Phase 3 metrics for cell proportions
        proportions_path_phase3 = os.path.join(output_folder, f"{sample_name}_cell_prop_phase3_results.csv")
        test_spots_df_phase3 = pd.read_csv(proportions_path_phase3, index_col=0).sort_index().sort_index(axis=1)

        # Verify columns match
        if not np.array_equal(test_spots_df_phase3.columns, spot_composition_df.columns):
            logging.warning(f"Phase 3 test_spots_df columns: {test_spots_df_phase3.columns}, spot_composition_df columns: {spot_composition_df.columns}")
            raise ValueError("ERROR: The column names in the Phase 3 CSV files do not match or are not in the same order!")

        # Convert DataFrames to numpy arrays for Phase 3 cell proportions
        test_spots_metadata_mtrx_phase3 = test_spots_df_phase3.values
        results_phase3 = benchmark_cell_proportions(test_spots_metadata_mtrx_phase3, spot_composition_mtrx, column_names)

        logging.info(f"Phase 3 cell proportion benchmarking results: {results_phase3}")

        # Convert to DataFrame and save Phase 3 cell proportion results
        prop_results_df_phase3 = pd.DataFrame([results_phase3])
        prop_results_name_phase3 = os.path.join(output_folder, f'{sample_name}_cellprop_results_summary_phase3_{suffix}_{radius}.csv')
        prop_results_df_phase3.to_csv(prop_results_name_phase3, index=False)

        print("Phase 3 cell proportion results:")
        print(prop_results_df_phase3)
        logging.info(f"Phase 3 cell proportion results summary: \n{prop_results_df_phase3}")

        # Calculate Phase 3 metrics for gene expression
        metrics_phase3 = calculate_expression_metrics(ground_truth_dir, layer_dir_phase3, normalize="range")

        average_rmse_phase3  = metrics_phase3.get('average_rmse')
        median_rmse_phase3   = metrics_phase3.get('median_rmse')
        average_nrmse_phase3 = metrics_phase3.get('average_nrmse')
        median_nrmse_phase3  = metrics_phase3.get('median_nrmse')
        average_mae_phase3   = metrics_phase3.get('average_mae')
        median_mae_phase3    = metrics_phase3.get('median_mae')
    
        # Create the metrics dictionary for Phase 3
        metrics_values_phase3 = {
            'Average RMSE': average_rmse_phase3,
            'Median RMSE': median_rmse_phase3,
            'Average NRMSE': average_nrmse_phase3,
            'Median NRMSE': median_nrmse_phase3,
            'Average MAE': average_mae_phase3,
            'Median MAE': median_mae_phase3
        }
    
        # Convert the Phase 3 metrics dictionary to a DataFrame
        gex_results_phase3 = pd.DataFrame(list(metrics_values_phase3.items()), columns=['Metric', 'Value'])

        print("Phase 3 gene expression metrics:")
        print(gex_results_phase3)
        logging.info(f"Phase 3 gene expression metrics: \n{gex_results_phase3}")

        # Save Phase 3 gene expression metrics
        gex_results_phase3.to_csv(os.path.join(output_folder,f'{sample_name}_gex_metrics_summary_phase3_{suffix}_{radius}.csv'), index=False)

        # Print comparison between original and Phase 3 results
        print("\nComparison of results:")
        print("Cell proportions:")
        comparison_df = pd.concat([
            prop_results_df.add_prefix('Original_'),
            prop_results_df_phase3.add_prefix('Phase3_')
        ], axis=1)
        print(comparison_df)
        logging.info(f"Comparison of cell proportions:\n{comparison_df}")

        print("\nGene expression:")
        gex_comparison = pd.DataFrame({
            'Metric': gex_results['Metric'],
            'Original': gex_results['Value'],
            'Phase3': gex_results_phase3['Value']
        })
        print(gex_comparison)
        logging.info(f"Comparison of gene expression metrics:\n{gex_comparison}")

        # Save comparisons
        comparison_df.to_csv(os.path.join(output_folder,f'{sample_name}_comparison_cellprop_{suffix}_{radius}.csv'), index=False)
        gex_comparison.to_csv(os.path.join(output_folder,f'{sample_name}_comparison_gex_{suffix}_{radius}.csv'), index=False)

        logging.info(f"Finished processing sample: {sample_name}")

        # ---------------------------------------------------------------------
        # Cleanup this sample to free memory
        # ---------------------------------------------------------------------
        del adata_cite
        del adata_gex
        del model
        gc.collect()

    logging.info("All samples processed successfully.")


if __name__ == "__main__":
    main()