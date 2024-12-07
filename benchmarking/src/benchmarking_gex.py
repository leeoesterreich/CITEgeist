import os
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
import argparse

def calculate_rmse(ground_truth_dir, predictions_dir, normalize='range'):
    # Dictionary to store RMSE, NRMSE and MAE per cell type
    metrics_per_cell_type = {}

    # Iterate through each ground-truth file
    for gt_filename in os.listdir(ground_truth_dir):
        if gt_filename.endswith("_GT.csv"):
            # Determine the cell type name by removing the "_GT.csv" suffix
            cell_type = gt_filename.replace("_GT.csv", "")
            
            # Construct the file paths for ground truth and prediction
            gt_filepath   = os.path.join(ground_truth_dir, gt_filename)
            pred_filepath = os.path.join(predictions_dir, f"{cell_type}_layer.csv")
            
            # Check if the prediction file exists
            if not os.path.exists(pred_filepath):
                print(f"Prediction file for {cell_type} not found. Skipping.")
                continue
            
            # Load the ground-truth and prediction data
            gt_df   = pd.read_csv(gt_filepath, index_col=0)
            pred_df = pd.read_csv(pred_filepath, index_col=0)
            
            # Find common genes (index) between GT and predictions
            common_genes = gt_df.index.intersection(pred_df.index)

            # Subset both DataFrames to include only common genes
            gt_df   = gt_df.loc[common_genes]
            pred_df = pred_df.loc[common_genes]

            # Ensure the columns (spots) are in the same order
            common_spots = gt_df.columns.intersection(pred_df.columns)
            gt_df   = gt_df[common_spots]
            pred_df = pred_df[common_spots]

            # Check for any empty DataFrames after filtering
            if gt_df.empty or pred_df.empty:
                print(f"No common genes or spots for {cell_type}. Skipping.")
                continue
            
            # Log1p normalization
            gt_df   = np.log1p(gt_df)    # Apply log1p to ground truth
            pred_df = np.log1p(pred_df)  # Apply log1p to predictions

            # Calculate RMSE for this cell type
            mse  = mean_squared_error(gt_df.values, pred_df.values)
            rmse = np.sqrt(mse)
            mae  = mean_absolute_error(gt_df.values, pred_df.values)
        
            # Calculate NRMSE
            if normalize == 'range':
                range_gt = gt_df.values.max() - gt_df.values.min()
                nrmse    = rmse / range_gt if range_gt != 0 else np.nan
            elif normalize == 'mean':
                mean_gt = gt_df.values.mean()
                nrmse   = rmse / mean_gt if mean_gt != 0 else np.nan
            else:
                raise ValueError("Normalization type not recognized. Use 'range' or 'mean'.")

            metrics_per_cell_type[cell_type] = {'RMSE': rmse, 'NRMSE': nrmse, 'MAE': mae}
    
    # Calculate overall RMSE and NRMSE statistics
    all_rmse_values  = [metrics['RMSE'] for metrics in metrics_per_cell_type.values()]
    all_nrmse_values = [metrics['NRMSE'] for metrics in metrics_per_cell_type.values()]
    all_mae_values   = [metrics['MAE'] for metrics in metrics_per_cell_type.values()]
    average_rmse = np.mean(all_rmse_values)
    median_rmse  = np.median(all_rmse_values)

    average_nrmse = np.nanmean(all_nrmse_values)
    median_nrmse  = np.nanmedian(all_nrmse_values)
    
    average_mae = np.nanmean(all_mae_values)
    median_mae  = np.nanmedian(all_mae_values)

    # Print RMSE and NRMSE per cell type and overall statistics
    print("RMSE, NRMSE, and MAE per cell type:")
    for cell_type, metrics in metrics_per_cell_type.items():
        print(f"\t{cell_type}: RMSE: {metrics['RMSE']:.4f}, NRMSE: {metrics['NRMSE']:.4f}, MAE; {metrics['MAE']:.4f}")
    
    print("\nOverall RMSE statistics:")
    print(f"\tAverage RMSE: {average_rmse:.4f}")
    print(f"\tMedian RMSE:  {median_rmse:.4f}")

    print("\nOverall NRMSE statistics:")
    print(f"\tAverage NRMSE: {average_nrmse:.4f}")
    print(f"\tMedian NRMSE:  {median_nrmse:.4f}")

    print("\nOverall MAE statistics:")
    print(f"\tAverage MAE: {average_mae:.4f}")
    print(f"\tMedian MAE:  {median_mae:.4f}")

    return { "metrics_per_cell_type" : metrics_per_cell_type, "average_rmse" : average_rmse, "median_rmse" : median_rmse, "average_nrmse" : average_nrmse, "median_nrmse" : median_nrmse, "average_mae" : average_mae, "median_mae" : median_mae }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate RMSE and NRMSE for gene count predictions.')
    parser.add_argument('ground_truth_dir', type=str, help='Directory containing ground truth CSV files.')
    parser.add_argument('predictions_dir',  type=str, help='Directory containing predicted CSV files.')
    parser.add_argument('--normalize', type=str, choices=['range', 'mean'], default='range', 
                        help='Normalization method for NRMSE calculation. Choose "range" or "mean".')

    args = parser.parse_args()
    
    calculate_rmse(args.ground_truth_dir, args.predictions_dir, args.normalize)
