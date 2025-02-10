import os
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error
import argparse

def calculate_rmse(ground_truth_dir, predictions_dir, normalize='range'):
    metrics_per_cell_type = {}

    for gt_filename in os.listdir(ground_truth_dir):
        if gt_filename.endswith("_GT.csv"):
            cell_type = gt_filename.replace("_GT.csv", "")
            gt_filepath = os.path.join(ground_truth_dir, gt_filename)
            pred_filepath = os.path.join(predictions_dir, f"{cell_type}_layer.csv")
            
            if not os.path.exists(pred_filepath):
                cell_type_underscore = cell_type.replace(" ", "_")
                alt_pred_filepath = os.path.join(predictions_dir, f"{cell_type_underscore}_layer.csv")
                if not os.path.exists(alt_pred_filepath):
                    print(f"Prediction file for {cell_type} not found. Skipping.")
                    continue
                pred_filepath = alt_pred_filepath

            gt_df = pd.read_csv(gt_filepath, index_col=0)
            pred_df = pd.read_csv(pred_filepath, index_col=0)
            
            common_genes = gt_df.index.intersection(pred_df.index)
            gt_df, pred_df = gt_df.loc[common_genes], pred_df.loc[common_genes]
            
            common_spots = gt_df.columns.intersection(pred_df.columns)
            gt_df, pred_df = gt_df[common_spots], pred_df[common_spots]
            
            if gt_df.empty or pred_df.empty:
                print(f"No common genes or spots for {cell_type}. Skipping.")
                continue
            
            gt_df, pred_df = np.log1p(gt_df), np.log1p(pred_df)
            
            mse = mean_squared_error(gt_df.values, pred_df.values)
            rmse = np.sqrt(mse)
            mae = mean_absolute_error(gt_df.values, pred_df.values)
            
            if normalize == 'range':
                range_gt = gt_df.values.max() - gt_df.values.min()
                nrmse = rmse / range_gt if range_gt != 0 else np.nan
            elif normalize == 'mean':
                mean_gt = gt_df.values.mean()
                nrmse = rmse / mean_gt if mean_gt != 0 else np.nan
            else:
                raise ValueError("Normalization type not recognized. Use 'range' or 'mean'.")
            
            metrics_per_cell_type[cell_type] = {'RMSE': rmse, 'NRMSE': nrmse, 'MAE': mae}
    
    metrics_df = pd.DataFrame.from_dict(metrics_per_cell_type, orient='index')
    metrics_df.to_csv(os.path.join(predictions_dir, "per_celltype_metrics.csv"))
    
    all_rmse_values  = [metrics['RMSE'] for metrics in metrics_per_cell_type.values()]
    all_nrmse_values = [metrics['NRMSE'] for metrics in metrics_per_cell_type.values()]
    all_mae_values   = [metrics['MAE'] for metrics in metrics_per_cell_type.values()]
    average_rmse = np.mean(all_rmse_values)
    median_rmse  = np.median(all_rmse_values)

    average_nrmse = np.nanmean(all_nrmse_values)
    median_nrmse  = np.nanmedian(all_nrmse_values)

    average_mae = np.nanmean(all_mae_values)
    median_mae  = np.nanmedian(all_mae_values)

    print("RMSE, NRMSE, and MAE per cell type:")
    for cell_type, metrics in metrics_per_cell_type.items():
        print(f"\t{cell_type}: RMSE: {metrics['RMSE']:.4f}, NRMSE: {metrics['NRMSE']:.4f}, MAE: {metrics['MAE']:.4f}")

    print("\nOverall RMSE statistics:")
    print(f"\tAverage RMSE: {average_rmse:.4f}")
    print(f"\tMedian RMSE:  {median_rmse:.4f}")

    print("\nOverall NRMSE statistics:")
    print(f"\tAverage NRMSE: {average_nrmse:.4f}")
    print(f"\tMedian NRMSE:  {median_nrmse:.4f}")

    print("\nOverall MAE statistics:")
    print(f"\tAverage MAE: {average_mae:.4f}")
    print(f"\tMedian MAE:  {median_mae:.4f}")
    
    print("Per-cell type metrics saved to per_celltype_metrics.csv")
    return { "metrics_per_cell_type" : metrics_per_cell_type, "average_rmse" : average_rmse, "median_rmse" : median_rmse, "average_nrmse" : average_nrmse, "median_nrmse" : median_nrmse, "average_mae" : average_mae, "median_mae" : median_mae }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate RMSE and NRMSE for gene count predictions.')
    parser.add_argument('ground_truth_dir', type=str, help='Directory containing ground truth CSV files.')
    parser.add_argument('predictions_dir', type=str, help='Directory containing predicted CSV files.')
    parser.add_argument('--normalize', type=str, choices=['range', 'mean'], default='range', help='Normalization method for NRMSE calculation.')
    
    args = parser.parse_args()
    calculate_rmse(args.ground_truth_dir, args.predictions_dir, args.normalize)
