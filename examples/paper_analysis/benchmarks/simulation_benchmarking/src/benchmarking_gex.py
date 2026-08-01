import argparse
import os

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error


def calculate_rmse(ground_truth_dir, predictions_dir, normalize="range"):
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

            gt_df = pd.read_csv(gt_filepath, index_col=0, on_bad_lines="skip")
            pred_df = pd.read_csv(pred_filepath, index_col=0, on_bad_lines="skip")

            # Normalize gene names to uppercase for case-insensitive matching
            gt_df.columns = gt_df.columns.str.upper()
            pred_df.columns = pred_df.columns.str.upper()

            # index = spot names, columns = gene names
            common_spots = gt_df.index.intersection(pred_df.index)
            gt_df, pred_df = gt_df.loc[common_spots], pred_df.loc[common_spots]

            common_genes = gt_df.columns.intersection(pred_df.columns)
            gt_df, pred_df = gt_df[common_genes], pred_df[common_genes]

            if gt_df.empty or pred_df.empty:
                print(f"No common genes or spots for {cell_type}. Skipping.")
                continue

            # Filter to GT-expressed genes (nonzero total before log1p)
            gt_totals = gt_df.sum(axis=0)
            expressed_mask = gt_totals > 0
            n_total = len(expressed_mask)
            n_expressed = int(expressed_mask.sum())
            gt_df = gt_df.loc[:, expressed_mask]
            pred_df = pred_df.loc[:, expressed_mask]
            print(
                f"  {cell_type}: {n_expressed}/{n_total} genes expressed in GT (filtered {n_total - n_expressed} zero-sum genes)"
            )

            if gt_df.empty:
                print(f"No expressed genes for {cell_type}. Skipping.")
                continue

            gt_df, pred_df = np.log1p(gt_df), np.log1p(pred_df)

            mse = mean_squared_error(gt_df.values, pred_df.values)
            rmse = np.sqrt(mse)
            mae = mean_absolute_error(gt_df.values, pred_df.values)

            if normalize == "range":
                range_gt = gt_df.values.max() - gt_df.values.min()
                nrmse = rmse / range_gt if range_gt != 0 else np.nan
            elif normalize == "mean":
                mean_gt = gt_df.values.mean()
                nrmse = rmse / mean_gt if mean_gt != 0 else np.nan
            else:
                raise ValueError("Normalization type not recognized. Use 'range' or 'mean'.")

            # Flattened Pearson r (all spots x genes concatenated)
            gt_flat = gt_df.values.ravel()
            pred_flat = pred_df.values.ravel()
            pearson_r = np.corrcoef(gt_flat, pred_flat)[0, 1]

            # Spot-layer r: per-gene Pearson r across spots, then averaged
            gene_rs = []
            for gene in gt_df.columns:
                gt_col = gt_df[gene].values
                pred_col = pred_df[gene].values
                if gt_col.std() > 0 and pred_col.std() > 0:
                    gene_rs.append(np.corrcoef(gt_col, pred_col)[0, 1])
            spot_layer_r = float(np.mean(gene_rs)) if gene_rs else np.nan

            metrics_per_cell_type[cell_type] = {
                "RMSE": rmse,
                "NRMSE": nrmse,
                "MAE": mae,
                "pearson_r": pearson_r,
                "spot_layer_r": spot_layer_r,
            }

    metrics_df = pd.DataFrame.from_dict(metrics_per_cell_type, orient="index")
    metrics_df.to_csv(os.path.join(predictions_dir, "per_celltype_metrics.csv"))

    all_rmse_values = [metrics["RMSE"] for metrics in metrics_per_cell_type.values()]
    all_nrmse_values = [metrics["NRMSE"] for metrics in metrics_per_cell_type.values()]
    all_mae_values = [metrics["MAE"] for metrics in metrics_per_cell_type.values()]
    all_pearson_r_values = [metrics["pearson_r"] for metrics in metrics_per_cell_type.values()]
    all_spot_layer_values = [metrics["spot_layer_r"] for metrics in metrics_per_cell_type.values()]

    average_rmse = np.mean(all_rmse_values)
    median_rmse = np.median(all_rmse_values)

    average_nrmse = np.nanmean(all_nrmse_values)
    median_nrmse = np.nanmedian(all_nrmse_values)

    average_mae = np.nanmean(all_mae_values)
    median_mae = np.nanmedian(all_mae_values)

    mean_pearson_r = float(np.nanmean(all_pearson_r_values))
    median_pearson_r = float(np.nanmedian(all_pearson_r_values))

    mean_spot_layer_r = float(np.nanmean(all_spot_layer_values))
    median_spot_layer_r = float(np.nanmedian(all_spot_layer_values))

    print("RMSE, NRMSE, MAE, Pearson r, and spot-layer r per cell type:")
    for cell_type, metrics in metrics_per_cell_type.items():
        print(
            f"\t{cell_type}: RMSE: {metrics['RMSE']:.4f}, NRMSE: {metrics['NRMSE']:.4f}, "
            f"MAE: {metrics['MAE']:.4f}, pearson_r: {metrics['pearson_r']:.4f}, "
            f"spot_layer_r: {metrics['spot_layer_r']:.4f}"
        )

    print("\nOverall RMSE statistics:")
    print(f"\tAverage RMSE: {average_rmse:.4f}")
    print(f"\tMedian RMSE:  {median_rmse:.4f}")

    print("\nOverall NRMSE statistics:")
    print(f"\tAverage NRMSE: {average_nrmse:.4f}")
    print(f"\tMedian NRMSE:  {median_nrmse:.4f}")

    print("\nOverall MAE statistics:")
    print(f"\tAverage MAE: {average_mae:.4f}")
    print(f"\tMedian MAE:  {median_mae:.4f}")

    print("\nOverall Pearson r statistics (flattened):")
    print(f"\tMean Pearson r:   {mean_pearson_r:.4f}")
    print(f"\tMedian Pearson r: {median_pearson_r:.4f}")

    print("\nOverall spot-layer r statistics (per-gene r, averaged):")
    print(f"\tMean spot-layer r:   {mean_spot_layer_r:.4f}")
    print(f"\tMedian spot-layer r: {median_spot_layer_r:.4f}")

    print("Per-cell type metrics saved to per_celltype_metrics.csv")
    return {
        "metrics_per_cell_type": metrics_per_cell_type,
        "average_rmse": average_rmse,
        "median_rmse": median_rmse,
        "average_nrmse": average_nrmse,
        "median_nrmse": median_nrmse,
        "average_mae": average_mae,
        "median_mae": median_mae,
        "mean_pearson_r": mean_pearson_r,
        "median_pearson_r": median_pearson_r,
        "mean_spot_layer_r": mean_spot_layer_r,
        "median_spot_layer_r": median_spot_layer_r,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RMSE and NRMSE for gene count predictions.")
    parser.add_argument("ground_truth_dir", type=str, help="Directory containing ground truth CSV files.")
    parser.add_argument("predictions_dir", type=str, help="Directory containing predicted CSV files.")
    parser.add_argument(
        "--normalize",
        type=str,
        choices=["range", "mean"],
        default="range",
        help="Normalization method for NRMSE calculation.",
    )

    args = parser.parse_args()
    calculate_rmse(args.ground_truth_dir, args.predictions_dir, args.normalize)
