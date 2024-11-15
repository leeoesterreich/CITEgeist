import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

def benchmark_performance(test_spots_metadata_mtrx, spot_composition_mtrx, column_names):
    # Check if inputs are numpy arrays (similar to matrix in R)
    if not isinstance(test_spots_metadata_mtrx, np.ndarray):
        raise ValueError("ERROR: test_spots_metadata_mtrx must be a numpy array!")

    if not isinstance(spot_composition_mtrx, np.ndarray):
        raise ValueError("ERROR: spot_composition_mtrx must be a numpy array!")

    # Initialize JSD matrix
    true_jsd_mtrx = np.zeros((test_spots_metadata_mtrx.shape[0], 1))

    # Calculate Jensen-Shannon Divergence (JSD)
    for i in range(test_spots_metadata_mtrx.shape[0]):
        x = np.vstack([test_spots_metadata_mtrx[i, :], spot_composition_mtrx[i, :]])

        if np.sum(spot_composition_mtrx[i, :]) > 0:
            true_jsd_mtrx[i, 0] = jensenshannon(x[0], x[1], base=2)  # JSD calculation
        else:
            true_jsd_mtrx[i, 0] = 1  # If composition is zero, assign 1 to JSD

    # Calculate RMSE and MAE for each cell type
    RMSE = {}
    MAE = {}
    all_rmse = 0
    all_mae = 0

    for i in range(test_spots_metadata_mtrx.shape[1]):
        # Sum of squared differences for the current cell type
        mse = np.sum((test_spots_metadata_mtrx[:, i] - spot_composition_mtrx[:, i]) ** 2)
        all_rmse += mse  # Accumulate the MSE for the overall RMSE calculation
        RMSE[column_names[i]] = np.sqrt(mse / test_spots_metadata_mtrx.shape[0])  # RMSE for each cell type

        mae = mean_absolute_error(test_spots_metadata_mtrx[:, i], spot_composition_mtrx[:, i])
        all_mae += mae
        MAE[column_names[i]] = mae

    # Calculate overall RMSE and MAE
    all_rmse = np.sqrt(all_rmse / (test_spots_metadata_mtrx.shape[0] * test_spots_metadata_mtrx.shape[1]))
    all_mae = all_mae / test_spots_metadata_mtrx.shape[1]

    # Quantiles of the JSD
    quants_jsd = np.quantile(np.min(true_jsd_mtrx, axis=1), [0.25, 0.5, 0.75])

    # Calculate correlation between matrices
    flat_test = test_spots_metadata_mtrx.flatten()
    flat_spot = spot_composition_mtrx.flatten()
    corr, _ = pearsonr(flat_test, flat_spot)

    return {
        'JSD': quants_jsd[1],  # Median JSD
        'RMSE': RMSE,
        'Sum_RMSE': all_rmse,
        'MAE': MAE,
        'Sum_MAE': all_mae,
        'corr': corr
    }

def main(test_spots_file, spot_composition_file):
    # Load data from CSV files & sort row names (spot IDs)
    test_spots_df = pd.read_csv(test_spots_file, index_col=0).sort_index()
    spot_composition_df = pd.read_csv(spot_composition_file, index_col=0).sort_index()

    # Drop "spot_x" and "spot_y" columns from both dataframes (if present)
    test_spots_df = test_spots_df.drop(columns=['spot_x', 'spot_y'], errors='ignore')
    spot_composition_df = spot_composition_df.drop(columns=['spot_x', 'spot_y'], errors='ignore')

    # Sort column names (cell types)
    test_spots_df = test_spots_df.sort_index(axis=1)
    spot_composition_df = spot_composition_df.sort_index(axis=1)

    # Check if columns match
    if not np.array_equal(test_spots_df.columns, spot_composition_df.columns):
        raise ValueError("ERROR: The column names in the input CSV files do not match or are not in the same order!")

    # Check if rows match
    if not np.array_equal(test_spots_df.index, spot_composition_df.index):
        raise ValueError("ERROR: The row names in the input CSV files do not match or are not in the same order!")

    # Convert DataFrames to numpy arrays
    test_spots_metadata_mtrx = test_spots_df.values
    spot_composition_mtrx = spot_composition_df.values
    column_names = test_spots_df.columns.tolist()  # Get column names for RMSE and MAE

    # Benchmark performance
    results = benchmark_performance(test_spots_metadata_mtrx, spot_composition_mtrx, column_names)

    # Print results
    print("Benchmarking Results:")
    print(f"Median JSD: {results['JSD']}")
    print("RMSE by Cell Type:")
    for cell_type, rmse_value in results['RMSE'].items():
        print(f"  {cell_type}: {rmse_value:.4f}")
    print(f"Sum RMSE: {results['Sum_RMSE']:.4f}")
    print("MAE by Cell Type:")
    for cell_type, mae_value in results['MAE'].items():
        print(f"  {cell_type}: {mae_value:.4f}")
    print(f"Sum MAE: {results['Sum_MAE']:.4f}")
    print(f"Correlation: {results['corr']:.4f}")

    return(results)
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python benchmarking.py <test_spots_file> <spot_composition_file>")
        sys.exit(1)

    test_spots_file = sys.argv[1]
    spot_composition_file = sys.argv[2]
    main(test_spots_file, spot_composition_file)
