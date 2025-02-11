import os
import argparse
import pandas as pd
from benchmarking_gex import calculate_rmse as run_gex_benchmarking_function  # Import GEX benchmarking
from benchmarking_spot_deconv import main as run_proportions_benchmarking_function  # Import proportions benchmarking

def run_gex_benchmarking(ground_truth_gex_dir, predicted_gex_dir, output_csv, normalize="range"):
    """Run the GEX benchmarking function and return results as a DataFrame."""
    print(f"Running GEX benchmarking for predicted layers in: {predicted_gex_dir}")
    
    # Call the GEX benchmarking function
    metrics = run_gex_benchmarking_function(ground_truth_gex_dir, predicted_gex_dir, normalize=normalize)
    
    # Extracting required metrics from the function result
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
    
    # Save the DataFrame to CSV
    gex_results.to_csv(output_csv, index=False)
    print(f"GEX benchmarking results saved to: {output_csv}")
    return gex_results


def run_proportions_benchmarking(ground_truth_proportions_file, predicted_proportions_file, output_csv):
    """Run the proportions benchmarking function and return results as a DataFrame."""
    print(f"Running proportions benchmarking for: {ground_truth_proportions_file} and {predicted_proportions_file}")
    results_dict = run_proportions_benchmarking_function(predicted_proportions_file, ground_truth_proportions_file)

    # Convert the dictionary into a DataFrame
    prop_results = pd.DataFrame([results_dict])

    # Save the results to CSV
    prop_results.to_csv(output_csv, index=False)
    print(f"Proportions benchmarking results saved to: {output_csv}")
    return prop_results


def main(replicates_dir, ground_truth_gex_dir, ground_truth_proportions_dir, output_dir, folder_pattern, normalize="range"):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for replicate_folder in os.listdir(replicates_dir):
        
        # Match the folder pattern
        if not replicate_folder.startswith(folder_pattern):
            continue
        print(f"Current folder: {replicate_folder}")
        
        replicate_path = os.path.join(replicates_dir, replicate_folder)
        predicted_gex_dir = os.path.join(replicate_path, "layers")
        
        # Match the replicate folder number with the ground truth sample folder
        replicate_number = replicate_folder.split("_")[-1]  # Extract the number from cell2location_map_<i>
        ground_truth_folder = f"sample_{replicate_number}"  # Corresponding sample folder in the ground truth directory
        
        # Path for the ground truth GEX folder
        ground_truth_gex_sample_folder = os.path.join(ground_truth_gex_dir, ground_truth_folder, "layers")
        
        # Ensure the GEX folder exists before proceeding
        if not os.path.isdir(ground_truth_gex_sample_folder):
            print(f"Skipping {replicate_folder}, ground truth GEX folder missing.")
            continue
        
        # Ground truth proportions file path (in ST_sim)
        ground_truth_proportions_file = os.path.join(ground_truth_proportions_dir, f"Wu_ST_{replicate_number}_prop.csv")
        
        # Predicted proportions (deconvolution) file path
        predicted_proportions_file = os.path.join(replicate_path, "cell2loc_deconv_predictions.csv")
        
        # Check if required files are present
        if not os.path.isdir(predicted_gex_dir) or not os.path.exists(predicted_proportions_file) or not os.path.exists(ground_truth_proportions_file):
            print(f"Skipping {replicate_folder}, required data missing.")
            continue

        # Output file paths for benchmarking results
        gex_output_csv  = os.path.join(output_dir, f"GEX_metrics_{replicate_folder}.csv")
        prop_output_csv = os.path.join(output_dir, f"Prop_metrics_{replicate_folder}.csv")

        # Print the current comparison being made
        print(f"\nComparing predicted layers from replicate: {replicate_folder} with ground truth: {ground_truth_gex_sample_folder}")

        # Run GEX benchmarking and save results
        run_gex_benchmarking(ground_truth_gex_sample_folder, predicted_gex_dir, gex_output_csv, normalize)

        # Run proportions benchmarking and save results
        run_proportions_benchmarking(ground_truth_proportions_file, predicted_proportions_file, prop_output_csv)

    print("All benchmarking completed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run full benchmarking for GEX and cell-type proportions.")
    parser.add_argument("replicates_dir", type=str, help="Path to the directory containing replicate subfolders.")
    parser.add_argument("ground_truth_gex_dir", type=str, help="Path to the directory containing ground truth GEX files.")
    parser.add_argument("ground_truth_proportions_dir", type=str, help="Path to the directory containing ground truth proportion files.")
    parser.add_argument("output_dir", type=str, help="Path to save the benchmarking results.")
    parser.add_argument("--folder_pattern", type=str, default="cell2location_map_",
                        help="Pattern to match replicate folders (default: 'cell2location_map_').")
    parser.add_argument("--normalize", type=str, choices=['range', 'mean'], default='range',
                        help="Normalization method for NRMSE calculation (default: 'range').")

    args = parser.parse_args()
    main(args.replicates_dir, args.ground_truth_gex_dir, args.ground_truth_proportions_dir, args.output_dir, args.folder_pattern, args.normalize)
