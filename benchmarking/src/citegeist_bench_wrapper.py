import os
import argparse
import pandas as pd
from benchmarking_gex import calculate_rmse as run_rmse_benchmarking_function  # Import RMSE benchmarking

def run_rmse_benchmarking(ground_truth_gex_dir, predictions_dir, output_csv, normalize="range"):
    """Run the RMSE benchmarking function and return results as a DataFrame."""
    print(f"\nRunning RMSE benchmarking:")
    print(f"Ground truth directory: {ground_truth_gex_dir}")
    print(f"Predictions directory: {predictions_dir}")
    print(f"Output CSV: {output_csv}")
    print(f"Normalize method: {normalize}")
    
    try:
        # Call the RMSE benchmarking function
        metrics = run_rmse_benchmarking_function(ground_truth_gex_dir, predictions_dir, normalize=normalize)
        
        if not metrics:
            print("Warning: No metrics returned from benchmarking function")
            return None
            
        # Print available keys in metrics for debugging
        print(f"Available metrics keys: {metrics.keys()}")
        
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

        # Print metrics for debugging
        print("\nCalculated metrics:")
        for metric, value in metrics_values.items():
            print(f"{metric}: {value}")

        # Convert the metrics dictionary to a DataFrame
        gex_results = pd.DataFrame(list(metrics_values.items()), columns=['Metric', 'Value'])
        
        # Save the DataFrame to CSV
        gex_results.to_csv(output_csv, index=False)
        print(f"GEX benchmarking results saved to: {output_csv}")
        return gex_results
        
    except Exception as e:
        print(f"Error in run_rmse_benchmarking: {str(e)}")
        return None

def main(predictions_dir, ground_truth_gex_dir, output_dir, folder_pattern, normalize="range"):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nStarting benchmarking process:")
    print(f"Predictions directory: {predictions_dir}")
    print(f"Ground truth base directory: {ground_truth_gex_dir}")
    print(f"Output directory: {output_dir}")

    # Debug: Check if directories exist
    print("\nChecking directory existence:")
    print(f"Predictions directory exists: {os.path.exists(predictions_dir)}")
    print(f"Ground truth directory exists: {os.path.exists(ground_truth_gex_dir)}")

    # Get list of sample folders from ground truth directory
    try:
        ground_truth_samples = [f for f in os.listdir(ground_truth_gex_dir) 
                              if os.path.isdir(os.path.join(ground_truth_gex_dir, f)) and f.startswith("sample_")]
        print(f"\nFound ground truth samples: {ground_truth_samples}")
    except Exception as e:
        print(f"Error listing ground truth directory: {str(e)}")
        return

    if not ground_truth_samples:
        print("No sample_X folders found in ground truth directory")
        return

    # Process each sample
    for sample_folder in ground_truth_samples:
        try:
            # Extract sample number
            sample_number = ''.join(c for c in sample_folder if c.isdigit())
            if not sample_number:
                print(f"Could not extract number from sample folder: {sample_folder}")
                continue

            print(f"\nProcessing sample {sample_number}")

            # Find corresponding prediction folder
            pred_folder = f"Wu_rep_{sample_number}_pass1"
            pred_path = os.path.join(predictions_dir, pred_folder)
            
            if not os.path.isdir(pred_path):
                print(f"Prediction folder not found: {pred_folder}")
                continue

            # Set up paths
            ground_truth_path = os.path.join(ground_truth_gex_dir, sample_folder, "layers")
            predicted_gex_dir = os.path.join(pred_path, "layers", "pass1")

            print(f"Ground truth path: {ground_truth_path}")
            print(f"Predictions path: {predicted_gex_dir}")

            if not os.path.isdir(ground_truth_path) or not os.path.isdir(predicted_gex_dir):
                print("Missing required directories, skipping...")
                continue

            # Create a temporary directory for renamed prediction files
            temp_dir = os.path.join(predicted_gex_dir, "temp_renamed")
            try:
                os.makedirs(temp_dir, exist_ok=True)
            except Exception as e:
                print(f"Error creating temporary directory: {str(e)}")
                continue

            # Copy and rename CITEgeist output files to match expected format
            renamed_files = []
            for filename in os.listdir(predicted_gex_dir):
                if filename.endswith("_layer_pass1.csv"):
                    old_path = os.path.join(predicted_gex_dir, filename)
                    # Extract cell type name and handle spaces
                    cell_type = filename.replace("_layer_pass1.csv", "")
                    new_filename = f"{cell_type}_layer.csv"  # Changed from _GT.csv to _layer.csv
                    new_path = os.path.join(temp_dir, new_filename)
                    print(f"Copying: {filename} -> {new_filename}")
                    # Copy instead of rename to preserve originals
                    with open(old_path, 'r') as src, open(new_path, 'w') as dst:
                        dst.write(src.read())
                    renamed_files.append(new_filename)

            print(f"Copied {len(renamed_files)} files to temporary directory")

            # Output file path for metrics
            output_csv = os.path.join(output_dir, f"GEX_metrics_{pred_folder}.csv")
            
            try:
                # Run RMSE benchmarking and get the results as a DataFrame
                rmse_results = run_rmse_benchmarking(ground_truth_path, temp_dir, output_csv, normalize)
                
                if rmse_results is None:
                    print(f"Warning: No results generated for {pred_folder}")
            except Exception as e:
                print(f"Error processing {pred_folder}: {str(e)}")
            finally:
                # Clean up temporary directory
                print("\nCleaning up temporary files...")
                for filename in renamed_files:
                    try:
                        os.remove(os.path.join(temp_dir, filename))
                    except Exception as e:
                        print(f"Error removing temporary file {filename}: {str(e)}")
                try:
                    os.rmdir(temp_dir)
                    print("Temporary directory removed successfully")
                except Exception as e:
                    print(f"Error removing temporary directory: {str(e)}")
        except Exception as e:
            print(f"Error processing sample {sample_folder}: {str(e)}")
            continue

    print("\nAll benchmarking completed!")
    print(f"Check output directory: {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run RMSE benchmarking for predicted layers.")
    parser.add_argument("predictions_dir", type=str, help="Path to the directory containing predicted CSV files (with Wu_rep_X folders).")
    parser.add_argument("ground_truth_gex_dir", type=str, help="Path to the directory containing ground truth GEX files (with sample_X folders).")
    parser.add_argument("output_dir", type=str, help="Path to save the benchmarking metrics.")
    parser.add_argument("--normalize", type=str, choices=['range', 'mean'], default='range',
                        help="Normalization method for NRMSE calculation (default: 'range').")

    args = parser.parse_args()
    # Note: folder_pattern is no longer needed as we're matching sample_X with Wu_rep_X directly
    main(args.predictions_dir, args.ground_truth_gex_dir, args.output_dir, "Wu_rep_", args.normalize)

