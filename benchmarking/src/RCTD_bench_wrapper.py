import os
import pandas as pd
import argparse
from benchmarking_spot_deconv import main as benchmark

def parse_args():
    # Setup argument parser to accept input paths
    parser = argparse.ArgumentParser(description="Benchmark RCTD predictions against ground truth.")
    parser.add_argument('--rctd_results_dir', type=str, required=True, help="Directory containing RCTD result files.")
    parser.add_argument('--ground_truth_dir', type=str, required=True, help="Directory containing ground truth files.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save the output benchmarking metrics.")
    
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_args()

    # Directories and file patterns from arguments
    rctd_results_dir = args.rctd_results_dir
    ground_truth_dir = args.ground_truth_dir
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

    # File patterns for predictions and ground truth
    prediction_pattern = "Wu_rep_*_RCTD_predictions.csv"
    ground_truth_pattern = "Wu_ST_*_prop.csv"

    # Gather files
    rctd_prediction_files = sorted(
        [os.path.join(rctd_results_dir, f) for f in os.listdir(rctd_results_dir) if f.endswith("RCTD_deconv_predictions.csv")]
    )
    ground_truth_files = sorted(
        [os.path.join(ground_truth_dir, f) for f in os.listdir(ground_truth_dir) if f.endswith("_prop.csv")]
    )

    # Sanity check
    if len(rctd_prediction_files) != len(ground_truth_files):
        raise ValueError("Mismatch between RCTD prediction files and ground-truth files!")

    # Run benchmarking for each pair of files
    for pred_file, gt_file in zip(rctd_prediction_files, ground_truth_files):
        replicate_name = os.path.basename(pred_file).split("_")[2]  # Extract replicate number
        output_file = os.path.join(output_dir, f"RCTD_benchmarking_metrics_rep_{replicate_name}.csv")
        print(f"Processing: {pred_file} vs {gt_file}, saving to {output_file}")

        try:
            # Run benchmarking
            results = benchmark(test_spots_file=gt_file, spot_composition_file=pred_file)

            # Convert results to DataFrame for saving
            results_df = pd.DataFrame([results])  # Wrap in a list to create a single-row DataFrame
            results_df.to_csv(output_file, index=False)
        except Exception as e:
            print(f"Error processing {pred_file} vs {gt_file}: {e}")

if __name__ == "__main__":
    main()
