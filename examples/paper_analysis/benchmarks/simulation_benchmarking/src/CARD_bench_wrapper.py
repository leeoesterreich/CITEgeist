"""
Benchmark CARD predictions against ground truth for simulation data.

Supports both reference and reference-free mode outputs.

Usage:
    python CARD_bench_wrapper.py \
        --card_results_dir /path/to/CARD/high_seg \
        --ground_truth_dir /path/to/ST_sim \
        --output_dir /path/to/output \
        --mode reference  # or reference_free
"""

import os
import argparse

import pandas as pd
from benchmarking_spot_deconv import main as benchmark


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark CARD predictions against ground truth.")
    parser.add_argument(
        "--card_results_dir",
        type=str,
        required=True,
        help="Directory containing CARD result files.",
    )
    parser.add_argument(
        "--ground_truth_dir",
        type=str,
        required=True,
        help="Directory containing ground truth files.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the output benchmarking metrics.",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="reference",
        choices=["reference", "reference_free"],
        help="CARD mode: 'reference' or 'reference_free'",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    card_results_dir = args.card_results_dir
    ground_truth_dir = args.ground_truth_dir
    output_dir = args.output_dir
    mode = args.mode

    os.makedirs(output_dir, exist_ok=True)

    # Determine prediction file suffix based on mode
    if mode == "reference":
        suffix = "_CARD_deconv_predictions.csv"
        method_name = "CARD"
    else:
        suffix = "_CARD_reffree_deconv_predictions.csv"
        method_name = "CARD_reffree"

    # Gather prediction files
    card_prediction_files = sorted(
        [os.path.join(card_results_dir, f) for f in os.listdir(card_results_dir) if f.endswith(suffix)]
    )

    # Gather ground truth files
    ground_truth_files = sorted(
        [os.path.join(ground_truth_dir, f) for f in os.listdir(ground_truth_dir) if f.endswith("_prop.csv")]
    )

    # Sanity check
    if len(card_prediction_files) != len(ground_truth_files):
        print(
            f"Warning: Found {len(card_prediction_files)} prediction files "
            f"and {len(ground_truth_files)} ground truth files"
        )

    # Match files by replicate number
    for pred_file in card_prediction_files:
        # Extract replicate number (e.g., "Wu_ST_0" -> "0")
        basename = os.path.basename(pred_file)
        parts = basename.replace(suffix, "").split("_")
        rep_num = parts[-1]  # Last part should be the replicate number

        # Find matching ground truth file
        gt_file = None
        for gt in ground_truth_files:
            if f"Wu_ST_{rep_num}_prop.csv" in gt:
                gt_file = gt
                break

        if gt_file is None:
            print(f"Warning: No ground truth found for {basename}, skipping")
            continue

        output_file = os.path.join(output_dir, f"{method_name}_benchmarking_metrics_rep_{rep_num}.csv")
        print(f"Processing: {pred_file} vs {gt_file}")

        try:
            results = benchmark(test_spots_file=gt_file, spot_composition_file=pred_file)
            results_df = pd.DataFrame([results])
            results_df.to_csv(output_file, index=False)
            print(f"  Saved to: {output_file}")
        except Exception as e:
            print(f"  Error: {e}")


if __name__ == "__main__":
    main()
