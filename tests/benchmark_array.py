#!/usr/bin/env python
"""
Benchmark Array Results Analyzer

Aggregates and compares results from CITEgeist array benchmarking jobs.
Compares manual vs auto profile modes across different hyperparameter settings.

Usage:
    python tests/benchmark_array.py --results-dir test_results/
    python tests/benchmark_array.py --results-dir test_results/ --output-dir benchmark_summary/
"""

import argparse
import glob
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_output_folder_name(folder_name: str) -> Optional[Dict[str, float]]:
    """
    Parse hyperparameters from output folder name.

    Expected format: radius_4_lambda_1_alpha_0.7_max_y_change_0.2.

    Returns:
        Dictionary with parsed parameters or None if parsing fails.
    """
    try:
        # Remove trailing dot and suffix
        clean_name = folder_name.rstrip('.')
        if 'FilteredRadii' in clean_name:
            clean_name = clean_name.split('.')[0]

        parts = clean_name.split('_')
        params = {}

        i = 0
        while i < len(parts):
            if parts[i] == 'radius' and i + 1 < len(parts):
                params['radius'] = float(parts[i + 1])
                i += 2
            elif parts[i] == 'lambda' and i + 1 < len(parts):
                params['lambda_reg'] = float(parts[i + 1])
                i += 2
            elif parts[i] == 'alpha' and i + 1 < len(parts):
                params['alpha_elastic'] = float(parts[i + 1])
                i += 2
            elif parts[i] == 'max' and i + 3 < len(parts) and parts[i + 1] == 'y' and parts[i + 2] == 'change':
                params['max_y_change'] = float(parts[i + 3])
                i += 4
            else:
                i += 1

        return params if params else None
    except (ValueError, IndexError) as e:
        logger.debug(f"Could not parse folder name '{folder_name}': {e}")
        return None


def find_result_files(results_dir: str) -> Dict[str, List[Dict]]:
    """
    Find all result CSV files in the results directory.

    Returns:
        Dictionary mapping test_set/profile_mode to list of result file info.
    """
    results = {
        'cell_proportions': [],
        'gex_metrics': [],
        'profile_discovery': [],
    }

    results_path = Path(results_dir)

    # Search pattern: test_results/{test_set}/{profile_mode}_profiles/...
    for test_set in ['mixed', 'high_seg']:
        for profile_mode in ['manual', 'auto']:
            mode_dir = results_path / test_set / f"{profile_mode}_profiles"

            if not mode_dir.exists():
                logger.debug(f"Directory not found: {mode_dir}")
                continue

            # Find all output subdirectories
            for output_dir in mode_dir.glob("radius_*"):
                if not output_dir.is_dir():
                    continue

                params = parse_output_folder_name(output_dir.name)
                if params is None:
                    continue

                base_info = {
                    'test_set': test_set,
                    'profile_mode': profile_mode,
                    'output_dir': str(output_dir),
                    **params
                }

                # Find cell proportion results
                for prop_file in output_dir.glob("*_cellprop_results_summary_*.csv"):
                    result_type = 'global' if '_global_' in prop_file.name else 'finetune'
                    results['cell_proportions'].append({
                        **base_info,
                        'file_path': str(prop_file),
                        'result_type': result_type,
                    })

                # Find GEX metrics
                for gex_file in output_dir.glob("*_gex_metrics_pass*.csv"):
                    pass_num = 1 if 'pass1' in gex_file.name else 2
                    results['gex_metrics'].append({
                        **base_info,
                        'file_path': str(gex_file),
                        'pass_number': pass_num,
                    })

                # Find profile discovery metrics (auto mode only)
                for disc_file in output_dir.glob("*_profile_discovery_metrics_*.csv"):
                    results['profile_discovery'].append({
                        **base_info,
                        'file_path': str(disc_file),
                    })

    return results


def load_cell_proportion_results(file_info_list: List[Dict]) -> pd.DataFrame:
    """Load and aggregate cell proportion benchmark results."""
    rows = []

    for info in file_info_list:
        try:
            df = pd.read_csv(info['file_path'])

            # Extract metrics from the single-row DataFrame
            row = {
                'test_set': info['test_set'],
                'profile_mode': info['profile_mode'],
                'result_type': info['result_type'],
                'radius': info.get('radius'),
                'lambda_reg': info.get('lambda_reg'),
                'alpha_elastic': info.get('alpha_elastic'),
                'max_y_change': info.get('max_y_change'),
            }

            # Add metrics from CSV
            for col in df.columns:
                if col in ['JSD', 'RMSE', 'MAE', 'corr', 'pearson', 'spearman']:
                    row[col] = df[col].values[0]

            rows.append(row)

        except Exception as e:
            logger.warning(f"Could not load {info['file_path']}: {e}")

    return pd.DataFrame(rows)


def load_gex_metrics(file_info_list: List[Dict]) -> pd.DataFrame:
    """Load and aggregate gene expression metrics."""
    rows = []

    for info in file_info_list:
        try:
            df = pd.read_csv(info['file_path'])

            row = {
                'test_set': info['test_set'],
                'profile_mode': info['profile_mode'],
                'pass_number': info['pass_number'],
                'radius': info.get('radius'),
                'lambda_reg': info.get('lambda_reg'),
                'alpha_elastic': info.get('alpha_elastic'),
                'max_y_change': info.get('max_y_change'),
            }

            # Parse metrics from the long-form DataFrame
            # Expected columns: Pass, Metric, Value
            for _, metric_row in df.iterrows():
                metric_name = metric_row.get('Metric', '')
                value = metric_row.get('Value', np.nan)

                # Clean up metric name for column
                clean_name = metric_name.replace(' ', '_')
                row[clean_name] = value

            rows.append(row)

        except Exception as e:
            logger.warning(f"Could not load {info['file_path']}: {e}")

    return pd.DataFrame(rows)


def load_profile_discovery_metrics(file_info_list: List[Dict]) -> pd.DataFrame:
    """Load and aggregate profile discovery metrics (auto mode only)."""
    rows = []

    for info in file_info_list:
        try:
            df = pd.read_csv(info['file_path'])

            row = {
                'test_set': info['test_set'],
                'profile_mode': info['profile_mode'],
                'radius': info.get('radius'),
                'lambda_reg': info.get('lambda_reg'),
                'alpha_elastic': info.get('alpha_elastic'),
                'max_y_change': info.get('max_y_change'),
            }

            # Add metrics from CSV
            for col in df.columns:
                if col not in ['matched_pairs', 'unmatched_discovered', 'missing_ground_truth']:
                    row[col] = df[col].values[0]

            rows.append(row)

        except Exception as e:
            logger.warning(f"Could not load {info['file_path']}: {e}")

    return pd.DataFrame(rows)


def compute_summary_statistics(df: pd.DataFrame, group_cols: List[str],
                                metric_cols: List[str]) -> pd.DataFrame:
    """Compute mean and std for metrics grouped by specified columns."""
    if df.empty:
        return pd.DataFrame()

    # Filter to existing columns
    existing_group = [c for c in group_cols if c in df.columns]
    existing_metrics = [c for c in metric_cols if c in df.columns]

    if not existing_group or not existing_metrics:
        return pd.DataFrame()

    agg_funcs = {col: ['mean', 'std', 'count'] for col in existing_metrics}
    summary = df.groupby(existing_group).agg(agg_funcs)

    # Flatten column names
    summary.columns = ['_'.join(col).strip() for col in summary.columns.values]

    return summary.reset_index()


def compare_manual_vs_auto(prop_df: pd.DataFrame, gex_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a comparison table of manual vs auto profile modes.

    Returns DataFrame with columns for each mode and the difference.
    """
    comparisons = []

    # Cell proportion comparison
    if not prop_df.empty:
        for test_set in prop_df['test_set'].unique():
            subset = prop_df[(prop_df['test_set'] == test_set) &
                            (prop_df['result_type'] == 'global')]

            manual = subset[subset['profile_mode'] == 'manual']
            auto = subset[subset['profile_mode'] == 'auto']

            if not manual.empty and not auto.empty:
                for metric in ['JSD', 'RMSE', 'MAE', 'corr']:
                    if metric in manual.columns:
                        manual_val = manual[metric].mean()
                        auto_val = auto[metric].mean()

                        # For JSD/RMSE/MAE, lower is better; for corr, higher is better
                        if metric == 'corr':
                            diff_pct = ((auto_val - manual_val) / abs(manual_val)) * 100 if manual_val != 0 else 0
                            better = 'auto' if auto_val > manual_val else 'manual'
                        else:
                            diff_pct = ((manual_val - auto_val) / abs(manual_val)) * 100 if manual_val != 0 else 0
                            better = 'auto' if auto_val < manual_val else 'manual'

                        comparisons.append({
                            'test_set': test_set,
                            'metric_type': 'cell_proportion',
                            'metric': metric,
                            'manual': manual_val,
                            'auto': auto_val,
                            'diff_pct': diff_pct,
                            'better': better,
                        })

    # GEX comparison
    if not gex_df.empty:
        for test_set in gex_df['test_set'].unique():
            for pass_num in [1, 2]:
                subset = gex_df[(gex_df['test_set'] == test_set) &
                               (gex_df['pass_number'] == pass_num)]

                manual = subset[subset['profile_mode'] == 'manual']
                auto = subset[subset['profile_mode'] == 'auto']

                if not manual.empty and not auto.empty:
                    for metric in ['Holistic_RMSE', 'Matched_RMSE', 'Holistic_MAE']:
                        if metric in manual.columns:
                            manual_val = manual[metric].mean()
                            auto_val = auto[metric].mean()

                            # Lower is better for all these metrics
                            diff_pct = ((manual_val - auto_val) / abs(manual_val)) * 100 if manual_val != 0 else 0
                            better = 'auto' if auto_val < manual_val else 'manual'

                            comparisons.append({
                                'test_set': test_set,
                                'metric_type': f'gex_pass{pass_num}',
                                'metric': metric,
                                'manual': manual_val,
                                'auto': auto_val,
                                'diff_pct': diff_pct,
                                'better': better,
                            })

    return pd.DataFrame(comparisons)


def print_summary_report(prop_df: pd.DataFrame, gex_df: pd.DataFrame,
                         discovery_df: pd.DataFrame, comparison_df: pd.DataFrame):
    """Print a formatted summary report to stdout."""

    print("\n" + "=" * 80)
    print("CITEGEIST BENCHMARK ARRAY RESULTS SUMMARY")
    print("=" * 80)

    # Cell Proportion Results
    print("\n" + "-" * 40)
    print("CELL PROPORTION BENCHMARKS")
    print("-" * 40)

    if not prop_df.empty:
        for test_set in sorted(prop_df['test_set'].unique()):
            print(f"\n[{test_set.upper()}]")
            subset = prop_df[(prop_df['test_set'] == test_set) &
                            (prop_df['result_type'] == 'global')]

            for mode in ['manual', 'auto']:
                mode_data = subset[subset['profile_mode'] == mode]
                if not mode_data.empty:
                    print(f"\n  {mode.upper()} profiles (n={len(mode_data)}):")
                    for metric in ['JSD', 'RMSE', 'MAE', 'corr']:
                        if metric in mode_data.columns:
                            mean_val = mode_data[metric].mean()
                            std_val = mode_data[metric].std()
                            print(f"    {metric:8s}: {mean_val:.4f} ± {std_val:.4f}")
    else:
        print("  No cell proportion results found.")

    # Gene Expression Results
    print("\n" + "-" * 40)
    print("GENE EXPRESSION BENCHMARKS")
    print("-" * 40)

    if not gex_df.empty:
        for test_set in sorted(gex_df['test_set'].unique()):
            print(f"\n[{test_set.upper()}]")

            for pass_num in [1]:  # Focus on Pass 1 for now
                print(f"\n  Pass {pass_num}:")
                subset = gex_df[(gex_df['test_set'] == test_set) &
                               (gex_df['pass_number'] == pass_num)]

                for mode in ['manual', 'auto']:
                    mode_data = subset[subset['profile_mode'] == mode]
                    if not mode_data.empty:
                        print(f"\n    {mode.upper()} profiles (n={len(mode_data)}):")

                        for metric in ['Holistic_RMSE', 'Matched_RMSE', 'Holistic_MAE',
                                      'N_Matched', 'N_Spurious', 'N_Missed']:
                            if metric in mode_data.columns:
                                mean_val = mode_data[metric].mean()
                                std_val = mode_data[metric].std()
                                if metric.startswith('N_'):
                                    print(f"      {metric:15s}: {mean_val:.1f} ± {std_val:.1f}")
                                else:
                                    print(f"      {metric:15s}: {mean_val:.4f} ± {std_val:.4f}")
    else:
        print("  No gene expression results found.")

    # Profile Discovery Results (auto mode only)
    print("\n" + "-" * 40)
    print("PROFILE DISCOVERY METRICS (AUTO MODE)")
    print("-" * 40)

    if not discovery_df.empty:
        for test_set in sorted(discovery_df['test_set'].unique()):
            print(f"\n[{test_set.upper()}] (n={len(discovery_df[discovery_df['test_set'] == test_set])})")
            subset = discovery_df[discovery_df['test_set'] == test_set]

            for metric in ['profile_recovery_rate', 'marker_precision', 'marker_recall',
                          'false_discovery_rate', 'n_matched', 'n_discovered']:
                if metric in subset.columns:
                    mean_val = subset[metric].mean()
                    std_val = subset[metric].std()
                    if metric.startswith('n_'):
                        print(f"    {metric:25s}: {mean_val:.1f} ± {std_val:.1f}")
                    else:
                        print(f"    {metric:25s}: {mean_val:.2%} ± {std_val:.2%}")
    else:
        print("  No profile discovery results found.")

    # Manual vs Auto Comparison
    print("\n" + "-" * 40)
    print("MANUAL vs AUTO COMPARISON")
    print("-" * 40)

    if not comparison_df.empty:
        for test_set in sorted(comparison_df['test_set'].unique()):
            print(f"\n[{test_set.upper()}]")
            subset = comparison_df[comparison_df['test_set'] == test_set]

            for metric_type in sorted(subset['metric_type'].unique()):
                print(f"\n  {metric_type}:")
                type_subset = subset[subset['metric_type'] == metric_type]

                for _, row in type_subset.iterrows():
                    arrow = "↓" if row['better'] == 'auto' else "↑"
                    sign = "+" if row['diff_pct'] > 0 else ""
                    print(f"    {row['metric']:15s}: manual={row['manual']:.4f}, auto={row['auto']:.4f} "
                          f"({sign}{row['diff_pct']:.1f}% {arrow} {row['better']})")
    else:
        print("  Insufficient data for comparison.")

    print("\n" + "=" * 80)
    print("END OF REPORT")
    print("=" * 80 + "\n")


def save_results(prop_df: pd.DataFrame, gex_df: pd.DataFrame,
                 discovery_df: pd.DataFrame, comparison_df: pd.DataFrame,
                 output_dir: str):
    """Save aggregated results to CSV files."""

    os.makedirs(output_dir, exist_ok=True)

    if not prop_df.empty:
        prop_path = os.path.join(output_dir, 'cell_proportion_results.csv')
        prop_df.to_csv(prop_path, index=False)
        logger.info(f"Saved cell proportion results to {prop_path}")

    if not gex_df.empty:
        gex_path = os.path.join(output_dir, 'gex_metrics_results.csv')
        gex_df.to_csv(gex_path, index=False)
        logger.info(f"Saved GEX metrics to {gex_path}")

    if not discovery_df.empty:
        disc_path = os.path.join(output_dir, 'profile_discovery_results.csv')
        discovery_df.to_csv(disc_path, index=False)
        logger.info(f"Saved profile discovery results to {disc_path}")

    if not comparison_df.empty:
        comp_path = os.path.join(output_dir, 'manual_vs_auto_comparison.csv')
        comparison_df.to_csv(comp_path, index=False)
        logger.info(f"Saved comparison results to {comp_path}")

    # Save summary statistics
    if not prop_df.empty:
        prop_summary = compute_summary_statistics(
            prop_df[prop_df['result_type'] == 'global'],
            ['test_set', 'profile_mode'],
            ['JSD', 'RMSE', 'MAE', 'corr']
        )
        if not prop_summary.empty:
            summary_path = os.path.join(output_dir, 'cell_proportion_summary.csv')
            prop_summary.to_csv(summary_path, index=False)
            logger.info(f"Saved cell proportion summary to {summary_path}")

    if not gex_df.empty:
        gex_summary = compute_summary_statistics(
            gex_df,
            ['test_set', 'profile_mode', 'pass_number'],
            ['Holistic_RMSE', 'Matched_RMSE', 'Holistic_MAE', 'N_Matched', 'N_Spurious', 'N_Missed']
        )
        if not gex_summary.empty:
            summary_path = os.path.join(output_dir, 'gex_metrics_summary.csv')
            gex_summary.to_csv(summary_path, index=False)
            logger.info(f"Saved GEX metrics summary to {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze CITEgeist benchmark array results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Analyze results in default location
    python tests/benchmark_array.py --results-dir test_results/

    # Save summary to specific output directory
    python tests/benchmark_array.py --results-dir test_results/ --output-dir benchmark_summary/

    # Quiet mode (only save files, minimal output)
    python tests/benchmark_array.py --results-dir test_results/ --quiet
        """
    )

    parser.add_argument('--results-dir', type=str, default='test_results/',
                        help='Directory containing benchmark results (default: test_results/)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Directory to save aggregated results (default: {results-dir}/summary/)')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress detailed output, only save files')

    args = parser.parse_args()

    # Set output directory
    if args.output_dir is None:
        args.output_dir = os.path.join(args.results_dir, 'summary')

    # Find all result files
    logger.info(f"Searching for results in: {args.results_dir}")
    result_files = find_result_files(args.results_dir)

    n_prop = len(result_files['cell_proportions'])
    n_gex = len(result_files['gex_metrics'])
    n_disc = len(result_files['profile_discovery'])

    logger.info(f"Found {n_prop} cell proportion files, {n_gex} GEX files, {n_disc} profile discovery files")

    if n_prop == 0 and n_gex == 0:
        logger.warning("No result files found. Check that jobs have completed.")
        print("\nNo results found. Make sure the benchmark jobs have completed.")
        print(f"Expected results in: {args.results_dir}/{{mixed,high_seg}}/{{manual,auto}}_profiles/")
        return 1

    # Load results
    prop_df = load_cell_proportion_results(result_files['cell_proportions'])
    gex_df = load_gex_metrics(result_files['gex_metrics'])
    discovery_df = load_profile_discovery_metrics(result_files['profile_discovery'])

    # Compare manual vs auto
    comparison_df = compare_manual_vs_auto(prop_df, gex_df)

    # Print summary report
    if not args.quiet:
        print_summary_report(prop_df, gex_df, discovery_df, comparison_df)

    # Save results
    save_results(prop_df, gex_df, discovery_df, comparison_df, args.output_dir)

    logger.info(f"Results saved to: {args.output_dir}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
