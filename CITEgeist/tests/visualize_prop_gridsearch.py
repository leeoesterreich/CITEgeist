import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple, Optional
import ast

def plot_prop_parameter_heatmaps(base_folder: str, output_dir: str = 'figures', label: Optional[str] = None) -> None:
    """
    Create heatmaps visualizing cell proportion metrics across different parameter combinations.
    
    Args:
        base_folder: Base directory containing the parameter sweep results
        output_dir: Directory to save output figures
        label: Optional label to append to output filename
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize dictionaries to store metrics for each parameter combination
    metrics_data: Dict[str, Dict[Tuple[float, float], List[float]]] = {
        'JSD': {},
        'Sum_RMSE': {},
        'Sum_MAE': {},
        'corr': {}
    }
    
    # Walk through the directory structure
    for root, dirs, files in os.walk(base_folder):
        for file in files:
            if '_summary_finetune_' in file:
                # Extract parameters from path
                path_parts = root.split('/')
                for part in path_parts:
                    if 'lambda_' in part and 'alpha_' in part:
                        try:
                            # Extract lambda and alpha values
                            lambda_val = float(part.split('lambda_')[1].split('_alpha_')[0])
                            alpha_val = float(part.split('alpha_')[1].split('_max_y')[0])
                            param_key = (lambda_val, alpha_val)
                            
                            # Read metrics file
                            df = pd.read_csv(os.path.join(root, file))
                            
                            # Initialize lists for this parameter combination if not exists
                            for metric in metrics_data:
                                if param_key not in metrics_data[metric]:
                                    metrics_data[metric][param_key] = []
                            
                            # Store metrics
                            metrics_data['JSD'][param_key].append(df['JSD'].values[0])
                            metrics_data['Sum_RMSE'][param_key].append(df['Sum_RMSE'].values[0])
                            metrics_data['Sum_MAE'][param_key].append(df['Sum_MAE'].values[0])
                            metrics_data['corr'][param_key].append(df['corr'].values[0])
                            
                        except (IndexError, ValueError) as e:
                            print(f"Error parsing parameters from {part}: {e}")
                            continue
    
    # Verify we have data before plotting
    if not any(metrics_data[metric] for metric in metrics_data):
        raise ValueError("No valid data was collected. Check the file structure and paths.")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    fig.suptitle('Cell Proportion Parameter Sweep Results', fontsize=16)
    
    # Function to create annotation text
    def format_annotation(mean: float, std: float) -> str:
        return f'{mean:.3f}\n(±{std:.3f})'
    
    # Plot heatmaps
    metrics = list(metrics_data.keys())
    titles = ['Jensen-Shannon Divergence', 'Sum RMSE', 'Sum MAE', 'Correlation']
    
    for idx, (metric, title) in enumerate(zip(metrics, titles)):
        i, j = divmod(idx, 2)
        
        # Create DataFrames for means and stds
        means = {}
        stds = {}
        for param_key in metrics_data[metric]:
            means[param_key] = np.mean(metrics_data[metric][param_key])
            stds[param_key] = np.std(metrics_data[metric][param_key])
        
        # Convert to matrices
        lambdas = sorted(list(set(k[0] for k in means.keys())))
        alphas = sorted(list(set(k[1] for k in means.keys())))
        
        mean_matrix = np.zeros((len(lambdas), len(alphas)))
        annotation_texts = np.empty((len(lambdas), len(alphas)), dtype=object)
        
        for idx_l, l in enumerate(lambdas):
            for idx_a, a in enumerate(alphas):
                if (l, a) in means:
                    mean_matrix[idx_l, idx_a] = means[(l, a)]
                    annotation_texts[idx_l, idx_a] = format_annotation(
                        means[(l, a)],
                        stds[(l, a)]
                    )
                else:
                    mean_matrix[idx_l, idx_a] = np.nan
                    annotation_texts[idx_l, idx_a] = 'N/A'
        
        # Use viridis colormap for all plots
        cmap = 'viridis' if metric != 'corr' else 'viridis_r'
        sns.heatmap(mean_matrix,
                   ax=axes[i, j],
                   cmap=cmap,
                   xticklabels=alphas,
                   yticklabels=lambdas,
                   annot=annotation_texts,
                   fmt='',
                   annot_kws={'size': 8})
        
        axes[i, j].set_title(title)
        axes[i, j].set_xlabel('Alpha')
        axes[i, j].set_ylabel('Lambda')
    
    plt.tight_layout()
    
    # Construct filename with optional label
    filename = 'prop_parameter_sweep_heatmaps'
    if label:
        filename += f'_{label}'
    filename += '.png'
    
    output_path = os.path.join(output_dir, filename)
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved heatmaps to: {output_path}")

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate cell proportion parameter sweep heatmaps')
    parser.add_argument('base_folder', type=str, help='Base folder containing the benchmark results')
    parser.add_argument('--output-dir', type=str, default='figures',
                        help='Directory to save the output figures (default: figures)')
    parser.add_argument('--label', type=str, help='Optional label to append to output filename')
    
    args = parser.parse_args()
    
    plot_prop_parameter_heatmaps(args.base_folder, args.output_dir, args.label)

    # Example: python tests/visualize_prop_gridsearch.py path/to/results/folder --output-dir figures