import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def plot_parameter_heatmaps(base_folder, output_dir='figures', label=None):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize dictionaries to store metrics for each parameter combination
    metrics_data = {
        'pass1': {'rmse': {}, 'nrmse': {}, 'mae': {}},
        'pass2': {'rmse': {}, 'nrmse': {}, 'mae': {}},
        'improvements': {'rmse': {}, 'nrmse': {}, 'mae': {}}
    }
    
    # Walk through the directory structure
    for root, dirs, files in os.walk(base_folder):
        for file in files:
            if file.endswith('_metrics_combined.csv') and 'Wu_rep_' in file:
                # Initialize parameters
                radius = None
                lambda_prior = None
                
                # Extract parameters from path
                path_parts = root.split('/')
                for part in path_parts:
                    if 'radius_' in part and 'lambda_prior_weight_' in part:
                        try:
                            part = part.rstrip('.')
                            radius_part = part.split('_lambda_prior')[0]
                            radius = float(radius_part.split('radius_')[1])
                            
                            lambda_part = part.split('lambda_prior_weight_')[1]
                            lambda_prior = float(lambda_part)
                            break
                        except (IndexError, ValueError) as e:
                            print(f"Error parsing combined parameters from {part}: {e}")
                            continue
                
                if radius is None or lambda_prior is None:
                    print(f"Warning: Could not extract parameters from path: {root}")
                    continue
                
                param_key = (radius, lambda_prior)
                
                try:
                    # Read combined metrics file
                    df = pd.read_csv(os.path.join(root, file))
                    
                    # Initialize lists for this parameter combination if not exists
                    for pass_type in metrics_data:
                        for metric in metrics_data[pass_type]:
                            if param_key not in metrics_data[pass_type][metric]:
                                metrics_data[pass_type][metric][param_key] = []
                    
                    # Extract pass1 metrics
                    pass1_metrics = df[df['Pass'] == 'Pass 1']
                    pass1_rmse = pass1_metrics[pass1_metrics['Metric'] == 'Average RMSE']['Value'].values[0]
                    pass1_nrmse = pass1_metrics[pass1_metrics['Metric'] == 'Average NRMSE']['Value'].values[0]
                    pass1_mae = pass1_metrics[pass1_metrics['Metric'] == 'Average MAE']['Value'].values[0]
                    
                    # Extract pass2 metrics
                    pass2_metrics = df[df['Pass'] == 'Pass 2']
                    pass2_rmse = pass2_metrics[pass2_metrics['Metric'] == 'Average RMSE']['Value'].values[0]
                    pass2_nrmse = pass2_metrics[pass2_metrics['Metric'] == 'Average NRMSE']['Value'].values[0]
                    pass2_mae = pass2_metrics[pass2_metrics['Metric'] == 'Average MAE']['Value'].values[0]
                    
                    # Store pass1 metrics
                    metrics_data['pass1']['rmse'][param_key].append(pass1_rmse)
                    metrics_data['pass1']['nrmse'][param_key].append(pass1_nrmse)
                    metrics_data['pass1']['mae'][param_key].append(pass1_mae)
                    
                    # Store pass2 metrics
                    metrics_data['pass2']['rmse'][param_key].append(pass2_rmse)
                    metrics_data['pass2']['nrmse'][param_key].append(pass2_nrmse)
                    metrics_data['pass2']['mae'][param_key].append(pass2_mae)
                    
                    # Calculate and store improvements (as percentages)
                    metrics_data['improvements']['rmse'][param_key].append(
                        ((pass1_rmse - pass2_rmse) / pass1_rmse) * 100)
                    metrics_data['improvements']['nrmse'][param_key].append(
                        ((pass1_nrmse - pass2_nrmse) / pass1_nrmse) * 100)
                    metrics_data['improvements']['mae'][param_key].append(
                        ((pass1_mae - pass2_mae) / pass1_mae) * 100)
                    
                except Exception as e:
                    print(f"Error processing file for parameters ({radius}, {lambda_prior}): {e}")
                    continue
    
    # Verify we have data before plotting
    if not any(metrics_data[pass_type][metric] for pass_type in metrics_data for metric in metrics_data[pass_type]):
        raise ValueError("No valid data was collected. Check the file structure and paths.")
    
    # Create figure with subplots
    fig, axes = plt.subplots(3, 3, figsize=(20, 20))
    fig.suptitle('Parameter Sweep Results', fontsize=16)
    
    # Function to create annotation text
    def format_annotation(mean, std):
        return f'{mean:.3f}\n(±{std:.3f})'
    
    # Plot heatmaps
    metrics = ['rmse', 'nrmse', 'mae']
    titles = ['RMSE', 'NRMSE', 'MAE']
    passes = ['pass1', 'pass2', 'improvements']
    
    for i, (metric, title) in enumerate(zip(metrics, titles)):
        for j, pass_type in enumerate(passes):
            # Create DataFrames for means and stds
            means = {}
            stds = {}
            for param_key in metrics_data[pass_type][metric]:
                means[param_key] = np.mean(metrics_data[pass_type][metric][param_key])
                stds[param_key] = np.std(metrics_data[pass_type][metric][param_key])
            
            # Convert to DataFrames
            radii = sorted(list(set(k[0] for k in means.keys())))
            lambdas = sorted(list(set(k[1] for k in means.keys())))
            
            if not radii or not lambdas:
                print(f"No data for {pass_type} {metric}")
                continue
                
            mean_matrix = np.zeros((len(lambdas), len(radii)))
            annotation_texts = np.empty((len(lambdas), len(radii)), dtype=object)
            
            for idx_l, l in enumerate(lambdas):
                for idx_r, r in enumerate(radii):
                    if (r, l) in means:
                        mean_matrix[idx_l, idx_r] = means[(r, l)]
                        annotation_texts[idx_l, idx_r] = format_annotation(
                            means[(r, l)],
                            stds[(r, l)]
                        )
                    else:
                        mean_matrix[idx_l, idx_r] = np.nan
                        annotation_texts[idx_l, idx_r] = 'N/A'
            
            # Create heatmap
            cmap = 'RdYlGn' if pass_type == 'improvements' else 'viridis'
            center = 0 if pass_type == 'improvements' else None
            
            sns.heatmap(mean_matrix, 
                       ax=axes[j,i], 
                       cmap=cmap,
                       center=center,
                       xticklabels=radii,
                       yticklabels=lambdas,
                       annot=annotation_texts,
                       fmt='',
                       annot_kws={'size': 8})
            
            row_title = 'Pass 1' if pass_type == 'pass1' else 'Pass 2' if pass_type == 'pass2' else 'Improvement'
            axes[j,i].set_title(f'{row_title} {title}')
            axes[j,i].set_xlabel('Radius')
            axes[j,i].set_ylabel('Lambda Prior Weight')
    
    plt.tight_layout()
    
    # Construct filename with optional label
    filename = 'parameter_sweep_heatmaps'
    if label:
        filename += f'_{label}'
    filename += '.png'
    
    output_path = os.path.join(output_dir, filename)
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved heatmaps to: {output_path}")

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate parameter sweep heatmaps')
    parser.add_argument('base_folder', type=str, help='Base folder containing the benchmark results')
    parser.add_argument('--output-dir', type=str, default='figures',
                        help='Directory to save the output figures (default: figures)')
    parser.add_argument('--label', type=str, help='Optional label to append to output filename')
    
    args = parser.parse_args()
    
    plot_parameter_heatmaps(args.base_folder, args.output_dir, args.label)