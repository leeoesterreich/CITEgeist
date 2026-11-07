import os
import gc
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error, mean_absolute_error

def validate_cell_profile_dict(cell_profile_dict):
    """
    Validate the structure of a cell profile dictionary.
    """
    if not isinstance(cell_profile_dict, dict):
        return False
    return all(isinstance(k, str) and isinstance(v, dict) for k, v in cell_profile_dict.items())

def save_results_to_output(results, filepath):
    """
    Save results as a CSV file.
    """
    df = pd.DataFrame(results)
    df.to_csv(filepath)

def cleanup_memory():
    """
    Force garbage collection to free memory.
    """
    gc.collect()

def setup_logging(output_folder, sample_name):
    """
    Set up dynamic logging.
    """
    log_file = os.path.join(output_folder, f"{sample_name}_CITEgeist.log")
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )
    logging.info('Logging initialized.')

### 📏 **Spatial Neighbor Functions**



def find_fixed_radius_neighbors(spot_index, adata, radius=50):
    """
    Find neighbors within a fixed radius of a given spot.

    Args:
        spot_index (int): Index of the central spot in the AnnData object.
        adata (AnnData): Spatial transcriptomics dataset with obsm['spatial'] containing spot coordinates.
        radius (float): Fixed radius to identify neighboring spots.

    Returns:
        tuple: (central_spot_name, list of neighbor spot names)
    """
    coordinates = adata.obsm['spatial']
    central_coord = coordinates[spot_index]
    
    # Identify all spots within the given radius, excluding the central spot
    neighbors = [idx for idx, coord in enumerate(coordinates)
                 if idx != spot_index and np.linalg.norm(coord - central_coord) <= radius]
    
    # Convert indices to spot names
    neighbor_names = adata.obs_names[neighbors].tolist()
    central_spot_name = adata.obs_names[spot_index]
    
    return central_spot_name, neighbor_names, spot_index, neighbors


def get_neighbors_with_fixed_radius(spot_index, adata, radius=50, include_center=True): # In Visium, 2 rings gives you the adjacent 6 units
    """
    Get indices of neighboring spots based on a fixed radius around the central spot.

    Parameters:
    - spot_index (int): The index of the central spot.
    - adata (AnnData): The AnnData object with spatial coordinates in obsm['spatial'].
    - radius (float): Fixed radius for finding neighbors.
    - include_center (bool): Whether to include the central spot in the neighbor list.

    Returns:
    - List of indices representing the central spot and its neighbors.
    """
    # Find neighbors within the given radius
    central_spot_names, neighbor_spots_names, spot_index, neighbors = find_fixed_radius_neighbors(spot_index, adata, radius)
    
    # Optionally include the central spot itself
    if include_center:
        neighbors = [spot_index] + neighbors
    
    logging.debug(f"Total neighbors for spot {spot_index} within radius {radius}: {neighbors}")
    return neighbors
def plot_neighbors_with_fixed_radius(adata, radius=50, num_spots=5):
    """
    Plot neighbors for multiple random central spots using `sc.pl.spatial`.

    Args:
        adata (AnnData): Spatial transcriptomics dataset with obsm['spatial'].
        radius (float): Fixed radius to identify neighboring spots.
        num_spots (int): Number of random spots to visualize.

    Returns:
        None: Displays a series of spatial plots showing neighbors.
    """
    import random
    
    # Select random spots
    random_spots = random.sample(range(adata.shape[0]), min(num_spots, adata.shape[0]))
    
    # Define colorblind-friendly colors that contrast well with orange background
    # Using dark blue, white, and black for maximum contrast
    color_dict = {
        'Other spots': '#00FF00',  # Green
        'Neighbor': '#40E0D0',     # Turquoise
        'Central spot': '#0000FF'  # Dark blue
    }
    
    for spot_index in random_spots:
        # Find neighbors within the given radius
        central_spot_names, neighbor_spots_names, spot_index, neighbors = find_fixed_radius_neighbors(spot_index, adata, radius)

        # Create a temporary column to highlight spots
        adata.obs['highlight'] = 'Other spots'
        adata.obs.loc[neighbor_spots_names, 'highlight'] = 'Neighbor'
        adata.obs.loc[central_spot_names, 'highlight'] = 'Central spot'
        
        # Plot using `sc.pl.spatial` with custom colors
        sc.pl.spatial(
            adata,
            color='highlight',
            title=f"Neighbors within {radius} units for Spot {central_spot_names}",
            spot_size = 75,
            frameon=False,
            palette=color_dict
        )
        
        # Clean up temporary column after each plot
        adata.obs.drop(columns=['highlight'], inplace=True)
        
def assert_neighborhood_size(adata, cell_profile_dict, radius=50, num_spots=5):
    """

    """
    import random
    
    # Select random spots
    random_spots = random.sample(range(adata.shape[0]), min(num_spots, adata.shape[0]))
    
    neighborhood_sizes = []
    
    for spot_index in random_spots:
        # Find neighbors within the given radius
        central_spot_names, neighbor_spots_names, spot_index, neighbors = find_fixed_radius_neighbors(spot_index, adata, radius)

    central_spot_names = list(central_spot_names) if not isinstance(central_spot_names, list) else central_spot_names
    neighbor_spots_names = list(neighbor_spots_names) if not isinstance(neighbor_spots_names, list) else neighbor_spots_names

    neighborhood_size = len(central_spot_names + neighbor_spots_names)
    
    
    assert all(x <= len(cell_profile_dict) for x in neighborhood_sizes), f"Some neighborhood values in the list are less than {len(cell_profile_dict)} celltypes being deconvoluted"

def benchmark_cell_proportions(true_proportions, predicted_proportions, cell_type_names):
    """
    Calculate performance metrics for cell type proportion predictions.
    
    Args:
        true_proportions (np.ndarray): Ground truth cell type proportions matrix
        predicted_proportions (np.ndarray): Predicted cell type proportions matrix
        cell_type_names (list): Names of cell types corresponding to matrix columns
        
    Returns:
        dict: Dictionary containing various performance metrics
    """
    if not isinstance(true_proportions, np.ndarray) or not isinstance(predicted_proportions, np.ndarray):
        raise ValueError("Input proportions must be numpy arrays")

    # Initialize JSD matrix
    true_jsd_mtrx = np.zeros((true_proportions.shape[0], 1))

    # Calculate Jensen-Shannon Divergence
    for i in range(true_proportions.shape[0]):
        x = np.vstack([true_proportions[i, :], predicted_proportions[i, :]])
        if np.sum(predicted_proportions[i, :]) > 0:
            true_jsd_mtrx[i, 0] = jensenshannon(x[0], x[1], base=2)
        else:
            true_jsd_mtrx[i, 0] = 1

    # Calculate per-celltype and overall metrics
    RMSE = {}
    MAE = {}
    all_rmse = 0
    all_mae = 0

    for i in range(true_proportions.shape[1]):
        mse = np.sum((true_proportions[:, i] - predicted_proportions[:, i]) ** 2)
        all_rmse += mse
        RMSE[cell_type_names[i]] = np.sqrt(mse / true_proportions.shape[0])

        mae = mean_absolute_error(true_proportions[:, i], predicted_proportions[:, i])
        all_mae += mae
        MAE[cell_type_names[i]] = mae

    # Calculate overall metrics
    all_rmse = np.sqrt(all_rmse / (true_proportions.shape[0] * true_proportions.shape[1]))
    all_mae = all_mae / true_proportions.shape[1]

    # Calculate JSD quantiles and correlation
    quants_jsd = np.quantile(np.min(true_jsd_mtrx, axis=1), [0.25, 0.5, 0.75])
    corr, _ = pearsonr(true_proportions.flatten(), predicted_proportions.flatten())

    return {
        'JSD': quants_jsd[1],
        'RMSE': RMSE,
        'Sum_RMSE': all_rmse,
        'MAE': MAE,
        'Sum_MAE': all_mae,
        'corr': corr
    }

def export_anndata_layers(adata, output_dir, pass_number=None):
    """
    Export all layers of an AnnData object to separate CSV files.
    Creates separate folders for different passes.
    
    Args:
        adata (AnnData): AnnData object containing the layers to export
        output_dir (str): Base directory where CSV files will be saved
        pass_number (int, optional): If specified, only export layers from this pass
    """
    # Create base output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create pass-specific directory if needed
    if pass_number is not None:
        target_dir = os.path.join(output_dir, f"pass{pass_number}")
        os.makedirs(target_dir, exist_ok=True)
    else:
        target_dir = output_dir

    # Filter layers for specific pass if requested
    layer_pattern = f"_pass{pass_number}" if pass_number is not None else None
    
    for layer_name in adata.layers.keys():
        # Skip if not matching pass number
        if layer_pattern is not None and layer_pattern not in layer_name:
            continue
        
        # Extract data and ensure it's dense
        layer_data = adata.layers[layer_name]
        dense_data = layer_data.toarray() if hasattr(layer_data, 'toarray') else layer_data
        
        # Create DataFrame
        df = pd.DataFrame(dense_data, index=adata.obs.index, columns=adata.var.index)
        
        # Extract cell type name from layer name consistently
        cell_type = layer_name.split('_genes_pass')[0]
        
        # Save with standardized naming including pass number
        if pass_number is not None:
            output_file = os.path.join(target_dir, f"{cell_type}_layer_pass{pass_number}.csv")
        else:
            output_file = os.path.join(target_dir, f"{cell_type}_layer.csv")
            
        df.to_csv(output_file)
        logging.info(f"Exported layer '{layer_name}' to '{output_file}'")

def calculate_expression_metrics(ground_truth_dir, predictions_dir, normalize='range', pass_number=None):
    """
    Calculate performance metrics for gene expression predictions.
    
    Args:
        ground_truth_dir (str): Directory containing ground truth CSV files
        predictions_dir (str): Directory containing prediction CSV files
        normalize (str): Normalization method for NRMSE ('range' or 'mean')
        pass_number (int, optional): If specified, look for predictions in pass-specific subdirectory
        
    Returns:
        dict: Dictionary containing performance metrics per cell type and overall statistics
    """
    metrics_per_cell_type = {}
    
    # Adjust predictions directory if pass number specified
    if pass_number is not None:
        predictions_dir = os.path.join(predictions_dir, f"pass{pass_number}")
    
    logging.info(f"Ground truth directory: {ground_truth_dir}")
    logging.info(f"Ground truth files: {sorted(os.listdir(ground_truth_dir))}")
    logging.info(f"Layer directory: {predictions_dir}")
    logging.info(f"Layer files: {sorted(os.listdir(predictions_dir))}")

    # Get sorted lists of files with pass number handling
    gt_files = sorted([f for f in os.listdir(ground_truth_dir) if f.endswith('_GT.csv')])
    if pass_number is not None:
        pred_files = sorted([f for f in os.listdir(predictions_dir) if f.endswith(f'_layer_pass{pass_number}.csv')])
    else:
        pred_files = sorted([f for f in os.listdir(predictions_dir) if f.endswith('_layer.csv')])

    print("GT files: ", gt_files)
    print("Pred files: ", pred_files)

    assert len(gt_files) == len(pred_files), "Number of ground truth files and prediction files do not match"

    # Create a dictionary to map cell types to their ground truth files
    gt_file_map = {f.replace('_GT.csv', ''): f for f in gt_files}

    # Create a list to store matched prediction and ground truth file pairs
    matched_files = []

    for pred_file in pred_files:
        # Remove pass number suffix if present
        base_pred_file = pred_file.replace(f'_pass{pass_number}', '') if pass_number is not None else pred_file
        cell_type = base_pred_file.replace('_layer.csv', '').split('_')[0]
        if cell_type in gt_file_map:
            matched_files.append((pred_file, gt_file_map[cell_type]))

    # Sort matched files by the cell type name
    matched_files.sort(key=lambda x: x[0])

    # Calculate metrics
    for pred_filename, gt_filename in matched_files:
        cell_type = gt_filename.replace("_GT.csv", "")
        gt_filepath = os.path.join(ground_truth_dir, gt_filename)
        pred_filepath = os.path.join(predictions_dir, pred_filename)

        if not os.path.exists(pred_filepath):
            logging.warning(f"Prediction file for {cell_type} not found. Skipping.")
            continue

        # Load and preprocess data
        gt_df = pd.read_csv(gt_filepath, index_col=0)
        pred_df = pd.read_csv(pred_filepath, index_col=0)

        # Find common genes and spots
        common_genes = gt_df.index.intersection(pred_df.index)
        common_spots = gt_df.columns.intersection(pred_df.columns)

        if len(common_genes) == 0 or len(common_spots) == 0:
            logging.warning(f"No common genes or spots for {cell_type}. Skipping.")
            continue

        # Subset and normalize data
        gt_subset = gt_df.reindex(index=common_genes, columns=common_spots)
        pred_subset = pred_df.reindex(index=common_genes, columns=common_spots)
        
        gt_df = pd.DataFrame(np.log1p(gt_subset.values), index=common_genes, columns=common_spots)
        pred_df = pd.DataFrame(np.log1p(pred_subset.values), index=common_genes, columns=common_spots)

        # Calculate metrics
        mse = mean_squared_error(gt_df.values, pred_df.values)
        rmse = np.sqrt(mse)
        mae = mean_absolute_error(gt_df.values, pred_df.values)

        # Calculate NRMSE
        if normalize == 'range':
            range_gt = gt_df.values.max() - gt_df.values.min()
            nrmse = rmse / range_gt if range_gt != 0 else np.nan
        elif normalize == 'mean':
            mean_gt = gt_df.values.mean()
            nrmse = rmse / mean_gt if mean_gt != 0 else np.nan
        else:
            raise ValueError("Normalization type must be 'range' or 'mean'")

        # Assertions that also print the celltype
        assert nrmse is not None and nrmse != np.nan, f"NRMSE is None for {cell_type}"
        assert rmse is not None and rmse != np.nan, f"RMSE is None for {cell_type}"
        assert mae is not None and mae != np.nan, f"MAE is None for {cell_type}"

        metrics_per_cell_type[cell_type] = {'RMSE': rmse, 'NRMSE': nrmse, 'MAE': mae}
        logging.info(f"Metrics for {cell_type}: RMSE={rmse:.4f}, NRMSE={nrmse:.4f}, MAE={mae:.4f}")
        
    return metrics_per_cell_type


