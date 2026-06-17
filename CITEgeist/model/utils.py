"""
Utility functions for CITEgeist including neighbor detection and validation.
"""

import gc
import logging
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error


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
    df.to_csv(filepath, index=False)


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
        filename=log_file, filemode="w", format="%(asctime)s - %(levelname)s - %(message)s", level=logging.INFO
    )
    logging.info("Logging initialized.")


# Spatial Neighbor Functions


def compute_optimal_radius(
    adata: sc.AnnData,
    n_rings: float = 3.0,
    fallback_spacing: float = 100.0,
    min_spots: int = 10,
) -> float:
    """
    Compute optimal spatial radius for neighborhood-based operations.

    Uses median nearest-neighbor distance from spatial coordinates
    to determine spot spacing, then computes radius for n_rings rings.

    The default of 3 rings was empirically validated via radius sweep
    benchmarking (see docs/plans/2026-02-07-radius-sweep-benchmark-design.md).

    Args:
        adata: AnnData with obsm['spatial'] coordinates.
        n_rings: Number of neighborhood rings (default 3.0, empirically optimal).
        fallback_spacing: Fallback spot spacing if auto-detection fails (default 100µm).
        min_spots: Minimum spots required for reliable estimation.

    Returns:
        Optimal radius in same units as spatial coordinates.

    Example:
        >>> radius = compute_optimal_radius(adata)
        >>> model.run_cell_proportion_model(radius=radius)
    """
    logger = logging.getLogger(__name__)

    coords = adata.obsm.get("spatial", None)

    if coords is None:
        logger.warning(
            "No spatial coordinates in adata.obsm['spatial']. " "Using fallback spacing: %s", fallback_spacing
        )
        spot_spacing = fallback_spacing
    elif coords.shape[0] < min_spots:
        logger.warning("Only %s spots (< %s). Using fallback spacing: %s", coords.shape[0], min_spots, fallback_spacing)
        spot_spacing = fallback_spacing
    else:
        from scipy.spatial import cKDTree  # pylint: disable=import-outside-toplevel

        # Filter out non-finite coordinates
        valid_mask = np.all(np.isfinite(coords), axis=1)
        valid_coords = coords[valid_mask]

        if valid_coords.shape[0] < min_spots:
            logger.warning(
                "Only %s valid coordinates. Using fallback spacing: %s", valid_coords.shape[0], fallback_spacing
            )
            spot_spacing = fallback_spacing
        else:
            tree = cKDTree(valid_coords)
            # Query 2 nearest neighbors (self + nearest)
            distances, _ = tree.query(valid_coords, k=2)
            # distances[:, 1] is distance to nearest neighbor ([:, 0] is self = 0)
            nn_distances = distances[:, 1]
            spot_spacing = float(np.median(nn_distances))

            logger.info(
                "Detected spot spacing: %s units (median of %s nearest-neighbor distances)",
                round(spot_spacing, 2),
                len(nn_distances),
            )

    # Compute radius for n rings
    # Formula: radius = spot_spacing * (n_rings + 0.05)
    # The 0.05 buffer ensures neighbors at exactly n_rings distance are included
    radius = spot_spacing * (n_rings + 0.05)

    logger.info(
        "Computed optimal radius: %s (%s rings, spot_spacing=%s)", round(radius, 2), n_rings, round(spot_spacing, 2)
    )

    return radius


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
    coordinates = adata.obsm["spatial"]
    central_coord = coordinates[spot_index]

    # Identify all spots within the given radius, excluding the central spot
    neighbors = [
        idx
        for idx, coord in enumerate(coordinates)
        if idx != spot_index and np.linalg.norm(coord - central_coord) <= radius
    ]

    # Convert indices to spot names
    neighbor_names = adata.obs_names[neighbors].tolist()
    central_spot_name = adata.obs_names[spot_index]

    return central_spot_name, neighbor_names, spot_index, neighbors


def get_neighbors_with_fixed_radius(
    spot_index, adata, radius=50, include_center=True
):  # In Visium, 2 rings gives you the adjacent 6 units
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
    _, _, spot_index, neighbors = find_fixed_radius_neighbors(spot_index, adata, radius)

    # Optionally include the central spot itself
    if include_center:
        neighbors = [spot_index] + neighbors

    logging.debug("Total neighbors for spot %s within radius %s: %s", spot_index, radius, neighbors)
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
    import random  # pylint: disable=import-outside-toplevel

    # Select random spots
    random_spots = random.sample(range(adata.shape[0]), min(num_spots, adata.shape[0]))

    # Define colorblind-friendly colors that contrast well with orange background
    # Using dark blue, white, and black for maximum contrast
    color_dict = {
        "Other spots": "#00FF00",  # Green
        "Neighbor": "#40E0D0",  # Turquoise
        "Central spot": "#0000FF",  # Dark blue
    }

    for spot_index in random_spots:
        # Find neighbors within the given radius
        central_spot_names, neighbor_spots_names, spot_index, _ = find_fixed_radius_neighbors(spot_index, adata, radius)

        # Create a temporary column to highlight spots
        adata.obs["highlight"] = "Other spots"
        adata.obs.loc[neighbor_spots_names, "highlight"] = "Neighbor"
        adata.obs.loc[central_spot_names, "highlight"] = "Central spot"

        # Plot using `sc.pl.spatial` with custom colors
        sc.pl.spatial(
            adata,
            color="highlight",
            title=f"Neighbors within {radius} units for Spot {central_spot_names}",
            spot_size=75,
            frameon=False,
            palette=color_dict,
        )

        # Clean up temporary column after each plot
        adata.obs.drop(columns=["highlight"], inplace=True)


def assert_neighborhood_size(adata, cell_profile_dict, radius=50, num_spots=5):
    """ """
    import random  # pylint: disable=import-outside-toplevel

    # Select random spots
    random_spots = random.sample(range(adata.shape[0]), min(num_spots, adata.shape[0]))

    neighborhood_sizes = []

    for spot_index in random_spots:
        # Find neighbors within the given radius
        central_spot_names, neighbor_spots_names, spot_index, _ = find_fixed_radius_neighbors(spot_index, adata, radius)

    central_spot_names = list(central_spot_names) if not isinstance(central_spot_names, list) else central_spot_names
    neighbor_spots_names = (
        list(neighbor_spots_names) if not isinstance(neighbor_spots_names, list) else neighbor_spots_names
    )

    assert all(
        x <= len(cell_profile_dict) for x in neighborhood_sizes
    ), f"Some neighborhood values in the list are less than {len(cell_profile_dict)} celltypes being deconvoluted"


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

    return {"JSD": quants_jsd[1], "RMSE": RMSE, "Sum_RMSE": all_rmse, "MAE": MAE, "Sum_MAE": all_mae, "corr": corr}


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
        dense_data = layer_data.toarray() if hasattr(layer_data, "toarray") else layer_data

        # Create DataFrame
        df = pd.DataFrame(dense_data, index=adata.obs.index, columns=adata.var.index)

        # Extract cell type name from layer name consistently
        cell_type = layer_name.split("_genes_pass")[0]

        # Save with standardized naming including pass number
        if pass_number is not None:
            output_file = os.path.join(target_dir, f"{cell_type}_layer_pass{pass_number}.csv")
        else:
            output_file = os.path.join(target_dir, f"{cell_type}_layer.csv")

        df.to_csv(output_file)
        logging.info("Exported layer '%s' to '%s'", layer_name, output_file)


def calculate_expression_metrics(ground_truth_dir, predictions_dir, normalize="range", pass_number=None):
    """
    Calculate performance metrics for gene expression predictions.

    Handles mismatched profiles between ground truth and predictions:
    - Profiles in predictions but not ground truth (e.g., Nonspecific, Unknown) are
      compared against zero (any allocation is wrong)
    - Profiles in ground truth but not predictions are compared against zero predictions
      (missed cell types)

    Args:
        ground_truth_dir (str): Directory containing ground truth CSV files
        predictions_dir (str): Directory containing prediction CSV files
        normalize (str): Normalization method for NRMSE ('range' or 'mean')
        pass_number (int, optional): If specified, look for predictions in pass-specific subdirectory

    Returns:
        dict: Dictionary containing performance metrics per cell type and overall statistics.
              Includes special keys '_spurious_profiles' and '_missed_profiles' for tracking.
    """
    if normalize not in ("range", "mean"):
        raise ValueError("Normalization type must be 'range' or 'mean'")

    metrics_per_cell_type = {}

    # Adjust predictions directory if pass number specified
    if pass_number is not None:
        predictions_dir = os.path.join(predictions_dir, f"pass{pass_number}")

    logging.info("Ground truth directory: %s", ground_truth_dir)
    logging.info("Ground truth files: %s", sorted(os.listdir(ground_truth_dir)))
    logging.info("Layer directory: %s", predictions_dir)
    logging.info("Layer files: %s", sorted(os.listdir(predictions_dir)))

    # Get sorted lists of files with pass number handling
    gt_files = sorted([f for f in os.listdir(ground_truth_dir) if f.endswith("_GT.csv")])
    if pass_number is not None:
        pred_files = sorted([f for f in os.listdir(predictions_dir) if f.endswith(f"_layer_pass{pass_number}.csv")])
    else:
        pred_files = sorted([f for f in os.listdir(predictions_dir) if f.endswith("_layer.csv")])

    print("GT files: ", gt_files)
    print("Pred files: ", pred_files)

    # Extract cell type names from filenames
    # GT: "B-cells_GT.csv" -> "B-cells"
    # Pred: "B-cells_Protein_1_B-cells_Protein_2_layer_pass1.csv" -> need to match to "B-cells"
    gt_cell_types = {f.replace("_GT.csv", ""): f for f in gt_files}

    # For predictions, we need to extract the base cell type
    # Handle formats like: "B-cells_Protein_1_B-cells_Protein_2_layer_pass1.csv"
    # or "Cancer_Epithelial_Protein_1_Cancer_Epithelial_Protein_2_layer_pass1.csv"
    pred_cell_types = {}
    for pred_file in pred_files:
        # Remove the layer suffix
        if pass_number is not None:
            base = pred_file.replace(f"_layer_pass{pass_number}.csv", "")
        else:
            base = pred_file.replace("_layer.csv", "")
        pred_cell_types[base] = pred_file

    # Create a mapping from prediction names to ground truth names
    # Try to match by finding GT name at the start of the prediction name
    pred_to_gt_map = {}
    for pred_name in pred_cell_types:
        matched_gt = None
        # Try exact match first (for simple names like "Unknown")
        if pred_name in gt_cell_types:
            matched_gt = pred_name
        else:
            # Try to find GT cell type that matches the start of pred name
            # Sort by length descending to match longest first (e.g., "Cancer Epithelial" before "Cancer")
            for gt_name in sorted(gt_cell_types.keys(), key=len, reverse=True):
                # Normalize names for comparison (replace spaces with underscores)
                gt_normalized = gt_name.replace(" ", "_")
                pred_normalized = pred_name.replace(" ", "_")
                if pred_normalized.startswith(gt_normalized + "_") or pred_normalized == gt_normalized:
                    matched_gt = gt_name
                    break
        pred_to_gt_map[pred_name] = matched_gt

    # Track matched, spurious (in pred but not gt), and missed (in gt but not pred) profiles
    matched_pred_names = set()
    spurious_profiles = []

    for pred_name, gt_name in pred_to_gt_map.items():
        if gt_name is None:
            spurious_profiles.append(pred_name)
        else:
            matched_pred_names.add(pred_name)

    # Find GT profiles that weren't matched by any prediction
    matched_gt_names = set(gt for gt in pred_to_gt_map.values() if gt is not None)
    missed_profiles = [gt for gt in gt_cell_types if gt not in matched_gt_names]

    logging.info("Matched profiles: %s", len(matched_pred_names))
    logging.info("Spurious profiles (pred only): %s", spurious_profiles)
    logging.info("Missed profiles (gt only): %s", missed_profiles)

    # Helper function to calculate metrics for a pair of dataframes
    def _calculate_pair_metrics(gt_df, pred_df, cell_type, normalize):
        # Find common genes and spots
        common_genes = gt_df.index.intersection(pred_df.index)
        common_spots = gt_df.columns.intersection(pred_df.columns)

        if len(common_genes) == 0 or len(common_spots) == 0:
            logging.warning("No common genes or spots for %s. Skipping.", cell_type)
            return None

        # Subset and normalize data
        gt_subset = gt_df.reindex(index=common_genes, columns=common_spots)
        pred_subset = pred_df.reindex(index=common_genes, columns=common_spots)

        gt_log = pd.DataFrame(np.log1p(gt_subset.values), index=common_genes, columns=common_spots)
        pred_log = pd.DataFrame(np.log1p(pred_subset.values), index=common_genes, columns=common_spots)

        # Calculate metrics
        mse = mean_squared_error(gt_log.values, pred_log.values)
        rmse = np.sqrt(mse)
        mae = mean_absolute_error(gt_log.values, pred_log.values)

        # Calculate NRMSE
        if normalize == "range":
            range_gt = gt_log.values.max() - gt_log.values.min()
            nrmse = rmse / range_gt if range_gt != 0 else np.nan
        elif normalize == "mean":
            mean_gt = gt_log.values.mean()
            nrmse = rmse / mean_gt if mean_gt != 0 else np.nan
        else:
            raise ValueError("Normalization type must be 'range' or 'mean'")

        return {"RMSE": rmse, "NRMSE": nrmse, "MAE": mae}

    # Process matched profiles
    for pred_name, gt_name in pred_to_gt_map.items():
        if gt_name is None:
            continue  # Handle spurious profiles separately

        gt_filepath = os.path.join(ground_truth_dir, gt_cell_types[gt_name])
        pred_filepath = os.path.join(predictions_dir, pred_cell_types[pred_name])

        if not os.path.exists(pred_filepath):
            logging.warning("Prediction file for %s not found. Skipping.", gt_name)
            continue

        # Load data
        gt_df = pd.read_csv(gt_filepath, index_col=0)
        pred_df = pd.read_csv(pred_filepath, index_col=0)

        metrics = _calculate_pair_metrics(gt_df, pred_df, gt_name, normalize)
        if metrics is not None:
            metrics_per_cell_type[gt_name] = metrics
            logging.info(
                "Metrics for %s: RMSE=%s, NRMSE=%s, MAE=%s", gt_name, metrics["RMSE"], metrics["NRMSE"], metrics["MAE"]
            )

    # Process spurious profiles (predictions without ground truth)
    # These should be compared against zero - any allocation is wrong
    for pred_name in spurious_profiles:
        pred_filepath = os.path.join(predictions_dir, pred_cell_types[pred_name])

        if not os.path.exists(pred_filepath):
            continue

        pred_df = pd.read_csv(pred_filepath, index_col=0)

        # Create zero ground truth with same shape
        gt_df = pd.DataFrame(0.0, index=pred_df.index, columns=pred_df.columns)

        metrics = _calculate_pair_metrics(gt_df, pred_df, pred_name, normalize)
        if metrics is not None:
            # Mark as spurious in the key
            metrics_per_cell_type[f"[SPURIOUS] {pred_name}"] = metrics
            logging.warning(
                "Spurious profile %s: RMSE=%s, NRMSE=%s, MAE=%s (compared against zero)",
                pred_name,
                metrics["RMSE"],
                metrics["NRMSE"],
                metrics["MAE"],
            )

    # Process missed profiles (ground truth without predictions)
    # These are compared against zero predictions - missed cell types
    for gt_name in missed_profiles:
        gt_filepath = os.path.join(ground_truth_dir, gt_cell_types[gt_name])

        if not os.path.exists(gt_filepath):
            continue

        gt_df = pd.read_csv(gt_filepath, index_col=0)

        # Create zero predictions with same shape
        pred_df = pd.DataFrame(0.0, index=gt_df.index, columns=gt_df.columns)

        metrics = _calculate_pair_metrics(gt_df, pred_df, gt_name, normalize)
        if metrics is not None:
            # Mark as missed in the key
            metrics_per_cell_type[f"[MISSED] {gt_name}"] = metrics
            logging.warning(
                "Missed profile %s: RMSE=%s, NRMSE=%s, MAE=%s (predicted as zero)",
                gt_name,
                metrics["RMSE"],
                metrics["NRMSE"],
                metrics["MAE"],
            )

    # Store tracking info as special keys
    metrics_per_cell_type["_spurious_profiles"] = spurious_profiles
    metrics_per_cell_type["_missed_profiles"] = missed_profiles

    return metrics_per_cell_type


# Spatial Diagnostic Plotting Functions
