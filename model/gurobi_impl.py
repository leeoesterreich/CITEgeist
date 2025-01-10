import numpy as np
import gurobipy as gp
from gurobipy import Model, GRB, quicksum
import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import logging
import traceback
import gc
from .utils import get_neighbors_with_fixed_radius
import os
from pathlib import Path
from scipy.stats import spearmanr
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

def map_antibodies_to_profiles(adata, cell_profile_dict):
    """
    Map antibody capture data to predefined cell type profiles.

    Args:
        adata (AnnData): Antibody capture AnnData object.
        cell_profile_dict (dict): Dictionary mapping cell types to antibody markers.

    Returns:
        np.ndarray: Profile-based antibody data matrix (N_spots x T_cell_types).
        list: List of cell type names (to ensure column order).
    """
    # Step 1: Subset data to relevant markers
    all_markers = [marker for profile in cell_profile_dict.values() for marker in profile['Major']]
    existing_markers = [marker for marker in all_markers if marker in adata.var_names]

    if len(existing_markers) == 0:
        print("Adata variables: ", adata.var_names)
        print("Antibody markers: ", all_markers)
        raise ValueError("No matching antibody markers found in adata.var_names.")
    
    adata.var_names_make_unique()
    
    adata = adata[:, existing_markers]

    # Step 2: Extract and prepare antibody capture data
    antibody_capture_data = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    antibody_capture_var_names = np.array(adata.var_names)

    cell_type_names = list(cell_profile_dict.keys())  # Cell type order
    N = antibody_capture_data.shape[0]  # Number of spots
    T = len(cell_type_names)  # Number of cell types

    profile_based_antibody_data = np.zeros((N, T))

    for profile_idx, (profile_name, profile_markers) in enumerate(cell_profile_dict.items()):
        major_markers = profile_markers.get("Major", [])
        relevant_marker_indices = [
            np.where(antibody_capture_var_names == marker)[0][0]
            for marker in major_markers if marker in antibody_capture_var_names
        ]
        if relevant_marker_indices:
            profile_based_antibody_data[:, profile_idx] = antibody_capture_data[:, relevant_marker_indices].mean(axis=1)
        else:
            print(f"⚠️ No valid markers found for profile '{profile_name}'.")

    # Step 3: Normalize columns to prevent zero-division
    column_max = np.max(profile_based_antibody_data, axis=0)
    zero_columns = column_max == 0
    if np.any(zero_columns):
        print(f"⚠️ Warning: Zero columns detected. Adding epsilon to prevent NaNs.")
        column_max[zero_columns] = 1e-6

    profile_based_antibody_data /= column_max

    # Validate final data
    if np.isnan(profile_based_antibody_data).any():
        raise ValueError("NaN values detected in `profile_based_antibody_data` after mapping.")

    return profile_based_antibody_data, cell_type_names


def optimize_cell_proportions(profile_based_antibody_data, cell_type_names, tolerance=1e-4, max_iterations=50, lambda_reg=1, alpha=0.5):
    """
    Perform EM-based optimization for cell type proportions using Gurobi.

    Args:
        profile_based_antibody_data (np.ndarray): N x T matrix of mapped antibody data.
        cell_type_names (list): List of cell type names.
        tolerance (float): Convergence tolerance for EM algorithm.
        max_iterations (int): Maximum number of iterations.
        lambda_reg (float): Regularization strength.
        alpha (float): L1-L2 tradeoff.

    Returns:
        pd.DataFrame: Cell type proportions per spot.
    """
    N, T = profile_based_antibody_data.shape
    
    # Initialize beta estimates
    beta_estimates = {ct: 1.0 for ct in cell_type_names}
    beta_prev = np.zeros(T)
    Y_prev = np.zeros((N, T))
    iteration = 0

    while iteration < max_iterations:
        print(f"\nIteration {iteration + 1}")
        model = gp.Model("EM_Cell_Proportions")
        model.setParam('OutputFlag', 0)
        
        # Define variables Y[i, j]
        Y = model.addVars(N, T, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="Y")
        
        # Objective: Total squared error + Elastic Net regularization
        error_terms = []
        for i in range(N):
            for j in range(T):
                S_ij = profile_based_antibody_data[i, j]
                beta_j = beta_estimates[cell_type_names[j]]
                Y_ij = Y[i, j]
                error_terms.append((S_ij - beta_j * Y_ij) * (S_ij - beta_j * Y_ij))

        total_error = gp.quicksum(error_terms)
        l1_term = gp.quicksum(Y[i, j] for i in range(N) for j in range(T))
        l2_term = gp.quicksum(Y[i, j] * Y[i, j] for i in range(N) for j in range(T))
        regularization_term = lambda_reg * (alpha * l1_term + (1 - alpha) * l2_term)
        
        model.setObjective(total_error + regularization_term, GRB.MINIMIZE)
        
        for i in range(N):
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) >= 0.95)
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) <= 1.05)
        
        model.optimize()
        
        if model.status == GRB.OPTIMAL:
            Y_values = np.array([[Y[i, j].X for j in range(T)] for i in range(N)])
        else:
            raise ValueError("Gurobi optimization failed to converge.")
        
        beta_new = np.array([
            np.dot(profile_based_antibody_data[:, j], Y_values[:, j]) / np.dot(Y_values[:, j], Y_values[:, j])
            if np.dot(Y_values[:, j], Y_values[:, j]) > 0 else 0.0
            for j in range(T)
        ])
        
        beta_diff = np.linalg.norm(beta_new - beta_prev)
        Y_diff = np.linalg.norm(Y_values - Y_prev)
        
        print(f"Change in beta: {beta_diff:.6f}, Change in Y: {Y_diff:.6f}")
        if beta_diff < tolerance and Y_diff < tolerance:
            print("Convergence achieved.")
            break
        
        beta_prev = beta_new
        Y_prev = Y_values
        iteration += 1

    return Y_values


def deconvolute_spot_with_neighbors(
    i, 
    filtered_adata, 
    cell_type_numbers_array, 
    radius, 
    alpha=0.5, 
    lambda_reg_gex=0.0001,
    spatial_smoothing_weight=0.1
):
    """
    Deconvolute gene expression incorporating negative binomial characteristics
    and spatial smoothing between neighboring spots.
    """
    # Validate parameters
    if not (0 <= alpha <= 1):
        raise ValueError("alpha must be between 0 and 1")
    if lambda_reg_gex < 0:
        raise ValueError("lambda_reg_gex must be non-negative")
    if spatial_smoothing_weight < 0:
        raise ValueError("spatial_smoothing_weight must be non-negative")

    model = None
    try:
        logging.info(f"Starting deconvolution for spot {i}")

        # Get Neighborhood Indices
        neighborhood_indices = get_neighbors_with_fixed_radius(i, filtered_adata, radius=radius, include_center=True)
        if not neighborhood_indices:
            logging.error(f"No valid neighbors found for spot {i}.")
            return None

        # Ensure indices are integers
        neighborhood_indices = [int(idx) for idx in neighborhood_indices if isinstance(idx, (int, np.integer))]
        neighborhood_indices = np.array(neighborhood_indices, dtype=int)
        
        logging.debug(f"Neighborhood indices for spot {i}: {neighborhood_indices}")

        # Extract Gene Expression Data
        deconvolution_expression_data = filtered_adata.X
        if hasattr(deconvolution_expression_data, 'toarray'):
            deconvolution_expression_data = deconvolution_expression_data.toarray()
        
        # Get neighborhood data
        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]
        
        # Calculate size factors and normalize
        size_factors = np.sum(neighborhood_expression_data, axis=1, keepdims=True)
        epsilon = 1e-10
        size_factors = np.maximum(size_factors, epsilon)
        local_median_umi = np.median(size_factors)
        if local_median_umi < epsilon:
            local_median_umi = epsilon
        size_factor_normalized = neighborhood_expression_data / size_factors * local_median_umi

        # Calculate spatial weights based on distances
        spot_coordinates = filtered_adata.obsm['spatial'][neighborhood_indices]
        distances = np.zeros((len(neighborhood_indices), len(neighborhood_indices)))
        for idx1 in range(len(neighborhood_indices)):
            for idx2 in range(idx1 + 1, len(neighborhood_indices)):
                dist = np.linalg.norm(spot_coordinates[idx1] - spot_coordinates[idx2])
                # Gaussian kernel for distance weighting
                weight = np.exp(-dist / (radius/2))  # radius/2 as length scale
                distances[idx1, idx2] = weight
                distances[idx2, idx1] = weight

        # Calculate negative binomial statistics with safety checks
        neighborhood_mean = np.mean(size_factor_normalized, axis=0)
        neighborhood_var = np.var(size_factor_normalized, axis=0)
        
        # Avoid division by zero and handle edge cases in dispersion calculation
        epsilon_mean = 1e-10  # Small constant to prevent division by zero
        
        # Calculate dispersion only where mean is sufficiently large
        mean_mask = neighborhood_mean > epsilon_mean
        dispersion = np.ones_like(neighborhood_mean)  # Initialize with ones
        
        # Calculate dispersion only for genes with valid means
        valid_means = neighborhood_mean[mean_mask]
        valid_vars = neighborhood_var[mean_mask]
        
        # Calculate dispersion for valid genes
        dispersion[mean_mask] = np.maximum(
            0, 
            (valid_vars - valid_means) / (valid_means ** 2 + epsilon_mean)
        )
        
        # Handle edge cases
        dispersion[np.isnan(dispersion)] = 1.0  # fallback for numerical issues
        dispersion[np.isinf(dispersion)] = 1.0   # handle infinite values
        dispersion = np.minimum(dispersion, 10.0) # cap maximum dispersion
        
        # Add logging for debugging
        # Only log NaN issues if they occur
        if np.any(np.isnan(size_factor_normalized)):
            logging.error(f"NaN values detected after normalization for spot {i}")
            return None
        
        # Scale each gene to its maximum expression in the neighborhood
        gene_max = np.max(size_factor_normalized, axis=0, keepdims=True)
        gene_max[gene_max == 0] = 1  # Avoid division by zero
        neighborhood_expression_data = size_factor_normalized / gene_max
        
        # Normalize cell type proportions within each spot
        spot_total_cells = np.sum(neighborhood_cell_type_numbers, axis=1, keepdims=True)
        spot_total_cells[spot_total_cells == 0] = 1  # Avoid division by zero
        neighborhood_cell_type_numbers = neighborhood_cell_type_numbers / spot_total_cells

        # Build Gurobi Model
        model = gp.Model(f"gene_expression_proportion_deconvolution_spot_{i}")
        model.setParam('OutputFlag', 0)
        model.setParam('Threads', 1)
        model.setParam('NodefileStart', 0.5)
        model.setParam('MIPGap', 0.01)
        model.setParam('TimeLimit', 600)
        model.setParam('NodeLimit', 1000000)
        
        T = neighborhood_cell_type_numbers.shape[1]  # Number of cell types
        M = neighborhood_expression_data.shape[1]    # Number of genes
        
        # Define variables
        P = model.addVars(T, M, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cell_type_gene_proportion")
        error = model.addVars(len(neighborhood_indices), M, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="error")
        
        # # New variables for spatial smoothing
        # spatial_diff = model.addVars(
        #     len(neighborhood_indices), 
        #     len(neighborhood_indices), 
        #     T, M, 
        #     lb=0, 
        #     ub=GRB.INFINITY, 
        #     vtype=GRB.CONTINUOUS, 
        #     name="spatial_diff"
        # )

        # Proportion constraint for each gene
        wiggle_room = 0.1
        for k in range(M):
            model.addConstr(
                gp.quicksum(P[j, k] for j in range(T)) >= 1 - wiggle_room,
                name=f"proportion_lower_bound_gene_{k}"
            )
            model.addConstr(
                gp.quicksum(P[j, k] for j in range(T)) <= 1 + wiggle_room,
                name=f"proportion_upper_bound_gene_{k}"
            )

        # Add NB-weighted error terms
        for n, spot_idx in enumerate(neighborhood_indices):
            for k in range(M):
                total_proportion = np.sum(neighborhood_cell_type_numbers[n, :])
                predicted_gene_expression = gp.quicksum(
                    (neighborhood_cell_type_numbers[n, j] / total_proportion) * P[j, k] 
                    for j in range(T)
                )
                
                # Weight errors by inverse dispersion (higher dispersion = lower weight)
                nb_weight = 1.0 / (1.0 + dispersion[k])
                
                model.addConstr(
                    error[n, k] >= nb_weight * (neighborhood_expression_data[n, k] - predicted_gene_expression)
                )
                model.addConstr(
                    error[n, k] >= nb_weight * (predicted_gene_expression - neighborhood_expression_data[n, k])
                )

        # # Add spatial smoothing constraints
        # for n1 in range(len(neighborhood_indices)):
        #     for n2 in range(n1 + 1, len(neighborhood_indices)):
        #         if distances[n1, n2] > 0.01:  # Only consider meaningful distances
        #             for j in range(T):
        #                 for k in range(M):
        #                     # Get predicted expressions for both spots
        #                     pred1 = gp.quicksum(
        #                         (neighborhood_cell_type_numbers[n1, j] / np.sum(neighborhood_cell_type_numbers[n1, :])) 
        #                         * P[j, k]
        #                     )
        #                     pred2 = gp.quicksum(
        #                         (neighborhood_cell_type_numbers[n2, j] / np.sum(neighborhood_cell_type_numbers[n2, :])) 
        #                         * P[j, k]
        #                     )
                            
        #                     # Add weighted difference constraints
        #                     model.addConstr(
        #                         spatial_diff[n1, n2, j, k] >= distances[n1, n2] * (pred1 - pred2)
        #                     )
        #                     model.addConstr(
        #                         spatial_diff[n1, n2, j, k] >= distances[n1, n2] * (pred2 - pred1)
        #                     )

        # # Enhanced objective function with spatial smoothing
        # spatial_smoothing_term = gp.quicksum(
        #     spatial_diff[n1, n2, j, k] 
        #     for n1 in range(len(neighborhood_indices)) 
        #     for n2 in range(n1 + 1, len(neighborhood_indices))
        #     for j in range(T)
        #     for k in range(M)
        # )

        # L1 and L2 terms for elastic net regularization
        l1_term = gp.quicksum(P[j, k] for j in range(T) for k in range(M))
        l2_term = gp.quicksum(P[j, k] * P[j, k] for j in range(T) for k in range(M))

        # Enhanced objective function with spatial smoothing AND elastic net regularization
        model.setObjective(
            gp.quicksum(error[n, k] for n in range(len(neighborhood_indices)) for k in range(M)) +
            # spatial_smoothing_weight * spatial_smoothing_term +
            lambda_reg_gex * (alpha * l1_term + (1 - alpha) * l2_term),
            GRB.MINIMIZE
        )

        model.optimize()

        if model.SolCount > 0:
            logging.info(f"Solution found for spot {i}")
            gene_expression_profile_solution = {key: P[key].X for key in P}
            return np.array([[gene_expression_profile_solution[j, k] for k in range(M)] for j in range(T)])
        else:
            logging.error(f"No feasible solution found for spot {i}.")
            return None

    except Exception as e:
        logging.error(f"Error during deconvolution of spot {i}: {str(e)}")
        logging.error(traceback.format_exc())
        if model:
            del model
            del neighborhood_expression_data
            del size_factor_normalized
            del spatial_diff
        gc.collect()
        return None


def optimize_gene_expression(
    sample_name,
    deconvolution_expression_data,
    cell_type_numbers_array,
    filtered_adata,
    radius=2,
    alpha=0.5,
    lambda_reg_gex=0.0001,
    max_workers=None,
    checkpoint_interval=100,
    output_dir="checkpoints",
    rerun=False
):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    N = deconvolution_expression_data.shape[0]
    M = deconvolution_expression_data.shape[1]
    T = cell_type_numbers_array.shape[1]

    # If rerun is True, delete all existing files for this sample
    if rerun:
        existing_files = [
            f for f in os.listdir(output_dir) 
            if f.startswith(f"{sample_name}_gene_expression") and f.endswith(".npz")
        ]
        print("Found existing files: ", existing_files)
        for file in existing_files:
            try:
                os.remove(os.path.join(output_dir, file))
                logging.info(f"Deleted existing file: {file}")
            except Exception as e:
                logging.warning(f"Failed to delete {file}: {e}")
        logging.info(f"Starting fresh analysis for {sample_name} (rerun=True)")
        completed_spots = set()
        spotwise_gene_expression_profiles = {}
    
    # If not rerunning, check for completed run
    else:
        complete_file = os.path.join(output_dir, f"{sample_name}_gene_expression_complete.npz")
        if os.path.exists(complete_file):
            try:
                # Try to load the complete file to verify it's valid
                complete_data = np.load(complete_file)
                if 'profiles' in complete_data and 'completed_spots' in complete_data:
                    logging.info(f"Found completed analysis for {sample_name}")
                    profiles = complete_data['profiles']
                    if profiles.shape == (N, T, M):  # Verify correct dimensions
                        return {i: profiles[i] for i in range(N)}
            except Exception as e:
                logging.error(f"Error loading complete file: {e}")
                # If complete file is corrupt, delete all checkpoints
                existing_files = [
                    f for f in os.listdir(output_dir) 
                    if f.startswith(f"{sample_name}_gene_expression") and f.endswith(".npz")
                ]
                for file in existing_files:
                    try:
                        os.remove(os.path.join(output_dir, file))
                        logging.info(f"Deleted corrupted checkpoint: {file}")
                    except Exception as e:
                        logging.warning(f"Failed to delete {file}: {e}")
                    completed_spots = set()
                    spotwise_gene_expression_profiles = {}
        else:
            # Look for the latest checkpoint
            checkpoints = [
                f for f in os.listdir(output_dir)
                if f.startswith(f"{sample_name}_gene_expression_checkpoint_") and f.endswith(".npz")
            ]
            
            if checkpoints:
                checkpoint_numbers = [
                    int(f.replace(f"{sample_name}_gene_expression_checkpoint_", "").replace(".npz", ""))
                    for f in checkpoints
                ]
                latest_number = max(checkpoint_numbers)
                latest_checkpoint = os.path.join(
                    output_dir, 
                    f"{sample_name}_gene_expression_checkpoint_{latest_number}.npz"
                )
                
                try:
                    checkpoint_data = np.load(latest_checkpoint)
                    if 'profiles' in checkpoint_data and 'completed_spots' in checkpoint_data:
                        logging.info(f"Found latest checkpoint with {latest_number} completed spots")
                        profiles = checkpoint_data['profiles']
                        completed_spots = set(checkpoint_data['completed_spots'].tolist())
                        
                        if profiles.shape == (N, T, M):
                            logging.info(f"Loading checkpoint data with shape {profiles.shape}")
                            # Initialize with checkpoint data
                            spotwise_gene_expression_profiles = {
                                i: profiles[i] 
                                for i in completed_spots 
                                if not np.any(np.isnan(profiles[i]))
                            }
                            # Update completed_spots to only include valid profiles
                            completed_spots = set(spotwise_gene_expression_profiles.keys())
                            logging.info(f"Loaded {len(completed_spots)} valid profiles from checkpoint")
                        else:
                            logging.warning(f"Invalid profile shape in checkpoint: {profiles.shape}, expected ({N}, {T}, {M})")
                            completed_spots = set()
                            spotwise_gene_expression_profiles = {}
                except Exception as e:
                    logging.error(f"Error loading checkpoint file: {e}")
                    # If checkpoint is corrupt, delete all checkpoints
                    for file in checkpoints:
                        try:
                            os.remove(os.path.join(output_dir, file))
                            logging.info(f"Deleted corrupted checkpoint: {file}")
                        except Exception as e:
                            logging.warning(f"Failed to delete {file}: {e}")
                    completed_spots = set()
                    spotwise_gene_expression_profiles = {}
            else:
                completed_spots = set()
                spotwise_gene_expression_profiles = {}

    logging.info(f"Starting analysis for {sample_name}")
    logging.info(f"Already completed spots: {len(completed_spots)}")




    try:
        workers = max_workers if max_workers is not None else os.cpu_count()
        logging.info(f"Using {workers} workers")
        
        # Only process spots that haven't been completed
        remaining_spots = [i for i in range(N) if i not in completed_spots]
        logging.info(f"Processing {len(remaining_spots)} remaining spots")
        
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    deconvolute_spot_with_neighbors,
                    i,
                    filtered_adata,
                    cell_type_numbers_array,
                    radius,
                    alpha,
                    lambda_reg_gex,

                ): i for i in remaining_spots
            }

            with tqdm(total=len(remaining_spots), desc="Deconvoluting Remaining Spots") as pbar:
                spots_since_last_save = 0
                
                for future in as_completed(futures):
                    i = futures[future]
                    try:
                        result = future.result()
                        if result is not None:
                            # Verify result shape before storing
                            if result.ndim != 2:
                                logging.error(f"Unexpected result shape for spot {i}: {result.shape}")
                                continue
                                
                            spotwise_gene_expression_profiles[i] = result.copy()
                            completed_spots.add(i)
                            spots_since_last_save += 1

                            # Save checkpoint every checkpoint_interval spots
                            if spots_since_last_save >= checkpoint_interval:
                                # Use number of completed spots instead of current spot index
                                n_completed = len(completed_spots)
                                checkpoint_path = os.path.join(
                                    output_dir, 
                                    f"{sample_name}_gene_expression_checkpoint_{n_completed}.npz"
                                )
                                
                                # Convert dictionary to numpy array for saving
                                max_spot = max(spotwise_gene_expression_profiles.keys())
                                profiles_array = np.full((max_spot + 1, T, M), np.nan)
                                
                                for spot_idx, profile in spotwise_gene_expression_profiles.items():
                                    if profile.shape != (T, M):
                                        logging.error(f"Invalid profile shape at spot {spot_idx}: {profile.shape}")
                                        continue
                                    profiles_array[spot_idx] = profile
                                
                                try:
                                    # Save as compressed numpy array with completed spots info
                                    np.savez_compressed(
                                        checkpoint_path,
                                        profiles=profiles_array,
                                        completed_spots=np.array(list(completed_spots)),
                                        n_completed=n_completed
                                    )
                                    
                                    # Only delete old checkpoints if new one saved successfully
                                    existing_checkpoints = [
                                        f for f in os.listdir(output_dir) 
                                        if f.startswith(f"{sample_name}_gene_expression_checkpoint_") 
                                        and f.endswith(".npz")
                                        and f != os.path.basename(checkpoint_path)  # Don't delete the one we just created
                                    ]
                                    for old_checkpoint in existing_checkpoints:
                                        try:
                                            os.remove(os.path.join(output_dir, old_checkpoint))
                                        except Exception as e:
                                            logging.warning(f"Failed to delete old checkpoint {old_checkpoint}: {e}")
                                    
                                    logging.info(f"Saved checkpoint after {n_completed} completed spots")
                                    spots_since_last_save = 0
                                except Exception as e:
                                    logging.error(f"Failed to save checkpoint: {e}")

                    except Exception as e:
                        logging.error(f"Error in spot {i}: {str(e)}")
                        logging.error(traceback.format_exc())
                    pbar.update(1)

    except Exception as e:
        logging.error(f"Error during parallel processing: {str(e)}")
        logging.error(traceback.format_exc())
    finally:
        futures.clear()
        gc.collect()
        
        # Save final results
        final_path = os.path.join(output_dir, f"{sample_name}_gene_expression_complete.npz")
        max_spot = max(spotwise_gene_expression_profiles.keys())
        final_profiles = np.full((max_spot + 1, T, M), np.nan)
        for spot_idx, profile in spotwise_gene_expression_profiles.items():
            final_profiles[spot_idx] = profile
            
        np.savez_compressed(
            final_path, 
            profiles=final_profiles, 
            completed_spots=np.array(list(completed_spots))
        )
        logging.info(f"Saved final results with {len(completed_spots)} completed spots")

    return spotwise_gene_expression_profiles