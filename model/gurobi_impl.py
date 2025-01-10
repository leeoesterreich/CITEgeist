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
    global_correlations=None,
    global_confidence=None
):
    """
    Deconvolute gene expression for a single spot using neighborhood data and Gurobi optimization.
    Now uses pre-calculated global correlations with p-value-based weighting.
    """
    model = None
    try:
        logging.info(f"Starting deconvolution for spot {i}")

        # ✅ Step 1: Get Neighborhood Indices
        try:
            neighborhood_indices = get_neighbors_with_fixed_radius(i, filtered_adata, radius=radius, include_center=True)
            if not neighborhood_indices:
                logging.error(f"No valid neighbors found for spot {i}.")
                return None

            # Ensure indices are integers
            neighborhood_indices = [
                int(idx) for idx in neighborhood_indices if isinstance(idx, (int, np.integer))
            ]
            neighborhood_indices = np.array(neighborhood_indices, dtype=int)
        except Exception as e:
            logging.error(f"Failed to retrieve or validate neighborhood indices for spot {i}: {str(e)}")
            return None

        logging.debug(f"Neighborhood indices for spot {i}: {neighborhood_indices}")

        # ✅ Step 2: Extract Gene Expression Data
        try:
            deconvolution_expression_data = filtered_adata.X

            
            if hasattr(deconvolution_expression_data, 'toarray'):
                deconvolution_expression_data = deconvolution_expression_data.toarray()
            
            # 1. Get neighborhood data
            neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
            neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]
            
            # 2. Calculate size factors with safety checks
            size_factors = np.sum(neighborhood_expression_data, axis=1, keepdims=True)
            
            # Add small epsilon to avoid division by zero
            epsilon = 1e-10
            size_factors = np.maximum(size_factors, epsilon)
            
            # 3. Safe normalization
            local_median_umi = np.median(size_factors)
            if local_median_umi < epsilon:
                local_median_umi = epsilon
                
            size_factor_normalized = neighborhood_expression_data / size_factors * local_median_umi
            
            # 4. Check for invalid values
            if np.any(np.isnan(size_factor_normalized)) or np.any(np.isinf(size_factor_normalized)):
                logging.warning(f"Invalid values detected after normalization for spot {i}. Using fallback.")
                # Fallback to original data if normalization fails
                size_factor_normalized = neighborhood_expression_data
            
            # 5. Scale each gene to its maximum expression in the neighborhood
            # This preserves relative differences while making genes comparable
            gene_max = np.max(size_factor_normalized, axis=0, keepdims=True)
            gene_max[gene_max == 0] = 1  # Avoid division by zero
            neighborhood_expression_data = size_factor_normalized / gene_max
            
            # 6. Normalize cell type proportions within each spot
            spot_total_cells = np.sum(neighborhood_cell_type_numbers, axis=1, keepdims=True)
            spot_total_cells[spot_total_cells == 0] = 1  # Avoid division by zero
            neighborhood_cell_type_numbers = neighborhood_cell_type_numbers / spot_total_cells
            
            logging.debug(f"Neighborhood data shape for spot {i}: {neighborhood_expression_data.shape}")
            logging.debug(f"Mean expression after normalization: {np.mean(neighborhood_expression_data):.3f}")
            logging.debug(f"Std expression after normalization: {np.std(neighborhood_expression_data):.3f}")


        except IndexError as e:
            logging.error(f"IndexError for spot {i} during data extraction: {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected error during data extraction for spot {i}: {e}")
            return None

        # ✅ Step 3: Build Gurobi Model
        try:
            model = gp.Model(f"gene_expression_proportion_deconvolution_spot_{i}")
            model.setParam('OutputFlag', 0)
            model.setParam('Threads', 1)
            model.setParam('NodefileStart', 0.5)
            model.setParam('MIPGap', 0.01)
            model.setParam('TimeLimit', 600)
            model.setParam('NodeLimit', 1000000)
            
            T = neighborhood_cell_type_numbers.shape[1]  # Number of cell types
            M = neighborhood_expression_data.shape[1]  # Number of genes
            
            P = model.addVars(T, M, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cell_type_gene_proportion")
            error = model.addVars(len(neighborhood_indices), M, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="error")

            # ✅ Proportion constraint for each gene
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

            # ✅ Error terms and Elastic Net Regularization
            for n, spot_idx in enumerate(neighborhood_indices):
                for k in range(M):
                    total_proportion = np.sum(neighborhood_cell_type_numbers[n, :])
                    predicted_gene_expression = gp.quicksum(
                        (neighborhood_cell_type_numbers[n, j] / total_proportion) * P[j, k] for j in range(T)
                    )
                    model.addConstr(error[n, k] >= neighborhood_expression_data[n, k] - predicted_gene_expression)
                    model.addConstr(error[n, k] >= predicted_gene_expression - neighborhood_expression_data[n, k])

            lambda_reg = lambda_reg_gex
            l1_term = gp.quicksum(P[j, k] for j in range(T) for k in range(M))
            l2_term = gp.quicksum(P[j, k] * P[j, k] for j in range(T) for k in range(M))

            correlation_weight = 1
            min_proportion_threshold = 0.2  # Minimum proportion per spot
            min_expression_threshold = 5.0   # Minimum expression level
            
            # Modified correlation terms that:
            # 1. Uses all correlations (not just significant ones)
            # 2. Weights impact by statistical significance
            # 3. Handles positive and negative correlations differently
            correlation_terms = gp.quicksum(
                global_confidence[j, k] * (  # global_confidence already includes (1 - adj_pval)
                    # For positive correlations
                    (max(0, global_correlations[j, k]) * P[j, k] * -1) +  # Encourage assignment
                    # For negative correlations (with reduced penalty for less significant ones)
                    (min(0, global_correlations[j, k]) * P[j, k])  # Discourage assignment
                )
                for j in range(T) for k in range(M)
            )

            alpha = 0.8
            
            # Update objective function with weighted correlation term
            model.setObjective(
                gp.quicksum(error[n, k] for n in range(len(neighborhood_indices)) for k in range(M)) +
                lambda_reg * (alpha * l1_term + (1 - alpha) * l2_term) +
                correlation_weight * correlation_terms,
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
        finally:
            if model:
                del model
                gc.collect()
    except Exception as e:
        logging.error(f"Error during deconvolution of spot {i}: {str(e)}")
        logging.error(traceback.format_exc())
        if model:
            del model
        gc.collect()
        return None

def calculate_global_correlations(filtered_adata, cell_type_numbers_array, min_proportion_threshold=0.2, min_expression_threshold=5.0):
        """
        Calculate global correlations with confidence weighted by adjusted p-values.
        """
        T = cell_type_numbers_array.shape[1]
        M = filtered_adata.shape[1]
        
        correlations = np.zeros((T, M))
        correlation_confidence = np.zeros((T, M))
        
        expression_data = filtered_adata.X.toarray() if sp.issparse(filtered_adata.X) else filtered_adata.X
        
        for j in range(T):
            # Store all p-values for this cell type
            all_pvals = []
            all_corrs = []
            valid_indices = []
            
            # First pass: calculate all correlations and p-values
            for k in range(M):
                cell_type_props = cell_type_numbers_array[:, j]
                gene_expr = expression_data[:, k]
                
                valid_mask = (
                    ~(np.isnan(cell_type_props) | np.isnan(gene_expr)) & 
                    (cell_type_props >= min_proportion_threshold) &
                    (gene_expr >= min_expression_threshold)
                )
                filtered_props = cell_type_props[valid_mask]
                filtered_expr = gene_expr[valid_mask]
                
                if len(filtered_props) < 3 or np.std(filtered_expr) == 0:
                    correlations[j, k] = 0
                    correlation_confidence[j, k] = 0
                else:
                    corr, pval = spearmanr(filtered_props, filtered_expr)
                    if not np.isnan(corr) and not np.isnan(pval):
                        all_corrs.append(corr)
                        all_pvals.append(pval)
                        valid_indices.append(k)
                    else:
                        correlations[j, k] = 0
                        correlation_confidence[j, k] = 0
            
            # Adjust p-values using Benjamini-Hochberg FDR
            if all_pvals:
                adjusted_pvals = multipletests(all_pvals, method='fdr_bh')[1]
                
                # Second pass: store results using adjusted p-values
                for idx, k in enumerate(valid_indices):
                    corr = all_corrs[idx]
                    adj_pval = adjusted_pvals[idx]
                    
                    # Store correlation as is
                    correlations[j, k] = corr
                    
                    # Calculate confidence based on:
                    # 1. Statistical significance (1 - adj_pval)
                    # 2. Number of valid observations
                    # 3. Strength of correlation
                    confidence = (
                        #(len(filtered_props) / len(cell_type_props)) *  # Data coverage
                        (1 - adj_pval) *                                # Statistical significance
                        abs(corr)                                       # Correlation strength
                    )
                    correlation_confidence[j, k] = confidence
                    
                    # Log correlations with their weights
                    if abs(corr) > 0.3:  # Lower threshold to see more correlations
                        gene_name = filtered_adata.var_names[k]
                        logging.info(f"Global correlation found:")
                        logging.info(f"Cell type {j} with gene {gene_name}:")
                        logging.info(f"Correlation: {corr:.3f}, adj. p-value: {adj_pval:.3e}")
                        logging.info(f"Confidence weight: {confidence:.3f}")
                        logging.info(f"Valid spots: {len(filtered_props)} of {len(cell_type_props)}")
        
        return correlations, correlation_confidence

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

    # Calculate global correlations once before parallel processing
    logging.info("Calculating global correlations...")
    global_correlations, global_confidence = calculate_global_correlations(filtered_adata, cell_type_numbers_array)
    logging.info("Global correlations calculated.")

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
                    global_correlations,
                    global_confidence
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