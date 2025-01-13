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
import numpy as np
import math
from sklearn.metrics import mutual_info_score


def compute_global_prior(
    spotwise_gene_expression_profiles,
    cell_type_numbers_array,
    lambda_prior=1.0,
    num_bins=5
):
    """
    Compute a global prior for each (celltype, gene) pair based on mutual information
    between the gene's assigned usage (from the first pass) and the cell type's proportion
    across all spots.

    Args:
        spotwise_gene_expression_profiles (dict[int -> (T x M) array]):
            A dictionary from spot index to a (T x M) matrix with the *first-pass* assignment
            of each gene among T cell types.
        cell_type_numbers_array (np.ndarray):
            N x T array of cell-type proportions (or counts) across N spots.
        lambda_prior (float):
            A scaling factor used to transform raw MI into a sharper or flatter distribution.
        num_bins (int):
            Number of bins to use when discretizing continuous values for MI calculation.

    Returns:
        global_prior (np.ndarray):
            A (T x M) matrix in [0,1], where global_prior[t, m] indicates how likely
            gene m belongs to cell type t based on global mutual information.
        mi_matrix (np.ndarray):
            The raw mutual information values (T x M) before transformation.
    """
    # Gather spot indices and shapes
    spot_indices = sorted(spotwise_gene_expression_profiles.keys())
    if not spot_indices:
        raise ValueError("No spotwise profiles found. Did the first pass run correctly?")
    T, M = spotwise_gene_expression_profiles[spot_indices[0]].shape
    N = len(spot_indices)

    # usage_array[n, t, m] = fraction of gene m assigned to cell type t at spot n
    usage_array = np.zeros((N, T, M), dtype=float)
    for idx, spot_id in enumerate(spot_indices):
        usage_array[idx, :, :] = spotwise_gene_expression_profiles[spot_id]

    # Helper for discretization
    def discretize(values, bins):
        return np.digitize(values, bins=np.histogram_bin_edges(values, bins=bins))

    mi_matrix = np.zeros((T, M), dtype=float)

    # For each cell type t, discretize the cell_type_numbers_array[:, t]
    # and compute MI with usage_array[:, t, m] for each gene m
    for t_idx in range(T):
        ct_disc = discretize(cell_type_numbers_array[:, t_idx], bins=num_bins)

        for m_idx in range(M):
            usage_vec = usage_array[:, t_idx, m_idx]
            usage_disc = discretize(usage_vec, bins=num_bins)
            mi_val = mutual_info_score(usage_disc, ct_disc)
            mi_matrix[t_idx, m_idx] = mi_val

    # Convert MI to [0,1] prior via softmax across t for each gene m
    global_prior = np.zeros((T, M), dtype=float)
    for m_idx in range(M):
        mis = mi_matrix[:, m_idx]
        scaled = mis * lambda_prior
        # subtract max(scaled) for numerical stability
        exps = np.exp(scaled - np.max(scaled))
        denom = np.sum(exps)
        if denom < 1e-12:
            # If everything is near 0, fallback to uniform
            global_prior[:, m_idx] = 1.0 / T
        else:
            global_prior[:, m_idx] = exps / denom

    return global_prior, mi_matrix


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

def deconvolute_spot_with_neighbors_with_prior(
    i, 
    filtered_adata, 
    cell_type_numbers_array, 
    radius,
    global_prior,
    lambda_prior_weight=0.5,
    alpha=0.5, 
    lambda_reg_gex=0.0001,
    local_weight=0.5,
    global_weight=0.5
):
    """
    Same structure as `deconvolute_spot_with_neighbors`, but includes a penalty term
    that discourages assigning gene k to cell type j if global_prior[j,k] is small.
    """
    model = None
    try:
        if not (0 <= alpha <= 1):
            raise ValueError("alpha must be between 0 and 1")
        if lambda_reg_gex < 0:
            raise ValueError("lambda_reg_gex must be non-negative")
        if local_weight < 0 or global_weight < 0:
            raise ValueError("local_weight and global_weight must be non-negative")

        # Neighborhood
        neighborhood_indices = get_neighbors_with_fixed_radius(i, filtered_adata, radius=radius, include_center=True)
        if not neighborhood_indices:
            logging.error(f"No valid neighbors found for spot {i}.")
            return None

        neighborhood_indices = np.array([
            int(idx) for idx in neighborhood_indices
            if isinstance(idx, (int, np.integer))
        ], dtype=int)

        # Extract local expression
        deconvolution_expression_data = filtered_adata.X
        if hasattr(deconvolution_expression_data, 'toarray'):
            deconvolution_expression_data = deconvolution_expression_data.toarray()

        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]

        epsilon = 1e-10
        size_factors = np.maximum(
            np.sum(neighborhood_expression_data, axis=1, keepdims=True),
            epsilon
        )
        local_median_umi = max(np.median(size_factors), epsilon)
        size_factor_normalized = neighborhood_expression_data / size_factors * local_median_umi

        # NB-like dispersion
        neighborhood_mean = np.mean(size_factor_normalized, axis=0)
        neighborhood_var = np.var(size_factor_normalized, axis=0)
        mean_mask = neighborhood_mean > epsilon
        dispersion = np.ones_like(neighborhood_mean)
        dispersion[mean_mask] = np.maximum(
            0,
            (neighborhood_var[mean_mask] - neighborhood_mean[mean_mask])
            / (neighborhood_mean[mean_mask] ** 2 + epsilon)
        )
        dispersion = np.clip(dispersion, 0, 10.0)

        gene_max = np.maximum(np.max(size_factor_normalized, axis=0, keepdims=True), epsilon)
        neighborhood_expression_data = size_factor_normalized / gene_max

        # Local & Global Enrichment
        local_proportions = np.mean(neighborhood_cell_type_numbers, axis=0)
        global_proportions = np.mean(cell_type_numbers_array, axis=0)

        local_enrichment_scores = local_proportions / (global_proportions + epsilon)
        local_enrichment_scores /= (np.max(local_enrichment_scores) + epsilon)

        # For demonstration, we reuse the same logic for "global_enrichment_scores"
        global_enrichment_scores = local_enrichment_scores.copy()

        # Build Gurobi model
        model = gp.Model(f"gene_expression_deconvolution_spot_{i}_with_prior")
        model.setParam('OutputFlag', 0)
        model.setParam('Threads', 1)
        model.setParam('NodefileStart', 0.5)
        model.setParam('MIPGap', 0.01)
        model.setParam('TimeLimit', 600)
        model.setParam('NodeLimit', 1000000)

        T = neighborhood_cell_type_numbers.shape[1]
        M = neighborhood_expression_data.shape[1]

        # Decision variables
        P = model.addVars(T, M, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="P")
        error = model.addVars(len(neighborhood_indices), M, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="error")

        # Proportion constraints
        wiggle_room = 0.05
        for k in range(M):
            model.addConstr(
                gp.quicksum(P[j, k] for j in range(T)) >= 1 - wiggle_room,
                name=f"prop_lower_gene_{k}"
            )
            model.addConstr(
                gp.quicksum(P[j, k] for j in range(T)) <= 1 + wiggle_room,
                name=f"prop_upper_gene_{k}"
            )

        # Predicted expression
        for n, spot_idx in enumerate(neighborhood_indices):
            for k in range(M):
                predicted_expression = gp.quicksum(
                    (
                        local_weight * local_enrichment_scores[j]
                        + global_weight * global_enrichment_scores[j]
                    )
                    * neighborhood_cell_type_numbers[n, j]
                    * P[j, k]
                    for j in range(T)
                )
                nb_weight = 1.0 / (1.0 + dispersion[k])
                model.addConstr(error[n, k] >= nb_weight * (neighborhood_expression_data[n, k] - predicted_expression))
                model.addConstr(error[n, k] >= nb_weight * (predicted_expression - neighborhood_expression_data[n, k]))

        # Elastic Net terms
        l1_term = gp.quicksum(P[j, k] for j in range(T) for k in range(M))
        l2_term = gp.quicksum(P[j, k] * P[j, k] for j in range(T) for k in range(M))
        reconstruction_error = gp.quicksum(error[n, k] for n in range(len(neighborhood_indices)) for k in range(M))

        # Global Prior Penalty
        #   sum_{j,k} [lambda_prior_weight * (1 - global_prior[j,k]) * P[j,k]]
        prior_penalty = gp.quicksum(
            lambda_prior_weight * (1.0 - global_prior[j, k]) * P[j, k]
            for j in range(T) for k in range(M)
        )

        # Final objective
        model.setObjective(
            reconstruction_error + 
            lambda_reg_gex * (alpha * l1_term + (1 - alpha) * l2_term) + 
            lambda_prior_weight * prior_penalty,
            GRB.MINIMIZE
        )

        model.optimize()
        if model.SolCount > 0:
            gene_expression_profile_solution = {key: P[key].X for key in P}
            return np.array([
                [gene_expression_profile_solution[j, k] for k in range(M)]
                for j in range(T)
            ])
        else:
            logging.error(f"No feasible solution found for spot {i}.")
            return None

    except Exception as e:
        logging.error(f"Error during deconvolution of spot {i} with prior: {str(e)}")
        logging.error(traceback.format_exc())
        return None
    finally:
        if model:
            del model
        gc.collect()


def deconvolute_spot_with_neighbors_wrapper(args):
    """
    Wrapper function for parallel processing that unpacks arguments and calls the main function.
    Args should be a tuple of (i, filtered_adata, cell_type_numbers_array, radius, global_prior, 
    lambda_prior_weight, alpha, lambda_reg_gex)
    """
    (i, filtered_adata, cell_type_numbers_array, radius, global_prior, 
     lambda_prior_weight, alpha, lambda_reg_gex) = args
    
    return deconvolute_spot_with_neighbors_with_prior(
        i,
        filtered_adata,
        cell_type_numbers_array,
        radius,
        global_prior,
        lambda_prior_weight,
        alpha,
        lambda_reg_gex
    )

def two_pass_optimize_gene_expression(
    sample_name,
    deconvolution_expression_data,
    cell_type_numbers_array,
    filtered_adata,
    radius=2,
    alpha=0.5,
    lambda_reg_gex=0.0001,
    lambda_prior_weight=0.5,
    prior_change_threshold=0.01,
    max_workers=None,
    checkpoint_interval=100,
    output_dir="checkpoints",
    rerun=False
):
    """
    1) First pass with normal optimization (no global prior).
    2) Compute global prior from first pass using mutual information.
    3) Second pass with prior penalty in the objective.
    4) Optionally check how much the new prior changed from first to second pass,
       and decide if we want to stop or do more passes.
    """

    #########################
    # === PASS 1 (No Prior)
    #########################
    logging.info("=== Pass 1: standard local optimization ===")
    pass1_result = optimize_gene_expression(
        sample_name=sample_name+"_pass1",
        deconvolution_expression_data=deconvolution_expression_data,
        cell_type_numbers_array=cell_type_numbers_array,
        filtered_adata=filtered_adata,
        radius=radius,
        alpha=alpha,
        lambda_reg_gex=lambda_reg_gex,
        max_workers=max_workers,
        checkpoint_interval=checkpoint_interval,
        output_dir=output_dir,
        rerun=rerun
    )

    # pass1_result is a dict: spot_index -> (T x M) array

    #########################
    # Compute Global Prior from Pass 1
    #########################
    logging.info("=== Computing global prior from pass 1 ===")
    global_prior_pass1, mi_matrix_pass1 = compute_global_prior(
        pass1_result,
        cell_type_numbers_array,
        lambda_prior=1.0,  # or some "temperature" factor for the MI-based prior
        num_bins=5
    )

    #########################
    # === PASS 2 (With Prior)
    #########################
    logging.info("=== Pass 2: re-optimization with global prior penalty ===")

    # Create argument tuples for each spot
    spot_args = [
        (i, filtered_adata, cell_type_numbers_array, radius, global_prior_pass1,
         lambda_prior_weight, alpha, lambda_reg_gex)
        for i in range(deconvolution_expression_data.shape[0])
    ]

    # Use the wrapper function with ProcessPoolExecutor
    pass2_result = optimize_gene_expression_with_custom_spot_fn(
        sample_name=sample_name+"_pass2",
        deconvolution_expression_data=deconvolution_expression_data,
        cell_type_numbers_array=cell_type_numbers_array,
        filtered_adata=filtered_adata,
        radius=radius,
        alpha=alpha,
        lambda_reg_gex=lambda_reg_gex,
        spot_function=deconvolute_spot_with_neighbors_wrapper,
        spot_args=spot_args,  # Pass the prepared arguments
        max_workers=max_workers,
        checkpoint_interval=checkpoint_interval,
        output_dir=output_dir,
        rerun=rerun
    )

    # Compute a new prior from pass 2 result, to see if it changed significantly
    global_prior_pass2, mi_matrix_pass2 = compute_global_prior(
        pass2_result,
        cell_type_numbers_array,
        lambda_prior=1.0,
        num_bins=5
    )

    # Calculate how much the prior changed (L1 norm or relative difference)
    prior_diff = np.mean(np.abs(global_prior_pass2 - global_prior_pass1))
    logging.info(f"Prior difference between pass1 -> pass2: {prior_diff:.5f}")
    logging.info(f"Prior convergence threshold: {prior_change_threshold}")

    if prior_diff < prior_change_threshold:
        logging.info("Global prior converged sufficiently; stopping after pass 2.")
        final_result = pass2_result
    else:
        logging.warning(f"Global prior changed significantly ({prior_diff:.5f} > {prior_change_threshold})")
        logging.info("Accepting pass 2 results, but you may want to consider additional passes.")
        final_result = pass2_result

    return final_result

def optimize_gene_expression_with_custom_spot_fn(
    sample_name,
    deconvolution_expression_data,
    cell_type_numbers_array,
    filtered_adata,
    radius=2,
    alpha=0.5,
    lambda_reg_gex=0.0001,
    spot_function=None,
    spot_args=None,
    max_workers=None,
    checkpoint_interval=100,
    output_dir="checkpoints",
    rerun=False
):
    """
    Similar to `optimize_gene_expression`, but calls `spot_function(...)` instead of
    `deconvolute_spot_with_neighbors(...)` for each spot. The rest of the checkpoint
    and parallel logic remains the same.

    Args:
        sample_name (str): Unique name for this sample (used in checkpoint files).
        deconvolution_expression_data (np.ndarray): (N x M) gene expression data.
        cell_type_numbers_array (np.ndarray): (N x T) cell type proportions or counts.
        filtered_adata (AnnData): Anndata with .X for expression and .obsm['spatial'] for coords.
        radius (float): Neighborhood radius.
        alpha (float): L1-L2 regularization mixing param.
        lambda_reg_gex (float): Regularization strength for the Gurobi model.
        spot_function (callable): The function used to deconvolute a single spot. 
            Must accept: (i, filtered_adata, cell_type_numbers_array, radius, alpha, lambda_reg_gex)
            and return a (T x M) np.array or None.
        max_workers (int): Number of parallel workers.
        checkpoint_interval (int): Save frequency in # of spots processed.
        output_dir (str): Directory for checkpoint files.
        rerun (bool): If True, delete all existing checkpoint files.

    Returns:
        dict: spot_index -> (T x M) array of gene expression profiles per spot.
    """
    import os
    from pathlib import Path
    import logging
    import traceback
    import gc
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from tqdm import tqdm
    import numpy as np

    if spot_function is None:
        raise ValueError("`spot_function` must be provided.")

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
    else:
        # If not rerunning, check for completed run
        complete_file = os.path.join(output_dir, f"{sample_name}_gene_expression_complete.npz")
        if os.path.exists(complete_file):
            try:
                # Try to load the complete file to verify it's valid
                complete_data = np.load(complete_file)
                if 'profiles' in complete_data and 'completed_spots' in complete_data:
                    logging.info(f"Found completed analysis for {sample_name}")
                    profiles = complete_data['profiles']
                    if profiles.shape == (N, T, M):  # Verify correct dims
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
                            completed_spots = set(spotwise_gene_expression_profiles.keys())
                            logging.info(f"Loaded {len(completed_spots)} valid profiles from checkpoint")
                        else:
                            logging.warning(f"Invalid profile shape in checkpoint: {profiles.shape}, expected ({N}, {T}, {M})")
                            completed_spots = set()
                            spotwise_gene_expression_profiles = {}
                    else:
                        logging.warning("Checkpoint data is missing 'profiles' or 'completed_spots'")
                        completed_spots = set()
                        spotwise_gene_expression_profiles = {}
                except Exception as e:
                    logging.error(f"Error loading checkpoint file: {e}")
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
        
        remaining_spots = [i for i in range(N) if i not in completed_spots]
        logging.info(f"Processing {len(remaining_spots)} remaining spots")

        with ProcessPoolExecutor(max_workers=workers) as executor:
            # Submit jobs using the args directly
            futures = {
                executor.submit(spot_function, args): i 
                for i, args in enumerate(spot_args)
            }

            with tqdm(total=len(remaining_spots), desc="Deconvoluting Remaining Spots") as pbar:
                spots_since_last_save = 0
                for future in as_completed(futures):
                    i = futures[future]
                    try:
                        result = future.result()
                        if result is not None:
                            if result.ndim != 2:
                                logging.error(f"Unexpected result shape for spot {i}: {result.shape}")
                                continue
                            if result.shape != (T, M):
                                logging.error(f"Result shape {result.shape} doesn't match (T={T}, M={M})")
                                continue
                            spotwise_gene_expression_profiles[i] = result.copy()
                            completed_spots.add(i)
                            spots_since_last_save += 1

                            # Save checkpoint every checkpoint_interval spots
                            if spots_since_last_save >= checkpoint_interval:
                                n_completed = len(completed_spots)
                                checkpoint_path = os.path.join(
                                    output_dir, 
                                    f"{sample_name}_gene_expression_checkpoint_{n_completed}.npz"
                                )
                                max_spot = max(spotwise_gene_expression_profiles.keys())
                                profiles_array = np.full((max_spot + 1, T, M), np.nan)
                                for spot_idx, profile in spotwise_gene_expression_profiles.items():
                                    profiles_array[spot_idx] = profile
                                
                                try:
                                    np.savez_compressed(
                                        checkpoint_path,
                                        profiles=profiles_array,
                                        completed_spots=np.array(list(completed_spots)),
                                        n_completed=n_completed
                                    )
                                    # Clean up old checkpoints
                                    existing_checkpoints = [
                                        f for f in os.listdir(output_dir) 
                                        if f.startswith(f"{sample_name}_gene_expression_checkpoint_") 
                                        and f.endswith(".npz")
                                        and f != os.path.basename(checkpoint_path)
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
        # End of parallel region
        gc.collect()

        final_path = os.path.join(output_dir, f"{sample_name}_gene_expression_complete.npz")
        if len(spotwise_gene_expression_profiles) > 0:
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
        else:
            logging.info("No profiles to save (empty result).")

    return spotwise_gene_expression_profiles


def deconvolute_spot_with_neighbors(
    i, 
    filtered_adata, 
    cell_type_numbers_array, 
    radius, 
    alpha=0.5, 
    lambda_reg_gex=0.0001,
    local_weight=0.5,   # Weight for local neighborhood enrichment
    global_weight=0.5   # Weight for global dataset enrichment
):
    """
    Deconvolution incorporating local vs. global cell type enrichment patterns.
    Removing the concept of a separate 'broader' region entirely:
      - 'local' means immediate neighborhood
      - 'global' means entire dataset
    If local_weight + global_weight < 1, you could treat the leftover as a baseline factor.
    """
    model = None
    try:
        # ------------------------------------------------------------------
        # Validate user parameters
        # ------------------------------------------------------------------
        if not (0 <= alpha <= 1):
            raise ValueError("alpha must be between 0 and 1")
        if lambda_reg_gex < 0:
            raise ValueError("lambda_reg_gex must be non-negative")
        if local_weight < 0 or global_weight < 0:
            raise ValueError("local_weight and global_weight must be non-negative")

        logging.info(f"Starting deconvolution for spot {i}")

        # ------------------------------------------------------------------
        # 1) Identify the Local Neighborhood
        # ------------------------------------------------------------------
        neighborhood_indices = get_neighbors_with_fixed_radius(
            i, filtered_adata, radius=radius, include_center=True
        )
        if not neighborhood_indices:
            logging.error(f"No valid neighbors found for spot {i}.")
            return None

        neighborhood_indices = np.array([
            int(idx) for idx in neighborhood_indices 
            if isinstance(idx, (int, np.integer))
        ], dtype=int)

        # ------------------------------------------------------------------
        # 2) Extract Data for the Neighborhood
        # ------------------------------------------------------------------
        deconvolution_expression_data = filtered_adata.X
        if hasattr(deconvolution_expression_data, 'toarray'):
            deconvolution_expression_data = deconvolution_expression_data.toarray()

        # Local expression and cell-type counts
        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]

        # ------------------------------------------------------------------
        # 3) Size-Factor Normalization (local)
        # ------------------------------------------------------------------
        epsilon = 1e-10
        size_factors = np.maximum(
            np.sum(neighborhood_expression_data, axis=1, keepdims=True),
            epsilon
        )
        local_median_umi = max(np.median(size_factors), epsilon)
        size_factor_normalized = neighborhood_expression_data / size_factors * local_median_umi

        # ------------------------------------------------------------------
        # 4) Compute NB-Like Dispersion & Scale Expression
        # ------------------------------------------------------------------
        neighborhood_mean = np.mean(size_factor_normalized, axis=0)
        neighborhood_var = np.var(size_factor_normalized, axis=0)
        mean_mask = neighborhood_mean > epsilon
        dispersion = np.ones_like(neighborhood_mean)
        dispersion[mean_mask] = np.maximum(
            0,
            (neighborhood_var[mean_mask] - neighborhood_mean[mean_mask])
            / (neighborhood_mean[mean_mask] ** 2 + epsilon)
        )
        dispersion = np.clip(dispersion, 0, 10.0)

        # Scale neighborhood gene expression to [0, 1]
        gene_max = np.maximum(np.max(size_factor_normalized, axis=0, keepdims=True), epsilon)
        neighborhood_expression_data = size_factor_normalized / gene_max

        # ------------------------------------------------------------------
        # 5) Local & Global Cell-Type Enrichment
        # ------------------------------------------------------------------
        # Local proportions for the neighborhood
        local_proportions = np.mean(neighborhood_cell_type_numbers, axis=0)

        # Global proportions for entire dataset
        global_proportions = np.mean(cell_type_numbers_array, axis=0)

        # Convert local/global proportions into "enrichment" scores
        local_enrichment_scores = local_proportions / (global_proportions + epsilon)
        local_enrichment_scores /= (np.max(local_enrichment_scores) + epsilon)

        # You might define another measure for "global" as well. 
        # Here, a direct approach: 
        global_enrichment_scores = local_proportions / (global_proportions + epsilon)
        global_enrichment_scores /= (np.max(global_enrichment_scores) + epsilon)

        # ------------------------------------------------------------------
        # 6) Build the Gurobi Model
        # ------------------------------------------------------------------
        model = gp.Model(f"gene_expression_deconvolution_spot_{i}")
        model.setParam('OutputFlag', 0)
        model.setParam('Threads', 1)
        model.setParam('NodefileStart', 0.5)
        model.setParam('MIPGap', 0.01)
        model.setParam('TimeLimit', 600)
        model.setParam('NodeLimit', 1000000)

        T = neighborhood_cell_type_numbers.shape[1]
        M = neighborhood_expression_data.shape[1]

        # Define variables: P[j, k] for cell type j, gene k
        P = model.addVars(T, M, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="P")

        # Error variables: absolute deviation from observed
        error = model.addVars(
            len(neighborhood_indices), M, 
            lb=0, ub=GRB.INFINITY,
            vtype=GRB.CONTINUOUS,
            name="error"
        )

        # ------------------------------------------------------------------
        # 7) Proportion Constraints per Gene
        # ------------------------------------------------------------------
        wiggle_room = 0.05
        for k in range(M):
            model.addConstr(
                gp.quicksum(P[j, k] for j in range(T)) >= 1 - wiggle_room,
                name=f"prop_lower_gene_{k}"
            )
            model.addConstr(
                gp.quicksum(P[j, k] for j in range(T)) <= 1 + wiggle_room,
                name=f"prop_upper_gene_{k}"
            )

        # ------------------------------------------------------------------
        # 8) Predicted Expression Incorporating Local & Global Context
        # ------------------------------------------------------------------
        for n, spot_idx in enumerate(neighborhood_indices):
            for k in range(M):
                # Weighted sum of local + global enrichment for each cell type
                # If local_weight + global_weight < 1, leftover can be baseline
                predicted_expression = gp.quicksum(
                    (
                        local_weight * local_enrichment_scores[j]
                        + global_weight * global_enrichment_scores[j]
                    )
                    * neighborhood_cell_type_numbers[n, j]
                    * P[j, k]
                    for j in range(T)
                )

                # NB-weight: penalize large deviations more strongly
                nb_weight = 1.0 / (1.0 + dispersion[k])

                # Constrain error[n, k] to be >= the absolute difference
                model.addConstr(error[n, k] >= nb_weight * (
                    neighborhood_expression_data[n, k] - predicted_expression
                ))
                model.addConstr(error[n, k] >= nb_weight * (
                    predicted_expression - neighborhood_expression_data[n, k]
                ))

        # ------------------------------------------------------------------
        # 9) Elastic Net Regularization
        # ------------------------------------------------------------------
        l1_term = gp.quicksum(P[j, k] for j in range(T) for k in range(M))
        l2_term = gp.quicksum(P[j, k] * P[j, k] for j in range(T) for k in range(M))
        model.setObjective(
            gp.quicksum(error[n, k] for n in range(len(neighborhood_indices)) for k in range(M))
            + lambda_reg_gex * (alpha * l1_term + (1 - alpha) * l2_term),
            GRB.MINIMIZE
        )

        model.optimize()

        if model.SolCount > 0:
            logging.info(f"Solution found for spot {i}")
            gene_expression_profile_solution = {key: P[key].X for key in P}
            return np.array([
                [gene_expression_profile_solution[j, k] for k in range(M)]
                for j in range(T)
            ])
        else:
            logging.error(f"No feasible solution found for spot {i}.")
            return None

    except Exception as e:
        logging.error(f"Error during deconvolution of spot {i}: {str(e)}")
        logging.error(traceback.format_exc())

        # Debug statements to aid troubleshooting
        print("Debug Info:")
        print(f"Spot index: {i}")
        print("Local proportions:", local_proportions)
        print("Global proportions:", global_proportions)
        print("Local enrichment:", local_enrichment_scores)
        print("Global enrichment:", global_enrichment_scores)
        return None

    finally:
        if model:
            del model
        gc.collect()


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