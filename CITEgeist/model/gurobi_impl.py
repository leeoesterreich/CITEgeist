# Standard library imports
import os
import logging
import traceback
import gc
import concurrent
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Any, Optional, List, Tuple, Union
import json
import psutil
import time

# Third-party imports
import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import Model, GRB, quicksum
import scanpy as sc
import scipy
from scipy.stats import spearmanr
from scipy.optimize import minimize
from scipy.special import loggamma, digamma
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# Local imports
from .utils import get_neighbors_with_fixed_radius
from .checkpoints import CheckpointManager


def compute_global_prior(
    spotwise_gene_expression_profiles: Dict[int, np.ndarray],
    cell_type_numbers_array: np.ndarray,
    lambda_prior: float = 1.0,
    min_expression_threshold: float = 0.1
) -> Dict[str, Any]:
    """
    Compute global prior from pass 1 results using normalized expression patterns.
    
    Args:
        spotwise_gene_expression_profiles: Dictionary mapping spot indices to profile matrices
        cell_type_numbers_array: Array of cell type proportions (N_spots × T_celltypes)
        lambda_prior: Strength of prior (default: 1.0)
        min_expression_threshold: Minimum expression to consider "active" (default: 0.1)
    
    Returns:
        Dict containing:
            - global_prior: Prior matrix (T_celltypes × M_genes)
            - confidence_scores: Confidence in each prior value
            - expression_patterns: Summary of expression patterns
    """
    # Validate inputs
    N = cell_type_numbers_array.shape[0]
    T = cell_type_numbers_array.shape[1]
    
    spot_keys = sorted(spotwise_gene_expression_profiles.keys())
    if len(spot_keys) != N:
        raise ValueError(f"Mismatch in number of spots: {len(spot_keys)} vs {N}")
    
    # Get dimensions from first profile
    example_profile = spotwise_gene_expression_profiles[spot_keys[0]]
    M = example_profile.shape[1]  # number of genes
    
    # Initialize arrays
    usage_array = np.zeros((N, T, M))
    for i, profile in spotwise_gene_expression_profiles.items():
        usage_array[i] = profile
    
    # Calculate expression statistics per cell type
    mean_expression = np.zeros((T, M))
    expression_frequency = np.zeros((T, M))
    expression_consistency = np.zeros((T, M))
    
    for t in range(T):
        # Weight profiles by cell type abundance
        weights = cell_type_numbers_array[:, t]  # Now 1D array of shape (N,)
        
        # Calculate weighted statistics
        active_expression = usage_array[:, t, :] > min_expression_threshold
        weighted_expression = usage_array[:, t, :]  # Shape: (N, M)
        
        # Mean expression when the cell type is present
        present_mask = weights > 0
        if np.any(present_mask):
            # Ensure weights match the data shape for averaging
            weights_for_average = weights[present_mask]  # 1D array of length n_present
            expression_for_average = weighted_expression[present_mask, :]  # (n_present, M)
            
            mean_expression[t] = np.average(
                expression_for_average,
                weights=weights_for_average,
                axis=0
            )
        
            # Expression consistency (coefficient of variation, inverse)
            # Calculate weighted std dev properly
            diff_squared = (expression_for_average - mean_expression[t]) ** 2  # (n_present, M)
            weighted_var = np.average(diff_squared, weights=weights_for_average, axis=0)  # (M,)
            std = np.sqrt(weighted_var)  # (M,)
            expression_consistency[t] = 1 / (1 + std / (mean_expression[t] + 1e-6))
        
        # Frequency of expression (properly weighted)
        total_weight = np.sum(weights) + 1e-6
        expression_frequency[t] = np.sum(active_expression * weights[:, np.newaxis], axis=0) / total_weight
    
    # Combine metrics into confidence scores
    confidence_scores = expression_frequency * expression_consistency
    
    # Generate prior probabilities
    # Scale mean expression to [0,1] per gene
    scaled_expression = mean_expression / (np.max(mean_expression, axis=0) + 1e-6)
    
    # Weight by confidence and apply prior strength
    weighted_scores = scaled_expression * np.power(confidence_scores, lambda_prior)
    
    # Convert to probabilities via softmax
    global_prior = np.zeros((T, M))
    for m in range(M):
        scores = weighted_scores[:, m]
        exp_scores = np.exp(scores - np.max(scores))  # Numerical stability
        global_prior[:, m] = exp_scores / (np.sum(exp_scores) + 1e-6)
    
    # Log statistics
    logging.info("Prior computation statistics:")
    logging.info(f" - Mean confidence score: {np.mean(confidence_scores):.4f}")
    logging.info(f" - Mean prior strength: {np.mean(global_prior):.4f}")
    logging.info(f" - % Strong signals (>0.5): {100 * np.mean(global_prior > 0.5):.2f}%")
    
    return {
        'global_prior': global_prior,
        'confidence_scores': confidence_scores,
        'expression_patterns': {
            'mean_expression': mean_expression,
            'expression_frequency': expression_frequency,
            'expression_consistency': expression_consistency
        }
    }

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
        logging.info("Adata variables: %s", adata.var_names)
        logging.info("Antibody markers: %s", all_markers)
        raise ValueError("No matching antibody markers found in adata.var_names.")
    
    adata.var_names_make_unique()
    adata = adata[:, existing_markers]

    # Step 2: Extract and prepare antibody capture data
    antibody_capture_data = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X.X
    antibody_capture_var_names = np.array(adata.var_names)

    cell_type_names = list(cell_profile_dict.keys())
    N = antibody_capture_data.shape[0]
    T = len(cell_type_names)

    profile_based_antibody_data = np.zeros((N, T))

    # Step 3: Map antibodies to profiles
    for profile_idx, (profile_name, profile_markers) in enumerate(cell_profile_dict.items()):
        major_markers = profile_markers.get("Major", [])
        try:
            relevant_marker_indices = [
                np.where(antibody_capture_var_names == marker)[0][0]
                for marker in major_markers if marker in antibody_capture_var_names
            ]
            if relevant_marker_indices:
                profile_based_antibody_data[:, profile_idx] = antibody_capture_data[:, relevant_marker_indices].mean(axis=1)
            else:
                logging.warning(f"No valid markers found for profile '{profile_name}'")
        except IndexError as e:
            logging.warning(f"Error processing markers for profile '{profile_name}': {str(e)}")

    # Step 4: Normalize with safety checks
    column_max = np.max(profile_based_antibody_data, axis=0)
    zero_columns = column_max == 0
    if np.any(zero_columns):
        logging.warning("Zero columns detected. Adding epsilon to prevent NaNs.")
        column_max[zero_columns] = 1e-6
    
    profile_based_antibody_data /= column_max

    if np.isnan(profile_based_antibody_data).any():
        raise ValueError("NaN values detected in profile_based_antibody_data after mapping.")

    return profile_based_antibody_data, cell_type_names

def optimize_cell_proportions(
    profile_based_antibody_data: np.ndarray,
    cell_type_names: List[str],
    tolerance: float = 1e-4,
    max_iterations: int = 50,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    normalize_beta: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform EM-based optimization for cell type proportions using Gurobi.

    Args:
        profile_based_antibody_data: N x T matrix of mapped antibody data
        cell_type_names: List of cell type names
        tolerance: Convergence tolerance for EM algorithm
        max_iterations: Maximum number of iterations
        lambda_reg: Regularization strength
        alpha: L1-L2 tradeoff factor (0 = L2, 1 = L1)
        normalize_beta: Whether to normalize beta values

    Returns:
        Tuple[np.ndarray, np.ndarray]: Y_values (N x T), beta_values (T,)
    """
    import gurobipy as gp
    from gurobipy import GRB

    N, T = profile_based_antibody_data.shape

    # Initialize beta estimates
    beta_estimates = {ct: 1.0 for ct in cell_type_names}
    beta_prev = np.zeros(T)
    Y_prev = np.zeros((N, T))
    iteration = 0

    while iteration < max_iterations:
        logging.info(f"\nIteration {iteration + 1}")
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

        # Sum of proportions constraints
        for i in range(N):
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) >= 0.9)
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) <= 1.2)

        try:
            model.optimize()
        except Exception as e:
            logging.error(f"Optimization error: {str(e)}")
            raise ValueError("Gurobi optimization failed") from e

        if model.status == GRB.OPTIMAL:
            Y_values = np.array([[Y[i, j].X for j in range(T)] for i in range(N)])
        else:
            raise ValueError("Gurobi optimization failed to converge")

        # Update beta
        beta_new = np.zeros(T)
        for j in range(T):
            Y_j = Y_values[:, j]
            S_j = profile_based_antibody_data[:, j]
            denominator = np.dot(Y_j, Y_j)
            
            if denominator > 0:
                beta_new[j] = np.dot(S_j, Y_j) / denominator
            beta_new[j] = max(beta_new[j], 0.0)  # Ensure non-negative

        # Optionally normalize beta values
        if normalize_beta:
            max_beta = np.max(beta_new)
            if max_beta > 0:
                beta_new = beta_new / max_beta

        # Convergence checks
        beta_diff = np.linalg.norm(beta_new - beta_prev)
        Y_diff = np.linalg.norm(Y_values - Y_prev)

        logging.info(f"Change in beta: {beta_diff:.6f}, Change in Y: {Y_diff:.6f}")
        if beta_diff < tolerance and Y_diff < tolerance:
            logging.info("Convergence achieved.")
            break

        # Update estimates for next iteration
        for j, ct_name in enumerate(cell_type_names):
            beta_estimates[ct_name] = beta_new[j]

        # Assert that beta_new is within the range [0, 1]
        assert np.all(beta_new >= 0) and np.all(beta_new <= 1), "Beta values must be within the range [0, 1]"

        beta_prev = beta_new.copy()
        Y_prev = Y_values.copy()
        iteration += 1

    return Y_values, beta_new



def finetune_cell_proportions(
    profile_based_antibody_data: np.ndarray,
    cell_type_names: List[str],
    initial_Y_values: np.ndarray,
    initial_beta_values: np.ndarray,
    adata: sc.AnnData,
    radius: float = 4.0,
    tolerance: float = 1e-4,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    max_iterations: int = 20,
    max_y_change: float = 0.4,
    beta_vary: bool = True,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 100,
    output_dir: str = "checkpoints",

    rerun: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Refine cell proportions using local neighborhood optimization with parallelization.

    Args:
        profile_based_antibody_data (np.ndarray):
            (N x T) array of mapped antibody intensities.
        cell_type_names (List[str]):
            Ordered list of length T specifying cell type names.
        initial_Y_values (np.ndarray):
            Initial cell proportion matrix of shape (N, T).
        initial_beta_values (np.ndarray):
            Initial beta estimates (shape T,). This is passed for consistency,
            but local solver decides how to use or ignore it based on beta_vary.
        adata (sc.AnnData):
            AnnData object with spot-level spatial coordinates in obsm['spatial'].
        radius (float):
            Radius for neighborhood-based local refinement.
        tolerance (float):
            Convergence tolerance for the local optimization loops.
        lambda_reg (float):
            Elastic net regularization strength for local solver.
        alpha (float):
            L1-L2 tradeoff (0 = purely L2, 1 = purely L1) in local solver.
        max_iterations (int):
            Maximum number of iterations for local solver.
        max_y_change (float):
            Maximum allowed change in Y values between iterations (default: 0.2).
            Values are constrained to vary by at most this amount while staying in [0,1].
        beta_vary (bool):
            If True, each spot's local solver is allowed to update betas;
            if False, betas remain fixed at the values passed in initial_beta_values.
        max_workers (int, optional):
            Maximum number of parallel workers. If None, uses os.cpu_count().
        checkpoint_interval (int):
            Number of spots between checkpoints.
        output_dir (str):
            Directory for checkpoints.
        rerun (bool):
            Whether to rerun if results exist.

    Returns:
        Tuple[np.ndarray, np.ndarray]:
            - A new (N x T) array of refined Y-values, obtained from a single pass
              of local refinements.
            - The original beta_values array, returned for interface consistency.
    """
    import os
    import gc
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool
    from tqdm import tqdm

    if initial_Y_values.ndim != 2:
        raise ValueError("initial_Y_values must be a 2D array (N x T).")
    if initial_beta_values.ndim != 1:
        raise ValueError("initial_beta_values must be a 1D array of length T.")

    N, T = profile_based_antibody_data.shape
    if initial_Y_values.shape != (N, T):
        raise ValueError("Mismatch between profile_based_antibody_data and initial_Y_values shapes.")
    if len(cell_type_names) != T:
        raise ValueError("cell_type_names length must match the number of columns in profile_based_antibody_data.")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Make a copy of the initial Y-values to store final refinements
    Y_refined = initial_Y_values.copy()

    # Calculate number of workers
    workers = max_workers if max_workers is not None else os.cpu_count()
    logging.info(f"Using {workers} workers for cell proportion refinement")

    # Process all spots in parallel
    futures = {}
    retry_count = 0
    max_retries = 3

    logging.info("Starting local cell proportion refinement")
    logging.info("Lambda reg: %s", lambda_reg)
    logging.info("Alpha: %s", alpha)

    while retry_count < max_retries:
        try:
            with ProcessPoolExecutor(max_workers=workers) as executor:
                futures.clear()
                for spot_idx in range(N):
                    future = executor.submit(
                        deconvolute_local_cell_proportions,
                        spot_idx=spot_idx,
                        adata=adata,
                        profile_based_antibody_data=profile_based_antibody_data,
                        radius=radius,
                        tolerance=tolerance,
                        lambda_reg=lambda_reg,
                        alpha=alpha,
                        beta_values=initial_beta_values,
                        beta_vary=beta_vary,
                        max_iterations=max_iterations,
                        max_y_change=max_y_change
                    )
                    futures[future] = spot_idx

                spots_processed = 0
                with tqdm(total=N, desc="Refining Cell Proportions") as pbar:
                    for future in as_completed(futures):
                        spot_idx = futures[future]
                        try:
                            result = future.result(timeout=300)
                            if result is not None:
                                Y_refined[spot_idx, :] = result
                                spots_processed += 1
                                pbar.update(1)

                                if spots_processed % checkpoint_interval == 0:
                                    # Save checkpoint
                                    checkpoint_path = os.path.join(output_dir, f"cell_prop_refinement_checkpoint_{spots_processed}.npy")
                                    np.save(checkpoint_path, Y_refined)
                                    logging.info(f"Saved checkpoint after {spots_processed} spots")

                        except TimeoutError:
                            logging.error(f"Timeout processing spot {spot_idx}")
                            continue
                        except Exception as e:
                            logging.error(f"Error processing spot {spot_idx}: {str(e)}")
                            continue

                break

        except BrokenProcessPool:
            retry_count += 1
            logging.warning(f"Process pool broken, retry {retry_count}/{max_retries}")
            if retry_count == max_retries:
                logging.error("Max retries reached, saving current progress")
            import time
            time.sleep(5)

    # Cleanup
    if futures:
        futures.clear()
    gc.collect()

    # Save final results
    final_path = os.path.join(output_dir, "cell_prop_refinement_final.npy")
    np.save(final_path, Y_refined)
    logging.info("Saved final refined cell proportions")

    return Y_refined, initial_beta_values



def deconvolute_local_cell_proportions(
    spot_idx: int,
    adata: sc.AnnData,
    profile_based_antibody_data: np.ndarray,
    radius: float = 2.0,
    tolerance: float = 1e-4,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    beta_values: Optional[np.ndarray] = None,
    beta_vary: bool = True,
    normalize_beta: bool = True,
    max_iterations: int = 20,
    max_y_change: float = 0.4
) -> Optional[np.ndarray]:
    """
    Refine cell proportions for a single spot via local neighborhood optimization.

    Args:
        spot_idx (int):
            Index of the spot to refine in the AnnData object.
        adata (sc.AnnData):
            AnnData containing spot-level spatial coordinates in obsm['spatial'].
        profile_based_antibody_data (np.ndarray):
            (N x T) global antibody intensities for N spots, T cell types.
        radius (float):
            Neighborhood radius for identifying neighbors.
        tolerance (float):
            Convergence threshold for Y- and beta-updates (if beta_vary=True).
        lambda_reg (float):
            Strength of elastic net regularization.
        alpha (float):
            L1-L2 tradeoff for the elastic net (0 = L2, 1 = L1).
        beta_values (Optional[np.ndarray]):
            Global or initial local beta values (length T). If None and beta_vary=True,
            local betas initialize at 1.0 each.
        beta_vary (bool):
            If True, local betas are iteratively updated.
            If False, beta_values remain fixed throughout optimization.
        normalize_beta (bool):
            Whether to normalize beta values after updates.
        max_iterations (int):
            Maximum iterations allowed for EM-like steps within this local function.
        max_y_change (float):
            Maximum allowed change in Y values between iterations (default: 0.2).
            Values are constrained to vary by at most this amount while staying in [0,1].

    Returns:
        Optional[np.ndarray]:
            Refined proportions (T,) for the specified spot, or None on failure.
    """
    import gurobipy as gp
    from gurobipy import GRB

    # Identify indices of spot's local neighborhood
    neighbor_indices = get_neighbors_with_fixed_radius(spot_idx, adata, radius=int(radius), include_center=True)
    if not neighbor_indices:
        logging.error(f"[Local Cell Props] No valid neighbors for spot {spot_idx}.")
        return None
    neighbor_indices = np.array(neighbor_indices, dtype=int)

    local_antibody_data = profile_based_antibody_data[neighbor_indices, :]
    local_N, T = local_antibody_data.shape

    if local_N == 0:
        logging.error(f"[Local Cell Props] Spot {spot_idx} has empty local antibody data.")
        return None

    # Identify center spot's position in neighbor list
    try:
        center_local_idx = np.where(neighbor_indices == spot_idx)[0][0]
    except IndexError:
        logging.error(f"[Local Cell Props] Could not find spot {spot_idx} in neighbor list.")
        return None

    # Initialize local betas
    if beta_values is not None and len(beta_values) == T:
        local_beta = beta_values.copy()
    else:
        local_beta = np.ones(T, dtype=float)

    beta_prev = local_beta.copy()

    # Initialize local Y to something uniform
    Y_prev = np.full((local_N, T), 1.0 / T)

    iteration = 0
    while iteration < max_iterations:
        try:
            model = gp.Model(f"Local_Cell_Props_spot_{spot_idx}")
            model.setParam('OutputFlag', 0)
            model.setParam('TimeLimit', 60)
            model.setParam('MIPGap', 0.01)

            # Build Y variables in [0, 1]
            Y_vars = model.addVars(local_N, T, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="Y")

            # Summation constraints on each row
            for i in range(local_N):
                model.addConstr(gp.quicksum(Y_vars[i, j] for j in range(T)) >= 0.9)
                model.addConstr(gp.quicksum(Y_vars[i, j] for j in range(T)) <= 1.2)

            # Add constraints to limit Y value changes from previous iteration
            if iteration > 0:
                for i in range(local_N):
                    for j in range(T):
                        prev_value = Y_prev[i, j]
                        # Lower bound: max(0, prev_value - max_y_change)
                        # Upper bound: min(1, prev_value + max_y_change)
                        lb = max(0.0, prev_value - max_y_change)
                        ub = min(1.0, prev_value + max_y_change)
                        model.addConstr(Y_vars[i, j] >= lb)
                        model.addConstr(Y_vars[i, j] <= ub)

            # Objective: sum of squared differences + elastic net
            error_terms = []
            for i in range(local_N):
                for j in range(T):
                    S_ij = local_antibody_data[i, j]
                    error_terms.append((S_ij - local_beta[j] * Y_vars[i, j]) ** 2)

            total_error = gp.quicksum(error_terms)
            l1 = gp.quicksum(Y_vars[i, j] for i in range(local_N) for j in range(T))
            l2 = gp.quicksum(Y_vars[i, j] * Y_vars[i, j] for i in range(local_N) for j in range(T))
            reg_term = lambda_reg * (alpha * l1 + (1.0 - alpha) * l2)
            model.setObjective(total_error + reg_term, GRB.MINIMIZE)

            model.optimize()

            if model.status != GRB.OPTIMAL:
                logging.warning(f"[Local Cell Props] Spot {spot_idx} local optimization not optimal (status: {model.status}).")
                return None

            # Extract current Y solution
            Y_values = np.array([[Y_vars[i, j].X for j in range(T)] for i in range(local_N)])

            # Update local beta if allowed
            if beta_vary:
                new_beta = np.zeros(T, dtype=float)
                for j in range(T):
                    Y_j = Y_values[:, j]
                    S_j = local_antibody_data[:, j]
                    denominator = np.dot(Y_j, Y_j)
                    
                    if denominator > 1e-15:
                        new_beta[j] = np.dot(S_j, Y_j) / denominator
                    new_beta[j] = max(new_beta[j], 0.0)  # Ensure non-negative

                # Optionally normalize beta values
                if normalize_beta:
                    max_beta = np.max(new_beta)
                    if max_beta > 0:
                        new_beta = new_beta / max_beta
            else:
                new_beta = local_beta.copy()

            # Check convergence
            beta_diff = np.linalg.norm(new_beta - beta_prev) if beta_vary else 0.0
            Y_diff = np.linalg.norm(Y_values - Y_prev)

            logging.debug(f"Spot {spot_idx} - Iteration {iteration + 1}: "
                        f"beta_diff={beta_diff:.6f}, Y_diff={Y_diff:.6f}")

            if beta_diff < tolerance and Y_diff < tolerance:
                logging.debug(f"Spot {spot_idx} converged after {iteration + 1} iterations")
                Y_prev = Y_values
                local_beta = new_beta
                break

            # Prepare for next iteration
            Y_prev = Y_values.copy()
            local_beta = new_beta.copy()
            beta_prev = new_beta.copy()
            iteration += 1

        except Exception as e:
            logging.error(f"Error in local optimization for spot {spot_idx}: {str(e)}")
            return None

        finally:
            if 'model' in locals():
                del model
            gc.collect()

    # Return just the center row of Y for this spot
    return Y_prev[center_local_idx, :]


################################################################################
# === DECONVOLUTION FOR GENES ===
################################################################################
def deconvolute_spot_with_neighbors_with_prior(
    spot_idx: int,
    adata: sc.AnnData,
    cell_type_numbers_array: np.ndarray,
    radius: float,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    local_enrichment_weight: float = 0.5,
    global_enrichment_weight: float = 0.5,
) -> Optional[np.ndarray]:
    """
    Deconvolute a spot with its neighbors, using both enrichment weights and optional prior.
    """
    model = None
    try:
        neighborhood_indices = get_neighbors_with_fixed_radius(
            spot_idx, adata, radius=int(radius), include_center=True
        )
        if not neighborhood_indices:
            logging.error(f"No valid neighbors found for spot {spot_idx}.")
            return None

        neighborhood_indices = np.array([
            int(idx) for idx in neighborhood_indices
            if isinstance(idx, (int, np.integer))
        ], dtype=int)

        # Extract expression data
        deconvolution_expression_data = adata.X
        if scipy.sparse.issparse(deconvolution_expression_data):
            deconvolution_expression_data = deconvolution_expression_data.toarray()  # type: ignore
        elif not isinstance(deconvolution_expression_data, np.ndarray):
            deconvolution_expression_data = np.array(deconvolution_expression_data)

        # Dimensions
        T = cell_type_numbers_array.shape[1]  # number of cell types
        M = deconvolution_expression_data.shape[1]  # number of genes

        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]

        # Compute normalized cell type weights to avoid abundance bias
        total_celltype_counts = np.sum(cell_type_numbers_array, axis=0) + 1e-10
        celltype_frequencies = total_celltype_counts / np.sum(total_celltype_counts)
        inverse_frequency_weights = 1.0 / (celltype_frequencies + 1e-10)
        normalized_weights = inverse_frequency_weights / np.max(inverse_frequency_weights)

        # Modified enrichment calculation
        def compute_expression_aware_enrichment(expression_data, cell_type_props, gene_idx):
            """
            Compute expression-aware enrichment scores.
            
            Args:
                expression_data (np.ndarray): Expression matrix
                cell_type_props (np.ndarray): Cell type proportions
                gene_idx (int): Gene index
                
            Returns:
                np.ndarray: Enrichment scores for each cell type
            """
            gene_expr = expression_data[:, gene_idx]
            expr_threshold = np.percentile(gene_expr[gene_expr > 0], 50) if np.any(gene_expr > 0) else 0
            high_expr_spots = gene_expr >= expr_threshold

            if not np.any(high_expr_spots):
                return np.ones(cell_type_props.shape[1]) / cell_type_props.shape[1]

            # Normalize cell type proportions by their global frequency
            normalized_props = cell_type_props / (celltype_frequencies + 1e-10)
            
            high_expr_props = np.mean(normalized_props[high_expr_spots], axis=0)
            background_props = np.mean(normalized_props, axis=0)

            epsilon = 1e-10
            enrichment = high_expr_props / (background_props + epsilon)
            
            # Apply smoothing to avoid extreme values
            smoothed_enrichment = 0.8 * enrichment + 0.2 * np.ones_like(enrichment)
            return smoothed_enrichment / (np.sum(smoothed_enrichment) + epsilon)

        # Compute expression-aware enrichment for each gene
        gene_specific_enrichment = np.zeros((M, T))

        for k in range(M):
            local_enrich = compute_expression_aware_enrichment(
                neighborhood_expression_data,
                neighborhood_cell_type_numbers,
                k
            )
            global_enrich = compute_expression_aware_enrichment(
                deconvolution_expression_data,
                cell_type_numbers_array,
                k
            )
            gene_specific_enrichment[k] = (
                local_enrichment_weight * local_enrich +
                global_enrichment_weight * global_enrich
            )

        # Build Gurobi model
        model = gp.Model(f"discrete_gene_expression_spot_{spot_idx}")
        model.setParam('OutputFlag', 0)
        model.setParam('Threads', 1)
        model.setParam('NodefileStart', 0.5)
        model.setParam('MIPGap', 0.01)
        model.setParam('TimeLimit', 600)
        model.setParam('NodeLimit', 1000000)

        # Variables for count assignment
        X = {}
        center_counts = deconvolution_expression_data[spot_idx, :]

        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    X[j, k] = model.addVar(
                        vtype=GRB.INTEGER,
                        lb=0,
                        ub=total_counts,
                        name=f"X_{j}_{k}"
                    )
                # Count conservation constraint
                model.addConstr(
                    gp.quicksum(X[j, k] for j in range(T)) == total_counts,
                    name=f"count_conservation_{k}"
                )

        # Validate prior if asked
        if global_prior is not None:
            if lambda_prior_weight > 0:
                if global_prior is None:
                    raise ValueError("lambda_prior_weight > 0 but no global_prior provided")
            if not isinstance(global_prior, np.ndarray):
                raise ValueError("global_prior must be a numpy array")
            if global_prior.shape != (T, M):
                raise ValueError(f"Prior matrix shape {global_prior.shape} does not match expected shape ({T}, {M})")

        # Modify objective terms to include frequency normalization
        obj_terms = []
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    # Get normalized weights
                    enrichment_weight = gene_specific_enrichment[k, j]
                    cell_type_weight = neighborhood_cell_type_numbers[len(neighborhood_indices) // 2, j]
                    
                    # Apply frequency normalization
                    normalized_weight = cell_type_weight * normalized_weights[j]
                    
                    # Add slight randomness to break ties using seeded RNG for reproducibility
                    rng = np.random.default_rng(42)  # Fixed seed for reproducibility
                    randomness = 0.9 + 0.2 * rng.random()
                    base_term = enrichment_weight * normalized_weight * randomness * X[j,k]
                    obj_terms.append(base_term)

                    # Prior terms remain unchanged
                    if global_prior is not None and lambda_prior_weight > 0:
                        try:
                            prior_value = float(global_prior[j, k])
                            prior_penalty = lambda_prior_weight * (1 - prior_value) * X[j, k]
                            obj_terms.append(-prior_penalty)
                        except Exception as e:
                            logging.warning(f"Error accessing prior at [{j}, {k}]: {str(e)}")
                            continue

        # Maximize the sum of all terms
        model.setObjective(
            gp.quicksum(obj_terms),
            GRB.MAXIMIZE
        )

        model.optimize()

        if model.status == GRB.OPTIMAL:
            logging.info(f"Solution found for spot {spot_idx}")
            result = np.zeros((T, M))
            for k in range(M):
                total_counts = int(center_counts[k])
                if total_counts > 0:
                    for j in range(T):
                        result[j, k] = X[j, k].X
            return result
        else:
            logging.error(f"No feasible solution found for spot {spot_idx}.")
            return None

    except Exception as e:
        logging.error(f"Error during deconvolution of spot {spot_idx}: {str(e)}")
        logging.error(traceback.format_exc())
        return None

    finally:
        if model:
            del model
        gc.collect()

def log_marker_gene_patterns(zero_patterns, marker_genes):
    """
    Log detailed patterns for marker genes.
    """
    for gene in marker_genes:
        logging.info(f"\nPatterns for {gene}:")
        for ct, genes_data in zero_patterns.items():
            if gene in genes_data:
                stats = genes_data[gene]
                logging.info(f"  {ct}:")
                logging.info(f"    Zero proportion: {stats['zero_proportion']:.3f}")
                if stats['n_nonzero'] > 0:
                    logging.info(f"    Mean nonzero expression: {stats['mean_nonzero_expression']:.3f}")
                else:
                    logging.info(f"    Mean nonzero expression: 0.0 (no nonzero values)")
                logging.info(f"    Number of spots: {stats['n_spots']}")
                logging.info(f"    Number of nonzero spots: {stats['n_nonzero']}")

def scale_genes(expression_matrix):
    """Scale each gene independently to [0,1] range.
    
    Args:
        expression_matrix (np.ndarray): Spots x Genes matrix
        
    Returns:
        tuple: (scaled_matrix, gene_mins, gene_maxs)
    """
    gene_mins = np.min(expression_matrix, axis=0)
    gene_maxs = np.max(expression_matrix, axis=0)
    
    # Avoid division by zero
    gene_ranges = np.maximum(gene_maxs - gene_mins, 1e-10)
    scaled_matrix = (expression_matrix - gene_mins) / gene_ranges
    
    return scaled_matrix, gene_mins, gene_maxs

def unscale_genes(scaled_matrix, gene_mins, gene_maxs):
    """Reverse the gene-wise scaling transformation.
    
    Args:
        scaled_matrix (np.ndarray): Scaled matrix
        gene_mins (np.ndarray): Original minimum values per gene
        gene_maxs (np.ndarray): Original maximum values per gene
        
    Returns:
        np.ndarray: Unscaled matrix
    """
    gene_ranges = np.maximum(gene_maxs - gene_mins, 1e-10)
    return (scaled_matrix * gene_ranges) + gene_mins

def optimize_gene_expression(
    sample_name: str,
    deconvolution_expression_data: np.ndarray,
    cell_type_numbers_array: np.ndarray,
    filtered_adata: sc.AnnData,
    radius: float = 2,
    global_enrichment_weight: float = 0.5,
    local_enrichment_weight: float = 0.5,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 100,
    output_dir: str = "checkpoints",
    rerun: bool = False
) -> Dict[str, Any]:
    """
    Optimize gene expression with enrichment weights and prior guidance.
    
    Args:
        sample_name (str): Name of the sample
        deconvolution_expression_data (np.ndarray): Gene expression data (N_spots x M_genes)
        cell_type_numbers_array (np.ndarray): Cell type proportions (N_spots x T_celltypes)
        filtered_adata (sc.AnnData): Filtered AnnData object containing gene expression data
        radius (float): Radius for neighbor detection
        alpha (float): Weight for spatial regularization
        lambda_reg_gex (float): Weight for gene expression regularization
        global_enrichment_weight (float): Weight for global expression enrichment (0-1)
        local_enrichment_weight (float): Weight for local expression enrichment (0-1)
        global_prior (np.ndarray, optional): Global prior matrix for guidance
        lambda_prior_weight (float): Weight for prior guidance
        max_workers (int, optional): Maximum number of parallel workers
        checkpoint_interval (int): Number of spots between checkpoints
        output_dir (str): Directory for checkpoints
        rerun (bool): Whether to rerun if results exist
        
    Returns:
        Dict[str, Any]: {
            'spotwise_profiles': Dict[int, np.ndarray],
            'dimensions': Tuple[int, int, int]
        }
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    N = deconvolution_expression_data.shape[0]
    M = deconvolution_expression_data.shape[1]
    T = cell_type_numbers_array.shape[1]

    # Initialize checkpoint manager
    checkpoint_mgr = CheckpointManager(output_dir, sample_name)
    
    # If rerun=False, check for completed run
    if not rerun:
        complete_results = checkpoint_mgr.check_complete_run(N, T, M)
        if complete_results is not None:
            return complete_results  # type: ignore
            
        # Load latest checkpoint if available
        completed_spots, spotwise_gene_expression_profiles = checkpoint_mgr.load_latest_checkpoint(N, T, M)
    else:
        completed_spots = set()
        spotwise_gene_expression_profiles = {}

    logging.info(f"Starting analysis for {sample_name}")
    logging.info(f"Already completed spots: {len(completed_spots)}")
    
    # Log configuration
    if global_enrichment_weight + local_enrichment_weight > 0:
        logging.info(f"Using enrichment weights - Global: {global_enrichment_weight}, Local: {local_enrichment_weight}")
    if global_prior is not None and lambda_prior_weight > 0:
        logging.info("Using prior-guided deconvolution")

    # Initialize futures as empty dict before try block
    futures = {}
    
    try:
        # Calculate number of workers (ensure it's an integer)
        workers = max_workers if max_workers is not None else os.cpu_count()
        logging.info(f"Using {workers} workers")
        
        # Only process spots that haven't been completed
        remaining_spots = [i for i in range(N) if i not in completed_spots]
        logging.info(f"Processing {len(remaining_spots)} remaining spots")
        
        retry_count = 0
        max_retries = 3
        while retry_count < max_retries:
            try:
                with ProcessPoolExecutor(max_workers=workers) as executor:
                    futures.clear()
                    for spot_idx in remaining_spots:
                        # Always use the same function with consistent args
                        future = executor.submit(
                            deconvolute_spot_with_neighbors_with_prior,
                            spot_idx,
                            filtered_adata,
                            cell_type_numbers_array,
                            radius,
                            global_prior,
                            lambda_prior_weight,
                            local_enrichment_weight,
                            global_enrichment_weight,
                        )
                        futures[future] = spot_idx

                    with tqdm(total=len(remaining_spots), desc="Deconvoluting Remaining Spots") as pbar:
                        spots_since_last_save = 0
                        
                        for future in as_completed(futures):
                            i = futures[future]
                            try:
                                result = future.result(timeout=300)
                                if result is not None and result.ndim == 2:
                                    spotwise_gene_expression_profiles[i] = result.copy()
                                    completed_spots.add(i)
                                    spots_since_last_save += 1
                                    pbar.update(1)

                                    if spots_since_last_save >= checkpoint_interval:
                                        checkpoint_mgr.save_checkpoint(
                                            completed_spots,
                                            spotwise_gene_expression_profiles,
                                            N, T, M
                                        )
                                        spots_since_last_save = 0
                            except TimeoutError:
                                logging.error(f"Timeout processing spot {i}")
                                continue
                            except Exception as e:
                                logging.error(f"Error processing spot {i}: {str(e)}")
                                logging.error(traceback.format_exc())
                                continue
                
                break
                
            except concurrent.futures.process.BrokenProcessPool:  # type: ignore
                retry_count += 1
                logging.warning(f"Process pool broken, retry {retry_count}/{max_retries}")
                if retry_count == max_retries:
                    logging.error("Max retries reached, saving current progress")
                import time
                time.sleep(5)

    finally:
        if futures:
            futures.clear()
        gc.collect()
        
        if spotwise_gene_expression_profiles:
            checkpoint_mgr.save_final_results(
                spotwise_gene_expression_profiles,
                completed_spots,
                N, T, M
            )

    return spotwise_gene_expression_profiles


def normalize_counts(adata, target_sum=10000, exclude_highly_expressed=False, max_fraction=0.05):
    """
    Normalize counts while preserving integer values and relative proportions.
    
    Args:
        adata: AnnData object
        target_sum: Target sum for each cell/spot
        exclude_highly_expressed: Whether to exclude highly expressed genes
        max_fraction: Maximum fraction for highly expressed genes
    """
    # Get matrix
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X.copy()
    
    # Handle highly expressed genes if requested
    if exclude_highly_expressed:
        counts_per_cell = X.sum(axis=1)
        gene_subset = ~(X > counts_per_cell[:, None] * max_fraction).any(axis=0)
        size_factors = X[:, gene_subset].sum(axis=1)
    else:
        size_factors = X.sum(axis=1)
    
    # Ensure positive values
    size_factors = np.maximum(size_factors, 1)
    median_size = max(np.median(size_factors), 1)
    
    # Calculate bounded scaling factors
    scaling_factors = np.clip(size_factors / median_size, 0.1, 10.0)
    scaled_factors = (target_sum / size_factors)
    
    # Scale and round to integers
    X_scaled = np.round(X * scaled_factors[:, None]).astype(int)
    
    # Safety bounds
    max_allowed = target_sum * 2
    X_scaled = np.clip(X_scaled, 0, max_allowed)
    
    # Create new AnnData with scaled counts
    adata_norm = adata.copy()
    adata_norm.X = X_scaled
    
    # Store normalization info
    adata_norm.obs['size_factors'] = scaling_factors
    adata_norm.obs['original_total'] = size_factors
    adata_norm.obs['scaled_total'] = X_scaled.sum(axis=1)
    
    # Log statistics
    logging.info(f"Normalization stats:")
    logging.info(f"Original median total: {median_size:.2f}")
    logging.info(f"Mean scaled total: {X_scaled.sum(axis=1).mean():.2f}")
    logging.info(f"Max scaled value: {X_scaled.max():.2f}")
    
    return adata_norm

def validate_prior_effect(spotwise_profiles_pass1, spotwise_profiles_pass2, global_prior):
    """
    Compare pass1 and pass2 results to verify prior influence.
    
    Args:
        spotwise_profiles_pass1 (dict): First pass results {spot_idx: profile_matrix}
        spotwise_profiles_pass2 (dict): Second pass results {spot_idx: profile_matrix}
        global_prior (np.ndarray): Global prior matrix (T x M)
        
    Returns:
        dict: Dictionary containing validation metrics
    """
    # Validate shapes
    if not spotwise_profiles_pass1 or not spotwise_profiles_pass2:
        raise ValueError("Empty profile dictionaries provided")
        
    # Get shapes from first profile
    first_profile1 = next(iter(spotwise_profiles_pass1.values()))
    T, M = first_profile1.shape
    
    if global_prior.shape != (T, M):
        raise ValueError(f"Prior shape {global_prior.shape} does not match profiles shape ({T}, {M})")
    
    prior_guided_changes = []
    spot_metrics = {}
    
    # Ensure we have matching spots
    common_spots = set(spotwise_profiles_pass1.keys()) & set(spotwise_profiles_pass2.keys())
    
    if not common_spots:
        logging.error("No matching spots found between pass1 and pass2 results")
        return None
        
    for spot in common_spots:
        profile1 = spotwise_profiles_pass1[spot]
        profile2 = spotwise_profiles_pass2[spot]
        
        # Calculate absolute changes between passes
        profile_diff = np.abs(profile2 - profile1)
        total_diff = np.sum(profile_diff)
        
        # Calculate prior influence on pass2 assignment
        prior_alignment = np.sum(global_prior * profile2)
        
        prior_guided_changes.append((total_diff, prior_alignment))
        
        # Store per-spot metrics
        spot_metrics[spot] = {
            'total_change': total_diff,
            'prior_alignment': prior_alignment,
            'mean_change': np.mean(profile_diff),
            'max_change': np.max(profile_diff)
        }
    
    # Calculate correlation between changes and prior influence
    changes = np.array([x[0] for x in prior_guided_changes])
    influences = np.array([x[1] for x in prior_guided_changes])
    
    correlation = np.corrcoef(changes, influences)[0,1]
    
    # Calculate summary statistics
    validation_metrics = {
        'prior_correlation': correlation,
        'mean_total_change': np.mean(changes),
        'mean_prior_influence': np.mean(influences),
        'std_total_change': np.std(changes),
        'std_prior_influence': np.std(influences),
        'n_spots_analyzed': len(common_spots),
        'spot_metrics': spot_metrics
    }
    
    # Log summary statistics
    logging.info("Prior Effect Validation:")
    logging.info(f"Prior-Change Correlation: {correlation:.4f}")
    logging.info(f"Mean Total Change: {validation_metrics['mean_total_change']:.4f}")
    logging.info(f"Mean Prior Influence: {validation_metrics['mean_prior_influence']:.4f}")
    logging.info(f"Number of Spots Analyzed: {validation_metrics['n_spots_analyzed']}")
    
    return validation_metrics

