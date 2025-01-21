# Standard library imports
import os
import logging
import traceback
import gc
import concurrent
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Any, Optional, List, Tuple, Union
import json

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

def fit_zinb(x, p):
    """
    Fit ZINB model for a gene's expression (x) against cell type proportion (p)
    
    Args:
        x (np.ndarray): Gene expression values
        p (np.ndarray): Cell type proportions
    
    Returns:
        tuple: (π, μ, θ) - zero-inflation probability, mean, and dispersion
    """
    def zinb_nll(params):
        pi, mu, theta = params
        # Bound parameters to valid ranges more strictly
        pi = np.clip(pi, 0.001, 0.999)
        mu = np.clip(mu, 0.001, None)
        theta = np.clip(theta, 0.001, None)
        
        # Zero and non-zero indices
        zero_idx = x == 0
        nonzero_idx = ~zero_idx
        
        if not np.any(nonzero_idx):
            return 1e10  # Penalize all-zero cases
        
        # Log-likelihood calculation
        ll_zeros = np.sum(np.log(pi + (1-pi) * np.power(theta/(theta+mu), theta))[zero_idx])
        
        ll_nonzeros = np.sum(
            np.log(1-pi) + loggamma(x[nonzero_idx] + theta) - loggamma(theta) - 
            loggamma(x[nonzero_idx] + 1) + theta * np.log(theta) + 
            x[nonzero_idx] * np.log(mu) - (x[nonzero_idx] + theta) * np.log(theta + mu)
        )
        
        return -(ll_zeros + ll_nonzeros)
    
    # Better initial guesses
    zero_prop = np.mean(x == 0)
    pi_init = min(0.9, zero_prop)  # Cap initial pi
    mu_init = max(0.1, np.mean(x[x > 0])) if np.any(x > 0) else 0.1
    theta_init = 1.0
    
    try:
        # Optimize with stricter bounds on pi
        result = minimize(
            zinb_nll, 
            [pi_init, mu_init, theta_init],
            bounds=[(0.001, 0.95), (0.001, None), (0.001, None)],  # Cap maximum pi at 0.95
            method='L-BFGS-B'
        )
        
        # Validate results
        if not result.success:
            return (zero_prop, mu_init, theta_init)
        
        return result.x
    except:
        return (zero_prop, mu_init, theta_init)

def analyze_zero_inflation_patterns(
    usage_array,
    cell_type_numbers,
    gene_names,
    cell_type_names,
    min_proportion=0.01,
    expression_threshold=1e-6
):
    """Analyze zero-inflation patterns in gene expression data."""
    patterns = {}
    
    # Ensure we're working with numpy arrays
    if isinstance(usage_array, dict):
        # Convert dictionary to array with shape (T x M) for cell types x genes
        T = len(cell_type_names)    # number of cell types
        M = len(gene_names)         # number of genes
        
        # Take the first spot's profile as template for shape
        first_profile = next(iter(usage_array.values()))
        if first_profile.shape != (T, M):
            raise ValueError(f"Expected profile shape (T={T}, M={M}), got {first_profile.shape}")
            
        usage_array = first_profile  # Use first spot's profile as reference
    
    # Convert to dense if sparse
    if hasattr(usage_array, 'toarray'):
        usage_array = usage_array.toarray()
    
    if usage_array.ndim != 2:
        raise ValueError(f"Expected 2D array (cell_types x genes), got shape {usage_array.shape}")
    
    T, M = usage_array.shape  # T cell types, M genes
    N = cell_type_numbers.shape[0]  # number of spots
    
    for m_idx, gene in enumerate(gene_names):
        patterns[gene] = {}
        
        for t_idx, cell_type in enumerate(cell_type_names):
            # Get spots where this cell type is present
            ct_spots = cell_type_numbers[:, t_idx] > min_proportion
            
            if np.any(ct_spots):
                # Get expression value for this cell type and gene
                expr_value = usage_array[t_idx, m_idx]
                n_spots = np.sum(ct_spots)
                
                # For a cell type's gene expression, we only have one value
                zero_prop = 1.0 if expr_value <= expression_threshold else 0.0
                patterns[gene][cell_type] = {
                    'zero_proportion': zero_prop,
                    'mean_nonzero': expr_value if expr_value > expression_threshold else 0,
                    'n_spots': n_spots
                }
            else:
                patterns[gene][cell_type] = {
                    'zero_proportion': np.nan,
                    'mean_nonzero': 0,
                    'n_spots': 0
                }
    
    return patterns

def suggest_zero_inflation_threshold(patterns, quantile=0.75):
    """
    Suggest a zero-inflation threshold based on the analysis results.
    """
    all_zero_props = []
    for gene_patterns in patterns.values():
        for ct_data in gene_patterns.values():
            if not np.isnan(ct_data['zero_proportion']):
                all_zero_props.append(ct_data['zero_proportion'])
    
    if all_zero_props:
        return np.quantile(all_zero_props, quantile)
    return 0.5  # Default fallback

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
        logging.info("Adata variables: ", adata.var_names)
        logging.info("Antibody markers: ", all_markers)
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
            logging.warning(f"⚠️ No valid markers found for profile '{profile_name}'.")

    # Step 3: Normalize columns to prevent zero-division
    column_max = np.max(profile_based_antibody_data, axis=0)
    zero_columns = column_max == 0
    if np.any(zero_columns):
        logging.warning(f"⚠️ Warning: Zero columns detected. Adding epsilon to prevent NaNs.")
        column_max[zero_columns] = 1e-6

    profile_based_antibody_data /= column_max

    # Validate final data
    if np.isnan(profile_based_antibody_data).any():
        raise ValueError("NaN values detected in `profile_based_antibody_data` after mapping.")

    return profile_based_antibody_data, cell_type_names

def optimize_cell_proportions(profile_based_antibody_data, cell_type_names, tolerance=1e-4, max_iterations=50, lambda_reg=1, alpha=0.7):
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
        
        for i in range(N):
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) >= 0.9)
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) <= 1.2)
        
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
        
        logging.info(f"Change in beta: {beta_diff:.6f}, Change in Y: {Y_diff:.6f}")
        if beta_diff < tolerance and Y_diff < tolerance:
            logging.info("Convergence achieved.")
            break
        
        beta_prev = beta_new
        Y_prev = Y_values
        iteration += 1

    return Y_values

################################################################################
# === DECONVOLUTION FOR GENES ===
################################################################################
def deconvolute_spot_with_neighbors_with_prior(
    spot_idx: int,
    adata: sc.AnnData,
    cell_type_numbers_array: np.ndarray,
    radius: float,
    alpha: float = 0.5,
    lambda_reg_gex: float = 0.001,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    local_weight: float = 0.5,
    global_weight: float = 0.5
) -> Optional[np.ndarray]:
    """
    Deconvolute a spot with its neighbors, optionally using a prior.
    
    Args:
        spot_idx: Index of spot to deconvolute
        adata: AnnData object containing gene expression
        cell_type_numbers_array: Array of cell type numbers
        radius: Radius for neighbor search
        alpha: Weight for spatial term
        lambda_reg_gex: L1/L2 regularization weight
        global_prior: Optional prior matrix (T x M)
        lambda_prior_weight: Weight for prior guidance (0.0 means no prior)
        
    Returns:
        Optional[np.ndarray]: Deconvoluted profile matrix or None if error
    """
    model = None
    try:
        # Get neighborhood data first to establish dimensions
        neighborhood_indices = get_neighbors_with_fixed_radius(
            spot_idx, adata, radius=radius, include_center=True
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
        if hasattr(deconvolution_expression_data, 'toarray'):
            deconvolution_expression_data = deconvolution_expression_data.toarray()
        
        # Get dimensions from the data
        T = cell_type_numbers_array.shape[1]  # number of cell types
        M = deconvolution_expression_data.shape[1]  # number of genes
        
        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]

        # Compute expression-aware enrichment for each gene
        gene_specific_enrichment = np.zeros((M, T))

        def compute_expression_aware_enrichment(expression_data, cell_type_props, gene_idx):
            gene_expr = expression_data[:, gene_idx]
            expr_threshold = np.percentile(gene_expr[gene_expr > 0], 50) if np.any(gene_expr > 0) else 0
            high_expr_spots = gene_expr >= expr_threshold
            
            if not np.any(high_expr_spots):
                return np.ones(cell_type_props.shape[1]) / cell_type_props.shape[1]
            
            high_expr_props = np.mean(cell_type_props[high_expr_spots], axis=0)
            background_props = np.mean(cell_type_props, axis=0)
            
            epsilon = 1e-10
            enrichment = high_expr_props / (background_props + epsilon)
            smoothed_enrichment = 0.8 * enrichment + 0.2 * np.ones_like(enrichment)
            
            return smoothed_enrichment / (np.sum(smoothed_enrichment) + epsilon)

        # Compute enrichment scores
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
                local_weight * local_enrich +
                global_weight * global_enrich
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
                    X[j,k] = model.addVar(
                        vtype=GRB.INTEGER,
                        lb=0,
                        ub=total_counts,
                        name=f"X_{j}_{k}"
                    )

                # Count conservation constraint
                model.addConstr(
                    gp.quicksum(X[j,k] for j in range(T)) == total_counts,
                    name=f"count_conservation_{k}"
                )

        # Only validate prior if both prior and weight are provided

        if global_prior is not None:
            if lambda_prior_weight > 0:  # First check if we're using prior guidance at all
                if global_prior is None:
                    raise ValueError("lambda_prior_weight > 0 but no global_prior provided")
            if not isinstance(global_prior, np.ndarray):
                raise ValueError("global_prior must be a numpy array")
            if global_prior.shape != (T, M):
                raise ValueError(f"Prior matrix shape {global_prior.shape} does not match expected shape ({T}, {M})")

        # Objective terms
        obj_terms = []
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    # Base enrichment term
                    enrichment_weight = gene_specific_enrichment[k, j]
                    cell_type_weight = neighborhood_cell_type_numbers[
                        len(neighborhood_indices)//2, j
                    ]
                    randomness = 0.9 + 0.2 * np.random.random()
                    
                    base_term = enrichment_weight * cell_type_weight * randomness * X[j,k]
                    obj_terms.append(base_term)
                    
                    # Add prior-based penalty only if we have both prior and positive weight
                    if global_prior is not None and lambda_prior_weight > 0:
                        try:
                            prior_value = float(global_prior[j, k])
                            prior_penalty = lambda_prior_weight * (1 - prior_value) * X[j,k]
                            obj_terms.append(-prior_penalty)  # Subtract penalty
                        except Exception as e:
                            logging.warning(f"Error accessing prior at [{j}, {k}]: {str(e)}")
                            continue

        model.setObjective(
            gp.quicksum(obj_terms),
            GRB.MAXIMIZE
        )

        # Optimize
        model.optimize()

        if model.status == GRB.OPTIMAL:
            logging.info(f"Solution found for spot {spot_idx}")
            result = np.zeros((T, M))
            for k in range(M):
                total_counts = int(center_counts[k])
                if total_counts > 0:
                    for j in range(T):
                        result[j,k] = X[j,k].X
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
    alpha: float = 0.5,
    lambda_reg_gex: float = 0.001,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 100,
    output_dir: str = "checkpoints",
    rerun: bool = False
) -> Dict[str, Any]:
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
            return complete_results
            
        # Load latest checkpoint if available
        completed_spots, spotwise_gene_expression_profiles = checkpoint_mgr.load_latest_checkpoint(N, T, M)
    else:
        completed_spots = set()
        spotwise_gene_expression_profiles = {}

    logging.info(f"Starting analysis for {sample_name}")
    logging.info(f"Already completed spots: {len(completed_spots)}")
    
    # Log whether using prior-guided deconvolution
    if global_prior is not None and lambda_prior_weight > 0:
        logging.info("Using prior-guided deconvolution")
    else:
        logging.info("Using standard deconvolution")

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
                            alpha,
                            lambda_reg_gex,
                            global_prior,  # Will be None if not using prior
                            lambda_prior_weight  # Will be 0.0 if not using prior
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
                
            except concurrent.futures.process.BrokenProcessPool:
                retry_count += 1
                logging.warning(f"Process pool broken, retry {retry_count}/{max_retries}")
                if retry_count == max_retries:
                    logging.error("Max retries reached, saving current progress")
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

def approximate_wgcna(adata_gex, max_clusters=10, correlation_method='pearson'):
    """
    Approximate WGCNA-like approach in Python:
      1) Compute gene-gene correlations
      2) Convert to distance matrix
      3) Cluster using hierarchical clustering
      4) Define modules from clusters
      5) Build a (genes x modules) binary membership matrix

    Args:
        adata_gex (AnnData): Spots × Genes
        max_clusters (int): Maximum number of modules to cut from the tree
        correlation_method (str): 'pearson' or 'spearman'

    Returns:
        module_matrix (np.ndarray): shape=(G, K), where G=number of genes, K=actual modules found
        module_labels (List[str]): names or IDs for the modules
    """
    # 1) Extract gene expression matrix
    gex_data = adata_gex.X.toarray() if hasattr(adata_gex.X, 'toarray') else adata_gex.X
    # shape: N (spots) × G (genes)
    # We want gene-gene correlation, so transpose to G (genes) × N (spots)
    gex_data_t = gex_data.T  # shape: G × N

    # Handle potential NaN or inf values
    gex_data_t = np.nan_to_num(gex_data_t, nan=0.0, posinf=0.0, neginf=0.0)

    # 2) Compute gene-gene correlation with proper handling of NaN values
    if correlation_method == 'spearman':
        corr = np.zeros((gex_data_t.shape[0], gex_data_t.shape[0]))
        for i in range(gex_data_t.shape[0]):
            for j in range(i+1):
                if i == j:
                    corr[i,j] = 1.0
                else:
                    # Calculate correlation only on finite values
                    mask = np.isfinite(gex_data_t[i]) & np.isfinite(gex_data_t[j])
                    if np.sum(mask) > 1:  # Need at least 2 points for correlation
                        rho, _ = spearmanr(gex_data_t[i][mask], gex_data_t[j][mask])
                        corr[i,j] = rho if np.isfinite(rho) else 0.0
                        corr[j,i] = corr[i,j]
                    else:
                        corr[i,j] = corr[j,i] = 0.0
    else:
        # Use numpy corrcoef with proper handling
        corr = np.corrcoef(gex_data_t)
        corr = np.nan_to_num(corr, nan=0.0)  # Replace NaN with 0

    # 3) Convert correlation to distance for clustering
    distance_matrix = 1.0 - np.abs(corr)  # Use absolute correlation
    distance_matrix = np.clip(distance_matrix, 0, 1)  # Ensure distances are in [0,1]

    # 4) Convert square distance to condensed form
    # Ensure the matrix is symmetric
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    dist_condensed = squareform(distance_matrix, checks=False)

    # Verify no non-finite values
    if not np.all(np.isfinite(dist_condensed)):
        logging.warning("Non-finite values found in distance matrix. Replacing with max distance.")
        dist_condensed = np.nan_to_num(dist_condensed, nan=1.0, posinf=1.0, neginf=0.0)

    # 5) Perform hierarchical clustering
    Z = linkage(dist_condensed, method='average')

    # 6) Cut into clusters and create module matrix
    cluster_ids = fcluster(Z, t=max_clusters, criterion='maxclust')
    unique_clusters = np.unique(cluster_ids)
    K = len(unique_clusters)

    # 6) Build module_matrix (G × K), where entry=1 if gene in that module
    module_matrix = np.zeros((gex_data.shape[1], K), dtype=float)
    for i_gene, cluster_id in enumerate(cluster_ids):
        module_matrix[i_gene, cluster_id - 1] = 1.0

    module_labels = [f"Module_{c}" for c in unique_clusters]
    return module_matrix, module_labels

def optimize_multimodal_phase_3_celltype_wnn(
    adata_gex,
    adata_antibody,
    cell_prop_array,
    cell_profiles,
    radius: float = 2.0,
    n_neighbors: int = 30,
    alpha_rna: float = 0.5,
    alpha_protein: float = 0.5,
    lambda_smooth: float = 0.1,
    min_cells_per_cluster: int = 5
):
    """
    Integrate RNA and protein data using cell-type-level WNN approach.
    
    Args:
        adata_gex (AnnData): Gene expression data with cell-type-specific layers
        adata_antibody (AnnData): Protein/antibody data
        cell_prop_array (np.ndarray): Current cell type proportions (N × T)
        cell_profiles (dict): Cell type marker profiles
    """
    from sklearn.neighbors import NearestNeighbors
    from sklearn.cluster import KMeans
    import numpy as np
    
    N = adata_gex.shape[0]  # number of spots
    T = cell_prop_array.shape[1]  # number of cell types
    cell_type_names = list(cell_profiles.keys())
    
    # 1. Extract features for clustering
    
    # Cell type proportions features
    prop_features = cell_prop_array  # (N x T)
    
    # Cell type-specific gene expression features
    gex_features = []
    for t, ct_name in enumerate(cell_type_names):
        # Get layer name for this cell type
        layer_key = f"{ct_name.replace(' ', '_')}_genes_pass2"
        if layer_key not in adata_gex.layers:
            raise ValueError(f"Layer {layer_key} not found in adata_gex")
        
        # Extract and normalize cell type-specific expression
        ct_expression = adata_gex.layers[layer_key]
        if hasattr(ct_expression, 'toarray'):
            ct_expression = ct_expression.toarray()
        
        # Weight by cell type proportion
        ct_prop = cell_prop_array[:, t].reshape(-1, 1)
        weighted_expression = ct_expression * ct_prop
        gex_features.append(weighted_expression)
    
    # Combine features
    gex_features = np.hstack(gex_features)  # (N x (T*G))
    
    # 2. Cluster spots based on both proportions and expression
    n_clusters = max(int(N / min_cells_per_cluster), T)
    
    # Normalize features before clustering
    prop_features = (prop_features - prop_features.mean(axis=0)) / (prop_features.std(axis=0) + 1e-6)
    gex_features = (gex_features - gex_features.mean(axis=0)) / (gex_features.std(axis=0) + 1e-6)
    
    # Cluster using proportions
    prop_clusters = KMeans(n_clusters=n_clusters).fit_predict(prop_features)
    
    # Cluster using gene expression
    gex_clusters = KMeans(n_clusters=n_clusters).fit_predict(gex_features)
    
    # 3. Compute cluster centroids
    def get_cluster_centroids(data, clusters):
        centroids = np.zeros((n_clusters, data.shape[1]))
        for k in range(n_clusters):
            mask = clusters == k
            if np.any(mask):
                centroids[k] = np.mean(data[mask], axis=0)
        return centroids
    
    prop_centroids = get_cluster_centroids(prop_features, prop_clusters)
    gex_centroids = get_cluster_centroids(gex_features, gex_clusters)
    
    # 4. Find nearest neighbors between clusters
    nbrs_prop = NearestNeighbors(n_neighbors=n_neighbors).fit(prop_centroids)
    nbrs_gex = NearestNeighbors(n_neighbors=n_neighbors).fit(gex_centroids)
    
    dists_prop, idx_prop = nbrs_prop.kneighbors(prop_centroids)
    dists_gex, idx_gex = nbrs_gex.kneighbors(gex_centroids)
    
    # 5. Create Gurobi model for refinement
    model = gp.Model("Phase3_CellType_WNN")
    model.setParam('OutputFlag', 0)
    
    # Variables for refined proportions and expression
    Y_refined = {}  # Refined proportions
    X_refined = {}  # Refined expression
    
    for i in range(N):
        for t in range(T):
            Y_refined[i,t] = model.addVar(lb=0, ub=1, name=f'Y_{i}_{t}')
            for g in range(adata_gex.shape[1]):
                X_refined[i,t,g] = model.addVar(lb=0, name=f'X_{i}_{t}_{g}')
    
    # 6. Objective terms
    obj_terms = []
    
    # Proportion-based consistency
    for i in range(N):
        prop_cluster = prop_clusters[i]
        neighbor_clusters = idx_prop[prop_cluster]
        weights = np.exp(-dists_prop[prop_cluster] / dists_prop[prop_cluster].mean())
        
        for k, n_cluster in enumerate(neighbor_clusters):
            neighbor_spots = np.where(prop_clusters == n_cluster)[0]
            if len(neighbor_spots) > 0:
                w = weights[k] * alpha_protein
                for t in range(T):
                    for j in neighbor_spots:
                        diff = Y_refined[i,t] - Y_refined[j,t]
                        obj_terms.append(w * diff * diff)
    
    # Expression-based consistency
    for i in range(N):
        gex_cluster = gex_clusters[i]
        neighbor_clusters = idx_gex[gex_cluster]
        weights = np.exp(-dists_gex[gex_cluster] / dists_gex[gex_cluster].mean())
        
        for k, n_cluster in enumerate(neighbor_clusters):
            neighbor_spots = np.where(gex_clusters == n_cluster)[0]
            if len(neighbor_spots) > 0:
                w = weights[k] * alpha_rna
                for t in range(T):
                    for g in range(adata_gex.shape[1]):
                        for j in neighbor_spots:
                            diff = X_refined[i,t,g] - X_refined[j,t,g]
                            obj_terms.append(w * diff * diff)
    
    # Spatial smoothing if available
    if 'spatial' in adata_gex.obsm:
        spatial_coords = adata_gex.obsm['spatial']
        nbrs_spatial = NearestNeighbors(radius=radius).fit(spatial_coords)
        spatial_graph = nbrs_spatial.radius_neighbors_graph(spatial_coords)
        
        spatial_indices = spatial_graph.nonzero()
        for idx in range(len(spatial_indices[0])):
            i, j = spatial_indices[0][idx], spatial_indices[1][idx]
            # Smooth proportions
            for t in range(T):
                diff_prop = Y_refined[i,t] - Y_refined[j,t]
                obj_terms.append(lambda_smooth * diff_prop * diff_prop)
                # Smooth expression
                for g in range(adata_gex.shape[1]):
                    diff_expr = X_refined[i,t,g] - X_refined[j,t,g]
                    obj_terms.append(lambda_smooth * diff_expr * diff_expr)
    
    # Add objective
    model.setObjective(gp.quicksum(obj_terms), GRB.MINIMIZE)
    
    # 7. Constraints
    # Sum-to-one for proportions
    for i in range(N):
        model.addConstr(gp.quicksum(Y_refined[i,t] for t in range(T)) == 1)
    
    # Expression must be proportional to cell type abundance
    for i in range(N):
        for t in range(T):
            for g in range(adata_gex.shape[1]):
                model.addConstr(X_refined[i,t,g] <= adata_gex.X[i,g] * Y_refined[i,t])
    
    # 8. Optimize
    model.optimize()
    
    if model.status != GRB.OPTIMAL:
        raise RuntimeError("Failed to find optimal solution")
    
    # 9. Extract results
    refined_props = np.zeros((N, T))
    refined_profiles = np.zeros((N, T, adata_gex.shape[1]))
    
    for i in range(N):
        for t in range(T):
            refined_props[i,t] = Y_refined[i,t].X
            for g in range(adata_gex.shape[1]):
                refined_profiles[i,t,g] = X_refined[i,t,g].X
    
    return {
        'refined_proportions': refined_props,
        'refined_profiles': refined_profiles,
        'cluster_assignments': {
            'proportion_clusters': prop_clusters,
            'expression_clusters': gex_clusters
        }
    }

def normalize_counts(adata, target_sum=10000):
    """
    Normalize counts while preserving integer values and relative proportions.
    
    Args:
        adata: AnnData object
        target_sum: Target sum for each cell/spot
    
    Returns:
        Normalized AnnData object with counts scaled to target_sum
    """
    # Get matrix
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    
    # Calculate scaling factors
    size_factors = np.maximum(X.sum(axis=1), 1)  # Avoid division by zero
    median_size = max(np.median(size_factors), 1)  # Ensure positive median
    
    # Calculate scaling factors with bounds
    scaling_factors = np.clip(size_factors / median_size, 0.1, 10.0)  # Limit scaling range
    
    # Scale to target sum while preserving integers
    # Use a more controlled scaling approach
    scaled_factors = (target_sum / np.maximum(size_factors, 1))
    X_scaled = np.round(X * scaled_factors[:, None]).astype(int)
    
    # Add safety check for extreme values
    max_allowed = target_sum * 2  # Set reasonable maximum
    X_scaled = np.clip(X_scaled, 0, max_allowed)
    
    # Create new AnnData with scaled counts
    adata_norm = adata.copy()
    adata_norm.X = X_scaled
    
    # Store size factors and scaling info
    adata_norm.obs['size_factors'] = scaling_factors
    adata_norm.obs['original_total'] = size_factors
    adata_norm.obs['scaled_total'] = X_scaled.sum(axis=1)
    
    # Log statistics for validation
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

