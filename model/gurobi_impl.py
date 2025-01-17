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
from scipy.optimize import minimize
from scipy.special import loggamma, digamma
import os
from pathlib import Path
import logging
import traceback
import gc
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import time
import concurrent.futures.process
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
        usage_array = usage_array[0]  # Take first pass data if dict
    
    # Convert to dense if sparse
    if hasattr(usage_array, 'toarray'):
        usage_array = usage_array.toarray()
    
    if usage_array.ndim != 2:
        raise ValueError(f"Expected 2D array, got shape {usage_array.shape}")
    
    N, M = usage_array.shape  # N spots, M genes
    
    for m_idx, gene in enumerate(gene_names):
        patterns[gene] = {}
        
        for t_idx, cell_type in enumerate(cell_type_names):
            # Get spots where this cell type is present
            ct_spots = cell_type_numbers[:, t_idx] > min_proportion
            
            if np.any(ct_spots):
                # Get expression values for these spots
                expr_values = usage_array[ct_spots, m_idx]
                
                if len(expr_values) > 0:
                    zero_prop = np.mean(expr_values <= expression_threshold)
                    # Store as dictionary with additional info
                    patterns[gene][cell_type] = {
                        'zero_proportion': zero_prop,
                        'mean_nonzero': np.mean(expr_values[expr_values > expression_threshold]) if np.any(expr_values > expression_threshold) else 0,
                        'n_spots': len(expr_values)
                    }
                else:
                    patterns[gene][cell_type] = {
                        'zero_proportion': np.nan,
                        'mean_nonzero': 0,
                        'n_spots': 0
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
    spotwise_gene_expression_profiles,
    cell_type_numbers_array,
    lambda_prior=1.0,
    zero_inflation_threshold=0.9
):
    """
    Compute global prior using ZINB fitting for discrete count data.
    
    Args:
        spotwise_gene_expression_profiles: Dict of spot profiles from first pass (now contains discrete counts)
        cell_type_numbers_array: Array of cell type proportions
        lambda_prior: Scaling factor for prior
        zero_inflation_threshold: Threshold for considering a gene zero-inflated
    """
    spot_indices = sorted(spotwise_gene_expression_profiles.keys())
    if not spot_indices:
        raise ValueError("No spotwise profiles found. Did the first pass run correctly?")
    T, M = spotwise_gene_expression_profiles[spot_indices[0]].shape
    N = len(spot_indices)

    # Construct usage array (now contains discrete counts)
    usage_array = np.zeros((N, T, M), dtype=float)
    for idx, spot_id in enumerate(spot_indices):
        usage_array[idx, :, :] = spotwise_gene_expression_profiles[spot_id]

    # Initialize matrices
    zinb_matrix = np.zeros((T, M), dtype=float)
    zero_inflation_probs = np.zeros((T, M), dtype=float)

    # For each cell type and gene
    for t_idx in range(T):
        ct_proportions = cell_type_numbers_array[:, t_idx]

        for m_idx in range(M):
            usage_vec = usage_array[:, t_idx, m_idx]
            
            try:
                # Fit ZINB and compute score
                pi, mu, theta = fit_zinb(usage_vec, ct_proportions)
                zero_inflation_probs[t_idx, m_idx] = pi
                
                if pi > zero_inflation_threshold:
                    zinb_score = 0.0  # Strong penalty for high zero-inflation
                else:
                    # Modified scoring for discrete counts
                    # Consider both frequency of expression and magnitude
                    nonzero_freq = np.mean(usage_vec > 0)
                    mean_magnitude = np.mean(usage_vec[usage_vec > 0]) if np.any(usage_vec > 0) else 0
                    
                    # Combined score considering both aspects
                    zinb_score = (1 - pi) * nonzero_freq * (mean_magnitude / (mean_magnitude + theta))
                    
                    # Additional boost for consistent expression
                    if nonzero_freq > 0.3:  # If expressed in >30% of spots
                        zinb_score *= 1.5
                        
                zinb_matrix[t_idx, m_idx] = zinb_score
            except:
                zinb_matrix[t_idx, m_idx] = 0.0
                zero_inflation_probs[t_idx, m_idx] = 1.0

    # Normalize matrix to [0,1] range per gene
    max_scores = np.maximum(np.max(zinb_matrix, axis=0), 1e-10)
    zinb_matrix = zinb_matrix / max_scores

    # Convert to prior via softmax
    global_prior = np.zeros((T, M), dtype=float)
    for m_idx in range(M):
        scores = zinb_matrix[:, m_idx]
        scaled = scores * lambda_prior * 2  # Stronger scaling for discrete data
        exps = np.exp(scaled - np.max(scaled))
        denom = np.sum(exps)
        if denom < 1e-12:
            global_prior[:, m_idx] = 1.0 / T
        else:
            global_prior[:, m_idx] = exps / denom

    return global_prior, zinb_matrix, zero_inflation_probs


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
    Error-minimizing discrete deconvolution with optional prior guidance.
    """
    model = None
    try:
        # Get neighborhood data first to establish dimensions
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

        # Extract expression data
        deconvolution_expression_data = filtered_adata.X
        if hasattr(deconvolution_expression_data, 'toarray'):
            deconvolution_expression_data = deconvolution_expression_data.toarray()
        
        # Get dimensions from the data
        T = cell_type_numbers_array.shape[1]  # number of cell types
        M = deconvolution_expression_data.shape[1]  # number of genes
        
        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]

        # Now we can safely create gene_specific_enrichment
        gene_specific_enrichment = np.zeros((M, T))

        def compute_expression_aware_enrichment(expression_data, cell_type_props, gene_idx):
            """
            Compute enrichment scores that account for both cell type presence
            and their gene expression patterns.
            """
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

        # Calculate enrichment scores for each gene
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

        # Build Gurobi model with relaxed parameters
        model = gp.Model(f"discrete_gene_expression_spot_{i}")
        model.setParam('OutputFlag', 0)
        model.setParam('Threads', 1)
        model.setParam('NodefileStart', 0.5)
        model.setParam('MIPGap', 0.1)  # Relaxed MIPGap for phase 2
        model.setParam('TimeLimit', 60)  # Shorter time limit per spot
        model.setParam('NodeLimit', 100000)
        model.setParam('FeasibilityTol', 1e-5)  # Slightly relaxed tolerance
        model.setParam('IntFeasTol', 1e-4)  # Slightly relaxed integer tolerance

        # Variables for count assignment
        X = {}
        center_counts = deconvolution_expression_data[i, :]

        # 1. Create error variables for each neighbor and gene
        error_vars = model.addVars(
            len(neighborhood_indices),
            M,
            lb=0, 
            ub=GRB.INFINITY,
            vtype=GRB.CONTINUOUS,
            name="error"
        )

        # 2. Define predicted expression for neighbors
        for n, spot_idx in enumerate(neighborhood_indices):
            total_prop_n = np.sum(neighborhood_cell_type_numbers[n, :])
            if total_prop_n < 1e-10:
                total_prop_n = 1e-10

            for k in range(M):
                if int(center_counts[k]) > 0:  # Only process non-zero genes
                    predicted_expr = gp.quicksum(
                        (neighborhood_cell_type_numbers[n, j] / total_prop_n) * X[j, k]
                        for j in range(T)
                    )
                    
                    observed_expr_nk = neighborhood_expression_data[n, k]
                    
                    # 3. Add absolute difference constraints
                    model.addConstr(error_vars[n, k] >= observed_expr_nk - predicted_expr)
                    model.addConstr(error_vars[n, k] >= predicted_expr - observed_expr_nk)

        # 4. Build minimization objective
        total_error = gp.quicksum(error_vars[n, k] 
                                 for n in range(len(neighborhood_indices)) 
                                 for k in range(M))

        # Optional regularization
        l1_term = gp.quicksum(X[j,k] for j in range(T) for k in range(M) 
                             if int(center_counts[k]) > 0)
        l2_term = gp.quicksum(X[j,k]*X[j,k] for j in range(T) for k in range(M) 
                             if int(center_counts[k]) > 0)
        regularization_term = lambda_reg_gex * (alpha*l1_term + (1-alpha)*l2_term)

        # Optional prior penalty (only used in pass 2)
        if global_prior is not None and lambda_prior_weight > 0:
            prior_penalty = lambda_prior_weight * gp.quicksum(
                (1 - global_prior[j, k]) * X[j, k]
                for j in range(T) for k in range(M)
                if int(center_counts[k]) > 0
            )
        else:
            prior_penalty = 0

        # Set final objective
        model.setObjective(
            total_error + regularization_term + prior_penalty,
            GRB.MINIMIZE
        )

        # Create variables and basic constraints
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

        # Add ZINB constraints with relaxed thresholds
        zero_inflation_threshold = 0.95  # Relaxed threshold
        max_allocation_fraction = 0.4  # Relaxed allocation limit
        
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    if global_prior[j, k] < zero_inflation_threshold:
                        # Add as a soft constraint using penalty in objective instead
                        penalty_weight = 1.0 - global_prior[j, k]
                        model.setObjective(
                            model.getObjective() - penalty_weight * X[j,k],
                            GRB.MAXIMIZE
                        )

        # Optimize with try-catch for numerical issues
        try:
            model.optimize()
        except Exception as e:
            logging.warning(f"Optimization error for spot {i}: {str(e)}")
            model.setParam('NumericFocus', 3)  # Try again with higher numeric focus
            model.optimize()

        if model.status == GRB.OPTIMAL or model.status == GRB.SUBOPTIMAL:
            result = np.zeros((T, M))
            for k in range(M):
                total_counts = int(center_counts[k])
                if total_counts > 0:
                    for j in range(T):
                        result[j,k] = X[j,k].X
            return result
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
    try:
        i, filtered_adata, cell_type_numbers_array, radius, global_prior, lambda_prior_weight, alpha, lambda_reg_gex = args
        return deconvolute_spot_with_neighbors_with_prior(
            i=i,
            filtered_adata=filtered_adata,
            cell_type_numbers_array=cell_type_numbers_array,
            radius=radius,
            global_prior=global_prior,
            lambda_prior_weight=lambda_prior_weight,
            alpha=alpha,
            lambda_reg_gex=lambda_reg_gex
        )
    except Exception as e:
        logging.error(f"Error in deconvolute_spot_with_neighbors_wrapper: {str(e)}")
        logging.error(traceback.format_exc())
        return None



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
                        future = executor.submit(
                            deconvolute_spot_with_neighbors,
                            spot_idx,
                            filtered_adata,
                            cell_type_numbers_array,
                            radius,
                            alpha,
                            lambda_reg_gex
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

def deconvolute_spot_with_neighbors(
    i, 
    filtered_adata, 
    cell_type_numbers_array, 
    radius, 
    alpha=0.5, 
    lambda_reg_gex=0.0001,
    local_weight=0.5,
    global_weight=0.5
):
    """
    Deconvolution with robustness to gene expression heterogeneity within cell types.
    """
    model = None
    try:
        # Get neighborhood data first to establish dimensions
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

        # Extract expression data
        deconvolution_expression_data = filtered_adata.X
        if hasattr(deconvolution_expression_data, 'toarray'):
            deconvolution_expression_data = deconvolution_expression_data.toarray()
        
        # Get dimensions from the data
        T = cell_type_numbers_array.shape[1]  # number of cell types
        M = deconvolution_expression_data.shape[1]  # number of genes
        
        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]

        # Now we can safely create gene_specific_enrichment
        gene_specific_enrichment = np.zeros((M, T))

        # Compute expression-aware enrichment for each gene
        def compute_expression_aware_enrichment(expression_data, cell_type_props, gene_idx):
            """
            Compute enrichment scores that account for both cell type presence
            and their gene expression patterns.
            """
            gene_expr = expression_data[:, gene_idx]
            expr_threshold = np.percentile(gene_expr[gene_expr > 0], 50) if np.any(gene_expr > 0) else 0
            high_expr_spots = gene_expr >= expr_threshold
            
            if not np.any(high_expr_spots):
                return np.ones(cell_type_props.shape[1]) / cell_type_props.shape[1]
            
            # Calculate cell type proportions in high-expression spots
            high_expr_props = np.mean(cell_type_props[high_expr_spots], axis=0)
            background_props = np.mean(cell_type_props, axis=0)
            
            # Compute enrichment scores
            epsilon = 1e-10
            enrichment = high_expr_props / (background_props + epsilon)
            
            # Smooth the enrichment scores
            smoothed_enrichment = 0.8 * enrichment + 0.2 * np.ones_like(enrichment)
            
            # Normalize
            return smoothed_enrichment / (np.sum(smoothed_enrichment) + epsilon)

        # Compute expression-aware enrichment for each gene
        gene_specific_enrichment = np.zeros((M, T))

        for k in range(M):
            # Local enrichment (in neighborhood)
            local_enrich = compute_expression_aware_enrichment(
                neighborhood_expression_data,
                neighborhood_cell_type_numbers,
                k
            )
            
            # Global enrichment (across all spots)
            global_enrich = compute_expression_aware_enrichment(
                deconvolution_expression_data,
                cell_type_numbers_array,
                k
            )
            
            # Combine local and global enrichment
            gene_specific_enrichment[k] = (
                local_weight * local_enrich +
                global_weight * global_enrich
            )

        # Build Gurobi model (similar to before)
        model = gp.Model(f"discrete_gene_expression_spot_{i}")
        model.setParam('OutputFlag', 0)
        model.setParam('Threads', 1)
        model.setParam('NodefileStart', 0.5)
        model.setParam('MIPGap', 0.01)
        model.setParam('TimeLimit', 600)
        model.setParam('NodeLimit', 1000000)

        T = neighborhood_cell_type_numbers.shape[1]
        M = neighborhood_expression_data.shape[1]

        # Variables for count assignment
        X = {}
        center_counts = deconvolution_expression_data[i, :]

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

        # Modified objective using gene-specific enrichment
        obj_terms = []
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    # Use gene-specific enrichment scores
                    enrichment_weight = gene_specific_enrichment[k, j]
                    
                    # Current spot's cell type proportion
                    cell_type_weight = neighborhood_cell_type_numbers[
                        len(neighborhood_indices)//2, j
                    ]
                    
                    # Add stochasticity to allow for heterogeneity
                    randomness = 0.9 + 0.2 * np.random.random()  # Random factor between 0.9 and 1.1
                    
                    obj_terms.append(
                        enrichment_weight * cell_type_weight * randomness * X[j,k]
                    )

        model.setObjective(
            gp.quicksum(obj_terms),
            GRB.MAXIMIZE
        )

        # Optimize
        model.optimize()

        if model.status == GRB.OPTIMAL:
            logging.info(f"Solution found for spot {i}")
            # Convert solution to array format
            result = np.zeros((T, M))
            for k in range(M):
                total_counts = int(center_counts[k])
                if total_counts > 0:
                    for j in range(T):
                        result[j,k] = X[j,k].X
            return result
        else:
            logging.error(f"No feasible solution found for spot {i}.")
            return None

    except Exception as e:
        logging.error(f"Error during discrete deconvolution of spot {i}: {str(e)}")
        logging.error(traceback.format_exc())
        return None

    finally:
        if model:
            del model
        gc.collect()


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
    Similar to optimize_gene_expression but with custom spot function.
    """
    if spot_function is None:
        raise ValueError("`spot_function` must be provided.")

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

    futures = {}
    
    try:
        workers = max_workers if max_workers is not None else max(1, os.cpu_count() // 2)
        logging.info(f"Using {workers} workers")
        
        remaining_spots = [i for i in range(N) if i not in completed_spots]
        logging.info(f"Processing {len(remaining_spots)} remaining spots")
        
        retry_count = 0
        max_retries = 3
        while retry_count < max_retries:
            try:
                with ProcessPoolExecutor(max_workers=workers) as executor:
                    futures.clear()
                    for spot_idx in remaining_spots:
                        if spot_idx < len(spot_args):
                            future = executor.submit(spot_function, spot_args[spot_idx])
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

    except Exception as e:
        logging.error(f"Error during parallel processing: {str(e)}")
        logging.error(traceback.format_exc())
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

    # 2) Compute gene-gene correlation
    if correlation_method == 'spearman':
        from scipy.stats import spearmanr
        corr, _ = spearmanr(gex_data_t, axis=1)
    else:
        corr = np.corrcoef(gex_data_t)

    # 3) Convert correlation to distance for clustering
    #    distance = 1 - corr
    distance_matrix = 1.0 - corr

    # 4) Perform hierarchical clustering on distance_matrix
    #    Convert square distance to condensed form
    dist_condensed = squareform(distance_matrix, checks=False)
    Z = linkage(dist_condensed, method='average')  # average linkage, or "ward"/"complete"

    # 5) Cut into clusters (modules)
    cluster_ids = fcluster(Z, t=max_clusters, criterion='maxclust')  # 1..max_clusters
    unique_clusters = np.unique(cluster_ids)
    K = len(unique_clusters)  # actual number of modules

    # 6) Build module_matrix (G × K), where entry=1 if gene in that module
    module_matrix = np.zeros((gex_data.shape[1], K), dtype=float)
    for i_gene, cluster_id in enumerate(cluster_ids):
        module_matrix[i_gene, cluster_id - 1] = 1.0

    module_labels = [f"Module_{c}" for c in unique_clusters]
    return module_matrix, module_labels

def optimize_multimodal_phase_3_wgcna(
    adata_gex,
    adata_antibody,
    cell_prop_array,
    cell_profiles,
    max_clusters=10,
    alpha_gex=1.0,
    alpha_antibody=1.0,
    alpha_spatial=0.0,
    lambda_reg_module=0.1,
    spatial_adjacency=None,
    random_seed=42
):
    """
    Integrate approximate WGCNA modules with an antibody-based term referencing 'cell_profiles'.

    Args:
        adata_gex (AnnData): Spots × Genes
        adata_antibody (AnnData): Spots × Markers
        cell_prop_array (np.ndarray): (N × T) existing cell-type proportions
        cell_profiles (dict): e.g. your Major/Minor marker dictionary
        max_clusters (int): upper bound on the number of gene modules
        alpha_gex, alpha_antibody (float): weights for GEX vs. antibody alignment
        alpha_spatial (float): weight for spatial term
        lambda_reg_module (float): L2 or L1 penalty for module usage
        spatial_adjacency (np.ndarray or None): (N × N) adjacency for spots
        random_seed (int): random seed

    Returns:
        results (dict):
          - modules: (G × K) matrix of module memberships
          - updated_cell_prop: (N × T) updated cell-type proportions
          - layers: (N × T × G) final refined expression per (spot, celltype, gene)
    """
    np.random.seed(random_seed)

    # 1) Approximate WGCNA to get (G × K) module membership matrix
    module_matrix, module_labels = approximate_wgcna(adata_gex, max_clusters=max_clusters)
    G, K = module_matrix.shape  # G = #genes, K = #modules

    # Data shapes
    gex_data = adata_gex.X.toarray() if hasattr(adata_gex.X, 'toarray') else adata_gex.X
    antibody_data = adata_antibody.X.toarray() if hasattr(adata_antibody.X, 'toarray') else adata_antibody.X
    N = gex_data.shape[0]  # spots
    T = cell_prop_array.shape[1]  # cell types
    M = antibody_data.shape[1]    # #markers

    logging.info(f"[Phase 3 (WGCNA)] N={N}, G={G}, K={K}, T={T}, M={M}")

    # 2) Optionally build a marker→celltype weighting from cell_profiles
    #    For simplicity, we create W[m, c] = 1 if marker m belongs to cell c, else 0.
    marker_var_names = list(adata_antibody.var_names)
    cell_type_names = list(cell_profiles.keys())
    marker_celltype_map = np.zeros((M, T), dtype=float)

    for t_idx, ct_name in enumerate(cell_type_names):
        major_markers = cell_profiles[ct_name].get("Major", [])
        # Convert your marker naming to var_names
        for mk in major_markers:
            # example mk = "CD68-1"
            if mk in marker_var_names:
                m_idx = marker_var_names.index(mk)
                marker_celltype_map[m_idx, t_idx] = 1.0

    # 3) Create Gurobi model
    model = gp.Model("Phase3_WGCNA")
    model.setParam('OutputFlag', 0)
    model.setParam('Seed', random_seed)

    # 4) Create variables:
    #    - X_module[i, k]: usage of module k in spot i
    #    - Y[i, c]: refined cell-type proportions
    #    (We also rely on module_matrix[G×K] for membership of gene→module.)
    X_module = {}
    for i in range(N):
        for k_ in range(K):
            X_module[(i, k_)] = model.addVar(lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                             name=f"X_module_{i}_{k_}")

    Y_refined = {}
    for i in range(N):
        for c in range(T):
            Y_refined[(i, c)] = model.addVar(lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS,
                                             name=f"Y_{i}_{c}")

    model.update()

    # 5) Objective Terms
    # a) GEX reconstruction error
    gex_residuals = []
    for i in range(N):
        for g_ in range(G):
            # predicted expression for gene g_ in spot i = sum_k( X_module[i,k] * module_matrix[g_,k] )
            pred_gex = gp.quicksum(X_module[(i, k_)] * module_matrix[g_, k_] for k_ in range(K))
            diff = gex_data[i, g_] - pred_gex
            gex_residuals.append(diff * diff)

    # b) Antibody-based constraints
    #    predicted antibody for marker m in spot i = sum_c( Y[i,c] * marker_celltype_map[m,c] ).
    #    We keep it simple, but you might do something more advanced.
    ab_residuals = []
    for i in range(N):
        for m_ in range(M):
            pred_ab = gp.quicksum(Y_refined[(i, c_)] * marker_celltype_map[m_, c_] for c_ in range(T))
            diff = antibody_data[i, m_] - pred_ab
            ab_residuals.append(diff * diff)

    # c) Spatial smoothing for Y
    spatial_terms = []
    if spatial_adjacency is not None and alpha_spatial > 0.0:
        # adjacency is NxN
        for i in range(N):
            for j in range(N):
                if spatial_adjacency[i, j] > 0:
                    for c_ in range(T):
                        delta = Y_refined[(i, c_)] - Y_refined[(j, c_)]
                        spatial_terms.append(delta * delta * spatial_adjacency[i, j])

    # d) Regularization on X_module usage
    #    L2 penalty: sum(X_module[i,k]^2)
    reg_module = []
    for i in range(N):
        for k_ in range(K):
            reg_module.append(X_module[(i, k_)] * X_module[(i, k_)])

    # Combine objectives
    obj = (
        alpha_gex * gp.quicksum(gex_residuals) +
        alpha_antibody * gp.quicksum(ab_residuals) +
        alpha_spatial * gp.quicksum(spatial_terms) +
        lambda_reg_module * gp.quicksum(reg_module)
    )
    model.setObjective(obj, GRB.MINIMIZE)

    # 6) Constraints
    # Force sum of Y[i,c] ~ 1
    for i in range(N):
        model.addConstr(
            gp.quicksum(Y_refined[(i, c_)] for c_ in range(T)) == 1.0,
            name=f"SumToOne_spot_{i}"
        )
        # Initialize Y with existing cell_prop_array
        for c_ in range(T):
            Y_refined[(i, c_)].start = float(cell_prop_array[i, c_])

    model.update()

    # 7) Optimize
    model.optimize()
    if model.status != GRB.OPTIMAL:
        logging.warning(f"WGCNA Phase 3: Model ended with status {model.status}")

    # 8) Extract results
    X_mod_sol = np.zeros((N, K), dtype=float)
    for i in range(N):
        for k_ in range(K):
            X_mod_sol[i, k_] = X_module[(i, k_)].X

    Y_sol = np.zeros((N, T), dtype=float)
    for i in range(N):
        for c_ in range(T):
            Y_sol[i, c_] = Y_refined[(i, c_)].X

    # 9) Construct updated layers: shape (N x T x G)
    #    e.g. each cell type's expression for gene g in spot i
    #    could be cell proportion * fraction of X_module usage for that gene's module.
    layers = np.zeros((N, T, G), dtype=float)
    for i in range(N):
        for c_ in range(T):
            for g_ in range(G):
                # sum across modules that gene belongs to
                # gene g_ is in module k_ if module_matrix[g_, k_] > 0
                val_mod = sum(X_mod_sol[i, k_] * module_matrix[g_, k_] for k_ in range(K))
                layers[i, c_, g_] = Y_sol[i, c_] * float(val_mod)

    results = {
        "modules": module_matrix,       # G x K
        "updated_cell_prop": Y_sol,     # N x T
        "layers": layers                # N x T x G
    }
    return results