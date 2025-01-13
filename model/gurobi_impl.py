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
    """Suggest threshold based on zero-inflation patterns."""
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

        # Modified objective using gene-specific enrichment and prior
        obj_terms = []
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    # Combine enrichment and prior with more weight on prior
                    enrichment_weight = gene_specific_enrichment[k, j]
                    prior_weight = np.clip(global_prior[j, k], 0.1, 1.0)  # Clip to avoid zero weights
                    cell_type_weight = neighborhood_cell_type_numbers[
                        len(neighborhood_indices)//2, j
                    ]
                    
                    combined_weight = (
                        (1 - lambda_prior_weight) * enrichment_weight +
                        lambda_prior_weight * prior_weight
                    )
                    
                    obj_terms.append(
                        combined_weight * cell_type_weight * X[j,k]
                    )

        model.setObjective(
            gp.quicksum(obj_terms),
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



def suggest_zero_inflation_threshold(zero_inflation_patterns, quantile=0.75):
    """
    Suggest a zero-inflation threshold based on the analysis results.
    """
    all_zero_props = []
    for ct_data in zero_inflation_patterns.values():
        for gene_data in ct_data.values():
            all_zero_props.append(gene_data['zero_proportion'])
    
    return np.quantile(all_zero_props, quantile)

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
    rerun=False,
    max_iterations=5  # Maximum number of iterations to prevent infinite loops
):
    """
    Iterative optimization that runs until the prior converges.
    """
    #########################
    # === Analyze Zero-Inflation Patterns
    #########################
    logging.info("=== Analyzing zero-inflation patterns from original data ===")
    
    # Get original expression matrix (raw counts)
    original_expression = filtered_adata.X
    if hasattr(original_expression, 'toarray'):
        original_expression = original_expression.toarray()
    
    # Analyze patterns
    zero_patterns = analyze_zero_inflation_patterns(
        original_expression,
        cell_type_numbers_array,
        filtered_adata.var_names,
        [f"CellType_{i}" for i in range(cell_type_numbers_array.shape[1])],
        min_proportion=0.01,
        expression_threshold=0.5
    )
    
    suggested_threshold = suggest_zero_inflation_threshold(zero_patterns, quantile=0.75)
    logging.info(f"Suggested zero-inflation threshold: {suggested_threshold:.3f}")

    #########################
    # === Iterative Optimization
    #########################
    iteration = 0
    prior_change = float('inf')
    previous_prior = None
    current_profiles = None

    while prior_change > prior_change_threshold and iteration < max_iterations:
        iteration += 1
        logging.info(f"\n=== Starting iteration {iteration} ===")

        #########################
        # === PASS 1 (No Prior or Using Previous Prior)
        #########################
        if iteration == 1:
            logging.info("=== Pass 1: Initial optimization without prior ===")
            spotwise_profiles_pass1 = optimize_gene_expression(
                sample_name=f"{sample_name}_iter{iteration}_pass1",
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
        else:
            logging.info(f"=== Pass 1 of iteration {iteration}: Using previous prior ===")
            # Create argument tuples for each spot using previous prior
            spot_args = [
                (i, filtered_adata, cell_type_numbers_array, radius, previous_prior,
                 lambda_prior_weight, alpha, lambda_reg_gex)
                for i in range(deconvolution_expression_data.shape[0])
            ]
            
            spotwise_profiles_pass1 = optimize_gene_expression_with_custom_spot_fn(
                sample_name=f"{sample_name}_iter{iteration}_pass1",
                deconvolution_expression_data=deconvolution_expression_data,
                cell_type_numbers_array=cell_type_numbers_array,
                filtered_adata=filtered_adata,
                radius=radius,
                alpha=alpha,
                lambda_reg_gex=lambda_reg_gex,
                spot_function=deconvolute_spot_with_neighbors_wrapper,
                spot_args=spot_args,
                max_workers=max_workers,
                checkpoint_interval=checkpoint_interval,
                output_dir=output_dir,
                rerun=rerun
            )

        #########################
        # === Compute New Prior
        #########################
        logging.info(f"=== Computing prior from iteration {iteration} ===")
        current_prior, zinb_matrix, zero_inflation_probs = compute_global_prior(
            spotwise_gene_expression_profiles=spotwise_profiles_pass1,
            cell_type_numbers_array=cell_type_numbers_array,
            lambda_prior=1.0,
            zero_inflation_threshold=suggested_threshold
        )

        # Calculate prior change if not first iteration
        if previous_prior is not None:
            prior_change = np.mean(np.abs(current_prior - previous_prior))
            logging.info(f"Prior change from previous iteration: {prior_change:.6f}")
        
        previous_prior = current_prior.copy()
        current_profiles = spotwise_profiles_pass1

        #########################
        # === PASS 2
        #########################
        logging.info(f"=== Pass 2 of iteration {iteration}: Optimization with current prior ===")
        
        spot_args = [
            (i, filtered_adata, cell_type_numbers_array, radius, current_prior,
             lambda_prior_weight, alpha, lambda_reg_gex)
            for i in range(deconvolution_expression_data.shape[0])
        ]

        current_profiles = optimize_gene_expression_with_custom_spot_fn(
            sample_name=f"{sample_name}_iter{iteration}_pass2",
            deconvolution_expression_data=deconvolution_expression_data,
            cell_type_numbers_array=cell_type_numbers_array,
            filtered_adata=filtered_adata,
            radius=radius,
            alpha=alpha,
            lambda_reg_gex=lambda_reg_gex,
            spot_function=deconvolute_spot_with_neighbors_wrapper,
            spot_args=spot_args,
            max_workers=max_workers,
            checkpoint_interval=checkpoint_interval,
            output_dir=output_dir,
            rerun=rerun
        )

    logging.info(f"\n=== Optimization completed after {iteration} iterations ===")
    if prior_change <= prior_change_threshold:
        logging.info("Converged: Prior change below threshold")
    else:
        logging.info("Stopped: Maximum iterations reached")

    return current_profiles

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