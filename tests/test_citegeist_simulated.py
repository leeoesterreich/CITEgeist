# tests/test_citegeist_simulated.py

import os
import sys
import argparse
import logging
import gc
import json
from datetime import datetime
from typing import Dict, Any, List, Optional
import time

import numpy as np
import scanpy as sc
import pandas as pd
import scipy.sparse
from scipy.stats import kurtosis as scipy_kurtosis
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Spatial statistics imports
try:
    from esda.moran import Moran
    from libpysal.weights import KNN as LibPySAL_KNN
    HAS_ESDA = True
except ImportError:
    HAS_ESDA = False
    Moran = None
    LibPySAL_KNN = None

# Add the CITEgeist package directory to the system path
# The model directory is in CITEgeist/model/, not at the repo root
repo_root = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.join(repo_root, 'CITEgeist'))

# Now import using the full package path
from model.citegeist_model import CitegeistModel
from model.utils import benchmark_cell_proportions, calculate_expression_metrics, export_anndata_layers
from model.profile_matching import (
    match_profiles_to_ground_truth,
    create_remapped_proportions,
    benchmark_profile_discovery,
)

# New spatial colocalization pipeline (Module 1 + 2a + 2b + 2c)
from model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,  # Default: spatial variance-based selection
)


##############################################################################
# Spatial Visualization Functions for SMM Debugging
##############################################################################

def compute_morans_i(
    values: np.ndarray,
    coords: np.ndarray,
    k: int = 8,
) -> float:
    """
    Compute global Moran's I for a single variable.

    Args:
        values: 1D array of values (n_spots,)
        coords: Spatial coordinates (n_spots, 2)
        k: Number of nearest neighbors for spatial weights

    Returns:
        Moran's I statistic (float). Returns NaN if computation fails.
    """
    if not HAS_ESDA:
        logging.warning("esda/libpysal not available for Moran's I computation")
        return np.nan

    # Handle edge cases
    if len(values) < k + 1:
        return np.nan
    if np.var(values) < 1e-10:
        return np.nan  # No variance = undefined Moran's I

    try:
        # Create spatial weights
        w = LibPySAL_KNN.from_array(coords, k=k)
        w.transform = 'r'  # Row-standardize

        # Compute Moran's I
        mi = Moran(values, w)
        return mi.I
    except Exception as e:
        logging.warning(f"Moran's I computation failed: {e}")
        return np.nan


def visualize_marker_spatial_patterns(
    marker_names: List[str],
    coords: np.ndarray,
    output_dir: str,
    discovery_result,
    prefix: str = "marker_spatial",
    morans_k: int = 8,
    morans_i_threshold: float = 0.0,
    alpha: float = 0.1,
    spot_size: float = 30.0,
):
    """
    Generate 5-panel spatial diagnostic plots for markers showing all SMM pipeline stages.

    Pipeline: GMM(raw) -> Soft-scale(X * P(signal)) -> Smooth -> Z-score -> Moran's I

    Creates a 5-panel figure per marker showing:
    1. Raw expression
    2. Signal probability P(signal | x) from GMM
    3. Background-corrected (X * P(signal))
    4. Smoothed + Z-scored values
    5. Verdict with Moran's I and pass/fail status

    Args:
        marker_names: List of marker names to visualize
        coords: Spatial coordinates (n_spots, 2)
        output_dir: Directory to save output plots
        discovery_result: ProfileDiscoveryResult with SMM data
        prefix: Filename prefix for output plots
        morans_k: Number of neighbors for Moran's I computation
        morans_i_threshold: Threshold for passing Moran's I filter
        alpha: Significance level for Moran's I p-value
        spot_size: Size of spots in scatter plot
    """
    from model.utils import plot_marker_processing_stages

    # Create output directory for spatial diagnostics
    diag_dir = os.path.join(output_dir, "spatial_diagnostics")
    os.makedirs(diag_dir, exist_ok=True)

    # Check if SMM data is available
    if not discovery_result.smm_applied or discovery_result.smm_raw_matrix is None:
        logging.warning("SMM data not available in discovery_result, skipping visualization")
        print("  Warning: SMM data not available, skipping spatial diagnostic plots")
        return

    # Get all marker names from the original data
    all_markers = list(discovery_result.smm_snr_values.keys()) if discovery_result.smm_snr_values else []

    for marker in marker_names:
        if marker not in all_markers:
            logging.warning(f"Marker {marker} not in SMM results, skipping visualization")
            continue

        # Get marker index
        try:
            m_idx = all_markers.index(marker)
        except ValueError:
            logging.warning(f"Marker {marker} index not found, skipping")
            continue

        # Extract data for this marker
        raw_values = discovery_result.smm_raw_matrix[:, m_idx]
        signal_prob = discovery_result.smm_signal_posteriors[:, m_idx]
        corrected_values = discovery_result.smm_corrected_matrix[:, m_idx]
        smoothed_values = discovery_result.smm_smoothed_matrix[:, m_idx]

        # Z-score the smoothed values for Moran's I
        zscore_values = (smoothed_values - smoothed_values.mean()) / (smoothed_values.std() + 1e-10)

        # Compute Moran's I on Z-scored data
        morans_i = compute_morans_i(zscore_values, coords, k=morans_k)

        # Simple permutation test for p-value (quick approximation)
        n_perm = 99
        rng = np.random.default_rng(42)
        null_i = np.zeros(n_perm)
        for b in range(n_perm):
            perm_z = rng.permutation(zscore_values)
            null_i[b] = compute_morans_i(perm_z, coords, k=morans_k)
        p_value = (1 + np.sum(null_i >= morans_i)) / (1 + n_perm)

        # Determine if marker passed
        passed = morans_i >= morans_i_threshold and p_value < alpha

        # Get SNR and signal fraction for display
        snr = discovery_result.smm_snr_values.get(marker)
        signal_fraction = discovery_result.smm_signal_fractions.get(marker) if discovery_result.smm_signal_fractions else None

        # Generate plot
        safe_marker_name = marker.replace(' ', '_').replace('/', '_')
        output_path = os.path.join(diag_dir, f"{prefix}_{safe_marker_name}_stages.png")

        plot_marker_processing_stages(
            marker_name=marker,
            coords=coords,
            raw_values=raw_values,
            signal_prob=signal_prob,
            corrected_values=corrected_values,
            smoothed_values=smoothed_values,
            zscore_values=zscore_values,
            morans_i=morans_i,
            p_value=p_value,
            passed=passed,
            output_path=output_path,
            spot_size=spot_size,
            snr=snr,
            signal_fraction=signal_fraction,
        )

        status = "PASSED" if passed else "FILTERED"
        print(f"  Saved: {output_path} (Moran's I={morans_i:.3f}, p={p_value:.3f}, {status})")


# Representative GT markers for visualization (one per cell type)
GT_MARKERS_FOR_VISUALIZATION = {
    "Cancer Epithelial": "Cancer Epithelial_Protein_1",
    "Normal Epithelial": "Normal Epithelial_Protein_1",
    "CAFs": "CAFs_Protein_1",  # Fibroblasts
    "T-cells": "T-cells_Protein_1",
    "B-cells": "B-cells_Protein_1",
}

# Nonspecific markers for diagnostic visualization (to compare filtering behavior)
# These should have low Moran's I and help calibrate filtering thresholds
NONSPECIFIC_MARKERS_FOR_VISUALIZATION = [
    "Nonspecific_Protein_1",
    "Nonspecific_Protein_10",
    "Nonspecific_Protein_25",
    "Nonspecific_Protein_50",
    "Nonspecific_Protein_75",
]


def calculate_gex_metrics(ground_truth_dir, layer_dir, pass_number=None):
    """
    Calculate gene expression metrics and format them into a DataFrame.
    
    Args:
        ground_truth_dir (str): Directory containing ground truth files
        layer_dir (str): Directory containing prediction layers
        pass_number (int, optional): Pass number for logging
        
    Returns:
        pd.DataFrame: DataFrame containing metrics
    """
    metrics = calculate_expression_metrics(ground_truth_dir, layer_dir, normalize="range", pass_number=pass_number)

    # Filter out special tracking keys (not metric dictionaries)
    metric_keys = [k for k in metrics.keys() if not k.startswith('_')]

    # Report spurious and missed profiles if any
    spurious = metrics.get('_spurious_profiles', [])
    missed = metrics.get('_missed_profiles', [])
    if spurious:
        print(f"  Spurious profiles (not in ground truth): {spurious}")
    if missed:
        print(f"  Missed profiles (not predicted): {missed}")

    # Check if all metrics are not None or Nan, if they are, print the celltype
    for celltype in metric_keys:
        metric = metrics[celltype]
        if metric['RMSE'] is None or metric['NRMSE'] is None or metric['MAE'] is None or np.isnan(metric['RMSE']) or np.isnan(metric['NRMSE']) or np.isnan(metric['MAE']):
            print(f"Cell type {celltype} has None or NaN metrics")

    # Create DataFrame with metrics while excluding None or NaN values
    # Only use actual metric entries (not special keys)
    metric_values_list = [metrics[k] for k in metric_keys]

    metrics_values = {
        'Pass': [f"Pass {pass_number}" if pass_number else "Unknown"] * 6,
        'Metric': [
            'Average RMSE', 'Median RMSE', 'Average NRMSE',
            'Median NRMSE', 'Average MAE', 'Median MAE'
        ],
        'Value': [
            np.nanmean([m['RMSE'] for m in metric_values_list if m['RMSE'] is not None]),
            np.nanmedian([m['RMSE'] for m in metric_values_list if m['RMSE'] is not None]),
            np.nanmean([m['NRMSE'] for m in metric_values_list if m['NRMSE'] is not None]),
            np.nanmedian([m['NRMSE'] for m in metric_values_list if m['NRMSE'] is not None]),
            np.nanmean([m['MAE'] for m in metric_values_list if m['MAE'] is not None]),
            np.nanmedian([m['MAE'] for m in metric_values_list if m['MAE'] is not None])
        ]
    }
    return pd.DataFrame(metrics_values)

def calculate_improvements(pass1_metrics, pass2_metrics):
    """
    Calculate improvement percentages between passes.
    
    Args:
        pass1_metrics (pd.DataFrame): Metrics from pass 1
        pass2_metrics (pd.DataFrame): Metrics from pass 2
        
    Returns:
        pd.DataFrame: DataFrame containing improvement percentages
    """
    improvements = {}
    for metric in pass1_metrics['Metric'].unique():
        pass1_value = pass1_metrics[pass1_metrics['Metric'] == metric]['Value'].values[0]
        pass2_value = pass2_metrics[pass2_metrics['Metric'] == metric]['Value'].values[0]
        improvement = ((pass1_value - pass2_value) / pass1_value) * 100
        improvements[metric] = improvement
    
    return pd.DataFrame({
        'Metric': list(improvements.keys()),
        'Improvement_Percentage': list(improvements.values())
    })

##############################################################################
# Example cell-type profile dictionary for demonstration (adjust as needed).
##############################################################################
cell_type_profiles = {
    "B-cells": {
        "Major": ["B-cells_Protein_1", "B-cells_Protein_2"]
    },
    "CAFs": {
        "Major": ["CAFs_Protein_1", "CAFs_Protein_2"]
    },
    "Cancer Epithelial": {
        "Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]
    },
    "Endothelial": {
        "Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]
    },
    "Myeloid": {
        "Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]
    },
    "Normal Epithelial": {
        "Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]
    },
    "PVL": {
        "Major": ["PVL_Protein_1", "PVL_Protein_2"]
    },
    "Plasmablasts": {
        "Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]
    },
    "T-cells": {
        "Major": ["T-cells_Protein_1", "T-cells_Protein_2"]
    }
}

def main():
    """
    Example script to run CITEgeist on simulated data using the newer model-based
    implementation. This script demonstrates how to:
      1) Read in parameter arguments
      2) Load data
      3) Subset and prepare antibody capture
      4) Map antibodies to cell-type profiles
      5) Run the CITEgeist model to obtain cell-type proportions and gene-expression deconvolution
      6) Save the outputs
    """

    parser = argparse.ArgumentParser(description='Run CITEgeist on simulated data.')
    parser.add_argument('--radius', type=float, required=True, help='Radius for neighbor detection')
    parser.add_argument('--lambda_reg', type=float, required=True, 
                       help='Regularization strength for elastic net')
    parser.add_argument('--alpha_elastic', type=float, required=True, 
                       help='Elastic net mixing parameter (0=L2, 1=L1)')
    parser.add_argument('--max_y_change', type=float, required=True, 
                       help='Maximum allowed change in Y values between iterations (0,1)')
    parser.add_argument('--input_folder', type=str, default='.', help='Folder all requisite samples and ground truth')
    parser.add_argument('--output_folder', type=str, default='citegeist_output', help='Output folder')
    parser.add_argument('--sample_prefix', type=str, default='Wu_rep', help='Prefix to filter sample files')
    parser.add_argument('--profiling_only', action='store_true', default=False, 
                        help='If set, only compute cell-type proportions (no gene expression deconvolution).')
    parser.add_argument('--skip_pass2', action='store_true', default=False,
                        help='If set, skip pass 2 and only run pass 1.')
    parser.add_argument('--auto-profiles', action='store_true', default=False,
                        help='Use auto-discovered profiles via spatial colocalization pipeline')
    parser.add_argument('--top-k', type=int, default=3,
                        help='Mutual top-k for profile pair sparsification (default: 3, tested value)')
    parser.add_argument('--discovery-seed', type=int, default=1234,
                        help='Random seed for profile discovery reproducibility')
    # Spatial colocalization pipeline parameters (Modules 1 + 2a + 2b + 2c)
    parser.add_argument('--morans-k', type=int, default=8,
                        help='Number of neighbors for spatial statistics (default: 8)')
    parser.add_argument('--smooth-k', type=int, default=6,
                        help='Number of neighbors for spatial smoothing before Moran\'s I (default: 6)')
    parser.add_argument('--n-permutations', type=int, default=999,
                        help='Number of permutations for significance testing (default: 999)')
    parser.add_argument('--fdr-threshold', type=float, default=0.05,
                        help='FDR threshold for significant markers/pairs (default: 0.05)')
    # Module 1 marker filtering thresholds
    parser.add_argument('--gmm-snr-threshold', type=float, default=0.5,
                        help='Minimum GMM SNR for marker to be interesting (default: 0.5)')
    parser.add_argument('--kurtosis-threshold', type=float, default=None,
                        help='Fixed kurtosis threshold (default: None = adaptive)')
    parser.add_argument('--morans-threshold', type=float, default=None,
                        help='Fixed Moran\'s I threshold (default: None = adaptive)')
    # Module 2c: Spatial variance-based profile selection
    parser.add_argument('--variance-target', type=float, default=0.90,
                        help='Target fraction of spatial variance to explain (default: 0.90)')
    parser.add_argument('--min-marginal-gain', type=float, default=0.005,
                        help='Minimum marginal variance gain to add profile (default: 0.005)')
    # Laplacian smoothing parameters for proportion optimization
    parser.add_argument('--lambda-laplacian', type=float, default=0.1,
                        help='Laplacian smoothing weight for spatial coherence (default: 0.1, 0 to disable)')
    parser.add_argument('--laplacian-k', type=int, default=8,
                        help='Number of neighbors for Laplacian graph (default: 8)')

    # Joint optimization parameters (replaces sequential profile discovery + proportion estimation)
    parser.add_argument('--joint', action='store_true', default=False,
                        help='Use joint optimization for profile discovery + proportions (recommended)')
    parser.add_argument('--joint-min-K', type=int, default=2,
                        help='Minimum number of cell types for joint optimization (default: 2)')
    parser.add_argument('--joint-max-K', type=int, default=12,
                        help='Maximum number of cell types for joint optimization (default: 12)')
    parser.add_argument('--joint-lambda-spatial', type=float, default=0.1,
                        help='Laplacian spatial smoothing weight for joint optimization (default: 0.1)')
    parser.add_argument('--joint-lambda-sparsity', type=float, default=0.1,
                        help='L1 sparsity weight on profile weights W (default: 0.1)')
    parser.add_argument('--joint-lambda-distinct', type=float, default=0.5,
                        help='Profile distinctness penalty weight (default: 0.5)')
    parser.add_argument('--joint-max-markers', type=int, default=3,
                        help='Maximum markers per cell type profile (default: 3)')
    parser.add_argument('--joint-max-iterations', type=int, default=50,
                        help='Maximum alternating minimization iterations (default: 50)')
    parser.add_argument('--joint-n-restarts', type=int, default=3,
                        help='Number of random restarts per K (default: 3)')
    parser.add_argument('--joint-profile-threshold', type=float, default=0.3,
                        help='Minimum W weight to include marker in profile (default: 0.3)')

    args = parser.parse_args()

    radius = args.radius
    lambda_reg = args.lambda_reg
    alpha_elastic = args.alpha_elastic
    max_y_change = args.max_y_change
    variables = f"radius_{radius}_lambda_{lambda_reg}_alpha_{alpha_elastic}_max_y_change_{max_y_change}"

    input_folder = args.input_folder
    output_folder = args.output_folder

    suffix = "FilteredRadiiArrayWinsorCLRDiscreteErrorMinimizing"

    output_folder = os.path.join(output_folder, f'{variables}.', suffix + "CITEgeistOutput")

    # Create an output directory
    os.makedirs(output_folder, exist_ok=True)

    # Initialize logging
    log_file = f"Simulated_CITEgeist_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{radius}_{lambda_reg}_{alpha_elastic}.log"
    log_path = os.path.join(args.output_folder, log_file)
    logging.basicConfig(
        filename=log_path,
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG
    )
    logging.info(f"Starting CITEgeist run with parameters: radius={radius}, lambda_reg={lambda_reg}, alpha_elastic={alpha_elastic}")

    # Find all unique sample numbers for Wu_rep_{number} pairs
    h5ad_dir = os.path.join(args.input_folder, "h5ad_objects")
    sample_numbers = []
    for f in os.listdir(h5ad_dir):
        if f.startswith(args.sample_prefix):
            number = f.split('_')[2].split('.')[0]
            if number not in sample_numbers:
                sample_numbers.append(number)
    
    if len(sample_numbers) == 0:
        logging.error(f"No samples found in {h5ad_dir} matching prefix: {args.sample_prefix}")
        print(f"No samples found in {h5ad_dir} matching prefix: {args.sample_prefix}")
        sys.exit(1)

    for number in sample_numbers:
        start_time = time.time()  # Start timing for this sample
        sample_name = f"{args.sample_prefix}_{number}"
        logging.info(f"Processing sample: {sample_name}")

        adata_cite_path = os.path.join(args.input_folder, "h5ad_objects", f"{sample_name}_CITE.h5ad")
        adata_gex_path = os.path.join(args.input_folder, "h5ad_objects", f"{sample_name}_GEX.h5ad")
        
        # Verify files exist before loading
        if not (os.path.exists(adata_cite_path) and os.path.exists(adata_gex_path)):
            logging.error(f"Missing required files for {sample_name}")
            continue
            
        adata_cite = sc.read_h5ad(adata_cite_path)
        logging.info(f"Loaded {adata_cite_path} with shape {adata_cite.shape}")
        
        adata_gex = sc.read_h5ad(adata_gex_path)
        logging.info(f"Loaded {adata_gex_path} with shape {adata_gex.shape}")

        X = adata_gex.X
        if X is None:
            logging.warning("adata_gex.X is None")
            return
        if scipy.sparse.issparse(X):
            X = X.toarray()  # type: ignore
        if isinstance(X, np.ndarray) and np.any(X < 0):
            logging.warning("adata_gex.X contains negative values, which is not expected in count data.")
        

        ##############################################################################
        # Initialize the model
        ##############################################################################
        
        model = CitegeistModel(sample_name=sample_name, output_folder=output_folder,
                               simulation=True,
                               gene_expression_adata=adata_gex, antibody_capture_adata=adata_cite)

        model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)

        # Save raw antibody data BEFORE preprocessing (Module 1 needs raw data)
        X_antibody_raw = model.antibody_capture_adata.X.copy()
        if scipy.sparse.issparse(X_antibody_raw):
            X_antibody_raw = X_antibody_raw.toarray()

        # Preprocess datasets
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()  # Applies Winsorizing + CLR transformation

        # Load or discover cell profiles
        if args.joint:
            # Use joint optimization (recommended - replaces sequential approach)
            logging.info("Running JOINT optimization for profile discovery + proportions...")
            logging.info(f"  K range: [{args.joint_min_K}, {args.joint_max_K}]")
            logging.info(f"  lambda_spatial={args.joint_lambda_spatial}, lambda_sparsity={args.joint_lambda_sparsity}")
            logging.info(f"  lambda_distinct={args.joint_lambda_distinct}, max_markers={args.joint_max_markers}")

            joint_result = model.run_joint_optimization(
                min_K=args.joint_min_K,
                max_K=args.joint_max_K,
                max_markers_per_type=args.joint_max_markers,
                lambda_spatial=args.joint_lambda_spatial,
                lambda_sparsity=args.joint_lambda_sparsity,
                lambda_distinct=args.joint_lambda_distinct,
                laplacian_k=args.laplacian_k,
                max_iterations=args.joint_max_iterations,
                tolerance=1e-4,
                n_restarts=args.joint_n_restarts,
                seed=args.discovery_seed,
                profile_threshold=args.joint_profile_threshold,
                verbose=True,
            )

            logging.info(f"Joint optimization selected K={joint_result.K} (BIC={joint_result.bic:.2f})")
            logging.info(f"Discovered {len(joint_result.profiles) - 1} profiles: {list(joint_result.profiles.keys())}")

            # Store joint result proportions for benchmarking
            # The joint optimization already computed proportions, so we don't need run_cell_proportion_model
            spot_names = model.antibody_capture_adata.obs_names

            # Use index_to_name mapping for column names (matches profile names)
            cell_type_names = [joint_result.index_to_name.get(k, f"CellType_{k}") for k in range(joint_result.K)]

            # Create proportions DataFrame with proper cell type names
            global_cell_type_proportions_df = pd.DataFrame(
                joint_result.Y,
                index=spot_names,
                columns=cell_type_names,
            )
            finetuned_cell_type_proportions_df = global_cell_type_proportions_df.copy()

            # For benchmarking: store the discovery result structure
            # Create a minimal discovery_result-like object for compatibility
            class JointDiscoveryCompat:
                def __init__(self, joint_result):
                    self.profiles = joint_result.profiles
                    self.proportions = joint_result.Y
                    self.beta = {name: joint_result.beta[i] for i, name in enumerate(joint_result.marker_names)}
                    self.smm_applied = False
                    self.metadata = joint_result.metadata
                    self.index_to_name = joint_result.index_to_name

            discovery_result = JointDiscoveryCompat(joint_result)

        elif args.auto_profiles:
            # ============================================================
            # NEW: Spatial Colocalization Pipeline (Modules 1 + 2a + 2b + 2c)
            # ============================================================
            # This replaces the old integrate_with_model() approach with
            # a cleaner, modular pipeline based on spatial statistics.

            # Get spatial coordinates
            coords = model.antibody_capture_adata.obsm.get('spatial', None)
            if coords is None:
                coords = model.gene_expression_adata.obsm.get('spatial', None)
            if coords is None:
                raise ValueError(
                    "Spatial coordinates required for auto-profile discovery. "
                    "Ensure obsm['spatial'] is populated in the AnnData object."
                )

            logging.info("Running spatial colocalization pipeline for auto profile discovery...")

            # Get marker names
            marker_names = list(model.antibody_capture_adata.var_names)

            # NOTE: Use RAW antibody data (X_antibody_raw) for Module 1-2c
            # The CLR transformation changes kurtosis/Moran's I characteristics
            # significantly and causes adaptive thresholds to fail.
            # The preprocessed CLR data is still used for proportion optimization.

            # -----------------------------------------------------------
            # Module 1: Identify spatially interesting markers (RAW DATA)
            # -----------------------------------------------------------
            logging.info("Module 1: Identifying interesting markers (using raw data)...")
            marker_result = identify_interesting_markers(
                X=X_antibody_raw,
                coords=coords,
                marker_names=marker_names,
                morans_k=args.morans_k,
                smooth_k=args.smooth_k,  # Spatial smoothing before Moran's I
                morans_n_perm=args.n_permutations,
                gmm_snr_threshold=args.gmm_snr_threshold,
                kurtosis_threshold=args.kurtosis_threshold,
                morans_threshold=args.morans_threshold,
                seed=args.discovery_seed,
                verbose=True,
            )

            interesting_markers = marker_result.interesting_markers
            logging.info(f"Module 1: Found {len(interesting_markers)} interesting markers")
            print(f"\nModule 1: {len(interesting_markers)} spatially interesting markers identified")

            # Print learned thresholds and detected markers for debugging
            print(f"  Learned thresholds: kurtosis={marker_result.kurtosis_threshold:.2f}, morans_i={marker_result.morans_threshold:.3f}")
            print(f"  Interesting markers: {interesting_markers}")

            if len(interesting_markers) < 2:
                logging.warning("Not enough interesting markers found. Falling back to manual profiles.")
                print("  Warning: Not enough interesting markers. Using manual profiles.")
                model.load_cell_profile_dict(cell_type_profiles)
            else:
                # -----------------------------------------------------------
                # Module 2a: Analyze marker colocalization (RAW DATA)
                # -----------------------------------------------------------
                logging.info("Module 2a: Analyzing marker colocalization (using raw data)...")
                coloc_result = analyze_marker_colocalization(
                    X=X_antibody_raw,
                    coords=coords,
                    marker_names=marker_names,
                    markers_to_analyze=interesting_markers,
                    neighbor_k=args.morans_k,
                    smooth_k=args.smooth_k,  # Spatial smoothing before bivariate Moran's I
                    n_permutations=args.n_permutations,
                    seed=args.discovery_seed,
                    verbose=True,
                )

                logging.info(f"Module 2a: Found {len(coloc_result.pairs)} significant marker pairs")
                print(f"Module 2a: {len(coloc_result.pairs)} significant colocalization pairs")

                # -----------------------------------------------------------
                # Module 2b: Discover profiles via hierarchical clustering
                # -----------------------------------------------------------
                logging.info("Module 2b: Discovering profiles...")
                discovery_result = discover_profiles(
                    colocalization_result=coloc_result,
                    fdr_alpha=args.fdr_threshold,
                    top_k=args.top_k,
                    seed=args.discovery_seed,
                    verbose=True,
                )

                logging.info(f"Module 2b: Discovered {len(discovery_result.profiles)} candidate profiles")
                print(f"Module 2b: {len(discovery_result.profiles)} candidate profiles discovered")
                for i, profile in enumerate(discovery_result.profiles):
                    print(f"  {i+1}. {profile}")

                # -----------------------------------------------------------
                # Module 2c: Select profiles by spatial variance (RAW DATA)
                # -----------------------------------------------------------
                logging.info("Module 2c: Selecting profiles by spatial variance (using raw data)...")
                selection_result = select_profiles(
                    X=X_antibody_raw,
                    coords=coords,
                    marker_names=marker_names,
                    profiles=discovery_result.profiles,
                    interesting_markers=interesting_markers,
                    colocalization_result=coloc_result,
                    min_spatial_explained=args.variance_target,
                    min_marginal_gain=args.min_marginal_gain,
                    verbose=True,
                )

                selected_profiles = selection_result.selected_profiles
                n_selected = selection_result.optimal_n
                # Get metrics at the selected n (arrays are indexed from 0, n=1 is index 0)
                total_ve = float(selection_result.variance_explained[n_selected - 1]) if n_selected > 0 else 0.0
                ps = float(selection_result.proportion_smoothness[n_selected - 1]) if n_selected > 0 else 0.0
                logging.info(f"Module 2c: Selected {len(selected_profiles)} profiles "
                           f"(explains {total_ve:.1%} spatial variance)")
                print(f"\nModule 2c: Selected {len(selected_profiles)} profiles")
                print(f"  Spatial variance explained: {total_ve:.1%}")
                print(f"  Proportion smoothness: {ps:.3f}")
                print(f"  Stopping reason: {selection_result.stopping_reason}")

                for i, profile in enumerate(selected_profiles):
                    print(f"  {i+1}. {profile}")

                # -----------------------------------------------------------
                # Convert to cell_profile_dict format for CitegeistModel
                # -----------------------------------------------------------
                cell_profile_dict = {}
                for i, profile in enumerate(selected_profiles):
                    profile_name = f"Profile_{i+1}"
                    markers_list = list(profile) if not isinstance(profile, list) else profile
                    cell_profile_dict[profile_name] = {"Major": markers_list}

                model.load_cell_profile_dict(cell_profile_dict)
                logging.info(f"Loaded {len(cell_profile_dict)} profiles into model")

        else:
            # Use manual cell_type_profiles (existing behavior)
            model.load_cell_profile_dict(cell_type_profiles)

        # Skip explicit Gurobi registration - module load sets GRB_LICENSE_FILE env var
        # model.register_gurobi("/ihome/crc/install/gurobi/gurobi1102/linux64/lic/gurobi.lic")

        ##############################################################################
        # 1) Cell Proportion Inference
        ##############################################################################
        if args.joint:
            # Joint optimization already computed proportions - skip this step
            logging.info(f"Using proportions from joint optimization for {sample_name}.")
            # global_cell_type_proportions_df and finetuned_cell_type_proportions_df
            # were already set in the joint optimization block above
        else:
            # Sequential approach: run proportion optimization after profile discovery/loading
            logging.info(f"Running cell proportion model for {sample_name} ...")

            global_cell_type_proportions_df, finetuned_cell_type_proportions_df = model.run_cell_proportion_model(
                radius=radius,
                tolerance=1e-4,
                max_iterations=20,
                lambda_reg=lambda_reg,
                alpha=alpha_elastic,
                max_workers=None,
                checkpoint_interval=100,
                max_y_change=max_y_change,
                validation_warn_only=args.auto_profiles,  # Warnings only when using auto-discovered profiles
                skip_finetuning=True,  # Disable finetuning for benchmarking (incompatible with auto-profiles)
                # Laplacian smoothing parameters
                lambda_laplacian=args.lambda_laplacian,
                laplacian_k=args.laplacian_k,
            )

            logging.info(f"Completed cell proportion inference for {sample_name}.")

        # Benchmarking Cell Proportions
        st_folder = os.path.join(input_folder, "ST_sim")

        spot_composition_df = pd.read_csv(os.path.join(st_folder, f"Wu_ST_{number}_prop.csv"), index_col=0).sort_index().sort_index(axis=1)
        spot_composition_df = spot_composition_df.iloc[:, :-2]

        # Sort indices numerically by spot number
        def sort_spot_indices(df):
            # Extract numbers from spot names and sort
            df.index = pd.Index(df.index, name='spot')
            return df.reindex(sorted(df.index, key=lambda x: int(x.split('_')[1]) if '_' in x else float('inf')))

        spot_composition_df = sort_spot_indices(spot_composition_df)
        gt_cell_types = list(cell_type_profiles.keys())

        # If using auto-profiles or joint optimization, calculate discovery metrics and remap proportions
        if args.auto_profiles or args.joint:
            # Calculate profile discovery accuracy metrics
            discovery_metrics = benchmark_profile_discovery(
                model.cell_profile_dict,  # Discovered profiles
                cell_type_profiles,       # Ground truth profiles
            )
            logging.info(f"Profile Discovery Metrics: {discovery_metrics}")
            print(f"\nProfile Discovery Metrics:")
            print(f"  Profile Recovery Rate: {discovery_metrics['profile_recovery_rate']:.2%}")
            print(f"  Marker Precision: {discovery_metrics['marker_precision']:.2%}")
            print(f"  Marker Recall: {discovery_metrics['marker_recall']:.2%}")
            print(f"  False Discovery Rate: {discovery_metrics['false_discovery_rate']:.2%}")
            print(f"  Matched: {discovery_metrics['n_matched']}/{discovery_metrics['n_ground_truth']} cell types")
            if discovery_metrics['matched_pairs']:
                print(f"  Matched pairs: {discovery_metrics['matched_pairs']}")
            if discovery_metrics['missing_ground_truth']:
                print(f"  Missing GT types: {discovery_metrics['missing_ground_truth']}")

            # Save profile discovery metrics
            discovery_metrics_df = pd.DataFrame([{
                k: str(v) if isinstance(v, list) else v
                for k, v in discovery_metrics.items()
            }])
            discovery_metrics_path = os.path.join(
                output_folder,
                f'{sample_name}_profile_discovery_metrics_{suffix}.csv'
            )
            discovery_metrics_df.to_csv(discovery_metrics_path, index=False)
            logging.info(f"Saved profile discovery metrics to {discovery_metrics_path}")

            # Remap discovered profile names to ground truth cell type names
            match_result = match_profiles_to_ground_truth(
                model.cell_profile_dict,
                cell_type_profiles,
            )

            # Remap proportion DataFrames
            global_cell_type_proportions_df = create_remapped_proportions(
                global_cell_type_proportions_df,
                match_result,
                gt_cell_types,
            )
            finetuned_cell_type_proportions_df = create_remapped_proportions(
                finetuned_cell_type_proportions_df,
                match_result,
                gt_cell_types,
            )
            logging.info(f"Remapped proportions to GT cell types. Columns: {list(global_cell_type_proportions_df.columns)}")

        results_dict = {
            'global': global_cell_type_proportions_df,
            'finetune': finetuned_cell_type_proportions_df
        }

        for key, test_spots_df in results_dict.items():
            # Sort both DataFrames: rows numerically by spot number, columns alphabetically
            test_spots_df = sort_spot_indices(test_spots_df).sort_index(axis=1)
            spot_composition_df = sort_spot_indices(spot_composition_df).sort_index(axis=1)

            # Normalize column names (strip whitespace, convert to string)
            test_spots_df.columns = pd.Index([str(c).strip() for c in test_spots_df.columns])
            spot_composition_df.columns = pd.Index([str(c).strip() for c in spot_composition_df.columns])

            print(test_spots_df.index)
            print(spot_composition_df.index)

            # Verify that indices match
            if not np.array_equal(test_spots_df.index.astype(str), spot_composition_df.index.astype(str)):
                logging.warning(f"test_spots_df indices: {test_spots_df.index}, spot_composition_df indices: {spot_composition_df.index}")
                raise ValueError("ERROR: The row indices in the input CSV files do not match or are not in the same order!")

            # For auto-profiles or joint mode, we've already remapped columns, so align with GT columns
            if args.auto_profiles or args.joint:
                # Get common columns (GT cell types that were matched)
                common_cols = [c for c in spot_composition_df.columns if c in test_spots_df.columns]
                if not common_cols:
                    logging.error("No common columns between predicted and ground truth after remapping!")
                    continue
                # Subset both DataFrames to common columns
                test_spots_df = test_spots_df[common_cols]
                spot_composition_df_subset = spot_composition_df[common_cols]
            else:
                # For manual profiles, model adds "Unknown" cell type - exclude it for benchmark comparison
                test_cols = set(test_spots_df.columns) - {"Unknown"}
                gt_cols = set(spot_composition_df.columns)
                common_cols = sorted(test_cols & gt_cols)

                # Check if GT columns are covered (ignore "Unknown" from predictions)
                missing_in_test = gt_cols - test_cols
                if missing_in_test:
                    print(f"Column mismatch - Missing in predictions: {missing_in_test}")
                    print(f"test_spots_df columns ({len(test_spots_df.columns)}): {list(test_spots_df.columns)}")
                    print(f"spot_composition_df columns ({len(spot_composition_df.columns)}): {list(spot_composition_df.columns)}")
                    logging.warning(f"Column mismatch - Missing in predictions: {missing_in_test}")
                    raise ValueError("ERROR: The column names in the input CSV files do not match!")

                # Reorder both to common sorted order (excludes "Unknown" from predictions)
                test_spots_df = test_spots_df[common_cols]
                spot_composition_df_subset = spot_composition_df[common_cols]

            # Convert DataFrames to numpy arrays
            test_spots_metadata_mtrx = test_spots_df.values
            spot_composition_mtrx = spot_composition_df_subset.values


            column_names = test_spots_df.columns.tolist()

            results = benchmark_cell_proportions(test_spots_metadata_mtrx, spot_composition_mtrx, column_names)
            logging.info(f"{key.capitalize()} Cell proportion benchmarking results: {results}")

            # Save cell proportion results
            prop_results_df = pd.DataFrame([results])
            prop_results_name = os.path.join(output_folder, f'{sample_name}_cellprop_results_summary_{key}_{suffix}_{radius}_.csv')
            
            prop_results_df.to_csv(prop_results_name, index=False)
            logging.info(f"{key.capitalize()} Cell proportion results summary: \n{prop_results_df}")
            print(f"{key.capitalize()} Cell proportion results summary: \n{prop_results_df}")


        # Try to append finetuned proportions if available (may be disabled)
        try:
            model.append_proportions_to_adata(key='finetuned')
        except FileNotFoundError:
            logging.info("Finetuned proportions file not found - skipping (finetuning may be disabled)")
            # Fall back to global proportions
            try:
                model.append_proportions_to_adata(key='global')
            except FileNotFoundError:
                logging.warning("No proportions file found to append to adata")


        if args.profiling_only:
            logging.info("Skipping gene-expression deconvolution (profiling_only=True).")
            continue

        ##############################################################################
        # 2) Gene Expression Deconvolution - Pass 1
        ##############################################################################
        logging.info(f"Running pass 1 gene expression optimization for {sample_name} ...")

        # Run first pass with weight parameters
        pass1_results = model.run_cell_expression_pass1(
            radius=radius,
            max_workers=None, 
            checkpoint_interval=100, 
            output_dir="checkpoints", 
            rerun=True
        )

        # Calculate pass 1 metrics
        ground_truth_folder = os.path.join(input_folder, "ST_GEX_sim")
        ground_truth_dir = os.path.join(ground_truth_folder, f"sample_{number}", "layers")
        layer_dir_pass1 = os.path.join(output_folder, f"{sample_name}_pass1/layers")

        if os.path.exists(ground_truth_dir):
            logging.info("Calculating metrics for pass 1...")
            pass1_metrics = calculate_gex_metrics(ground_truth_dir, layer_dir_pass1, pass_number=1)

            
            assert pass1_metrics is not None, "Pass 1 metrics are None"
            
            print(f"Pass 1 metrics:\n{pass1_metrics}")
            print(f"Pass 1 metrics:\n{pass1_metrics}")
            logging.info(f"Pass 1 metrics:\n{pass1_metrics}")
            
            # Save pass 1 metrics
            metrics_path_pass1 = os.path.join(output_folder, f"{sample_name}_gex_metrics_pass1.csv")
            pass1_metrics.to_csv(metrics_path_pass1, index=False)

            ##############################################################################
            # 3) Compute Prior and Run Pass 2
            ##############################################################################
            logging.info("Computing prior from pass 1 results...")
            prior_info = model.compute_expression_prior(
                spotwise_profiles_pass1=pass1_results['spotwise_profiles'],
                cell_type_numbers_array=model.results['cell_prop'].values
            )

            
        end_time = time.time()
        runtime = end_time - start_time
        runtime_message = f"Runtime for sample {sample_name}: {runtime:.2f} seconds ({runtime/60:.2f} minutes)"
        print(runtime_message)
        logging.info(runtime_message)

if __name__ == "__main__":
    main()