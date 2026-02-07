"""
Main CitegeistModel class for spatial transcriptomics deconvolution.

This module implements the two-pass CITEgeist algorithm for deconvolving
spatial transcriptomics data using CITE-seq profiles.
"""
# Standard library imports
import logging
import os
from typing import Any, Dict, List, Optional, Tuple, Union

# Third-party imports
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc
from scipy.ndimage import gaussian_filter

# Local imports
from .gurobi_impl import (
    compute_global_prior,
    compute_marker_exclusivity,
    estimate_true_expression_cell,
    finetune_cell_proportions,
    finetune_cell_proportions_per_marker,
    map_antibodies_to_profiles,
    map_antibodies_to_profiles_v2,
    normalize_counts,
    optimize_cell_proportions,
    optimize_cell_proportions_per_marker,
    optimize_gene_expression,
)
from .utils import (
    assert_neighborhood_size,
    cleanup_memory,
    compute_optimal_radius,
    export_anndata_layers,
    get_neighbors_with_fixed_radius,
    setup_logging,
    validate_cell_profile_dict,
)


RESOLUTION_DEFAULTS = {
    "spot": {
        "neighbor_k": 8,
        "morans_k": 8,
        "smooth_k": 6,
        "coloc_neighbor_k": 6,
        "coloc_multi_scale_k": [6, 12, 24, 48, 64],
        "laplacian_k": 8,
        "lambda_spatial": 0.1,
        "lambda_sparse": 0.0,
        "pass2_library_slack": 1.0,
    },
    "cell": {
        "neighbor_k": 50,
        "morans_k": 50,
        "smooth_k": 20,
        "coloc_neighbor_k": 30,
        "coloc_multi_scale_k": [20, 40, 60, 80, 100],
        "laplacian_k": 50,
        "lambda_spatial": 0.01,
        "lambda_sparse": 0.1,
        "pass2_library_slack": 1.5,
    },
}


class CitegeistModel:
    def __init__(
        self,
        sample_name,
        adata=None,
        output_folder=None,
        simulation=False,
        gene_expression_adata=None,
        antibody_capture_adata=None,
        resolution="spot",
        resolution_overrides=None,
    ):
        """
        Initialize the CitegeistModel with an AnnData object and output folder.

        Args:
            adata (AnnData, optional): Spatial transcriptomics data object.
            output_folder (str): Path to save results and outputs.
            simulation (bool): Flag indicating if the data comes from a simulation framework.
            gene_expression_adata (AnnData, optional): Gene expression AnnData object (for simulations).
            antibody_capture_adata (AnnData, optional): Antibody capture AnnData object (for simulations).
        """
        if simulation:
            if gene_expression_adata is None or antibody_capture_adata is None:
                raise ValueError(
                    "In simulation mode, both `gene_expression_adata` and `antibody_capture_adata` must be provided."
                )
            self.gene_expression_adata = gene_expression_adata
            self.antibody_capture_adata = antibody_capture_adata
            self.adata = None  # Clear `adata` since separate datasets are provided
        else:
            if adata is None:
                raise ValueError("In non-simulation mode, `adata` must be provided.")
            self.adata = adata
            self.gene_expression_adata = None
            self.antibody_capture_adata = None

        self.sample_name = sample_name

        if output_folder is None:
            raise ValueError("output_folder must be provided")
        self.output_folder = str(output_folder)  # Ensure string type

        os.makedirs(self.output_folder, exist_ok=True)
        setup_logging(self.output_folder, self.sample_name)

        self.results = {}
        self.cell_profile_dict = None
        self.preprocessed_gex = False
        self.preprocessed_antibody = False

        # Resolution mode and parameter presets
        if resolution not in RESOLUTION_DEFAULTS:
            raise ValueError(
                f"resolution must be one of {list(RESOLUTION_DEFAULTS.keys())}, got '{resolution}'"
            )
        self.resolution = resolution
        self.resolution_params = dict(RESOLUTION_DEFAULTS[resolution])
        if resolution_overrides:
            for key, val in resolution_overrides.items():
                if key not in self.resolution_params:
                    raise ValueError(f"Unknown resolution parameter: '{key}'")
                self.resolution_params[key] = val

        print("CitegeistModel initialized successfully.")

    def __repr__(self):
        """
        Developer-friendly representation of the CitegeistModel.
        """
        return (
            f"<CitegeistModel(adata={'Loaded' if self.adata else 'Not Loaded'}, "
            f"gene_expression_adata={'Loaded' if self.gene_expression_adata else 'Not Loaded'}, "
            f"antibody_capture_adata={'Loaded' if self.antibody_capture_adata else 'Not Loaded'}, "
            f"cell_profile_dict={'Loaded' if self.cell_profile_dict else 'Not Loaded'}, "
            f"preprocessed_gex={'Yes' if self.preprocessed_gex else 'No'}, "
            f"preprocessed_antibody={'Yes' if self.preprocessed_antibody else 'No'}, "
            f"output_folder='{self.output_folder}')>"
        )

    def __str__(self):
        """
        User-friendly representation of the CitegeistModel.
        """
        details = [
            "CitegeistModel Summary:",
            f"- Output Folder: {self.output_folder}",
            f"- Main AnnData Loaded: {'Yes' if self.adata else 'No'}",
            f"- Gene Expression AnnData Loaded: {'Yes' if self.gene_expression_adata else 'No'}",
            f"- Antibody Capture AnnData Loaded: {'Yes' if self.antibody_capture_adata else 'No'}",
            f"- Cell Profile Dictionary Loaded: {'Yes' if self.cell_profile_dict else 'No'}",
            f"- Gene Expression Preprocessed: {'Yes' if self.preprocessed_gex else 'No'}",
            f"- Antibody Capture Preprocessed: {'Yes' if self.preprocessed_antibody else 'No'}",
        ]
        return "\n".join(details)

    def register_gurobi(self, license_file_path):
        """
        Configure Gurobi by setting only the license file path.

        Args:
            license_file_path (str): Path to the Gurobi license file.
        """
        if not os.path.isfile(license_file_path):
            raise FileNotFoundError(f"❌ License file not found at: {license_file_path}")

        # Set only the license file environment variable
        os.environ["GRB_LICENSE_FILE"] = license_file_path

        print("✅ Gurobi license file has been successfully configured.")
        print(f" - GRB_LICENSE_FILE: {os.environ['GRB_LICENSE_FILE']}")

    # -----------------------------------------
    # Data Splitting
    # -----------------------------------------

    def split_adata(self):
        """
        Split the AnnData object into separate gene expression and antibody capture sub-objects
        based on 'feature_types' in `adata.var`.

        Returns:
            None
        """
        if self.adata is None:
            raise ValueError("No valid data loaded. Ensure `adata` or split datasets are loaded properly.")

        if "feature_types" not in self.adata.var.columns:
            raise ValueError("The 'feature_types' column is missing in `adata.var`. Cannot split data.")

        if self.adata is None:
            raise ValueError("No valid data loaded. Ensure `adata` or split datasets are loaded properly.")

        self.adata.var_names_make_unique()

        if self.gene_expression_adata or self.antibody_capture_adata:
            raise ValueError("Data seems to already be split")

        # Identify indices for Gene Expression and Antibody Capture
        gene_expression_idx = np.where(self.adata.var["feature_types"] == "Gene Expression")[0]
        antibody_capture_idx = np.where(self.adata.var["feature_types"] == "Antibody Capture")[0]

        if len(gene_expression_idx) == 0:
            raise ValueError("No 'Gene Expression' features found in `adata.var['feature_types']`.")
        if len(antibody_capture_idx) == 0:
            raise ValueError("No 'Antibody Capture' features found in `adata.var['feature_types']`.")

        # Split AnnData object
        self.gene_expression_adata = self.adata[:, gene_expression_idx].copy()
        self.antibody_capture_adata = self.adata[:, antibody_capture_idx].copy()

        print("AnnData has been successfully split into 'gene_expression_adata' and 'antibody_capture_adata'.")

    # -----------------------------------------
    # Utility Functions
    # -----------------------------------------
    @staticmethod
    def winsorize(matrix, lower_percentile=5, upper_percentile=95):
        """Winsorize a 2D NumPy array."""
        lower_bound = np.percentile(matrix, lower_percentile)
        upper_bound = np.percentile(matrix, upper_percentile)
        return np.clip(matrix, lower_bound, upper_bound)

    @staticmethod
    def row_normalize(matrix, target_sum=1e4):
        """Row normalize a 2D NumPy array to a fixed target sum."""
        row_sums = matrix.sum(axis=1, keepdims=True)
        normalized = (matrix / row_sums) * target_sum
        return normalized

    @staticmethod
    def global_clr(matrix, epsilon=1e-6):
        """
        Apply margin=2 CLR normalization (global geometric mean per marker).
        Args:
            matrix (numpy.ndarray): Input matrix (spots x markers).
            epsilon (float): Small constant to avoid division by zero.
        Returns:
            numpy.ndarray: CLR-normalized matrix.
        """
        matrix = matrix + epsilon  # Avoid division by zero
        geom_mean = np.exp(np.mean(np.log(matrix), axis=0))
        normalized_matrix = matrix / geom_mean
        return normalized_matrix

    def load_cell_profile_dict(self, cell_profile_dict):
        """
        Load and validate the cell profile dictionary.

        Args:
            cell_profile_dict (dict): Dictionary of cell type profiles.
        """
        if validate_cell_profile_dict(cell_profile_dict):
            self.cell_profile_dict = cell_profile_dict
        else:
            raise ValueError("Invalid cell_profile_dict format.")

    # -----------------------------------------
    # Preprocessing Functions
    # -----------------------------------------

    def filter_gex(self, nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=10):
        """
        Filter genes in the gene expression AnnData object based on user-defined criteria.

        Filters genes that have:
        1. A count > 0 in at least `nonzero_percentage` of spots
        2. Mean expression > `mean_expression_threshold` in nonzero spots

        Args:
            nonzero_percentage (float): Minimum percentage of spots where a gene must have a count > 0 (default: 1%)
            mean_expression_threshold (float): Minimum mean expression value in nonzero spots
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        # Extract the data matrix
        matrix = (
            self.gene_expression_adata.X.toarray()
            if hasattr(self.gene_expression_adata.X, "toarray")
            else self.gene_expression_adata.X
        )
        matrix = np.asarray(matrix)  # Ensure dense matrix

        # Calculate the number of spots
        num_spots = matrix.shape[0]

        # First filter: minimum percentage of nonzero spots
        count_filter = (matrix > 0).sum(axis=0) >= (nonzero_percentage * num_spots)

        # Calculate mean expression in nonzero spots for each gene
        nonzero_means = np.zeros(matrix.shape[1])
        for j in range(matrix.shape[1]):
            nonzero_vals = matrix[:, j][matrix[:, j] > 0]
            nonzero_means[j] = np.mean(nonzero_vals) if len(nonzero_vals) > 0 else 0

        # Second filter: mean expression in nonzero spots
        mean_filter = nonzero_means > mean_expression_threshold

        # Combine filters
        col_filter = count_filter & mean_filter

        # Apply the filter and subset the AnnData object
        filtered_gene_count = np.sum(col_filter)
        initial_gene_count = self.gene_expression_adata.shape[1]

        self.gene_expression_adata = self.gene_expression_adata[:, col_filter].copy()

        initial_spot_count = self.gene_expression_adata.shape[0]

        sc.pp.filter_cells(self.gene_expression_adata, min_counts=min_counts)

        print(
            f"Filtered gene expression data: {initial_gene_count} → {filtered_gene_count} genes "
            f"(count > 0 in at least {nonzero_percentage*100}% of spots, mean expression > {mean_expression_threshold} "
            f"in nonzero spots). Remaining spots: {self.gene_expression_adata.shape[0]} "
            f"Filtered spots: {initial_spot_count} to {self.gene_expression_adata.shape[0]}"
        )

    def copy_gex_to_protein_adata(self):
        """
        Copy the number of spots in the gene expression AnnData object to the antibody capture AnnData object.
        """
        if self.antibody_capture_adata is None:
            raise ValueError("Antibody capture data has not been split. Run `split_adata` first.")
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        # Get the spot names from gene expression data
        gex_spots = set(self.gene_expression_adata.obs_names)

        # Filter antibody capture data to keep only spots present in gene expression data
        filtered_spots = [spot for spot in self.antibody_capture_adata.obs_names if spot in gex_spots]

        if not filtered_spots:
            raise ValueError("No matching spots found between gene expression and antibody capture data.")

        self.antibody_capture_adata = self.antibody_capture_adata[filtered_spots, :].copy()

        logging.info(f"Filtered antibody capture data to {len(filtered_spots)} spots present in gene expression data.")

    def preprocess_gex(self, target_sum=10000):
        """
        Preprocess gene expression data with count-preserving normalization.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        # Normalize while preserving counts
        self.gene_expression_adata = normalize_counts(self.gene_expression_adata, target_sum=target_sum)

        # Validate integer format
        matrix = self.gene_expression_adata.X
        if hasattr(matrix, "toarray"):
            matrix = matrix.toarray()

        if not np.all(np.equal(np.mod(matrix, 1), 0)):
            raise ValueError("Gene expression data contains non-integer values after normalization.")

        self.preprocessed_gex = True
        logging.info(
            f"Gene expression data normalized to {target_sum} counts per spot and validated for discrete count analysis."
        )

    def preprocess_antibody(self):
        """
        Preprocess antibody capture data:
        - Winsorize extreme values.
        - Apply Gaussian smoothing for local background correction.
        - Apply global CLR normalization.
        - Raise an error if NaNs or Infs are detected in the processed data.

        """
        if self.antibody_capture_adata is None:
            raise ValueError("Antibody capture data has not been split. Run `split_adata` first.")

        # Step 1: Extract and ensure matrix is dense
        matrix = (
            self.antibody_capture_adata.X.toarray()
            if hasattr(self.antibody_capture_adata.X, "toarray")
            else self.antibody_capture_adata.X
        )
        matrix = np.asarray(matrix)

        # Step 2: Validate initial data (no NaNs or Infs at start)
        if np.isnan(matrix).any() or np.isinf(matrix).any():
            raise ValueError("Antibody capture matrix contains NaN or Inf values before preprocessing.")

        # Step 3: Winsorize to cap extreme values
        matrix = self.winsorize(matrix, lower_percentile=5, upper_percentile=95)

        column_sums = matrix.sum(axis=0)
        zero_columns = column_sums == 0
        if np.any(zero_columns):
            matrix[:, zero_columns] += 1e-6

        # Step 5: Apply CLR Normalization
        matrix = self.global_clr(matrix)

        # Step 6: Final Validation for NaNs or Infs
        if np.isnan(matrix).any() or np.isinf(matrix).any():
            raise ValueError("NaN or Inf values detected in antibody capture matrix after preprocessing.")

        # Step 7: Reassign processed matrix to AnnData object
        self.antibody_capture_adata.X = matrix

        # Update status flag
        self.preprocessed_antibody = True

        print("Antibody capture data preprocessing completed: Winsorized, CLR applied, no NaNs detected.")

    def run_cell_proportion_model(
        self,
        radius=None,
        tolerance=1e-4,
        max_iterations=20,
        lambda_reg=1,
        alpha=0.5,
        max_y_change=0.4,
        max_workers=None,
        checkpoint_interval=100,
        unknown_threshold=0.05,
        min_celltype_threshold=0.01,
        redundancy_threshold=0.1,
        validation_warn_only=False,
        skip_finetuning=False,
        # Laplacian smoothing parameters
        lambda_laplacian=0.1,
        laplacian_k=8,
        # Per-marker beta parameters
        per_marker_beta=True,
        beta_min=0.1,
        beta_max=2.0,
        # Cell classification parameters (cell resolution only)
        use_gating=None,
        priority_dict=None,
        threshold_method="auto",
        use_negative_gates=False,
    ):
        """
        Orchestrates the cell proportion optimization workflow.
        Delegates optimization to `optimize_cell_proportions` and `finetune_cell_proportions` in `gurobi_impl.py`.

        Args:
            tolerance (float): Convergence tolerance for EM algorithm
            max_iterations (int): Maximum number of iterations
            lambda_reg (float): Regularization strength
            alpha (float): L1-L2 tradeoff factor (0 = L2, 1 = L1)
            max_workers (int, optional): Maximum number of parallel workers for finetuning
            checkpoint_interval (int): Number of spots between checkpoints during finetuning
            unknown_threshold (float): Maximum allowed mean proportion for Unknown cell type (default: 0.05 = 5%)
            min_celltype_threshold (float): Minimum required mean proportion for defined cell types (default: 0.01 = 1%)
            redundancy_threshold (float): Maximum allowed fraction of redundant cell types (default: 0.1 = 10%)
            skip_finetuning (bool): If True, skip the finetuning step and return global proportions only
            lambda_laplacian (float): Weight for Laplacian spatial smoothing (default: 0.1, 0 to disable)
            laplacian_k (int): Number of neighbors for Laplacian graph (default: 8)
            per_marker_beta (bool): If True (default), use per-marker beta learning which preserves
                marker-level signal variation. If False, use legacy per-celltype beta with marker averaging.
            beta_min (float): Minimum allowed beta value for per-marker optimization (default: 0.1)
            beta_max (float): Maximum allowed beta value for per-marker optimization (default: 2.0)
        """

        # Use resolution preset for laplacian_k if caller used default
        if laplacian_k == 8 and self.resolution == "cell":
            laplacian_k = self.resolution_params["laplacian_k"]

        # Use resolution preset for lambda_laplacian if caller used default
        if lambda_laplacian == 0.1 and self.resolution == "cell":
            lambda_laplacian = self.resolution_params["lambda_spatial"]

        if radius is None:
            # Auto-detect optimal radius from spatial coordinates
            source_adata = self.antibody_capture_adata or self.gene_expression_adata
            if source_adata is None:
                raise ValueError("No AnnData available for radius auto-detection")
            radius = compute_optimal_radius(source_adata)
            logging.info(f"Auto-detected radius: {radius:.2f} (3 rings)")

        if self.adata is None and (self.gene_expression_adata is None or self.antibody_capture_adata is None):
            raise ValueError("No valid data loaded. Ensure `adata` or split datasets are loaded properly.")

        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary has not been loaded. Run `load_cell_profile_dict` first.")

        # Extract spatial coordinates for Laplacian smoothing
        coords = None
        if lambda_laplacian > 0:
            coords = self.antibody_capture_adata.obsm.get('spatial', None)
            if coords is None and self.gene_expression_adata is not None:
                coords = self.gene_expression_adata.obsm.get('spatial', None)
            if coords is not None:
                logging.info(f"Using Laplacian smoothing with lambda={lambda_laplacian}, k={laplacian_k}")
            else:
                logging.warning("No spatial coordinates found for Laplacian smoothing - disabling")
                lambda_laplacian = 0

        # Dispatch to gating-based classification for cell resolution
        if use_gating is None:
            use_gating = (self.resolution == "cell")
        if use_gating:
            logging.info("Cell resolution: dispatching to gating-based classification")
            return self._run_cell_classification(
                threshold_method=threshold_method,
                priority_dict=priority_dict,
                use_negative_gates=use_negative_gates,
                coords=coords,
            )

        spot_names = self.antibody_capture_adata.obs_names

        if per_marker_beta:
            # ====== NEW PER-MARKER BETA APPROACH ======
            logging.info("Using per-marker beta optimization (preserves marker-level signal variation)")

            # Get marker-level data instead of averaged profiles
            marker_level_data, marker_names, assignment_matrix, cell_type_names = map_antibodies_to_profiles_v2(
                self.antibody_capture_adata, self.cell_profile_dict
            )

            try:
                logging.info(f"Running Stage 1 cell proportion optimization with validation thresholds: "
                            f"Unknown<{unknown_threshold*100:.1f}%, CellTypes>{min_celltype_threshold*100:.1f}%, Redundancy<{redundancy_threshold*100:.0f}%")

                Y_values, beta_values, marker_beta_dict, alpha_values = optimize_cell_proportions_per_marker(
                    marker_level_data=marker_level_data,
                    marker_names=marker_names,
                    assignment_matrix=assignment_matrix,
                    cell_type_names=cell_type_names,
                    tolerance=tolerance,
                    max_iterations=max_iterations,
                    lambda_reg=lambda_reg,
                    alpha=alpha,
                    beta_min=beta_min,
                    beta_max=beta_max,
                    unknown_threshold=unknown_threshold,
                    min_celltype_threshold=min_celltype_threshold,
                    redundancy_threshold=redundancy_threshold,
                    warn_only=validation_warn_only,
                    lambda_laplacian=lambda_laplacian,
                    coords=coords,
                    laplacian_k=laplacian_k,
                    lambda_sparse=self.resolution_params.get("lambda_sparse", 0.0),
                )

                # Store marker betas and baselines for downstream analysis
                self.results["marker_beta"] = marker_beta_dict
                marker_alpha_dict = {marker_names[i]: alpha_values[i] for i in range(len(marker_names))}
                self.results["marker_alpha"] = marker_alpha_dict
                for m_idx, m_name in enumerate(marker_names):
                    if alpha_values[m_idx] > 0.05:
                        logging.info(f"  Marker baseline: {m_name} = {alpha_values[m_idx]:.3f}")

                # Compute marker exclusivity scores for finetuning
                marker_owners = []
                for m_idx in range(assignment_matrix.shape[0]):
                    owners = [j for j in range(assignment_matrix.shape[1]) if assignment_matrix[m_idx, j] > 0]
                    marker_owners.append(owners)

                marker_exclusivity = compute_marker_exclusivity(
                    marker_level_data=marker_level_data,
                    Y_values=Y_values,
                    marker_owners=marker_owners,
                    assignment_matrix=assignment_matrix,
                )

                # Log exclusivity scores
                for m_idx, m_name in enumerate(marker_names):
                    if marker_owners[m_idx]:
                        logging.info(f"  Marker exclusivity: {m_name} = {marker_exclusivity[m_idx]:.3f}")
                self.results["marker_exclusivity"] = {
                    marker_names[i]: marker_exclusivity[i] for i in range(len(marker_names))
                }

            except ValueError as e:
                error_msg = f"Cell proportion validation failed for sample '{self.sample_name}': {str(e)}"
                logging.error(error_msg)
                raise ValueError(error_msg) from e

            # Optionally skip finetuning
            if skip_finetuning:
                logging.info("Skipping finetuning step (skip_finetuning=True)")
                global_cell_type_proportions_df = pd.DataFrame(Y_values, index=spot_names, columns=cell_type_names)
                finetuned_cell_type_proportions_df = global_cell_type_proportions_df.copy()
            else:
                # Create finetuning output directory
                finetune_output_dir = os.path.join(self.output_folder, "cell_prop_finetuning")
                os.makedirs(finetune_output_dir, exist_ok=True)

                if self.antibody_capture_adata is None:
                    raise ValueError("Antibody capture data has not been split. Run `split_adata` first.")

                Y_prev, beta_prev = finetune_cell_proportions_per_marker(
                    marker_level_data=marker_level_data,
                    marker_names=marker_names,
                    assignment_matrix=assignment_matrix,
                    cell_type_names=cell_type_names,
                    initial_Y_values=Y_values,
                    initial_beta_values=beta_values,
                    adata=self.antibody_capture_adata,
                    radius=radius,
                    tolerance=tolerance,
                    lambda_reg=lambda_reg,
                    alpha=alpha,
                    max_iterations=max_iterations,
                    max_y_change=max_y_change,
                    beta_vary=True,
                    beta_min=beta_min,
                    beta_max=beta_max,
                    marker_exclusivity=marker_exclusivity,
                    marker_alpha=alpha_values,
                    max_workers=max_workers,
                    checkpoint_interval=checkpoint_interval,
                    output_dir=finetune_output_dir,
                    rerun=True,
                )

                global_cell_type_proportions_df = pd.DataFrame(Y_values, index=spot_names, columns=cell_type_names)
                finetuned_cell_type_proportions_df = pd.DataFrame(Y_prev, index=spot_names, columns=cell_type_names)

        else:
            # ====== LEGACY PER-CELLTYPE BETA APPROACH ======
            logging.info("Using legacy per-celltype beta optimization (averages markers per cell type)")

            profile_based_antibody_data, cell_type_names = map_antibodies_to_profiles(
                self.antibody_capture_adata, self.cell_profile_dict
            )

            try:
                logging.info(f"Running Stage 1 cell proportion optimization with validation thresholds: "
                            f"Unknown<{unknown_threshold*100:.1f}%, CellTypes>{min_celltype_threshold*100:.1f}%, Redundancy<{redundancy_threshold*100:.0f}%")

                Y_values, beta_values = optimize_cell_proportions(
                    profile_based_antibody_data,
                    cell_type_names,
                    unknown_threshold=unknown_threshold,
                    min_celltype_threshold=min_celltype_threshold,
                    redundancy_threshold=redundancy_threshold,
                    warn_only=validation_warn_only,
                    lambda_laplacian=lambda_laplacian,
                    coords=coords,
                    laplacian_k=laplacian_k,
                )
            except ValueError as e:
                error_msg = f"Cell proportion validation failed for sample '{self.sample_name}': {str(e)}"
                logging.error(error_msg)
                raise ValueError(error_msg) from e

            # Optionally skip finetuning
            if skip_finetuning:
                logging.info("Skipping finetuning step (skip_finetuning=True)")
                global_cell_type_proportions_df = pd.DataFrame(Y_values, index=spot_names, columns=cell_type_names)
                finetuned_cell_type_proportions_df = global_cell_type_proportions_df.copy()
            else:
                # Create finetuning output directory
                finetune_output_dir = os.path.join(self.output_folder, "cell_prop_finetuning")
                os.makedirs(finetune_output_dir, exist_ok=True)

                if self.antibody_capture_adata is None:
                    raise ValueError("Antibody capture data has not been split. Run `split_adata` first.")

                Y_prev, beta_prev = finetune_cell_proportions(
                    profile_based_antibody_data,
                    cell_type_names,
                    Y_values,
                    beta_values,
                    self.antibody_capture_adata,
                    radius=radius,
                    max_workers=max_workers,
                    checkpoint_interval=checkpoint_interval,
                    output_dir=finetune_output_dir,
                    rerun=True,
                    beta_vary=True,
                    tolerance=tolerance,
                    max_iterations=max_iterations,
                    lambda_reg=lambda_reg,
                    alpha=alpha,
                    max_y_change=max_y_change,
                )

                global_cell_type_proportions_df = pd.DataFrame(Y_values, index=spot_names, columns=cell_type_names)
                finetuned_cell_type_proportions_df = pd.DataFrame(Y_prev, index=spot_names, columns=cell_type_names)

        global_cell_type_proportions_df = global_cell_type_proportions_df.sort_index()
        finetuned_cell_type_proportions_df = finetuned_cell_type_proportions_df.sort_index()

        global_cell_type_proportions_df.to_csv(
            os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_global_results.csv")
        )
        finetuned_cell_type_proportions_df.to_csv(
            os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_finetuned_results.csv")
        )

        # Store to self.results for use by run_cell_expression_pass1
        self.results["cell_prop"] = finetuned_cell_type_proportions_df

        return global_cell_type_proportions_df, finetuned_cell_type_proportions_df

    def _run_cell_classification(self, **kwargs):
        """Archived — gating-based cell classification module removed."""
        raise NotImplementedError(
            "Gating-based cell classification has been archived. "
            "Use spot-resolution QP deconvolution (use_gating=False) instead."
        )

    def run_cell_expression_pass1(
        self,
        radius=None,
        alpha=0.5,
        global_enrichment_weight=0.5,
        local_enrichment_weight=0.5,
        max_workers=None,
        checkpoint_interval=100,
        output_dir="checkpoints",
        rerun=True,
        continuous_relaxation=True,
        lambda_gex_reg=0.01,
    ):
        """
        Run first pass of gene expression deconvolution.

        Args:
            radius (float): Radius for neighbor detection
            alpha (float): Weight for spatial regularization
            global_enrichment_weight (float): Weight for global expression enrichment (0-1)
            local_enrichment_weight (float): Weight for local expression enrichment (0-1)
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
        if not self.preprocessed_gex:
            raise ValueError("Gene expression data not preprocessed. Run preprocess_gex() first.")

        # Auto-detect radius if not specified
        if radius is None:
            source_adata = self.gene_expression_adata or self.antibody_capture_adata
            if source_adata is None:
                raise ValueError("No AnnData available for radius auto-detection")
            radius = compute_optimal_radius(source_adata)
            logging.info(f"Auto-detected radius: {radius:.2f} (3 rings)")

        logging.info("Starting Pass 1: Error minimization with enrichment weights...")

        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        if "cell_prop" not in self.results or self.results["cell_prop"] is None:
            raise ValueError("Cell proportions not computed. Run cell proportion model first.")

        cell_props_values = self.results["cell_prop"].values

        # Diagnostic check for low or zero cell proportions
        total_props_per_spot = cell_props_values.sum(axis=1)
        zero_prop_spots = np.where(total_props_per_spot < 1e-9)[0]

        if len(zero_prop_spots) > 0:
            logging.warning(
                f"Found {len(zero_prop_spots)} spots with zero or negligible total cell proportions. "
                f"Deconvolution for these spots may fail and their profiles will be imputed as zeros. "
                f"Example spot indices: {zero_prop_spots[:5]}"
            )

        if self.resolution == "cell":
            # Cell-level: estimate true expression per cell (not deconvolution)
            logging.info("Cell-level mode: using true count estimation instead of deconvolution")

            # Get gene expression data
            gex_data = self.gene_expression_adata.X
            if hasattr(gex_data, 'toarray'):
                gex_data = gex_data.toarray()
            gex_data = np.asarray(gex_data, dtype=np.float64)

            # Get spatial coordinates
            cell_coords = self.gene_expression_adata.obsm.get('spatial', None)
            if cell_coords is None:
                cell_coords = self.antibody_capture_adata.obsm.get('spatial', None)
            if cell_coords is None:
                raise ValueError("No spatial coordinates found for cell-level expression estimation")

            # Build enrichment weights from cell type assignments
            n_types = cell_props_values.shape[1]
            n_genes = gex_data.shape[1]
            dominant_type = np.argmax(cell_props_values, axis=1)

            enrichment = np.ones((n_types, n_genes)) * 0.1
            for ct_idx in range(n_types):
                ct_cells = np.where(dominant_type == ct_idx)[0]
                if len(ct_cells) > 0:
                    ct_mean = gex_data[ct_cells].mean(axis=0)
                    global_mean = gex_data.mean(axis=0) + 1e-10
                    enrichment[ct_idx] = ct_mean / global_mean

            X_true = estimate_true_expression_cell(
                X_obs=gex_data,
                Y_assignments=cell_props_values,
                coords=cell_coords,
                enrichment_weights=enrichment,
                library_slack=self.resolution_params["pass2_library_slack"],
                lambda_spatial=self.resolution_params["lambda_spatial"],
                spatial_k=self.resolution_params["neighbor_k"],
                max_workers=max_workers,
            )

            # Store results in same format as spot-level for downstream compatibility
            N = X_true.shape[0]
            T = n_types
            M = n_genes
            spotwise_gene_expression_profiles = {}
            for i in range(N):
                # For cell-level: all expression assigned to dominant type
                profile = np.zeros((T, M))
                dt = dominant_type[i]
                profile[dt, :] = X_true[i, :]
                spotwise_gene_expression_profiles[i] = profile

            self.results["gene_expression"] = spotwise_gene_expression_profiles
            logging.info(f"Cell-level expression estimation complete for {N} cells")
            return

        spotwise_profiles = optimize_gene_expression(
            sample_name=self.sample_name,
            deconvolution_expression_data=self.gene_expression_adata.X,
            cell_type_numbers_array=cell_props_values,
            filtered_adata=self.gene_expression_adata,
            radius=radius,
            global_enrichment_weight=global_enrichment_weight,
            local_enrichment_weight=local_enrichment_weight,
            global_prior=None,  # No prior in pass 1
            lambda_prior_weight=0.0,  # No prior weight in pass 1
            max_workers=max_workers,
            checkpoint_interval=checkpoint_interval,
            output_dir=output_dir,
            rerun=rerun,
            continuous_relaxation=continuous_relaxation,
            lambda_gex_reg=lambda_gex_reg,
        )

        # Get dimensions for NaN imputation and consistency checks
        N = self.gene_expression_adata.shape[0]  # number of spots
        T = cell_props_values.shape[1]  # number of cell types
        M = self.gene_expression_adata.shape[1]  # number of genes
        dimensions = (N, T, M)

        # Impute spots that failed to converge (NaN spots)
        nan_spots = [i for i in range(N) if i not in spotwise_profiles]
        if nan_spots:
            logging.info(f"Found {len(nan_spots)} spots that failed to converge. Starting imputation...")

            imputed_count = 0
            zero_profile = np.zeros((T, M), dtype=float)

            for spot_idx in nan_spots:
                # Prioritize imputing with zeros if cell proportions are negligible
                if total_props_per_spot[spot_idx] < 1e-9:
                    spotwise_profiles[spot_idx] = zero_profile
                    logging.info(f"Imputed spot {spot_idx} with a zero profile due to negligible cell proportions.")
                    imputed_count += 1
                    continue

                # Otherwise, use neighbor-based imputation
                neighbor_indices = get_neighbors_with_fixed_radius(
                    spot_idx, self.gene_expression_adata, radius=radius, include_center=False
                )

                # Corrected key usage: use integer keys consistently
                neighbor_profiles = [spotwise_profiles[i] for i in neighbor_indices if i in spotwise_profiles]

                if neighbor_profiles:
                    imputed_profile = np.nanmean(neighbor_profiles, axis=0)
                    spotwise_profiles[spot_idx] = imputed_profile
                    logging.info(f"Imputed spot {spot_idx} using {len(neighbor_profiles)} neighbors.")
                    imputed_count += 1
                else:
                    logging.warning(
                        f"No valid neighbors found to impute spot {spot_idx}. It will be filled with zeros as a fallback."
                    )
                    spotwise_profiles[spot_idx] = zero_profile  # Fallback to prevent downstream errors
                    imputed_count += 1

            logging.info(f"Finished imputation. {imputed_count}/{len(nan_spots)} failed spots were imputed.")

        # Final check for any remaining NaNs, though the logic above should prevent this
        final_nan_spots = [i for i in range(N) if i not in spotwise_profiles]
        if final_nan_spots:
            logging.error(
                f"FATAL: {len(final_nan_spots)} spots could not be imputed: {final_nan_spots[:10]}. Check imputation logic."
            )
            # This case should ideally not be reached. Raising an error might be appropriate.
            raise RuntimeError("Failed to impute all necessary spot profiles.")

        # Store first pass results
        self.results["gene_expression_pass1"] = spotwise_profiles

        # Save and evaluate results
        parquet_path = os.path.join(self.output_folder, f"{self.sample_name}_gene_expression_pass1.parquet")
        self._save_profiles_to_parquet(spotwise_profiles, parquet_path)

        self.append_gex_to_adata(pass_number=1)

        layer_dir = os.path.join(self.output_folder, f"{self.sample_name}_pass1/layers")
        export_anndata_layers(self.gene_expression_adata, layer_dir, pass_number=1)

        return {"spotwise_profiles": spotwise_profiles, "dimensions": dimensions}

    def compute_expression_prior(
        self,
        spotwise_profiles_pass1: Dict[int, np.ndarray],
        cell_type_numbers_array: np.ndarray,
        lambda_prior: float = 1.0,
        min_expression_threshold: float = 0.1,
    ) -> Dict[str, Any]:
        """
        Compute global prior from pass 1 results.

        Args:
            spotwise_profiles_pass1: Dictionary mapping spot indices to profile matrices
            cell_type_numbers_array: Array of cell type proportions (N_spots × T_celltypes)
            lambda_prior: Strength of prior influence (default: 1.0)
            min_expression_threshold: Minimum expression to consider "active" (default: 0.1)

        Returns:
            Dict[str, Any]: {
                'global_prior': np.ndarray,  # shape (T_celltypes, M_genes)
                'confidence_scores': np.ndarray,  # shape (T_celltypes, M_genes)
                'expression_patterns': Dict containing detailed statistics
            }
        """
        if not self.preprocessed_gex:
            raise ValueError("Gene expression data not preprocessed. Run preprocess_gex() first.")

        logging.info("Computing prior from pass 1 results...")

        # Get gene and cell type names for validation
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data not available")
        gene_names = self.gene_expression_adata.var_names

        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary not loaded. Run load_cell_profile_dict() first.")
        cell_type_names = list(self.cell_profile_dict.keys())

        # Compute global prior with new approach
        prior_info = compute_global_prior(
            spotwise_profiles_pass1,
            cell_type_numbers_array,
            lambda_prior=lambda_prior,
            min_expression_threshold=min_expression_threshold,
        )

        # Validate prior shape
        T = cell_type_numbers_array.shape[1]  # num cell types
        M = self.gene_expression_adata.shape[1]  # num genes

        if prior_info["global_prior"].shape != (T, M):
            raise ValueError(f"Prior shape {prior_info['global_prior'].shape} does not match expected ({T}, {M})")

        # Log detailed statistics about the prior
        logging.info("\nPrior computation details:")
        logging.info(f"Number of cell types: {T}")
        logging.info(f"Number of genes: {M}")

        # Per cell-type statistics
        for t, cell_type in enumerate(cell_type_names):
            mean_conf = np.mean(prior_info["confidence_scores"][t])
            strong_signals = np.mean(prior_info["global_prior"][t] > 0.5)
            logging.info(f"\n{cell_type}:")
            logging.info(f" - Mean confidence score: {mean_conf:.4f}")
            logging.info(f" - % Strong signals: {100 * strong_signals:.2f}%")

            # Expression pattern summary
            mean_exp = np.mean(prior_info["expression_patterns"]["mean_expression"][t])
            freq = np.mean(prior_info["expression_patterns"]["expression_frequency"][t])
            cons = np.mean(prior_info["expression_patterns"]["expression_consistency"][t])
            logging.info(f" - Mean expression: {mean_exp:.4f}")
            logging.info(f" - Mean expression frequency: {freq:.4f}")
            logging.info(f" - Mean expression consistency: {cons:.4f}")

        return prior_info

    def _save_profiles_to_parquet(self, profiles, path):
        """Helper method to save profiles to parquet format with consistent naming."""
        if not profiles:
            logging.warning("No profiles to save.")
            return

        N = max(profiles.keys()) + 1
        T = profiles[0].shape[0]
        M = profiles[0].shape[1]

        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary not loaded. Run load_cell_profile_dict() first.")

        # Get cell type names from the dictionary
        cell_type_names = list(self.cell_profile_dict.keys())

        # Create combined matrix with proper cell type names and spot formatting
        spot_celltype_indices = []

        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data not available")

        # Use ::: as delimiter to avoid conflicts with underscores in cell type names
        # (e.g., auto-discovered profiles like "Cancer Epithelial_Protein_1")
        for i in range(N):
            spot_name = self.gene_expression_adata.obs_names[i]  # Use actual spot names from AnnData
            for cell_type in cell_type_names:
                spot_celltype_indices.append(f"{spot_name}:::{cell_type}")

        gene_names = self.gene_expression_adata.var_names
        nan_matrix = np.full((T, M), np.nan)
        data_combined = np.vstack([profiles.get(i, nan_matrix) for i in range(N)])

        # Create DataFrame
        df = pd.DataFrame(data_combined, index=spot_celltype_indices, columns=gene_names)

        df.to_parquet(path, compression="gzip")
        logging.info(f"Saved profiles to {path} with cell types: {cell_type_names}")

    def append_proportions_to_adata(self, proportions_path=None, key="finetuned"):
        """Append cell type proportions to AnnData object."""
        if proportions_path is None:
            proportions_path = os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_{key}_results.csv")

        # Load proportions CSV
        spot_by_celltype_df = pd.read_csv(proportions_path, index_col=0)

        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data not available")

        # Debug prints before sorting
        print("\nBefore sorting:")
        print("CSV spots 1-10:", list(spot_by_celltype_df.index[:10]))
        print("AnnData spots 1-10:", list(self.gene_expression_adata.obs_names[:10]))

        if "spot_" in str(spot_by_celltype_df.index[0]):
            # Sort both numerically by the spot number
            def get_spot_number(x):
                return int(x.split("spot_")[1])

            # Sort using reindex instead of sort_index
            sorted_csv_idx = sorted(spot_by_celltype_df.index, key=get_spot_number)
            sorted_adata_idx = sorted(self.gene_expression_adata.obs_names, key=get_spot_number)

            spot_by_celltype_df = spot_by_celltype_df.reindex(sorted_csv_idx)
            self.gene_expression_adata = self.gene_expression_adata[sorted_adata_idx].copy()

            # Debug prints after sorting
            print("\nAfter sorting:")
            print("CSV spots 1-10:", list(spot_by_celltype_df.index[:10]))
            print("AnnData spots 1-10:", list(self.gene_expression_adata.obs_names[:10]))

        # Check if indices match after sorting
        if not all(spot_by_celltype_df.index == self.gene_expression_adata.obs_names):
            raise ValueError("Spot indices still don't match after sorting. Please verify your data.")

        # Add cell type proportions to adata.obs
        for cell_type in spot_by_celltype_df.columns:
            self.gene_expression_adata.obs[cell_type] = spot_by_celltype_df[cell_type]

        self.results["cell_prop"] = spot_by_celltype_df

        print("✅ Cell type proportions have been appended to adata.obs and results['cell_prop']")

    def append_gex_to_adata(self, parquet_path=None, pass_number=1):
        """
        Append gene expression layers from a Parquet file back into the gene_expression_adata object.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        if parquet_path is None:
            parquet_path = os.path.join(
                self.output_folder, f"{self.sample_name}_gene_expression_pass{pass_number}.parquet"
            )

        # Step 1: Read the Parquet file into a pandas DataFrame
        table = pq.read_table(parquet_path)
        df = table.to_pandas()
        print(f"Parquet file for pass {pass_number} loaded successfully.")

        # Step 2: Reset the index to extract 'Spot' and 'CellType'
        df = df.reset_index()

        # Check which delimiter format is used (new ::: or legacy _)
        sample_index = df["index"].iloc[0]
        if ":::" in sample_index:
            # New format: spot_1:::CellType
            df[["Spot", "CellType"]] = df["index"].str.split(":::", n=1, expand=True)
            print("Spot and CellType successfully split (using ::: delimiter).")
        else:
            # Legacy format: spot_1_CellType - need to match against known cell types
            # Get cell type names from the dictionary for parsing
            if self.cell_profile_dict is None:
                raise ValueError("Cell profile dictionary not loaded. Run load_cell_profile_dict() first.")
            known_cell_types = list(self.cell_profile_dict.keys())

            # Parse each index by finding which cell type suffix matches
            spots = []
            cell_types = []
            for idx_val in df["index"]:
                matched = False
                for ct in known_cell_types:
                    suffix = f"_{ct}"
                    if idx_val.endswith(suffix):
                        spots.append(idx_val[:-len(suffix)])
                        cell_types.append(ct)
                        matched = True
                        break
                if not matched:
                    raise ValueError(f"Could not parse index '{idx_val}' - no matching cell type found in {known_cell_types}")

            df["Spot"] = spots
            df["CellType"] = cell_types
            print("Spot and CellType successfully split (using legacy _ delimiter with cell type matching).")

        df = df.drop(columns=["index"])

        # Debug print spot names
        print("\nSpot name formats:")
        print("AnnData spot names format:", self.gene_expression_adata.obs_names[:5])
        print("Parquet spot names format:", df["Spot"].unique()[:5])

        # Get cell type names from the dictionary for validation
        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary not loaded. Run load_cell_profile_dict() first.")

        expected_cell_types = set(self.cell_profile_dict.keys())
        found_cell_types = set(df["CellType"].unique())

        if not found_cell_types.issubset(expected_cell_types):
            logging.warning(f"Found unexpected cell types: {found_cell_types - expected_cell_types}")
            logging.warning(f"Expected cell types: {expected_cell_types}")
            raise ValueError("Cell type mismatch in loaded data")

        # Step 3: Process each cell type
        for cell_type in found_cell_types:
            # Filter data for this cell type
            celltype_data = df[df["CellType"] == cell_type].copy()
            celltype_data = celltype_data.drop(columns=["CellType"])

            # Ensure spot names match AnnData format
            if (
                "spot_" in str(self.gene_expression_adata.obs_names[0])
                and not celltype_data["Spot"].str.contains("spot_").all()
            ):
                celltype_data["Spot"] = "spot_" + celltype_data["Spot"].astype(str)
            elif celltype_data["Spot"].str.contains("spot_").all() and not "spot_" in str(
                self.gene_expression_adata.obs_names[0]
            ):
                celltype_data["Spot"] = celltype_data["Spot"].str.replace("spot_", "")

            # Set Spot as index
            celltype_data = celltype_data.set_index("Spot")

            # Verify all spots exist in AnnData
            missing_spots = set(celltype_data.index) - set(self.gene_expression_adata.obs_names)
            if missing_spots:
                raise ValueError(f"Found spots in parquet that don't exist in AnnData: {missing_spots}")

            # Create matrix with proper spot ordering
            celltype_matrix = np.zeros(
                (len(self.gene_expression_adata.obs_names), len(self.gene_expression_adata.var_names))
            )
            for spot in self.gene_expression_adata.obs_names:
                if spot in celltype_data.index:
                    idx = self.gene_expression_adata.obs_names.get_loc(spot)
                    celltype_matrix[idx] = celltype_data.loc[spot].values

            # Add as layer with consistent naming
            layer_name = f"{cell_type.replace(' ', '_')}_genes_pass{pass_number}"
            self.gene_expression_adata.layers[layer_name] = celltype_matrix
            print(f"Added layer: {layer_name} (Shape: {celltype_matrix.shape})")

            # After adding each layer, verify it was added correctly
            if layer_name not in self.gene_expression_adata.layers:
                logging.error(f"Failed to add layer: {layer_name}")
            else:
                logging.info(f"Successfully added layer: {layer_name}")

    def get_adata(self):
        """
        Retrieve the internal AnnData object for downstream analysis.

        Returns:
            AnnData: The internal `adata` object.
        """
        if self.gene_expression_adata is None:
            raise ValueError("AnnData object is not initialized in the model.")

        print("✅ Returning the internal AnnData object.")
        return self.gene_expression_adata

    def cleanup(self):
        """Free memory and clean up temporary data."""
        cleanup_memory()

    def validate_neighborhood_size(self, radius):
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dict has not been loaded. Run 'load_cell_profile_dict' first.")
        assert_neighborhood_size(self.gene_expression_adata, self.cell_profile_dict, radius=radius, num_spots=5)
