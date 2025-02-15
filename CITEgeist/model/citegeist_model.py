# Standard library imports
import os
import logging
from typing import Dict, Any, Optional, Tuple, List, Union


# Third-party imports
import numpy as np
import pandas as pd
import scanpy as sc
import pyarrow.parquet as pq
from scipy.ndimage import gaussian_filter

# Local imports
from .gurobi_impl import (
    map_antibodies_to_profiles, 
    optimize_cell_proportions, 
    optimize_gene_expression, 
    compute_global_prior,
    normalize_counts,
    finetune_cell_proportions
)
from .utils import (
    validate_cell_profile_dict,  
    cleanup_memory, 
    setup_logging, 
    get_neighbors_with_fixed_radius, 
    assert_neighborhood_size,
    export_anndata_layers,
)


class CitegeistModel:
    def __init__(self, sample_name, adata=None, output_folder=None, simulation=False, gene_expression_adata=None, antibody_capture_adata=None):
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
        
        if 'feature_types' not in self.adata.var.columns:
            raise ValueError("The 'feature_types' column is missing in `adata.var`. Cannot split data.")
        
        if self.adata is None:
            raise ValueError("No valid data loaded. Ensure `adata` or split datasets are loaded properly.")

        self.adata.var_names_make_unique()
        
        if self.gene_expression_adata or self.antibody_capture_adata :
                raise ValueError(
                    "Data seems to already be split"
                )

        # Identify indices for Gene Expression and Antibody Capture
        gene_expression_idx = np.where(self.adata.var['feature_types'] == 'Gene Expression')[0]
        antibody_capture_idx = np.where(self.adata.var['feature_types'] == 'Antibody Capture')[0]

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
        matrix = self.gene_expression_adata.X.toarray() if hasattr(self.gene_expression_adata.X, 'toarray') else self.gene_expression_adata.X
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
        

        print(f"Filtered gene expression data: {initial_gene_count} → {filtered_gene_count} genes "
              f"(count > 0 in at least {nonzero_percentage*100}% of spots, mean expression > {mean_expression_threshold} "
              f"in nonzero spots). Remaining spots: {self.gene_expression_adata.shape[0]} "
              f"Filtered spots: {initial_spot_count} to {self.gene_expression_adata.shape[0]}")

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
        if hasattr(matrix, 'toarray'):
            matrix = matrix.toarray()
        
        if not np.all(np.equal(np.mod(matrix, 1), 0)):
            raise ValueError("Gene expression data contains non-integer values after normalization.")
        
        self.preprocessed_gex = True
        logging.info(f"Gene expression data normalized to {target_sum} counts per spot and validated for discrete count analysis.")

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
        matrix = self.antibody_capture_adata.X.toarray() if hasattr(self.antibody_capture_adata.X, 'toarray') else self.antibody_capture_adata.X
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

    def run_cell_proportion_model(self, radius=None, tolerance=1e-4, max_iterations=20, lambda_reg=1, alpha=0.5, max_y_change=0.4, max_workers=None, checkpoint_interval=100):
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
        """

        if radius is None:
            raise ValueError("Radius must be provided. Run `run_cell_proportion_model` with a radius argument.")

        if self.adata is None and (self.gene_expression_adata is None or self.antibody_capture_adata is None):
            raise ValueError("No valid data loaded. Ensure `adata` or split datasets are loaded properly.")

        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary has not been loaded. Run `load_cell_profile_dict` first.")

        profile_based_antibody_data, cell_type_names = map_antibodies_to_profiles(self.antibody_capture_adata, self.cell_profile_dict)
        
        Y_values, beta_values = optimize_cell_proportions(profile_based_antibody_data, cell_type_names)
        
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
            max_y_change=max_y_change
        )

        spot_names = self.antibody_capture_adata.obs_names
        
        global_cell_type_proportions_df = pd.DataFrame(Y_values, index=spot_names, columns=cell_type_names)
        finetuned_cell_type_proportions_df = pd.DataFrame(Y_prev, index=spot_names, columns=cell_type_names)

        global_cell_type_proportions_df = global_cell_type_proportions_df.sort_index()
        finetuned_cell_type_proportions_df = finetuned_cell_type_proportions_df.sort_index()
        
        global_cell_type_proportions_df.to_csv(os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_global_results.csv"))
        finetuned_cell_type_proportions_df.to_csv(os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_finetuned_results.csv"))

        return global_cell_type_proportions_df, finetuned_cell_type_proportions_df

    def run_cell_expression_pass1(self, radius, alpha=0.5, lambda_reg_gex=0.001,
                            global_enrichment_weight=0.5, local_enrichment_weight=0.5,
                            max_workers=None, checkpoint_interval=100, 
                            output_dir="checkpoints", rerun=True):
        """
        Run first pass of gene expression deconvolution.
        
        Args:
            radius (float): Radius for neighbor detection
            alpha (float): Weight for spatial regularization
            lambda_reg_gex (float): Weight for gene expression regularization
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

        logging.info("Starting Pass 1: Error minimization with enrichment weights...")
        
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")
        
        spotwise_profiles = optimize_gene_expression(
            sample_name=self.sample_name,
            deconvolution_expression_data=self.gene_expression_adata.X,
            cell_type_numbers_array=self.results['cell_prop'].values,
            filtered_adata=self.gene_expression_adata,
            radius=radius,
            global_enrichment_weight=global_enrichment_weight,
            local_enrichment_weight=local_enrichment_weight,
            global_prior=None,  # No prior in pass 1
            lambda_prior_weight=0.0,  # No prior weight in pass 1
            max_workers=max_workers,
            checkpoint_interval=checkpoint_interval,
            output_dir=output_dir,
            rerun=rerun
        )

       
       
        # Get dimensions for NaN imputation
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data not available")
        if 'cell_prop' not in self.results or self.results['cell_prop'] is None:
            raise ValueError("Cell proportions not computed. Run cell proportion model first.")
            
        N = self.gene_expression_adata.shape[0]  # number of spots
        T = self.results['cell_prop'].values.shape[1]  # number of cell types
        M = self.gene_expression_adata.shape[1]  # number of genes
        
        # Impute NaN spots for first pass
        nan_spots = [i for i in range(N) if i not in spotwise_profiles]
        for nan_spot in nan_spots:
            neighbor_indices = get_neighbors_with_fixed_radius(
                nan_spot, 
                self.gene_expression_adata, 
                radius=radius, 
                include_center=False
            )

            if spotwise_profiles is None:
                raise ValueError("Spotwise profiles not computed. Run cell expression pass 1 first.")

            neighbor_profiles = [
                spotwise_profiles[str(i)]
                for i in neighbor_indices
                if i in spotwise_profiles
            ]

            if neighbor_profiles:
                # Round to nearest integer for count data
                imputed_profile = np.round(np.nanmean(neighbor_profiles, axis=0)).astype(int)
                spotwise_profiles[str(nan_spot)] = imputed_profile
                logging.info(f"Imputed spot {nan_spot} using neighbors at radius {radius} (Pass 1).")
            else:
                logging.warning(f"No valid neighbors found to impute spot {nan_spot} (Pass 1). Leaving as NaN.")

        # Store first pass results
        self.results['gene_expression_pass1'] = spotwise_profiles
        
        # Get dimensions for consistency checks
        N = self.gene_expression_adata.shape[0]  # spots
        T = self.results['cell_prop'].values.shape[1]  # cell types
        M = self.gene_expression_adata.shape[1]  # genes
        dimensions = (N, T, M)

        # Save and evaluate results
        parquet_path = os.path.join(self.output_folder, f"{self.sample_name}_gene_expression_pass1.parquet")
        self._save_profiles_to_parquet(spotwise_profiles, parquet_path)
        
        self.append_gex_to_adata(pass_number=1)
        
        layer_dir = os.path.join(self.output_folder, f"{self.sample_name}_pass1/layers")
        export_anndata_layers(self.gene_expression_adata, layer_dir, pass_number=1)

        return {
            'spotwise_profiles': spotwise_profiles,
            'dimensions': dimensions
        }

    def compute_expression_prior(
        self, 
        spotwise_profiles_pass1: Dict[int, np.ndarray],
        cell_type_numbers_array: np.ndarray,
        lambda_prior: float = 1.0,
        min_expression_threshold: float = 0.1
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
            min_expression_threshold=min_expression_threshold
        )
        
        # Validate prior shape
        T = cell_type_numbers_array.shape[1]  # num cell types
        M = self.gene_expression_adata.shape[1]  # num genes
        
        if prior_info['global_prior'].shape != (T, M):
            raise ValueError(f"Prior shape {prior_info['global_prior'].shape} does not match expected ({T}, {M})")
        
        # Log detailed statistics about the prior
        logging.info("\nPrior computation details:")
        logging.info(f"Number of cell types: {T}")
        logging.info(f"Number of genes: {M}")
        
        # Per cell-type statistics
        for t, cell_type in enumerate(cell_type_names):
            mean_conf = np.mean(prior_info['confidence_scores'][t])
            strong_signals = np.mean(prior_info['global_prior'][t] > 0.5)
            logging.info(f"\n{cell_type}:")
            logging.info(f" - Mean confidence score: {mean_conf:.4f}")
            logging.info(f" - % Strong signals: {100 * strong_signals:.2f}%")
            
            # Expression pattern summary
            mean_exp = np.mean(prior_info['expression_patterns']['mean_expression'][t])
            freq = np.mean(prior_info['expression_patterns']['expression_frequency'][t])
            cons = np.mean(prior_info['expression_patterns']['expression_consistency'][t])
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


        for i in range(N):
            spot_name = self.gene_expression_adata.obs_names[i]  # Use actual spot names from AnnData
            for cell_type in cell_type_names:
                spot_celltype_indices.append(f"{spot_name}_{cell_type}")
        
        gene_names = self.gene_expression_adata.var_names
        nan_matrix = np.full((T, M), np.nan)
        data_combined = np.vstack([
            profiles.get(i, nan_matrix) 
            for i in range(N)
        ])
        
        # Create DataFrame
        df = pd.DataFrame(
            data_combined, 
            index=spot_celltype_indices, 
            columns=gene_names
        )
        
        df.to_parquet(path, compression="gzip")
        logging.info(f"Saved profiles to {path} with cell types: {cell_type_names}")

    def append_proportions_to_adata(self, proportions_path=None, key='finetuned'):
        """Append cell type proportions to AnnData object."""
        if proportions_path is None:
            proportions_path = os.path.join(self.output_folder, f'{self.sample_name}_cell_prop_{key}_results.csv')

        # Load proportions CSV
        spot_by_celltype_df = pd.read_csv(proportions_path, index_col=0)


        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data not available")
        
        # Debug prints before sorting
        print("\nBefore sorting:")
        print("CSV spots 1-10:", list(spot_by_celltype_df.index[:10]))
        print("AnnData spots 1-10:", list(self.gene_expression_adata.obs_names[:10]))
        
        if 'spot_' in str(spot_by_celltype_df.index[0]):
            # Sort both numerically by the spot number
            def get_spot_number(x):
                return int(x.split('spot_')[1])
            
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

        self.results['cell_prop'] = spot_by_celltype_df
        
        print("✅ Cell type proportions have been appended to adata.obs and results['cell_prop']")
        
        
    def append_gex_to_adata(self, parquet_path=None, pass_number=1):
        """
        Append gene expression layers from a Parquet file back into the gene_expression_adata object.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")
        
        if parquet_path is None:
            parquet_path = os.path.join(
                self.output_folder, 
                f"{self.sample_name}_gene_expression_pass{pass_number}.parquet"
            )

        # Step 1: Read the Parquet file into a pandas DataFrame
        table = pq.read_table(parquet_path)
        df = table.to_pandas()
        print(f"Parquet file for pass {pass_number} loaded successfully.")

        # Step 2: Reset the index to extract 'Spot' and 'CellType'
        df = df.reset_index()
        df[['Spot', 'CellType']] = df['index'].str.rsplit('_', n=1, expand=True)
        df = df.drop(columns=['index'])
        print("Spot and CellType successfully split.")

        # Debug print spot names
        print("\nSpot name formats:")
        print("AnnData spot names format:", self.gene_expression_adata.obs_names[:5])
        print("Parquet spot names format:", df['Spot'].unique()[:5])

        # Get cell type names from the dictionary for validation
        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary not loaded. Run load_cell_profile_dict() first.")
        
        expected_cell_types = set(self.cell_profile_dict.keys())
        found_cell_types = set(df['CellType'].unique())
        
        if not found_cell_types.issubset(expected_cell_types):
            logging.warning(f"Found unexpected cell types: {found_cell_types - expected_cell_types}")
            logging.warning(f"Expected cell types: {expected_cell_types}")
            raise ValueError("Cell type mismatch in loaded data")

        # Step 3: Process each cell type
        for cell_type in found_cell_types:
            # Filter data for this cell type
            celltype_data = df[df['CellType'] == cell_type].copy()
            celltype_data = celltype_data.drop(columns=['CellType'])
            
            # Ensure spot names match AnnData format
            if 'spot_' in str(self.gene_expression_adata.obs_names[0]) and not celltype_data['Spot'].str.contains('spot_').all():
                celltype_data['Spot'] = 'spot_' + celltype_data['Spot'].astype(str)
            elif celltype_data['Spot'].str.contains('spot_').all() and not 'spot_' in str(self.gene_expression_adata.obs_names[0]):
                celltype_data['Spot'] = celltype_data['Spot'].str.replace('spot_', '')
            
            # Set Spot as index
            celltype_data = celltype_data.set_index('Spot')
            
            # Verify all spots exist in AnnData
            missing_spots = set(celltype_data.index) - set(self.gene_expression_adata.obs_names)
            if missing_spots:
                raise ValueError(f"Found spots in parquet that don't exist in AnnData: {missing_spots}")
            
            # Create matrix with proper spot ordering
            celltype_matrix = np.zeros((len(self.gene_expression_adata.obs_names), len(self.gene_expression_adata.var_names)))
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
        