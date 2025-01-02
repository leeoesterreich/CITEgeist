import os
import scanpy as sc
from .gurobi_impl import map_antibodies_to_profiles, optimize_cell_proportions, optimize_gene_expression
from .utils import validate_cell_profile_dict, save_results_to_output, cleanup_memory, setup_logging, get_neighbors_with_fixed_radius, assert_neighborhood_size
import numpy as np
import pandas as pd
import logging
import pyarrow.parquet as pq


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
            
        self.output_folder = output_folder
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
        if 'feature_types' not in self.adata.var.columns:
            raise ValueError("The 'feature_types' column is missing in `adata.var`. Cannot split data.")
        
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
    
    def filter_gex(self, row_cutoff=10, value_cutoff=1):
        """
        Filter genes in the gene expression AnnData object based on user-defined criteria.

        Args:
            row_cutoff (int): Minimum number of spots where a gene must meet the value threshold.
            value_cutoff (float): Minimum value threshold for gene expression.

        Raises:
            ValueError: If gene expression data has not been split or initialized.

        Returns:
            None: Updates `self.gene_expression_adata` in place.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        # Extract the data matrix
        matrix = self.gene_expression_adata.X.toarray() if hasattr(self.gene_expression_adata.X, 'toarray') else self.gene_expression_adata.X
        matrix = np.asarray(matrix)  # Ensure dense matrix

        # Compute boolean filter: Genes with at least `row_cutoff` spots having `value_cutoff` expression
        col_filter = (matrix >= value_cutoff).sum(axis=0) >= row_cutoff

        # Convert to 1D array if it's sparse
        col_filter = col_filter.A1 if hasattr(col_filter, "A1") else col_filter

        # Apply the filter and subset the AnnData object
        filtered_gene_count = np.sum(col_filter)
        initial_gene_count = self.gene_expression_adata.shape[1]

        self.gene_expression_adata = self.gene_expression_adata[:, col_filter].copy()

        print(f"Filtered gene expression data: {initial_gene_count} → {filtered_gene_count} genes "
              f"(row_cutoff={row_cutoff}, value_cutoff={value_cutoff}).")
    
    def preprocess_gex(self):
        """
        Preprocess gene expression data:
        - Normalize to target sum.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        sc.pp.normalize_total(self.gene_expression_adata, target_sum=1e4)
        
        self.preprocessed_gex = True

        print("Gene expression data preprocessing completed: Normalized.")

    def preprocess_antibody(self):
        """
        Preprocess antibody capture data:
        - Winsorize extreme values.
        - Apply global CLR normalization.
        - Raise an error if NaNs or Infs are detected in the processed data.
        """
        if self.antibody_capture_adata is None:
            raise ValueError("Antibody capture data has not been split. Run `split_adata` first.")

        # Step 1: Extract and ensure matrix is dense
        matrix = self.antibody_capture_adata.X.toarray() if hasattr(self.antibody_capture_adata.X, 'toarray') else self.antibody_capture_adata.X
        matrix = np.asarray(matrix)  # Ensure it's a NumPy array

        # Step 2: Validate initial data (no NaNs or Infs at start)
        if np.isnan(matrix).any() or np.isinf(matrix).any():
            raise ValueError("Antibody capture matrix contains NaN or Inf values before preprocessing.")

        # Step 3: Winsorize to cap extreme values
        matrix = self.winsorize(matrix, lower_percentile=5, upper_percentile=95)

        # Step 4: Prevent zero-only columns (add small epsilon if column sum is zero)
        column_sums = matrix.sum(axis=0)
        zero_columns = column_sums == 0
        if np.any(zero_columns):
            matrix[:, zero_columns] += 1e-6  # Add epsilon to zero columns

        # Step 5: Apply CLR Normalization
        matrix = self.global_clr(matrix)

        # Step 6: Final Validation for NaNs or Infs
        if np.isnan(matrix).any() or np.isinf(matrix).any():
            raise ValueError("NaN or Inf values detected in antibody capture matrix after preprocessing. "
                             "Check input data and preprocessing steps.")

        # Step 7: Reassign processed matrix to AnnData object
        self.antibody_capture_adata.X = matrix

        # Update status flag
        self.preprocessed_antibody = True

        print("Antibody capture data preprocessing completed: Winsorized, CLR applied, no NaNs detected.")

    def run_cell_proportion_model(self, tolerance=1e-4, max_iterations=50, lambda_reg=1, alpha=0.5):
            """
            Orchestrates the cell proportion optimization workflow.
            Delegates optimization to `optimize_cell_proportions` in `gurobi_impl.py`.
            """
            if self.adata is None and (self.gene_expression_adata is None or self.antibody_capture_adata is None):
                raise ValueError("No valid data loaded. Ensure `adata` or split datasets are loaded properly.")

            adata = self.antibody_capture_adata if self.antibody_capture_adata is not None else self.adata

            if self.cell_profile_dict is None:
                raise ValueError("Cell profile dictionary has not been loaded. Run `load_cell_profile_dict` first.")

            profile_based_antibody_data, cell_type_names = map_antibodies_to_profiles(self.antibody_capture_adata, self.cell_profile_dict)
            Y_values = optimize_cell_proportions(profile_based_antibody_data, cell_type_names)
            spot_names = self.adata.obs_names
            cell_type_proportions_df = pd.DataFrame(Y_values, index=spot_names, columns=cell_type_names)
            # Store and save results
            self.results['cell_prop'] = cell_type_proportions_df
            output_file = os.path.join(self.output_folder, f'{self.sample_name}_cell_prop_results.csv')
            save_results_to_output(cell_type_proportions_df, output_file)
            
            print(f"Cell type proportions saved to '{output_file}'.")
            
    def run_cell_expression_model(self, radius=2, alpha=0.5, lambda_reg_gex=0.0001, max_workers=16):
        """
        Run the gene expression deconvolution model using neighborhood-based Gurobi optimization.

        Args:
            radius (int): Neighborhood radius for neighboring spot selection.
            alpha (float): L1-L2 Elastic Net regularization parameter.
            lambda_reg_gex (float): Regularization strength for gene expression model.
            max_workers (int): Number of parallel workers.

        Raises:
            ValueError: If gene expression data is not initialized.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        self.validate_neighborhood_size(radius)
        
        # Extract required data
        filtered_data = self.gene_expression_adata.X.toarray() if hasattr(self.gene_expression_adata.X, 'toarray') else self.gene_expression_adata.X
        deconvolution_expression_data = filtered_data
        cell_type_numbers_array = self.results.get('cell_prop').values
        filtered_adata = self.gene_expression_adata
        # Run optimization
        spotwise_gene_expression_profiles = optimize_gene_expression(
            deconvolution_expression_data=deconvolution_expression_data,
            cell_type_numbers_array=cell_type_numbers_array,
            filtered_adata=self.gene_expression_adata,
            radius=radius,
            alpha=alpha,
            lambda_reg_gex=lambda_reg_gex,
            max_workers=max_workers
        )

        # Save results
        spot_names = self.gene_expression_adata.obs_names
        gene_names = self.gene_expression_adata.var_names
        cell_type_names = list(self.cell_profile_dict.keys())
        N = deconvolution_expression_data.shape[0]
        M = deconvolution_expression_data.shape[1]
        T = len(cell_type_names)

        # Imputation for NaN spots
        nan_spots = [i for i in range(N) if i not in spotwise_gene_expression_profiles]
        for nan_spot in nan_spots:
            # Get the indices of neighbors at radius 2
            neighbor_indices = get_neighbors_with_fixed_radius(nan_spot, filtered_adata, radius=2, include_center=False)

            # Collect non-NaN neighbors
            neighbor_profiles = [
                spotwise_gene_expression_profiles[i]
                for i in neighbor_indices
                if i in spotwise_gene_expression_profiles
            ]

            if neighbor_profiles:
                # Compute the average proportions from neighbors
                imputed_profile = np.nanmean(neighbor_profiles, axis=0)
                spotwise_gene_expression_profiles[nan_spot] = imputed_profile
                logging.info(f"Imputed spot {nan_spot} using neighbors at radius 2.")
            else:
                logging.warning(f"No valid neighbors found to impute spot {nan_spot}. Leaving as NaN.")
            
    
        spot_celltype_indices = [f"{spot_names[i]}_{cell_type_names[j]}" for i in range(N) for j in range(T)]
        nan_matrix = np.full((T, M), np.nan)
        gene_expression_data_combined = np.vstack(
            [spotwise_gene_expression_profiles.get(i, nan_matrix) for i in range(N)]
        )

        # Validate shape
        expected_shape = (N * T, M)
        if gene_expression_data_combined.shape != expected_shape:
            raise ValueError(
                f"Invalid shape in combined results. Expected {expected_shape}, got {gene_expression_data_combined.shape}."
            )

        gene_expression_df = pd.DataFrame(gene_expression_data_combined, index=spot_celltype_indices, columns=gene_names)

        # Save to Parquet
        parquet_path = os.path.join(self.output_folder, f"{self.sample_name}_gene_expression.parquet")
        gene_expression_df.to_parquet(parquet_path, compression="gzip")
        print(f"✅ Gene expression data saved to '{parquet_path}'.")

    def append_proportions_to_adata(self, proportions_path = None):
        """
        Append cell type proportion results back into the AnnData object.

        Args:
            proportions_path (str): Path to the CSV file containing cell type proportions.

        Raises:
            ValueError: If `adata` is not initialized or spot indices do not match.

        Returns:
            None: Updates `self.adata` in place.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression AnnData object is not initialized in the model.")
            
        if proportions_path is None:
            proportions_path = os.path.join(self.output_folder, f'{self.sample_name}_cell_prop_results.csv')

        # Load proportions CSV
        spot_by_celltype_df = pd.read_csv(proportions_path, index_col=0)

        # Check index alignment
        if not all(spot_by_celltype_df.index == self.gene_expression_adata.obs_names):
            print(spot_by_celltype_df.index)
            print(self.gene_expression_adata.obs_names)
            raise ValueError("Spot indices do not match between the CSV and AnnData object. Please verify your data.")
        else:
            print("✅ Spot indices match between CSV and adata.obs.")

        # Append cell type proportions to `adata.obs`
        for cell_type in spot_by_celltype_df.columns:
            self.gene_expression_adata.obs[cell_type] = spot_by_celltype_df[cell_type]

        print("✅ Cell type proportions have been appended to adata.obs.")
        
        
    def append_gex_to_adata(self, parquet_path = None):
        """
        Append gene expression layers from a Parquet file back into the gene_expression_adata object.

        Args:
            parquet_path (str): Path to the Parquet file containing spot-by-celltype gene expression data.

        Returns:
            None: Updates the `gene_expression_adata` object in place.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")
        if parquet_path is None:
            parquet_path = os.path.join(self.output_folder, f"{self.sample_name}_gene_expression.parquet")

        # Step 1: Read the Parquet file into a pandas DataFrame
        table = pq.read_table(parquet_path)
        df = table.to_pandas()
        print("Parquet file loaded successfully.")

        # Step 2: Reset the index to extract 'Spot' and 'CellType'
        df = df.reset_index()
        df[['Spot', 'CellType']] = df['index'].str.rsplit('_', n=1, expand=True)
        df = df.drop(columns=['index'])
        print("Spot and CellType successfully split.")

        # Step 3: Process each unique CellType
        for celltype in df['CellType'].unique():
            # Subset DataFrame to this cell type
            celltype_df = df[df['CellType'] == celltype].drop(columns=['Spot', 'CellType'])
            celltype_df.index = df[df['CellType'] == celltype]['Spot']

            # Convert the DataFrame to a numpy array
            celltype_matrix = celltype_df.values

            # Ensure gene_expression_adata.X is dense
            X_dense = self.gene_expression_adata.X.toarray() if hasattr(self.gene_expression_adata.X, 'toarray') else self.gene_expression_adata.X

            if X_dense.shape != celltype_matrix.shape:
                raise ValueError(
                    f"Shape mismatch: CellType matrix {celltype_matrix.shape} does not match adata.X {X_dense.shape}"
                )

            # Calculate gene contribution for this cell type
            celltype_gene_matrix = celltype_matrix * X_dense

            # Add matrices to layers
            layer_name = celltype.replace(' ', '_')  # Sanitize layer name
            self.gene_expression_adata.layers[layer_name + "_contribution"] = celltype_matrix
            self.gene_expression_adata.layers[layer_name + "_genes"] = celltype_gene_matrix

            print(f"Added layers: {layer_name}_contribution, {layer_name}_genes (Shape: {celltype_matrix.shape})")

        print("All layers added successfully.")
        print("Available layers:", self.gene_expression_adata.layers.keys())

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
        