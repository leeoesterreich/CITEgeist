"""
Main CitegeistModel class for spatial transcriptomics deconvolution.

This module implements the two-pass CITEgeist algorithm for deconvolving
spatial transcriptomics data using CITE-seq profiles.
"""

# Standard library imports
import logging
import os
from typing import Any, Dict, List, Optional, Tuple

# Third-party imports
import numpy as np
import pandas as pd
import scanpy as sc

try:
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover - exercised only in reduced environments
    pq = None

# Local imports
# Production-path imports (cuOPT backend)
from .deconvolution.qp_solver import (  # pylint: disable=unused-import
    compute_global_prior,
    compute_marker_exclusivity,
    map_antibodies_to_profiles_v2,
    normalize_counts,
    optimize_cell_proportions_per_marker,
)
from .marker_utils import strip_antibody_suffix
from .morphology.segmentation import (
    compute_spot_nuclei_counts_from_adata,
    normalize_nuclei_counts_for_prior,
    save_segmentation_artifacts,
)
from .utils import (
    assert_neighborhood_size,
    cleanup_memory,
    compute_optimal_radius,
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
        "lambda_spatial": 0.0,
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
        "lambda_spatial": 0.0,
        "lambda_sparse": 0.1,
        "pass2_library_slack": 1.5,
    },
}


def _require_pyarrow(feature_name: str) -> None:
    """Raise a clear error when parquet support is requested without pyarrow."""
    if pq is None:
        raise ImportError(
            f"{feature_name} requires optional parquet support. Install `pyarrow` "
            "or use the project's conda environment before using parquet I/O."
        )


class CitegeistModel:
    """Orchestrates the CITEgeist pipeline for one sample.

    Typical sequence: ``split_adata()`` -> ``run_cell_proportion_model()`` (M3 cuOPT QP)
    -> ``run_sace_allocation()`` (M3-gex SACE). See module docstrings for details.

    Attributes:
        results (dict): Pipeline outputs keyed by stage; ``results["cell_prop"]`` holds
            the finetuned per-spot proportion DataFrame consumed by SACE allocation.
        resolution (str): Active resolution mode ("spot" or "cell").
        resolution_params (dict): Parameter preset for the active resolution.
        cell_profile_dict (dict | None): Antibody -> cell-type profiles once loaded.
        output_folder (str): Directory for logs and saved outputs.
    """

    def __init__(
        self,
        sample_name,
        adata=None,
        output_folder=None,
        *,
        simulation=False,
        gene_expression_adata=None,
        antibody_capture_adata=None,
        resolution="spot",
        resolution_overrides=None,
    ):
        """
        Initialize the CitegeistModel with an AnnData object and output folder.

        Args:
            sample_name (str): Required sample identifier; used for logging setup and output
                file/log naming.
            adata (AnnData, optional): Spatial transcriptomics data object.
            output_folder (str): Path to save results and outputs.
            simulation (bool): Flag indicating if the data comes from a simulation framework.
            gene_expression_adata (AnnData, optional): Gene expression AnnData object (for simulations).
            antibody_capture_adata (AnnData, optional): Antibody capture AnnData object (for simulations).
            resolution (str): Resolution mode selecting parameter presets from
                ``RESOLUTION_DEFAULTS`` — "spot" (default) or "cell". Governs spatial-smoothing
                defaults (e.g. laplacian_k, lambda_spatial, lambda_sparse) and the classification
                dispatch path used by ``run_cell_proportion_model``.
            resolution_overrides (dict, optional): Per-key overrides applied on top of the
                selected resolution preset; each key must already exist in that preset, or a
                ValueError is raised.
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
            raise ValueError(f"resolution must be one of {list(RESOLUTION_DEFAULTS.keys())}, got '{resolution}'")
        self.resolution = resolution
        self.resolution_params = dict(RESOLUTION_DEFAULTS[resolution])
        if resolution_overrides:
            for key, val in resolution_overrides.items():
                if key not in self.resolution_params:
                    raise ValueError(f"Unknown resolution parameter: '{key}'")
                self.resolution_params[key] = val

        logging.info("CitegeistModel initialized successfully.")

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

        logging.info("AnnData has been successfully split into 'gene_expression_adata' and 'antibody_capture_adata'.")

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

    def filter_gex(self, nonzero_percentage=0.01, mean_expression_threshold=1.0, min_counts=10):
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

        logging.info(
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

        logging.info("Filtered antibody capture data to %s spots present in gene expression data.", len(filtered_spots))

    def preprocess_gex(self, target_sum=10000):
        """
        Preprocess gene expression data with count-preserving normalization.
        """
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        # Preserve raw counts for SACE (requires unnormalized data)
        if "raw_counts" not in self.gene_expression_adata.layers:
            self.gene_expression_adata.layers["raw_counts"] = self.gene_expression_adata.X.copy()

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
            "Gene expression data normalized to %s counts per spot and validated for discrete count analysis.",
            target_sum,
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

        # Step 2: Store raw counts before any transformation (for cellularity QP)
        self.antibody_capture_adata.layers["raw_counts"] = matrix.copy()

        # Step 3: Validate initial data (no NaNs or Infs at start)
        if np.isnan(matrix).any() or np.isinf(matrix).any():
            raise ValueError("Antibody capture matrix contains NaN or Inf values before preprocessing.")

        # Step 4: Winsorize to cap extreme values
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

        logging.info("Antibody capture data preprocessing completed: Winsorized, CLR applied, no NaNs detected.")

    def preprocess_antibody_discrete(
        self,
        winsorize_lower: int = 5,
        winsorize_upper: int = 95,
        scale_per_marker: bool = True,
        scale_mode: str = "per_marker",
    ) -> None:
        """
        Preprocess antibody data for discrete cell assignment.

        Unlike preprocess_antibody(), this method does NOT apply per-spot normalization
        (CLR or row normalization), preserving the cellularity signal where more cells
        in a spot produce higher total antibody signal.

        Args:
            winsorize_lower: Lower percentile for winsorization (default: 5)
            winsorize_upper: Upper percentile for winsorization (default: 95)
            scale_per_marker: DEPRECATED - use scale_mode instead. If True and
                scale_mode not specified, uses 'per_marker' scaling.
            scale_mode: Scaling mode after winsorization. Options:
                - 'per_marker': Scale each marker independently to [0, 1] (default)
                - 'global': Scale all markers together using global min/max.
                    Better for datasets with rare cell types, as it preserves
                    relative signal magnitude across markers.
                - 'none': No scaling, just winsorization.

        Note:
            This preprocessing is designed for discrete cell assignment where
            the total signal intensity correlates with cell count. The standard
            preprocess_antibody() normalizes away this information.

            For datasets with rare cell types, use scale_mode='global' to avoid
            inflating signal from low-abundance markers.
        """
        if self.antibody_capture_adata is None:
            raise ValueError("Antibody capture data has not been split. Run `split_adata` first.")

        # Step 1: Extract and ensure matrix is dense
        matrix = (
            self.antibody_capture_adata.X.toarray()
            if hasattr(self.antibody_capture_adata.X, "toarray")
            else self.antibody_capture_adata.X
        )
        matrix = np.asarray(matrix, dtype=np.float64)

        # Step 2: Validate initial data
        if np.isnan(matrix).any() or np.isinf(matrix).any():
            raise ValueError("Antibody capture matrix contains NaN or Inf values before preprocessing.")

        # Step 3: Winsorize per marker (column-wise) to cap extreme values
        for col_idx in range(matrix.shape[1]):
            col = matrix[:, col_idx]
            lower_bound = np.percentile(col, winsorize_lower)
            upper_bound = np.percentile(col, winsorize_upper)
            matrix[:, col_idx] = np.clip(col, lower_bound, upper_bound)

        # Step 4: Scaling based on scale_mode
        # Handle deprecated scale_per_marker parameter
        if not scale_per_marker and scale_mode == "per_marker":
            scale_mode = "none"

        if scale_mode == "per_marker":
            # Scale each marker independently to [0, 1]
            for col_idx in range(matrix.shape[1]):
                col = matrix[:, col_idx]
                col_min = col.min()
                col_max = col.max()
                if col_max > col_min:
                    matrix[:, col_idx] = (col - col_min) / (col_max - col_min)
                else:
                    # Constant column - set to 0
                    matrix[:, col_idx] = 0.0
        elif scale_mode == "global":
            # Scale all markers together using global min/max
            # This preserves relative signal magnitude across markers,
            # preventing rare cell type markers from being inflated
            global_min = matrix.min()
            global_max = matrix.max()
            if global_max > global_min:
                matrix = (matrix - global_min) / (global_max - global_min)
            else:
                matrix = np.zeros_like(matrix)
        elif scale_mode == "none":
            # No scaling, just use winsorized values
            pass
        else:
            raise ValueError(f"Invalid scale_mode: {scale_mode}. Use 'per_marker', 'global', or 'none'.")

        # Step 5: Final validation
        if np.isnan(matrix).any() or np.isinf(matrix).any():
            raise ValueError("NaN or Inf values detected in antibody capture matrix after preprocessing.")

        # Step 6: Reassign processed matrix
        self.antibody_capture_adata.X = matrix

        # Update status flag
        self.preprocessed_antibody = True

        # Log row sum statistics to confirm cellularity signal is preserved
        row_sums = matrix.sum(axis=1)
        logging.info(
            "Discrete antibody preprocessing complete: " "Row sums range [%s, %s], mean=%s, std=%s",
            row_sums.min(),
            row_sums.max(),
            row_sums.mean(),
            row_sums.std(),
        )
        logging.info(
            f"Antibody capture data preprocessing completed for discrete mode: "
            f"Winsorized [{winsorize_lower}%, {winsorize_upper}%], "
            f"scale_mode={scale_mode}, no per-spot normalization."
        )

    def compute_spot_nuclei_counts(
        self,
        resolution_mode: str = "hires",
        max_fullres_side: int = 9000,
        save_masks: bool = True,
        modality: str = "he",
        **stardist_kwargs,
    ) -> pd.Series:
        """
        Compute per-spot nuclei counts from Visium histology using StarDist.

        Writes the following columns to both gene/protein AnnData .obs (when present):
            - nuclei_count_raw
            - nuclei_count_target

        Args:
            resolution_mode: Image resolution to use ('lowres', 'hires', 'fullres').
            max_fullres_side: Max image dimension for fullres fallback.
            save_masks: Whether to save StarDist mask arrays to disk.
            modality: StarDist modality - ``"he"`` for H&E or ``"dapi"`` for fluorescence.
            **stardist_kwargs: Extra arguments forwarded to ``StarDistSegmenter.segment()``
                (e.g. ``prob_thresh``, ``nms_thresh``, ``scale``).

        Returns:
            pd.Series: Raw nuclei counts indexed by spot.
        """
        source_adata = self.antibody_capture_adata
        if source_adata is None:
            source_adata = self.gene_expression_adata
        if source_adata is None:
            source_adata = self.adata
        if source_adata is None:
            raise ValueError("No AnnData available to run nuclei segmentation.")

        seg_result = compute_spot_nuclei_counts_from_adata(
            adata=source_adata,
            resolution_mode=resolution_mode,
            max_fullres_side=max_fullres_side,
            modality=modality,
            **stardist_kwargs,
        )
        target = normalize_nuclei_counts_for_prior(seg_result.nuclei_count_raw)

        # Attach to available AnnData objects using aligned index assignment
        for ad in (self.gene_expression_adata, self.antibody_capture_adata, self.adata):
            if ad is None:
                continue
            common = ad.obs_names.intersection(seg_result.nuclei_count_raw.index)
            if len(common) == 0:
                continue
            ad.obs.loc[common, "nuclei_count_raw"] = seg_result.nuclei_count_raw.loc[common].astype(float).values
            ad.obs.loc[common, "nuclei_count_target"] = target.loc[common].astype(float).values

        output_paths = save_segmentation_artifacts(
            output_folder=self.output_folder,
            sample_name=self.sample_name,
            result=seg_result,
            save_masks=save_masks,
        )
        self.results["nuclei_counts"] = {
            "resolution_mode": resolution_mode,
            "n_nuclei": int(seg_result.nuclei_count_raw.sum()),
            "n_spots": int(seg_result.nuclei_count_raw.shape[0]),
            "image_shape": seg_result.image_shape,
            "outputs": output_paths,
        }
        logging.info(
            "Computed nuclei counts from StarDist (%s): nuclei=%d, spots=%d.",
            resolution_mode,
            int(seg_result.nuclei_count_raw.sum()),
            int(seg_result.nuclei_count_raw.shape[0]),
        )
        return seg_result.nuclei_count_raw

    def run_cell_proportion_model(
        self,
        *,
        radius=None,
        tolerance=1e-4,
        max_iterations=20,
        lambda_reg=1,
        alpha=0.5,
        _max_y_change=0.4,
        _max_workers=None,
        _checkpoint_interval=100,
        unknown_threshold=0.05,
        min_celltype_threshold=0.01,
        redundancy_threshold=0.1,
        validation_warn_only=False,
        # Laplacian smoothing: small weight reduces spatial error clustering (+0.002 r)
        lambda_laplacian=0.03,
        laplacian_k=8,
        # Per-marker beta parameters
        beta_min=0.1,
        beta_max=2.0,
        lambda_coverage=1.0,
        # Optional nuclei abundance prior (spot mode)
        use_nuclei_prior=False,
        nuclei_prior_lambda=1.0,
        nuclei_target_col="nuclei_count_target",
        nuclei_prior_bounds=None,  # Override row-sum bounds, e.g. (0.1, 5.0)
        # Cellularity-aware Laplacian (bilateral smoothing)
        use_cellularity_laplacian=False,
        cellularity_col=None,  # Column in antibody_capture_adata.obs with raw nuclei counts
        cellularity_sigma=0.5,  # Bandwidth for cellularity similarity kernel
        # IQP→QP cascade: sparsity mask from discrete assignment
        sparsity_mask=None,  # (N, T) array of per-spot-per-type upper bounds
        # GMM detection gating: run per-type GMM to identify absent types (r=0.766)
        use_detection_gating=True,
        detection_gate_ub=0.01,  # UB for gated-out types (lowered from 0.05: structural FP reduction)
        # GEX-informed detection (requires use_detection_gating=True)
        use_gex_detection=True,
        gex_detection_k=10,
        gex_detection_min_corr=0.15,
        fusion_mode="union",
        # Iterative sparsity refinement (opt-in)
        refine_sparsity=False,
        refine_suppress_threshold=0.02,
        refine_rescue_threshold=0.08,
        # Entropy-weighted markers (opt-in)
        use_entropy_weights=False,
        entropy_weight_alpha=1.0,
        # Direct per-marker weight override (dict: marker_name -> weight, default 1.0 for unlisted)
        marker_weight_dict=None,
        # Confusion penalty (opt-in): penalize confused type pairs
        lambda_confusion=0.0,
        confusion_pairs_manual=None,  # List of (type_name_a, type_name_b) tuples, or None for auto-detect
        confusion_corr_threshold=-0.25,
        # Cellularity-scaled QP (count-space solver)
        nuclei_counts=None,  # pd.Series of per-spot nuclei counts
        _cellularity_slack=0.3,  # δ fraction for hard constraints
        _lambda_cellularity=1.0,  # soft penalty weight for count-sum deviation
        # Cell classification parameters (cell resolution only)
        use_gating=None,
        priority_dict=None,
        threshold_method="auto",
        use_negative_gates=False,
        # Method selection
        method="qp",  # "qp" (cuOPT, default) or "per_type_beta"
        # Per-type beta params (ignored when method != "per_type_beta")
        sigma_scale=1.0,
        # NB-specific params (ignored when method="qp")
        nb_device="cpu",
        nb_n_iter=100,
        nb_use_detection=True,
        nb_detection_threshold=0.5,
        nb_gpu_adam_steps=200,
        # Morphology prior (refinement)
        morphology_device="cpu",
    ):
        """
        Orchestrates the cell proportion optimization workflow.

        Args:
            tolerance (float): Convergence tolerance for EM algorithm
            max_iterations (int): Maximum number of iterations
            lambda_reg (float): Regularization strength
            alpha (float): L1-L2 tradeoff factor (0 = L2, 1 = L1)
            unknown_threshold (float): Maximum allowed mean proportion for Unknown cell type (default: 0.05 = 5%)
            min_celltype_threshold (float): Minimum required mean proportion for defined cell types (default: 0.01 = 1%)
            redundancy_threshold (float): Maximum allowed fraction of redundant cell types (default: 0.1 = 10%)
            lambda_laplacian (float): Weight for Laplacian spatial smoothing (default: 0.0, disabled)
            laplacian_k (int): Number of neighbors for Laplacian graph (default: 8)
            beta_min (float): Minimum allowed beta value for per-marker optimization (default: 0.1)
            beta_max (float): Maximum allowed beta value for per-marker optimization (default: 2.0)
            lambda_coverage (float): Exponent for marker-count asymmetric loss scaling. Default: 1.0

        Args (production-relevant; see signature for the full set):
            method: proportion solver variant ("qp" — per-marker cuOPT QP — is production;
                also "per_type_beta").
            fusion_mode: detection-mask fusion mode ("union" default).
            use_gex_detection / gex_detection_k / gex_detection_min_corr: GEX-informed detection fusion.
            use_detection_gating / detection_gate_ub: GMM detection gating controls.
            radius: neighborhood radius for spatial terms.
            nuclei_counts: optional per-spot nuclei counts (pd.Series) enabling the
                cellularity-scaled count-space QP; None disables it.
            use_nuclei_prior / nuclei_prior_lambda / nuclei_target_col: optional
                nuclei-abundance row-sum prior (spot mode).
            use_cellularity_laplacian / cellularity_col / cellularity_sigma: bilateral
                spatial smoothing weighted by per-spot nuclei similarity.

        Returns:
            tuple[pandas.DataFrame, pandas.DataFrame]: (global, finetuned) per-spot cell-type
            proportions (rows = spot barcodes, columns = cell types); the finetuned frame is also
            stored on the model (``self.results["cell_prop"]``) for downstream SACE allocation.
            Other dispatch paths (cell-resolution gating) return their own
            result type instead.
        """

        # --- Detection gating guard ---
        # GMM detection gating is critical for accurate proportions (r=0.766 vs 0.36 without).
        # Disabling it is almost always a mistake in benchmarks and production.
        if not use_detection_gating:
            import warnings  # pylint: disable=import-outside-toplevel

            warnings.warn(
                "use_detection_gating=False: GMM detection gating is disabled. "
                "This typically degrades proportion accuracy (r drops ~0.40). "
                "Only disable for ablation studies or debugging.",
                UserWarning,
                stacklevel=2,
            )

        # --- Method dispatch ---
        if method not in ("qp", "per_type_beta"):
            raise ValueError(f"Unknown method '{method}'. Use 'qp' or 'per_type_beta'.")

        # --- QP path (cuOPT backend) ---

        # Use resolution preset for laplacian_k if caller used default
        if laplacian_k == 8 and self.resolution == "cell":
            laplacian_k = self.resolution_params["laplacian_k"]

        # Use resolution preset for lambda_laplacian if caller used old default
        if lambda_laplacian == 0.1 and self.resolution == "cell":
            lambda_laplacian = self.resolution_params["lambda_spatial"]  # now 0.0

        if radius is None:
            # Auto-detect optimal radius from spatial coordinates
            source_adata = self.antibody_capture_adata or self.gene_expression_adata
            if source_adata is None:
                raise ValueError("No AnnData available for radius auto-detection")
            radius = compute_optimal_radius(source_adata)
            logging.info("Auto-detected radius: %s (3 rings)", radius)

        if self.adata is None and (self.gene_expression_adata is None or self.antibody_capture_adata is None):
            raise ValueError("No valid data loaded. Ensure `adata` or split datasets are loaded properly.")

        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary has not been loaded. Run `load_cell_profile_dict` first.")

        # Extract spatial coordinates for Laplacian smoothing
        coords = None
        if lambda_laplacian > 0:
            coords = self.antibody_capture_adata.obsm.get("spatial", None)
            if coords is None and self.gene_expression_adata is not None:
                coords = self.gene_expression_adata.obsm.get("spatial", None)
            if coords is not None:
                logging.info("Using Laplacian smoothing with lambda=%s, k=%s", lambda_laplacian, laplacian_k)
            else:
                logging.warning("No spatial coordinates found for Laplacian smoothing - disabling")
                lambda_laplacian = 0

        # Dispatch to gating-based classification for cell resolution
        if use_gating is None:
            use_gating = self.resolution == "cell"
        if use_gating:
            logging.info("Cell resolution: dispatching to gating-based classification")
            return self._run_cell_classification(
                threshold_method=threshold_method,
                priority_dict=priority_dict,
                use_negative_gates=use_negative_gates,
                coords=coords,
            )

        # === Cellularity-informed spot weighting ===
        # Compute spot weights from nuclei counts: sqrt(N/median(N))
        # Dense spots have more reliable protein signal → weight reconstruction more
        cellularity_spot_weights = None
        cellularity_array = None
        if nuclei_counts is not None:
            from .assignment.cellularity_utils import prepare_cellularity  # pylint: disable=import-outside-toplevel

            spot_names_for_cell = self.antibody_capture_adata.obs_names
            cellularity_array = prepare_cellularity(nuclei_counts, spot_names_for_cell)
            median_N = np.median(cellularity_array)
            cellularity_spot_weights = np.sqrt(cellularity_array / max(median_N, 1.0))
            logging.info(
                "Cellularity spot weights: median N=%.1f, weight range=[%.2f, %.2f]",
                median_N,
                cellularity_spot_weights.min(),
                cellularity_spot_weights.max(),
            )

        spot_names = self.antibody_capture_adata.obs_names
        spot_abundance_target = None
        lambda_abundance_prior = 0.0
        row_sum_bounds = None  # default (0.9, 1.2) inside solvers
        if use_nuclei_prior:
            import warnings  # pylint: disable=import-outside-toplevel

            warnings.warn(
                "use_nuclei_prior is deprecated and will be removed in a future version. "
                "Use nuclei_counts= parameter instead for cellularity-scaled QP deconvolution.",
                DeprecationWarning,
                stacklevel=2,
            )
            if self.antibody_capture_adata is None:
                raise ValueError("Antibody capture data must be loaded to use nuclei prior.")
            if nuclei_target_col not in self.antibody_capture_adata.obs.columns:
                raise ValueError(
                    f"Nuclei prior requested but column '{nuclei_target_col}' not found in antibody_capture_adata.obs. "
                    "Run compute_spot_nuclei_counts() first."
                )
            target_series = self.antibody_capture_adata.obs.loc[spot_names, nuclei_target_col]
            spot_abundance_target = target_series.to_numpy(dtype=float)
            if np.isnan(spot_abundance_target).any() or np.isinf(spot_abundance_target).any():
                raise ValueError(f"Invalid values in nuclei target column '{nuclei_target_col}'.")
            lambda_abundance_prior = float(nuclei_prior_lambda)
            row_sum_bounds = nuclei_prior_bounds if nuclei_prior_bounds is not None else (0.25, 2.5)
            logging.info(
                "Using nuclei abundance prior: lambda=%.3f, target_col='%s', median=%.3f, "
                "row_sum_bounds=(%.2f, %.2f)",
                lambda_abundance_prior,
                nuclei_target_col,
                float(np.median(spot_abundance_target)),
                row_sum_bounds[0],
                row_sum_bounds[1],
            )

        # Extract cellularity for Laplacian edge weighting (independent of abundance prior)
        cellularity_laplacian_array = None
        if use_cellularity_laplacian:
            col = cellularity_col or nuclei_target_col
            if col in self.antibody_capture_adata.obs.columns:
                cellularity_laplacian_array = self.antibody_capture_adata.obs.loc[spot_names, col].to_numpy(dtype=float)
                logging.info(
                    "Cellularity-aware Laplacian: col='%s', sigma=%.2f, " "median=%.1f, range=[%.0f, %.0f]",
                    col,
                    cellularity_sigma,
                    float(np.median(cellularity_laplacian_array)),
                    float(np.min(cellularity_laplacian_array)),
                    float(np.max(cellularity_laplacian_array)),
                )
            else:
                logging.warning(
                    "use_cellularity_laplacian=True but column '%s' not found. " "Falling back to standard Laplacian.",
                    col,
                )

        recon_error = None  # Will be set by per-marker optimization path

        # Per-marker beta optimization (preserves marker-level signal variation)
        logging.info("Using per-marker beta optimization (preserves marker-level signal variation)")

        # Get marker-level data instead of averaged profiles
        marker_level_data, marker_names, assignment_matrix, cell_type_names = map_antibodies_to_profiles_v2(
            self.antibody_capture_adata, self.cell_profile_dict
        )

        # For per_type_beta: expand to all 17 markers from MARKER_TYPE_TABLE.
        # Keep the 7-Major assignment_matrix for detection gating (unchanged).
        if method == "per_type_beta":
            from .deconvolution.emission_init import (  # pylint: disable=import-outside-toplevel
                _strip_suffix,
                build_marker_config,
            )

            available = [_strip_suffix(v) for v in self.antibody_capture_adata.var_names]
            ptb_markers, _, ptb_type_names = build_marker_config(available)

            if len(ptb_markers) > len(marker_names):
                # Extract full 17-marker data from antibody adata
                ab_data = self.antibody_capture_adata.X
                ab_data = ab_data.toarray() if hasattr(ab_data, "toarray") else np.asarray(ab_data)
                var_canonical = [_strip_suffix(v) for v in self.antibody_capture_adata.var_names]
                canon_to_idx = {}
                for idx, cn in enumerate(var_canonical):
                    if cn not in canon_to_idx:
                        canon_to_idx[cn] = idx

                ptb_col_idx = [canon_to_idx[m] for m in ptb_markers if m in canon_to_idx]
                ptb_markers_found = [m for m in ptb_markers if m in canon_to_idx]

                ptb_data = ab_data[:, ptb_col_idx].astype(np.float64)
                # Column-max normalize (same as map_antibodies_to_profiles_v2)
                col_max = np.max(ptb_data, axis=0)
                col_max[col_max == 0] = 1e-6
                ptb_data = ptb_data / col_max

                logging.info(
                    "Per-type beta: expanded markers %d -> %d (%s)",
                    len(marker_names),
                    len(ptb_markers_found),
                    ", ".join(set(ptb_markers_found) - set(marker_names)),
                )
                # Replace marker data for the optimizer (detection still uses original assignment_matrix)
                marker_level_data = ptb_data
                marker_names = ptb_markers_found
                cell_type_names = ptb_type_names

        # Cellularity correction: divide marker data by log(N+1)/log(median(N)+1)
        # so dense spots need proportionally more signal. Proven +0.020 r on Xenium.
        if cellularity_array is not None:
            median_N = float(np.median(cellularity_array))
            log_scale = np.log(cellularity_array + 1.0) / np.log(median_N + 1.0)
            log_scale = np.maximum(log_scale, 0.3)  # floor to avoid division by tiny N
            marker_level_data = marker_level_data / log_scale[:, np.newaxis]
            logging.info(
                "Cellularity correction (S/logN): median_N=%.1f, " "scale range=[%.3f, %.3f]",
                median_N,
                log_scale.min(),
                log_scale.max(),
            )

        # GMM detection gating: identify absent types per spot
        if use_detection_gating and sparsity_mask is None:
            from .deconvolution.detection import (  # pylint: disable=import-outside-toplevel
                build_detection_marker_groups,
                detect_cell_types,
            )

            active_mask = assignment_matrix.T > 0  # (T, M)
            marker_groups = build_detection_marker_groups(active_mask, cell_type_names)

            # Map marker indices to raw antibody columns
            raw_ab = self.antibody_capture_adata.layers.get("raw_counts", self.antibody_capture_adata.X)
            raw_ab = raw_ab.toarray() if hasattr(raw_ab, "toarray") else np.asarray(raw_ab)

            def _strip_suffix(name: str) -> str:
                """Strip -1 suffix from antibody names (e.g., CD68-1 -> CD68)."""
                return strip_antibody_suffix(name)

            var_names = list(self.antibody_capture_adata.var_names)
            m2raw = {}
            for idx, vn in enumerate(var_names):
                c = _strip_suffix(vn)
                if c not in m2raw:
                    m2raw[c] = idx

            mapped_groups = {}
            for ct, m_indices in marker_groups.items():
                raw_idx = [
                    m2raw[_strip_suffix(marker_names[mi])]
                    for mi in m_indices
                    if _strip_suffix(marker_names[mi]) in m2raw
                ]
                if raw_idx:
                    mapped_groups[ct] = raw_idx

            detected = detect_cell_types(
                raw_ab,
                mapped_groups,
                threshold=0.5,
                log_transform=True,
                adaptive_threshold=True,
            )

            # Persist raw boolean detection mask for downstream Bayesian assignment
            detection_mask_bool = np.ones((marker_level_data.shape[0], len(cell_type_names)), dtype=bool)
            detected_types_list = list(mapped_groups.keys())
            for k, ct in enumerate(detected_types_list):
                if ct in cell_type_names:
                    t_idx = cell_type_names.index(ct)
                    detection_mask_bool[:, t_idx] = detected[:, k]
            # Ensure at least 2 types detected per spot (same rescue as sparsity_mask)
            for i in range(marker_level_data.shape[0]):
                if detection_mask_bool[i].sum() < 2:
                    detection_mask_bool[i] = True
            self.results["detection_mask_bool"] = detection_mask_bool
            logging.info(
                "Persisted detection_mask_bool: %d spots, %d types, %.1f%% gated",
                marker_level_data.shape[0],
                len(cell_type_names),
                100.0 * (~detection_mask_bool).sum() / (marker_level_data.shape[0] * len(cell_type_names)),
            )

            N_spots = marker_level_data.shape[0]
            T = len(cell_type_names)
            sparsity_mask = np.ones((N_spots, T), dtype=np.float64)
            detected_types = list(mapped_groups.keys())
            for k, ct in enumerate(detected_types):
                if ct in cell_type_names:
                    t_idx = cell_type_names.index(ct)
                    if cellularity_array is not None:
                        # Cellularity-informed gate: can't have less than 1 cell
                        # gate_ub_i = min(detection_gate_ub, 1/N_i)
                        per_spot_ub = np.minimum(detection_gate_ub, 1.0 / np.maximum(cellularity_array, 1.0))
                        sparsity_mask[~detected[:, k], t_idx] = per_spot_ub[~detected[:, k]]
                    else:
                        sparsity_mask[~detected[:, k], t_idx] = detection_gate_ub

            # Ensure at least 2 types open per spot
            for i in range(N_spots):
                if (sparsity_mask[i] > detection_gate_ub).sum() < 2:
                    sparsity_mask[i] = 1.0

            n_gated = (sparsity_mask < 1.0).sum()
            logging.info(
                "Detection gating: %d/%d (%.1f%%) spot-type pairs gated (ub=%.2f)",
                n_gated,
                N_spots * T,
                100.0 * n_gated / (N_spots * T),
                detection_gate_ub,
            )

        # GEX-informed detection fusion
        if use_gex_detection and use_detection_gating and sparsity_mask is not None:
            from .deconvolution.detection_refinement import (  # pylint: disable=import-outside-toplevel
                compute_gene_type_correlations,
                detect_cell_types_gex,
                fuse_detection_masks,
            )

            gex_X = self.gene_expression_adata.X
            gex_X = gex_X.toarray() if hasattr(gex_X, "toarray") else np.asarray(gex_X)
            gene_names_list = list(self.gene_expression_adata.var_names)

            ab_names_stripped = [_strip_suffix(n) for n in var_names]
            H = compute_gene_type_correlations(
                gex_X,
                raw_ab,
                ab_names_stripped,
                self.cell_profile_dict,
                cell_type_names,
            )
            self.results["gex_detection_corr"] = H

            gex_detected = detect_cell_types_gex(
                gex_X,
                H,
                gene_names_list,
                cell_type_names,
                k=gex_detection_k,
                min_corr=gex_detection_min_corr,
            )
            self.results["detection_protein"] = detection_mask_bool.copy()
            self.results["detection_gex"] = gex_detected

            # Fuse protein + GEX detection
            fused_detected = fuse_detection_masks(
                detection_mask_bool,
                gex_detected,
                assignment_matrix,
                mode=fusion_mode,
            )
            self.results["detection_fused"] = fused_detected

            # Rebuild sparsity mask from fused detection
            sparsity_mask = np.ones((N_spots, T), dtype=np.float64)
            detected_types = list(mapped_groups.keys())
            for _, ct in enumerate(detected_types):
                if ct in cell_type_names:
                    t_idx = cell_type_names.index(ct)
                    if cellularity_array is not None:
                        per_spot_ub = np.minimum(detection_gate_ub, 1.0 / np.maximum(cellularity_array, 1.0))
                        sparsity_mask[~fused_detected[:, t_idx], t_idx] = per_spot_ub[~fused_detected[:, t_idx]]
                    else:
                        sparsity_mask[~fused_detected[:, t_idx], t_idx] = detection_gate_ub

            for i in range(N_spots):
                if (sparsity_mask[i] > detection_gate_ub).sum() < 2:
                    sparsity_mask[i] = 1.0

            # Update detection_mask_bool to reflect fusion
            detection_mask_bool = fused_detected
            self.results["detection_mask_bool"] = detection_mask_bool

            n_gated_fused = (sparsity_mask < 1.0).sum()
            logging.info(
                "GEX-fused detection: %d/%d (%.1f%%) spot-type pairs gated",
                n_gated_fused,
                N_spots * T,
                100.0 * n_gated_fused / (N_spots * T),
            )

        # Marker weighting: explicit dict > entropy weights > uniform (None)
        marker_weight = None
        if marker_weight_dict is not None:
            marker_weight = np.array([marker_weight_dict.get(m, 1.0) for m in marker_names], dtype=np.float64)
            self.results["marker_weight_source"] = "explicit"
        elif use_entropy_weights:
            from .deconvolution.detection_refinement import (  # pylint: disable=import-outside-toplevel
                compute_marker_entropy_weights,
            )

            marker_weight = compute_marker_entropy_weights(
                marker_level_data,
                marker_names,
                alpha=entropy_weight_alpha,
            )
            self.results["marker_weight_source"] = "entropy"
        if marker_weight is not None:
            self.results["marker_weights"] = dict(zip(marker_names, marker_weight.tolist()))

        try:
            logging.info(
                "Running Stage 1 cell proportion optimization with validation thresholds: "
                "Unknown<%s%%, CellTypes>%s%%, Redundancy<%s%%",
                round(unknown_threshold * 100, 1),
                round(min_celltype_threshold * 100, 1),
                round(redundancy_threshold * 100),
            )

            if method == "per_type_beta":
                # --- Per-type beta EM path ---
                from .deconvolution.emission_init import (  # pylint: disable=import-outside-toplevel
                    build_beta_prior_sigma,
                    initialize_beta_matrix,
                )
                from .deconvolution.qp_solver import (  # pylint: disable=import-outside-toplevel
                    optimize_cell_proportions_per_type_beta,
                )

                # Build beta init from marker-level data
                raw_for_init = marker_level_data  # already cellularity-corrected if applicable
                median_N = float(np.median(cellularity_array)) if cellularity_array is not None else 1.0
                beta_init = initialize_beta_matrix(
                    raw_counts=raw_for_init,
                    markers=marker_names,
                    type_names=cell_type_names,
                    median_N=median_N,
                )
                prior_sigma = build_beta_prior_sigma(
                    markers=marker_names,
                    type_names=cell_type_names,
                    sigma_scale=sigma_scale,
                )
                logging.info(
                    "Per-type beta: beta_init shape=%s, prior_sigma shape=%s, sigma_scale=%.2f",
                    beta_init.shape,
                    prior_sigma.shape,
                    sigma_scale,
                )

                (
                    Y_values,
                    _,
                    marker_beta_matrix_dict,
                    alpha_values,
                    recon_error,
                    objective_history,
                ) = optimize_cell_proportions_per_type_beta(
                    marker_level_data=marker_level_data,
                    marker_names=marker_names,
                    cell_type_names=cell_type_names,
                    beta_init=beta_init,
                    prior_sigma=prior_sigma,
                    tolerance=tolerance,
                    max_iterations=max_iterations,
                    lambda_reg=lambda_reg,
                    alpha_elastic=alpha,
                    beta_max=beta_max,
                    lambda_laplacian=lambda_laplacian,
                    coords=coords,
                    laplacian_k=laplacian_k,
                    row_sum_bounds=row_sum_bounds,
                    sparsity_mask=sparsity_mask,
                    spot_weights=cellularity_spot_weights,
                    marker_weight=marker_weight,
                )

                # Store per_type_beta-specific results
                self.results["marker_beta_matrix"] = marker_beta_matrix_dict
                self.results["objective_history"] = objective_history

                # Create scalar marker_beta_dict for backward compat (max across types)
                marker_beta_dict = {m: max(marker_beta_matrix_dict[m].values()) for m in marker_beta_matrix_dict}
                # beta_values as 1-D array for downstream compat (kept for logging symmetry)
                _beta_values = np.array([marker_beta_dict[m] for m in marker_names], dtype=np.float64)  # noqa: F841

                logging.info(
                    "Per-type beta EM converged in %d iterations, final obj=%.4f",
                    len(objective_history),
                    objective_history[-1] if objective_history else float("nan"),
                )

            else:
                Y_values, _, marker_beta_dict, alpha_values, recon_error = optimize_cell_proportions_per_marker(
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
                    lambda_coverage=lambda_coverage,
                    spot_abundance_target=spot_abundance_target,
                    lambda_abundance_prior=lambda_abundance_prior,
                    row_sum_bounds=row_sum_bounds,
                    cellularity=cellularity_laplacian_array,
                    cellularity_sigma=cellularity_sigma,
                    sparsity_mask=sparsity_mask,
                    spot_weights=cellularity_spot_weights,
                    marker_weight=marker_weight,
                )

            # Post-QP subcell cleanup: zero out GMM-gated types below 1/(2*N_i)
            # Only applies to types the GMM flagged as absent — detected types are protected
            if cellularity_array is not None and sparsity_mask is not None:
                min_prop = 0.5 / np.maximum(cellularity_array, 1.0)  # (N_spots,)
                gated = sparsity_mask < 1.0  # (N_spots, T) — True where GMM said absent
                subcell_mask = gated & (Y_values < min_prop[:, np.newaxis])
                n_zeroed = subcell_mask.sum()
                if n_zeroed > 0:
                    Y_values[subcell_mask] = 0.0
                    row_sums = Y_values.sum(axis=1, keepdims=True)
                    row_sums = np.maximum(row_sums, 1e-10)
                    Y_values = Y_values / row_sums
                    logging.info(
                        "Subcell cleanup: zeroed %d/%d (%.1f%%) entries below 1/(2N)",
                        n_zeroed,
                        Y_values.size,
                        100.0 * n_zeroed / Y_values.size,
                    )

            # Iterative sparsity refinement: re-solve with QP-informed mask
            if refine_sparsity and sparsity_mask is not None:
                from .deconvolution.detection_refinement import (  # pylint: disable=import-outside-toplevel
                    refine_sparsity_from_proportions,
                )

                sparsity_mask_pass1 = sparsity_mask.copy()
                refined_mask = refine_sparsity_from_proportions(
                    Y_values,
                    sparsity_mask,
                    cellularity_array,
                    suppress_threshold=refine_suppress_threshold,
                    rescue_threshold=refine_rescue_threshold,
                    detection_gate_ub=detection_gate_ub,
                )

                n_changes = int((refined_mask != sparsity_mask).sum())
                logging.info(
                    "Sparsity refinement: %d/%d mask entries changed (%.1f%%)",
                    n_changes,
                    N_spots * T,
                    100.0 * n_changes / (N_spots * T),
                )

                if n_changes > 0:
                    logging.info("Re-solving QP with refined sparsity mask (Pass 2)...")
                    Y_values, _, marker_beta_dict, alpha_values, recon_error = optimize_cell_proportions_per_marker(
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
                        lambda_coverage=lambda_coverage,
                        spot_abundance_target=spot_abundance_target,
                        lambda_abundance_prior=lambda_abundance_prior,
                        row_sum_bounds=row_sum_bounds,
                        cellularity=cellularity_laplacian_array,
                        cellularity_sigma=cellularity_sigma,
                        sparsity_mask=refined_mask,
                        spot_weights=cellularity_spot_weights,
                        marker_weight=marker_weight,
                    )

                    # Re-run subcell cleanup with refined mask
                    if cellularity_array is not None:
                        min_prop = 0.5 / np.maximum(cellularity_array, 1.0)
                        gated = refined_mask < 1.0
                        subcell_mask = gated & (Y_values < min_prop[:, np.newaxis])
                        n_zeroed = subcell_mask.sum()
                        if n_zeroed > 0:
                            Y_values[subcell_mask] = 0.0
                            row_sums = Y_values.sum(axis=1, keepdims=True)
                            row_sums = np.maximum(row_sums, 1e-10)
                            Y_values = Y_values / row_sums

                    sparsity_mask = refined_mask

                self.results["sparsity_mask_pass1"] = sparsity_mask_pass1
                self.results["sparsity_mask_refined"] = refined_mask

            # Confusion penalty: use manual pairs or auto-detect, then re-solve
            if lambda_confusion > 0:
                confusion_pairs = []

                if confusion_pairs_manual:
                    # Manual pairs: convert type names to indices
                    for name_a, name_b in confusion_pairs_manual:
                        if name_a in cell_type_names and name_b in cell_type_names:
                            idx_a = cell_type_names.index(name_a)
                            idx_b = cell_type_names.index(name_b)
                            confusion_pairs.append((idx_a, idx_b))
                            logging.info("Manual confusion pair: %s ↔ %s", name_a, name_b)
                else:
                    # Auto-detect from proportion anti-correlations
                    # Use proportion values directly — types that compete for the same signal
                    corr_matrix = np.corrcoef(Y_values.T)  # (T, T)
                    for i in range(T):
                        for j in range(i + 1, T):
                            if corr_matrix[i, j] < confusion_corr_threshold:
                                confusion_pairs.append((i, j))
                                logging.info(
                                    "Confusion pair detected: %s ↔ %s (corr=%.3f)",
                                    cell_type_names[i],
                                    cell_type_names[j],
                                    corr_matrix[i, j],
                                )

                if confusion_pairs:
                    logging.info(
                        "Re-solving QP with confusion penalty (λ=%.3f, %d pairs)...",
                        lambda_confusion,
                        len(confusion_pairs),
                    )
                    Y_values, _, marker_beta_dict, alpha_values, recon_error = optimize_cell_proportions_per_marker(
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
                        lambda_coverage=lambda_coverage,
                        spot_abundance_target=spot_abundance_target,
                        lambda_abundance_prior=lambda_abundance_prior,
                        row_sum_bounds=row_sum_bounds,
                        cellularity=cellularity_laplacian_array,
                        cellularity_sigma=cellularity_sigma,
                        sparsity_mask=sparsity_mask,
                        spot_weights=cellularity_spot_weights,
                        marker_weight=marker_weight,
                        confusion_pairs=confusion_pairs,
                        lambda_confusion=lambda_confusion,
                    )

                    # Re-run subcell cleanup
                    if cellularity_array is not None and sparsity_mask is not None:
                        min_prop = 0.5 / np.maximum(cellularity_array, 1.0)
                        gated = sparsity_mask < 1.0
                        subcell_mask = gated & (Y_values < min_prop[:, np.newaxis])
                        n_zeroed = subcell_mask.sum()
                        if n_zeroed > 0:
                            Y_values[subcell_mask] = 0.0
                            row_sums = Y_values.sum(axis=1, keepdims=True)
                            row_sums = np.maximum(row_sums, 1e-10)
                            Y_values = Y_values / row_sums

                    self.results["confusion_pairs"] = [
                        (cell_type_names[a], cell_type_names[b]) for a, b in confusion_pairs
                    ]
                    if not confusion_pairs_manual:
                        self.results["confusion_corr_matrix"] = corr_matrix
                else:
                    logging.info("No confusion pairs found (threshold=%.2f)", confusion_corr_threshold)

            # Store cellularity results (post-hoc conversion)
            if cellularity_array is not None:
                from .assignment.cellularity_utils import (  # pylint: disable=import-outside-toplevel
                    round_counts_largest_remainder,
                )

                spot_names = self.antibody_capture_adata.obs_names
                c_continuous = Y_values * cellularity_array[:, np.newaxis]
                c_df = pd.DataFrame(c_continuous, index=spot_names, columns=cell_type_names)
                self.results["cell_counts_continuous"] = c_df

                int_counts = np.zeros_like(c_continuous, dtype=np.int64)
                for i in range(len(spot_names)):
                    int_counts[i] = round_counts_largest_remainder(c_continuous[i], round(cellularity_array[i]))
                int_df = pd.DataFrame(int_counts, index=spot_names, columns=cell_type_names)
                self.results["cell_counts_integer"] = int_df

            # Store marker betas and baselines for downstream analysis
            self.results["marker_beta"] = marker_beta_dict
            marker_alpha_dict = {marker_names[i]: alpha_values[i] for i in range(len(marker_names))}
            self.results["marker_alpha"] = marker_alpha_dict
            for m_idx, m_name in enumerate(marker_names):
                if alpha_values[m_idx] > 0.05:
                    logging.info("  Marker baseline: %s = %s", m_name, alpha_values[m_idx])

            # Compute marker exclusivity scores for finetuning
            # (skip for per_type_beta — exclusivity is encoded in the beta matrix itself)
            if method != "per_type_beta":
                marker_owners = []
                for m_idx in range(assignment_matrix.shape[0]):
                    owners = [j for j in range(assignment_matrix.shape[1]) if assignment_matrix[m_idx, j] > 0]
                    marker_owners.append(owners)

                marker_exclusivity = compute_marker_exclusivity(
                    marker_level_data=marker_level_data,
                    Y_values=Y_values,
                    marker_owners=marker_owners,
                    _assignment_matrix=assignment_matrix,
                )

                # Log exclusivity scores
                for m_idx, m_name in enumerate(marker_names):
                    if marker_owners[m_idx]:
                        logging.info("  Marker exclusivity: %s = %s", m_name, marker_exclusivity[m_idx])
                self.results["marker_exclusivity"] = {
                    marker_names[i]: marker_exclusivity[i] for i in range(len(marker_names))
                }

        except ValueError as e:
            error_msg = f"Cell proportion validation failed for sample '{self.sample_name}': {str(e)}"
            logging.error(error_msg)
            raise ValueError(error_msg) from e

        global_cell_type_proportions_df = pd.DataFrame(Y_values, index=spot_names, columns=cell_type_names)
        finetuned_cell_type_proportions_df = global_cell_type_proportions_df.copy()

        global_cell_type_proportions_df = global_cell_type_proportions_df.sort_index()
        finetuned_cell_type_proportions_df = finetuned_cell_type_proportions_df.sort_index()

        # Add per-spot reconstruction error if available (from per-marker beta optimization)
        if recon_error is not None:
            recon_series = pd.Series(recon_error, index=spot_names, name="recon_error")
            global_cell_type_proportions_df["recon_error"] = recon_series

        global_cell_type_proportions_df.to_csv(
            os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_global_results.csv")
        )
        finetuned_cell_type_proportions_df.to_csv(
            os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_finetuned_results.csv")
        )

        # Store to self.results for use by run_cell_expression_pass1
        self.results["cell_prop"] = finetuned_cell_type_proportions_df

        # Store base proportions for Bayesian assignment prior (pre-morphology)
        if "cell_prop_base" not in self.results:
            self.results["cell_prop_base"] = self.results["cell_prop"].copy()

        return global_cell_type_proportions_df, finetuned_cell_type_proportions_df

    def discretize_proportions(
        self,
        proportions_df: pd.DataFrame,
        nuclei_counts: pd.Series,
    ) -> pd.DataFrame:
        """
        Convert continuous cell type proportions to integer cell counts.

        Uses the largest remainder method (Hamilton/Vinton method) to allocate
        integer cell counts while ensuring sum(counts) = nuclei_count for each spot.
        This is the recommended approach when discrete cell counts are needed,
        as it achieves ~94% of continuous model performance vs ~69% for direct IQP.

        Args:
            proportions_df: DataFrame with cell type columns and proportion values (0-1).
                Output from run_cell_proportion_model().
            nuclei_counts: Series with nuclei count per spot (from StarDist or similar).

        Returns:
            DataFrame with cell type columns and integer count values per spot.
            Guaranteed: sum of counts per spot == nuclei_count for that spot.

        Example:
            >>> # Run continuous model first
            >>> global_props, finetuned_props = model.run_cell_proportion_model()
            >>> # Get nuclei counts from StarDist segmentation
            >>> nuclei_counts = pd.Series({'spot_1': 10, 'spot_2': 8, ...})
            >>> # Discretize
            >>> cell_counts = model.discretize_proportions(finetuned_props, nuclei_counts)

        Note:
            The largest remainder method:
            1. Computes raw_counts = proportions * nuclei_count
            2. Takes floor of each value
            3. Allocates remaining cells to cell types with largest fractional parts
            This guarantees exact sum while minimizing total rounding error.
        """
        from .assignment.cellularity_utils import (  # pylint: disable=import-outside-toplevel
            round_counts_largest_remainder,
        )

        # Align indices
        common_spots = proportions_df.index.intersection(nuclei_counts.index)
        if len(common_spots) == 0:
            raise ValueError("No common spots between proportions_df and nuclei_counts")

        if len(common_spots) < len(proportions_df):
            logging.warning(
                "Only %s/%s spots have nuclei counts. Missing spots will be excluded.",
                len(common_spots),
                len(proportions_df),
            )

        props = proportions_df.loc[common_spots]
        nuclei = nuclei_counts.loc[common_spots]

        cell_counts = pd.DataFrame(index=common_spots, columns=props.columns, dtype=int)

        for spot in common_spots:
            N = int(nuclei[spot])
            p = props.loc[spot].values.astype(float)

            if N == 0:
                cell_counts.loc[spot] = 0
                continue

            floor_counts = round_counts_largest_remainder(p, N)

            cell_counts.loc[spot] = floor_counts

        logging.info(
            "Discretized %s spots: total cells=%s, mean per spot=%s",
            len(common_spots),
            cell_counts.values.sum(),
            round(cell_counts.sum(axis=1).mean(), 1),
        )

        return cell_counts

    def run_cell_expression_pass1(self):
        """Run gene expression deconvolution using SACE (single-pass allocation).

        Convenience wrapper that calls run_sace_allocation(output_mode="layers")
        and stores results under the canonical 'gene_expression_pass1' key.

        Returns:
            Dict[str, Any]: {
                'spotwise_profiles': Dict[int, np.ndarray],
                'dimensions': None
            }
        """
        if not self.preprocessed_gex:
            raise ValueError("Gene expression data not preprocessed. Run preprocess_gex() first.")

        logging.info("Using SACE GEX deconvolution (single-pass, no solver required)")

        spot_type_gex, _ = self.run_sace_allocation(
            output_mode="layers",
        )

        self.results["gene_expression_pass1"] = spot_type_gex

        return {"spotwise_profiles": spot_type_gex, "dimensions": None}

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
        _gene_names = self.gene_expression_adata.var_names  # noqa: F841

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
        logging.info("Number of cell types: %s", T)
        logging.info("Number of genes: %s", M)

        # Per cell-type statistics
        for t, cell_type in enumerate(cell_type_names):
            mean_conf = np.mean(prior_info["confidence_scores"][t])
            strong_signals = np.mean(prior_info["global_prior"][t] > 0.5)
            logging.info("\n%s:", cell_type)
            logging.info(" - Mean confidence score: %s", mean_conf)
            logging.info(" - % Strong signals: %s%", 100 * strong_signals)

            # Expression pattern summary
            mean_exp = np.mean(prior_info["expression_patterns"]["mean_expression"][t])
            freq = np.mean(prior_info["expression_patterns"]["expression_frequency"][t])
            cons = np.mean(prior_info["expression_patterns"]["expression_consistency"][t])
            logging.info(" - Mean expression: %s", mean_exp)
            logging.info(" - Mean expression frequency: %s", freq)
            logging.info(" - Mean expression consistency: %s", cons)

        return prior_info

    def _save_profiles_to_parquet(self, profiles, path):
        """Helper method to save profiles to parquet format with consistent naming."""
        _require_pyarrow("Saving profiles to parquet")
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
        logging.info("Saved profiles to %s with cell types: %s", path, cell_type_names)

    def append_proportions_to_adata(self, proportions_path=None, key="finetuned"):
        """Append cell type proportions to AnnData object."""
        if proportions_path is None:
            proportions_path = os.path.join(self.output_folder, f"{self.sample_name}_cell_prop_{key}_results.csv")

        # Load proportions CSV
        spot_by_celltype_df = pd.read_csv(proportions_path, index_col=0)
        spot_by_celltype_df.columns = spot_by_celltype_df.columns.str.replace(" ", "_")

        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data not available")

        if "spot_" in str(spot_by_celltype_df.index[0]):
            # Sort both numerically by the spot number
            def get_spot_number(x):
                return int(x.split("spot_")[1])

            # Sort using reindex instead of sort_index
            sorted_csv_idx = sorted(spot_by_celltype_df.index, key=get_spot_number)
            sorted_adata_idx = sorted(self.gene_expression_adata.obs_names, key=get_spot_number)

            spot_by_celltype_df = spot_by_celltype_df.reindex(sorted_csv_idx)
            self.gene_expression_adata = self.gene_expression_adata[sorted_adata_idx].copy()

        # Check if indices match after sorting
        if not all(spot_by_celltype_df.index == self.gene_expression_adata.obs_names):
            raise ValueError("Spot indices still don't match after sorting. Please verify your data.")

        # Add cell type proportions to adata.obs
        for cell_type in spot_by_celltype_df.columns:
            self.gene_expression_adata.obs[cell_type] = spot_by_celltype_df[cell_type]

        self.results["cell_prop"] = spot_by_celltype_df

        logging.info("✅ Cell type proportions have been appended to adata.obs and results['cell_prop']")

    def append_gex_to_adata(self, parquet_path=None, pass_number=1):
        """
        Append gene expression layers from a Parquet file back into the gene_expression_adata object.
        """
        _require_pyarrow("Appending parquet-backed gene expression")
        if self.gene_expression_adata is None:
            raise ValueError("Gene expression data has not been split. Run `split_adata` first.")

        if parquet_path is None:
            parquet_path = os.path.join(
                self.output_folder, f"{self.sample_name}_gene_expression_pass{pass_number}.parquet"
            )

        # Step 1: Read the Parquet file into a pandas DataFrame
        table = pq.read_table(parquet_path)
        df = table.to_pandas()
        logging.info(f"Parquet file for pass {pass_number} loaded successfully.")

        # Step 2: Reset the index to extract 'Spot' and 'CellType'
        df = df.reset_index()

        # Check which delimiter format is used (new ::: or legacy _)
        sample_index = df["index"].iloc[0]
        if ":::" in sample_index:
            # New format: spot_1:::CellType
            df[["Spot", "CellType"]] = df["index"].str.split(":::", n=1, expand=True)
            logging.info("Spot and CellType successfully split (using ::: delimiter).")
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
                        spots.append(idx_val[: -len(suffix)])
                        cell_types.append(ct)
                        matched = True
                        break
                if not matched:
                    raise ValueError(
                        f"Could not parse index '{idx_val}' - no matching cell type found in {known_cell_types}"
                    )

            df["Spot"] = spots
            df["CellType"] = cell_types
            logging.info("Spot and CellType successfully split (using legacy _ delimiter with cell type matching).")

        df = df.drop(columns=["index"])

        # Get cell type names from the dictionary for validation
        if self.cell_profile_dict is None:
            raise ValueError("Cell profile dictionary not loaded. Run load_cell_profile_dict() first.")

        expected_cell_types = set(self.cell_profile_dict.keys())
        found_cell_types = set(df["CellType"].unique())

        if not found_cell_types.issubset(expected_cell_types):
            logging.warning("Found unexpected cell types: %s", found_cell_types - expected_cell_types)
            logging.warning("Expected cell types: %s", expected_cell_types)
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
            elif celltype_data["Spot"].str.contains("spot_").all() and "spot_" not in str(
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
            logging.info(f"Added layer: {layer_name} (Shape: {celltype_matrix.shape})")

            # After adding each layer, verify it was added correctly
            if layer_name not in self.gene_expression_adata.layers:
                logging.error("Failed to add layer: %s", layer_name)
            else:
                logging.info("Successfully added layer: %s", layer_name)

    def get_adata(self):
        """
        Retrieve the internal AnnData object for downstream analysis.

        Returns:
            AnnData: The internal `adata` object.
        """
        if self.gene_expression_adata is None:
            raise ValueError("AnnData object is not initialized in the model.")

        logging.info("✅ Returning the internal AnnData object.")
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

    def assign_cells(
        self,
        nuclei_counts: "pd.Series",
        cell_to_spot: np.ndarray,
        *,
        cell_ids: np.ndarray = None,
        random_state: int = 42,
        assignment_method: str = "hungarian",
        detection_mask: np.ndarray = None,
        proportion_prior: np.ndarray = None,
        morph_scores_precomputed: np.ndarray = None,
        morphology_weight: float = 0.5,
    ) -> "pd.DataFrame":
        """Assign individual cells to types using proportions.

        Post-processing step that runs after run_cell_proportion_model().
        Discretizes spot-level proportions into per-cell type assignments via
        Hungarian matching, or Bayesian posterior assignment when
        morph_scores_precomputed is supplied.

        Args:
            nuclei_counts: Per-spot nuclei counts.
            cell_to_spot: (C,) int array mapping each cell to a spot index
                (positional into self.results["cell_prop"]).
            cell_ids: (C,) cell/nucleus identifiers. Defaults to integer indices.
            random_state: Seed for deterministic no-morphology assignment.
            assignment_method: "hungarian" or "bayesian". Default "hungarian".
            detection_mask: (I, T) boolean mask required for bayesian assignment.
            proportion_prior: (I, T) alternative prior for bayesian assignment.
            morph_scores_precomputed: (C, T) precomputed morphology scores for bayesian.
            morphology_weight: Hungarian morphology nudge in [0,1]; ignored by bayesian.

        Returns:
            DataFrame (C rows): spot_id, cell_id, per-type scores, assigned_type, confidence.
        """
        from .assignment.cell_assignment import assign_cells as _assign_cells  # pylint: disable=import-outside-toplevel

        if "cell_prop" not in self.results:
            raise RuntimeError(
                "run_cell_proportion_model() must be called before assign_cells(). " "No 'cell_prop' found in results."
            )

        # Auto-populate for Bayesian assignment
        if assignment_method == "bayesian":
            if detection_mask is None:
                detection_mask = self.results.get("detection_mask_bool")
                if detection_mask is None:
                    raise RuntimeError(
                        "No detection mask available. Run run_cell_proportion_model() "
                        "with use_detection_gating=True first."
                    )
            if proportion_prior is None:
                base = self.results.get("cell_prop_base", self.results.get("cell_prop"))
                if base is not None:
                    proportion_prior = base.values if hasattr(base, "values") else base

        existing_int = self.results.get("cell_counts_integer", None)

        result = _assign_cells(
            cell_prop=self.results["cell_prop"],
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
            cell_ids=cell_ids,
            existing_integer_counts=existing_int,
            output_folder=getattr(self, "output_folder", None),
            sample_name=getattr(self, "sample_name", "sample"),
            random_state=random_state,
            assignment_method=assignment_method,
            detection_mask=detection_mask,
            proportion_prior=proportion_prior,
            morph_scores_precomputed=morph_scores_precomputed,
            morphology_weight=morphology_weight,
        )

        self.results["cell_assignments"] = result
        self.results["cell_types_hard"] = result["assigned_type"].values

        logging.info("Assigned %d cells to %d types", len(result), len(self.results["cell_prop"].columns))

        return result

    @staticmethod
    def _validate_he_resolution(source_adata, resolution_mode: str):
        """Validate H&E image resolution for single-cell segmentation."""
        from .morphology.segmentation import (
            _get_first_library_payload,
            estimate_pixel_size_um,
            get_resolution_image_and_scale,
        )

        MIN_SEGMENTATION_UM_PER_PX = 1.0

        try:
            _, scale = get_resolution_image_and_scale(source_adata, resolution_mode)
        except ValueError as e:
            raise ValueError(
                f"Single-cell output requires a high-resolution H&E image for "
                f"nuclei segmentation. No '{resolution_mode}' image found in "
                f"uns['spatial']. Load data with sq.read.visium(load_images=True) "
                f"or provide a fullres image."
            ) from e

        if resolution_mode == "fullres":
            lib = _get_first_library_payload(source_adata)
            if "fullres" not in lib.get("images", {}):
                raise ValueError(
                    "No true fullres image found. The 'fullres' resolution mode "
                    "falls back to upscaling hires, which does not improve optical "
                    "resolution. Use resolution_mode='hires' or provide a real "
                    "fullres image."
                )

        pixel_size = estimate_pixel_size_um(source_adata)
        if pixel_size is None:
            raise ValueError(
                "Cannot determine image resolution from scalefactors. "
                "Single-cell output requires known pixel size. "
                "Ensure scalefactors['spot_diameter_fullres'] is present."
            )

        effective_um_per_px = pixel_size / scale
        if effective_um_per_px > MIN_SEGMENTATION_UM_PER_PX:
            raise ValueError(
                f"Image resolution too coarse for nuclei segmentation "
                f"({effective_um_per_px:.2f} µm/px in {resolution_mode} frame, "
                f"need ≤{MIN_SEGMENTATION_UM_PER_PX} µm/px). "
                f"Provide a hires or fullres H&E image."
            )

    def run_sace_allocation(
        self,
        *,
        output_mode: str = "layers",
        init_method: str = "marker",
        integer_allocation: bool = False,
        resolution_mode: str = "hires",
        modality: str = "he",
        assignment_method: str = "hungarian",
        stardist_kwargs: dict = None,
        use_morphology: bool = True,
        morphology_weight: float = 0.5,
    ):
        """Run SACE per-cell GEX allocation.

        Args:
            output_mode: "layers" (spot-level type compartments, no nuclei
                needed) or "single_cell" (auto-runs StarDist + assignment +
                per-cell projection; requires high-res H&E).
            init_method: "marker" (default) or "prop".
            integer_allocation: If True, apply largest-remainder rounding.
            resolution_mode: Image resolution for StarDist (single_cell only).
            modality: StarDist modality "he" or "dapi" (single_cell only).
            assignment_method: "hungarian" (default) or "bayesian" (single_cell only,
                requires morph_scores_precomputed via assign_cells).
            stardist_kwargs: Extra args for StarDist (single_cell only).
            use_morphology: compute non-DL morphology scores and assign with the
                Bayesian posterior (the benchmarked winner: 0.434 pooled vs 0.336
                Hungarian on Xenium RCC); default True. If morphology scoring
                fails, falls back to count-constrained Hungarian. Set False for
                byte-identical, proportion-consistent Hungarian assignment.
            morphology_weight: Hungarian morphology nudge in [0,1], forwarded to
                assign_cells; ignored by bayesian.

        Returns:
            layers mode: (spot_type_gex, diagnostics)
            single_cell mode: (spot_type_gex, cell_adata, diagnostics)
        """
        if output_mode not in ("layers", "single_cell"):
            raise ValueError(f"output_mode must be 'layers' or 'single_cell', got '{output_mode}'")

        from .gex.sace_gex import integerize_spot_type_gex, project_sace_to_cells, run_sace_layers

        if "cell_prop" not in self.results:
            raise ValueError("Run run_cell_proportion_model() first.")

        # Fail fast: validate H&E before any computation in single_cell mode
        if output_mode == "single_cell":
            source_adata = self.antibody_capture_adata or self.gene_expression_adata or self.adata
            self._validate_he_resolution(source_adata, resolution_mode)

        adata = self.gene_expression_adata
        if "raw_counts" in adata.layers:
            raw_X = adata.layers["raw_counts"]
        else:
            logging.warning("No raw_counts layer found — using .X (may be normalized)")
            raw_X = adata.X

        ab_data = None
        ab_names = None
        if hasattr(self, "antibody_capture_adata") and self.antibody_capture_adata is not None:
            ab_adata = self.antibody_capture_adata
            ab_names = list(ab_adata.var_names)
            gex_spots = set(adata.obs_names)
            ab_mask = ab_adata.obs_names.isin(gex_spots)
            if ab_mask.sum() < len(adata):
                logging.warning(
                    "Only %d/%d GEX spots found in antibody data — " "skipping marker-guided SACE init",
                    ab_mask.sum(),
                    len(adata),
                )
            else:
                ab_sub = ab_adata[ab_mask]
                ab_sub = ab_sub[adata.obs_names]
                ab_raw = ab_sub.layers.get("raw_counts", ab_sub.X)
                ab_data = np.asarray(ab_raw)
                logging.info(
                    "Aligned antibody data: %d spots × %d markers",
                    ab_data.shape[0],
                    ab_data.shape[1],
                )

        spot_type_gex, E_x, internals, diagnostics = run_sace_layers(
            spot_counts=raw_X,
            proportions=self.results["cell_prop"],
            spot_coords=adata.obsm["spatial"],
            gene_names=list(adata.var_names),
            spotwise_profiles_init=self.results.get("gene_expression_pass1"),
            antibody_data=ab_data,
            antibody_names=ab_names,
            cell_profile_dict=getattr(self, "cell_profile_dict", None),
            init_method=init_method,
        )

        if integer_allocation:
            spot_type_gex = integerize_spot_type_gex(spot_type_gex, raw_X)
            logging.info("Applied LR integer allocation (conserving spot totals)")

        self.results["sace_spot_type_gex"] = spot_type_gex
        self.results["sace_diagnostics"] = diagnostics

        if output_mode == "layers":
            return spot_type_gex, diagnostics

        # --- single_cell mode ---
        source_adata = self.antibody_capture_adata or self.gene_expression_adata or self.adata

        from .morphology.segmentation import (
            compute_spot_nuclei_counts_from_adata,
            detect_spot_diameter_pixels,
            get_resolution_image_and_scale,
            normalize_nuclei_counts_for_prior,
            save_segmentation_artifacts,
        )
        from .morphology.segmentation_qc import generate_segmentation_qc_pdf

        seg_result = compute_spot_nuclei_counts_from_adata(
            adata=source_adata,
            resolution_mode=resolution_mode,
            modality=modality,
            **(stardist_kwargs or {}),
        )
        self.results["segmentation_result"] = seg_result

        save_segmentation_artifacts(
            output_folder=self.output_folder,
            sample_name=self.sample_name,
            result=seg_result,
            save_masks=True,
        )

        image, scale = get_resolution_image_and_scale(source_adata, resolution_mode)
        spot_centers_fullres = np.asarray(source_adata.obsm["spatial"], dtype=np.float64)
        spot_centers_image = spot_centers_fullres * float(scale)
        spot_diam_fullres = detect_spot_diameter_pixels(adata=source_adata)
        spot_radius_image = (spot_diam_fullres * float(scale)) / 2.0

        generate_segmentation_qc_pdf(
            output_folder=self.output_folder,
            sample_name=self.sample_name,
            image=image,
            masks=seg_result.masks,
            centroids_xy=seg_result.centroids_xy,
            spot_centers_xy=spot_centers_image,
            spot_radius_px=spot_radius_image,
            nuclei_count_raw=seg_result.nuclei_count_raw,
            qc=seg_result.qc,
        )

        nucleus_spot_map = seg_result.nucleus_spot_map
        if nucleus_spot_map is None or len(nucleus_spot_map) == 0:
            raise ValueError(
                "StarDist segmentation found no nuclei assigned to spots. "
                "Check that the H&E image covers the tissue area."
            )

        nuclei_counts = seg_result.nuclei_count_raw
        target_counts = normalize_nuclei_counts_for_prior(nuclei_counts)

        for ad in (self.gene_expression_adata, self.antibody_capture_adata, self.adata):
            if ad is None:
                continue
            common = ad.obs_names.intersection(nuclei_counts.index)
            if len(common) == 0:
                continue
            ad.obs.loc[common, "nuclei_count_raw"] = nuclei_counts.loc[common].astype(float).values
            ad.obs.loc[common, "nuclei_count_target"] = target_counts.loc[common].astype(float).values

        spot_names = list(adata.obs_names)
        spot_to_idx = {s: i for i, s in enumerate(spot_names)}
        nsm = nucleus_spot_map.copy()
        nsm["spot_idx"] = nsm["spot_id"].map(spot_to_idx)
        nsm = nsm.dropna(subset=["spot_idx"])
        nsm["spot_idx"] = nsm["spot_idx"].astype(int)
        cell_to_spot = nsm["spot_idx"].values
        cell_ids = nsm["nucleus_id"].values.astype(str)

        morph_scores_precomputed = None
        sc_method = assignment_method
        if use_morphology:
            from .morphology.nucleus_scores import compute_morphology_scores_safe

            morph_scores_precomputed = compute_morphology_scores_safe(
                labeled_mask=seg_result.masks,
                nuclei_spot_map=cell_to_spot,
                spot_proportions=self.results["cell_prop"].values,
                cell_type_names=list(self.results["cell_prop"].columns),
                nucleus_ids=nsm["nucleus_id"].values.astype(int),
            )
            if morph_scores_precomputed is not None:
                sc_method = "bayesian"
            else:
                logging.warning(
                    "Morphology-informed Bayesian assignment unavailable; using "
                    "count-constrained Hungarian assignment instead."
                )

        assignment_result = self.assign_cells(
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
            cell_ids=cell_ids,
            assignment_method=sc_method,
            morph_scores_precomputed=morph_scores_precomputed,
            morphology_weight=morphology_weight,
        )

        cell_assignments_dict = dict(
            zip(
                assignment_result["cell_id"].astype(str),
                assignment_result["assigned_type"],
            )
        )

        cell_spot_map = pd.DataFrame(
            {
                "cell_id": cell_ids,
                "spot_barcode": nsm["spot_id"].values,
                "spot_idx": cell_to_spot,
                "x": nsm["x_pixel"].values.astype(float),
                "y": nsm["y_pixel"].values.astype(float),
            }
        )

        type_names = list(self.results["cell_prop"].columns)
        cell_adata = project_sace_to_cells(
            E_x,
            internals,
            cell_assignments_dict,
            cell_spot_map,
            list(adata.var_names),
            type_names,
        )

        self.results["sace_cell_adata"] = cell_adata
        logging.info(
            "Single-cell SACE: %d cells × %d genes",
            cell_adata.n_obs,
            cell_adata.n_vars,
        )

        return spot_type_gex, cell_adata, diagnostics

    def run_module3_5_functional_annotation(
        self,
        *,
        functional_marker_table=None,
        max_iter=200,
        lr=0.01,
        early_stopping_patience=20,
        holdout_fraction=0.1,
        device="cpu",
        lambda_sigma=2.0,
        gmm_min_proportion=0.05,
        min_spots=20,
        m1_result=None,
        m2_result=None,
        profile_discovery_result=None,
        coverage_threshold=0.75,
    ):
        """Module 3.5: Functional protein annotation.

        Learns per-type emission rates for functional markers via NB MLE
        (frozen QP proportions), then GMM-gates to call positive/negative
        states using observed-to-expected ratio.

        Requires run_cell_proportion_model() to have been called first.

        Args:
            functional_marker_table: Dict of marker -> {function, active_types}.
                Defaults to DEFAULT_FUNCTIONAL_TABLE.
            max_iter: Max NB optimization iterations.
            lr: Adam learning rate.
            early_stopping_patience: Early stopping patience.
            holdout_fraction: Fraction held out for early stopping.
            device: PyTorch device.
            lambda_sigma: Prior sigma for active lambda pairs.
            gmm_min_proportion: Min cell-type proportion for GMM inclusion.
            min_spots: Min spots for a (type,marker) pair to be estimated.

        Returns:
            Dict with functional_lambda, functional_intensity,
            functional_gates, functional_summary.
        """
        from .annotation.functional_annotation import (  # pylint: disable=import-outside-toplevel
            DEFAULT_FUNCTIONAL_TABLE,
            build_active_mask,
            compute_spatial_statistics,
            gate_functional_markers,
            learn_functional_emissions,
        )

        if "cell_prop" not in self.results:
            raise ValueError("Run run_cell_proportion_model() first.")

        if functional_marker_table is None:
            functional_marker_table = DEFAULT_FUNCTIONAL_TABLE

        proportions = self.results["cell_prop"]
        cell_types = list(proportions.columns)

        # Find functional markers present in the antibody panel
        ab_adata = self.antibody_capture_adata
        panel_markers = list(ab_adata.var_names)

        # Try direct match first, then try stripping -1 suffix
        functional_markers = [m for m in functional_marker_table if m in panel_markers]
        marker_to_panel = {m: m for m in functional_markers}

        if len(functional_markers) == 0:
            panel_stripped = {}
            for pm in panel_markers:
                stripped = strip_antibody_suffix(pm)
                panel_stripped[stripped] = pm
            functional_markers = [m for m in functional_marker_table if m in panel_stripped]
            marker_to_panel = {m: panel_stripped[m] for m in functional_markers}

        logging.info(
            "Module 3.5: %d/%d functional markers found in panel", len(functional_markers), len(functional_marker_table)
        )

        if len(functional_markers) == 0:
            raise ValueError(
                "No functional markers found in antibody panel. " f"Panel markers: {panel_markers[:10]}..."
            )

        # Extract raw counts for functional markers (NB needs integer counts)
        panel_names = [marker_to_panel[m] for m in functional_markers]
        if "raw_counts" in ab_adata.layers:
            ab_raw = ab_adata[:, panel_names].layers["raw_counts"]
        else:
            logging.warning("No raw_counts layer — using .X (may be CLR-normalized)")
            ab_raw = ab_adata[:, panel_names].X
        if hasattr(ab_raw, "toarray"):
            ab_raw = ab_raw.toarray()
        observed = np.asarray(ab_raw, dtype=np.float32)

        # Build active mask
        active_mask = build_active_mask(functional_markers, cell_types, functional_marker_table)

        # Size factors: fixed from antibody total counts
        total_ab = (
            np.asarray(ab_adata.X.sum(axis=1)).flatten() if hasattr(ab_adata.X, "toarray") else ab_adata.X.sum(axis=1)
        )
        total_ab = np.asarray(total_ab, dtype=np.float64).flatten()
        median_ab = np.median(total_ab[total_ab > 0]) if (total_ab > 0).any() else 1.0
        size_factors = np.clip(total_ab / max(median_ab, 1.0), 0.1, 10.0).astype(np.float32)

        # Learn emission rates
        nb_results = learn_functional_emissions(
            observed=observed,
            proportions=proportions.values,
            active_mask=active_mask,
            size_factors=size_factors,
            max_iter=max_iter,
            lr=lr,
            early_stopping_patience=early_stopping_patience,
            holdout_fraction=holdout_fraction,
            lambda_sigma=lambda_sigma,
            device=device,
        )

        # GMM gating (ratio-based)
        intensity_df, gates_df, summary = gate_functional_markers(
            observed=observed,
            proportions=proportions.values,
            lam=nb_results["lambda"],
            background=nb_results["background"],
            size_factors=size_factors,
            active_mask=active_mask,
            cell_types=cell_types,
            functional_markers=functional_markers,
            gmm_min_proportion=gmm_min_proportion,
            min_spots=min_spots,
        )
        # Assign spot barcodes as index so projection can join on spot_id
        gates_df.index = proportions.index
        intensity_df.index = proportions.index

        # Spatial statistics
        spot_coords = self.gene_expression_adata.obsm["spatial"]
        active_pairs = [
            (ct, m)
            for t_idx, ct in enumerate(cell_types)
            for m_idx, m in enumerate(functional_markers)
            if active_mask[t_idx, m_idx] > 0.5
        ]
        spatial_stats = compute_spatial_statistics(gates_df, spot_coords, active_pairs)

        # Merge spatial stats into summary
        for pair, stats in spatial_stats.items():
            if pair in summary:
                summary[pair].update(stats)

        # Store results
        lam_df = pd.DataFrame(nb_results["lambda"], index=cell_types, columns=functional_markers)
        self.results["functional_lambda"] = lam_df
        self.results["functional_intensity"] = intensity_df
        self.results["functional_gates"] = gates_df
        self.results["functional_summary"] = summary

        # Save to output folder
        if self.output_folder:
            lam_df.to_csv(os.path.join(self.output_folder, f"{self.sample_name}_functional_lambda.csv"))
            intensity_df.to_csv(
                os.path.join(self.output_folder, f"{self.sample_name}_functional_intensity.csv"), index=False
            )
            gates_df.to_csv(os.path.join(self.output_folder, f"{self.sample_name}_functional_gates.csv"), index=False)

            import json  # pylint: disable=import-outside-toplevel

            summary_ser = {f"{k[0]}:{k[1]}": v for k, v in summary.items()}
            with open(os.path.join(self.output_folder, f"{self.sample_name}_functional_summary.json"), "w") as f:
                json.dump(summary_ser, f, indent=2)

        n_sig = sum(1 for s in spatial_stats.values() if s.get("morans_p", 1.0) < 0.05)
        logging.info(
            "Module 3.5 complete: %d markers, %d active pairs, %d spatially significant",
            len(functional_markers),
            len(active_pairs),
            n_sig,
        )

        results = {
            "functional_lambda": lam_df,
            "functional_intensity": intensity_df,
            "functional_gates": gates_df,
            "functional_summary": summary,
        }

        # Optional coverage gap analysis (fires only when M1/M2 results are provided)
        if m1_result is not None or m2_result is not None:
            from .annotation.coverage_check import check_module_coverage  # pylint: disable=import-outside-toplevel

            coverage = check_module_coverage(
                m1_result=m1_result,
                m2_result=m2_result,
                cell_profile_dict=getattr(self, "cell_profile_dict", {}),
                functional_marker_table=functional_marker_table or {},
                profile_discovery_result=profile_discovery_result,
                colocalization_threshold=coverage_threshold,
            )
            for line in coverage.warning_lines:
                logging.warning(line)
            if coverage.n_warnings > 0 and self.output_folder:
                coverage.to_csv(self.output_folder)
                logging.info("Coverage check saved to %s/coverage_check_*.csv", self.output_folder)
            results["coverage_check"] = coverage

        return results

    def run_functional_annotation(self, *args, **kwargs):
        """Compatibility wrapper for the Module 3.5 functional annotation entrypoint."""
        return self.run_module3_5_functional_annotation(*args, **kwargs)

    def run_sace_protein(
        self,
        cell_assignments: Dict[str, str],
        cell_spot_map: "pd.DataFrame",
        *,
        module3_5_candidates_df: Optional["pd.DataFrame"] = None,
        functional_table: Optional[Dict] = None,
        bimodality_threshold: float = 1.5,
        posterior_threshold: float = 0.5,
        min_high_component_log_mean: float = 1.0,
    ):
        """Per-cell functional protein deconvolution via SACE.

        Allocates spot-level functional protein counts to individual cells
        using the same SACE machinery as GEX, then GMM-gates each cell.

        Preconditions:
            - run_cell_proportion_model() called (proportions available)
            - assign_cells() called (cell_assignments and cell_spot_map available)

        Args:
            cell_assignments: Dict[cell_id -> type_name] from assign_cells().
            cell_spot_map: DataFrame with cell_id, spot_barcode, spot_idx, x, y.
            module3_5_candidates_df: Optional DataFrame with parent_type and
                functional_marker columns for custom active mask.
            functional_table: Optional dict in DEFAULT_FUNCTIONAL_TABLE format
                (marker -> {"active_types": [...]}). If provided, replaces
                DEFAULT_FUNCTIONAL_TABLE. Use to pass M1-filtered per-sample tables.
            bimodality_threshold: Forwarded to gmm_gate_cells(); separation
                multiplier on pooled_std (default 1.5).
            posterior_threshold: Forwarded to gmm_gate_cells(); min posterior
                probability to call a cell positive (default 0.5).
            min_high_component_log_mean: Forwarded to gmm_gate_cells(); high
                component log-mean must exceed this or gate is suppressed
                (default 1.0, calibrated on Xenium RCC Region 0).

        Returns:
            Dict with:
                cell_protein: (N_cells, M_func) ndarray
                protein_names: list of marker names
                protein_gates_df: DataFrame with gate columns
                gmm_summary: per-pair gating summary
                sace_diagnostics: SACE convergence info
        """
        from .annotation.functional_annotation import (  # pylint: disable=import-outside-toplevel
            DEFAULT_FUNCTIONAL_TABLE,
            build_active_mask,
            gmm_gate_cells,
        )
        from .gex.sace_gex import run_sace  # pylint: disable=import-outside-toplevel

        if "cell_prop" not in self.results:
            raise ValueError("run_cell_proportion_model() must be called before run_sace_protein().")

        # --- 1. Identify unused functional markers ---
        if not hasattr(self, "cell_profile_dict") or self.cell_profile_dict is None:
            raise ValueError("cell_profile_dict must be loaded before run_sace_protein().")

        consumed = set()
        for ct_info in self.cell_profile_dict.values():
            for key in ("Major", "Soft"):
                consumed.update(ct_info.get(key, []))

        panel_markers = list(self.antibody_capture_adata.var_names)

        # Handle -1 suffix: "CD68-1" -> "CD68" for matching, but preserve
        # names that naturally end in "-1" like "PD-1"
        def _canon(name):
            return strip_antibody_suffix(name)

        panel_stripped = {_canon(m): m for m in panel_markers}

        unused_panel_names = []
        for canon, original in panel_stripped.items():
            if canon not in consumed:
                unused_panel_names.append(original)

        unused_canon = [_canon(m) for m in unused_panel_names]

        logging.info(
            "SACE protein: %d unused markers out of %d panel markers",
            len(unused_canon),
            len(panel_markers),
        )

        if len(unused_canon) == 0:
            logging.warning("All panel markers consumed by QP; no markers for SACE protein.")
            return {}

        # --- 2. Build active mask ---
        cell_types = list(self.results["cell_prop"].columns)

        if module3_5_candidates_df is not None:
            type_marker_pairs = []
            for _, row in module3_5_candidates_df.iterrows():
                ct = row["parent_type"]
                marker = row["functional_marker"]
                if _canon(marker) in unused_canon or marker in unused_canon:
                    type_marker_pairs.append((ct, _canon(marker)))
            active_mask = build_active_mask(
                unused_canon,
                cell_types,
                functional_table={
                    m: {"active_types": [ct for ct2, m2 in type_marker_pairs if m2 == m for ct in [ct2]]}
                    for m in unused_canon
                },
            )
        else:
            active_mask = build_active_mask(
                unused_canon,
                cell_types,
                functional_table if functional_table is not None else DEFAULT_FUNCTIONAL_TABLE,
            )

        logging.info(
            "SACE protein active mask: %d active pairs out of %d possible",
            int(active_mask.sum()),
            active_mask.size,
        )

        # --- 3. Extract raw protein counts ---
        ab_adata = self.antibody_capture_adata
        gex_spots = set(self.gene_expression_adata.obs_names)

        # Subset and reorder antibody to match GEX spots
        ab_mask = ab_adata.obs_names.isin(gex_spots)
        ab_sub = ab_adata[ab_mask]
        ab_sub = ab_sub[self.gene_expression_adata.obs_names]

        # Extract raw counts for unused markers
        try:
            col_idx = [list(ab_sub.var_names).index(m) for m in unused_panel_names]
        except ValueError as e:
            logging.error("SACE protein: marker not found in antibody panel — %s. " "Cannot extract protein counts.", e)
            return {}

        raw_layer = ab_sub.layers.get("raw_counts", ab_sub.X)
        raw = raw_layer[:, col_idx] if hasattr(raw_layer, "__getitem__") else raw_layer
        spot_protein = np.asarray(raw, dtype=np.float64)

        # --- 4. Adapt cell_spot_map columns ---
        csm = cell_spot_map.copy()
        if "nucleus_id" in csm.columns and "cell_id" not in csm.columns:
            csm = csm.rename(columns={"nucleus_id": "cell_id"})
        if "spot_id" in csm.columns and "spot_barcode" not in csm.columns:
            csm = csm.rename(columns={"spot_id": "spot_barcode"})
        if "x_pixel" in csm.columns and "x" not in csm.columns:
            csm = csm.rename(columns={"x_pixel": "x", "y_pixel": "y"})
        if "spot_idx" not in csm.columns:
            spot_names = list(self.gene_expression_adata.obs_names)
            spot_to_idx = {s: i for i, s in enumerate(spot_names)}
            csm["spot_idx"] = csm["spot_barcode"].map(spot_to_idx)
            csm = csm.dropna(subset=["spot_idx"])
            csm["spot_idx"] = csm["spot_idx"].astype(int)

        # --- 5. Build warm-start from functional_lambda if available ---
        spotwise_init = None
        prev_lam = self.results.get("functional_lambda")
        if prev_lam is not None:
            M = len(unused_canon)
            init_profile = np.zeros((len(cell_types), M), dtype=np.float64)
            for m_idx, canon in enumerate(unused_canon):
                if canon in prev_lam.columns:
                    init_profile[:, m_idx] = prev_lam[canon].values
            row_sums = init_profile.sum(axis=1, keepdims=True)
            row_sums = np.where(row_sums > 0, row_sums, 1.0)
            init_profile_norm = init_profile / row_sums

            N_spots = spot_protein.shape[0]
            spotwise_init = {s: init_profile_norm.copy() for s in range(N_spots)}
            logging.info("SACE protein: warm-starting from functional_lambda")

        # --- 6. Call run_sace() ---
        _, cell_adata, diagnostics = run_sace(
            spot_counts=spot_protein,
            proportions=self.results["cell_prop"],
            cell_assignments=cell_assignments,
            cell_spot_map=csm,
            spot_coords=self.gene_expression_adata.obsm["spatial"],
            gene_names=unused_canon,
            spotwise_profiles_init=spotwise_init,
        )

        # Extract per-cell protein from cell_adata.X
        cell_protein = np.asarray(cell_adata.X, dtype=np.float32)
        cell_type_arr = np.array(cell_adata.obs["cell_type"].values)

        # --- 7. Post-hoc active mask: zero out inactive (type, marker) pairs ---
        type_to_idx = {ct: i for i, ct in enumerate(cell_types)}
        type_indices = np.array([type_to_idx.get(ct, -1) for ct in cell_type_arr])
        valid = type_indices >= 0
        # Build per-cell mask: (N_cells, M) — active_mask is (T, M)
        row_idx = np.where(valid, type_indices, 0)
        cell_active = active_mask[row_idx]  # (N_cells, M)
        cell_active[~valid] = 1  # unknown types: preserve all values
        cell_protein = np.where(cell_active.astype(bool), cell_protein, np.nan)

        # --- 8. GMM gating ---
        gates_df, gmm_summary = gmm_gate_cells(
            cell_protein=cell_protein,
            cell_types=cell_type_arr,
            type_names=cell_types,
            marker_names=unused_canon,
            active_mask=active_mask,
            bimodality_threshold=bimodality_threshold,
            posterior_threshold=posterior_threshold,
            min_high_component_log_mean=min_high_component_log_mean,
        )
        # Align gate index to cell_adata index
        gates_df.index = cell_adata.obs.index

        # --- 9. Store results ---
        self.results["sace_cell_protein"] = cell_protein
        self.results["sace_protein_names"] = unused_canon
        self.results["sace_protein_gates"] = gates_df
        self.results["sace_protein_gmm_summary"] = gmm_summary
        self.results["sace_protein_diagnostics"] = diagnostics

        logging.info(
            "SACE protein complete: %d cells x %d markers, %d active pairs gated",
            cell_protein.shape[0],
            cell_protein.shape[1],
            len(gmm_summary),
        )

        # Save outputs if output_folder configured
        if self.output_folder:
            np.save(
                os.path.join(self.output_folder, f"{self.sample_name}_cell_protein.npy"),
                cell_protein,
            )
            gates_df.to_csv(
                os.path.join(self.output_folder, f"{self.sample_name}_protein_gates.csv"),
            )
            logging.info("SACE protein outputs saved to %s", self.output_folder)

        return {
            "cell_protein": cell_protein,
            "protein_names": unused_canon,
            "protein_gates_df": gates_df,
            "gmm_summary": gmm_summary,
            "sace_diagnostics": diagnostics,
        }

    def run_protein_subtype_split(
        self,
        cell_assignments: Dict[str, str],
        cell_spot_map: "pd.DataFrame",
        validated_pairs: Optional[List[Tuple[str, str]]] = None,
        min_subtype_cells: int = 50,
    ) -> Tuple[Dict[str, str], "pd.DataFrame"]:
        """Split cell types into protein-gate-defined subtypes for SACE GEX.

        Intended to run after ``run_sace_protein()`` and before ``run_sace()``
        (GEX deconvolution).  Splits cell assignments and spot-level proportions
        using per-cell binary gate calls from SACE protein allocation.

        Only pairs whose gates were fired as ``gmm_bimodal`` (not suppressed as
        weak_signal or insufficient) are eligible for splitting.

        Args:
            cell_assignments: Dict[cell_id → type_name] from assign_cells().
            cell_spot_map: DataFrame with cell_id and spot_id columns.
            validated_pairs: List of (type, marker) pairs to split.  Defaults to
                the biologically motivated exhaustion/activation pairs:
                T cells×{LAG-3,PD-1}, Macrophages×PD-L1.
                Only pairs with gmm_bimodal gates are actually applied.
            min_subtype_cells: Minimum cells in each subtype to proceed (default 50).

        Returns:
            Tuple of (updated_cell_assignments, updated_proportions_df).
            Also stores results in self.results["subtype_assignments"] and
            self.results["subtype_proportions"].

        Raises:
            ValueError: If run_sace_protein() has not been called.
        """
        from .annotation.subtype_splitting import split_by_protein_gates  # pylint: disable=import-outside-toplevel

        if "sace_protein_gates" not in self.results:
            raise ValueError("run_sace_protein() must be called before run_protein_subtype_split().")
        if "cell_prop" not in self.results:
            raise ValueError("run_cell_proportion_model() must be called before run_protein_subtype_split().")

        protein_gates_df = self.results["sace_protein_gates"]
        gmm_summary = self.results.get("sace_protein_gmm_summary", {})
        proportions = self.results["cell_prop"].copy()

        if validated_pairs is None:
            validated_pairs = [
                ("T cells", "LAG-3"),
                ("T cells", "PD-1"),
                ("Macrophages", "PD-L1"),
            ]

        # Filter to pairs that actually fired gmm_bimodal — don't split on weak signal
        bimodal_pairs = [
            (ct, mk)
            for ct, mk in validated_pairs
            if gmm_summary.get((ct, mk), {}).get("gating_method") == "gmm_bimodal"
        ]
        if not bimodal_pairs:
            logging.info(
                "run_protein_subtype_split: no gmm_bimodal pairs found among %s; "
                "returning original assignments unchanged.",
                validated_pairs,
            )
            return cell_assignments, proportions

        logging.info(
            "run_protein_subtype_split: splitting on %d bimodal pairs: %s",
            len(bimodal_pairs),
            bimodal_pairs,
        )

        updated_assignments, updated_proportions = split_by_protein_gates(
            cell_assignments=cell_assignments,
            protein_gates_df=protein_gates_df,
            proportions=proportions,
            cell_spot_map=cell_spot_map,
            validated_pairs=bimodal_pairs,
            min_subtype_cells=min_subtype_cells,
        )

        self.results["subtype_assignments"] = updated_assignments
        self.results["subtype_proportions"] = updated_proportions

        logging.info(
            "run_protein_subtype_split: %d types now; original %d",
            len(updated_proportions.columns),
            len(proportions.columns),
        )

        return updated_assignments, updated_proportions

    def build_validated_module3_5_annotations(self, _assignments_df=None, _benchmark_summary=None):
        """Build validated Module 3.5 annotations from SACE protein output.

        If SACE protein has been run, annotations are already per-cell in
        self.results['sace_protein_gates']. This method exists for backward
        compatibility and to apply benchmark validation gating.
        """
        gates = self.results.get("sace_protein_gates")
        if gates is not None:
            self.results["projected_functional_calls"] = gates
            return gates

        logging.warning(
            "build_validated_module3_5_annotations: no SACE protein output found. " "Run run_sace_protein() first."
        )
        return None

    def build_validated_functional_annotations(self, assignments_df, benchmark_summary):
        """Compatibility wrapper for the Module 3.5 projection entrypoint."""
        return self.build_validated_module3_5_annotations(assignments_df, benchmark_summary)
