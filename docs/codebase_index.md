<!-- codebase-index git:cfdae96b -->
# CITEgeist Codebase Index

## CITEgeist/model/

### CITEgeist/model/__init__.py
*CITEgeist model package - spatial transcriptomics deconvolution using CITE-seq.*

### CITEgeist/model/annotation/__init__.py

### CITEgeist/model/annotation/coverage_check.py
*Coverage gap analysis for Module 3 and Module 3.5 inputs.*
- `check_module_coverage(m1_result, m2_result, cell_profile_dict, functional_marker_table, profile_discovery_result=None, colocalization_threshold=0.75) -> 'CoverageCheckResult'` — Check for coverage gaps between M1/M2 outputs and Module 3/3.5 input dicts.
- `class CoverageCheckResult` — Coverage gap analysis results.
  - `to_csv(output_dir) -> None` — Write coverage_check_markers.csv and coverage_check_pairs.csv to output_dir.

### CITEgeist/model/annotation/functional_annotation.py
*Module 3.5: Functional protein annotation for CITEgeist.*
- `build_active_mask(functional_markers, cell_types, functional_table=None) -> np.ndarray` — Build binary (T, M_func) active mask from the functional table.
- `learn_functional_emissions(observed, proportions, active_mask, size_factors, max_iter=200, lr=0.01, early_stopping_patience=20, holdout_fraction=0.1, lambda_sigma=2.0, device='cpu', seed=42) -> Dict` — Learn per-type NB emission rates for functional markers via PyTorch.
- `gate_functional_markers(observed, proportions, lam, background, size_factors, active_mask, cell_types, functional_markers, gmm_min_proportion=0.05, min_spots=20) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]` — Gate spots as functionally active using observed-to-expected ratio GMM.
- `compute_spatial_statistics(gates_df, spot_coords, active_pairs) -> Dict[Tuple[str, str], Dict]` — Compute Moran's I spatial autocorrelation on gate calls.
- `gmm_gate_cells(cell_protein, cell_types, type_names, marker_names, active_mask, min_cells=20, bimodality_threshold=1.5, posterior_threshold=0.5, min_high_component_log_mean=1.0) -> Tuple[pd.DataFrame, Dict]` — Gate per-cell deconvolved protein levels via 2-component GMM.

### CITEgeist/model/annotation/module3_5_benchmark.py
- `build_gt_binary_calls(cell_table, cell_type, marker) -> pd.Series` — Build binary ground-truth calls using a 2-component GMM on log1p protein values.
- `score_pair_predictions(pair_id, cell_type, marker, gt_calls, pred_scores, pred_calls, n_supporting_spots, headline) -> PairBenchmarkResult` — Score a functional pair against ground-truth and predicted calls.
- `aggregate_module3_5_results(pairs, thresholds) -> dict[str, Any]`
- `score_spot_attribution(cell_type, marker, lambda_df, proportions_df, protein_df, mapping) -> dict | None` — Compute SACE-protein spot attribution accuracy in GT-positive spots.
- `class PairBenchmarkResult`

### CITEgeist/model/annotation/module3_5_projection.py
- `should_enrich_module3_5_output(summary) -> bool`
- `normalized_validated_pairs(summary) -> set[str]`

### CITEgeist/model/annotation/subtype_splitting.py
*Gate identifiability audit, subtype proportion splitting, and per-cell protein-gate splitting.*
- `audit_gate_identifiability(parent_props, gates, type_names, func_marker_names, config, collinearity_threshold=0.8, min_effective_rank=8) -> IdentifiabilityReport` — Run Phase 0 identifiability checks on gate-split subtype proportions.
- `build_subtype_proportions(parent_props, gates, type_names, func_marker_names, config) -> pd.DataFrame` — Split parent proportions into subtypes using functional gates.
- `permute_gates_within_type(gates, rng) -> np.ndarray` — Shuffle gate values across spots independently per (type, marker) pair.
- `split_by_protein_gates(cell_assignments, protein_gates_df, proportions, cell_spot_map, validated_pairs, min_subtype_cells=50) -> Tuple[Dict[str, str], pd.DataFrame]` — Split cell types into functional subtypes using per-cell protein gates.
- `class IdentifiabilityReport`

### CITEgeist/model/assignment/__init__.py

### CITEgeist/model/assignment/cell_assignment.py
*Discrete cell assignment with optional morphology nudge.*
- `discretize_proportions(cell_prop, nuclei_counts, existing_integer_counts=None) -> pd.DataFrame` — Discretize continuous proportions to integer cell counts per spot.
- `hungarian_assign_spot(n_cells, integer_counts, cell_types, morph_scores, morphology_weight) -> List[str]` — Assign n_cells to types for a single spot using Hungarian matching.
- `assign_cells_to_types(integer_counts, cell_to_spot, cell_types, morph_scores, morphology_weight, spot_proportions, cell_ids=None) -> pd.DataFrame` — Assign every cell to a cell type via per-spot Hungarian matching.
- `fit_morphology_scores(embeddings, cell_to_spot, oracle_props, num_types, n_epochs=100, device='cpu', seed=42) -> np.ndarray` — Fit prototype-contrastive LLP on embeddings and return per-cell scores.
- `bayesian_assign_cells(morph_scores, cell_to_spot, proportion_prior, detection_mask, cell_ids=None, spot_ids=None, cell_types=None, eps=1e-10) -> pd.DataFrame` — Assign cells to types via Bayesian posterior with protein detection mask.
- `assign_cells(cell_prop, nuclei_counts, cell_to_spot, cell_ids=None, morphology_embeddings=None, patches=None, encoder_checkpoint=None, morphology_weight=0.5, existing_integer_counts=None, output_folder=None, sample_name='sample', n_morph_epochs=100, device='cpu', random_state=42, assignment_method='hungarian', detection_mask=None, proportion_prior=None, morph_scores_precomputed=None) -> pd.DataFrame` — Assign individual cells to types using proportions + optional morphology.
- `extract_embeddings(patches, encoder_checkpoint, output_folder, sample_name, device='cpu') -> np.ndarray` — Extract morphology embeddings from DAPI patches with hash-based caching.

### CITEgeist/model/assignment/cellularity_utils.py
*Utilities for cellularity-scaled QP deconvolution.*
- `round_counts_largest_remainder(c, N) -> np.ndarray` — Round continuous cell counts to integers summing to N.
- `winsorize_cellularity(N_raw, percentile=95) -> np.ndarray` — Winsorize nuclei counts to cap outliers and clamp zeros.
- `prepare_cellularity(nuclei_counts, spot_names, percentile=95) -> np.ndarray` — Prepare cellularity array from optional nuclei counts.

### CITEgeist/model/assignment/hungarian_assignment.py
*Hungarian algorithm for optimal nucleus-to-celltype assignment.*
- `assign_nuclei_to_types(probs, counts, nucleus_ids, lambda_prior=0.0, proportions=None) -> Dict[int, int]` — Assign nuclei to cell types using Hungarian algorithm.

### CITEgeist/model/assignment/module3b_nucleus_assignment.py
*Module 3b: Per-Nucleus Cell Type Assignment.*
- `random_assign_nuclei_to_types(nucleus_ids, counts, rng=None) -> Dict[int, int]` — Randomly assign nuclei to cell types respecting count constraints.
- `run_nucleus_assignment(mask, nuclei_spot_map, proportions, nuclei_counts, cell_types, cell_mask=None, use_morphology=False, random_seed=None) -> NucleusAssignmentResult` — Run full nucleus assignment pipeline.
- `class NucleusAssignmentResult` — Result of nucleus assignment.

### CITEgeist/model/assignment/single_cell_output.py
*Create single-cell AnnData output.*
- `create_single_cell_adata(cell_gex, morphology_features, assignments, sample_name, classifier=None, gene_metadata=None, assignment_method='hungarian', functional_annotations=None, functional_metadata=None, cell_protein=None, protein_names=None, protein_gates=None) -> ad.AnnData` — Create single-cell AnnData from cell-level expression and metadata.

### CITEgeist/model/checkpoints.py
*Checkpoint management for saving and resuming optimization state.*
- `class CheckpointManager` — Manages loading and saving of optimization checkpoints.
  - `check_complete_run(N, T, M)` — Check if a complete run exists.
  - `load_latest_checkpoint(N, T, M)` — Load the latest valid checkpoint.
  - `save_checkpoint(completed_spots, spotwise_profiles, N, T, M)` — Save current progress as checkpoint.
  - `save_final_results(spotwise_profiles, completed_spots, N, T, M)` — Save final results.

### CITEgeist/model/citegeist_model.py
*Main CitegeistModel class for spatial transcriptomics deconvolution.*
- `class CitegeistModel`
  - `register_gurobi(license_file_path)` — Configure Gurobi by setting only the license file path.
  - `split_adata()` — Split the AnnData object into separate gene expression and antibody capture sub-objects
  - `[static] winsorize(matrix, lower_percentile=5, upper_percentile=95)` — Winsorize a 2D NumPy array.
  - `[static] row_normalize(matrix, target_sum=10000.0)` — Row normalize a 2D NumPy array to a fixed target sum.
  - `[static] global_clr(matrix, epsilon=1e-06)` — Apply margin=2 CLR normalization (global geometric mean per marker).
  - `load_cell_profile_dict(cell_profile_dict)` — Load and validate the cell profile dictionary.
  - `filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=10)` — Filter genes in the gene expression AnnData object based on user-defined criteria.
  - `copy_gex_to_protein_adata()` — Copy the number of spots in the gene expression AnnData object to the antibody capture AnnData object.
  - `preprocess_gex(target_sum=10000)` — Preprocess gene expression data with count-preserving normalization.
  - `preprocess_antibody()` — Preprocess antibody capture data:
  - `preprocess_antibody_discrete(winsorize_lower=5, winsorize_upper=95, scale_per_marker=True, scale_mode='per_marker') -> None` — Preprocess antibody data for discrete cell assignment.
  - `compute_spot_nuclei_counts(resolution_mode='hires', max_fullres_side=9000, save_masks=True, modality='he', **stardist_kwargs) -> pd.Series` — Compute per-spot nuclei counts from Visium histology using StarDist.
  - `train_reference(reference_adata, cell_type_col='cell_type', n_markers_per_type=20, de_method='wilcoxon', type_mapping=None) -> 'ReferenceProfile'` — Train NB reference profiles from annotated scRNA-seq.
  - `run_cell_proportion_model(radius=None, tolerance=0.0001, max_iterations=20, lambda_reg=1, alpha=0.5, _max_y_change=0.4, _max_workers=None, _checkpoint_interval=100, unknown_threshold=0.05, min_celltype_threshold=0.01, redundancy_threshold=0.1, validation_warn_only=False, lambda_laplacian=0.03, laplacian_k=8, beta_min=0.1, beta_max=2.0, lambda_coverage=1.0, use_nuclei_prior=False, nuclei_prior_lambda=1.0, nuclei_target_col='nuclei_count_target', nuclei_prior_bounds=None, use_cellularity_laplacian=False, cellularity_col=None, cellularity_sigma=0.5, sparsity_mask=None, use_detection_gating=True, detection_gate_ub=0.01, use_gex_detection=True, gex_detection_k=10, gex_detection_min_corr=0.15, fusion_mode='union', refine_sparsity=False, refine_suppress_threshold=0.02, refine_rescue_threshold=0.08, use_entropy_weights=False, entropy_weight_alpha=1.0, marker_weight_dict=None, lambda_confusion=0.0, confusion_pairs_manual=None, confusion_corr_threshold=-0.25, nuclei_counts=None, _cellularity_slack=0.3, _lambda_cellularity=1.0, use_gating=None, priority_dict=None, threshold_method='auto', use_negative_gates=False, method='qp', sigma_scale=1.0, nb_device='cpu', nb_n_iter=100, nb_use_detection=True, nb_detection_threshold=0.5, nb_gpu_adam_steps=200, morphology_prior=None, morphology_patches=None, morphology_cell_to_spot=None, lambda_morphology=0.1, morphology_device='cpu', reference=None, kappa=10.0)` — Orchestrates the cell proportion optimization workflow.
  - `discretize_proportions(proportions_df, nuclei_counts) -> pd.DataFrame` — Convert continuous cell type proportions to integer cell counts.
  - `run_cell_expression_pass1(cell_assignments=None, cell_spot_map=None, sace_max_iter=1, sace_n_0=10.0, sace_bandwidth=None)` — Run gene expression deconvolution using SACE (Spatially-Adaptive Compositional EM).
  - `compute_expression_prior(spotwise_profiles_pass1, cell_type_numbers_array, lambda_prior=1.0, min_expression_threshold=0.1) -> Dict[str, Any]` — Compute global prior from pass 1 results.
  - `append_proportions_to_adata(proportions_path=None, key='finetuned')` — Append cell type proportions to AnnData object.
  - `append_gex_to_adata(parquet_path=None, pass_number=1)` — Append gene expression layers from a Parquet file back into the gene_expression_adata object.
  - `get_adata()` — Retrieve the internal AnnData object for downstream analysis.
  - `cleanup()` — Free memory and clean up temporary data.
  - `validate_neighborhood_size(radius)`
  - `assign_cells(nuclei_counts, cell_to_spot, cell_ids=None, morphology_embeddings=None, patches=None, encoder_checkpoint=None, morphology_weight=0.5, random_state=42, device='cpu', assignment_method='hungarian', detection_mask=None, proportion_prior=None, morph_scores_precomputed=None) -> 'pd.DataFrame'` — Assign individual cells to types using proportions + optional morphology.
  - `run_sace_allocation(cell_assignments, cell_spot_map, n_0=10.0, bandwidth=None, max_iter=1, tol=0.0001)` — Run SACE per-cell GEX allocation.
  - `run_module3_5_functional_annotation(functional_marker_table=None, max_iter=200, lr=0.01, early_stopping_patience=20, holdout_fraction=0.1, device='cpu', lambda_sigma=2.0, gmm_min_proportion=0.05, min_spots=20, m1_result=None, m2_result=None, profile_discovery_result=None, coverage_threshold=0.75)` — Module 3.5: Functional protein annotation.
  - `run_functional_annotation(*args, **kwargs)` — Compatibility wrapper for the Module 3.5 functional annotation entrypoint.
  - `run_sace_protein(cell_assignments, cell_spot_map, module3_5_candidates_df=None, functional_table=None, max_iter=1, n_0=10.0, bandwidth=None, bimodality_threshold=1.5, posterior_threshold=0.5, min_high_component_log_mean=1.0)` — Per-cell functional protein deconvolution via SACE.
  - `run_protein_subtype_split(cell_assignments, cell_spot_map, validated_pairs=None, min_subtype_cells=50) -> Tuple[Dict[str, str], 'pd.DataFrame']` — Split cell types into protein-gate-defined subtypes for SACE GEX.
  - `build_validated_module3_5_annotations(_assignments_df=None, _benchmark_summary=None)` — Build validated Module 3.5 annotations from SACE protein output.
  - `build_validated_functional_annotations(assignments_df, benchmark_summary)` — Compatibility wrapper for the Module 3.5 projection entrypoint.

### CITEgeist/model/deconvolution/__init__.py

### CITEgeist/model/deconvolution/cuopt_impl.py
*cuOPT-based optimization backend for CITEgeist.*
- `build_spatial_laplacian(coords, k=8, normed=True, cellularity=None, cellularity_sigma=0.5) -> scipy.sparse.spmatrix` — Build graph Laplacian from spatial coordinates using k-NN.
- `compute_global_prior(spotwise_gene_expression_profiles, cell_type_numbers_array, lambda_prior=1.0, min_expression_threshold=0.1) -> Dict[str, Any]` — Compute global prior from pass 1 results using normalized expression patterns.
- `map_antibodies_to_profiles(adata, cell_profile_dict)` — Map antibody capture data to predefined cell type profiles.
- `map_antibodies_to_profiles_v2(adata, cell_profile_dict)` — Map antibody capture data while preserving marker-level granularity.
- `map_antibodies_to_profiles_cellularity(adata, cell_profile_dict, clip_percentile=99, scale='median')` — Map antibody data for cellularity-scaled QP, preserving between-spot amplitude.
- `map_antibodies_raw_counts(adata, cell_profile_dict, winsorize_lower=1.0, winsorize_upper=99.0, scale='median', eps=1e-06)` — Map raw antibody counts for cellularity-scaled QP (no CLR).
- `compute_marker_exclusivity(marker_level_data, Y_values, marker_owners, _assignment_matrix, floor=0.3, epsilon=1e-09) -> np.ndarray` — Compute per-marker exclusivity scores measuring discriminative power.
- `validate_cell_proportions(Y_values, cell_type_names, profile_based_antibody_data=None, unknown_threshold=0.05, min_celltype_threshold=0.01, redundancy_threshold=0.2, warn_only=False) -> None` — Validate cell type proportions after Stage 1 optimization.
- `compute_proportion_enrichment(gene_expr, cell_type_props, celltype_frequencies=None) -> np.ndarray` — Compute proportion-based enrichment WITHOUT smoothing.
- `compute_marker_enrichment(gene_expr, anchor_expr, anchor_weights) -> np.ndarray` — Compute marker-guided enrichment via correlation with anchor genes.
- `compute_adaptive_enrichment(prop_enrichment, marker_enrichment, max_variance=0.1) -> np.ndarray` — Adaptively blend proportion and marker enrichment based on proportion variance.
- `precompute_anchor_expression(gene_expression, anchor_genes, anchor_weights) -> Tuple[np.ndarray, np.ndarray]` — Precompute weighted mean anchor expression per cell type.
- `estimate_true_expression_cell(X_obs, Y_assignments, coords, enrichment_weights, library_slack=1.5, lambda_enrich=1.0, lambda_spatial=0.0, spatial_k=50, max_workers=None, _checkpoint_interval=500) -> np.ndarray` — Estimate true gene expression per cell using optimization.
- `normalize_counts(adata, target_sum=10000, exclude_highly_expressed=False, max_fraction=0.05)` — Normalize counts while preserving integer values and relative proportions.
- `deconvolute_spot_with_neighbors_with_prior(spot_idx, adata, cell_type_numbers_array, radius, global_prior=None, lambda_prior_weight=0.0, local_enrichment_weight=0.5, global_enrichment_weight=0.5, continuous_relaxation=True, lambda_gex_reg=0.01, enrichment_smoothing=0.2, use_kl_regularization=True, kl_temperature=0.3, lambda_kl=0.1) -> Optional[np.ndarray]` — Deconvolute a spot with its neighbors, using enrichment weights and optional prior.
- `optimize_gene_expression(sample_name, deconvolution_expression_data, cell_type_numbers_array, filtered_adata, radius=2, global_enrichment_weight=0.5, local_enrichment_weight=0.5, global_prior=None, lambda_prior_weight=0.0, max_workers=None, checkpoint_interval=100, output_dir='checkpoints', rerun=False, continuous_relaxation=True, lambda_gex_reg=0.01, enrichment_smoothing=0.2, _anchor_genes=None, _anchor_weights=None, _module_weight=0.5, use_kl_regularization=True, kl_temperature=0.3, lambda_kl=0.1) -> Dict[str, Any]` — Optimize gene expression: prep-cook pattern (CPU-parallel prep + sequential GPU solve).
- `optimize_cell_proportions_per_marker(marker_level_data, marker_names, assignment_matrix, cell_type_names, tolerance=0.0001, max_iterations=50, lambda_reg=1.0, alpha=0.5, normalize_beta=True, beta_min=0.1, beta_max=2.0, unknown_threshold=0.05, min_celltype_threshold=0.01, redundancy_threshold=0.2, warn_only=False, lambda_laplacian=0.0, coords=None, laplacian_k=8, lambda_sparse=0.0, alpha_max=0.8, lambda_alpha=1.0, lambda_coverage=1.0, spot_abundance_target=None, lambda_abundance_prior=0.0, row_sum_bounds=None, cellularity=None, cellularity_sigma=0.5, sparsity_mask=None, spot_weights=None, morphology_prior=None, lambda_morphology=0.0, freeze_globals=False, marker_weight=None, _confusion_pairs=None, _lambda_confusion=0.0) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray, np.ndarray]` — Perform EM-based optimization for cell type proportions with per-marker beta.
- `optimize_cell_proportions_per_type_beta(marker_level_data, marker_names, cell_type_names, beta_init, prior_sigma, tolerance=0.0001, max_iterations=50, lambda_reg=1.0, alpha_elastic=0.5, beta_max=2.0, alpha_max=0.8, lambda_alpha=1.0, lambda_laplacian=0.0, coords=None, laplacian_k=8, row_sum_bounds=None, sparsity_mask=None, spot_weights=None, marker_weight=None) -> Tuple[np.ndarray, np.ndarray, Dict[str, Dict[str, float]], np.ndarray, np.ndarray, List[float]]` — EM-based optimization with per-type-per-marker beta matrix.
- `optimize_cell_proportions_per_marker_matrix(marker_level_data, marker_names, assignment_matrix, cell_type_names, tolerance=0.0001, max_iterations=50, lambda_reg=1.0, alpha=0.5, normalize_beta=True, beta_min=0.1, beta_max=2.0, unknown_threshold=0.05, min_celltype_threshold=0.01, redundancy_threshold=0.2, warn_only=False, lambda_laplacian=0.0, coords=None, laplacian_k=8, lambda_sparse=0.0, alpha_max=0.8, lambda_alpha=1.0, lambda_coverage=1.0, spot_abundance_target=None, lambda_abundance_prior=0.0, row_sum_bounds=None, cellularity=None, cellularity_sigma=0.5, sparsity_mask=None, spot_weights=None, marker_weight=None, confusion_pairs=None, lambda_confusion=0.0) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray, np.ndarray]` — Matrix-formulated version of optimize_cell_proportions_per_marker.

### CITEgeist/model/deconvolution/detection.py
*Cell type detection using multivariate Gaussian Mixture Models.*
- `detect_cell_types(X, marker_groups, threshold=0.5, random_state=42, log_transform=True, adaptive_threshold=True) -> np.ndarray` — Binary detection per cell type using multivariate GMM.
- `build_detection_marker_groups(active_mask, type_names) -> Dict[str, List[int]]` — Build marker groups for GMM detection using curated strategy.

### CITEgeist/model/deconvolution/detection_refinement.py
*Detection refinement: GEX-based detection and iterative sparsity tuning.*
- `compute_gene_type_correlations(Y, antibody_data, antibody_names, cell_profile_dict, type_names) -> np.ndarray` — Compute correlation between each gene and each cell type's protein signal.
- `detect_cell_types_gex(Y, H, _gene_names, type_names, k=10, min_corr=0.15, threshold=0.5) -> np.ndarray` — Detect cell type presence per spot using gene expression signatures.
- `fuse_detection_masks(protein_detected, gex_detected, assignment_matrix, mode='adaptive') -> np.ndarray` — Fuse protein-based and GEX-based detection masks.
- `refine_sparsity_from_proportions(Y, sparsity_mask, cellularity=None, suppress_threshold=0.02, rescue_threshold=0.08, detection_gate_ub=0.01) -> np.ndarray` — Refine detection upper bounds using Pass 1 proportion estimates.
- `compute_marker_entropy_weights(marker_level_data, marker_names, alpha=1.0, eps=1e-10, weight_floor=0.1) -> np.ndarray` — Compute per-marker weights based on expression entropy.

### CITEgeist/model/deconvolution/emission_init.py
*Emission initialization for per-type beta QP.*
- `build_marker_config(available_markers, marker_subset=None) -> Tuple[List[str], np.ndarray, List[str]]` — Build marker panel config from available markers in the data.
- `initialize_beta_matrix(raw_counts, markers, type_names, median_N, soft_scale=0.1, inactive_val=0.001) -> np.ndarray` — Initialize per-type-per-marker beta matrix from raw counts.
- `build_beta_prior_sigma(markers, type_names, sigma_exclusive=5.0, sigma_shared=2.0, sigma_inactive=0.1, sigma_scale=1.0) -> np.ndarray` — Build per-pair regularization sigma matrix.

### CITEgeist/model/deconvolution/reference_model.py
*NB reference model for proportion refinement using scRNA-seq marker genes.*
- `train_reference(reference_adata, cell_type_col='cell_type', spatial_genes=None, n_markers_per_type=20, de_method='wilcoxon', type_mapping=None) -> ReferenceProfile` — Train NB reference profiles from annotated scRNA-seq.
- `refine_proportions_nb(spot_counts, pi_protein, reference, kappa=10.0, epsilon=1.0, lr=0.1, max_iter=200, convergence_tol=0.0001, device='cpu') -> np.ndarray` — Refine protein QP proportions using NB MAP inference on marker genes.
- `class ReferenceProfile` — Frozen NB reference profiles learned from scRNA-seq.
  - `save(path) -> None`
  - `[cls] load(path) -> 'ReferenceProfile'`

### CITEgeist/model/discovery/__init__.py

### CITEgeist/model/discovery/marker_interest.py
*Marker interest scoring for spatial transcriptomics data.*
- `identify_interesting_markers(X, coords, marker_names, kurtosis_threshold=None, morans_threshold=None, gmm_snr_threshold=1.0, morans_k=8, smooth_k=6, morans_n_perm=199, seed=1234, verbose=True) -> MarkerInterestResult` — Identify markers with interesting signal magnitude variation.
- `class MarkerInterest` — Interest metrics for a single marker.
- `class MarkerInterestResult` — Results container for marker interest analysis.
  - `to_dataframe() -> pd.DataFrame` — Convert results to ranked DataFrame sorted by interest_score descending.
  - `interesting_markers() -> List[str]` — Return names of markers passing EITHER kurtosis OR Moran's I gate (plus GMM).
  - `boring_markers() -> List[str]` — Return names of markers failing both kurtosis AND Moran's I gates (or GMM).

### CITEgeist/model/discovery/spatial_colocalization.py
*Spatial colocalization analysis for marker pairs and profile discovery.*
- `analyze_marker_colocalization(X, coords, marker_names, markers_to_analyze=None, neighbor_k=6, smooth_k=6, signal_threshold_percentile=75.0, n_permutations=999, seed=1234, verbose=True, multi_scale_k=[6, 12, 24, 48, 64], multi_scale_aggregation='max') -> ColocalizationResult` — Analyze pairwise spatial colocalization between markers.
- `discover_profiles_continuous(colocalization_result, top_k=5, distance_metric='colocalization_score', seed=1234, verbose=True) -> ProfileDiscoveryResult` — Discover cell type profiles using continuous edge weighting (no FDR gate).
- `rescue_singletons(profiles, signal_masks, signal_mask_marker_names, min_unique_coverage=0.3, min_signal_fraction=0.05, verbose=False) -> List[List[str]]` — Filter singletons by unique spatial coverage.
- `select_profiles(X, coords, marker_names, profiles, _interesting_markers, _colocalization_result, min_spatial_explained=0.9, min_protein_explained=0.9, max_profiles=15, _min_profiles=5, neighbor_k=8, verbose=False) -> ProfileSelectionResult` — Module 2c: Select profiles by dual variance contribution.
- `discover_hierarchical_profiles_continuous(coloc_result, antibody_expression, marker_names, _coords, _improvement_threshold=0.05, _sharing_ratio=0.5, _sharing_min_I=0.2, _max_depth=5, _neighbor_k=6, top_k=3, distance_metric='colocalization_score', verbose=True) -> HierarchicalProfileResult` — Discover hierarchical profiles using flat-first, hierarchy-second approach.
- `class MarkerPairColocalization` — Colocalization metrics for a pair of markers.
- `class ColocalizationResult` — Results container for colocalization analysis.
  - `to_dataframe() -> pd.DataFrame` — Convert to ranked DataFrame sorted by colocalization_score descending.
  - `get_pairs_for_marker(marker) -> List[MarkerPairColocalization]` — Get all pairs involving a specific marker.
  - `top_pairs(n=20) -> List[MarkerPairColocalization]` — Return top N pairs by colocalization score.
- `class LineageDendrogram` — Hierarchical clustering result for a single lineage (connected component).
  - `get_newick() -> str` — Convert to Newick format for visualization.
- `class ProfileDiscoveryResult` — Results from profile discovery.
  - `to_dataframe() -> pd.DataFrame` — Convert profiles to DataFrame.
  - `get_profile_for_marker(marker) -> Optional[List[str]]` — Get the profile containing a specific marker.
  - `summary() -> str` — Return a summary string.
- `class ProfileTreeNode` — A node in the hierarchical profile tree.
  - `is_leaf() -> bool` — True if this is a leaf node (no children).
  - `get_all_markers() -> List[str]` — Get all markers in this subtree (self + descendants).
- `class ProfileTree` — Hierarchical tree of marker profiles.
  - `get_leaves() -> List[ProfileTreeNode]` — Return all leaf nodes (final cell type profiles).
  - `get_depth() -> int` — Return maximum depth of tree.
- `class HierarchicalProfileResult` — Results from hierarchical profile discovery.
  - `to_profile_dict() -> Dict[str, List[str]]` — Convert to Module 3 compatible profile dictionary.
  - `summary() -> str` — Return a human-readable summary string.
- `class ProfileSelectionResult` — Results from Module 2c profile selection.
  - `summary() -> str` — Return a summary string.

### CITEgeist/model/gex/__init__.py

### CITEgeist/model/gex/cell_level_gex.py
*Distribute deconvolved GEX to individual cells.*
- `distribute_gex_to_cells(deconvolved_gex, assignments, nucleus_spot_map) -> pd.DataFrame` — Distribute spot-level deconvolved GEX to individual cells.
- `allocate_gex_type_reference(hard_labels, scores, type_names, barcodes, nucleus_ids, proportions, spot_gex) -> pd.DataFrame` — Allocate spot GEX to cells using type-specific reference profiles.

### CITEgeist/model/gex/gex_modules.py
*GEX module-aware enrichment and KL regularization for gene expression deconvolution.*
- `discover_anchor_genes(gene_expression, cell_proportions, min_anchors=5, max_anchors=10, _initial_min_correlation=0.3, min_expressing_spots=0.1) -> Tuple[Dict[int, List[int]], Dict[int, float], Dict[int, Dict[int, float]]]` — Discover anchor genes for each cell type based on correlation with proportions.
- `compute_module_aware_enrichment(_spot_idx, neighborhood_expression, base_enrichment, anchor_genes, module_weight=0.5, min_neighbors_for_corr=10) -> np.ndarray` — Compute module-aware enrichment by correlating genes with anchors in neighborhood.
- `compute_softmax_target(enrichment, temperature=0.3) -> np.ndarray` — Compute softmax target distribution from enrichment scores.
- `compute_kl_penalty_coefficients(target_distribution, total_counts, lambda_kl=0.1) -> Dict[str, np.ndarray]` — Compute coefficients for KL-divergence penalty in Gurobi objective.

### CITEgeist/model/gex/sace_gex.py
*Spatially-Adaptive Compositional EM (SACE) for per-cell GEX deconvolution.*
- `build_kernel_matrix(coords, bandwidth=None, truncate_at=3.0) -> sp.csr_matrix` — Build row-normalized Gaussian spatial kernel matrix.
- `run_sace(spot_counts, proportions, cell_assignments, cell_spot_map, spot_coords, gene_names, spotwise_profiles_init=None, n_0=10.0, bandwidth=None, max_iter=10, tol=0.0001, eps=1e-06, min_mass_threshold=1.0, damping_eta=1.0, antibody_data=None, antibody_names=None, cell_profile_dict=None) -> Tuple[Dict[int, np.ndarray], sc.AnnData, Dict]` — Run Spatially-Adaptive Compositional EM for per-cell GEX.

### CITEgeist/model/module2_proposal_builder.py
- `build_candidate_rank_lists(profiles, role_config, ontology_name) -> tuple[pd.DataFrame, pd.DataFrame]`
- `class MarkerRoleConfig`

### CITEgeist/model/morphology/__init__.py

### CITEgeist/model/morphology/morphology_backbone.py
*Morphology backbone encoders for single-cell assignment.*
- `class MorphologyBackbone` — Abstract base class for morphology feature extraction.
  - `extract(patches) -> torch.Tensor` — Extract embeddings from image patches.
  - `extract_numpy(patches, batch_size=64, device='cpu') -> np.ndarray` — Extract embeddings from numpy patches with batching.
  - `embed_dim() -> int`
- `class DAPIBackbone` — DAPI + boundary (+ optional 3rd) channel backbone using SimCLR ViT-Small encoder.
  - `extract(patches) -> torch.Tensor`
- `class HEBackbone` — H&E backbone using ImageNet-pretrained ViT-Small.
  - `extract(patches) -> torch.Tensor`

### CITEgeist/model/morphology/morphology_features.py
*Nuclear morphology feature extraction from segmentation masks and image patches.*
- `extract_nucleus_features(mask) -> pd.DataFrame` — Extract morphology features from a labeled nucleus mask.
- `extract_cell_features(nucleus_mask, cell_mask) -> pd.DataFrame` — Extract morphology features from paired nucleus and cell masks.
- `largest_remainder_discretize(proportions, n_total) -> np.ndarray` — Convert proportions to integer counts using largest remainder method.

### CITEgeist/model/morphology/morphology_prior.py
*Morphology-informed proportion prior.*
- `preprocess_morphology_prior(p_morph, detection_mask=None, eps=1e-06) -> np.ndarray` — Clip, mask, and renormalize morphology prior for solver use.
- `compute_morphology_prior(patches, cell_to_spot, oracle_props, num_types, num_spots, n_epochs=100, device='cpu', detection_mask=None) -> np.ndarray` — Train prototype-contrastive LLP and compute spot-level morphology prior.

### CITEgeist/model/morphology/patch_extraction.py
*Nucleus patch extraction for LLP backbone training and morphology-informed cell assignment.*
- `compute_global_stats(image, norm_method='percentile') -> Dict[str, Any]` — Compute normalization stats per channel for entire image.
- `extract_patch(image, bbox, expansion=0.75, output_size=96, global_stats=None) -> np.ndarray` — Extract expanded patch around a nucleus with global normalization.
- `extract_patch_with_size(image, bbox, expansion=0.75, output_size=96, global_stats=None) -> Tuple[np.ndarray, np.ndarray]` — Extract patch AND return size features.

### CITEgeist/model/morphology/prototype_contrastive.py
*Prototype-Contrastive LLP: Dual-head cell typing from spot-level proportion labels.*
- `proportion_kl_loss(oracle, predicted, eps=1e-08, type_weights=None) -> torch.Tensor` — KL divergence loss between oracle and predicted spot-level proportions.
- `consistency_loss(q_class, q_proto, eps=1e-08) -> torch.Tensor` — Symmetric KL consistency loss between classifier and prototype heads.
- `variance_covariance_loss(z, gamma=1.0) -> Tuple[torch.Tensor, torch.Tensor]` — VICReg-style variance and covariance regularisation on embeddings.
- `separation_loss(prototypes, margin=0.1) -> torch.Tensor` — Margin-based cosine hinge loss that pushes prototypes apart.
- `sharpening_loss(q, eps=1e-08) -> torch.Tensor` — Mean entropy sharpening loss (minimise to sharpen assignments).
- `tta_8x(patch) -> list` — Generate 8 deterministic test-time augmentation views of a patch.
- `train_prototype_contrastive(patches, cell_to_spot, oracle_props, num_types=6, embed_dim=128, in_channels=3, encoder_depth=12, n_epochs_2a=100, n_epochs_2b=50, encoder_warmup_epochs=10, finetune_layers=2, spots_per_batch=200, lr_2a=0.0001, lr_2b=1e-05, weight_decay=0.0001, tau_start=0.1, tau_end=0.05, lambda_consist=1.0, lambda_var=0.1, lambda_cov=0.01, lambda_sep=0.1, lambda_sharp=0.5, sep_margin=0.1, val_frac=0.2, spot_coords=None, device='cuda', simclr_checkpoint=None, prototype_refresh=True, use_type_weights=False, oracle_label_smoothing=0.0, encoder_img_size=96, augmentation_continuous_rotation=True, from_embeddings=False, encoder_embed_dim=384, external_backbone=None, warmup_freeze_backbone=True, early_stop_patience=20, seed=42) -> Dict` — Train the PrototypeContrastiveModel with two-stage LLP.
- `run_inference_tta(model, patches, cell_to_spot, num_spots, device='cuda', batch_size=256) -> Dict` — Run inference with 8-fold test-time augmentation.
- `class PrototypeContrastiveModel` — Dual-head prototype-contrastive model for cell-type proportion LLP.
  - `forward(patches) -> Dict[str, torch.Tensor]` — Forward pass.
  - `get_normalized_prototypes() -> torch.Tensor` — Return L2-normalized prototype matrix.
  - `freeze_encoder() -> None` — Freeze all encoder parameters (set requires_grad=False).
  - `unfreeze_last_n_blocks(n) -> None` — Unfreeze the last n transformer blocks and the final LayerNorm.
  - `load_simclr_encoder(checkpoint_path, device=None) -> None` — Load encoder weights from a SimCLR checkpoint.
  - `init_prototypes_from_kmeans(embeddings, num_types=None) -> None` — Initialise prototypes using K-Means on a batch of embeddings.
- `class CellPatchAugmentation` — Stochastic augmentation pipeline for cell-patch tensors.

### CITEgeist/model/morphology/segmentation.py
*StarDist-based nuclei segmentation utilities for Visium-style datasets.*
- `run_nuclei_segmentation(image, modality='he', **kwargs) -> Tuple[np.ndarray, pd.DataFrame]` — Convenience function: segment nuclei with StarDist.
- `run_cellpose_nuclei_segmentation(image_rgb_uint8, _use_gpu=False, _diameter=None, _flow_threshold=0.4, _cellprob_threshold=0.0, _model=None, _model_type='nuclei') -> Tuple[np.ndarray, np.ndarray]` — Deprecated: use ``run_nuclei_segmentation(image, modality='dapi')`` instead.
- `estimate_pixel_size_um(adata, spot_diameter_um=55.0) -> Optional[float]` — Estimate microns-per-pixel from Visium scalefactors.
- `detect_spot_diameter_pixels(adata, pixel_size_um=None, spot_diameter_um=None) -> float` — Auto-detect spot diameter in pixels for nuclei-to-spot assignment.
- `get_resolution_image_and_scale(adata, resolution_mode='hires', max_fullres_side=9000) -> Tuple[np.ndarray, float]` — Return segmentation image and fullres->image coordinate scale.
- `assign_nuclei_centroids_to_spots(centroids_xy, spot_centers_xy, spot_radius_px, spot_names) -> pd.Series` — Assign each nucleus centroid to nearest spot within spot radius.
- `normalize_nuclei_counts_for_prior(nuclei_count_raw, clip_min=0.25, clip_max=2.5) -> pd.Series` — Normalize raw nuclei counts into a soft abundance target around 1.0.
- `compute_spot_nuclei_counts_from_adata(adata, resolution_mode='hires', max_fullres_side=9000, spot_diameter_um=None, pixel_size_um=None, modality='he', **stardist_kwargs) -> SegmentationResult` — End-to-end nuclei segmentation + centroid assignment for an AnnData Visium object.
- `compute_spot_nuclei_counts_patchwise(fullres_image, spot_centers_fullres, spot_names, pixel_size_um, spot_diameter_um=55.0, patch_margin_um=25.0, modality='he', **stardist_kwargs) -> SegmentationResult` — Spot-patch-based nuclei segmentation for Visium on fullres images.
- `save_segmentation_artifacts(output_folder, sample_name, result, save_masks=True) -> Dict[str, str]` — Persist segmentation artifacts and return output paths.
- `class SegmentationQC` — Quality-control metrics for a segmentation run.
  - `to_dict() -> dict`
- `class SegmentationResult` — Container for segmentation and mapping outputs.
- `class StarDistSegmenter` — Wrapper around StarDist 2D pretrained models for nuclei segmentation.
  - `segment(image, pixel_size_um=None, min_area_um2=None, max_area_um2=None, **kwargs) -> Tuple[np.ndarray, pd.DataFrame]` — Run StarDist prediction on an image with area-based post-filtering.

### CITEgeist/model/morphology/simclr.py
*SimCLR (Simple Contrastive Learning) for nucleus morphology SSL.*
- `create_simclr_model(in_channels=2, img_size=96, patch_size=16, embed_dim=384, encoder_depth=12, encoder_num_heads=6, **kwargs) -> SimCLR` — Factory function for SimCLR model.
- `class SimCLRProjector` — SimCLR projection head (MLP).
  - `forward(x) -> torch.Tensor`
- `class SimCLR` — SimCLR self-supervised learning for nucleus morphology.
  - `forward(x1, x2) -> torch.Tensor` — Compute NT-Xent loss for two views.
  - `nt_xent_loss(z1, z2) -> torch.Tensor` — Normalized Temperature-scaled Cross Entropy Loss (NT-Xent).
  - `encode(x) -> torch.Tensor` — Encode images (for inference).

### CITEgeist/model/morphology/soft_label_classifier.py
*Soft-label classifier for morphology-to-celltype prediction.*
- `class SoftLabelClassifier` — Multinomial classifier trained on soft labels (spot proportions).
  - `fit(X, y_soft) -> 'SoftLabelClassifier'` — Fit classifier on soft labels.
  - `predict_proba(X) -> np.ndarray` — Predict cell type probabilities for each sample.
  - `feature_importances() -> np.ndarray` — Get feature importances (coefficient magnitudes per class).

### CITEgeist/model/morphology/vit_encoder.py
*Vision Transformer (ViT) encoder for nucleus morphology patches.*
- `create_vit_small(img_size=96, in_chans=2, **kwargs) -> ViTEncoder` — Factory function to create a ViT-Small encoder.
- `class PatchEmbed` — Patch embedding layer using Conv2d.
  - `forward(x) -> torch.Tensor` — Args:
- `class Attention` — Multi-head self-attention.
  - `forward(x) -> torch.Tensor` — Args:
- `class MLP` — MLP block with GELU activation.
  - `forward(x) -> torch.Tensor` — Args:
- `class TransformerBlock` — Transformer block with pre-norm architecture.
  - `forward(x) -> torch.Tensor` — Args:
- `class ViTEncoder` — Vision Transformer encoder for nucleus morphology patches.
  - `forward(x, return_patch_tokens=False) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]` — Forward pass through the ViT encoder.

### CITEgeist/model/morphology/vit_extractor.py
*Vision Transformer feature extractor for histopathology images.*
- `load_uni_extractor(weights_path, device='cuda' if torch.cuda.is_available() else 'cpu') -> ViTFeatureExtractor` — Load UNI foundation model.
- `class ViTFeatureExtractor` — ViT-based feature extractor for nucleus patches.
  - `normalize(x) -> torch.Tensor` — Apply ImageNet normalization.
  - `forward(x) -> torch.Tensor` — Extract features from patches.
  - `extract_numpy(patches, batch_size=32) -> np.ndarray` — Extract features from numpy array of patches.

### CITEgeist/model/programs/__init__.py

### CITEgeist/model/programs/anchored_program_discovery.py
*Module 4: Protein-Anchored Spatial Transcriptomic Program Discovery.*
- `analyze_program_relationships(result, adata, fdr_threshold=0.05, min_bivariate_i=0.1, n_permutations=199, neighbor_k=8, include_within_anchor=True, random_state=42) -> BivariateProgramResult` — Module 4b: Analyze bivariate spatial relationships between programs.
- `stack_deconvolved_layers(adata, layer_pattern='_genes_pass1', coord_key='spatial', cell_type_names=None) -> sc.AnnData` — Stack cell-type-specific deconvolved gene expression layers into a single AnnData.
- `unstack_program_results(stacked_H, stacked_adata, n_spots) -> Dict[str, NDArray[np.floating]]` — Unstack program activities back to per-cell-type matrices.
- `extract_celltype_expression(adata, cell_type, layer_pattern='_genes_pass1') -> Tuple[NDArray[np.floating], NDArray[np.floating]]` — Extract deconvolved expression for a single cell type from Module 3 layers.
- `anchored_spatial_nmf(X, weights, coords, K=5, lambda_spatial=0.0, lambda_sparsity=0.01, max_iter=200, random_state=42) -> Tuple[NDArray[np.floating], NDArray[np.floating], float]` — Perform weighted NMF with spatial Laplacian regularization.
- `validate_programs_with_proteins(H, protein_data, protein_names, anchor_proteins) -> pd.DataFrame` — Validate discovered programs against protein expression.
- `detect_spatial_subpopulations(H, coords, n_clusters=3, spatial_weight=0.3, min_spots_per_cluster=10) -> List[SpatialSubpopulation]` — Detect spatially distinct subpopulations within an anchor cell type.
- `discover_anchored_programs(adata, cell_type_proportions, profile_discovery_result=None, protein_adata=None, K_programs=5, lambda_spatial=0.0, lambda_sparsity=0.01, min_proportion_threshold=0.1, contrastive_strength=0.7, use_enriched_genes=True, enriched_gene_fc=1.2, max_enriched_genes=2000, validate_with_proteins=True, top_n_genes=50, detect_subpopulations=True, n_subpopulations=3, random_state=42) -> AnchoredProgramDiscoveryResult` — Discover ANCHOR-SPECIFIC protein-anchored spatial transcriptomic programs.
- `discover_joint_programs(adata, cell_type_proportions, K_programs=10, layer_pattern='_genes_pass1', lambda_spatial=0.0, lambda_sparsity=0.01, top_n_genes=50, random_state=42) -> JointDiscoveryResult` — Discover spatial programs jointly across all cell types.
- `store_results_in_adata(adata, result) -> None` — Store program discovery results in AnnData object.
- `contrastive_anchored_nmf(X, anchor_weights, background, coords, K=5, lambda_spatial=0.0, lambda_sparsity=0.01, contrastive_strength=0.8, max_iter=200, random_state=42) -> Tuple[NDArray[np.floating], NDArray[np.floating], float]` — Contrastive NMF for anchor-specific program discovery.
- `discover_programs_from_layers(adata, layer_pattern='_genes_pass1', cell_type_proportions=None, profile_discovery_result=None, protein_adata=None, K_programs=5, lambda_spatial=0.0, lambda_sparsity=0.01, min_expression_threshold=0.0, validate_with_proteins=True, top_n_genes=50, detect_subpopulations=True, n_subpopulations=3, random_state=42) -> AnchoredProgramDiscoveryResult` — Discover spatial transcriptomic programs from DECONVOLVED expression layers.
- `analyze_program_regions(result, adata, region_column, min_spots_per_region=20) -> AnchoredProgramResult` — Annotate programs with region enrichment statistics.
- `compare_programs_by_region(result, adata, region_column, region_a, region_b, top_n_genes=50) -> pd.DataFrame` — Compare program activities and gene loadings between two regions.
- `extract_program_context_genes(result, program_id, target_gene, top_n=50, exclude_target=True) -> List[Tuple[str, float]]` — Extract genes co-loaded with a target gene in a specific program.
- `class SpatialSubpopulation` — A spatially distinct subpopulation within an anchor cell type.
- `class SpatialProgram` — A single spatial transcriptomic program discovered within a cell type context.
- `class JointProgram` — A spatial program discovered from joint analysis of all cell types.
- `class JointDiscoveryResult` — Results from joint program discovery across all cell types.
  - `summary() -> str` — Return summary string.
  - `to_dataframe() -> pd.DataFrame` — Convert programs to DataFrame.
- `class AnchoredProgramResult` — Results from discovering programs anchored to a specific cell type.
  - `get_program_genes(program_idx, top_n=50) -> List[Tuple[str, float]]` — Get top genes for a specific program.
  - `to_dataframe() -> pd.DataFrame` — Convert program summaries to DataFrame.
- `class AnchoredProgramDiscoveryResult` — Complete results from protein-anchored program discovery across all cell types.
  - `summary() -> str` — Return a summary string.
  - `to_dataframe() -> pd.DataFrame` — Convert all program summaries to a single DataFrame.
- `class ProgramPairRelationship` — Bivariate spatial relationship between two programs.
- `class BivariateProgramResult` — Module 4b results: bivariate relationships between programs.
  - `summary() -> str` — Return a summary string.
  - `to_dataframe() -> pd.DataFrame` — Convert all pairs to DataFrame.

### CITEgeist/model/programs/cross_sample_integration.py
*Module 5: Cross-Sample Integration for Generalizable Program Relationships.*
- `load_module4_from_csv(csv_path, sample_name=None) -> Dict[str, AnchoredProgramResult]` — Load Module 4 results from CSV file.
- `load_module4b_from_csv(csv_path) -> BivariateProgramResult` — Load Module 4b results from CSV file.
- `load_multi_sample_results(sample_dirs, module4_pattern='*_module4*programs.csv', module4b_pattern='*_module4b*relationships.csv') -> Tuple[Dict[str, AnchoredProgramDiscoveryResult], Dict[str, BivariateProgramResult]]` — Load Module 4 and 4b results from multiple sample directories.
- `align_gene_sets(results_by_sample) -> Tuple[Dict[Tuple[str, str], NDArray], List[str], pd.DataFrame]` — Align gene sets across samples by creating a union of all genes.
- `integrate_programs_harmony(aligned_W, metadata, n_components=30, n_clusters=50, theta=2.0, max_iter=20, tol=0.0001, random_state=42) -> Tuple[NDArray, NDArray, Dict[str, Any]]` — Harmony-style integration of program signatures across samples.
- `match_programs_across_samples(Z, metadata, aligned_W, all_genes, similarity_threshold=0.7, top_n_genes=20) -> List[AlignedProgram]` — Match programs across samples based on embedding similarity.
- `compare_bivariate_relationships(bivariate_results, aligned_programs, min_samples=2, _significance_threshold=0.05, colocalization_threshold=0.1) -> List[ConservedRelationship]` — Compare bivariate relationships across samples using aligned program IDs.
- `build_similarity_network(aligned_programs, conserved_relationships, min_program_conservation=0.3, min_relationship_conservation=0.3) -> nx.Graph` — Build NetworkX graph of conserved programs and relationships.
- `integrate_samples(module4_results, module4b_results=None, n_components=30, n_clusters=50, theta=2.0, similarity_threshold=0.7, max_iter=20, random_state=42, build_network=True) -> IntegrationResult` — Main function: integrate programs across multiple samples.
- `save_integration_results(result, output_dir, prefix='module5') -> Dict[str, Path]` — Save integration results to files.
- `class AlignedProgram` — A program aligned across multiple samples.
- `class ConservedRelationship` — A bivariate relationship that appears across multiple samples.
- `class IntegrationResult` — Module 5 results: cross-sample integration.
  - `summary() -> str` — Return a summary string.
  - `to_programs_dataframe() -> pd.DataFrame` — Convert aligned programs to DataFrame.
  - `to_relationships_dataframe() -> pd.DataFrame` — Convert conserved relationships to DataFrame.

### CITEgeist/model/proposal_review_loader.py
- `build_module3_profile_dict_from_review(reviewed) -> dict[str, dict[str, list[str]]]`
- `build_module3_5_table_from_review(reviewed) -> dict[str, dict[str, list[str]]]`

### CITEgeist/model/qc/__init__.py
*QC module for validating CITEgeist single-cell outputs.*
- `run_qc(*args, **kwargs)` — Orchestrate all QC checks. See report.run_qc for full signature.

### CITEgeist/model/qc/_types.py
*Shared type definitions for the QC subpackage.*
- `class QCResult` — Result container for a single QC module.

### CITEgeist/model/qc/canonical_markers.py
*Canonical RNA markers for cell type validation.*
- `get_available_markers(cell_type, gene_names) -> list[str]` — Return canonical markers for a cell type that exist in the gene list.

### CITEgeist/model/qc/gex_qc.py
*GEX error characterization against ground truth.*
- `compute_gex_correlations(pred_layers, gt_layers) -> dict` — Compute per-type GEX Pearson r and Spearman rho on log1p values.
- `analyze_per_gene(pred_layer, gt_layer) -> pd.DataFrame` — Per-gene spatial r and magnitude ratio for one cell type.
- `compare_nmf_programs(pred_layer, gt_layer, k=5) -> np.ndarray` — Compare NMF programs between predicted and GT GEX.
- `run_gex_qc(pred_gex_layers, gt_gex_layers=None, reference_profiles=None, nmf_k=5, exemplar_type=None) -> QCResult` — Run GEX error characterization.

### CITEgeist/model/qc/marker_enrichment.py
*Canonical marker enrichment, cross-patient consistency, and internal coherence.*
- `compute_marker_enrichment(adata, patient_id=None) -> pd.DataFrame` — Compute marker enrichment per cell type: log2FC, Wilcoxon p, AUC.
- `check_cross_patient_consistency(enrichment_all_patients, min_positive_fraction=9 / 12, marker_fail_fraction=0.3) -> list[str]` — Flag inconsistent markers and anomalous patients.
- `check_internal_coherence(proportions, gex_layers, dominant_threshold=0.3) -> dict` — Check that dominant-type spots have matching GEX marker expression.
- `run_marker_enrichment(adata, proportions, gex_layers=None, multi_patient_enrichments=None, patient_id=None) -> QCResult` — Run canonical marker validation, cross-patient, and internal coherence.

### CITEgeist/model/qc/proportion_qc.py
*Proportion error characterization against ground truth.*
- `compute_proportion_correlations(pred, gt) -> dict` — Compute Pearson r, Spearman rho, RMSE, MAE per cell type and overall.
- `decompose_scaling_error(pred, gt) -> dict` — Decompose prediction error into scaling and residual components.
- `compute_discrete_round_trip(adata, proportions) -> dict` — Aggregate single-cell assignments to per-spot fractions.
- `compute_error_confusion(pred, gt) -> pd.DataFrame` — Build signed error confusion matrix.
- `run_proportion_qc(adata, proportions, gt_proportions, spatial_coords=None) -> QCResult` — Run proportion error characterization against ground truth.

### CITEgeist/model/qc/report.py
*QC orchestrator and figure composition.*
- `run_qc(adata_per_cell, proportions, mode='self_consistency', gt_proportions=None, gt_gex_layers=None, pred_gex_layers=None, _reference_adata=None, output_dir='./qc_output', empty_umi_threshold=50, empty_genes_threshold=25) -> dict[str, QCResult]` — Orchestrate all QC checks.

### CITEgeist/model/qc/single_cell_qc.py
*Standard single-cell QC metrics and empty cell detection.*
- `compute_cell_metrics(adata) -> pd.DataFrame` — Compute per-cell QC metrics from raw count matrix.
- `detect_empty_cells(metrics, umi_threshold=50, genes_threshold=25) -> np.ndarray` — Flag cells below UMI or gene detection thresholds.
- `check_compartment_emptiness(cell_types, is_empty, threshold=0.5) -> list[str]` — Flag cell types where >threshold fraction of cells are empty.
- `run_single_cell_qc(adata, umi_threshold=50, genes_threshold=25, patient_id=None) -> QCResult` — Run single-cell QC: metrics, empty detection, summary, figures.

### CITEgeist/model/unified_config.py
*Shared configuration for the unified PC-MIL pipeline.*

### CITEgeist/model/utils.py
*Utility functions for CITEgeist including neighbor detection and validation.*
- `validate_cell_profile_dict(cell_profile_dict)` — Validate the structure of a cell profile dictionary.
- `save_results_to_output(results, filepath)` — Save results as a CSV file.
- `cleanup_memory()` — Force garbage collection to free memory.
- `setup_logging(output_folder, sample_name)` — Set up dynamic logging.
- `compute_optimal_radius(adata, n_rings=3.0, fallback_spacing=100.0, min_spots=10) -> float` — Compute optimal spatial radius for neighborhood-based operations.
- `find_fixed_radius_neighbors(spot_index, adata, radius=50)` — Find neighbors within a fixed radius of a given spot.
- `get_neighbors_with_fixed_radius(spot_index, adata, radius=50, include_center=True)` — Get indices of neighboring spots based on a fixed radius around the central spot.
- `plot_neighbors_with_fixed_radius(adata, radius=50, num_spots=5)` — Plot neighbors for multiple random central spots using `sc.pl.spatial`.
- `assert_neighborhood_size(adata, cell_profile_dict, radius=50, num_spots=5)`
- `benchmark_cell_proportions(true_proportions, predicted_proportions, cell_type_names)` — Calculate performance metrics for cell type proportion predictions.
- `export_anndata_layers(adata, output_dir, pass_number=None)` — Export all layers of an AnnData object to separate CSV files.
- `calculate_expression_metrics(ground_truth_dir, predictions_dir, normalize='range', pass_number=None)` — Calculate performance metrics for gene expression predictions.
- `plot_marker_processing_stages(marker_name, coords, raw_values, signal_prob, corrected_values, _smoothed_values, zscore_values, morans_i, p_value, passed, output_path, spot_size=30.0, dpi=150, snr=None, signal_fraction=None)` — Generate 5-panel diagnostic plot showing processing stages for a marker.

## Benchmarking/

`Benchmarking/discrete_assignment/benchmark_discrete.py`
  → Benchmark discrete cell assignment against ground truth.

`Benchmarking/morphology_benchmarking/03_sctype_annotation.py`
  → scType cell type annotation for morphology audit.

`Benchmarking/morphology_benchmarking/04_pseudovisium.py`
  → Build pseudo-Visium from ENACT segmentation + scType labels.

`Benchmarking/morphology_benchmarking/04a_qc_plots.py`
  → QC plots for morphology audit: cell density, type distribution, example patches.

`Benchmarking/morphology_benchmarking/05_simclr_pretrain.py`
  → SimCLR self-supervised pretraining on H&E patches.

`Benchmarking/morphology_benchmarking/06_run_experiment.py`
  → Run morphology audit experiment (conditions 1-3).

`Benchmarking/morphology_benchmarking/07_evaluate.py`
  → Evaluate morphology audit results across all 3 conditions.

`Benchmarking/morphology_benchmarking/07_evaluate_xenium.py`
  → Xenium morphology audit — per-sample LLP evaluation.

`Benchmarking/morphology_benchmarking/__init__.py`

`Benchmarking/morphology_benchmarking/src/__init__.py`

`Benchmarking/morphology_benchmarking/src/patch_extractor.py`
  → H&E patch extraction from WSI at nucleus centroids.

`Benchmarking/morphology_benchmarking/src/pseudovisium.py`
  → Pseudo-Visium hex grid construction for morphology audit.

`Benchmarking/morphology_benchmarking/src/type_collapsing.py`
  → Deterministic 9→6 cell type collapsing for morphology audit.

`Benchmarking/scripts/audit_existing_results.py`
  → Audit existing benchmark results across the Benchmarking directory.

`Benchmarking/scripts/consolidate_results.py`
  → Consolidate method results into canonical_results.json.

`Benchmarking/simulation_benchmarking/CARD/convert_reference_to_csv.py`
  → Convert simulation reference h5ad to CSV format for CARD R scripts.

`Benchmarking/simulation_benchmarking/CARD/generate_marker_genes_simulation.py`
  → Generate marker gene lists for CARD reference-free mode on simulation data.

`Benchmarking/simulation_benchmarking/CITEgeist/debug/analyze_antibody_distribution.py`
  → Analyze antibody data distribution for rare cell types.

`Benchmarking/simulation_benchmarking/CITEgeist/debug/analyze_mixed_vs_highseg.py`
  → Analysis script to compare mixed vs high_seg simulation datasets

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py`
  → Discrete cell assignment benchmark on scCube simulation data.

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_gex_detection_sim.py`
  → Robustness check: does GEX-informed union detection help on simulation data?

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_gex_enrichment_diagnostic.py`
  → GEX enrichment diagnostic: magnitude-aware profile initialization.

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_nb_regression.py`
  → Regression benchmark: run cuOPT NB pipeline (with GMM gating, S/log(N+1) cellularity)

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_oracle_ceiling.py`
  → Oracle ceiling experiments for SACE GEX gap decomposition.

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_poisson_nmf_simulation.py`
  → Poisson NMF GEX deconvolution benchmark on scCube simulation data.

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_sace_simulation.py`
  → CITEgeist-SACE GEX deconvolution benchmark on scCube simulation data.

`Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_simulation_unified.py`
  → Unified simulation benchmark: StarDist + continuous proportions + SC evaluation.

`Benchmarking/simulation_benchmarking/CITEgeist/src/debug_gex_allocation.py`
  → Debug: trace GEX deconvolution enrichment + allocation for a CE-specific gene.

`Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py`
  → Generate synthetic nuclei images from scCube simulation data.

`Benchmarking/simulation_benchmarking/CITEgeist/src/generate_simulation_images.py`
  → Generate StarDist-tuned synthetic nuclei images from scCube simulation data.

`Benchmarking/simulation_benchmarking/CITEgeist/src/run_sace_simulation_benchmark.py`
  → Unified GEX benchmark: SACE vs QP vs NB on simulated data.

`Benchmarking/simulation_benchmarking/CITEgeist/src/summarize_simulation_results.py`
  → Aggregate simulation benchmark results into summary table.

`Benchmarking/simulation_benchmarking/CITEgeist/src/test_discretize_continuous.py`
  → Test hybrid approach: run continuous model, then discretize using nuclei counts.

`Benchmarking/simulation_benchmarking/CITEgeist/src/test_laplacian_ablation.py`
  → Test whether Laplacian spatial smoothing is the key difference maker.

`Benchmarking/simulation_benchmarking/Cell2Location/high_seg/high_seg_cell2loc.py`

`Benchmarking/simulation_benchmarking/Cell2Location/mixed/mixed_cell2loc.py`

`Benchmarking/simulation_benchmarking/Tangram/high_seg/Tangram.py`

`Benchmarking/simulation_benchmarking/Tangram/mixed/Tangram.py`

`Benchmarking/simulation_benchmarking/Tangram/src/run_benchmark.py`
  → Run Tangram deconvolution on Wu BRCA simulation data.

`Benchmarking/simulation_benchmarking/scResolve/src/generate_synthetic_images.py`
  → Generate synthetic images for scResolve benchmarking using actual scCube cell data.

`Benchmarking/simulation_benchmarking/scResolve/src/run_benchmark.py`
  → scResolve benchmark wrapper for simulated spatial transcriptomics data.

`Benchmarking/simulation_benchmarking/src/CARD_bench_wrapper.py`
  → Benchmark CARD predictions against ground truth for simulation data.

`Benchmarking/simulation_benchmarking/src/RCTD_bench_wrapper.py`

`Benchmarking/simulation_benchmarking/src/benchmarking_gex.py`

`Benchmarking/simulation_benchmarking/src/benchmarking_spot_deconv.py`

`Benchmarking/simulation_benchmarking/src/cell2location_bench_wrapper.py`

`Benchmarking/simulation_benchmarking/src/citegeist_bench_wrapper.py`

`Benchmarking/simulation_benchmarking/src/citegeist_hierarchical_wrapper.py`
  → Hierarchical Profile Discovery Benchmark on Simulated Data

`Benchmarking/simulation_benchmarking/src/compare_gex_methods.py`
  → Unified GEX method comparison for simulation benchmarks.

`Benchmarking/simulation_benchmarking/src/compare_prop_methods.py`
  → Unified proportion method comparison for simulation benchmarks.

`Benchmarking/simulation_benchmarking/src/seurat_bench_wrapper.py`

`Benchmarking/simulation_benchmarking/src/tangram_bench_wrapper.py`

`Benchmarking/visiumhd_benchmarking/slurm/download_uni.py`
  → Download UNI pathology foundation model from HuggingFace.

`Benchmarking/visiumhd_benchmarking/src/__init__.py`
  → Benchmarking framework for H&E morphology-based cell type assignment

`Benchmarking/visiumhd_benchmarking/src/create_pseudo_visium.py`
  → Create pseudo-Visium spots from Visium HD single-cell data.

`Benchmarking/visiumhd_benchmarking/src/evaluate_single_cell.py`
  → Single-cell assignment evaluation.

`Benchmarking/visiumhd_benchmarking/src/extract_patches_he.py`
  → Extract nucleus patches from H&E images.

`Benchmarking/visiumhd_benchmarking/src/run_benchmark.py`
  → Main benchmark script for Visium HD H&E morphology.

`Benchmarking/visiumhd_benchmarking/src/run_segmentation_he.py`
  → StarDist segmentation for H&E images.

`Benchmarking/visiumhd_benchmarking/src/test_create_pseudo_visium.py`
  → Tests for pseudo-Visium spot creation.

`Benchmarking/visiumhd_benchmarking/src/test_evaluate_single_cell.py`
  → Tests for single-cell evaluation.

`Benchmarking/visiumhd_benchmarking/src/test_extract_patches_he.py`
  → Tests for H&E patch extraction.

`Benchmarking/visiumhd_benchmarking/src/test_run_segmentation_he.py`
  → Tests for StarDist H&E segmentation.

`Benchmarking/visiumhd_benchmarking/src/test_train_mil.py`
  → Tests for MIL training pipeline.

`Benchmarking/visiumhd_benchmarking/src/train_mil.py`
  → Training pipeline for proportion-guided MIL.

`Benchmarking/visiumhd_segfree_poc/src/assign_transcripts_to_nuclei.py`
  → Assign Xenium transcripts to StarDist-detected nuclei.

`Benchmarking/visiumhd_segfree_poc/src/benchmark_gpu_scanning.py`
  → Benchmark: GPU-batched Adam vs CPU L-BFGS-B for scanning NB proportion solving.

`Benchmarking/visiumhd_segfree_poc/src/constants.py`
  → Shared constants for Visium HD segmentation-free POC.

`Benchmarking/visiumhd_segfree_poc/src/create_square_bins.py`
  → Create 16um square bins from Xenium transcripts and protein data.

`Benchmarking/visiumhd_segfree_poc/src/ensemble_and_evaluate.py`
  → Task 3: Confidence-weighted ensemble of CITEgeist + PC-MIL predictions.

`Benchmarking/visiumhd_segfree_poc/src/evaluate_and_plot.py`
  → Evaluate segmentation-free benchmark: Oracle vs SingleR vs CITEgeist (cellularity QP).

`Benchmarking/visiumhd_segfree_poc/src/evaluate_scanning_nb.py`
  → Evaluate saved scanning NB proportions against SingleR ground truth.

`Benchmarking/visiumhd_segfree_poc/src/extract_if_patches_and_features.py`
  → Extract 3-channel IF patches per StarDist nucleus and run frozen ViT-S for 384D features.

`Benchmarking/visiumhd_segfree_poc/src/extract_simclr_features.py`
  → Extract 384D SimCLR features using the trained encoder for all 648K nuclei.

`Benchmarking/visiumhd_segfree_poc/src/extract_simclr_patches.py`
  → Extract 96x96 2-channel (DAPI + boundary) patches for all StarDist nuclei.

`Benchmarking/visiumhd_segfree_poc/src/oracle_label_transfer.py`
  → Task 4 / Panel A: Oracle label transfer from Xenium SingleR labels to StarDist nuclei.

`Benchmarking/visiumhd_segfree_poc/src/prepare_singler_input.py`
  → Task 5a: Prepare SingleR input from StarDist nucleus transcript assignments.

`Benchmarking/visiumhd_segfree_poc/src/run_chunked_cellularity.py`
  → Chunked cellularity-aware CITEgeist deconvolution for Visium HD seg-free benchmark.

`Benchmarking/visiumhd_segfree_poc/src/run_citegeist_bins.py`
  → Task 6: Run CITEgeist deconvolution on binned pseudo-Visium data.

`Benchmarking/visiumhd_segfree_poc/src/run_scanning_nb.py`
  → End-to-end scanning NB deconvolution benchmark on Xenium pseudo-Visium HD.

`Benchmarking/visiumhd_segfree_poc/src/run_stardist_xenium.py`
  → Run StarDist nucleus segmentation on Xenium DAPI image.

`Benchmarking/visiumhd_segfree_poc/src/singler_proportions.py`
  → Task 5c: Compute bin proportions from StarDist+SingleR labels.

`Benchmarking/visiumhd_segfree_poc/src/sweep_detection.py`
  → Sweep detection strategies on saved scanning NB protein counts.

`Benchmarking/visiumhd_segfree_poc/src/train_pcmil_segfree.py`
  → Train PC-MIL on Visium HD segmentation-free POC.

`Benchmarking/visiumhd_segfree_poc/src/train_simclr_segfree.py`
  → Train SimCLR on 2-channel (DAPI + boundary) nucleus patches for segfree POC.

`Benchmarking/xenium_benchmarking/CARD/src/convert_h5ad.py`
  → Convert h5ad files to CSV format for RCTD (R).

`Benchmarking/xenium_benchmarking/CARD/src/generate_marker_genes.py`
  → Generate marker gene lists for CARD reference-free mode.

`Benchmarking/xenium_benchmarking/CITEgeist/src/__init__.py`
  → CITEgeist Benchmarking on Xenium Pseudo-Visium Data.

`Benchmarking/xenium_benchmarking/CITEgeist/src/add_18s_channel.py`
  → Add 18S channel to existing 2ch patches using same bounding boxes.

`Benchmarking/xenium_benchmarking/CITEgeist/src/aggregate_module3_5_results.py`
  → Aggregate Module 3.5 benchmark results across all 5 Xenium regions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/aggregate_sace_results.py`
  → Aggregate SACE benchmark results across all 5 Xenium pseudo-Visium regions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/aggregate_unified_results.py`
  → Aggregate unified GEX benchmark results across regions/replicates.

`Benchmarking/xenium_benchmarking/CITEgeist/src/analyze_18s_embeddings.py`
  → Analyze whether 18S channel adds discriminative information to SimCLR embeddings.

`Benchmarking/xenium_benchmarking/CITEgeist/src/analyze_concordance.py`
  → Concordance Analysis for Module 4/5 Validation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/analyze_radius_sweep.py`
  → Analyze radius sweep results and generate comparison metrics.

`Benchmarking/xenium_benchmarking/CITEgeist/src/analyze_staged_results.py`
  → Analyze Staged Evaluation Results

`Benchmarking/xenium_benchmarking/CITEgeist/src/analyze_unassigned_doublets.py`
  → Analyze unassigned rate and doublet rate in gating-based cell classification.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_assign_cells.py`
  → End-to-end benchmark for assign_cells() on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_assign_cells_morphology.py`
  → End-to-end benchmark for assign_cells() WITH morphology on Xenium Region 1.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_cell_morphology.py`
  → Benchmark cell morphology features for single-cell assignment.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_cellpose_resolution.py`
  → Benchmark nuclei-count fidelity across image resolutions on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_detection_estimation.py`
  → Benchmark detection + estimation model on Xenium pseudo-Visium data.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py`
  → DEPRECATED (2026-03-17): Discrete IQP benchmark archived.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_em_qp.py`
  → Benchmark: EM-wrapped QP proportion estimation on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_exp1a_exclusivity.py`
  → Exp 1a/1b: Exclusivity-weighted QP proportion estimation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_exp2_multivariate_gmm.py`
  → Exp 2: Joint multivariate GMM on Major markers.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_exp3_protein_clustering.py`
  → Exp 3: Protein profile clustering vs GT — information ceiling test.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py`
  → Hybrid cell assignment benchmark: continuous model + discretization via nuclei counts.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_marker_ablation.py`
  → Benchmark: graded marker ablation for cuOPT QP on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_mil_assignment.py`
  → Benchmark Module 3b MIL single-cell assignment on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_module2a_cpu_refactor.py`
  → Benchmark legacy/refactor Module 2a colocalization runtimes.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_morphology_proportions_qp.py`
  → Benchmark: morphology-informed QP proportion estimation on Xenium Region 1.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_multimodal_cellpose.py`
  → DEPRECATED (2026-03-17): Multimodal benchmark archived.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_nb_joint.py`
  → Joint tri-modal NB-GNN benchmark on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_nb_reference.py`
  → Benchmark NB reference proportion refinement on Xenium RCC.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py`
  → Benchmark PC-MIL on Xenium pseudo-Visium with 5-fold cross-validation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_per_type_beta.py`
  → Benchmark: per-type beta QP vs single-beta QP on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_prototype_llp.py`
  → Benchmark: Prototype-Contrastive LLP for single-cell typing.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_qp_fixes.py`
  → Benchmark 3 targeted QP fixes from the weakness audit.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_qp_fixes_v2.py`
  → Benchmark QP fixes v2: test Major-only vs Major+Minor markers,

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_sace_xenium.py`
  → CITEgeist-SACE GEX benchmark on Xenium pseudo-Visium (own proportions).

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_simplex_assignment.py`
  → Benchmark SimplexMIL vs gated-attention MIL for 2ch/3ch SimCLR on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_single_cell_resolution.py`
  → Single-cell resolution benchmark using Module 3b nucleus assignment pipeline.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_sparsity_refinement.py`
  → Benchmark: sparsity refinement evaluation on Xenium pseudo-Visium (5 regions).

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_spatial_gnn.py`
  → Spatial GNN for cell type classification from tissue graphs.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage.py`
  → Benchmark Two-Stage VAE-Guided Assignment.

`Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage_morphology.py`
  → Two-Stage Morphology Benchmark: Hybrid proportions + morphology-based single-cell assignment.

`Benchmarking/xenium_benchmarking/CITEgeist/src/check_panel_discrimination.py`
  → Quick check: which of the 27 proteins best discriminate Fibroblasts from Epithelial?

`Benchmarking/xenium_benchmarking/CITEgeist/src/compare_gex_cellularity.py`
  → Compare GEX deconvolution: baseline vs cellularity prior proportions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/compare_protein_rna_integration.py`
  → Compare protein-anchored vs RNA-only program discovery.

`Benchmarking/xenium_benchmarking/CITEgeist/src/compare_resolutions.py`
  → Compare CITEgeist results at single-cell vs spot-level resolution.

`Benchmarking/xenium_benchmarking/CITEgeist/src/debug_nb_standalone.py`
  → Debug NB standalone discrepancy: profile dict (10 markers) vs full panel (18 markers).

`Benchmarking/xenium_benchmarking/CITEgeist/src/diag_gt_calling_validity.py`
  → Diagnostic: GT calling validity and interaction with Module 2 / Module 3.5.

`Benchmarking/xenium_benchmarking/CITEgeist/src/diagnose_stardist_overcounting.py`
  → Diagnostic: compare StarDist nuclei counts with/without scale+prob_thresh tuning.

`Benchmarking/xenium_benchmarking/CITEgeist/src/enrich_single_cell_output_with_module3_5.py`
  → Benchmark-gated enrichment helpers for Module 3.5 single-cell outputs.

`Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_direct_softmax.py`
  → Evaluate DirectSoftmaxModel for single-cell type assignment.

`Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py`
  → Module-by-Module Pipeline Evaluation for CITEgeist

`Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_singlecell.py`
  → Evaluate single-cell demonstration results.

`Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_ssl_embeddings.py`
  → Evaluate SSL embeddings quality and classification performance.

`Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_vim_epithelial.py`
  → Evaluate VIM-Epithelial experiment results and compare with baseline.

`Benchmarking/xenium_benchmarking/CITEgeist/src/extract_ssl_embeddings.py`
  → Extract embeddings from trained SSL models (MAE or DINO).

`Benchmarking/xenium_benchmarking/CITEgeist/src/extract_vit_features.py`
  → Extract ViT-S features from Xenium nucleus patches.

`Benchmarking/xenium_benchmarking/CITEgeist/src/generate_module2_candidate_reviews.py`
  → Create ranked review artifacts from Module 2 proposals across Xenium regions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/generate_singlecell_figures.py`
  → Generate publication figures for single-cell demonstration.

`Benchmarking/xenium_benchmarking/CITEgeist/src/investigate_18s_ceiling.py`
  → Investigate the discriminative ceiling of the 18S channel for cell-type classification.

`Benchmarking/xenium_benchmarking/CITEgeist/src/load_xenium_singlecell.py`
  → Load Xenium single-cell data for CITEgeist analysis.

`Benchmarking/xenium_benchmarking/CITEgeist/src/merged_morphology_ceiling.py`
  → 4-class merged morphology ceiling: Lymphocyte, Macrophage, Stromal, Epithelial.

`Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_masked_patches.py`
  → Prepare cell-isolated masked patches for VAE training.

`Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py`
  → Prepare nucleus patches for VAE training.

`Benchmarking/xenium_benchmarking/CITEgeist/src/quick_dino_vs_simclr.py`
  → Quick DINO vs SimCLR comparison on region 1, ViT-Tiny, 50 epochs.

`Benchmarking/xenium_benchmarking/CITEgeist/src/rescore_module3_5_benchmark.py`
  → Re-score Module 3.5 benchmark using updated GT calling and pair list.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`
  → CITEgeist benchmark runner for Xenium pseudo-Visium data.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark_vim_epithelial.py`
  → Experiment: Add Vimentin to Epithelial profile in achievable-7 benchmark.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_gwr_benchmark.py`
  → GWR-NNLS GEX benchmark on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_module3_5_prod_pipeline.py`
  → Idempotent production wiring for benchmark-gated Module 3.5 enrichment.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_module3_5_xenium_benchmark.py`
  → Run the Module 3.5 benchmark on Xenium pseudo-Visium regions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_module4_validation.py`
  → Module 4/5 Validation on Xenium Data.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_nb_7type_benchmark.py`
  → NB emission 7-type benchmark for Xenium pseudo-Visium evaluation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_nb_canonical_benchmark.py`
  → NB emission canonical benchmark for Xenium pseudo-Visium evaluation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_panregion_module4.py`
  → Pan-Region Module 4: NMF Program Discovery with Consensus Profiles

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_phase1_segmentation.py`
  → Phase 1: StarDist nuclei segmentation for Xenium pseudo-Visium regions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_phase2_module3.py`
  → Phase 2: Module 3 cell-type proportion estimation for Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_phase3_features.py`
  → Phase 3: ViT feature extraction from nucleus patches for Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_phase4_mil_hungarian.py`
  → Phase 4: PC-MIL training + Hungarian assignment across all 5 Xenium regions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_phase5_evaluate.py`
  → Phase 5: Evaluate Phase 4 assignments against ground truth labels.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_preprocess.py`
  → Unified preprocessing for Xenium benchmark: StarDist segmentation,

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_qc_xenium_benchmark.py`
  → Run benchmark-mode QC on CITEgeist Xenium results vs SingleR ground truth.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_qp_7type_benchmark.py`
  → QP 7-type benchmark for Xenium pseudo-Visium evaluation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_qp_only_benchmark.py`
  → Gold standard QP benchmark for Xenium pseudo-Visium evaluation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_radius_sweep.py`
  → Radius sweep benchmark for CITEgeist Module 3.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_sace_benchmark.py`
  → SACE per-cell GEX benchmark on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_module12.py`
  → Run Modules 1 and 2 on Xenium single-cell data.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_module4.py`
  → Run Module 4 on Xenium single-cell data using discovered profiles.

`Benchmarking/xenium_benchmarking/CITEgeist/src/run_unified_gex_benchmark.py`
  → Unified GEX benchmark: SACE vs QP vs NB on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/supervised_cnn_ceiling.py`
  → Supervised CNN ceiling: train end-to-end on patches with per-cell GT labels.

`Benchmarking/xenium_benchmarking/CITEgeist/src/sweep_gex_lambda.py`
  → Sweep lambda_gex_reg values on region 0 to find optimal L2 regularization for GEX.

`Benchmarking/xenium_benchmarking/CITEgeist/src/sweep_gmm_threshold.py`
  → Sweep GMM bimodality and posterior thresholds on Region 0 to fix calibration.

`Benchmarking/xenium_benchmarking/CITEgeist/src/sweep_lambda_sparse.py`
  → Sweep lambda_sparse values for cell-resolution Pass 1 on a single region.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_blend_api.py`
  → Test the blend method via CitegeistModel API on Xenium 5-region.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cellularity_ceiling.py`
  → Cellularity ceiling ladder + noise robustness sweep.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cellularity_qp_benchmark.py`
  → Benchmark cellularity-scaled QP on Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_count_gated_qp.py`
  → Two-pass QP: use integer cell counts from pass 1 as sparsity gate for pass 2.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cuopt_adaptive_blend.py`
  → Adaptive blend: QP where strong, NB where QP is weak.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cuopt_matrix_perf.py`
  → A/B performance test: original vs matrix-formulated cuOPT QP solver.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cuopt_nb_fair.py`
  → Fair comparison: cuOPT QP vs NB (cuOPT warm) with IDENTICAL evaluation.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cuopt_nb_sweep.py`
  → Hyperparameter sweep: NB tuning for cuOPT QP warm start (Region 0 only).

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cuopt_parity.py`
  → Parity test: cuOPT QP vs Gurobi QP on Xenium 5-region.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_cuopt_warm_nb.py`
  → Test cuOPT QP warm start → NB model on Xenium 5-region.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_dual_objective_qp.py`
  → Dual-objective QP: combine CLR-normalized and raw-count data for deconvolution.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_gmm_multicomponent.py`
  → GMM multi-component elbow analysis + gated QP benchmark.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_imagenet_vit_3ch.py`
  → Quick test: ImageNet ViT-Small on 3ch Xenium patches vs 2ch padded.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_l1_sparsity_sweep.py`
  → Sweep L1 fraction (alpha) and lambda_reg with S/log(N+1) cellularity correction.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_nb_gated_qp.py`
  → NB-gated QP: use GMM detection as sparsity mask for cuOPT QP.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_nb_model_gated_qp.py`
  → NB-model-gated QP: use NB proportions as sparsity mask for cuOPT QP.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_qp_cellularity.py`
  → Quick experiment: cuOPT QP with cellularity-scaled reconstruction.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_qp_cellularity_v2.py`
  → Cellularity correction v2: push beyond r=0.739.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_sparse_nb_5region.py`
  → Sparse NB (CLR+Gaussian) on all 5 Xenium regions vs cuOPT QP baseline.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_spatial_residuals.py`
  → Test whether GEX deconvolution residuals have spatial structure.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_unified_prototypes.py`
  → Benchmark both unified solver prototypes on Xenium Region 0.

`Benchmarking/xenium_benchmarking/CITEgeist/src/test_weighted_assignment_qp.py`
  → Weighted assignment matrix QP: learn per-type marker weights from data.

`Benchmarking/xenium_benchmarking/CITEgeist/src/train_dino.py`
  → DINO training script for nucleus patch representation learning.

`Benchmarking/xenium_benchmarking/CITEgeist/src/train_mae.py`
  → Masked Autoencoder (MAE) Training Script for Nucleus Patches

`Benchmarking/xenium_benchmarking/CITEgeist/src/train_simclr.py`
  → SimCLR training script for nucleus patch representation learning.

`Benchmarking/xenium_benchmarking/CITEgeist/src/train_xgboost_combined.py`
  → Combined VAE + Comprehensive Morphology + XGBoost Pipeline.

`Benchmarking/xenium_benchmarking/CITEgeist/src/tune_sparse_nb_grid.py`
  → Grid search for Sparse NB tuning on Xenium 5-region benchmark.

`Benchmarking/xenium_benchmarking/CITEgeist/src/validate_module3_5_singlecell.py`
  → Validate Module 3.5 per-cell functional calls against Xenium GT protein expression.

`Benchmarking/xenium_benchmarking/CITEgeist/src/validate_morphology_fixes.py`
  → Validate morphology sprint + Epithelial fixes on all 5 Xenium regions.

`Benchmarking/xenium_benchmarking/CITEgeist/src/xenium_spatial_displacement.py`
  → Spatial displacement analysis for Xenium benchmark (region 1).

`Benchmarking/xenium_benchmarking/CITEgeist/tests/__init__.py`
  → Benchmarking tests

`Benchmarking/xenium_benchmarking/CITEgeist/tests/test_preprocessing_pipeline.py`
  → Integration test for new preprocessing pipeline.

`Benchmarking/xenium_benchmarking/Cell2Location/src/run_benchmark.py`
  → Run Cell2Location deconvolution on Xenium pseudo-Visium data.

`Benchmarking/xenium_benchmarking/Cell2Location/src/train_reference.py`
  → Train Cell2Location reference model on GSE156632 scRNA-seq data.

`Benchmarking/xenium_benchmarking/RCTD/src/convert_h5ad.py`
  → Convert h5ad files to CSV format for RCTD (R).

`Benchmarking/xenium_benchmarking/Seurat/src/convert_h5ad.py`
  → Convert h5ad files to CSV format for RCTD (R).

`Benchmarking/xenium_benchmarking/Tangram/src/run_benchmark.py`
  → Run Tangram deconvolution on Xenium pseudo-Visium data.

`Benchmarking/xenium_benchmarking/benchmark_constants.py`
  → Shared constants for Xenium benchmarking.

`Benchmarking/xenium_benchmarking/evaluation/consolidate_results.py`
  → Consolidate all benchmark results into a single canonical JSON.

`Benchmarking/xenium_benchmarking/evaluation/src/__init__.py`
  → Benchmark Evaluation for Xenium Deconvolution Methods.

`Benchmarking/xenium_benchmarking/evaluation/src/audit_qp_pipeline.py`
  → Comprehensive audit of the cuOPT QP gold standard pipeline.

`Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py`
  → Full method comparison against ground truth.

`Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py`
  → GEX deconvolution comparison across methods against protein-gated ground truth.

`Benchmarking/xenium_benchmarking/evaluation/src/compare_marker_genes_only.py`
  → Test scResolve performance on marker genes only.

`Benchmarking/xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py`
  → Cross-method GEX evaluation on Xenium pseudo-Visium (SingleR GT, 7 types).

`Benchmarking/xenium_benchmarking/evaluation/src/compute_heterogeneity_metrics.py`
  → compute_heterogeneity_metrics.py

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_benchmark.py`
  → Unified evaluation for CITEgeist Xenium benchmark.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_cell_classification.py`
  → Three-tier evaluation of Module 3 protein gating against RNA-based annotations.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_cell_morphology.py`
  → Evaluate cell morphology features for cell type classification.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_discovery_methods.py`
  → Evaluate Module 1-2 vs Leiden clustering on discovery quality.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_gex.py`
  → Evaluate CITEgeist gene expression deconvolution results against ground truth.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_gex_granular.py`
  → Evaluate GEX deconvolution across methods using granular 10-cell-type ground truth.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_gex_spatial.py`
  → Evaluate GEX deconvolution with spatial-aware metrics.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_metrics.py`
  → Evaluate benchmark metrics for CITEgeist predictions.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`
  → Evaluate morphology-guided cell type assignment accuracy.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_single_cell_resolution.py`
  → Evaluate single-cell resolution benchmark results.

`Benchmarking/xenium_benchmarking/evaluation/src/evaluate_vae_assignment.py`
  → Evaluate VAE + Sinkhorn single-cell assignment against Xenium ground truth.

`Benchmarking/xenium_benchmarking/evaluation/src/export_reference_for_seurat.py`
  → Export GSE156632 reference h5ad to 10x-style MTX format for Seurat loading.

`Benchmarking/xenium_benchmarking/evaluation/src/generate_heterogeneity_figures.py`
  → generate_heterogeneity_figures.py

`Benchmarking/xenium_benchmarking/evaluation/src/generate_singler_gex_gt.py`
  → Generate SingleR-based GEX ground truth for Xenium pseudo-Visium.

`Benchmarking/xenium_benchmarking/evaluation/src/leiden_baseline_comparison.py`
  → Leiden clustering baseline for Module 1-2 discovery comparison.

`Benchmarking/xenium_benchmarking/evaluation/src/module12_discovery_runner.py`
  → Module 1-2 discovery runner for comparison experiment.

`Benchmarking/xenium_benchmarking/evaluation/src/run_cell_evaluation_subset.py`
  → Run three-tier cell classification evaluation on the intersection of

`Benchmarking/xenium_benchmarking/evaluation/src/run_validated_benchmark.py`
  → Run benchmark with spatially-validated cell type definitions.

`Benchmarking/xenium_benchmarking/evaluation/src/sctype_annotation.py`
  → ScType-based per-cell RNA annotation for Xenium single-cell data.

`Benchmarking/xenium_benchmarking/evaluation/src/singlecell_validation.py`
  → Single-cell validation analysis for CITEgeist autodiscovery.

`Benchmarking/xenium_benchmarking/evaluation/src/spot_level_validation.py`
  → Spot-level validation: Compare autodiscovered profiles to actual cell counts.

`Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py`
  → Tests for morphology-guided cell type assignment evaluation.

`Benchmarking/xenium_benchmarking/figures/plot_benchmark_overview.py`
  → Generate publication-quality benchmark figures:

`Benchmarking/xenium_benchmarking/figures/plot_gex_comparison.py`
  → Fair GEX deconvolution comparison: CITEgeist vs scResolve (Region 0).

`Benchmarking/xenium_benchmarking/ground_truth_singler/src/aggregate_to_spots.py`
  → Aggregate per-cell SingleR labels to pseudo-Visium spot proportions.

`Benchmarking/xenium_benchmarking/ground_truth_singler/src/extract_xenium_data.py`
  → Extract cell-level GEX and protein data from Xenium cell_feature_matrix.h5.

`Benchmarking/xenium_benchmarking/ground_truth_singler/src/split_cd4_cd8.py`
  → Split SingleR 'T cells' into CD4+ and CD8+ using Xenium protein expression.

`Benchmarking/xenium_benchmarking/ground_truth_singler/src/validate_annotations.py`
  → Validate SingleR annotations and compute concordance with existing GT approaches.

`Benchmarking/xenium_benchmarking/reference_data/GSE156632/src/process_reference.py`
  → Process GSE156632 scRNA-seq reference data for Xenium benchmarking.

`Benchmarking/xenium_benchmarking/reference_data/GSE156632/src/split_tcells_to_protein7.py`
  → Derive the 7-type Xenium reference from the authoritative 6-type annotation.

`Benchmarking/xenium_benchmarking/scResolve/figures/generate_gex_investigation_figure.py`
  → Generate figure documenting scResolve GEX performance investigation.

`Benchmarking/xenium_benchmarking/scResolve/src/annotate_cells.py`
  → Annotate scResolve segmented cells using GSE156632 reference scRNA-seq.

`Benchmarking/xenium_benchmarking/scResolve/src/annotate_cells_balanced.py`
  → Balanced label transfer for scResolve cells using reference scRNA-seq.

`Benchmarking/xenium_benchmarking/scResolve/src/compare_multimodal_normalization.py`
  → Compare scResolve multimodal results with different protein normalization strategies.

`Benchmarking/xenium_benchmarking/scResolve/src/diagnose_rna_by_celltype.py`
  → Diagnose scResolve RNA expression by cell type.

`Benchmarking/xenium_benchmarking/scResolve/src/diagnose_rna_only_consistency.py`
  → Diagnose scResolve RNA-only mode internal consistency.

`Benchmarking/xenium_benchmarking/scResolve/src/diagnose_xenium_ground_truth.py`
  → Diagnose Xenium ground truth RNA-protein concordance.

`Benchmarking/xenium_benchmarking/scResolve/src/evaluate_cell_recovery.py`
  → Evaluate scResolve's ability to recover individual cells from pseudo-Visium spots.

`Benchmarking/xenium_benchmarking/scResolve/src/evaluate_gex_all_regions.py`
  → Compare GEX recovery across all deconvolution methods - ALL REGIONS.

`Benchmarking/xenium_benchmarking/scResolve/src/evaluate_gex_comparison.py`
  → Compare GEX recovery across all deconvolution methods.

`Benchmarking/xenium_benchmarking/scResolve/src/evaluate_gex_reconstruction.py`
  → Evaluate scResolve gene expression reconstruction quality.

`Benchmarking/xenium_benchmarking/scResolve/src/evaluate_protein_recovery.py`
  → Evaluate scResolve's ability to recover per-cell PROTEIN expression.

`Benchmarking/xenium_benchmarking/scResolve/src/extract_dapi_regions.py`
  → Extract high-resolution DAPI crops for each pseudo-Visium region.

`Benchmarking/xenium_benchmarking/scResolve/src/extract_morphology_regions.py`
  → Extract morphology image regions from Xenium OME-TIFF for scResolve benchmarking.

`Benchmarking/xenium_benchmarking/scResolve/src/fair_comparison.py`
  → Fair comparison framework for scResolve vs deconvolution methods.

`Benchmarking/xenium_benchmarking/scResolve/src/gate_and_aggregate_gex.py`
  → Gate scResolve cells by protein markers and aggregate GEX per cell type.

`Benchmarking/xenium_benchmarking/scResolve/src/generate_cell_types.py`
  → Generate cell_types.csv from Xenium RNA cluster annotations.

`Benchmarking/xenium_benchmarking/scResolve/src/generate_synthetic_images.py`
  → Generate synthetic images for scResolve benchmarking.

`Benchmarking/xenium_benchmarking/scResolve/src/inspect_cell_level.py`
  → Inspect scResolve cell-level outputs to evaluate segmentation quality

`Benchmarking/xenium_benchmarking/scResolve/src/monkey_patch_scresolve.py`
  → Monkey patch for scResolve to fix compatibility with newer SpaceRanger outputs.

`Benchmarking/xenium_benchmarking/scResolve/src/resume_segmentation.py`
  → Resume incomplete scResolve segmentation and run aggregation.

`Benchmarking/xenium_benchmarking/scResolve/src/run_benchmark.py`
  → scResolve benchmark wrapper for Xenium pseudo-Visium data with real morphology.

`Benchmarking/xenium_benchmarking/scResolve/src/run_benchmark_multimodal.py`
  → scResolve benchmark wrapper with MULTIMODAL support (GEX + Protein).

`Benchmarking/xenium_benchmarking/scResolve/src/run_benchmark_multimodal_v2.py`
  → scResolve benchmark wrapper with MULTIMODAL support (GEX + Protein) - Version 2.

`Benchmarking/xenium_benchmarking/scResolve/src/run_benchmark_multimodal_v3.py`
  → scResolve benchmark wrapper with MULTIMODAL support - Version 3.

`Benchmarking/xenium_pseudovisium/analyze_cluster_profiles.py`
  → Analyze protein marker profiles for each RNA cluster to create accurate cell type annotations.

`Benchmarking/xenium_pseudovisium/src/__init__.py`
  → Xenium Pseudo-Visium Data Generation.

`Benchmarking/xenium_pseudovisium/src/analyze_celltype_concordance.py`
  → Analyze concordance between protein-gated and RNA-defined cell types.

`Benchmarking/xenium_pseudovisium/src/analyze_consensus_gt.py`
  → Analyze what ground truth looks like when requiring BOTH protein AND RNA agreement.

`Benchmarking/xenium_pseudovisium/src/analyze_unknown_rna_distribution.py`
  → Analyze what RNA cell types the protein-Unknown cells belong to.

`Benchmarking/xenium_pseudovisium/src/analyze_vim_epithelial_overlap.py`
  → Analyze vimentin/epithelial overlap and fibroblast gate breakdown.

`Benchmarking/xenium_pseudovisium/src/create_protein_gt.py`
  → Create protein-based ground truth using single-cell protein gating.

`Benchmarking/xenium_pseudovisium/src/create_pseudo_spots.py`
  → Create pseudo-Visium spots from Xenium single-cell data.

`Benchmarking/xenium_pseudovisium/src/create_rna_gt.py`
  → Create RNA-based ground truth using Xenium's gene expression k-means clustering.

`Benchmarking/xenium_pseudovisium/src/define_cell_types.py`
  → Define cell types based on protein marker expression.

`Benchmarking/xenium_pseudovisium/src/define_cell_types_validated.py`
  → Spatially-validated cell type definitions based on single-cell co-expression analysis.

`Benchmarking/xenium_pseudovisium/src/diagnose_cd4_discrepancy.py`
  → Diagnose CD4+ T cell discrepancy between protein gating and RNA clustering.

`Benchmarking/xenium_pseudovisium/src/diagnose_unknown_cells.py`
  → Diagnose protein expression in "Unknown" cells to assess gating stringency.

`Benchmarking/xenium_pseudovisium/src/generate_expression_weighted_gt.py`
  → Generate expression-weighted ground truth proportions.

`Benchmarking/xenium_pseudovisium/src/generate_gex_ground_truth.py`
  → Generate gene expression ground truth for Xenium benchmarking.

`Benchmarking/xenium_pseudovisium/src/generate_granular_dataset.py`
  → Generate granular 10-cell-type pseudo-Visium dataset.

`Benchmarking/xenium_pseudovisium/src/generate_ground_truth.py`
  → Generate ground truth cell type proportions from single-cell data.

`Benchmarking/xenium_pseudovisium/src/generate_protein_gt_gex.py`
  → Generate GEX ground truth using protein-gated cell type assignments.

`Benchmarking/xenium_pseudovisium/src/load_xenium.py`
  → Load Xenium data into AnnData format.

`Benchmarking/xenium_pseudovisium/src/rna_cell_types.py`
  → RNA-based cell type classification using Xenium's built-in clustering.

`Benchmarking/xenium_pseudovisium/src/split_regions.py`
  → Split tissue into non-overlapping regions for pseudo-replicates.

`Benchmarking/xenium_pseudovisium/src/test_lower_thresholds.py`
  → Test concordance with lower protein gating thresholds.

`Benchmarking/xenium_pseudovisium/src/validate_ground_truth_umap.py`
  → UMAP Validation of Xenium 10-Cell-Type Ground Truth

## tests/

`tests/benchmark_array.py`
  → Benchmark Array Results Analyzer

`tests/benchmarking/__init__.py`

`tests/benchmarking/test_patch_extractor.py`
  → Tests for H&E patch extraction from WSI.

`tests/benchmarking/test_pseudovisium.py`
  → Tests for pseudo-Visium hex grid construction.

`tests/benchmarking/test_type_collapsing.py`
  → Tests for deterministic 9→6 cell type collapsing.

`tests/conftest.py`
  → Pytest configuration and fixtures for CITEgeist tests.

`tests/test_anchored_programs.py`
  → Test Module 4: Protein-Anchored Spatial Transcriptomic Program Discovery.

`tests/test_build_codebase_index.py`
  → tests/test_build_codebase_index.py

`tests/test_cell_assignment.py`
  → Tests for discrete cell assignment with morphology nudge.

`tests/test_cell_level_gex.py`
  → Tests for cell-level GEX distribution.

`tests/test_checkpoints.py`
  → Unit tests for CITEgeist checkpoint management module.

`tests/test_cross_sample_integration.py`
  → Test Module 5: Cross-Sample Integration.

`tests/test_cuopt_qp.py`
  → Smoke test for cuOPT QP proportion solver.

`tests/test_detection.py`
  → Tests for cell type detection module.

`tests/test_detection_estimation.py`
  → Tests for combined detection + estimation pipeline.

`tests/test_detection_refinement.py`
  → Tests for detection refinement module (GEX detection + sparsity refinement).

`tests/test_ensemble_proportions.py`
  → CITEgeist/tests/test_ensemble_proportions.py

`tests/test_figures/__init__.py`

`tests/test_figures/test_spatial_utils.py`
  → Tests for manuscript/figures/_shared/spatial_utils.py

`tests/test_functional_annotation.py`
  → Unit tests for CITEgeist/model/functional_annotation.py (Module 3.5).

`tests/test_generate_module2_candidate_reviews.py`

`tests/test_gex_initialization.py`
  → Tests for GEX initialization: anchor masks, size factors, Poisson IRLS.

`tests/test_gex_nb_likelihood.py`
  → Tests for GEX NB log-likelihood.

`tests/test_gwr_gex.py`
  → Tests for hierarchical shrinkage GEX deconvolution (gwr_gex.py).

`tests/test_hungarian_assignment.py`
  → Tests for Hungarian assignment algorithm.

`tests/test_hungarian_weighted.py`
  → Tests for proportion-weighted Hungarian assignment.

`tests/test_integration.py`
  → Integration tests for CITEgeist full pipeline.

`tests/test_joint_nb_optimizer.py`
  → Tests for joint NB optimizer (Stage 3).

`tests/test_mae.py`
  → Tests for Masked Autoencoder (MAE) module.

`tests/test_marker_interest.py`
  → Test harness for marker_interest module.

`tests/test_marker_validation.py`
  → Tests for single-cell marker gene validation.

`tests/test_masked_iqp.py`
  → Tests for masked IQP solver.

`tests/test_model.py`
  → Unit tests for CITEgeist model preprocessing and core functionality.

`tests/test_module2_profile_discovery.py`
  → Test script for Module 2: Profile Discovery from Spatial Colocalization.

`tests/test_module2_proposal_builder.py`

`tests/test_module2b_relaxed.py`
  → Quick diagnostic to test Module 2b with relaxed parameters.

`tests/test_module2c_profile_selection.py`
  → Test script for Module 2c: Reconstruction-Based Profile Selection.

`tests/test_module3_5_benchmark.py`

`tests/test_module3_5_migration.py`

`tests/test_module3_5_prod_pipeline.py`
  → Tests for benchmark-gated Module 3.5 orchestration.

`tests/test_module3_5_projection.py`
  → Tests for Module 3.5 projection helpers.

`tests/test_module3b_nucleus_assignment.py`
  → Tests for Module 3b: Per-Nucleus Assignment.

`tests/test_module4_regions.py`
  → Test Module 4c: Region-Aware Program Analysis.

`tests/test_module_one_marker_detection.py`
  → Test script for Module 1: Marker Interest Detection on real patient data.

`tests/test_module_robustness.py`
  → Test Module 1-2c robustness across mixed and high_seg simulated data.

`tests/test_morphology_backbone.py`
  → Tests for morphology backbone abstraction.

`tests/test_morphology_features.py`
  → Tests for nuclear morphology feature extraction.

`tests/test_morphology_prior.py`
  → Tests for morphology-informed proportion estimation.

`tests/test_multimodal_refinement.py`
  → Tests for multimodal refinement (Pass 1.5 + Pass 2 EM).

`tests/test_nb_cold_start.py`
  → Tests for cold-start NB deconvolution (Gurobi-free Module 3).

`tests/test_nb_functional.py`
  → Tests for nb_functional_attribution module.

`tests/test_nb_gated.py`
  → Tests for the gated NB likelihood function.

`tests/test_nb_joint_model.py`
  → Tests for the joint tri-modal model orchestrator (nb_joint_model.py).

`tests/test_nb_joint_utils.py`
  → Tests for CITEgeist.model.nb_joint_utils.

`tests/test_nb_spatial_gnn.py`
  → Tests for CITEgeist.model.nb_spatial_gnn.SpatialGNN.

`tests/test_patch_extraction.py`
  → Tests for nucleus patch extraction with global normalization.

`tests/test_pc_mil.py`
  → Tests for Protein-Conditioned MIL model.

`tests/test_per_cell_gex.py`
  → Tests for Phase 4 per-cell GEX allocation with module constraints.

`tests/test_profile_adapter.py`
  → Tests for profile dict adapter between Module 3 and PC-MIL formats.

`tests/test_projection_heads.py`
  → Tests for projection heads and prototypes.

`tests/test_proportion_mil.py`
  → Tests for proportion-guided MIL.

`tests/test_proposal_review_loader.py`

`tests/test_prototype_contrastive.py`
  → Tests for Prototype-Contrastive LLP loss functions and model.

`tests/test_prototype_learning.py`
  → Tests for prototype learning (Stage 2).

`tests/test_qc/__init__.py`

`tests/test_qc/conftest.py`
  → QC test fixtures.

`tests/test_qc/test_canonical_markers.py`
  → Tests for canonical_markers.py.

`tests/test_qc/test_gex_qc.py`
  → Tests for gex_qc.py.

`tests/test_qc/test_integration.py`
  → Integration test: full QC pipeline end-to-end.

`tests/test_qc/test_marker_enrichment.py`
  → Tests for marker_enrichment.py.

`tests/test_qc/test_proportion_qc.py`
  → Tests for proportion_qc.py.

`tests/test_qc/test_report.py`
  → Tests for report.py orchestrator.

`tests/test_qc/test_single_cell_qc.py`
  → Tests for single_cell_qc.py.

`tests/test_sace_gex.py`
  → Tests for SACE per-cell GEX deconvolution.

`tests/test_segmentation.py`
  → Unit tests for segmentation utilities (non-Cellpose dependent pieces).

`tests/test_simulation_benchmark.py`
  → Tests for simulation benchmark utility functions.

`tests/test_single_cell_e2e.py`
  → End-to-end test for single-cell resolution on synthetic data.

`tests/test_single_cell_integration.py`
  → Integration tests for single-cell resolution pipeline.

`tests/test_single_cell_output.py`
  → Tests for single-cell AnnData output.

`tests/test_sinkhorn.py`
  → Tests for Sinkhorn optimal transport.

`tests/test_smoke_benchmark_outputs.py`
  → Smoke tests validating canonical benchmark outputs exist and meet performance thresholds.

`tests/test_smoke_patient_pipeline.py`
  → End-to-end smoke test for the CITEgeist patient pipeline.

`tests/test_soft_assignment.py`
  → Tests for graded (Soft) marker assignment in map_antibodies_to_profiles_v2.

`tests/test_soft_label_classifier.py`
  → Tests for soft-label morphology classifier.

`tests/test_spatial_colocalization.py`
  → Test harness for spatial_colocalization module.

`tests/test_ssl_utils.py`
  → Tests for SSL utilities module.

`tests/test_stage2_high_purity.py`
  → Tests for high-purity spot detection.

`tests/test_stage2_integration.py`
  → Integration tests for Stage 2 pipeline.

`tests/test_stage2_model.py`
  → Tests for Stage 2 complete model.

`tests/test_stage2_morphology.py`
  → Tests for Stage 2 morphology-based assignment.

`tests/test_stage2_projection.py`
  → Tests for Stage 2 projection head.

`tests/test_stage2_prototypes.py`
  → Tests for Stage 2 type prototypes.

`tests/test_stage2_trainer.py`
  → Tests for Stage 2 trainer.

`tests/test_subtype_splitting.py`
  → Unit tests for subtype splitting via protein gates.

`tests/test_train_vae_vicreg.py`
  → Tests for VICReg-enhanced VAE training.

`tests/test_two_stage_pipeline.py`
  → Tests for two-stage pipeline orchestration.

`tests/test_utils.py`
  → Unit tests for CITEgeist utils module.

`tests/test_vae.py`
  → Tests for VAE architecture.

`tests/test_vicreg.py`
  → Tests for VICReg loss module.

`tests/test_vit_encoder.py`
  → Tests for ViT-Small encoder for nucleus morphology patches.

`tests/test_vit_extractor.py`
  → Tests for ViT feature extraction.

`tests/test_watershed_segmentation.py`
  → Tests for watershed cell segmentation.

`tests/visualize_prop_gridsearch.py`

`tests/visualize_simulated_benchmarks.py`

## manuscript/

`manuscript/figures/_shared/__init__.py`

`manuscript/figures/_shared/spatial_utils.py`
  → Shared spatial plotting utilities for manuscript figures.

`manuscript/figures/_shared/style.py`
  → Shared style configuration for CITEgeist manuscript figures.

`manuscript/figures/figure_style.py`
  → Shared style configuration for CITEgeist manuscript figures.

`manuscript/figures/midkine/__init__.py`

`manuscript/figures/midkine/generate.py`
  → Midkine (MDK) Figure — CITEgeist Patterns v7.

`manuscript/figures/midkine/test_generate.py`

`manuscript/figures/morphology_backbone/generate.py`
  → Supplemental Figure S-PF6: Morphology-Guided Cell Assignment Backbone.

`manuscript/figures/profile_discovery/__init__.py`

`manuscript/figures/profile_discovery/generate.py`
  → Figure: Profile Discovery (Module 2)

`manuscript/figures/profile_discovery/test_generate.py`

`manuscript/figures/simulated_benchmarking/__init__.py`

`manuscript/figures/simulated_benchmarking/generate.py`
  → Simulated Benchmarking Figure — CITEgeist Patterns v7.

`manuscript/figures/simulated_benchmarking/test_generate.py`
  → Smoke test: script runs and produces expected output files.

`manuscript/figures/supplementary/generate_supp_SPF5.py`
  → Supplementary Figure S-PF5 -- Proportion accuracy drives GEX quality.

`manuscript/figures/xenium_benchmarking/__init__.py`

`manuscript/figures/xenium_benchmarking/generate.py`
  → Xenium Benchmarking Figure — CITEgeist Patterns v8.

`manuscript/figures/xenium_benchmarking/generate_supp_SBV5.py`
  → Supplementary Figure S-BV5 — CARD Reference Imbalance Failure Mode.

`manuscript/figures/xenium_benchmarking/test_generate.py`
  → Smoke test: script runs and produces expected output files.

`manuscript/generate_figures_docx.py`
  → Generate CITEgeist_Patterns_v8_Figures.docx with embedded composite figures + legends.

## examples/

`examples/analyze_treatment_response.py`
  → Analyze Treatment Response: Biopsy (S1) vs Surgical (S2) Comparison.

`examples/compare_profiles.py`
  → Compare Discovered vs Curated Profiles.

`examples/compute_sample.py`

`examples/evaluate_sc_markers.py`
  → Evaluate top marker genes per cell type from single-cell GEX AnnData.

`examples/run_cuopt_qp_patient.py`
  → Run cuOPT QP proportions for a single patient sample.

`examples/run_module12_discovery.py`
  → Module 1-2 Discovery Runner for Patient Samples.

`examples/run_module3_unified.py`
  → Module 3 Runner with Unified Profile.

`examples/run_module4_discovery.py`
  → Module 4 Program Discovery Runner for Patient Samples.

`examples/run_module4b_bivariate.py`
  → Module 4b Bivariate Program Relationship Analysis.

`examples/run_module5_integration.py`
  → Module 5 Cross-Sample Integration Runner.

`examples/run_module5_pydeseq.py`
  → Module 5 PyDESeq2 Analysis

`examples/run_module5_pydeseq_pseudobulk.py`
  → Module 5 PyDESeq2 Analysis - Pseudo-bulk Approach

`examples/run_morphology_assignment.py`
  → H&E Morphology Single-Cell Assignment Pipeline for Patient Data.

`examples/summarize_morphology.py`
  → Generate summary report across all 12 patient samples.

`examples/train_pooled_mil.py`
  → Pooled multi-sample MIL training.

`examples/visualize_sc_spatial.py`
  → Visualize single-cell spatial outputs from CITEgeist 12-patient pipeline.
