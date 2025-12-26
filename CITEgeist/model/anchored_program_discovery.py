"""
Module 4: Protein-Anchored Spatial Transcriptomic Program Discovery.

Discovers gene expression programs anchored by cell type proportions from Module 3
and validated against protein profiles from Module 2.

This fills a gap in existing methods (NSF, STAMP, SpaTM) which work on
transcriptomics alone - CITEgeist leverages spatial protein data for
interpretable, validated program discovery.

Pipeline integration:
    Module 1 (marker_interest) -> Filter spatially-variable proteins
    Module 2 (spatial_colocalization) -> Discover protein profiles
    Module 3 (gurobi_impl) -> Cell type proportions (Y_refined)
    Module 4 (THIS MODULE) -> Protein-anchored gene programs
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scipy.sparse
from numpy.typing import NDArray
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import NMF

try:
    from .gurobi_impl import build_spatial_laplacian
    from .spatial_colocalization import ProfileDiscoveryResult
except ImportError:
    from gurobi_impl import build_spatial_laplacian
    from spatial_colocalization import ProfileDiscoveryResult


logger = logging.getLogger(__name__)


# =============================================================================
# Data Structures
# =============================================================================


@dataclass
class SpatialSubpopulation:
    """A spatially distinct subpopulation within an anchor cell type."""

    subpop_id: int
    """Unique identifier for this subpopulation."""

    n_spots: int
    """Number of spots in this subpopulation."""

    spot_indices: List[int]
    """Indices of spots belonging to this subpopulation."""

    spatial_centroid: Tuple[float, float]
    """Mean (x, y) coordinates of this subpopulation."""

    dominant_program: int
    """Index of the most active program in this subpopulation."""

    program_activities: Dict[int, float]
    """Mean activity of each program in this subpopulation."""

    location_label: str
    """Descriptive label (e.g., 'tumor_core', 'stromal', 'interface')."""


@dataclass
class SpatialProgram:
    """A single spatial transcriptomic program discovered within a cell type context."""

    program_id: int
    """Unique identifier for this program within its anchor."""

    top_genes: List[str]
    """Top N genes by loading weight (most characteristic of this program)."""

    gene_loadings: Dict[str, float]
    """Gene name -> loading weight for all genes in the program."""

    variance_explained: float
    """Fraction of variance explained by this program."""

    spatial_moran_i: float
    """Moran's I of program activity across spots (spatial coherence)."""

    spatial_moran_pvalue: float
    """P-value for Moran's I (from permutation test)."""

    mean_activity: float
    """Mean program activity across spots."""

    active_spots_fraction: float
    """Fraction of spots with above-median program activity."""

    subpopulations: Optional[List[int]] = None
    """Subpopulation IDs where this program is dominant (if detected)."""


@dataclass
class AnchoredProgramResult:
    """Results from discovering programs anchored to a specific cell type."""

    anchor_name: str
    """Name of the cell type anchor (e.g., 'Macrophage', 'T-cell')."""

    anchor_proteins: List[str]
    """Protein markers defining this anchor (from Module 2 profile)."""

    programs: List[SpatialProgram]
    """List of K discovered programs for this anchor."""

    W: NDArray[np.floating]
    """Gene loadings matrix: (n_genes, K_programs)."""

    H: NDArray[np.floating]
    """Spot loadings matrix: (K_programs, n_spots)."""

    gene_names: List[str]
    """Gene names corresponding to rows of W."""

    protein_correlations: pd.DataFrame
    """Program x Protein correlation matrix for validation."""

    reconstruction_error: float
    """Frobenius norm of reconstruction error."""

    n_spots_used: int
    """Number of spots used (with anchor proportion > threshold)."""

    parameters: Dict[str, Any] = field(default_factory=dict)
    """Parameters used for discovery (K, lambda_spatial, etc.)."""

    subpopulations: List[SpatialSubpopulation] = field(default_factory=list)
    """Spatially distinct subpopulations detected within this anchor."""

    def get_program_genes(self, program_idx: int, top_n: int = 50) -> List[Tuple[str, float]]:
        """
        Get top genes for a specific program.

        Args:
            program_idx: Index of the program (0 to K-1).
            top_n: Number of top genes to return.

        Returns:
            List of (gene_name, loading) tuples sorted by loading.
        """
        if program_idx >= self.W.shape[1]:
            raise ValueError(f"Program index {program_idx} out of range (K={self.W.shape[1]})")

        loadings = self.W[:, program_idx]
        top_indices = np.argsort(loadings)[::-1][:top_n]

        return [(self.gene_names[i], loadings[i]) for i in top_indices]

    def to_dataframe(self) -> pd.DataFrame:
        """Convert program summaries to DataFrame."""
        records = []
        for p in self.programs:
            records.append({
                "anchor": self.anchor_name,
                "program_id": p.program_id,
                "n_top_genes": len(p.top_genes),
                "top_genes": ", ".join(p.top_genes[:10]),
                "variance_explained": p.variance_explained,
                "spatial_moran_i": p.spatial_moran_i,
                "spatial_moran_pvalue": p.spatial_moran_pvalue,
                "mean_activity": p.mean_activity,
                "active_spots_fraction": p.active_spots_fraction,
            })
        return pd.DataFrame(records)


@dataclass
class AnchoredProgramDiscoveryResult:
    """Complete results from protein-anchored program discovery across all cell types."""

    results_by_anchor: Dict[str, AnchoredProgramResult]
    """Results for each cell type anchor."""

    n_anchors: int
    """Number of cell type anchors analyzed."""

    total_programs: int
    """Total number of programs discovered across all anchors."""

    profile_discovery_result: Optional[ProfileDiscoveryResult]
    """Module 2 result used for protein validation (if provided)."""

    parameters: Dict[str, Any] = field(default_factory=dict)
    """Global parameters used for discovery."""

    def summary(self) -> str:
        """Return a summary string."""
        lines = [
            f"Anchored Program Discovery Results",
            f"=" * 40,
            f"Anchors analyzed: {self.n_anchors}",
            f"Total programs: {self.total_programs}",
            "",
        ]

        for anchor_name, result in self.results_by_anchor.items():
            lines.append(f"  {anchor_name}:")
            lines.append(f"    Proteins: {', '.join(result.anchor_proteins)}")
            lines.append(f"    Programs: {len(result.programs)}")
            lines.append(f"    Spots used: {result.n_spots_used}")
            lines.append("")

        return "\n".join(lines)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert all program summaries to a single DataFrame."""
        dfs = []
        for result in self.results_by_anchor.values():
            dfs.append(result.to_dataframe())
        if dfs:
            return pd.concat(dfs, ignore_index=True)
        return pd.DataFrame()


# =============================================================================
# Core Functions
# =============================================================================


def _get_celltype_weights(
    cell_type_proportions: pd.DataFrame,
    cell_type_name: str,
    min_proportion_threshold: float = 0.1,
) -> Tuple[NDArray[np.floating], NDArray[np.bool_]]:
    """
    Extract cell type weights from Module 3 deconvolution output.

    For Visium data, the anchor weights come directly from deconvolved proportions.

    Args:
        cell_type_proportions: Module 3 output (N_spots x T_cell_types).
        cell_type_name: Name of the cell type to get weights for.
        min_proportion_threshold: Minimum proportion to include a spot.

    Returns:
        Tuple of (weights array, mask of included spots).
    """
    if cell_type_name not in cell_type_proportions.columns:
        raise ValueError(
            f"Cell type '{cell_type_name}' not found in proportions. "
            f"Available: {list(cell_type_proportions.columns)}"
        )

    weights = cell_type_proportions[cell_type_name].values.astype(np.float64)

    # Create mask for spots with sufficient proportion
    mask = weights >= min_proportion_threshold

    logger.debug(
        f"Cell type '{cell_type_name}': {mask.sum()}/{len(mask)} spots "
        f"have proportion >= {min_proportion_threshold}"
    )

    return weights, mask


def _compute_spatial_moran_i(
    values: NDArray[np.floating],
    coords: NDArray[np.floating],
    k: int = 8,
    n_permutations: int = 199,
) -> Tuple[float, float]:
    """
    Compute Moran's I spatial autocorrelation statistic.

    Args:
        values: Values to test for spatial autocorrelation (n_spots,).
        coords: Spatial coordinates (n_spots, 2).
        k: Number of neighbors for spatial weights.
        n_permutations: Number of permutations for p-value.

    Returns:
        Tuple of (Moran's I, p-value).
    """
    from scipy.spatial import cKDTree

    n = len(values)
    if n < 10:
        return 0.0, 1.0

    # Build k-NN weights
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=min(k + 1, n))

    # Create weight matrix (binary adjacency)
    W = np.zeros((n, n))
    for i, neighbors in enumerate(indices):
        for j in neighbors[1:]:  # Skip self
            W[i, j] = 1
            W[j, i] = 1

    # Normalize weights
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    W = W / row_sums

    # Compute Moran's I
    z = values - values.mean()
    numerator = np.sum(W * np.outer(z, z))
    denominator = np.sum(z ** 2)

    if denominator == 0:
        return 0.0, 1.0

    I = (n / W.sum()) * (numerator / denominator)

    # Permutation test for p-value
    I_perm = np.zeros(n_permutations)
    for p in range(n_permutations):
        z_perm = np.random.permutation(z)
        num_perm = np.sum(W * np.outer(z_perm, z_perm))
        I_perm[p] = (n / W.sum()) * (num_perm / denominator)

    # Two-tailed p-value
    p_value = (np.sum(np.abs(I_perm) >= np.abs(I)) + 1) / (n_permutations + 1)

    return float(I), float(p_value)


def anchored_spatial_nmf(
    X: NDArray[np.floating],
    weights: NDArray[np.floating],
    coords: NDArray[np.floating],
    K: int = 5,
    lambda_spatial: float = 0.1,
    lambda_sparsity: float = 0.01,
    max_iter: int = 200,
    random_state: int = 42,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], float]:
    """
    Perform weighted NMF with spatial Laplacian regularization.

    Solves: min_{W,H} ||D(X - WH)||_F^2 + lambda_spatial * Tr(H L H^T) + lambda_sparsity * ||W||_1

    For simplicity, we use sklearn NMF on weighted data and apply spatial
    smoothing as a post-processing step. A more rigorous implementation
    would integrate the Laplacian into the optimization loop.

    Args:
        X: Gene expression matrix (n_spots, n_genes).
        weights: Cell type weights per spot from Module 3 (n_spots,).
        coords: Spatial coordinates (n_spots, 2).
        K: Number of programs to discover.
        lambda_spatial: Spatial smoothness weight (applied via Laplacian smoothing).
        lambda_sparsity: Gene sparsity weight (sklearn's alpha parameter).
        max_iter: Maximum iterations for NMF.
        random_state: Random seed for reproducibility.

    Returns:
        Tuple of (W gene loadings, H spot loadings, reconstruction error).
    """
    n_spots, n_genes = X.shape

    # Weight the expression matrix by cell type proportion
    # Spots with higher proportion contribute more to program discovery
    sqrt_weights = np.sqrt(weights).reshape(-1, 1)
    X_weighted = X * sqrt_weights

    # Handle sparse matrices
    if scipy.sparse.issparse(X_weighted):
        X_weighted = X_weighted.toarray()

    # Ensure non-negative (required for NMF)
    X_weighted = np.maximum(X_weighted, 0)

    # Run NMF
    nmf = NMF(
        n_components=K,
        init='nndsvda',
        max_iter=max_iter,
        random_state=random_state,
        alpha_W=lambda_sparsity,
        l1_ratio=0.5,  # Balance L1 and L2 for W
    )

    try:
        H = nmf.fit_transform(X_weighted)  # (n_spots, K)
        W = nmf.components_.T  # (n_genes, K) after transpose
    except ValueError as e:
        logger.warning(f"NMF failed: {e}. Using random initialization.")
        W = np.random.rand(n_genes, K)
        H = np.random.rand(n_spots, K)

    # Apply spatial smoothing to H using Laplacian
    if lambda_spatial > 0 and n_spots >= 10:
        L = build_spatial_laplacian(coords, k=8, normed=True)

        # Smooth each program's loadings
        # H_smooth = (I + lambda * L)^(-1) @ H
        I = scipy.sparse.eye(n_spots)
        smoothing_matrix = I + lambda_spatial * L

        try:
            from scipy.sparse.linalg import spsolve
            H_smooth = np.zeros_like(H)
            for k in range(K):
                H_smooth[:, k] = spsolve(smoothing_matrix.tocsc(), H[:, k])
            H = np.maximum(H_smooth, 0)  # Keep non-negative
        except Exception as e:
            logger.debug(f"Spatial smoothing failed: {e}. Using unsmoothed H.")

    # Compute reconstruction error
    X_reconstructed = H @ W.T
    reconstruction_error = np.linalg.norm(X_weighted - X_reconstructed * sqrt_weights, 'fro')

    return W, H.T, reconstruction_error  # H transposed to (K, n_spots)


def validate_programs_with_proteins(
    H: NDArray[np.floating],
    protein_data: NDArray[np.floating],
    protein_names: List[str],
    anchor_proteins: List[str],
) -> pd.DataFrame:
    """
    Validate discovered programs against protein expression.

    Computes correlations between program activities (H) and protein levels
    to ensure programs align with expected biological patterns.

    Args:
        H: Program activities (K_programs, n_spots).
        protein_data: Protein expression matrix (n_spots, n_proteins).
        protein_names: Names of proteins.
        anchor_proteins: Proteins defining the anchor (should correlate highly).

    Returns:
        DataFrame with program x protein correlations and validation flags.
    """
    K = H.shape[0]
    n_proteins = len(protein_names)

    correlations = np.zeros((K, n_proteins))
    pvalues = np.zeros((K, n_proteins))

    for k in range(K):
        h_k = H[k, :]
        for p, prot in enumerate(protein_names):
            prot_values = protein_data[:, p]

            # Handle constant values
            if np.std(h_k) < 1e-10 or np.std(prot_values) < 1e-10:
                correlations[k, p] = 0.0
                pvalues[k, p] = 1.0
            else:
                r, pval = pearsonr(h_k, prot_values)
                correlations[k, p] = r
                pvalues[k, p] = pval

    # Create DataFrame
    df = pd.DataFrame(correlations, columns=protein_names)
    df.index.name = 'program'

    # Add validation column: does program correlate with anchor proteins?
    anchor_mask = [p in anchor_proteins for p in protein_names]
    if any(anchor_mask):
        df['anchor_correlation'] = df.iloc[:, anchor_mask].max(axis=1)
        df['validated'] = df['anchor_correlation'] > 0.3
    else:
        df['anchor_correlation'] = np.nan
        df['validated'] = False

    return df


def detect_spatial_subpopulations(
    H: NDArray[np.floating],
    coords: NDArray[np.floating],
    n_clusters: int = 3,
    spatial_weight: float = 0.3,
    min_spots_per_cluster: int = 10,
) -> List[SpatialSubpopulation]:
    """
    Detect spatially distinct subpopulations within an anchor cell type.

    Clusters spots based on both program activity patterns AND spatial location,
    enabling discovery of e.g., tumor-infiltrating vs stromal immune cells.

    Args:
        H: Program loadings matrix (K_programs x n_spots)
        coords: Spatial coordinates (n_spots x 2)
        n_clusters: Number of subpopulations to detect
        spatial_weight: Weight for spatial coordinates vs program loadings (0-1)
        min_spots_per_cluster: Minimum spots required for valid subpopulation

    Returns:
        List of SpatialSubpopulation objects, one per detected cluster
    """
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    n_spots = H.shape[1]
    K = H.shape[0]

    if n_spots < n_clusters * min_spots_per_cluster:
        logger.warning(f"Too few spots ({n_spots}) for {n_clusters} subpopulations")
        return []

    # Normalize program loadings and coordinates separately
    H_norm = StandardScaler().fit_transform(H.T)  # (n_spots, K)
    coords_norm = StandardScaler().fit_transform(coords)  # (n_spots, 2)

    # Combine features with spatial weighting
    # Higher spatial_weight = more emphasis on spatial location
    features = np.hstack([
        (1 - spatial_weight) * H_norm,
        spatial_weight * coords_norm
    ])

    # Cluster spots
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(features)

    # Build subpopulation objects
    subpopulations = []
    for cluster_id in range(n_clusters):
        mask = labels == cluster_id
        spot_indices = np.where(mask)[0].tolist()

        if len(spot_indices) < min_spots_per_cluster:
            continue

        # Compute cluster statistics
        cluster_H = H[:, mask]  # (K, n_spots_in_cluster)
        mean_activities = {k: float(cluster_H[k, :].mean()) for k in range(K)}
        dominant_program = int(np.argmax([mean_activities[k] for k in range(K)]))

        cluster_coords = coords[mask, :]
        centroid = (float(cluster_coords[:, 0].mean()), float(cluster_coords[:, 1].mean()))

        # Generate location label based on centroid position relative to all spots
        # (This is a heuristic - could be refined with actual tissue annotation)
        all_centroid = coords.mean(axis=0)
        dist_from_center = np.sqrt(
            (centroid[0] - all_centroid[0])**2 + (centroid[1] - all_centroid[1])**2
        )
        max_dist = np.sqrt(
            ((coords[:, 0] - all_centroid[0])**2 + (coords[:, 1] - all_centroid[1])**2).max()
        )

        if dist_from_center < max_dist * 0.33:
            location_label = "core"
        elif dist_from_center < max_dist * 0.66:
            location_label = "intermediate"
        else:
            location_label = "peripheral"

        subpop = SpatialSubpopulation(
            subpop_id=cluster_id,
            n_spots=len(spot_indices),
            spot_indices=spot_indices,
            spatial_centroid=centroid,
            dominant_program=dominant_program,
            program_activities=mean_activities,
            location_label=location_label,
        )
        subpopulations.append(subpop)

    # Sort by size
    subpopulations.sort(key=lambda x: x.n_spots, reverse=True)

    return subpopulations


def discover_anchored_programs(
    adata,
    cell_type_proportions: pd.DataFrame,
    profile_discovery_result: Optional[ProfileDiscoveryResult] = None,
    protein_adata=None,
    K_programs: int = 5,
    lambda_spatial: float = 0.1,
    lambda_sparsity: float = 0.01,
    min_proportion_threshold: float = 0.1,
    contrastive_strength: float = 0.7,
    use_enriched_genes: bool = True,
    enriched_gene_fc: float = 1.2,
    max_enriched_genes: int = 2000,
    validate_with_proteins: bool = True,
    top_n_genes: int = 50,
    detect_subpopulations: bool = True,
    n_subpopulations: int = 3,
    random_state: int = 42,
) -> AnchoredProgramDiscoveryResult:
    """
    Discover ANCHOR-SPECIFIC protein-anchored spatial transcriptomic programs.

    Main entry point for Module 4. Uses CONTRASTIVE NMF to find programs that
    are specific to each anchor cell type, rather than shared global patterns.

    Key features:
    1. Subtracts background expression from other cell types (contrastive)
    2. Focuses on genes enriched in anchor spots (optional)
    3. Produces programs unique to each anchor

    Args:
        adata: AnnData with gene expression (X) and spatial coords (obsm['spatial']).
        cell_type_proportions: Module 3 output - cell type proportions per spot (N x T).
        profile_discovery_result: Module 2 output - protein profiles for validation.
        protein_adata: AnnData with protein data for validation (optional).
        K_programs: Number of programs to discover per cell type.
        lambda_spatial: Spatial smoothness regularization weight.
        lambda_sparsity: Gene loading sparsity weight.
        min_proportion_threshold: Minimum cell type proportion to include a spot.
        contrastive_strength: How much background to subtract (0=none, 1=full).
        use_enriched_genes: Pre-filter to genes enriched in anchor spots.
        enriched_gene_fc: Minimum fold-change to consider a gene enriched.
        max_enriched_genes: Maximum number of enriched genes to use per anchor.
        validate_with_proteins: Whether to validate programs against proteins.
        top_n_genes: Number of top genes to report per program.
        random_state: Random seed for reproducibility.

    Returns:
        AnchoredProgramDiscoveryResult with anchor-specific programs.
    """
    import scanpy as sc

    logger.info("Starting CONTRASTIVE anchor-specific program discovery (Module 4)")
    logger.info(f"  contrastive_strength={contrastive_strength}, use_enriched_genes={use_enriched_genes}")

    # Get gene expression data
    if scipy.sparse.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)

    gene_names = list(adata.var_names)
    coords = adata.obsm['spatial']
    n_spots, n_genes = X.shape

    logger.info(f"Input: {n_spots} spots, {n_genes} genes, {len(cell_type_proportions.columns)} cell types")

    # Get protein data if available
    protein_data = None
    protein_names = None
    if validate_with_proteins and protein_adata is not None:
        if scipy.sparse.issparse(protein_adata.X):
            protein_data = protein_adata.X.toarray()
        else:
            protein_data = np.array(protein_adata.X)
        protein_names = list(protein_adata.var_names)
        logger.info(f"Protein data: {protein_data.shape[1]} proteins for validation")

    # Build profile name -> proteins mapping from Module 2
    profile_to_proteins = {}
    if profile_discovery_result is not None:
        for i, profile in enumerate(profile_discovery_result.profiles):
            # Use profile markers as the anchor proteins
            profile_name = f"Profile_{i}"
            profile_to_proteins[profile_name] = list(profile)

    # Discover programs for each cell type using CONTRASTIVE approach
    results_by_anchor = {}
    total_programs = 0

    for cell_type in cell_type_proportions.columns:
        if cell_type == "Unknown":
            continue

        logger.info(f"Processing anchor: {cell_type}")

        # Get weights from Module 3 proportions
        weights, mask = _get_celltype_weights(
            cell_type_proportions, cell_type, min_proportion_threshold
        )

        n_included = mask.sum()
        if n_included < 20:
            logger.warning(
                f"Skipping {cell_type}: only {n_included} spots have "
                f"proportion >= {min_proportion_threshold}"
            )
            continue

        # Compute background expression from anchor-LOW spots (contrastive)
        background = _compute_background_expression(
            X, cell_type_proportions, cell_type,
            low_threshold=min_proportion_threshold / 2  # Use stricter threshold for background
        )

        # Optionally select anchor-enriched genes
        if use_enriched_genes:
            enriched_idx, enriched_names = _select_anchor_enriched_genes(
                X, weights, gene_names,
                high_threshold=0.3 if n_included > 100 else min_proportion_threshold * 2,
                low_threshold=min_proportion_threshold / 2,
                min_fold_change=enriched_gene_fc,
                max_genes=max_enriched_genes,
            )
            X_subset = X[mask, :][:, enriched_idx]
            background_subset = background[mask, :][:, enriched_idx]
            gene_names_used = enriched_names
            print(f"    Enriched genes: {len(enriched_names)} (FC >= {enriched_gene_fc})")
        else:
            X_subset = X[mask, :]
            background_subset = background[mask, :]
            gene_names_used = gene_names
            enriched_idx = list(range(n_genes))
            print(f"    Using all {len(gene_names)} genes")

        weights_subset = weights[mask]
        coords_subset = coords[mask, :]

        # Determine anchor proteins
        anchor_proteins = profile_to_proteins.get(cell_type, [])
        if not anchor_proteins:
            for pname, prots in profile_to_proteins.items():
                if cell_type.lower() in pname.lower() or any(
                    cell_type.lower() in p.lower() for p in prots
                ):
                    anchor_proteins = prots
                    break

        # Run CONTRASTIVE spatially-regularized NMF
        W, H, recon_error = contrastive_anchored_nmf(
            X_subset,
            weights_subset,
            background_subset,
            coords_subset,
            K=K_programs,
            lambda_spatial=lambda_spatial,
            lambda_sparsity=lambda_sparsity,
            contrastive_strength=contrastive_strength,
            random_state=random_state,
        )

        # Build full H matrix (zeros for excluded spots)
        H_full = np.zeros((K_programs, n_spots))
        H_full[:, mask] = H

        # Validate with proteins if available
        protein_correlations = pd.DataFrame()
        if validate_with_proteins and protein_data is not None and protein_names is not None:
            protein_subset = protein_data[mask, :]
            protein_correlations = validate_programs_with_proteins(
                H, protein_subset, protein_names, anchor_proteins
            )

        # Build SpatialProgram objects for each program
        programs = []
        total_var = np.var(X_subset)

        for k in range(K_programs):
            # Get top genes (from enriched subset if applicable)
            loadings = W[:, k]
            top_indices = np.argsort(loadings)[::-1][:top_n_genes]
            top_genes = [gene_names_used[i] for i in top_indices]
            gene_loadings = {gene_names_used[i]: float(loadings[i]) for i in top_indices}

            # Compute variance explained (approximate)
            program_var = np.var(H[k, :]) * np.sum(loadings ** 2)
            var_explained = program_var / total_var if total_var > 0 else 0

            # Compute spatial Moran's I for this program
            moran_i, moran_p = _compute_spatial_moran_i(
                H[k, :], coords_subset, k=8, n_permutations=99
            )

            # Program activity statistics
            mean_activity = float(np.mean(H[k, :]))
            median_activity = float(np.median(H[k, :]))
            active_fraction = float(np.mean(H[k, :] > median_activity))

            programs.append(SpatialProgram(
                program_id=k,
                top_genes=top_genes,
                gene_loadings=gene_loadings,
                variance_explained=var_explained,
                spatial_moran_i=moran_i,
                spatial_moran_pvalue=moran_p,
                mean_activity=mean_activity,
                active_spots_fraction=active_fraction,
            ))

        # Map W back to full gene space if using enriched genes
        if use_enriched_genes:
            W_full = np.zeros((n_genes, K_programs))
            for i, idx in enumerate(enriched_idx):
                W_full[idx, :] = W[i, :]
        else:
            W_full = W

        # Detect spatial subpopulations within this anchor
        subpopulations = []
        if detect_subpopulations and n_included >= n_subpopulations * 10:
            subpopulations = detect_spatial_subpopulations(
                H=H,  # Program loadings for spots with anchor
                coords=coords_subset,
                n_clusters=n_subpopulations,
                spatial_weight=0.3,
                min_spots_per_cluster=10,
            )
            if subpopulations:
                print(f"    Subpopulations: {len(subpopulations)} detected")
                for sp in subpopulations:
                    print(f"      {sp.location_label}: {sp.n_spots} spots, dominant program {sp.dominant_program}")

        # Build AnchoredProgramResult
        results_by_anchor[cell_type] = AnchoredProgramResult(
            anchor_name=cell_type,
            anchor_proteins=anchor_proteins,
            programs=programs,
            W=W_full,
            H=H_full,
            gene_names=gene_names,
            protein_correlations=protein_correlations,
            reconstruction_error=recon_error,
            n_spots_used=n_included,
            parameters={
                'K_programs': K_programs,
                'lambda_spatial': lambda_spatial,
                'contrastive_strength': contrastive_strength,
                'use_enriched_genes': use_enriched_genes,
                'n_enriched_genes': len(enriched_idx) if use_enriched_genes else n_genes,
            },
            subpopulations=subpopulations,
        )

        total_programs += K_programs
        logger.info(
            f"  {cell_type}: {K_programs} programs, {n_included} spots, "
            f"{len(enriched_idx) if use_enriched_genes else n_genes} genes"
        )

    # Build final result
    result = AnchoredProgramDiscoveryResult(
        results_by_anchor=results_by_anchor,
        n_anchors=len(results_by_anchor),
        total_programs=total_programs,
        profile_discovery_result=profile_discovery_result,
        parameters={
            'K_programs': K_programs,
            'lambda_spatial': lambda_spatial,
            'contrastive_strength': contrastive_strength,
            'use_enriched_genes': use_enriched_genes,
            'min_proportion_threshold': min_proportion_threshold,
            'random_state': random_state,
        },
    )

    logger.info(f"Completed: {result.n_anchors} anchors, {result.total_programs} programs (contrastive)")

    return result


def store_results_in_adata(
    adata,
    result: AnchoredProgramDiscoveryResult,
) -> None:
    """
    Store program discovery results in AnnData object.

    Creates:
    - adata.obsm['X_anchored_programs_{anchor}']: Per-spot program activities
    - adata.varm['anchored_program_loadings_{anchor}']: Gene loadings
    - adata.uns['anchored_programs']: Metadata and summaries

    Args:
        adata: AnnData to store results in.
        result: AnchoredProgramDiscoveryResult to store.
    """
    # Store program activities in obsm
    for anchor_name, anchor_result in result.results_by_anchor.items():
        key = f'X_anchored_programs_{anchor_name}'
        adata.obsm[key] = anchor_result.H.T  # (n_spots, K)

        # Store gene loadings in varm if genes match
        if list(adata.var_names) == anchor_result.gene_names:
            varm_key = f'anchored_program_loadings_{anchor_name}'
            adata.varm[varm_key] = anchor_result.W  # (n_genes, K)

    # Store metadata in uns
    adata.uns['anchored_programs'] = {
        'n_anchors': result.n_anchors,
        'total_programs': result.total_programs,
        'parameters': result.parameters,
        'anchors': list(result.results_by_anchor.keys()),
        'summary': result.summary(),
    }

    logger.info(f"Stored results in adata: {result.n_anchors} anchors in obsm/varm/uns")


# =============================================================================
# MODULE 4 v2: CONTRASTIVE ANCHOR-SPECIFIC PROGRAM DISCOVERY
# =============================================================================


def _compute_background_expression(
    X: NDArray[np.floating],
    all_proportions: pd.DataFrame,
    anchor_name: str,
    low_threshold: float = 0.1,
) -> NDArray[np.floating]:
    """
    Compute background expression from spots with LOW anchor proportion.

    Background = mean expression in spots where anchor cell type is rare.
    This represents tissue expression in the ABSENCE of the anchor cell type,
    giving us a proper baseline to subtract for contrastive analysis.

    Args:
        X: Gene expression matrix (n_spots, n_genes)
        all_proportions: Cell type proportions (n_spots, n_cell_types)
        anchor_name: Name of anchor cell type to EXCLUDE from background
        low_threshold: Proportion threshold for "anchor-low" spots

    Returns:
        Background expression matrix (n_spots, n_genes) - same value per gene,
        broadcast to all spots for subtraction.
    """
    n_spots, n_genes = X.shape

    # Find spots with LOW anchor proportion
    anchor_props = all_proportions[anchor_name].values
    low_mask = anchor_props < low_threshold

    n_low = low_mask.sum()
    if n_low < 10:
        # Fallback: use median across all spots
        print(f"    Background: only {n_low} low spots, using median")
        background_vector = np.median(X, axis=0)
    else:
        # Mean expression in anchor-low spots = background
        background_vector = np.mean(X[low_mask, :], axis=0)
        print(f"    Background: {n_low} low spots (prop < {low_threshold:.2f})")

    # Broadcast to all spots (same background for each spot)
    background = np.tile(background_vector, (n_spots, 1))

    return background


def _select_anchor_enriched_genes(
    X: NDArray[np.floating],
    anchor_proportions: NDArray[np.floating],
    gene_names: List[str],
    high_threshold: float = 0.3,
    low_threshold: float = 0.1,
    min_fold_change: float = 1.2,
    max_genes: int = 2000,
) -> Tuple[List[int], List[str]]:
    """
    Select genes enriched in anchor-high spots vs anchor-low spots.

    Uses a simple fold-change approach to identify genes more highly
    expressed where the anchor cell type is abundant.

    Args:
        X: Gene expression matrix (n_spots, n_genes)
        anchor_proportions: Anchor proportion per spot (n_spots,)
        gene_names: Names of genes
        high_threshold: Proportion threshold for "anchor-high" spots
        low_threshold: Proportion threshold for "anchor-low" spots
        min_fold_change: Minimum fold-change to consider enriched
        max_genes: Maximum number of enriched genes to return

    Returns:
        Tuple of (gene indices, gene names) for enriched genes
    """
    high_mask = anchor_proportions >= high_threshold
    low_mask = anchor_proportions < low_threshold

    n_high = high_mask.sum()
    n_low = low_mask.sum()

    if n_high < 10 or n_low < 10:
        logger.warning(f"Too few spots for enrichment (high={n_high}, low={n_low}), using all genes")
        return list(range(X.shape[1])), gene_names

    # Mean expression in high vs low spots
    mean_high = np.mean(X[high_mask, :], axis=0) + 1e-6
    mean_low = np.mean(X[low_mask, :], axis=0) + 1e-6

    fold_change = mean_high / mean_low

    # Select genes with sufficient fold-change
    enriched_mask = fold_change >= min_fold_change
    enriched_indices = np.where(enriched_mask)[0]

    # Sort by fold-change and take top max_genes
    if len(enriched_indices) > max_genes:
        fc_order = np.argsort(fold_change[enriched_indices])[::-1]
        enriched_indices = enriched_indices[fc_order[:max_genes]]

    enriched_names = [gene_names[i] for i in enriched_indices]

    logger.debug(f"Selected {len(enriched_indices)} enriched genes (FC >= {min_fold_change})")

    return list(enriched_indices), enriched_names


def contrastive_anchored_nmf(
    X: NDArray[np.floating],
    anchor_weights: NDArray[np.floating],
    background: NDArray[np.floating],
    coords: NDArray[np.floating],
    K: int = 5,
    lambda_spatial: float = 0.1,
    lambda_sparsity: float = 0.01,
    contrastive_strength: float = 0.8,
    max_iter: int = 200,
    random_state: int = 42,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], float]:
    """
    Contrastive NMF for anchor-specific program discovery.

    Subtracts a fraction of the background (expression from other cell types)
    before running NMF, allowing discovery of anchor-SPECIFIC patterns.

    Solves: NMF(X_anchor - contrastive_strength * X_background)

    Args:
        X: Gene expression matrix (n_spots, n_genes)
        anchor_weights: Cell type weights from Module 3 (n_spots,)
        background: Background expression from other cell types (n_spots, n_genes)
        coords: Spatial coordinates (n_spots, 2)
        K: Number of programs to discover
        lambda_spatial: Spatial smoothness weight
        lambda_sparsity: Gene sparsity weight
        contrastive_strength: How much background to subtract (0-1)
        max_iter: Maximum NMF iterations
        random_state: Random seed

    Returns:
        Tuple of (W gene loadings, H spot loadings, reconstruction error)
    """
    n_spots, n_genes = X.shape

    # Compute contrastive expression: anchor contribution minus background
    X_contrastive = X - contrastive_strength * background

    # Report contrastive effect
    var_orig = np.var(X)
    var_contrastive = np.var(X_contrastive)
    pct_removed = 100 * (1 - var_contrastive / var_orig) if var_orig > 0 else 0
    print(f"    Contrastive: removed {pct_removed:.1f}% variance (strength={contrastive_strength})")

    # Ensure non-negative (NMF requirement) - shift to positive range
    X_min = X_contrastive.min()
    if X_min < 0:
        X_contrastive = X_contrastive - X_min + 1e-6

    # Weight by anchor proportion (spots with more anchor contribute more)
    sqrt_weights = np.sqrt(anchor_weights).reshape(-1, 1)
    X_weighted = X_contrastive * sqrt_weights

    # Handle sparse matrices
    if scipy.sparse.issparse(X_weighted):
        X_weighted = X_weighted.toarray()

    X_weighted = np.maximum(X_weighted, 0)

    # Run NMF
    nmf = NMF(
        n_components=K,
        init='nndsvda',
        max_iter=max_iter,
        random_state=random_state,
        alpha_W=lambda_sparsity,
        l1_ratio=0.5,
    )

    try:
        H = nmf.fit_transform(X_weighted)  # (n_spots, K)
        W = nmf.components_.T  # (n_genes, K)
    except ValueError as e:
        logger.warning(f"Contrastive NMF failed: {e}. Falling back to standard NMF.")
        return anchored_spatial_nmf(
            X, anchor_weights, coords, K, lambda_spatial, lambda_sparsity,
            max_iter, random_state
        )

    # Apply spatial smoothing
    if lambda_spatial > 0 and n_spots >= 10:
        L = build_spatial_laplacian(coords, k=8, normed=True)
        I = scipy.sparse.eye(n_spots)
        smoothing_matrix = I + lambda_spatial * L

        try:
            from scipy.sparse.linalg import spsolve
            H_smooth = np.zeros_like(H)
            for k in range(K):
                H_smooth[:, k] = spsolve(smoothing_matrix.tocsc(), H[:, k])
            H = np.maximum(H_smooth, 0)
        except Exception as e:
            logger.debug(f"Spatial smoothing failed: {e}")

    # Reconstruction error (on original weighted data)
    X_reconstructed = H @ W.T
    reconstruction_error = np.linalg.norm(X_weighted - X_reconstructed * sqrt_weights, 'fro')

    return W, H.T, reconstruction_error
