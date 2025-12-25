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


def discover_anchored_programs(
    adata,
    cell_type_proportions: pd.DataFrame,
    profile_discovery_result: Optional[ProfileDiscoveryResult] = None,
    protein_adata=None,
    K_programs: int = 5,
    lambda_spatial: float = 0.1,
    lambda_sparsity: float = 0.01,
    min_proportion_threshold: float = 0.1,
    validate_with_proteins: bool = True,
    top_n_genes: int = 50,
    random_state: int = 42,
) -> AnchoredProgramDiscoveryResult:
    """
    Discover protein-anchored spatial transcriptomic programs.

    Main entry point for Module 4. For each cell type from Module 3 deconvolution,
    discovers K gene expression programs weighted by cell type proportion.
    Optionally validates against Module 2 protein profiles.

    Args:
        adata: AnnData with gene expression (X) and spatial coords (obsm['spatial']).
        cell_type_proportions: Module 3 output - cell type proportions per spot (N x T).
        profile_discovery_result: Module 2 output - protein profiles for validation.
        protein_adata: AnnData with protein data for validation (optional).
        K_programs: Number of programs to discover per cell type.
        lambda_spatial: Spatial smoothness regularization weight.
        lambda_sparsity: Gene loading sparsity weight.
        min_proportion_threshold: Minimum cell type proportion to include a spot.
        validate_with_proteins: Whether to validate programs against proteins.
        top_n_genes: Number of top genes to report per program.
        random_state: Random seed for reproducibility.

    Returns:
        AnchoredProgramDiscoveryResult with programs for each cell type anchor.
    """
    import scanpy as sc

    logger.info("Starting protein-anchored program discovery (Module 4)")

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

    # Discover programs for each cell type
    results_by_anchor = {}
    total_programs = 0

    for cell_type in cell_type_proportions.columns:
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

        # Subset to spots with sufficient proportion
        X_subset = X[mask, :]
        weights_subset = weights[mask]
        coords_subset = coords[mask, :]

        # Determine anchor proteins
        anchor_proteins = profile_to_proteins.get(cell_type, [])
        if not anchor_proteins:
            # Try to match cell type name to profile
            for pname, prots in profile_to_proteins.items():
                if cell_type.lower() in pname.lower() or any(
                    cell_type.lower() in p.lower() for p in prots
                ):
                    anchor_proteins = prots
                    break

        # Run spatially-regularized NMF
        W, H, recon_error = anchored_spatial_nmf(
            X_subset,
            weights_subset,
            coords_subset,
            K=K_programs,
            lambda_spatial=lambda_spatial,
            lambda_sparsity=lambda_sparsity,
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
            # Get top genes
            loadings = W[:, k]
            top_indices = np.argsort(loadings)[::-1][:top_n_genes]
            top_genes = [gene_names[i] for i in top_indices]
            gene_loadings = {gene_names[i]: float(loadings[i]) for i in top_indices}

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

        # Build AnchoredProgramResult
        results_by_anchor[cell_type] = AnchoredProgramResult(
            anchor_name=cell_type,
            anchor_proteins=anchor_proteins,
            programs=programs,
            W=W,
            H=H_full,
            gene_names=gene_names,
            protein_correlations=protein_correlations,
            reconstruction_error=recon_error,
            n_spots_used=n_included,
            parameters={
                'K_programs': K_programs,
                'lambda_spatial': lambda_spatial,
                'lambda_sparsity': lambda_sparsity,
                'min_proportion_threshold': min_proportion_threshold,
            },
        )

        total_programs += K_programs
        logger.info(
            f"  {cell_type}: {K_programs} programs from {n_included} spots, "
            f"reconstruction error = {recon_error:.2f}"
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
            'lambda_sparsity': lambda_sparsity,
            'min_proportion_threshold': min_proportion_threshold,
            'validate_with_proteins': validate_with_proteins,
            'random_state': random_state,
        },
    )

    logger.info(f"Completed: {result.n_anchors} anchors, {result.total_programs} total programs")

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
