"""
RNA marker selection for CITEgeist pipeline.

This module enables CITEgeist to work with RNA-only data by selecting
spatially-informative marker genes that can substitute for protein markers.

Supports three modes:
1. CURATED_ONLY: User provides a list of known marker genes
2. HYBRID: Start with curated markers, discover additional spatial markers
3. AUTODISCOVERY: Fully automated selection from HVGs (exploratory)

The selected markers feed into Modules 1-2 for profile discovery.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc
from numpy.typing import NDArray
from scipy.stats import spearmanr

from .marker_interest import identify_interesting_markers, MarkerInterestResult

logger = logging.getLogger(__name__)


class MarkerMode(Enum):
    """Marker selection strategy."""
    CURATED_ONLY = "curated"           # User provides marker list
    HYBRID = "hybrid"                   # Curated + auto-discovered
    AUTODISCOVERY = "autodiscovery"    # Fully automated from HVGs


@dataclass
class RNAMarkerSelectionResult:
    """Results container for RNA marker selection."""

    selected_markers: List[str]
    curated_markers: List[str]          # Curated markers found in data
    discovered_markers: List[str]       # Auto-discovered markers
    excluded_redundant: List[str]       # Markers excluded for redundancy
    mode: MarkerMode

    # Spatial interest results for discovered markers
    spatial_interest: Optional[MarkerInterestResult] = None

    # Per-marker metadata
    marker_metadata: Dict[str, Dict] = field(default_factory=dict)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame with marker info."""
        records = []
        for marker in self.selected_markers:
            meta = self.marker_metadata.get(marker, {})
            records.append({
                "marker": marker,
                "source": "curated" if marker in self.curated_markers else "discovered",
                "interest_score": meta.get("interest_score", np.nan),
                "morans_i": meta.get("morans_i", np.nan),
                "kurtosis": meta.get("kurtosis", np.nan),
            })
        return pd.DataFrame(records)

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            f"RNA Marker Selection ({self.mode.value} mode)",
            f"  Total selected: {len(self.selected_markers)}",
            f"  From curated list: {len(self.curated_markers)}",
            f"  Auto-discovered: {len(self.discovered_markers)}",
            f"  Excluded (redundant): {len(self.excluded_redundant)}",
        ]
        return "\n".join(lines)


# ============================================================================
# Default curated marker gene sets
# ============================================================================

# Common cell type markers (subset of PanglaoDB/CellMarker)
# These are canonical markers that tend to have good spatial signal
DEFAULT_CURATED_MARKERS = {
    # T cells
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
    # B cells
    "CD19", "MS4A1", "CD79A", "CD79B",
    # NK cells
    "NCAM1", "NKG7", "GNLY",
    # Monocytes/Macrophages
    "CD14", "CD68", "CD163", "FCGR3A", "CSF1R",
    # Dendritic cells
    "ITGAX", "CD1C", "CLEC9A",
    # Fibroblasts
    "COL1A1", "COL1A2", "DCN", "FAP", "PDGFRA",
    # Endothelial
    "PECAM1", "VWF", "CDH5", "KDR",
    # Epithelial
    "EPCAM", "KRT8", "KRT18", "KRT19",
    # Smooth muscle
    "ACTA2", "MYH11", "TAGLN",
    # Proliferation
    "MKI67", "TOP2A", "PCNA",
}


# ============================================================================
# Biological equivalence mapping: Protein antibodies -> RNA genes
# ============================================================================
# This mapping enables semantic comparison between protein-discovered and
# RNA-discovered profiles. Protein antibody names (left) map to their
# corresponding RNA gene names (right).

PROTEIN_TO_RNA_MAPPING: Dict[str, List[str]] = {
    # T cell markers
    "CD3E": ["CD3E", "CD3D", "CD3G"],
    "CD3": ["CD3E", "CD3D", "CD3G"],
    "CD4": ["CD4"],
    "CD8A": ["CD8A", "CD8B"],
    "CD8": ["CD8A", "CD8B"],
    "CD45": ["PTPRC"],
    "CD45RA": ["PTPRC"],  # Isoform - maps to same gene
    "CD45RO": ["PTPRC"],  # Isoform - maps to same gene

    # B cell markers
    "CD20": ["MS4A1"],  # CD20 = MS4A1 gene
    "CD19": ["CD19"],
    "CD79A": ["CD79A"],
    "CD79B": ["CD79B"],

    # Myeloid markers (Macrophages/Monocytes/DCs)
    "CD68": ["CD68"],
    "CD163": ["CD163"],
    "CD14": ["CD14"],
    "CD16": ["FCGR3A", "FCGR3B"],  # CD16 = FcγRIII
    "HLA-DR": ["HLA-DRA", "HLA-DRB1", "HLA-DRB5"],
    "CD11c": ["ITGAX"],  # Integrin alpha X
    "CD11b": ["ITGAM"],  # Integrin alpha M

    # NK cell markers
    "CD56": ["NCAM1"],
    "NKG7": ["NKG7"],

    # Proliferation markers
    "Ki-67": ["MKI67"],
    "Ki67": ["MKI67"],
    "PCNA": ["PCNA"],

    # Epithelial markers
    "E-Cadherin": ["CDH1"],
    "EPCAM": ["EPCAM"],
    "PanCK": ["KRT8", "KRT18", "KRT19", "KRT7"],  # Pan-cytokeratin

    # Mesenchymal/Stromal markers
    "Vimentin": ["VIM"],
    "alphaSMA": ["ACTA2"],  # Alpha smooth muscle actin
    "alpha-SMA": ["ACTA2"],
    "SMA": ["ACTA2"],

    # Endothelial markers
    "CD31": ["PECAM1"],
    "VWF": ["VWF"],

    # Immune checkpoint/activation
    "PD-1": ["PDCD1"],
    "PD-L1": ["CD274"],
    "LAG-3": ["LAG3"],
    "VISTA": ["VSIR"],
    "GranzymeB": ["GZMB"],
    "Granzyme-B": ["GZMB"],

    # Plasma cell markers
    "CD138": ["SDC1"],  # Syndecan-1

    # Other
    "PTEN": ["PTEN"],
    "Beta-catenin": ["CTNNB1"],
}

# Reverse mapping: RNA gene -> protein antibody names
RNA_TO_PROTEIN_MAPPING: Dict[str, List[str]] = {}
for protein, rna_genes in PROTEIN_TO_RNA_MAPPING.items():
    for gene in rna_genes:
        if gene not in RNA_TO_PROTEIN_MAPPING:
            RNA_TO_PROTEIN_MAPPING[gene] = []
        RNA_TO_PROTEIN_MAPPING[gene].append(protein)


def get_equivalent_markers(
    marker: str,
    source: str = "protein",
) -> Set[str]:
    """
    Get biologically equivalent markers in the other modality.

    Args:
        marker: Marker name to look up.
        source: "protein" (look up RNA equivalents) or "rna" (look up protein equivalents).

    Returns:
        Set of equivalent marker names.
    """
    if source == "protein":
        # Look up RNA equivalents for a protein
        return set(PROTEIN_TO_RNA_MAPPING.get(marker, [marker]))
    else:
        # Look up protein equivalents for an RNA gene
        return set(RNA_TO_PROTEIN_MAPPING.get(marker, [marker]))


def get_curated_markers(
    tissue_type: Optional[str] = None,
    additional_markers: Optional[List[str]] = None,
) -> Set[str]:
    """
    Get curated marker gene set.

    Args:
        tissue_type: Optional tissue type for tissue-specific markers.
            Currently supports: None (default pan-tissue set)
            Future: "tumor", "immune", "brain", etc.
        additional_markers: Additional markers to include.

    Returns:
        Set of marker gene names.
    """
    markers = DEFAULT_CURATED_MARKERS.copy()

    # TODO: Add tissue-specific marker sets
    if tissue_type is not None:
        logger.warning(f"Tissue-specific markers for '{tissue_type}' not yet implemented, using defaults")

    if additional_markers:
        markers.update(additional_markers)

    return markers


# ============================================================================
# Core selection functions
# ============================================================================

def _compute_redundancy_matrix(
    X: NDArray[np.floating],
    gene_names: List[str],
) -> pd.DataFrame:
    """
    Compute pairwise Spearman correlation matrix for redundancy detection.

    Args:
        X: Expression matrix (n_spots, n_genes).
        gene_names: Gene names corresponding to columns.

    Returns:
        DataFrame with pairwise correlations.
    """
    n_genes = len(gene_names)
    corr_matrix = np.zeros((n_genes, n_genes))

    for i in range(n_genes):
        for j in range(i, n_genes):
            if i == j:
                corr_matrix[i, j] = 1.0
            else:
                rho, _ = spearmanr(X[:, i], X[:, j])
                corr_matrix[i, j] = rho
                corr_matrix[j, i] = rho

    return pd.DataFrame(corr_matrix, index=gene_names, columns=gene_names)


def _filter_redundant_markers(
    candidates: List[str],
    reference_markers: List[str],
    X: NDArray[np.floating],
    gene_to_idx: Dict[str, int],
    threshold: float = 0.7,
) -> Tuple[List[str], List[str]]:
    """
    Filter candidates that are redundant with reference markers.

    A candidate is redundant if it has Spearman correlation > threshold
    with any reference marker.

    Args:
        candidates: Candidate markers to filter.
        reference_markers: Reference markers to compare against.
        X: Expression matrix.
        gene_to_idx: Mapping from gene name to column index.
        threshold: Correlation threshold for redundancy.

    Returns:
        Tuple of (kept_markers, excluded_markers).
    """
    if not reference_markers or not candidates:
        return candidates, []

    kept = []
    excluded = []

    # Get reference expression vectors
    ref_indices = [gene_to_idx[g] for g in reference_markers if g in gene_to_idx]
    if not ref_indices:
        return candidates, []

    X_ref = X[:, ref_indices]

    for candidate in candidates:
        if candidate not in gene_to_idx:
            continue

        candidate_expr = X[:, gene_to_idx[candidate]]

        # Check correlation with all reference markers
        is_redundant = False
        for i, ref_idx in enumerate(ref_indices):
            rho, _ = spearmanr(candidate_expr, X[:, ref_idx])
            if abs(rho) > threshold:
                is_redundant = True
                logger.debug(
                    f"  {candidate} redundant with {reference_markers[i]} (r={rho:.3f})"
                )
                break

        if is_redundant:
            excluded.append(candidate)
        else:
            kept.append(candidate)

    return kept, excluded


def _run_spatial_filtering(
    adata: sc.AnnData,
    genes: List[str],
    morans_k: int = 8,
    smooth_k: int = 6,
    n_permutations: int = 99,
    seed: int = 42,
) -> MarkerInterestResult:
    """
    Run Module 1 spatial interest detection on gene subset.

    Args:
        adata: AnnData with spatial coordinates in obsm['spatial'].
        genes: Genes to analyze.
        morans_k: Neighbors for Moran's I.
        smooth_k: Neighbors for smoothing.
        n_permutations: Permutations for p-value.
        seed: Random seed.

    Returns:
        MarkerInterestResult with spatial interest scores.
    """
    # Subset to genes present in data
    genes_present = [g for g in genes if g in adata.var_names]
    if not genes_present:
        raise ValueError("None of the specified genes found in data")

    adata_subset = adata[:, genes_present].copy()

    # Get expression matrix
    X = adata_subset.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    # Get coordinates
    coords = adata_subset.obsm["spatial"]

    logger.info(f"Running spatial interest detection on {len(genes_present)} genes...")

    result = identify_interesting_markers(
        X=X,
        coords=coords,
        marker_names=genes_present,
        morans_k=morans_k,
        smooth_k=smooth_k,
        morans_n_perm=n_permutations,
        seed=seed,
        verbose=True,
    )

    return result


# ============================================================================
# Main API
# ============================================================================

def select_rna_markers(
    adata: sc.AnnData,
    mode: Union[MarkerMode, str] = MarkerMode.HYBRID,
    curated_markers: Optional[List[str]] = None,
    n_hvgs: int = 2000,
    max_discovered: int = 50,
    max_total: int = 100,
    redundancy_threshold: float = 0.7,
    morans_k: int = 8,
    smooth_k: int = 6,
    n_permutations: int = 99,
    strict_spatial_threshold: bool = True,
    seed: int = 42,
) -> RNAMarkerSelectionResult:
    """
    Select RNA markers for CITEgeist pipeline.

    This function selects spatially-informative marker genes that can
    substitute for protein markers in the CITEgeist pipeline. It supports
    three modes:

    1. CURATED_ONLY: Use only user-provided marker genes
    2. HYBRID (recommended): Start with curated markers, discover additional
       spatially-interesting genes that are non-redundant
    3. AUTODISCOVERY: Fully automated selection from HVGs (exploratory)

    Args:
        adata: AnnData object with:
            - .X: Gene expression matrix
            - .var_names: Gene names
            - .obsm['spatial']: Spatial coordinates
        mode: Selection strategy (MarkerMode enum or string).
        curated_markers: User-provided marker gene list.
            If None, uses DEFAULT_CURATED_MARKERS.
        n_hvgs: Number of highly variable genes to consider for discovery.
        max_discovered: Maximum auto-discovered markers to add (hybrid/auto modes).
        max_total: Maximum total markers (hard cap for Module 2 tractability).
        redundancy_threshold: Spearman correlation threshold for redundancy.
            Genes with |r| > threshold with existing markers are excluded.
        morans_k: Neighbors for Moran's I in spatial filtering.
        smooth_k: Neighbors for smoothing in spatial filtering.
        n_permutations: Permutations for spatial interest p-values.
        strict_spatial_threshold: If True, use stricter thresholds for RNA
            (higher Moran's I, higher kurtosis requirements).
        seed: Random seed for reproducibility.

    Returns:
        RNAMarkerSelectionResult with selected markers and metadata.

    Example:
        >>> # Hybrid mode with custom markers
        >>> result = select_rna_markers(
        ...     adata,
        ...     mode="hybrid",
        ...     curated_markers=["CD3E", "CD68", "EPCAM", "COL1A1"],
        ...     max_discovered=30,
        ... )
        >>> print(result.summary())
        >>>
        >>> # Use selected markers in Module 1
        >>> adata_markers = adata[:, result.selected_markers]
    """
    # Parse mode
    if isinstance(mode, str):
        mode = MarkerMode(mode)

    logger.info(f"Selecting RNA markers (mode={mode.value})")

    # Get curated markers
    if curated_markers is None:
        curated_set = get_curated_markers()
    else:
        curated_set = set(curated_markers)

    # Find curated markers present in data
    curated_present = [g for g in curated_set if g in adata.var_names]
    curated_missing = curated_set - set(curated_present)

    if curated_missing:
        logger.info(f"  {len(curated_missing)} curated markers not in data: {list(curated_missing)[:5]}...")
    logger.info(f"  {len(curated_present)} curated markers found in data")

    # Build gene-to-index mapping
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

    # Get expression matrix
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    # Initialize result containers
    discovered_markers = []
    excluded_redundant = []
    spatial_interest = None
    marker_metadata = {}

    # Mode-specific logic
    if mode == MarkerMode.CURATED_ONLY:
        # Just use curated markers
        selected = curated_present[:max_total]

    elif mode in (MarkerMode.HYBRID, MarkerMode.AUTODISCOVERY):
        # Get HVGs for discovery
        logger.info(f"  Computing {n_hvgs} highly variable genes...")

        # Work on a copy to avoid modifying original
        adata_temp = adata.copy()

        # Compute HVGs if not already done
        if "highly_variable" not in adata_temp.var.columns:
            sc.pp.highly_variable_genes(
                adata_temp,
                n_top_genes=min(n_hvgs, adata_temp.n_vars),
                flavor="seurat_v3" if "counts" in adata_temp.layers else "seurat",
            )

        hvgs = list(adata_temp.var_names[adata_temp.var["highly_variable"]])
        logger.info(f"  Found {len(hvgs)} HVGs")

        # For hybrid mode, exclude curated markers from HVG discovery
        if mode == MarkerMode.HYBRID:
            hvgs = [g for g in hvgs if g not in curated_set]
            logger.info(f"  {len(hvgs)} HVGs after excluding curated markers")

        # Run spatial filtering on HVGs
        if hvgs:
            spatial_interest = _run_spatial_filtering(
                adata_temp,
                genes=hvgs,
                morans_k=morans_k,
                smooth_k=smooth_k,
                n_permutations=n_permutations,
                seed=seed,
            )

            # Get spatially interesting genes
            if strict_spatial_threshold:
                # Stricter filtering for RNA: require BOTH kurtosis AND Moran's I
                # (protein uses OR logic)
                interesting_genes = [
                    m.name for m in spatial_interest.markers
                    if m.passed_kurtosis and m.passed_morans and m.passed_gmm
                ]
                logger.info(
                    f"  {len(interesting_genes)} genes pass strict spatial filter "
                    "(kurtosis AND Moran's I)"
                )
            else:
                interesting_genes = spatial_interest.interesting_markers
                logger.info(f"  {len(interesting_genes)} genes pass spatial filter")

            # Sort by interest score
            gene_scores = {m.name: m.interest_score for m in spatial_interest.markers}
            interesting_genes = sorted(
                interesting_genes,
                key=lambda g: gene_scores.get(g, 0),
                reverse=True,
            )

            # Store metadata
            for m in spatial_interest.markers:
                marker_metadata[m.name] = {
                    "interest_score": m.interest_score,
                    "morans_i": m.morans_i,
                    "kurtosis": m.kurtosis,
                    "gmm_snr": m.gmm_snr,
                }

            # Filter redundant genes
            reference = curated_present if mode == MarkerMode.HYBRID else []

            candidates_to_check = interesting_genes[:max_discovered * 2]  # Check more than needed
            non_redundant, redundant = _filter_redundant_markers(
                candidates=candidates_to_check,
                reference_markers=reference + discovered_markers,
                X=X,
                gene_to_idx=gene_to_idx,
                threshold=redundancy_threshold,
            )

            # Take top non-redundant up to max_discovered
            discovered_markers = non_redundant[:max_discovered]
            excluded_redundant = redundant

            logger.info(f"  {len(discovered_markers)} non-redundant markers discovered")
            logger.info(f"  {len(excluded_redundant)} markers excluded for redundancy")

        # Combine curated and discovered
        if mode == MarkerMode.HYBRID:
            selected = curated_present + discovered_markers
        else:  # AUTODISCOVERY
            selected = discovered_markers

        # Apply total cap
        selected = selected[:max_total]

    else:
        raise ValueError(f"Unknown mode: {mode}")

    # Add metadata for curated markers (spatial scores if computed)
    for marker in curated_present:
        if marker not in marker_metadata:
            marker_metadata[marker] = {
                "interest_score": np.nan,
                "morans_i": np.nan,
                "kurtosis": np.nan,
            }

    result = RNAMarkerSelectionResult(
        selected_markers=selected,
        curated_markers=curated_present if mode != MarkerMode.AUTODISCOVERY else [],
        discovered_markers=discovered_markers,
        excluded_redundant=excluded_redundant,
        mode=mode,
        spatial_interest=spatial_interest,
        marker_metadata=marker_metadata,
    )

    logger.info(result.summary())

    return result


def select_rna_markers_from_adata(
    adata: sc.AnnData,
    mode: str = "hybrid",
    tissue_type: Optional[str] = None,
    custom_markers: Optional[List[str]] = None,
    **kwargs,
) -> RNAMarkerSelectionResult:
    """
    Convenience wrapper for select_rna_markers with sensible defaults.

    Args:
        adata: AnnData with gene expression and spatial coordinates.
        mode: "curated", "hybrid", or "autodiscovery".
        tissue_type: Optional tissue type for tissue-specific markers.
        custom_markers: Optional custom marker list (overrides tissue_type).
        **kwargs: Additional arguments passed to select_rna_markers.

    Returns:
        RNAMarkerSelectionResult.
    """
    if custom_markers is not None:
        curated = custom_markers
    else:
        curated = list(get_curated_markers(tissue_type=tissue_type))

    return select_rna_markers(
        adata=adata,
        mode=mode,
        curated_markers=curated,
        **kwargs,
    )


# ============================================================================
# Validation utilities
# ============================================================================

def validate_marker_spatial_quality(
    adata: sc.AnnData,
    markers: List[str],
    min_morans_i: float = 0.1,
    min_expressing_fraction: float = 0.05,
) -> pd.DataFrame:
    """
    Validate spatial quality of selected markers.

    Checks that markers have:
    1. Sufficient expression (not too sparse)
    2. Spatial structure (Moran's I > threshold)

    Args:
        adata: AnnData with expression and spatial coordinates.
        markers: Markers to validate.
        min_morans_i: Minimum Moran's I for "good" spatial signal.
        min_expressing_fraction: Minimum fraction of spots with expression.

    Returns:
        DataFrame with quality metrics and pass/fail status.
    """
    markers_present = [m for m in markers if m in adata.var_names]

    # Run spatial interest analysis
    spatial_result = _run_spatial_filtering(
        adata,
        genes=markers_present,
        n_permutations=99,
    )

    # Get expression fractions
    X = adata[:, markers_present].X
    if hasattr(X, "toarray"):
        X = X.toarray()

    expressing_fractions = (X > 0).mean(axis=0)

    # Build validation table
    records = []
    for i, marker in enumerate(markers_present):
        m_info = next((m for m in spatial_result.markers if m.name == marker), None)

        morans_i = m_info.morans_i if m_info else np.nan
        expr_frac = expressing_fractions[i]

        passed = (
            morans_i >= min_morans_i and
            expr_frac >= min_expressing_fraction
        )

        records.append({
            "marker": marker,
            "morans_i": morans_i,
            "expressing_fraction": expr_frac,
            "passed_morans": morans_i >= min_morans_i,
            "passed_expression": expr_frac >= min_expressing_fraction,
            "passed": passed,
        })

    df = pd.DataFrame(records)

    n_passed = df["passed"].sum()
    logger.info(f"Marker validation: {n_passed}/{len(df)} markers passed quality checks")

    return df


# ============================================================================
# Profile comparison / validation
# ============================================================================

@dataclass
class ProfileComparisonResult:
    """Results from comparing protein vs RNA discovered profiles."""

    # Profile-level metrics
    n_protein_profiles: int
    n_rna_profiles: int
    profile_jaccard_matrix: pd.DataFrame  # (n_protein × n_rna) raw Jaccard similarities
    best_matches: Dict[str, Tuple[str, float]]  # protein_profile -> (best_rna_profile, jaccard)

    # Semantic matching (using biological equivalence)
    semantic_jaccard_matrix: Optional[pd.DataFrame] = None  # Jaccard with equivalence mapping
    best_semantic_matches: Optional[Dict[str, Tuple[str, float]]] = None
    mean_semantic_jaccard: float = 0.0
    semantic_matched_profiles: int = 0

    # Spatial concordance
    spatial_correlation_matrix: Optional[pd.DataFrame] = None  # Correlation of profile scores
    mean_spatial_correlation: float = 0.0
    best_spatial_matches: Optional[Dict[str, Tuple[str, float]]] = None

    # Assignment agreement (if assignments provided)
    assignment_agreement: Optional[float] = None  # Fraction of spots with same top profile
    confusion_matrix: Optional[pd.DataFrame] = None

    # Summary metrics
    mean_best_jaccard: float = 0.0
    matched_profiles: int = 0  # Number of protein profiles with good RNA match (raw Jaccard)

    # Detailed match info
    match_details: Optional[List[Dict]] = None

    def summary(self) -> str:
        """Human-readable summary."""
        lines = [
            "Profile Comparison: Protein vs RNA",
            "=" * 50,
            f"Protein profiles: {self.n_protein_profiles}",
            f"RNA profiles: {self.n_rna_profiles}",
            "",
            "--- Raw Jaccard (marker name overlap) ---",
            f"Mean best Jaccard: {self.mean_best_jaccard:.3f}",
            f"Matched profiles (Jaccard > 0.3): {self.matched_profiles}",
        ]

        if self.semantic_jaccard_matrix is not None:
            lines.extend([
                "",
                "--- Semantic Jaccard (biological equivalence) ---",
                f"Mean semantic Jaccard: {self.mean_semantic_jaccard:.3f}",
                f"Matched profiles (semantic > 0.3): {self.semantic_matched_profiles}",
            ])

        if self.spatial_correlation_matrix is not None:
            lines.extend([
                "",
                "--- Spatial Concordance ---",
                f"Mean spatial correlation: {self.mean_spatial_correlation:.3f}",
            ])

        if self.assignment_agreement is not None:
            lines.append(f"Assignment agreement: {self.assignment_agreement:.1%}")

        # Show best matches with semantic info
        if self.match_details:
            lines.extend(["", "Best matches (by spatial correlation):"])
            for match in self.match_details[:10]:  # Top 10
                p_name = match["protein_profile"]
                r_name = match["rna_profile"]
                spatial = match.get("spatial_corr", 0)
                semantic = match.get("semantic_jaccard", 0)
                raw = match.get("raw_jaccard", 0)
                lines.append(
                    f"  {p_name} -> {r_name} "
                    f"(spatial={spatial:.3f}, semantic={semantic:.3f}, raw={raw:.3f})"
                )
                if match.get("equivalent_markers"):
                    lines.append(f"    Equivalent markers: {match['equivalent_markers']}")

        return "\n".join(lines)


def _compute_profile_scores(
    X: NDArray[np.floating],
    profiles: List[List[str]],
    marker_names: List[str],
    profile_names: Optional[List[str]] = None,
) -> NDArray[np.floating]:
    """
    Compute profile activity scores for each spot.

    Score = mean expression of profile markers.

    Args:
        X: Expression matrix (n_spots, n_markers).
        profiles: List of marker lists defining each profile.
        marker_names: Names corresponding to X columns.
        profile_names: Optional names for profiles.

    Returns:
        Score matrix (n_spots, n_profiles).
    """
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}
    n_spots = X.shape[0]
    n_profiles = len(profiles)

    scores = np.zeros((n_spots, n_profiles))

    for p_idx, profile_markers in enumerate(profiles):
        marker_indices = [marker_to_idx[m] for m in profile_markers if m in marker_to_idx]
        if marker_indices:
            scores[:, p_idx] = X[:, marker_indices].mean(axis=1)

    return scores


def _compute_semantic_jaccard(
    protein_markers: List[str],
    rna_markers: List[str],
) -> Tuple[float, List[str]]:
    """
    Compute Jaccard similarity using biological equivalence mapping.

    Instead of matching marker names directly, this maps protein antibody
    names to their corresponding RNA gene names for comparison.

    Args:
        protein_markers: List of protein marker names.
        rna_markers: List of RNA marker names.

    Returns:
        Tuple of (semantic_jaccard, list_of_equivalent_markers_found).
    """
    # Convert protein markers to their RNA equivalents
    protein_as_rna = set()
    for p_marker in protein_markers:
        equivalents = get_equivalent_markers(p_marker, source="protein")
        protein_as_rna.update(equivalents)

    rna_set = set(rna_markers)

    # Compute intersection and union
    intersection = protein_as_rna & rna_set
    union = protein_as_rna | rna_set

    jaccard = len(intersection) / len(union) if union else 0.0

    return jaccard, list(intersection)


def compare_protein_vs_rna_profiles(
    protein_profiles: List[List[str]],
    rna_profiles: List[List[str]],
    X_protein: Optional[NDArray[np.floating]] = None,
    X_rna: Optional[NDArray[np.floating]] = None,
    protein_marker_names: Optional[List[str]] = None,
    rna_marker_names: Optional[List[str]] = None,
    coords: Optional[NDArray[np.floating]] = None,
    protein_profile_names: Optional[List[str]] = None,
    rna_profile_names: Optional[List[str]] = None,
    jaccard_match_threshold: float = 0.3,
    use_semantic_matching: bool = True,
) -> ProfileComparisonResult:
    """
    Compare profiles discovered from protein vs RNA markers.

    This validates whether RNA-based profile discovery recovers
    similar structure to protein-based discovery.

    Computes:
    1. Raw Jaccard similarity (marker name overlap)
    2. Semantic Jaccard similarity (using biological equivalence mapping)
    3. Spatial concordance (correlation of profile scores across spots)
    4. Assignment agreement (optional, if expression data provided)

    The semantic Jaccard uses a biological equivalence mapping (e.g., CD20 -> MS4A1)
    to properly compare protein antibodies with RNA genes that encode the same protein.

    Args:
        protein_profiles: List of marker lists from protein analysis.
        rna_profiles: List of marker lists from RNA analysis.
        X_protein: Optional protein expression matrix for spatial analysis.
        X_rna: Optional RNA expression matrix for spatial analysis.
        protein_marker_names: Names for protein markers (columns of X_protein).
        rna_marker_names: Names for RNA markers (columns of X_rna).
        coords: Spatial coordinates (required for spatial concordance).
        protein_profile_names: Optional names for protein profiles.
        rna_profile_names: Optional names for RNA profiles.
        jaccard_match_threshold: Minimum Jaccard for a "good" match.
        use_semantic_matching: Whether to compute semantic Jaccard.

    Returns:
        ProfileComparisonResult with comparison metrics.

    Example:
        >>> comparison = compare_protein_vs_rna_profiles(
        ...     protein_profiles=[["CD20", "CD45RA"], ["CD68", "CD163"]],
        ...     rna_profiles=[["MS4A1", "PTPRC"], ["CD68", "CD163"]],
        ... )
        >>> # Semantic Jaccard will recognize CD20=MS4A1 and CD45RA=PTPRC
        >>> print(comparison.mean_semantic_jaccard)  # Higher than raw
    """
    n_protein = len(protein_profiles)
    n_rna = len(rna_profiles)

    # Generate profile names if not provided
    if protein_profile_names is None:
        protein_profile_names = [f"Protein_{i}" for i in range(n_protein)]
    if rna_profile_names is None:
        rna_profile_names = [f"RNA_{i}" for i in range(n_rna)]

    # 1. Compute raw Jaccard similarity (marker name overlap)
    jaccard_matrix = np.zeros((n_protein, n_rna))

    for i, p_markers in enumerate(protein_profiles):
        p_set = set(p_markers)
        for j, r_markers in enumerate(rna_profiles):
            r_set = set(r_markers)
            intersection = len(p_set & r_set)
            union = len(p_set | r_set)
            jaccard_matrix[i, j] = intersection / union if union > 0 else 0.0

    jaccard_df = pd.DataFrame(
        jaccard_matrix,
        index=protein_profile_names,
        columns=rna_profile_names,
    )

    # Find best raw matches
    best_matches = {}
    for i, p_name in enumerate(protein_profile_names):
        best_j = np.argmax(jaccard_matrix[i, :])
        best_jaccard = jaccard_matrix[i, best_j]
        best_matches[p_name] = (rna_profile_names[best_j], float(best_jaccard))

    mean_best_jaccard = np.mean([v[1] for v in best_matches.values()])
    matched_profiles = sum(1 for v in best_matches.values() if v[1] >= jaccard_match_threshold)

    # 2. Compute semantic Jaccard similarity (using biological equivalence)
    semantic_jaccard_df = None
    best_semantic_matches = None
    mean_semantic_jaccard = 0.0
    semantic_matched_profiles = 0
    equivalent_markers_found = {}  # Store which equivalents matched

    if use_semantic_matching:
        semantic_matrix = np.zeros((n_protein, n_rna))

        for i, p_markers in enumerate(protein_profiles):
            for j, r_markers in enumerate(rna_profiles):
                sem_jaccard, equiv_markers = _compute_semantic_jaccard(p_markers, r_markers)
                semantic_matrix[i, j] = sem_jaccard
                if equiv_markers:
                    key = (protein_profile_names[i], rna_profile_names[j])
                    equivalent_markers_found[key] = equiv_markers

        semantic_jaccard_df = pd.DataFrame(
            semantic_matrix,
            index=protein_profile_names,
            columns=rna_profile_names,
        )

        # Find best semantic matches
        best_semantic_matches = {}
        for i, p_name in enumerate(protein_profile_names):
            best_j = np.argmax(semantic_matrix[i, :])
            best_jaccard = semantic_matrix[i, best_j]
            best_semantic_matches[p_name] = (rna_profile_names[best_j], float(best_jaccard))

        mean_semantic_jaccard = np.mean([v[1] for v in best_semantic_matches.values()])
        semantic_matched_profiles = sum(
            1 for v in best_semantic_matches.values() if v[1] >= jaccard_match_threshold
        )

    # 3. Compute spatial concordance (if expression data provided)
    spatial_corr_matrix = None
    mean_spatial_corr = 0.0
    best_spatial_matches = None
    protein_scores = None
    rna_scores = None

    if X_protein is not None and X_rna is not None and coords is not None:
        # Ensure arrays
        if hasattr(X_protein, "toarray"):
            X_protein = X_protein.toarray()
        if hasattr(X_rna, "toarray"):
            X_rna = X_rna.toarray()

        X_protein = np.asarray(X_protein, dtype=np.float64)
        X_rna = np.asarray(X_rna, dtype=np.float64)

        # Compute profile scores
        protein_scores = _compute_profile_scores(
            X_protein, protein_profiles, protein_marker_names
        )
        rna_scores = _compute_profile_scores(
            X_rna, rna_profiles, rna_marker_names
        )

        # Compute correlation matrix
        spatial_corr = np.zeros((n_protein, n_rna))
        for i in range(n_protein):
            for j in range(n_rna):
                if protein_scores[:, i].std() > 1e-10 and rna_scores[:, j].std() > 1e-10:
                    corr, _ = spearmanr(protein_scores[:, i], rna_scores[:, j])
                    spatial_corr[i, j] = corr if not np.isnan(corr) else 0.0
                else:
                    spatial_corr[i, j] = 0.0

        spatial_corr_matrix = pd.DataFrame(
            spatial_corr,
            index=protein_profile_names,
            columns=rna_profile_names,
        )

        # Find best spatial matches
        best_spatial_matches = {}
        for i, p_name in enumerate(protein_profile_names):
            best_j = np.argmax(spatial_corr[i, :])
            best_corr = spatial_corr[i, best_j]
            best_spatial_matches[p_name] = (rna_profile_names[best_j], float(best_corr))

        # Mean of best spatial correlations
        mean_spatial_corr = float(np.mean([v[1] for v in best_spatial_matches.values()]))

    # 4. Assignment agreement (optional)
    assignment_agreement = None
    confusion_matrix = None

    if protein_scores is not None and rna_scores is not None:
        # Assign spots to profiles
        protein_assignments = np.argmax(protein_scores, axis=1)
        rna_assignments = np.argmax(rna_scores, axis=1)

        # Use spatial correlation to determine best match mapping
        match_source = best_spatial_matches if best_spatial_matches else best_matches
        corr_threshold = 0.3  # Correlation threshold for matching

        rna_to_protein_map = {}
        for i, p_name in enumerate(protein_profile_names):
            best_rna, score = match_source[p_name]
            if score >= corr_threshold:
                j = rna_profile_names.index(best_rna)
                rna_to_protein_map[j] = i

        # Compute agreement only for matched profiles
        if rna_to_protein_map:
            agreement_count = 0
            total_count = 0
            for spot_idx in range(len(protein_assignments)):
                p_assign = protein_assignments[spot_idx]
                r_assign = rna_assignments[spot_idx]
                if r_assign in rna_to_protein_map:
                    total_count += 1
                    if rna_to_protein_map[r_assign] == p_assign:
                        agreement_count += 1

            if total_count > 0:
                assignment_agreement = agreement_count / total_count

        # Build confusion matrix
        conf = np.zeros((n_protein, n_rna), dtype=int)
        for spot_idx in range(len(protein_assignments)):
            conf[protein_assignments[spot_idx], rna_assignments[spot_idx]] += 1

        confusion_matrix = pd.DataFrame(
            conf,
            index=protein_profile_names,
            columns=rna_profile_names,
        )

    # 5. Build detailed match info (sorted by spatial correlation)
    match_details = []
    for i, p_name in enumerate(protein_profile_names):
        # Get best RNA match by spatial correlation
        if best_spatial_matches:
            best_rna, spatial = best_spatial_matches[p_name]
        else:
            best_rna = best_matches[p_name][0]
            spatial = 0.0

        j = rna_profile_names.index(best_rna)

        detail = {
            "protein_profile": p_name,
            "protein_markers": protein_profiles[i],
            "rna_profile": best_rna,
            "rna_markers": rna_profiles[j],
            "raw_jaccard": jaccard_matrix[i, j],
            "spatial_corr": spatial,
        }

        if use_semantic_matching:
            detail["semantic_jaccard"] = semantic_matrix[i, j]
            key = (p_name, best_rna)
            detail["equivalent_markers"] = equivalent_markers_found.get(key, [])

        match_details.append(detail)

    # Sort by spatial correlation (descending)
    match_details.sort(key=lambda x: x.get("spatial_corr", 0), reverse=True)

    return ProfileComparisonResult(
        n_protein_profiles=n_protein,
        n_rna_profiles=n_rna,
        profile_jaccard_matrix=jaccard_df,
        best_matches=best_matches,
        semantic_jaccard_matrix=semantic_jaccard_df,
        best_semantic_matches=best_semantic_matches,
        mean_semantic_jaccard=mean_semantic_jaccard,
        semantic_matched_profiles=semantic_matched_profiles,
        spatial_correlation_matrix=spatial_corr_matrix,
        mean_spatial_correlation=mean_spatial_corr,
        best_spatial_matches=best_spatial_matches,
        assignment_agreement=assignment_agreement,
        confusion_matrix=confusion_matrix,
        mean_best_jaccard=mean_best_jaccard,
        matched_profiles=matched_profiles,
        match_details=match_details,
    )


def validate_rna_profiles_against_protein(
    adata: sc.AnnData,
    protein_profiles: List[List[str]],
    rna_marker_result: RNAMarkerSelectionResult,
    rna_pipeline_config: Optional[Dict] = None,
) -> ProfileComparisonResult:
    """
    Run RNA pipeline and compare to known protein profiles.

    Convenience function for validation experiments.

    Args:
        adata: AnnData with both protein and RNA expression.
        protein_profiles: Ground truth profiles from protein analysis.
        rna_marker_result: Result from select_rna_markers().
        rna_pipeline_config: Optional config for RNA pipeline.

    Returns:
        ProfileComparisonResult.
    """
    from .spatial_colocalization import (
        analyze_marker_colocalization,
        discover_profiles,
    )

    # Get RNA marker expression
    rna_markers = rna_marker_result.selected_markers
    adata_rna = adata[:, [m for m in rna_markers if m in adata.var_names]].copy()

    X_rna = adata_rna.X
    if hasattr(X_rna, "toarray"):
        X_rna = X_rna.toarray()
    X_rna = np.asarray(X_rna, dtype=np.float64)

    coords = adata_rna.obsm["spatial"]
    marker_names = list(adata_rna.var_names)

    # Run Module 1 on RNA markers
    from .marker_interest import identify_interesting_markers

    interest_result = identify_interesting_markers(
        X=X_rna,
        coords=coords,
        marker_names=marker_names,
        verbose=False,
    )

    # Run Module 2
    coloc_result = analyze_marker_colocalization(
        X=X_rna,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interest_result.interesting_markers,
        verbose=False,
    )

    profile_result = discover_profiles(
        colocalization_result=coloc_result,
        verbose=False,
    )

    rna_profiles = profile_result.profiles

    # Get protein expression for comparison
    all_protein_markers = set()
    for profile in protein_profiles:
        all_protein_markers.update(profile)

    protein_markers_present = [m for m in all_protein_markers if m in adata.var_names]

    if protein_markers_present:
        adata_protein = adata[:, protein_markers_present].copy()
        X_protein = adata_protein.X
        if hasattr(X_protein, "toarray"):
            X_protein = X_protein.toarray()
        X_protein = np.asarray(X_protein, dtype=np.float64)
        protein_marker_names = list(adata_protein.var_names)
    else:
        X_protein = None
        protein_marker_names = None

    # Compare profiles
    comparison = compare_protein_vs_rna_profiles(
        protein_profiles=protein_profiles,
        rna_profiles=rna_profiles,
        X_protein=X_protein,
        X_rna=X_rna,
        protein_marker_names=protein_marker_names,
        rna_marker_names=marker_names,
        coords=coords,
    )

    return comparison


if __name__ == "__main__":
    # Simple test with synthetic data
    import logging
    logging.basicConfig(level=logging.INFO)

    print("=" * 60)
    print("RNA Marker Selection Module Test")
    print("=" * 60)

    # Create synthetic spatial data
    np.random.seed(42)
    n_spots = 500
    n_genes = 1000

    # Random expression
    X = np.random.exponential(1.0, (n_spots, n_genes))

    # Add spatial structure to some genes
    coords = np.random.rand(n_spots, 2) * 100

    # Make genes 0-4 spatially structured (higher in left half)
    left_mask = coords[:, 0] < 50
    X[left_mask, 0:5] += 5.0

    # Make genes 5-9 spatially structured (higher in top half)
    top_mask = coords[:, 1] > 50
    X[top_mask, 5:10] += 5.0

    # Create AnnData
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    # Include some "curated" marker names
    gene_names[0] = "CD3E"
    gene_names[1] = "CD68"
    gene_names[5] = "EPCAM"
    gene_names[6] = "COL1A1"

    adata = sc.AnnData(X)
    adata.var_names = gene_names
    adata.obsm["spatial"] = coords

    # Mark HVGs manually for this test
    adata.var["highly_variable"] = False
    adata.var.loc[adata.var_names[:100], "highly_variable"] = True

    print("\nTest 1: Curated-only mode")
    result1 = select_rna_markers(
        adata,
        mode="curated",
        curated_markers=["CD3E", "CD68", "EPCAM", "COL1A1", "MISSING_GENE"],
    )
    print(result1.summary())
    print(f"Selected: {result1.selected_markers}")

    print("\n" + "=" * 60)
    print("Test 2: Hybrid mode")
    result2 = select_rna_markers(
        adata,
        mode="hybrid",
        curated_markers=["CD3E", "CD68"],
        max_discovered=10,
        n_permutations=49,  # Faster for test
    )
    print(result2.summary())
    print(f"Selected: {result2.selected_markers}")

    print("\n" + "=" * 60)
    print("Test 3: Autodiscovery mode")
    result3 = select_rna_markers(
        adata,
        mode="autodiscovery",
        max_discovered=15,
        n_permutations=49,
    )
    print(result3.summary())
    print(f"Selected: {result3.selected_markers}")

    print("\nRNA marker selection module test complete!")
