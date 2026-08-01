"""
Module 4: Protein-Anchored Spatial Transcriptomic Program Discovery.

Discovers gene expression programs from DECONVOLVED cell-type-specific expression
layers output by Module 3, validated against protein profiles from Module 2.

This fills a gap in existing methods (NSF, STAMP, SpaTM) which work on
transcriptomics alone - CITEgeist leverages spatial protein data AND
cell-type-specific deconvolved expression for interpretable program discovery.

Pipeline integration:
    Module 1 (marker_interest) -> Filter spatially-variable proteins
    Module 2 (spatial_colocalization) -> Discover protein profiles
    Module 3 (qp_solver) -> Cell type proportions + DECONVOLVED GEX LAYERS
    Module 4 (THIS MODULE) -> Programs from deconvolved layers

Key architectural change (v3):
    - Module 4 now uses DECONVOLVED gene expression layers from Module 3
    - Each layer contains cell-type-specific expression (what that cell type expresses)
    - No need for contrastive subtraction - data is already cell-type separated
    - Helper function stacks layers for unified analysis
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from numpy.typing import NDArray
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import NMF

from ..deconvolution.qp_solver import build_spatial_laplacian
from ..discovery.spatial_colocalization import ProfileDiscoveryResult

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

    # Region analysis fields (populated by analyze_program_regions)
    region_enrichment: Optional[Dict[str, float]] = None
    """Mean program activity per region (region_value -> mean_activity)."""

    region_specificity: Optional[float] = None
    """Specificity score: 0 = uniform across regions, 1 = completely region-specific."""

    region_pvalue: Optional[float] = None
    """P-value for region enrichment (Mann-Whitney U or Kruskal-Wallis)."""

    enriched_region: Optional[str] = None
    """Region with highest mean activity (if significantly enriched)."""


@dataclass
class JointProgram:
    """A spatial program discovered from joint analysis of all cell types."""

    program_id: int
    """Unique identifier for this program."""

    top_genes: List[str]
    """Top N genes by loading magnitude."""

    gene_loadings: Dict[str, float]
    """Gene name -> loading value for top genes."""

    variance_explained: float
    """Fraction of variance explained by this program."""

    spatial_moran_i: float
    """Moran's I spatial autocorrelation of program activity."""

    spatial_moran_pvalue: float
    """P-value for Moran's I test."""

    mean_activity: float
    """Mean program activity across all spots."""

    active_spots_fraction: float
    """Fraction of spots with above-median activity."""

    cell_type_enrichments: Dict[str, float]
    """Cell type -> enrichment score (correlation with proportions)."""

    primary_cell_type: str
    """Cell type with highest enrichment."""

    secondary_cell_type: Optional[str]
    """Second highest cell type if interaction program."""

    interaction_score: float
    """Score 0-1 indicating how multi-cell-type the program is (0=single, 1=balanced)."""

    program_type: str
    """'single_celltype', 'interaction', or 'microenvironment'."""


@dataclass
class JointDiscoveryResult:
    """Results from joint program discovery across all cell types."""

    programs: List[JointProgram]
    """List of discovered programs with cell type assignments."""

    W: NDArray[np.floating]
    """Gene loadings matrix (n_genes, K_programs)."""

    H: NDArray[np.floating]
    """Program activities matrix (K_programs, n_spots)."""

    gene_names: List[str]
    """Gene names corresponding to W rows."""

    cell_type_names: List[str]
    """Cell types included in analysis."""

    n_spots: int
    """Number of spots analyzed."""

    reconstruction_error: float
    """NMF reconstruction error."""

    H_by_celltype: Optional[Dict[str, NDArray]] = None
    """Program activities split by cell type (from unstacking)."""

    parameters: Dict[str, Any] = field(default_factory=dict)
    """Parameters used for discovery."""

    def summary(self) -> str:
        """Return summary string."""
        n_single = sum(1 for p in self.programs if p.program_type == "single_celltype")
        n_interaction = sum(1 for p in self.programs if p.program_type == "interaction")
        n_micro = sum(1 for p in self.programs if p.program_type == "microenvironment")

        lines = [
            "Joint Program Discovery Results",
            "=" * 40,
            f"Total programs: {len(self.programs)}",
            f"  Single cell-type: {n_single}",
            f"  Interaction: {n_interaction}",
            f"  Microenvironment: {n_micro}",
            f"Spots: {self.n_spots}",
            f"Genes: {len(self.gene_names)}",
            f"Cell types: {', '.join(self.cell_type_names)}",
            f"Reconstruction error: {self.reconstruction_error:.4f}",
        ]
        return "\n".join(lines)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert programs to DataFrame."""
        records = []
        for prog in self.programs:
            records.append(
                {
                    "program_id": prog.program_id,
                    "top_genes": ", ".join(prog.top_genes[:10]),
                    "variance_explained": prog.variance_explained,
                    "spatial_moran_i": prog.spatial_moran_i,
                    "spatial_moran_pvalue": prog.spatial_moran_pvalue,
                    "mean_activity": prog.mean_activity,
                    "active_spots_fraction": prog.active_spots_fraction,
                    "primary_cell_type": prog.primary_cell_type,
                    "secondary_cell_type": prog.secondary_cell_type,
                    "interaction_score": prog.interaction_score,
                    "program_type": prog.program_type,
                }
            )
        return pd.DataFrame(records)


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
            records.append(
                {
                    "anchor": self.anchor_name,
                    "program_id": p.program_id,
                    "n_top_genes": len(p.top_genes),
                    "top_genes": ", ".join(p.top_genes[:10]),
                    "variance_explained": p.variance_explained,
                    "spatial_moran_i": p.spatial_moran_i,
                    "spatial_moran_pvalue": p.spatial_moran_pvalue,
                    "mean_activity": p.mean_activity,
                    "active_spots_fraction": p.active_spots_fraction,
                }
            )
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
            "Anchored Program Discovery Results",
            "=" * 40,
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
# Module 4b: Bivariate Program Relationships
# =============================================================================


@dataclass
class ProgramPairRelationship:
    """Bivariate spatial relationship between two programs."""

    anchor1: str
    """Cell type 1."""

    program1_id: int
    """Program index in anchor1."""

    anchor2: str
    """Cell type 2."""

    program2_id: int
    """Program index in anchor2."""

    bivariate_morans_i: float
    """Spatial cross-correlation [-1, 1]."""

    bivariate_pvalue: float
    """Permutation test p-value."""

    pearson_correlation: float
    """Non-spatial correlation."""

    pearson_pvalue: float
    """Pearson p-value."""

    n_spots_used: int
    """Spots where both programs are active."""

    top_genes_overlap: List[str]
    """Shared top genes (if any)."""

    relationship_type: str
    """'co-localized', 'exclusive', or 'independent'."""

    def __repr__(self) -> str:
        return (
            f"ProgramPairRelationship({self.anchor1}_prog{self.program1_id} ↔ "
            f"{self.anchor2}_prog{self.program2_id}: I={self.bivariate_morans_i:.3f}, "
            f"type={self.relationship_type})"
        )


@dataclass
class BivariateProgramResult:
    """Module 4b results: bivariate relationships between programs."""

    significant_pairs: List[ProgramPairRelationship]
    """Pairs with FDR-corrected significant spatial relationships."""

    all_pairs: List[ProgramPairRelationship]
    """All tested pairs."""

    n_programs_total: int
    """Total number of programs analyzed."""

    n_pairs_tested: int
    """Number of pairs tested."""

    n_significant: int
    """Number of significant pairs after FDR correction."""

    fdr_threshold: float
    """FDR threshold used."""

    parameters: Dict[str, Any] = field(default_factory=dict)
    """Parameters used for analysis."""

    def summary(self) -> str:
        """Return a summary string."""
        lines = [
            "Module 4b: Bivariate Program Relationships",
            "=" * 45,
            f"Programs analyzed: {self.n_programs_total}",
            f"Pairs tested: {self.n_pairs_tested}",
            f"Significant pairs (FDR < {self.fdr_threshold}): {self.n_significant}",
            "",
        ]

        # Summarize co-localized pairs
        colocalized = [p for p in self.significant_pairs if p.relationship_type == "co-localized"]
        exclusive = [p for p in self.significant_pairs if p.relationship_type == "exclusive"]

        if colocalized:
            lines.append(f"Co-localized pairs ({len(colocalized)}):")
            for pair in colocalized[:5]:  # Top 5
                lines.append(
                    f"  {pair.anchor1}_prog{pair.program1_id} ↔ "
                    f"{pair.anchor2}_prog{pair.program2_id}: "
                    f"I={pair.bivariate_morans_i:.3f} (p={pair.bivariate_pvalue:.3f})"
                )
            if len(colocalized) > 5:
                lines.append(f"  ... and {len(colocalized) - 5} more")
            lines.append("")

        if exclusive:
            lines.append(f"Mutually exclusive pairs ({len(exclusive)}):")
            for pair in exclusive[:5]:  # Top 5
                lines.append(
                    f"  {pair.anchor1}_prog{pair.program1_id} ↔ "
                    f"{pair.anchor2}_prog{pair.program2_id}: "
                    f"I={pair.bivariate_morans_i:.3f} (p={pair.bivariate_pvalue:.3f})"
                )
            if len(exclusive) > 5:
                lines.append(f"  ... and {len(exclusive) - 5} more")

        return "\n".join(lines)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert all pairs to DataFrame."""
        records = []
        for p in self.all_pairs:
            records.append(
                {
                    "anchor1": p.anchor1,
                    "program1_id": p.program1_id,
                    "anchor2": p.anchor2,
                    "program2_id": p.program2_id,
                    "bivariate_morans_i": p.bivariate_morans_i,
                    "bivariate_pvalue": p.bivariate_pvalue,
                    "pearson_correlation": p.pearson_correlation,
                    "pearson_pvalue": p.pearson_pvalue,
                    "n_spots_used": p.n_spots_used,
                    "n_genes_overlap": len(p.top_genes_overlap),
                    "top_genes_overlap": ", ".join(p.top_genes_overlap[:10]),
                    "relationship_type": p.relationship_type,
                }
            )
        return pd.DataFrame(records)


def _build_program_pair_list(
    result: AnchoredProgramDiscoveryResult,
    include_within_anchor: bool = True,
) -> List[Tuple[str, int, str, int]]:
    """
    Generate all (anchor1, prog1, anchor2, prog2) tuples to test.

    Args:
        result: Module 4 output with programs per cell type.
        include_within_anchor: Also compare programs within same cell type.

    Returns:
        List of (anchor1, prog1_id, anchor2, prog2_id) tuples.
    """
    pairs = []
    anchors = list(result.results_by_anchor.keys())

    for i, anchor1 in enumerate(anchors):
        result1 = result.results_by_anchor[anchor1]
        n_progs1 = len(result1.programs)

        # Within-anchor pairs
        if include_within_anchor:
            for p1 in range(n_progs1):
                for p2 in range(p1 + 1, n_progs1):
                    pairs.append((anchor1, p1, anchor1, p2))

        # Cross-anchor pairs
        for anchor2 in anchors[i + 1 :]:
            result2 = result.results_by_anchor[anchor2]
            n_progs2 = len(result2.programs)

            for p1 in range(n_progs1):
                for p2 in range(n_progs2):
                    pairs.append((anchor1, p1, anchor2, p2))

    return pairs


def _classify_relationship(
    bivariate_i: float,
    _pvalue: float,
    fdr_significant: bool,
    threshold: float = 0.1,
) -> str:
    """
    Classify relationship as co-localized, exclusive, or independent.

    Args:
        bivariate_i: Bivariate Moran's I value.
        pvalue: P-value for the bivariate Moran's I.
        fdr_significant: Whether this passes FDR correction.
        threshold: Minimum |I| to consider significant.

    Returns:
        Relationship type string.
    """
    if not fdr_significant:
        return "independent"

    if bivariate_i > threshold:  # pylint: disable=no-else-return
        return "co-localized"
    elif bivariate_i < -threshold:
        return "exclusive"
    else:
        return "independent"


def _compute_gene_overlap(
    program1: SpatialProgram,
    program2: SpatialProgram,
    top_n: int = 50,
) -> List[str]:
    """
    Find shared genes between program signatures.

    Args:
        program1: First program.
        program2: Second program.
        top_n: Number of top genes to compare.

    Returns:
        List of overlapping gene names.
    """
    genes1 = set(program1.top_genes[:top_n])
    genes2 = set(program2.top_genes[:top_n])
    return sorted(genes1 & genes2)


def analyze_program_relationships(
    result: AnchoredProgramDiscoveryResult,
    adata: sc.AnnData,
    *,
    fdr_threshold: float = 0.05,
    min_bivariate_i: float = 0.1,
    n_permutations: int = 199,
    neighbor_k: int = 8,
    include_within_anchor: bool = True,
    random_state: int = 42,
) -> BivariateProgramResult:
    """
    Module 4b: Analyze bivariate spatial relationships between programs.

    Computes bivariate Moran's I between all pairs of programs to identify
    spatially correlated transcriptomic programs. This reveals:
    - Co-localized programs: Different cell type programs that peak together spatially
    - Exclusive programs: Programs that avoid each other spatially
    - Independent programs: No significant spatial relationship

    Args:
        result: Module 4 output with H matrices per cell type.
        adata: AnnData with spatial coordinates.
        fdr_threshold: FDR threshold for significance.
        min_bivariate_i: Minimum ``|bivariate Moran's I|`` to report as strong.
        n_permutations: Permutations for p-value calculation.
        neighbor_k: k-NN for spatial weights.
        include_within_anchor: Also compare programs within same cell type.
        random_state: Random seed.

    Returns:
        BivariateProgramResult with all pairwise relationships.
    """
    from scipy.spatial import cKDTree  # pylint: disable=import-outside-toplevel
    from scipy.stats import false_discovery_control  # pylint: disable=import-outside-toplevel

    logger.info("Starting Module 4b: Bivariate Program Relationships")

    # Get spatial coordinates
    coords = adata.obsm.get("spatial", None)
    if coords is None:
        raise ValueError("No spatial coordinates found in adata.obsm['spatial']")

    n_spots = coords.shape[0]

    # Build k-NN neighbor graph
    logger.info("Building %s-NN spatial neighbor graph for %s spots", neighbor_k, n_spots)
    tree = cKDTree(coords)
    _, neighbor_indices = tree.query(coords, k=neighbor_k + 1)

    # Convert to list of neighbor lists (excluding self)
    neighbors: List[List[int]] = []
    for i in range(n_spots):
        # Exclude self (first element)
        neighbors.append(neighbor_indices[i, 1:].tolist())

    # Build list of program pairs to test
    pairs = _build_program_pair_list(result, include_within_anchor)
    n_pairs = len(pairs)

    # Count total programs
    n_programs = sum(len(r.programs) for r in result.results_by_anchor.values())

    logger.info(
        "Testing %s program pairs (%s programs, %s within-anchor)",
        n_pairs,
        n_programs,
        "including" if include_within_anchor else "excluding",
    )

    # Initialize random number generator
    rng = np.random.default_rng(random_state)

    # Compute bivariate Moran's I for each pair
    all_relationships: List[ProgramPairRelationship] = []
    pvalues_for_fdr: List[float] = []

    for pair_idx, (anchor1, prog1_id, anchor2, prog2_id) in enumerate(pairs):
        if (pair_idx + 1) % 50 == 0 or pair_idx == n_pairs - 1:
            logger.info("  Processing pair %s/%s", pair_idx + 1, n_pairs)

        # Get H vectors for both programs
        result1 = result.results_by_anchor[anchor1]
        result2 = result.results_by_anchor[anchor2]

        H1 = result1.H[prog1_id, :]  # (n_spots,)
        H2 = result2.H[prog2_id, :]  # (n_spots,)

        # Find spots where both programs have activity
        active_mask = (H1 > 0) & (H2 > 0)
        n_active = active_mask.sum()

        if n_active < 20:
            # Not enough active spots for reliable analysis
            relationship = ProgramPairRelationship(
                anchor1=anchor1,
                program1_id=prog1_id,
                anchor2=anchor2,
                program2_id=prog2_id,
                bivariate_morans_i=0.0,
                bivariate_pvalue=1.0,
                pearson_correlation=0.0,
                pearson_pvalue=1.0,
                n_spots_used=n_active,
                top_genes_overlap=[],
                relationship_type="independent",
            )
            all_relationships.append(relationship)
            pvalues_for_fdr.append(1.0)
            continue

        # Compute bivariate Moran's I with permutation test
        # Uses the full H vectors (not just active spots) to capture spatial structure
        bivariate_i = _compute_bivariate_morans_i_programs(H1, H2, neighbors)

        # Permutation test for p-value
        null_i = np.zeros(n_permutations)
        for p in range(n_permutations):
            shuffled_h2 = rng.permutation(H2)
            null_i[p] = _compute_bivariate_morans_i_programs(H1, shuffled_h2, neighbors)

        # Two-tailed p-value
        bivariate_pvalue = (np.sum(np.abs(null_i) >= np.abs(bivariate_i)) + 1) / (n_permutations + 1)

        # Compute Pearson correlation (non-spatial baseline)
        if np.std(H1) > 1e-10 and np.std(H2) > 1e-10:
            pearson_r, pearson_p = pearsonr(H1, H2)
        else:
            pearson_r, pearson_p = 0.0, 1.0

        # Compute gene overlap
        program1 = result1.programs[prog1_id]
        program2 = result2.programs[prog2_id]
        gene_overlap = _compute_gene_overlap(program1, program2)

        # Store relationship (classification done after FDR)
        relationship = ProgramPairRelationship(
            anchor1=anchor1,
            program1_id=prog1_id,
            anchor2=anchor2,
            program2_id=prog2_id,
            bivariate_morans_i=float(bivariate_i),
            bivariate_pvalue=float(bivariate_pvalue),
            pearson_correlation=float(pearson_r),
            pearson_pvalue=float(pearson_p),
            n_spots_used=n_active,
            top_genes_overlap=gene_overlap,
            relationship_type="pending",  # Set after FDR
        )
        all_relationships.append(relationship)
        pvalues_for_fdr.append(bivariate_pvalue)

    # Apply FDR correction (Benjamini-Hochberg)
    pvalues_array = np.array(pvalues_for_fdr)
    fdr_significant = false_discovery_control(pvalues_array, method="bh") < fdr_threshold

    # Classify relationships based on FDR-corrected significance
    significant_pairs: List[ProgramPairRelationship] = []
    for i, relationship in enumerate(all_relationships):
        is_significant = fdr_significant[i]
        relationship_type = _classify_relationship(
            relationship.bivariate_morans_i,
            relationship.bivariate_pvalue,
            is_significant,
            threshold=min_bivariate_i,
        )
        # Update relationship type (dataclass is mutable)
        all_relationships[i] = ProgramPairRelationship(
            anchor1=relationship.anchor1,
            program1_id=relationship.program1_id,
            anchor2=relationship.anchor2,
            program2_id=relationship.program2_id,
            bivariate_morans_i=relationship.bivariate_morans_i,
            bivariate_pvalue=relationship.bivariate_pvalue,
            pearson_correlation=relationship.pearson_correlation,
            pearson_pvalue=relationship.pearson_pvalue,
            n_spots_used=relationship.n_spots_used,
            top_genes_overlap=relationship.top_genes_overlap,
            relationship_type=relationship_type,
        )

        if is_significant and abs(relationship.bivariate_morans_i) >= min_bivariate_i:
            significant_pairs.append(all_relationships[i])

    # Sort significant pairs by |bivariate Moran's I|
    significant_pairs.sort(key=lambda x: abs(x.bivariate_morans_i), reverse=True)

    n_significant = len(significant_pairs)
    logger.info("Found %s significant pairs (FDR < %s, |I| >= %s)", n_significant, fdr_threshold, min_bivariate_i)

    # Build result
    bivariate_result = BivariateProgramResult(
        significant_pairs=significant_pairs,
        all_pairs=all_relationships,
        n_programs_total=n_programs,
        n_pairs_tested=n_pairs,
        n_significant=n_significant,
        fdr_threshold=fdr_threshold,
        parameters={
            "min_bivariate_i": min_bivariate_i,
            "n_permutations": n_permutations,
            "neighbor_k": neighbor_k,
            "include_within_anchor": include_within_anchor,
            "random_state": random_state,
        },
    )

    return bivariate_result


def _compute_bivariate_morans_i_programs(
    H1: NDArray[np.floating],
    H2: NDArray[np.floating],
    neighbors: List[List[int]],
) -> float:
    """
    Compute bivariate Moran's I between two program activity vectors.

    This measures spatial cross-correlation: when program 1 is active in a spot,
    are neighboring spots also active for program 2?

    Formula: I_AB = mean(a_centered * spatial_lag(b_centered))

    Args:
        H1: Program 1 activities (n_spots,).
        H2: Program 2 activities (n_spots,).
        neighbors: List of neighbor indices per spot.

    Returns:
        Bivariate Moran's I in range approximately [-1, 1].
    """
    n_spots = len(H1)

    # Center and standardize
    mean1 = np.mean(H1)
    mean2 = np.mean(H2)
    std1 = np.std(H1)
    std2 = np.std(H2)

    if std1 < 1e-10 or std2 < 1e-10:
        return 0.0

    H1_centered = (H1 - mean1) / std1
    H2_centered = (H2 - mean2) / std2

    # Compute spatially-lagged H2 (mean of H2 in neighborhood of each spot)
    H2_lagged = np.zeros(n_spots)
    for i in range(n_spots):
        if len(neighbors[i]) > 0:
            H2_lagged[i] = np.mean(H2_centered[neighbors[i]])

    # Bivariate Moran's I = correlation between H1 and spatially-lagged H2
    bivariate_i = np.mean(H1_centered * H2_lagged)

    return float(bivariate_i)


# =============================================================================
# Helper Functions for Deconvolved Layers
# =============================================================================


def stack_deconvolved_layers(
    adata: sc.AnnData,
    layer_pattern: str = "_genes_pass1",
    coord_key: str = "spatial",
    cell_type_names: Optional[List[str]] = None,
) -> sc.AnnData:
    """
    Stack cell-type-specific deconvolved gene expression layers into a single AnnData.

    Takes the deconvolved layers from Module 3 (one layer per cell type) and stacks
    them vertically so each (spot, cell_type) combination becomes a row. This allows
    NMF to discover programs from already-separated cell-type-specific expression.

    Example:
        Input layers:
            adata.layers['Profile_1_genes_pass1']: (N_spots, M_genes)
            adata.layers['Profile_2_genes_pass1']: (N_spots, M_genes)

        Output stacked AnnData:
            X: (N_spots * T_celltypes, M_genes)
            obs['original_spot']: spot barcode
            obs['cell_type']: which cell type this row represents
            obsm['spatial']: duplicated (x,y) coordinates

    Args:
        adata: AnnData with deconvolved layers from Module 3.
        layer_pattern: Pattern to identify deconvolved layers (e.g., "_genes_pass1").
        coord_key: Key in obsm for spatial coordinates.
        cell_type_names: Specific cell types to include. If None, auto-detect from layers.

    Returns:
        Stacked AnnData with shape (N_spots * T_celltypes, M_genes).
            - obs['original_spot']: Original spot barcode
            - obs['cell_type']: Cell type for this row
            - obs['original_spot_idx']: Integer index into original adata
            - obsm['spatial']: Spatial coordinates (duplicated per cell type)
    """
    # Find deconvolved layers
    if cell_type_names is None:
        deconv_layers = [layer_name for layer_name in adata.layers.keys() if layer_pattern in layer_name]
        if not deconv_layers:
            raise ValueError(
                f"No deconvolved layers found matching pattern '{layer_pattern}'. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        # Extract cell type names from layer names
        # e.g., 'Profile_1_genes_pass1' -> 'Profile_1'
        cell_type_names = [layer.replace(layer_pattern, "") for layer in deconv_layers]
    else:
        deconv_layers = [f"{ct}{layer_pattern}" for ct in cell_type_names]
        # Verify all layers exist
        missing = [lyr for lyr in deconv_layers if lyr not in adata.layers]
        if missing:
            raise ValueError(f"Missing layers: {missing}")

    n_spots = adata.shape[0]
    n_genes = adata.shape[1]
    n_celltypes = len(cell_type_names)

    logger.info("Stacking %s deconvolved layers: %s", n_celltypes, cell_type_names)
    logger.info("  Input: %s spots x %s genes", n_spots, n_genes)
    logger.info("  Output: %s rows x %s genes", n_spots * n_celltypes, n_genes)

    # Stack expression matrices
    stacked_X = []
    stacked_obs = []
    stacked_coords = []

    coords = adata.obsm.get(coord_key, None)
    if coords is None:
        logger.warning("No spatial coordinates found at obsm['%s']", coord_key)
        coords = np.zeros((n_spots, 2))

    for ct_idx, (cell_type, layer_name) in enumerate(zip(cell_type_names, deconv_layers)):
        # Get layer data
        layer_data = adata.layers[layer_name]
        if scipy.sparse.issparse(layer_data):
            layer_data = layer_data.toarray()

        stacked_X.append(layer_data)

        # Build obs for this cell type
        for spot_idx, spot_name in enumerate(adata.obs_names):
            stacked_obs.append(
                {
                    "original_spot": spot_name,
                    "cell_type": cell_type,
                    "original_spot_idx": spot_idx,
                    "cell_type_idx": ct_idx,
                }
            )

        # Duplicate coordinates
        stacked_coords.append(coords.copy())

    # Concatenate
    stacked_X = np.vstack(stacked_X)  # type: ignore[assignment]
    stacked_coords = np.vstack(stacked_coords)  # type: ignore[assignment]
    stacked_obs_df = pd.DataFrame(stacked_obs)

    # Create new row names: {spot}_{celltype}
    row_names = [f"{row['original_spot']}_{row['cell_type']}" for _, row in stacked_obs_df.iterrows()]
    stacked_obs_df.index = pd.Index(row_names)

    # Create stacked AnnData
    stacked_adata = sc.AnnData(
        X=stacked_X,
        obs=stacked_obs_df,
        var=adata.var.copy(),
    )
    stacked_adata.obsm[coord_key] = stacked_coords

    # Copy over any useful uns metadata
    if "anchored_programs" in adata.uns:
        stacked_adata.uns["anchored_programs"] = adata.uns["anchored_programs"]

    logger.info("  Created stacked AnnData: %s", stacked_adata.shape)

    return stacked_adata


def unstack_program_results(
    stacked_H: NDArray[np.floating],
    stacked_adata: sc.AnnData,
    n_spots: int,
) -> Dict[str, NDArray[np.floating]]:
    """
    Unstack program activities back to per-cell-type matrices.

    Takes the H matrix from NMF on stacked data and separates it back into
    per-cell-type program activity matrices.

    Args:
        stacked_H: Program activities from NMF (K_programs, N_spots * T_celltypes).
        stacked_adata: The stacked AnnData used for NMF.
        n_spots: Number of original spots.

    Returns:
        Dict mapping cell_type -> H matrix (K_programs, n_spots)
    """
    cell_types = stacked_adata.obs["cell_type"].unique()

    # H is (K, N*T), need to split into T matrices of (K, N)
    H_by_celltype = {}

    for ct_idx, cell_type in enumerate(cell_types):
        # Rows for this cell type are at indices [ct_idx*n_spots : (ct_idx+1)*n_spots]
        start_idx = ct_idx * n_spots
        end_idx = (ct_idx + 1) * n_spots
        H_by_celltype[cell_type] = stacked_H[:, start_idx:end_idx]

    return H_by_celltype


def extract_celltype_expression(
    adata: sc.AnnData,
    cell_type: str,
    layer_pattern: str = "_genes_pass1",
) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
    """
    Extract deconvolved expression for a single cell type from Module 3 layers.

    Args:
        adata: AnnData with deconvolved layers.
        cell_type: Name of cell type (e.g., 'Profile_1').
        layer_pattern: Pattern for layer names.

    Returns:
        Tuple of (expression matrix, spatial coordinates)
        Expression shape: (n_spots, n_genes)
    """
    layer_name = f"{cell_type}{layer_pattern}"

    if layer_name not in adata.layers:
        raise ValueError(f"Layer '{layer_name}' not found. Available: {list(adata.layers.keys())}")

    X = adata.layers[layer_name]
    if scipy.sparse.issparse(X):
        X = X.toarray()

    coords = adata.obsm.get("spatial", np.zeros((adata.shape[0], 2)))

    return np.array(X), np.array(coords)


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
        "Cell type '%s': %s/%s spots have proportion >= %s",
        cell_type_name,
        mask.sum(),
        len(mask),
        min_proportion_threshold,
    )

    return weights, mask


def _assign_program_cell_types(
    H: NDArray[np.floating],
    cell_type_proportions: pd.DataFrame,
    spot_names: Optional[List[str]] = None,
    single_threshold: float = 0.7,
    interaction_threshold: float = 0.25,
) -> List[Dict[str, Any]]:
    """
    Assign cell type labels to joint programs based on correlation with proportions.

    Args:
        H: Program activities (K_programs, n_spots).
        cell_type_proportions: Cell type proportions per spot (n_spots, n_celltypes).
        spot_names: Names of spots in H (used for alignment with proportions).
            If None, assumes H and proportions have same order.
        single_threshold: Min enrichment for single cell-type classification.
        interaction_threshold: Min enrichment for secondary cell type in interaction.

    Returns:
        List of dicts with cell type assignment info per program.
    """
    K_programs = H.shape[0]
    cell_types = list(cell_type_proportions.columns)

    # Handle spot alignment between H and proportions
    if spot_names is not None:
        # Find common spots between H and proportions
        h_spots = pd.Index(spot_names)
        prop_spots = cell_type_proportions.index
        common_spots = h_spots.intersection(prop_spots)

        if len(common_spots) == 0:
            logger.warning("No common spots between H and proportions. Using all spots.")
            h_indices = np.arange(H.shape[1])
            prop_aligned = cell_type_proportions
        else:
            logger.info(
                "Aligned %s spots (H has %s, proportions has %s)",
                len(common_spots),
                len(h_spots),
                len(prop_spots),
            )
            # Get indices in H for common spots
            h_indices = np.array([h_spots.get_loc(s) for s in common_spots])
            # Align proportions to same order
            prop_aligned = cell_type_proportions.loc[common_spots]
    else:
        # Assume same order
        h_indices = np.arange(H.shape[1])
        prop_aligned = cell_type_proportions

    results = []

    for k in range(K_programs):
        # Use only aligned spots
        h_k = H[k, h_indices]

        # Compute correlation with each cell type's proportions
        enrichments = {}
        for ct in cell_types:
            if ct == "Unknown":
                continue
            props = prop_aligned[ct].values
            # Use Spearman correlation (rank-based, more robust)
            if np.std(h_k) > 1e-10 and np.std(props) > 1e-10:
                corr, _ = spearmanr(h_k, props)
                enrichments[ct] = max(0, corr)  # Only positive correlations
            else:
                enrichments[ct] = 0.0

        # Normalize to sum to 1
        total = sum(enrichments.values())
        if total > 0:
            enrichments = {ct: v / total for ct, v in enrichments.items()}

        # Sort by enrichment
        sorted_cts = sorted(enrichments.items(), key=lambda x: x[1], reverse=True)

        primary_ct = sorted_cts[0][0] if sorted_cts else "Unknown"
        primary_score = sorted_cts[0][1] if sorted_cts else 0.0

        secondary_ct = sorted_cts[1][0] if len(sorted_cts) > 1 else None
        secondary_score = sorted_cts[1][1] if len(sorted_cts) > 1 else 0.0

        # Compute interaction score: how balanced is the distribution?
        # 0 = all in one cell type, 1 = perfectly balanced
        if len(enrichments) > 1:
            max_enrich = max(enrichments.values())
            # Gini-like: 1 - (max / ideal_balanced)
            ideal_balanced = 1.0 / len(enrichments)
            interaction_score = 1.0 - (max_enrich - ideal_balanced) / (1.0 - ideal_balanced)
            interaction_score = max(0, min(1, interaction_score))
        else:
            interaction_score = 0.0

        # Classify program type
        if primary_score >= single_threshold:
            program_type = "single_celltype"
        elif secondary_score >= interaction_threshold:
            program_type = "interaction"
        else:
            program_type = "microenvironment"

        results.append(
            {
                "cell_type_enrichments": enrichments,
                "primary_cell_type": primary_ct,
                "secondary_cell_type": secondary_ct if secondary_score >= interaction_threshold else None,
                "interaction_score": interaction_score,
                "program_type": program_type,
            }
        )

    return results


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
    from scipy.spatial import cKDTree  # pylint: disable=import-outside-toplevel

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
    denominator = np.sum(z**2)

    if denominator == 0:
        return 0.0, 1.0

    I = (n / W.sum()) * (numerator / denominator)  # noqa: E741 — Moran's I

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
    *,
    K: int = 5,
    lambda_spatial: float = 0.0,
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
        init="nndsvda",
        max_iter=max_iter,
        random_state=random_state,
        solver="cd",  # Frobenius norm, coordinate descent
        loss="frobenius",
        alpha_W=lambda_sparsity,
        l1_ratio=0.5,  # Balance L1 and L2 for W
    )

    try:
        H = nmf.fit_transform(X_weighted)  # (n_spots, K)
        W = nmf.components_.T  # (n_genes, K) after transpose
    except ValueError as e:
        logger.warning("NMF failed: %s. Using random initialization.", e)
        W = np.random.rand(n_genes, K)
        H = np.random.rand(n_spots, K)

    # Apply spatial smoothing to H using Laplacian
    if lambda_spatial > 0 and n_spots >= 10:
        L = build_spatial_laplacian(coords, k=8, normed=True)

        # Smooth each program's loadings
        # H_smooth = (I + lambda * L)^(-1) @ H
        I = scipy.sparse.eye(n_spots)  # noqa: E741 — identity matrix
        smoothing_matrix = I + lambda_spatial * L

        try:
            from scipy.sparse.linalg import spsolve  # pylint: disable=import-outside-toplevel

            H_smooth = np.zeros_like(H)
            for k in range(K):
                H_smooth[:, k] = spsolve(smoothing_matrix.tocsc(), H[:, k])
            H = np.maximum(H_smooth, 0)  # Keep non-negative
        except (ImportError, OSError) as e:
            logger.debug("Spatial smoothing failed: %s. Using unsmoothed H.", e)

    # Compute reconstruction error
    X_reconstructed = H @ W.T
    reconstruction_error = np.linalg.norm(X_weighted - X_reconstructed * sqrt_weights, "fro")

    return W, H.T, float(reconstruction_error)  # H transposed to (K, n_spots)


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
        for p, _ in enumerate(protein_names):
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
    df.index.name = "program"

    # Add validation column: does program correlate with anchor proteins?
    anchor_mask = [p in anchor_proteins for p in protein_names]
    if any(anchor_mask):
        df["anchor_correlation"] = df.iloc[:, anchor_mask].max(axis=1)
        df["validated"] = df["anchor_correlation"] > 0.3
    else:
        df["anchor_correlation"] = np.nan
        df["validated"] = False

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
    from sklearn.cluster import KMeans  # pylint: disable=import-outside-toplevel
    from sklearn.preprocessing import StandardScaler  # pylint: disable=import-outside-toplevel

    n_spots = H.shape[1]
    K = H.shape[0]

    if n_spots < n_clusters * min_spots_per_cluster:
        logger.warning("Too few spots (%s) for %s subpopulations", n_spots, n_clusters)
        return []

    # Normalize program loadings and coordinates separately
    H_norm = StandardScaler().fit_transform(H.T)  # (n_spots, K)
    coords_norm = StandardScaler().fit_transform(coords)  # (n_spots, 2)

    # Combine features with spatial weighting
    # Higher spatial_weight = more emphasis on spatial location
    features = np.hstack([(1 - spatial_weight) * H_norm, spatial_weight * coords_norm])

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
        dist_from_center = np.sqrt((centroid[0] - all_centroid[0]) ** 2 + (centroid[1] - all_centroid[1]) ** 2)
        max_dist = np.sqrt(((coords[:, 0] - all_centroid[0]) ** 2 + (coords[:, 1] - all_centroid[1]) ** 2).max())

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
    *,
    profile_discovery_result: Optional[ProfileDiscoveryResult] = None,
    protein_adata=None,
    K_programs: int = 5,
    lambda_spatial: float = 0.0,
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
    logger.info("Starting CONTRASTIVE anchor-specific program discovery (Module 4)")
    logger.info("  contrastive_strength=%s, use_enriched_genes=%s", contrastive_strength, use_enriched_genes)

    # Get gene expression data
    if scipy.sparse.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)

    gene_names = list(adata.var_names)
    coords = adata.obsm["spatial"]
    n_spots, n_genes = X.shape

    logger.info("Input: %s spots, %s genes, %s cell types", n_spots, n_genes, len(cell_type_proportions.columns))

    # Get protein data if available
    protein_data = None
    protein_names = None
    if validate_with_proteins and protein_adata is not None:
        if scipy.sparse.issparse(protein_adata.X):
            protein_data = protein_adata.X.toarray()
        else:
            protein_data = np.array(protein_adata.X)
        protein_names = list(protein_adata.var_names)
        logger.info("Protein data: %s proteins for validation", protein_data.shape[1])

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

        logger.info("Processing anchor: %s", cell_type)

        # Get weights from Module 3 proportions
        weights, mask = _get_celltype_weights(cell_type_proportions, cell_type, min_proportion_threshold)

        n_included = mask.sum()
        if n_included < 20:
            logger.warning(
                "Skipping %s: only %s spots have proportion >= %s", cell_type, n_included, min_proportion_threshold
            )
            continue

        # Compute background expression from anchor-LOW spots (contrastive)
        background = _compute_background_expression(
            X,
            cell_type_proportions,
            cell_type,
            low_threshold=min_proportion_threshold / 2,  # Use stricter threshold for background
        )

        # Optionally select anchor-enriched genes
        if use_enriched_genes:
            enriched_idx, enriched_names = _select_anchor_enriched_genes(
                X,
                weights,
                gene_names,
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
                if cell_type.lower() in pname.lower() or any(cell_type.lower() in p.lower() for p in prots):
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
            protein_correlations = validate_programs_with_proteins(H, protein_subset, protein_names, anchor_proteins)

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
            program_var = np.var(H[k, :]) * np.sum(loadings**2)
            var_explained = program_var / total_var if total_var > 0 else 0

            # Compute spatial Moran's I for this program
            moran_i, moran_p = _compute_spatial_moran_i(H[k, :], coords_subset, k=8, n_permutations=99)

            # Program activity statistics
            mean_activity = float(np.mean(H[k, :]))
            median_activity = float(np.median(H[k, :]))
            active_fraction = float(np.mean(H[k, :] > median_activity))

            programs.append(
                SpatialProgram(
                    program_id=k,
                    top_genes=top_genes,
                    gene_loadings=gene_loadings,
                    variance_explained=var_explained,
                    spatial_moran_i=moran_i,
                    spatial_moran_pvalue=moran_p,
                    mean_activity=mean_activity,
                    active_spots_fraction=active_fraction,
                )
            )

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
                "K_programs": K_programs,
                "lambda_spatial": lambda_spatial,
                "contrastive_strength": contrastive_strength,
                "use_enriched_genes": use_enriched_genes,
                "n_enriched_genes": len(enriched_idx) if use_enriched_genes else n_genes,
            },
            subpopulations=subpopulations,
        )

        total_programs += K_programs
        logger.info(
            "  %s: %s programs, %s spots, %s genes",
            cell_type,
            K_programs,
            n_included,
            len(enriched_idx) if use_enriched_genes else n_genes,
        )

    # Build final result
    result = AnchoredProgramDiscoveryResult(
        results_by_anchor=results_by_anchor,
        n_anchors=len(results_by_anchor),
        total_programs=total_programs,
        profile_discovery_result=profile_discovery_result,
        parameters={
            "K_programs": K_programs,
            "lambda_spatial": lambda_spatial,
            "contrastive_strength": contrastive_strength,
            "use_enriched_genes": use_enriched_genes,
            "min_proportion_threshold": min_proportion_threshold,
            "random_state": random_state,
        },
    )

    logger.info("Completed: %s anchors, %s programs (contrastive)", result.n_anchors, result.total_programs)

    return result


def discover_joint_programs(
    adata: sc.AnnData,
    cell_type_proportions: pd.DataFrame,
    *,
    K_programs: int = 10,
    layer_pattern: str = "_genes_pass1",
    lambda_spatial: float = 0.0,
    lambda_sparsity: float = 0.01,
    top_n_genes: int = 50,
    random_state: int = 42,
) -> JointDiscoveryResult:
    """
    Discover spatial programs jointly across all cell types.

    Unlike discover_anchored_programs which runs NMF per cell type, this function
    stacks all deconvolved layers and runs a single NMF to find programs that
    may span multiple cell types (e.g., tumor-immune interface programs).

    Args:
        adata: AnnData with deconvolved layers from Module 3.
        cell_type_proportions: Module 3 output - cell type proportions per spot.
        K_programs: Number of programs to discover.
        layer_pattern: Pattern to identify deconvolved layers.
        lambda_spatial: Spatial smoothness regularization (not yet implemented).
        lambda_sparsity: Sparsity regularization via NMF alpha.
        top_n_genes: Number of top genes to report per program.
        random_state: Random seed.

    Returns:
        JointDiscoveryResult with programs and cell type assignments.
    """
    logger.info("Starting JOINT program discovery across all cell types")

    # Stack deconvolved layers
    stacked_adata = stack_deconvolved_layers(
        adata,
        layer_pattern=layer_pattern,
    )

    # Get stacked expression matrix
    if scipy.sparse.issparse(stacked_adata.X):
        X_stacked = stacked_adata.X.toarray()
    else:
        X_stacked = np.array(stacked_adata.X)

    # Ensure non-negative
    X_stacked = np.maximum(X_stacked, 0)

    gene_names = list(stacked_adata.var_names)
    cell_type_names = list(stacked_adata.obs["cell_type"].unique())
    n_spots = adata.shape[0]

    logger.info("Stacked data: %s rows (%s spots x %s cell types)", X_stacked.shape[0], n_spots, len(cell_type_names))
    logger.info("Discovering %s joint programs", K_programs)

    # Run NMF on stacked data
    # X_stacked is (n_spots * n_celltypes, n_genes)
    # NMF: X ≈ W @ H where W is (n_samples, K), H is (K, n_features)
    nmf = NMF(
        n_components=K_programs,
        init="nndsvda",
        random_state=random_state,
        max_iter=500,
        alpha_W=lambda_sparsity,
        alpha_H=0.0,
        l1_ratio=0.5,
    )

    # W_stacked: (n_spots * n_celltypes, K) - activities per stacked row
    # H_nmf: (K, n_genes) - gene loadings
    W_stacked = nmf.fit_transform(X_stacked)
    H_nmf = nmf.components_

    # We want:
    # W = gene loadings (n_genes, K)
    # H = spot activities (K, n_spots)
    W = H_nmf.T  # (n_genes, K)

    # Unstack W_stacked to get per-cell-type activities
    # W_stacked is (n_spots * n_celltypes, K)
    H_by_celltype = {}
    for ct_idx, cell_type in enumerate(cell_type_names):
        start_idx = ct_idx * n_spots
        end_idx = (ct_idx + 1) * n_spots
        H_by_celltype[cell_type] = W_stacked[start_idx:end_idx, :].T  # (K, n_spots)

    # Average across cell types for overall spot activity
    H = np.zeros((K_programs, n_spots))
    for ct_H in H_by_celltype.values():
        H += ct_H
    H /= len(H_by_celltype)

    # Compute reconstruction error
    X_reconstructed = W_stacked @ H_nmf
    recon_error = np.mean((X_stacked - X_reconstructed) ** 2)

    # Assign cell types to programs
    # Pass spot names to handle alignment between H (from adata) and proportions
    spot_names = list(adata.obs_names)
    cell_type_assignments = _assign_program_cell_types(H, cell_type_proportions, spot_names=spot_names)

    # Get spatial coordinates
    coords = adata.obsm.get("spatial", np.zeros((n_spots, 2)))

    # Build JointProgram objects
    programs = []
    total_var = np.var(X_stacked)

    for k in range(K_programs):
        # Top genes
        loadings = W[:, k]
        top_indices = np.argsort(loadings)[::-1][:top_n_genes]
        top_genes = [gene_names[i] for i in top_indices]
        gene_loadings = {gene_names[i]: float(loadings[i]) for i in top_indices}

        # Variance explained
        program_var = np.var(H[k, :]) * np.sum(loadings**2)
        var_explained = program_var / total_var if total_var > 0 else 0

        # Spatial Moran's I
        moran_i, moran_p = _compute_spatial_moran_i(H[k, :], coords, k=8, n_permutations=99)

        # Activity stats
        mean_activity = float(np.mean(H[k, :]))
        median_activity = float(np.median(H[k, :]))
        active_fraction = float(np.mean(H[k, :] > median_activity))

        # Cell type assignment
        ct_info = cell_type_assignments[k]

        programs.append(
            JointProgram(
                program_id=k,
                top_genes=top_genes,
                gene_loadings=gene_loadings,
                variance_explained=var_explained,
                spatial_moran_i=moran_i,
                spatial_moran_pvalue=moran_p,
                mean_activity=mean_activity,
                active_spots_fraction=active_fraction,
                cell_type_enrichments=ct_info["cell_type_enrichments"],
                primary_cell_type=ct_info["primary_cell_type"],
                secondary_cell_type=ct_info["secondary_cell_type"],
                interaction_score=ct_info["interaction_score"],
                program_type=ct_info["program_type"],
            )
        )

    logger.info("Discovered %s joint programs", len(programs))
    logger.info("  Single cell-type: %s", sum(1 for p in programs if p.program_type == "single_celltype"))
    logger.info("  Interaction: %s", sum(1 for p in programs if p.program_type == "interaction"))
    logger.info("  Microenvironment: %s", sum(1 for p in programs if p.program_type == "microenvironment"))

    return JointDiscoveryResult(
        programs=programs,
        W=W,
        H=H,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        n_spots=n_spots,
        reconstruction_error=recon_error,
        H_by_celltype=H_by_celltype,
        parameters={
            "K_programs": K_programs,
            "layer_pattern": layer_pattern,
            "lambda_spatial": lambda_spatial,
            "lambda_sparsity": lambda_sparsity,
            "random_state": random_state,
        },
    )


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
        key = f"X_anchored_programs_{anchor_name}"
        adata.obsm[key] = anchor_result.H.T  # (n_spots, K)

        # Store gene loadings in varm if genes match
        if list(adata.var_names) == anchor_result.gene_names:
            varm_key = f"anchored_program_loadings_{anchor_name}"
            adata.varm[varm_key] = anchor_result.W  # (n_genes, K)

    # Store metadata in uns
    adata.uns["anchored_programs"] = {
        "n_anchors": result.n_anchors,
        "total_programs": result.total_programs,
        "parameters": result.parameters,
        "anchors": list(result.results_by_anchor.keys()),
        "summary": result.summary(),
    }

    logger.info("Stored results in adata: %s anchors in obsm/varm/uns", result.n_anchors)


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
    n_spots, _ = X.shape

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
    *,
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
        logger.warning("Too few spots for enrichment (high=%s, low=%s), using all genes", n_high, n_low)
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

    logger.debug("Selected %s enriched genes (FC >= %s)", len(enriched_indices), min_fold_change)

    return list(enriched_indices), enriched_names


def contrastive_anchored_nmf(
    X: NDArray[np.floating],
    anchor_weights: NDArray[np.floating],
    background: NDArray[np.floating],
    coords: NDArray[np.floating],
    *,
    K: int = 5,
    lambda_spatial: float = 0.0,
    lambda_sparsity: float = 0.01,
    contrastive_strength: float = 0.8,
    max_iter: int = 200,
    random_state: int = 42,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], float]:
    """
    Contrastive NMF for anchor-specific program discovery.

    DEPRECATED: Use discover_programs_from_layers() instead, which uses
    deconvolved expression from Module 3 rather than contrastive subtraction.

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
    n_spots, _ = X.shape

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
        init="nndsvda",
        max_iter=max_iter,
        random_state=random_state,
        alpha_W=lambda_sparsity,
        l1_ratio=0.5,
    )

    try:
        H = nmf.fit_transform(X_weighted)  # (n_spots, K)
        W = nmf.components_.T  # (n_genes, K)
    except ValueError as e:
        logger.warning("Contrastive NMF failed: %s. Falling back to standard NMF.", e)
        return anchored_spatial_nmf(
            X,
            anchor_weights,
            coords,
            K=K,
            lambda_spatial=lambda_spatial,
            lambda_sparsity=lambda_sparsity,
            max_iter=max_iter,
            random_state=random_state,
        )

    # Apply spatial smoothing
    if lambda_spatial > 0 and n_spots >= 10:
        L = build_spatial_laplacian(coords, k=8, normed=True)
        I = scipy.sparse.eye(n_spots)  # noqa: E741 — identity matrix
        smoothing_matrix = I + lambda_spatial * L

        try:
            from scipy.sparse.linalg import spsolve  # pylint: disable=import-outside-toplevel

            H_smooth = np.zeros_like(H)
            for k in range(K):
                H_smooth[:, k] = spsolve(smoothing_matrix.tocsc(), H[:, k])
            H = np.maximum(H_smooth, 0)
        except (ImportError, OSError) as e:
            logger.debug("Spatial smoothing failed: %s", e)

    # Reconstruction error (on original weighted data)
    X_reconstructed = H @ W.T
    reconstruction_error = np.linalg.norm(X_weighted - X_reconstructed * sqrt_weights, "fro")

    return W, H.T, float(reconstruction_error)


# =============================================================================
# MODULE 4 v3: DECONVOLVED LAYER-BASED PROGRAM DISCOVERY
# =============================================================================


def discover_programs_from_layers(
    adata: sc.AnnData,
    *,
    layer_pattern: str = "_genes_pass1",
    cell_type_proportions: Optional[pd.DataFrame] = None,
    profile_discovery_result: Optional[ProfileDiscoveryResult] = None,
    protein_adata: Optional[sc.AnnData] = None,
    K_programs: int = 5,
    lambda_spatial: float = 0.0,
    lambda_sparsity: float = 0.01,
    min_expression_threshold: float = 0.0,
    validate_with_proteins: bool = True,
    top_n_genes: int = 50,
    detect_subpopulations: bool = True,
    n_subpopulations: int = 3,
    random_state: int = 42,
) -> AnchoredProgramDiscoveryResult:
    """
    Discover spatial transcriptomic programs from DECONVOLVED expression layers.

    This is the recommended Module 4 entry point. Uses the cell-type-specific
    deconvolved expression layers from Module 3 rather than raw expression with
    contrastive subtraction.

    Key advantages over contrastive approach:
    1. Uses actual deconvolved expression (not approximated)
    2. No need to estimate/subtract background
    3. More accurate cell-type-specific programs
    4. Leverages Module 3's QP optimization results directly

    Args:
        adata: AnnData with deconvolved layers from Module 3.
            Expected layers: '{CellType}_genes_pass1' for each cell type.
        layer_pattern: Pattern to identify deconvolved layers (default: "_genes_pass1").
        cell_type_proportions: Module 3 output - cell type proportions (optional).
            Used to weight spots by cell type abundance.
        profile_discovery_result: Module 2 output - protein profiles for validation.
        protein_adata: AnnData with protein data for validation (optional).
        K_programs: Number of programs to discover per cell type.
        lambda_spatial: Spatial smoothness regularization weight.
        lambda_sparsity: Gene loading sparsity weight.
        min_expression_threshold: Minimum total expression to include a spot.
        validate_with_proteins: Whether to validate programs against proteins.
        top_n_genes: Number of top genes to report per program.
        detect_subpopulations: Whether to detect spatial subpopulations.
        n_subpopulations: Number of subpopulations to detect per cell type.
        random_state: Random seed for reproducibility.

    Returns:
        AnchoredProgramDiscoveryResult with programs for each cell type.
    """
    logger.info("Starting DECONVOLVED LAYER-BASED program discovery (Module 4 v3)")

    # Find deconvolved layers
    deconv_layers = [layer_name for layer_name in adata.layers.keys() if layer_pattern in layer_name]
    if not deconv_layers:
        raise ValueError(
            f"No deconvolved layers found matching pattern '{layer_pattern}'. "
            f"Available layers: {list(adata.layers.keys())}. "
            f"Run Module 3 (optimize_gene_expression) first."
        )

    # Extract cell type names from layer names
    cell_type_names = [layer.replace(layer_pattern, "") for layer in deconv_layers]

    logger.info("Found %s deconvolved layers: %s", len(cell_type_names), cell_type_names)

    # Get spatial coordinates
    coords = adata.obsm.get("spatial", None)
    if coords is None:
        raise ValueError("No spatial coordinates found in adata.obsm['spatial']")

    gene_names = list(adata.var_names)
    n_spots = adata.shape[0]
    n_genes = adata.shape[1]

    logger.info("Input: %s spots, %s genes", n_spots, n_genes)

    # Get protein data if available
    protein_data = None
    protein_names = None
    if validate_with_proteins and protein_adata is not None:
        if scipy.sparse.issparse(protein_adata.X):
            protein_data = protein_adata.X.toarray()
        else:
            protein_data = np.array(protein_adata.X)
        protein_names = list(protein_adata.var_names)
        logger.info("Protein data: %s proteins for validation", protein_data.shape[1])

    # Build profile name -> proteins mapping from Module 2
    profile_to_proteins = {}
    if profile_discovery_result is not None:
        for i, profile in enumerate(profile_discovery_result.profiles):
            profile_name = f"Profile_{i}"
            profile_to_proteins[profile_name] = list(profile)

    # Discover programs for each cell type
    results_by_anchor = {}
    total_programs = 0

    for cell_type in cell_type_names:
        logger.info("Processing cell type: %s", cell_type)

        # Extract deconvolved expression for this cell type
        X_celltype, _ = extract_celltype_expression(adata, cell_type, layer_pattern)

        # Get cell type proportions for weighting (if available)
        if cell_type_proportions is not None and cell_type in cell_type_proportions.columns:
            weights = cell_type_proportions[cell_type].values
        else:
            # Equal weighting if no proportions available
            weights = np.ones(n_spots)

        # Filter spots with sufficient expression
        spot_totals = np.sum(X_celltype, axis=1)
        active_mask = spot_totals > min_expression_threshold
        n_active = active_mask.sum()

        if n_active < 20:
            logger.warning(
                "Skipping %s: only %s spots have expression > %s", cell_type, n_active, min_expression_threshold
            )
            continue

        X_subset = X_celltype[active_mask, :]
        coords_subset = coords[active_mask, :]
        weights_subset = weights[active_mask]

        logger.info("  %s active spots (expression > %s)", n_active, min_expression_threshold)

        # Determine anchor proteins for validation
        anchor_proteins = profile_to_proteins.get(cell_type, [])

        # Run spatially-regularized NMF on deconvolved expression
        W, H, recon_error = _deconvolved_spatial_nmf(
            X_subset,
            weights_subset,
            coords_subset,
            K=K_programs,
            lambda_spatial=lambda_spatial,
            lambda_sparsity=lambda_sparsity,
            random_state=random_state,
        )

        # Build full H matrix (zeros for filtered spots)
        H_full = np.zeros((K_programs, n_spots))
        H_full[:, active_mask] = H

        # Validate with proteins if available
        protein_correlations = pd.DataFrame()
        if validate_with_proteins and protein_data is not None and protein_names is not None:
            protein_subset = protein_data[active_mask, :]
            protein_correlations = validate_programs_with_proteins(H, protein_subset, protein_names, anchor_proteins)

        # Build SpatialProgram objects
        programs = []

        # Compute total variance for normalization (proper calculation)
        X_reconstructed = H.T @ W.T  # (n_spots_subset, n_genes)
        total_ss = np.sum((X_subset - X_subset.mean()) ** 2)
        residual_ss = np.sum((X_subset - X_reconstructed) ** 2)
        total_var_explained = 1 - (residual_ss / total_ss) if total_ss > 0 else 0

        for k in range(K_programs):
            # Get top genes
            loadings = W[:, k]
            top_indices = np.argsort(loadings)[::-1][:top_n_genes]
            top_genes = [gene_names[i] for i in top_indices]
            gene_loadings = {gene_names[i]: float(loadings[i]) for i in top_indices}

            # Compute per-program variance explained (approximate)
            # Using the fraction of total W norm
            W_k_norm = np.sum(loadings**2)
            total_W_norm = np.sum(W**2)
            var_explained = (W_k_norm / total_W_norm) * total_var_explained if total_W_norm > 0 else 0

            # Compute spatial Moran's I for this program
            moran_i, moran_p = _compute_spatial_moran_i(H[k, :], coords_subset, k=8, n_permutations=199)

            # Program activity statistics
            mean_activity = float(np.mean(H[k, :]))
            median_activity = float(np.median(H[k, :]))
            active_fraction = float(np.mean(H[k, :] > median_activity))

            programs.append(
                SpatialProgram(
                    program_id=k,
                    top_genes=top_genes,
                    gene_loadings=gene_loadings,
                    variance_explained=var_explained * 100,  # As percentage
                    spatial_moran_i=moran_i,
                    spatial_moran_pvalue=moran_p,
                    mean_activity=mean_activity,
                    active_spots_fraction=active_fraction,
                )
            )

        # Detect spatial subpopulations
        subpopulations = []
        if detect_subpopulations and n_active >= n_subpopulations * 10:
            subpopulations = detect_spatial_subpopulations(
                H=H,
                coords=coords_subset,
                n_clusters=n_subpopulations,
                spatial_weight=0.3,
                min_spots_per_cluster=10,
            )
            if subpopulations:
                logger.info("  Subpopulations: %s detected", len(subpopulations))

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
            n_spots_used=n_active,
            parameters={
                "K_programs": K_programs,
                "lambda_spatial": lambda_spatial,
                "lambda_sparsity": lambda_sparsity,
                "min_expression_threshold": min_expression_threshold,
                "source": "deconvolved_layers",
            },
            subpopulations=subpopulations,
        )

        total_programs += K_programs
        logger.info("  %s: %s programs, %s spots", cell_type, K_programs, n_active)

    # Build final result
    result = AnchoredProgramDiscoveryResult(
        results_by_anchor=results_by_anchor,
        n_anchors=len(results_by_anchor),
        total_programs=total_programs,
        profile_discovery_result=profile_discovery_result,
        parameters={
            "K_programs": K_programs,
            "lambda_spatial": lambda_spatial,
            "lambda_sparsity": lambda_sparsity,
            "min_expression_threshold": min_expression_threshold,
            "layer_pattern": layer_pattern,
            "random_state": random_state,
            "source": "deconvolved_layers",
        },
    )

    logger.info("Completed: %s cell types, %s programs", result.n_anchors, result.total_programs)

    return result


def _deconvolved_spatial_nmf(
    X: NDArray[np.floating],
    weights: NDArray[np.floating],
    coords: NDArray[np.floating],
    *,
    K: int = 5,
    lambda_spatial: float = 0.0,
    lambda_sparsity: float = 0.01,
    max_iter: int = 200,
    random_state: int = 42,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], float]:
    """
    NMF on deconvolved cell-type-specific expression with spatial regularization.

    Unlike contrastive_anchored_nmf, this operates directly on the deconvolved
    expression from Module 3 - no background subtraction needed.

    Args:
        X: Deconvolved gene expression matrix (n_spots, n_genes).
        weights: Cell type proportions for weighting (n_spots,).
        coords: Spatial coordinates (n_spots, 2).
        K: Number of programs to discover.
        lambda_spatial: Spatial smoothness weight.
        lambda_sparsity: Gene sparsity weight.
        max_iter: Maximum NMF iterations.
        random_state: Random seed.

    Returns:
        Tuple of (W gene loadings, H spot loadings, reconstruction error).
    """
    n_spots, n_genes = X.shape

    # Weight expression by cell type proportion
    sqrt_weights = np.sqrt(weights + 1e-10).reshape(-1, 1)
    X_weighted = X * sqrt_weights

    # Ensure non-negative
    X_weighted = np.maximum(X_weighted, 0)

    # Run NMF
    nmf = NMF(
        n_components=K,
        init="nndsvda",
        max_iter=max_iter,
        random_state=random_state,
        alpha_W=lambda_sparsity,
        l1_ratio=0.5,
    )

    try:
        H = nmf.fit_transform(X_weighted)  # (n_spots, K)
        W = nmf.components_.T  # (n_genes, K)
    except ValueError as e:
        logger.warning("NMF failed: %s. Using random initialization.", e)
        W = np.abs(np.random.randn(n_genes, K))
        H = np.abs(np.random.randn(n_spots, K))

    # Apply spatial smoothing to H
    if lambda_spatial > 0 and n_spots >= 10:
        L = build_spatial_laplacian(coords, k=8, normed=True)
        I = scipy.sparse.eye(n_spots)  # noqa: E741 — identity matrix
        smoothing_matrix = I + lambda_spatial * L

        try:
            from scipy.sparse.linalg import spsolve  # pylint: disable=import-outside-toplevel

            H_smooth = np.zeros_like(H)
            for k in range(K):
                H_smooth[:, k] = spsolve(smoothing_matrix.tocsc(), H[:, k])
            H = np.maximum(H_smooth, 0)
        except (ImportError, OSError) as e:
            logger.debug("Spatial smoothing failed: %s", e)

    # Compute reconstruction error
    X_reconstructed = H @ W.T
    reconstruction_error = np.linalg.norm(X_weighted - X_reconstructed * sqrt_weights, "fro")

    return W, H.T, float(reconstruction_error)


# =============================================================================
# MODULE 4c: REGION-AWARE PROGRAM ANALYSIS
# =============================================================================


def analyze_program_regions(
    result: AnchoredProgramResult,
    adata: sc.AnnData,
    region_column: str,
    min_spots_per_region: int = 20,
) -> AnchoredProgramResult:
    """
    Annotate programs with region enrichment statistics.

    For each program, computes activity distribution across regions defined by
    a categorical column in adata.obs, and tests for significant enrichment.

    This is a post-hoc analysis step - programs are discovered first, then
    analyzed for regional variation.

    Args:
        result: AnchoredProgramResult from discover_programs_from_layers()
        adata: AnnData with region annotations in obs
        region_column: Column in adata.obs defining regions (e.g., "D538G Mutation")
        min_spots_per_region: Minimum spots required per region for analysis

    Returns:
        AnchoredProgramResult with region fields populated in each SpatialProgram

    Example:
        >>> result = discover_programs_from_layers(adata, ...)
        >>> result = analyze_program_regions(result, adata, "D538G Mutation")
        >>> for prog in result.programs:
        ...     print(f"Program {prog.program_id}: enriched in {prog.enriched_region}")
    """
    from scipy.stats import kruskal, mannwhitneyu  # pylint: disable=import-outside-toplevel

    if region_column not in adata.obs.columns:
        raise ValueError(f"Region column '{region_column}' not found in adata.obs")

    regions = adata.obs[region_column].unique()
    n_regions = len(regions)

    logger.info("Analyzing program regions using '%s' (%s regions)", region_column, n_regions)

    # Get program activities (H matrix)
    H = result.H  # Shape: (K_programs, n_spots)

    # Analyze each program
    for k, program in enumerate(result.programs):
        h_k = H[k, :]  # Activity for this program across spots

        # Compute mean activity per region
        region_means = {}
        region_activities = {}
        for region in regions:
            mask = (adata.obs[region_column] == region).values
            if mask.sum() >= min_spots_per_region:
                activities = h_k[mask]
                region_means[str(region)] = float(np.mean(activities))
                region_activities[str(region)] = activities

        program.region_enrichment = region_means

        # Skip statistical test if insufficient regions
        if len(region_activities) < 2:
            program.region_specificity = 0.0
            program.region_pvalue = 1.0
            program.enriched_region = None
            continue

        # Statistical test for region differences
        activity_lists = list(region_activities.values())
        if len(region_activities) == 2:
            # Mann-Whitney U for two regions
            _, pval = mannwhitneyu(activity_lists[0], activity_lists[1], alternative="two-sided")
        else:
            # Kruskal-Wallis for multiple regions
            _, pval = kruskal(*activity_lists)

        program.region_pvalue = float(pval)

        # Compute specificity score (coefficient of variation of region means)
        means = np.array(list(region_means.values()))
        if means.mean() > 0:
            program.region_specificity = float(means.std() / means.mean())
        else:
            program.region_specificity = 0.0

        # Identify enriched region
        if pval < 0.05 and len(region_means) > 0:
            program.enriched_region = max(region_means, key=region_means.get)
        else:
            program.enriched_region = None

    n_enriched = sum(1 for p in result.programs if p.enriched_region is not None)
    logger.info(
        "Region analysis complete: %s/%s programs show significant region enrichment", n_enriched, len(result.programs)
    )

    return result


def compare_programs_by_region(
    result: AnchoredProgramResult,
    adata: sc.AnnData,
    region_column: str,
    region_a: Any,
    region_b: Any,
    *,
    top_n_genes: int = 50,
) -> pd.DataFrame:
    """
    Compare program activities and gene loadings between two regions.

    For each program, computes:
    - Mean activity in region A vs B
    - Fold change and statistical significance
    - Top genes driving the program

    Useful for identifying which programs and genes are region-specific.

    Args:
        result: AnchoredProgramResult with discovered programs
        adata: AnnData with region annotations
        region_column: Column defining regions
        region_a: First region value (e.g., True for D538G+)
        region_b: Second region value (e.g., False for D538G-)
        top_n_genes: Number of top genes to include per program

    Returns:
        DataFrame with columns:
        - program_id: Program index
        - mean_activity_a, mean_activity_b: Mean activity per region
        - fold_change: region_a / region_b
        - pvalue: Mann-Whitney U test
        - top_genes: Comma-separated top genes
        - gene_loadings: Dict of gene -> loading

    Example:
        >>> df = compare_programs_by_region(result, adata, "D538G Mutation", True, False)
        >>> d538g_enriched = df[df['fold_change'] > 1.5]
    """
    from scipy.stats import mannwhitneyu  # pylint: disable=import-outside-toplevel

    mask_a = (adata.obs[region_column] == region_a).values
    mask_b = (adata.obs[region_column] == region_b).values

    H = result.H
    W = result.W

    records = []
    for k, _ in enumerate(result.programs):
        h_k = H[k, :]

        activities_a = h_k[mask_a]
        activities_b = h_k[mask_b]

        mean_a = float(np.mean(activities_a))
        mean_b = float(np.mean(activities_b))

        # Fold change (with pseudocount to avoid division by zero)
        fc = (mean_a + 1e-6) / (mean_b + 1e-6)

        # Statistical test
        if len(activities_a) >= 10 and len(activities_b) >= 10:
            _, pval = mannwhitneyu(activities_a, activities_b, alternative="two-sided")
        else:
            pval = 1.0

        # Top genes
        loadings = W[:, k]
        top_idx = np.argsort(loadings)[::-1][:top_n_genes]
        top_genes = [result.gene_names[i] for i in top_idx]
        gene_loadings = {result.gene_names[i]: float(loadings[i]) for i in top_idx}

        records.append(
            {
                "program_id": k,
                "mean_activity_a": mean_a,
                "mean_activity_b": mean_b,
                "fold_change": fc,
                "pvalue": float(pval),
                "n_spots_a": int(mask_a.sum()),
                "n_spots_b": int(mask_b.sum()),
                "top_genes": ", ".join(top_genes[:10]),
                "gene_loadings": gene_loadings,
            }
        )

    df = pd.DataFrame(records)
    df = df.sort_values("fold_change", ascending=False)

    return df


def extract_program_context_genes(
    result: AnchoredProgramResult,
    program_id: int,
    target_gene: str,
    top_n: int = 50,
    exclude_target: bool = True,
) -> List[Tuple[str, float]]:
    """
    Extract genes co-loaded with a target gene in a specific program.

    Useful for finding "contextual factors" - genes that are co-expressed
    with a gene of interest (e.g., MDK) in a spatial program.

    Args:
        result: AnchoredProgramResult with discovered programs
        program_id: Index of the program to analyze
        target_gene: Gene of interest (e.g., "MDK")
        top_n: Number of top co-loaded genes to return
        exclude_target: Whether to exclude target gene from results

    Returns:
        List of (gene_name, loading) tuples sorted by loading

    Example:
        >>> # Find genes co-expressed with MDK in program 2
        >>> context_genes = extract_program_context_genes(result, 2, "MDK")
        >>> print("Contextual factors:", [g[0] for g in context_genes[:10]])
    """
    if program_id >= len(result.programs):
        raise ValueError(f"Program {program_id} not found (max: {len(result.programs)-1})")

    W = result.W  # (n_genes, K_programs)
    loadings = W[:, program_id]

    # Check if target gene is in this program
    if target_gene in result.gene_names:
        target_idx = result.gene_names.index(target_gene)
        target_loading = loadings[target_idx]
        logger.info("Target gene '%s' loading in program %s: %s", target_gene, program_id, round(target_loading, 4))
    else:
        logger.warning("Target gene '%s' not found in gene list", target_gene)

    # Get top genes by loading
    top_idx = np.argsort(loadings)[::-1]

    context_genes = []
    for idx in top_idx:
        gene = result.gene_names[idx]
        if exclude_target and gene == target_gene:
            continue
        context_genes.append((gene, float(loadings[idx])))
        if len(context_genes) >= top_n:
            break

    return context_genes
