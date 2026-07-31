"""
Module 5: Cross-Sample Integration for Generalizable Program Relationships.

Integrates gene expression programs across multiple patient samples using
Harmony-style batch correction, producing a similarity network showing
which programs and bivariate relationships are conserved across patients.

Key concepts:
    - INTRA-sample: Bivariate Moran's I computed within each sample (Module 4b)
    - INTER-sample: We compare whether aligned program pairs show consistent
                   relationships across multiple patients
    - We do NOT recompute spatial statistics in integrated space

Pipeline integration:
    Module 4 (per sample)  -> W matrices, programs, H matrices
    Module 4b (per sample) -> Bivariate relationships
    Module 5 (THIS MODULE) -> Cross-sample integration + similarity network
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.special import softmax
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import cosine_similarity

try:
    from .anchored_program_discovery import (
        AnchoredProgramDiscoveryResult,
        AnchoredProgramResult,
        BivariateProgramResult,
        ProgramPairRelationship,
        SpatialProgram,
    )
except ImportError:
    from anchored_program_discovery import (  # type: ignore[no-redef]
        AnchoredProgramDiscoveryResult,
        AnchoredProgramResult,
        BivariateProgramResult,
        ProgramPairRelationship,
        SpatialProgram,
    )


logger = logging.getLogger(__name__)


# =============================================================================
# Data Structures
# =============================================================================


@dataclass
class AlignedProgram:
    """A program aligned across multiple samples."""

    program_id: str
    """Unique identifier for this aligned program (e.g., 'aligned_001')."""

    cell_type: str
    """Anchor cell type (e.g., 'Macrophage', 'T-cell')."""

    samples_present: List[str]
    """List of sample names where this program appears."""

    consensus_signature: Dict[str, float]
    """Meta-analyzed gene loadings (gene_name -> mean loading across samples)."""

    sample_signatures: Dict[str, NDArray]
    """Per-sample W vectors (sample_name -> (n_genes,) array)."""

    sample_program_ids: Dict[str, int]
    """Mapping of sample_name -> original program_id in that sample."""

    embedding: NDArray
    """Position in integrated latent space (n_components,)."""

    conservation_score: float
    """Fraction of samples where this program appears (0-1)."""

    top_genes: List[str]
    """Top N genes by consensus loading."""

    def __repr__(self) -> str:
        return (
            f"AlignedProgram({self.program_id}, cell_type={self.cell_type}, "
            f"samples={len(self.samples_present)}/{self.conservation_score:.0%})"
        )


@dataclass
class ConservedRelationship:
    """A bivariate relationship that appears across multiple samples."""

    program1_id: str
    """Aligned program 1 ID."""

    program2_id: str
    """Aligned program 2 ID."""

    samples_with_relationship: List[str]
    """List of samples where this relationship appears."""

    bivariate_i_values: List[float]
    """Per-sample bivariate Moran's I values."""

    mean_bivariate_i: float
    """Mean bivariate Moran's I across samples."""

    std_bivariate_i: float
    """Standard deviation of bivariate I across samples."""

    conservation_score: float
    """Fraction of samples with this relationship (0-1)."""

    relationship_type: str
    """'co-localized', 'exclusive', or 'variable'."""

    def __repr__(self) -> str:
        return (
            f"ConservedRelationship({self.program1_id} <-> {self.program2_id}: "
            f"I={self.mean_bivariate_i:.3f}+/-{self.std_bivariate_i:.3f}, "
            f"{self.relationship_type}, {len(self.samples_with_relationship)} samples)"
        )


@dataclass
class IntegrationResult:
    """Module 5 results: cross-sample integration."""

    aligned_programs: List[AlignedProgram]
    """List of all aligned programs across samples."""

    conserved_relationships: List[ConservedRelationship]
    """Bivariate relationships that persist across samples."""

    sample_names: List[str]
    """Names of all integrated samples."""

    n_programs_per_sample: Dict[str, int]
    """Number of programs per sample before alignment."""

    integration_embedding: NDArray
    """Full embedding matrix (n_total_programs, n_components)."""

    program_metadata: pd.DataFrame
    """Metadata for each program in the embedding."""

    similarity_graph: Optional[nx.Graph]
    """NetworkX graph of program similarities (if built)."""

    harmony_converged: bool
    """Whether Harmony iteration converged."""

    n_iterations: int
    """Number of Harmony iterations."""

    parameters: Dict[str, Any] = field(default_factory=dict)
    """Parameters used for integration."""

    def summary(self) -> str:
        """Return a summary string."""
        n_samples = len(self.sample_names)
        n_aligned = len(self.aligned_programs)
        n_conserved = len([p for p in self.aligned_programs if p.conservation_score > 0.5])
        n_relationships = len(self.conserved_relationships)
        n_strong_rel = len([r for r in self.conserved_relationships if r.conservation_score > 0.5])

        lines = [
            "Module 5: Cross-Sample Integration",
            "=" * 40,
            f"Samples integrated: {n_samples}",
            f"Aligned program groups: {n_aligned}",
            f"  Highly conserved (>50%): {n_conserved}",
            "",
            f"Conserved relationships: {n_relationships}",
            f"  Strong (>50% samples): {n_strong_rel}",
            "",
            f"Harmony converged: {self.harmony_converged}",
            f"Iterations: {self.n_iterations}",
        ]

        # Top conserved programs
        if self.aligned_programs:
            lines.append("")
            lines.append("Top conserved programs:")
            sorted_progs = sorted(self.aligned_programs, key=lambda p: p.conservation_score, reverse=True)
            for prog in sorted_progs[:5]:
                lines.append(
                    f"  {prog.program_id} ({prog.cell_type}): "
                    f"{len(prog.samples_present)}/{n_samples} samples, "
                    f"top genes: {', '.join(prog.top_genes[:3])}"
                )

        # Top conserved relationships
        if self.conserved_relationships:
            lines.append("")
            lines.append("Top conserved relationships:")
            sorted_rels = sorted(self.conserved_relationships, key=lambda r: r.conservation_score, reverse=True)
            for rel in sorted_rels[:5]:
                lines.append(
                    f"  {rel.program1_id} <-> {rel.program2_id}: "
                    f"I={rel.mean_bivariate_i:.3f}, {rel.relationship_type}, "
                    f"{len(rel.samples_with_relationship)} samples"
                )

        return "\n".join(lines)

    def to_programs_dataframe(self) -> pd.DataFrame:
        """Convert aligned programs to DataFrame."""
        records = []
        for prog in self.aligned_programs:
            records.append(
                {
                    "program_id": prog.program_id,
                    "cell_type": prog.cell_type,
                    "n_samples": len(prog.samples_present),
                    "conservation_score": prog.conservation_score,
                    "samples": ",".join(prog.samples_present),
                    "top_genes": ",".join(prog.top_genes[:10]),
                }
            )
        return pd.DataFrame(records)

    def to_relationships_dataframe(self) -> pd.DataFrame:
        """Convert conserved relationships to DataFrame."""
        records = []
        for rel in self.conserved_relationships:
            records.append(
                {
                    "program1_id": rel.program1_id,
                    "program2_id": rel.program2_id,
                    "n_samples": len(rel.samples_with_relationship),
                    "conservation_score": rel.conservation_score,
                    "mean_bivariate_i": rel.mean_bivariate_i,
                    "std_bivariate_i": rel.std_bivariate_i,
                    "relationship_type": rel.relationship_type,
                    "samples": ",".join(rel.samples_with_relationship),
                }
            )
        return pd.DataFrame(records)


# =============================================================================
# Loading Functions
# =============================================================================


def load_module4_from_csv(
    csv_path: Path,
    sample_name: Optional[str] = None,
) -> Dict[str, AnchoredProgramResult]:
    """
    Load Module 4 results from CSV file.

    Args:
        csv_path: Path to CSV file (from AnchoredProgramDiscoveryResult.to_dataframe())
        sample_name: Sample name to use (extracted from filename if not provided)

    Returns:
        Dictionary mapping anchor names to AnchoredProgramResult (partial - W only)
    """
    if sample_name is None:
        sample_name = csv_path.stem.replace("_module4_v3_programs", "").replace("_module4_programs", "")

    # Try loading with and without index column
    df = pd.read_csv(csv_path)

    # Check if first column is an unnamed index
    if df.columns[0] == "" or df.columns[0].startswith("Unnamed"):
        df = pd.read_csv(csv_path, index_col=0)

    # Ensure 'anchor' column exists
    if "anchor" not in df.columns:
        logger.warning("No 'anchor' column in %s", csv_path)
        return {}

    results = {}
    for anchor in df["anchor"].unique():
        anchor_df = df[df["anchor"] == anchor]

        # Parse gene loadings from top_genes column
        # This is a simplified reconstruction - full W matrix would need separate storage
        programs = []
        for _, row in anchor_df.iterrows():
            top_genes = row["top_genes"].split(", ") if pd.notna(row["top_genes"]) else []
            programs.append(
                SpatialProgram(
                    program_id=int(row["program_id"]),
                    top_genes=top_genes,
                    gene_loadings={g: 1.0 / (i + 1) for i, g in enumerate(top_genes)},  # Approximate
                    variance_explained=row["variance_explained"],
                    spatial_moran_i=row["spatial_moran_i"],
                    spatial_moran_pvalue=row["spatial_moran_pvalue"],
                    mean_activity=row["mean_activity"],
                    active_spots_fraction=row["active_spots_fraction"],
                )
            )

        # Create partial result (W matrix approximate from top genes)
        all_genes = list(set(g for p in programs for g in p.top_genes))
        n_progs = len(programs)
        W = np.zeros((len(all_genes), n_progs))
        for k, prog in enumerate(programs):
            for g, loading in prog.gene_loadings.items():
                if g in all_genes:
                    W[all_genes.index(g), k] = loading

        results[anchor] = AnchoredProgramResult(
            anchor_name=anchor,
            anchor_proteins=[],
            programs=programs,
            W=W,
            H=np.zeros((n_progs, 1)),  # Placeholder
            gene_names=all_genes,
            protein_correlations=pd.DataFrame(),
            reconstruction_error=0.0,
            n_spots_used=0,
            parameters={"sample_name": sample_name},
        )

    return results


def load_module4b_from_csv(csv_path: Path) -> BivariateProgramResult:
    """
    Load Module 4b results from CSV file.

    Args:
        csv_path: Path to CSV file (from BivariateProgramResult.to_dataframe())

    Returns:
        BivariateProgramResult
    """
    # Load without assuming index column
    df = pd.read_csv(csv_path)

    # Check if first column is an unnamed index
    if df.columns[0] == "" or df.columns[0].startswith("Unnamed"):
        df = pd.read_csv(csv_path, index_col=0)

    all_pairs = []
    for _, row in df.iterrows():
        pair = ProgramPairRelationship(
            anchor1=row["anchor1"],
            program1_id=int(row["program1_id"]),
            anchor2=row["anchor2"],
            program2_id=int(row["program2_id"]),
            bivariate_morans_i=row["bivariate_morans_i"],
            bivariate_pvalue=row.get("bivariate_pvalue", 0.05),
            pearson_correlation=row.get("pearson_correlation", 0.0),
            pearson_pvalue=row.get("pearson_pvalue", 1.0),
            n_spots_used=int(row.get("n_spots_used", 0)),
            top_genes_overlap=(
                row.get("top_genes_overlap", "").split(",") if pd.notna(row.get("top_genes_overlap", "")) else []
            ),
            relationship_type=row.get("relationship_type", "unknown"),
        )
        all_pairs.append(pair)

    significant = [p for p in all_pairs if p.bivariate_pvalue < 0.05]

    return BivariateProgramResult(
        significant_pairs=significant,
        all_pairs=all_pairs,
        n_programs_total=len(
            set([(p.anchor1, p.program1_id) for p in all_pairs] + [(p.anchor2, p.program2_id) for p in all_pairs])
        ),
        n_pairs_tested=len(all_pairs),
        n_significant=len(significant),
        fdr_threshold=0.05,
    )


def load_multi_sample_results(
    sample_dirs: List[Path],
    module4_pattern: str = "*_module4*programs.csv",
    module4b_pattern: str = "*_module4b*relationships.csv",
) -> Tuple[Dict[str, AnchoredProgramDiscoveryResult], Dict[str, BivariateProgramResult],]:
    """
    Load Module 4 and 4b results from multiple sample directories.

    Args:
        sample_dirs: List of directories containing sample outputs
        module4_pattern: Glob pattern for Module 4 CSV files
        module4b_pattern: Glob pattern for Module 4b CSV files

    Returns:
        Tuple of (module4_results, module4b_results) dictionaries keyed by sample name
    """
    module4_results = {}
    module4b_results = {}

    for sample_dir in sample_dirs:
        sample_dir = Path(sample_dir)

        # Find Module 4 files
        m4_files = list(sample_dir.glob(module4_pattern))
        if m4_files:
            m4_path = m4_files[0]
            sample_name = m4_path.stem.split("_module4")[0]
            results_by_anchor = load_module4_from_csv(m4_path, sample_name)

            module4_results[sample_name] = AnchoredProgramDiscoveryResult(
                results_by_anchor=results_by_anchor,
                n_anchors=len(results_by_anchor),
                total_programs=sum(len(r.programs) for r in results_by_anchor.values()),
                profile_discovery_result=None,
                parameters={"source_file": str(m4_path)},
            )
            logger.info("Loaded Module 4 results for %s: %s anchors", sample_name, len(results_by_anchor))

        # Find Module 4b files
        m4b_files = list(sample_dir.glob(module4b_pattern))
        if m4b_files:
            m4b_path = m4b_files[0]
            sample_name = m4b_path.stem.split("_module4b")[0]
            module4b_results[sample_name] = load_module4b_from_csv(m4b_path)
            logger.info(
                "Loaded Module 4b results for %s: %s pairs",
                sample_name,
                module4b_results[sample_name].n_pairs_tested,
            )

    return module4_results, module4b_results


# =============================================================================
# Gene Set Alignment
# =============================================================================


def align_gene_sets(
    results_by_sample: Dict[str, AnchoredProgramDiscoveryResult],
) -> Tuple[Dict[Tuple[str, str], NDArray], List[str], pd.DataFrame]:
    """
    Align gene sets across samples by creating a union of all genes.

    Args:
        results_by_sample: Module 4 results keyed by sample name

    Returns:
        Tuple of:
            - aligned_W: Dict[(sample, cell_type), W_aligned] with zero-padded matrices
            - all_genes: List of all genes in union
            - metadata: DataFrame with program metadata (sample, cell_type, program_id)
    """
    # Find union of all genes
    all_genes_set = set()
    for sample_name, result in results_by_sample.items():
        for anchor_result in result.results_by_anchor.values():
            all_genes_set.update(anchor_result.gene_names)

    all_genes = sorted(list(all_genes_set))
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}

    logger.info("Gene set alignment: %s genes across %s samples", len(all_genes), len(results_by_sample))

    # Create aligned W matrices
    aligned_W: Dict[Tuple[str, str], NDArray] = {}
    metadata_records = []

    for sample_name, result in results_by_sample.items():
        for cell_type, anchor_result in result.results_by_anchor.items():
            n_programs = anchor_result.W.shape[1]
            W_aligned = np.zeros((len(all_genes), n_programs))

            # Fill in known gene loadings
            for i, gene in enumerate(anchor_result.gene_names):
                if gene in gene_to_idx:
                    W_aligned[gene_to_idx[gene], :] = anchor_result.W[i, :]

            aligned_W[(sample_name, cell_type)] = W_aligned

            # Record metadata for each program
            for k in range(n_programs):
                metadata_records.append(
                    {
                        "sample": sample_name,
                        "cell_type": cell_type,
                        "program_id": k,
                        "program_key": f"{sample_name}_{cell_type}_prog{k}",
                    }
                )

    metadata = pd.DataFrame(metadata_records)

    logger.info("Created %s aligned W matrices with %s total programs", len(aligned_W), len(metadata))

    return aligned_W, all_genes, metadata


# =============================================================================
# Harmony-Style Integration
# =============================================================================


def integrate_programs_harmony(
    aligned_W: Dict[Tuple[str, str], NDArray],
    metadata: pd.DataFrame,  # pylint: disable=unused-argument
    *,
    n_components: int = 30,
    n_clusters: int = 50,
    theta: float = 2.0,
    max_iter: int = 20,
    tol: float = 1e-4,
    random_state: int = 42,
) -> Tuple[NDArray, NDArray, Dict[str, Any]]:
    """
    Harmony-style integration of program signatures across samples.

    Algorithm:
    1. Stack all W matrices (columns as programs)
    2. PCA to reduce dimensionality
    3. Iterative soft clustering + batch correction
    4. Return integrated embedding + cluster assignments

    Args:
        aligned_W: Dict[(sample, cell_type), W_aligned] from align_gene_sets
        metadata: DataFrame with sample, cell_type, program_id columns
        n_components: Number of PCA components
        n_clusters: Number of soft clusters for Harmony
        theta: Diversity penalty (higher = more mixing)
        max_iter: Maximum iterations
        tol: Convergence tolerance
        random_state: Random seed

    Returns:
        Tuple of:
            - Z: Integrated embedding (n_programs, n_components)
            - R: Soft cluster assignments (n_programs, n_clusters)
            - info: Dict with convergence info
    """
    np.random.seed(random_state)

    # Stack all W matrices (each program is a column)
    W_list: List[Any] = []
    sample_labels_list: List[str] = []

    for (sample, _), W in aligned_W.items():
        n_progs = W.shape[1]
        for k in range(n_progs):
            W_list.append(W[:, k])
            sample_labels_list.append(sample)

    W_stacked = np.column_stack(W_list)  # (n_genes, n_programs)
    sample_labels: NDArray[Any] = np.array(sample_labels_list)
    unique_samples = np.unique(sample_labels)

    n_programs = W_stacked.shape[1]
    n_genes = W_stacked.shape[0]

    logger.info("Harmony integration: %s programs, %s genes, %s samples", n_programs, n_genes, len(unique_samples))

    # Adjust n_components if needed
    actual_components = min(n_components, n_programs, n_genes)

    # PCA reduction
    pca = PCA(n_components=actual_components, random_state=random_state)
    Z = pca.fit_transform(W_stacked.T)  # (n_programs, n_components)

    logger.info(
        "PCA: %s%% variance explained with %s components",
        round(pca.explained_variance_ratio_.sum() * 100, 1),
        actual_components,
    )

    # L2 normalize
    norms = np.linalg.norm(Z, axis=1, keepdims=True)
    norms[norms == 0] = 1
    Z = Z / norms

    # Adjust n_clusters if needed
    actual_clusters = min(n_clusters, n_programs)

    # Initialize soft clustering with k-means
    kmeans = KMeans(n_clusters=actual_clusters, random_state=random_state, n_init=10)
    kmeans.fit(Z)
    C = kmeans.cluster_centers_  # (n_clusters, n_components)

    # Initialize soft assignments
    distances = pairwise_distances(Z, C)
    R = softmax(-theta * distances, axis=1)  # (n_programs, n_clusters)

    # Iterative Harmony correction
    converged = False
    prev_Z = Z.copy()

    for iteration in range(max_iter):
        # Update cluster centroids (global)
        R_sum = R.sum(axis=0, keepdims=True).T  # (n_clusters, 1)
        R_sum[R_sum == 0] = 1
        C = (R.T @ Z) / R_sum  # (n_clusters, n_components)

        # Sample-specific correction
        for sample in unique_samples:
            mask = sample_labels == sample

            if mask.sum() == 0:
                continue

            # Sample-specific centroids
            R_sample = R[mask]
            Z_sample = Z[mask]

            R_sample_sum = R_sample.sum(axis=0, keepdims=True).T
            R_sample_sum[R_sample_sum == 0] = 1
            C_sample = (R_sample.T @ Z_sample) / R_sample_sum

            # Correction factor
            correction = C - C_sample

            # Apply weighted correction
            Z[mask] += R_sample @ correction

        # Re-normalize
        norms = np.linalg.norm(Z, axis=1, keepdims=True)
        norms[norms == 0] = 1
        Z = Z / norms

        # Update cluster assignments
        distances = pairwise_distances(Z, C)
        R = softmax(-theta * distances, axis=1)

        # Check convergence
        delta = np.mean(np.abs(Z - prev_Z))
        if delta < tol:
            converged = True
            logger.info("Harmony converged at iteration %s (delta=%s)", iteration + 1, round(delta, 6))
            break

        prev_Z = Z.copy()

        if (iteration + 1) % 5 == 0:
            logger.info("Harmony iteration %s: delta=%s", iteration + 1, round(delta, 6))

    if not converged:
        logger.warning("Harmony did not converge after %s iterations", max_iter)

    info = {
        "converged": converged,
        "n_iterations": iteration + 1,
        "pca_variance_explained": pca.explained_variance_ratio_.sum(),
        "n_components": actual_components,
        "n_clusters": actual_clusters,
        "theta": theta,
    }

    return Z, R, info


# =============================================================================
# Program Matching
# =============================================================================


def match_programs_across_samples(
    Z: NDArray,
    metadata: pd.DataFrame,
    aligned_W: Dict[Tuple[str, str], NDArray],
    all_genes: List[str],
    *,
    similarity_threshold: float = 0.7,
    top_n_genes: int = 20,
) -> List[AlignedProgram]:
    """
    Match programs across samples based on embedding similarity.

    Args:
        Z: Integrated embedding (n_programs, n_components)
        metadata: DataFrame with sample, cell_type, program_id columns
        aligned_W: Dict of aligned W matrices
        all_genes: List of all genes
        similarity_threshold: Minimum cosine similarity to group programs
        top_n_genes: Number of top genes to include

    Returns:
        List of AlignedProgram objects
    """
    n_samples = len(metadata["sample"].unique())

    # Compute pairwise cosine similarity
    similarity = cosine_similarity(Z)

    # Convert similarity to distance for clustering
    distance_matrix = 1 - similarity
    np.fill_diagonal(distance_matrix, 0)

    # Hierarchical clustering
    try:
        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=1 - similarity_threshold,
            metric="precomputed",
            linkage="average",
        )
        labels = clustering.fit_predict(distance_matrix)
    except (ValueError, RuntimeError) as e:
        logger.warning("Clustering failed: %s. Using fallback method.", e)
        # Fallback: each program is its own cluster
        labels = np.arange(len(Z))

    logger.info("Program matching: %s aligned program groups from %s total programs", len(np.unique(labels)), len(Z))

    # Build aligned programs
    aligned_programs = []

    for cluster_id in np.unique(labels):
        mask = labels == cluster_id
        cluster_meta = metadata[mask]
        cluster_indices = np.where(mask)[0]

        # Use the most common cell type as the representative
        cell_type = cluster_meta["cell_type"].mode().iloc[0] if len(cluster_meta) > 0 else "Unknown"

        # Collect per-sample signatures
        sample_signatures = {}
        sample_program_ids = {}

        for _, (_, row) in zip(cluster_indices, cluster_meta.iterrows()):
            sample = row["sample"]
            ct = row["cell_type"]
            prog_id = row["program_id"]

            if (sample, ct) in aligned_W:
                W = aligned_W[(sample, ct)]
                if prog_id < W.shape[1]:
                    sample_signatures[sample] = W[:, prog_id]
                    sample_program_ids[sample] = prog_id

        samples_present = list(sample_signatures.keys())

        if len(samples_present) == 0:
            continue

        # Compute consensus signature (mean across samples)
        sig_matrix = np.column_stack(list(sample_signatures.values()))
        consensus = sig_matrix.mean(axis=1)

        # Top genes by consensus loading
        top_indices = np.argsort(np.abs(consensus))[::-1][:top_n_genes]
        top_genes = [all_genes[i] for i in top_indices]

        consensus_dict = {all_genes[i]: consensus[i] for i in range(len(all_genes))}

        # Mean embedding
        mean_embedding = Z[mask].mean(axis=0)

        aligned_programs.append(
            AlignedProgram(
                program_id=f"aligned_{cluster_id:03d}",
                cell_type=cell_type,
                samples_present=samples_present,
                consensus_signature=consensus_dict,
                sample_signatures=sample_signatures,
                sample_program_ids=sample_program_ids,
                embedding=mean_embedding,
                conservation_score=len(samples_present) / n_samples,
                top_genes=top_genes,
            )
        )

    # Sort by conservation score
    aligned_programs.sort(key=lambda p: p.conservation_score, reverse=True)

    logger.info("Created %s aligned programs", len(aligned_programs))

    return aligned_programs


# =============================================================================
# Bivariate Relationship Comparison
# =============================================================================


def compare_bivariate_relationships(
    bivariate_results: Dict[str, BivariateProgramResult],
    aligned_programs: List[AlignedProgram],
    min_samples: int = 2,
    _significance_threshold: float = 0.05,
    colocalization_threshold: float = 0.1,
) -> List[ConservedRelationship]:
    """
    Compare bivariate relationships across samples using aligned program IDs.

    IMPORTANT: This compares INTRA-sample relationships (computed by Module 4b)
    to assess their persistence INTER-sample (across patients).

    Args:
        bivariate_results: Module 4b results keyed by sample name
        aligned_programs: List of aligned programs with sample mappings
        min_samples: Minimum samples for a relationship to be considered
        significance_threshold: P-value threshold for significant relationships
        colocalization_threshold: Moran's I threshold for co-localization

    Returns:
        List of ConservedRelationship objects
    """
    # Build reverse mapping: (sample, cell_type, program_id) -> aligned_program_id
    program_to_aligned = {}
    for aligned_prog in aligned_programs:
        for sample, prog_id in aligned_prog.sample_program_ids.items():
            key = (sample, aligned_prog.cell_type, prog_id)
            program_to_aligned[key] = aligned_prog.program_id

    # Collect bivariate I values per aligned program pair
    relationship_data: Dict[Any, Any] = defaultdict(lambda: {"values": [], "samples": []})

    for sample_name, result in bivariate_results.items():
        for pair in result.all_pairs:
            # Map sample-specific programs to aligned IDs
            key1 = (sample_name, pair.anchor1, pair.program1_id)
            key2 = (sample_name, pair.anchor2, pair.program2_id)

            aligned1 = program_to_aligned.get(key1)
            aligned2 = program_to_aligned.get(key2)

            if aligned1 is None or aligned2 is None:
                continue

            # Use sorted tuple as key to ensure consistency
            pair_key = tuple(sorted([aligned1, aligned2]))

            relationship_data[pair_key]["values"].append(pair.bivariate_morans_i)
            relationship_data[pair_key]["samples"].append(sample_name)

    # Convert to ConservedRelationship objects
    n_samples = len(bivariate_results)
    conserved_relationships = []

    for (prog1, prog2), data in relationship_data.items():
        if len(data["samples"]) < min_samples:
            continue

        values = np.array(data["values"])
        mean_i = np.mean(values)
        std_i = np.std(values)

        # Classify relationship type
        if mean_i > colocalization_threshold:
            rel_type = "co-localized"
        elif mean_i < -colocalization_threshold:
            rel_type = "exclusive"
        else:
            rel_type = "independent"

        conserved_relationships.append(
            ConservedRelationship(
                program1_id=prog1,
                program2_id=prog2,
                samples_with_relationship=data["samples"],
                bivariate_i_values=data["values"],
                mean_bivariate_i=mean_i,
                std_bivariate_i=std_i,
                conservation_score=len(data["samples"]) / n_samples,
                relationship_type=rel_type,
            )
        )

    # Sort by conservation score
    conserved_relationships.sort(key=lambda r: r.conservation_score, reverse=True)

    logger.info(
        "Found %s relationships across samples, %s conserved in >50%% samples",
        len(conserved_relationships),
        sum(1 for r in conserved_relationships if r.conservation_score > 0.5),
    )

    return conserved_relationships


# =============================================================================
# Similarity Network
# =============================================================================


def build_similarity_network(
    aligned_programs: List[AlignedProgram],
    conserved_relationships: List[ConservedRelationship],
    min_program_conservation: float = 0.3,
    min_relationship_conservation: float = 0.3,
) -> nx.Graph:
    """
    Build NetworkX graph of conserved programs and relationships.

    Args:
        aligned_programs: List of aligned programs
        conserved_relationships: List of conserved relationships
        min_program_conservation: Minimum conservation to include program as node
        min_relationship_conservation: Minimum conservation to include edge

    Returns:
        NetworkX Graph with programs as nodes and relationships as edges
    """
    G = nx.Graph()

    # Add program nodes
    for prog in aligned_programs:
        if prog.conservation_score >= min_program_conservation:
            G.add_node(
                prog.program_id,
                cell_type=prog.cell_type,
                conservation=prog.conservation_score,
                n_samples=len(prog.samples_present),
                top_genes=",".join(prog.top_genes[:5]),
            )

    # Add relationship edges
    for rel in conserved_relationships:
        if rel.conservation_score >= min_relationship_conservation:
            if rel.program1_id in G.nodes and rel.program2_id in G.nodes:
                G.add_edge(
                    rel.program1_id,
                    rel.program2_id,
                    weight=rel.mean_bivariate_i,
                    conservation=rel.conservation_score,
                    std=rel.std_bivariate_i,
                    type=rel.relationship_type,
                    n_samples=len(rel.samples_with_relationship),
                )

    logger.info("Built similarity network: %s nodes, %s edges", G.number_of_nodes(), G.number_of_edges())

    return G


# =============================================================================
# Main Integration Function
# =============================================================================


def integrate_samples(
    module4_results: Dict[str, AnchoredProgramDiscoveryResult],
    module4b_results: Optional[Dict[str, BivariateProgramResult]] = None,
    *,
    n_components: int = 30,
    n_clusters: int = 50,
    theta: float = 2.0,
    similarity_threshold: float = 0.7,
    max_iter: int = 20,
    random_state: int = 42,
    build_network: bool = True,
) -> IntegrationResult:
    """
    Main function: integrate programs across multiple samples.

    Args:
        module4_results: Dict of Module 4 results keyed by sample name
        module4b_results: Optional dict of Module 4b results for relationship analysis
        n_components: Number of PCA components for Harmony
        n_clusters: Number of soft clusters for Harmony
        theta: Diversity penalty (higher = more mixing)
        similarity_threshold: Threshold for program matching
        max_iter: Maximum Harmony iterations
        random_state: Random seed
        build_network: Whether to build similarity network

    Returns:
        IntegrationResult with aligned programs, relationships, and network
    """
    logger.info("Starting cross-sample integration with %s samples", len(module4_results))

    # Step 1: Align gene sets
    aligned_W, all_genes, metadata = align_gene_sets(module4_results)

    # Step 2: Harmony integration
    Z, _, harmony_info = integrate_programs_harmony(
        aligned_W=aligned_W,
        metadata=metadata,
        n_components=n_components,
        n_clusters=n_clusters,
        theta=theta,
        max_iter=max_iter,
        random_state=random_state,
    )

    # Step 3: Match programs across samples
    aligned_programs = match_programs_across_samples(
        Z=Z,
        metadata=metadata,
        aligned_W=aligned_W,
        all_genes=all_genes,
        similarity_threshold=similarity_threshold,
    )

    # Step 4: Compare bivariate relationships (if available)
    conserved_relationships = []
    if module4b_results is not None and len(module4b_results) > 0:
        conserved_relationships = compare_bivariate_relationships(
            bivariate_results=module4b_results,
            aligned_programs=aligned_programs,
        )

    # Step 5: Build similarity network
    similarity_graph = None
    if build_network:
        similarity_graph = build_similarity_network(
            aligned_programs=aligned_programs,
            conserved_relationships=conserved_relationships,
        )

    # Compute n_programs per sample
    n_programs_per_sample = {sample: result.total_programs for sample, result in module4_results.items()}

    result = IntegrationResult(
        aligned_programs=aligned_programs,
        conserved_relationships=conserved_relationships,
        sample_names=list(module4_results.keys()),
        n_programs_per_sample=n_programs_per_sample,
        integration_embedding=Z,
        program_metadata=metadata,
        similarity_graph=similarity_graph,
        harmony_converged=harmony_info["converged"],
        n_iterations=harmony_info["n_iterations"],
        parameters={
            "n_components": n_components,
            "n_clusters": n_clusters,
            "theta": theta,
            "similarity_threshold": similarity_threshold,
            "pca_variance_explained": harmony_info["pca_variance_explained"],
        },
    )

    logger.info(result.summary())

    return result


def save_integration_results(
    result: IntegrationResult,
    output_dir: Path,
    prefix: str = "module5",
) -> Dict[str, Path]:
    """
    Save integration results to files.

    Args:
        result: IntegrationResult to save
        output_dir: Output directory
        prefix: Prefix for output files

    Returns:
        Dictionary of file types to paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    # Save aligned programs
    programs_file = output_dir / f"{prefix}_aligned_programs.csv"
    result.to_programs_dataframe().to_csv(programs_file, index=False)
    saved_files["programs"] = programs_file

    # Save conserved relationships
    relationships_file = output_dir / f"{prefix}_conserved_relationships.csv"
    result.to_relationships_dataframe().to_csv(relationships_file, index=False)
    saved_files["relationships"] = relationships_file

    # Save embedding
    embedding_file = output_dir / f"{prefix}_embedding.npy"
    np.save(embedding_file, result.integration_embedding)
    saved_files["embedding"] = embedding_file

    # Save metadata
    metadata_file = output_dir / f"{prefix}_program_metadata.csv"
    result.program_metadata.to_csv(metadata_file, index=False)
    saved_files["metadata"] = metadata_file

    # Save network if available
    if result.similarity_graph is not None:
        network_file = output_dir / f"{prefix}_similarity_network.graphml"
        nx.write_graphml(result.similarity_graph, network_file)
        saved_files["network"] = network_file

    logger.info("Saved integration results to %s", output_dir)

    return saved_files
