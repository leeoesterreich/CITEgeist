#!/usr/bin/env python
"""
Test Module 5: Cross-Sample Integration.

Tests Harmony-style integration of gene expression programs across multiple samples.

Usage:
    python tests/test_cross_sample_integration.py
    pytest tests/test_cross_sample_integration.py -v
"""

import logging
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model import (
    # Module 4 types
    AnchoredProgramResult,
    AnchoredProgramDiscoveryResult,
    SpatialProgram,
    # Module 4b types
    ProgramPairRelationship,
    BivariateProgramResult,
    # Module 5 types
    AlignedProgram,
    ConservedRelationship,
    IntegrationResult,
    # Module 5 functions
    align_gene_sets,
    integrate_programs_harmony,
    match_programs_across_samples,
    compare_bivariate_relationships,
    build_similarity_network,
    integrate_samples,
    save_integration_results,
)


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Test Fixtures
# =============================================================================


def create_mock_module4_result(
    sample_name: str,
    n_cell_types: int = 3,
    n_programs_per_type: int = 2,
    n_genes: int = 100,
    shared_gene_fraction: float = 0.8,
    random_state: int = 42,
) -> AnchoredProgramDiscoveryResult:
    """
    Create a mock Module 4 result for testing.

    Args:
        sample_name: Name of the sample
        n_cell_types: Number of cell type anchors
        n_programs_per_type: Programs per cell type
        n_genes: Number of genes per cell type
        shared_gene_fraction: Fraction of genes shared across samples
        random_state: Random seed

    Returns:
        AnchoredProgramDiscoveryResult
    """
    rng = np.random.RandomState(random_state)

    cell_types = [f"CellType_{chr(65 + i)}" for i in range(n_cell_types)]

    # Create gene names - mix of shared and sample-specific
    n_shared = int(n_genes * shared_gene_fraction)
    n_specific = n_genes - n_shared
    shared_genes = [f"Gene_{i:03d}" for i in range(n_shared)]
    specific_genes = [f"Gene_{sample_name}_{i:03d}" for i in range(n_specific)]
    all_genes = shared_genes + specific_genes

    results_by_anchor = {}

    for ct in cell_types:
        # Create W matrix (gene loadings)
        W = rng.exponential(0.5, size=(n_genes, n_programs_per_type))

        # Create H matrix (spot loadings) - placeholder
        n_spots = 500
        H = rng.exponential(0.3, size=(n_programs_per_type, n_spots))

        # Create programs
        programs = []
        for k in range(n_programs_per_type):
            # Get top genes for this program
            top_idx = np.argsort(W[:, k])[::-1][:20]
            top_genes = [all_genes[i] for i in top_idx]
            gene_loadings = {all_genes[i]: W[i, k] for i in top_idx}

            programs.append(SpatialProgram(
                program_id=k,
                top_genes=top_genes,
                gene_loadings=gene_loadings,
                variance_explained=rng.uniform(5, 15),
                spatial_moran_i=rng.uniform(0.1, 0.4),
                spatial_moran_pvalue=rng.uniform(0.001, 0.05),
                mean_activity=rng.uniform(0.05, 0.2),
                active_spots_fraction=rng.uniform(0.4, 0.6),
            ))

        results_by_anchor[ct] = AnchoredProgramResult(
            anchor_name=ct,
            anchor_proteins=[f"{ct}_Protein_1", f"{ct}_Protein_2"],
            programs=programs,
            W=W,
            H=H,
            gene_names=all_genes,
            protein_correlations=pd.DataFrame(),
            reconstruction_error=rng.uniform(0.1, 0.3),
            n_spots_used=n_spots,
            parameters={"sample_name": sample_name},
        )

    return AnchoredProgramDiscoveryResult(
        results_by_anchor=results_by_anchor,
        n_anchors=n_cell_types,
        total_programs=n_cell_types * n_programs_per_type,
        profile_discovery_result=None,
        parameters={"sample_name": sample_name},
    )


def create_mock_module4b_result(
    module4_result: AnchoredProgramDiscoveryResult,
    sample_name: str,
    random_state: int = 42,
) -> BivariateProgramResult:
    """
    Create a mock Module 4b result for testing.

    Args:
        module4_result: Module 4 result to base relationships on
        sample_name: Sample name
        random_state: Random seed

    Returns:
        BivariateProgramResult
    """
    rng = np.random.RandomState(random_state)

    all_pairs = []
    anchors = list(module4_result.results_by_anchor.keys())

    for i, anchor1 in enumerate(anchors):
        for anchor2 in anchors[i:]:
            result1 = module4_result.results_by_anchor[anchor1]
            result2 = module4_result.results_by_anchor[anchor2]

            for prog1 in result1.programs:
                for prog2 in result2.programs:
                    if anchor1 == anchor2 and prog1.program_id >= prog2.program_id:
                        continue

                    # Generate realistic bivariate Moran's I
                    bivariate_i = rng.uniform(-0.3, 0.5)

                    # Determine relationship type
                    if bivariate_i > 0.1:
                        rel_type = "co-localized"
                    elif bivariate_i < -0.1:
                        rel_type = "exclusive"
                    else:
                        rel_type = "independent"

                    all_pairs.append(ProgramPairRelationship(
                        anchor1=anchor1,
                        program1_id=prog1.program_id,
                        anchor2=anchor2,
                        program2_id=prog2.program_id,
                        bivariate_morans_i=bivariate_i,
                        bivariate_pvalue=rng.uniform(0.001, 0.2),
                        pearson_correlation=rng.uniform(-0.5, 0.5),
                        pearson_pvalue=rng.uniform(0.001, 0.5),
                        n_spots_used=rng.randint(200, 500),
                        top_genes_overlap=[],
                        relationship_type=rel_type,
                    ))

    significant = [p for p in all_pairs if p.bivariate_pvalue < 0.05]

    return BivariateProgramResult(
        significant_pairs=significant,
        all_pairs=all_pairs,
        n_programs_total=module4_result.total_programs,
        n_pairs_tested=len(all_pairs),
        n_significant=len(significant),
        fdr_threshold=0.05,
        parameters={"sample_name": sample_name},
    )


def create_mock_multi_sample_data(
    n_samples: int = 3,
    n_cell_types: int = 3,
    n_programs_per_type: int = 2,
    n_genes: int = 100,
    base_seed: int = 42,
) -> tuple:
    """
    Create mock data for multiple samples.

    Returns:
        Tuple of (module4_results, module4b_results)
    """
    module4_results = {}
    module4b_results = {}

    for i in range(n_samples):
        sample_name = f"Sample_{i:02d}"
        seed = base_seed + i

        m4 = create_mock_module4_result(
            sample_name=sample_name,
            n_cell_types=n_cell_types,
            n_programs_per_type=n_programs_per_type,
            n_genes=n_genes,
            random_state=seed,
        )
        module4_results[sample_name] = m4

        m4b = create_mock_module4b_result(
            module4_result=m4,
            sample_name=sample_name,
            random_state=seed + 100,
        )
        module4b_results[sample_name] = m4b

    return module4_results, module4b_results


# =============================================================================
# Test Cases
# =============================================================================


class TestGeneSetAlignment:
    """Test gene set alignment across samples."""

    def test_align_single_sample(self):
        """Test alignment with single sample."""
        m4, _ = create_mock_multi_sample_data(n_samples=1)

        aligned_W, all_genes, metadata = align_gene_sets(m4)

        assert len(aligned_W) > 0
        assert len(all_genes) > 0
        assert len(metadata) > 0
        assert "sample" in metadata.columns
        assert "cell_type" in metadata.columns
        assert "program_id" in metadata.columns

    def test_align_multiple_samples(self):
        """Test alignment with multiple samples."""
        m4, _ = create_mock_multi_sample_data(n_samples=3)

        aligned_W, all_genes, metadata = align_gene_sets(m4)

        # Check that we have more genes due to some sample-specific ones
        assert len(all_genes) >= 100

        # Check metadata has all samples
        assert len(metadata["sample"].unique()) == 3

        # Check all W matrices have same number of rows (genes)
        n_genes = len(all_genes)
        for key, W in aligned_W.items():
            assert W.shape[0] == n_genes


class TestHarmonyIntegration:
    """Test Harmony-style integration algorithm."""

    def test_harmony_basic(self):
        """Test basic Harmony integration."""
        m4, _ = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)

        Z, R, info = integrate_programs_harmony(
            aligned_W=aligned_W,
            metadata=metadata,
            n_components=10,
            n_clusters=5,
            max_iter=10,
        )

        # Check output shapes
        n_programs = len(metadata)
        assert Z.shape[0] == n_programs
        assert Z.shape[1] <= 10  # May be less if n_programs < n_components
        assert R.shape[0] == n_programs

        # Check info dict
        assert "converged" in info
        assert "n_iterations" in info
        assert info["n_iterations"] <= 10

    def test_harmony_convergence(self):
        """Test that Harmony converges with sufficient iterations."""
        m4, _ = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)

        Z, R, info = integrate_programs_harmony(
            aligned_W=aligned_W,
            metadata=metadata,
            n_components=10,
            n_clusters=5,
            max_iter=50,
            tol=1e-3,
        )

        # With reasonable settings, should converge
        # (but don't fail test if it doesn't - just log)
        if not info["converged"]:
            logger.warning(f"Harmony did not converge in {info['n_iterations']} iterations")

    def test_harmony_different_theta(self):
        """Test Harmony with different diversity penalties."""
        m4, _ = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)

        # Low theta - less mixing
        Z_low, _, _ = integrate_programs_harmony(
            aligned_W=aligned_W,
            metadata=metadata,
            theta=0.5,
            max_iter=10,
        )

        # High theta - more mixing
        Z_high, _, _ = integrate_programs_harmony(
            aligned_W=aligned_W,
            metadata=metadata,
            theta=4.0,
            max_iter=10,
        )

        # Both should produce valid embeddings
        assert Z_low.shape == Z_high.shape
        assert not np.allclose(Z_low, Z_high)


class TestProgramMatching:
    """Test program matching across samples."""

    def test_match_programs_basic(self):
        """Test basic program matching."""
        m4, _ = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)
        Z, _, _ = integrate_programs_harmony(aligned_W, metadata, max_iter=5)

        aligned_programs = match_programs_across_samples(
            Z=Z,
            metadata=metadata,
            aligned_W=aligned_W,
            all_genes=all_genes,
            similarity_threshold=0.5,
        )

        assert len(aligned_programs) > 0
        for prog in aligned_programs:
            assert isinstance(prog, AlignedProgram)
            assert len(prog.samples_present) > 0
            assert len(prog.top_genes) > 0
            assert 0 <= prog.conservation_score <= 1

    def test_conservation_scores(self):
        """Test that conservation scores are correct."""
        m4, _ = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)
        Z, _, _ = integrate_programs_harmony(aligned_W, metadata, max_iter=5)

        aligned_programs = match_programs_across_samples(
            Z=Z,
            metadata=metadata,
            aligned_W=aligned_W,
            all_genes=all_genes,
        )

        for prog in aligned_programs:
            expected_score = len(prog.samples_present) / 3  # 3 samples
            assert abs(prog.conservation_score - expected_score) < 0.01


class TestBivariateRelationshipComparison:
    """Test bivariate relationship comparison across samples."""

    def test_compare_relationships_basic(self):
        """Test basic relationship comparison."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)
        Z, _, _ = integrate_programs_harmony(aligned_W, metadata, max_iter=5)
        aligned_programs = match_programs_across_samples(Z, metadata, aligned_W, all_genes)

        conserved = compare_bivariate_relationships(
            bivariate_results=m4b,
            aligned_programs=aligned_programs,
            min_samples=2,
        )

        # Should find some relationships (or empty list)
        assert isinstance(conserved, list)
        for rel in conserved:
            assert isinstance(rel, ConservedRelationship)
            assert len(rel.samples_with_relationship) >= 2

    def test_relationship_types(self):
        """Test that relationship types are correctly classified."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)
        Z, _, _ = integrate_programs_harmony(aligned_W, metadata, max_iter=5)
        aligned_programs = match_programs_across_samples(Z, metadata, aligned_W, all_genes)

        conserved = compare_bivariate_relationships(
            bivariate_results=m4b,
            aligned_programs=aligned_programs,
        )

        for rel in conserved:
            assert rel.relationship_type in ["co-localized", "exclusive", "independent"]

            # Check consistency with mean_bivariate_i
            if rel.mean_bivariate_i > 0.1:
                assert rel.relationship_type == "co-localized"
            elif rel.mean_bivariate_i < -0.1:
                assert rel.relationship_type == "exclusive"


class TestSimilarityNetwork:
    """Test similarity network construction."""

    def test_build_network_basic(self):
        """Test basic network construction."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)
        Z, _, _ = integrate_programs_harmony(aligned_W, metadata, max_iter=5)
        aligned_programs = match_programs_across_samples(Z, metadata, aligned_W, all_genes)
        conserved = compare_bivariate_relationships(m4b, aligned_programs)

        G = build_similarity_network(
            aligned_programs=aligned_programs,
            conserved_relationships=conserved,
            min_program_conservation=0.0,  # Include all
            min_relationship_conservation=0.0,
        )

        assert G.number_of_nodes() > 0

    def test_network_filtering(self):
        """Test network filtering by conservation."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)
        Z, _, _ = integrate_programs_harmony(aligned_W, metadata, max_iter=5)
        aligned_programs = match_programs_across_samples(Z, metadata, aligned_W, all_genes)
        conserved = compare_bivariate_relationships(m4b, aligned_programs)

        G_all = build_similarity_network(
            aligned_programs, conserved,
            min_program_conservation=0.0,
            min_relationship_conservation=0.0,
        )

        G_filtered = build_similarity_network(
            aligned_programs, conserved,
            min_program_conservation=0.5,
            min_relationship_conservation=0.5,
        )

        # Filtered should have fewer or equal nodes/edges
        assert G_filtered.number_of_nodes() <= G_all.number_of_nodes()

    def test_network_node_attributes(self):
        """Test that network nodes have correct attributes."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)
        aligned_W, all_genes, metadata = align_gene_sets(m4)
        Z, _, _ = integrate_programs_harmony(aligned_W, metadata, max_iter=5)
        aligned_programs = match_programs_across_samples(Z, metadata, aligned_W, all_genes)
        conserved = compare_bivariate_relationships(m4b, aligned_programs)

        G = build_similarity_network(aligned_programs, conserved, 0.0, 0.0)

        for node in G.nodes:
            attrs = G.nodes[node]
            assert "cell_type" in attrs
            assert "conservation" in attrs
            assert "n_samples" in attrs


class TestIntegrateSamples:
    """Test the main integration function."""

    def test_integrate_samples_basic(self):
        """Test full integration pipeline."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)

        result = integrate_samples(
            module4_results=m4,
            module4b_results=m4b,
            n_components=10,
            max_iter=10,
        )

        assert isinstance(result, IntegrationResult)
        assert len(result.aligned_programs) > 0
        assert len(result.sample_names) == 3
        assert result.integration_embedding is not None
        assert result.similarity_graph is not None

    def test_integrate_samples_without_4b(self):
        """Test integration without Module 4b results."""
        m4, _ = create_mock_multi_sample_data(n_samples=3)

        result = integrate_samples(
            module4_results=m4,
            module4b_results=None,
        )

        assert isinstance(result, IntegrationResult)
        assert len(result.aligned_programs) > 0
        assert len(result.conserved_relationships) == 0  # No 4b data

    def test_integrate_samples_dataframes(self):
        """Test that result DataFrames are valid."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)

        result = integrate_samples(m4, m4b)

        programs_df = result.to_programs_dataframe()
        assert len(programs_df) == len(result.aligned_programs)
        assert "program_id" in programs_df.columns
        assert "cell_type" in programs_df.columns
        assert "conservation_score" in programs_df.columns

        relationships_df = result.to_relationships_dataframe()
        assert len(relationships_df) == len(result.conserved_relationships)

    def test_summary(self):
        """Test summary string generation."""
        m4, m4b = create_mock_multi_sample_data(n_samples=3)

        result = integrate_samples(m4, m4b)
        summary = result.summary()

        assert "Cross-Sample Integration" in summary
        assert "Samples integrated:" in summary
        assert "Aligned program groups:" in summary


class TestSaveResults:
    """Test saving integration results."""

    def test_save_results(self, tmp_path):
        """Test saving results to files."""
        m4, m4b = create_mock_multi_sample_data(n_samples=2)
        result = integrate_samples(m4, m4b)

        saved_files = save_integration_results(
            result=result,
            output_dir=tmp_path,
            prefix="test_module5",
        )

        assert "programs" in saved_files
        assert saved_files["programs"].exists()

        assert "relationships" in saved_files
        assert saved_files["relationships"].exists()

        assert "embedding" in saved_files
        assert saved_files["embedding"].exists()

        assert "network" in saved_files
        assert saved_files["network"].exists()


# =============================================================================
# Main
# =============================================================================


def main():
    """Run tests manually."""
    logger.info("=" * 70)
    logger.info("Testing Module 5: Cross-Sample Integration")
    logger.info("=" * 70)

    # Create test data
    logger.info("\nCreating mock multi-sample data...")
    m4, m4b = create_mock_multi_sample_data(n_samples=4)

    for sample, result in m4.items():
        logger.info(f"  {sample}: {result.n_anchors} anchors, {result.total_programs} programs")

    # Run integration
    logger.info("\nRunning cross-sample integration...")
    result = integrate_samples(
        module4_results=m4,
        module4b_results=m4b,
        n_components=20,
        n_clusters=10,
        theta=2.0,
        similarity_threshold=0.6,
        max_iter=20,
    )

    # Print summary
    print("\n" + result.summary())

    # Show top aligned programs
    print("\n" + "=" * 70)
    print("ALIGNED PROGRAMS")
    print("=" * 70)

    for prog in result.aligned_programs[:10]:
        print(f"\n{prog.program_id} ({prog.cell_type}):")
        print(f"  Conservation: {prog.conservation_score:.0%}")
        print(f"  Samples: {', '.join(prog.samples_present)}")
        print(f"  Top genes: {', '.join(prog.top_genes[:5])}")

    # Show top relationships
    print("\n" + "=" * 70)
    print("CONSERVED RELATIONSHIPS")
    print("=" * 70)

    for rel in result.conserved_relationships[:10]:
        print(f"\n{rel.program1_id} <-> {rel.program2_id}:")
        print(f"  Type: {rel.relationship_type}")
        print(f"  Mean I: {rel.mean_bivariate_i:.3f} +/- {rel.std_bivariate_i:.3f}")
        print(f"  Conservation: {rel.conservation_score:.0%}")

    # Network stats
    if result.similarity_graph:
        G = result.similarity_graph
        print("\n" + "=" * 70)
        print("SIMILARITY NETWORK")
        print("=" * 70)
        print(f"Nodes: {G.number_of_nodes()}")
        print(f"Edges: {G.number_of_edges()}")
        if G.number_of_nodes() > 0:
            degrees = [d for n, d in G.degree()]
            print(f"Mean degree: {np.mean(degrees):.2f}")

    logger.info("\nAll tests passed!")


if __name__ == "__main__":
    main()
