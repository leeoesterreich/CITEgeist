"""
Integration tests for CITEgeist full pipeline.

These tests verify that components work together correctly.
Marked as integration and may be slow.
"""

import os
import sys

import pytest
import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.checkpoints import CheckpointManager


# ==================== Full Pipeline Integration Tests ====================

@pytest.mark.integration
class TestFullPipeline:
    """Integration tests for complete CITEgeist pipeline."""

    def test_model_initialization_and_setup(self, temp_output_dir, sample_name,
                                           mock_gex_adata, mock_protein_adata,
                                           mock_cell_profile_dict):
        """Test full model initialization and setup workflow."""
        # Initialize model
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        # Load cell profiles
        model.load_cell_profile_dict(mock_cell_profile_dict)

        # Verify state
        assert model.cell_profile_dict is not None
        assert model.gene_expression_adata is not None
        assert model.antibody_capture_adata is not None

    def test_data_split_workflow(self, temp_output_dir, sample_name, mock_combined_adata,
                                 mock_cell_profile_dict):
        """Test workflow: initialize -> split -> load profiles."""
        # Initialize with combined data
        model = CitegeistModel(
            sample_name=sample_name,
            adata=mock_combined_adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        # Split data
        model.split_adata()

        # Load profiles
        model.load_cell_profile_dict(mock_cell_profile_dict)

        # Verify components are ready
        assert model.gene_expression_adata is not None
        assert model.antibody_capture_adata is not None
        assert model.cell_profile_dict is not None

    def test_preprocessing_workflow(self, temp_output_dir, sample_name,
                                    mock_gex_adata, mock_protein_adata,
                                    mock_cell_profile_dict):
        """Test preprocessing workflow (without optimization)."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        model.load_cell_profile_dict(mock_cell_profile_dict)

        # Note: Actual preprocessing methods would be called here
        # But we're testing the setup and state management

        assert model.preprocessed_gex is False
        assert model.preprocessed_antibody is False


# ==================== Checkpoint Integration Tests ====================

@pytest.mark.integration
class TestCheckpointIntegration:
    """Integration tests for checkpoint system with model."""

    def test_checkpoint_save_and_load_workflow(self, temp_output_dir, sample_name):
        """Test complete checkpoint save/load workflow."""
        # Setup checkpoint manager
        checkpoint_dir = os.path.join(temp_output_dir, "checkpoints")
        manager = CheckpointManager(checkpoint_dir, sample_name)

        N, T, M = 100, 9, 200

        # Simulate partial completion
        profiles = {i: np.random.rand(T, M) for i in range(50)}
        completed_spots = set(profiles.keys())

        # Save checkpoint
        manager.save_checkpoint(completed_spots, profiles, N, T, M)

        # Load checkpoint
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        # Verify recovery
        assert len(loaded_spots) == 50
        assert len(loaded_profiles) == 50

        # Simulate completion of remaining spots
        for i in range(50, N):
            profiles[i] = np.random.rand(T, M)
        completed_spots = set(profiles.keys())

        # Save final results
        manager.save_final_results(profiles, completed_spots, N, T, M)

        # Verify complete run exists
        complete_profiles = manager.check_complete_run(N, T, M)
        assert complete_profiles is not None
        assert len(complete_profiles) == N

    def test_checkpoint_recovery_after_interruption(self, temp_output_dir, sample_name):
        """Test recovery from interruption using checkpoints."""
        checkpoint_dir = os.path.join(temp_output_dir, "checkpoints")
        N, T, M = 100, 9, 200

        # First run: partial completion
        manager1 = CheckpointManager(checkpoint_dir, sample_name)
        profiles_partial = {i: np.random.rand(T, M) for i in range(30)}
        manager1.save_checkpoint(set(profiles_partial.keys()), profiles_partial, N, T, M)

        # Simulate restart: new manager instance
        manager2 = CheckpointManager(checkpoint_dir, sample_name)

        # Check for existing progress
        loaded_spots, loaded_profiles = manager2.load_latest_checkpoint(N, T, M)

        # Should recover partial work
        assert len(loaded_spots) == 30
        assert len(loaded_profiles) == 30

        # Continue from checkpoint
        remaining_spots = set(range(N)) - loaded_spots
        assert len(remaining_spots) == 70


# ==================== Data Flow Integration Tests ====================

@pytest.mark.integration
class TestDataFlow:
    """Integration tests for data flow through pipeline."""

    def test_data_consistency_through_split(self, temp_output_dir, sample_name,
                                           n_genes, n_proteins, n_spots):
        """Test that data remains consistent through splitting."""
        # Create combined dataset
        import scanpy as sc
        from anndata import AnnData
        import scipy.sparse as sp

        gene_data = sp.random(n_spots, n_genes, density=0.8, format='csr')
        protein_data = np.random.rand(n_spots, n_proteins)

        gene_names = [f"GENE_{i}" for i in range(n_genes)]
        protein_names = [f"PROT_{i}" for i in range(n_proteins)]
        spot_names = [f"spot_{i:04d}" for i in range(n_spots)]

        combined_data = np.hstack([gene_data.toarray(), protein_data])

        adata = AnnData(X=combined_data)
        adata.obs_names = spot_names
        adata.var_names = gene_names + protein_names
        adata.var['feature_types'] = ['Gene Expression'] * n_genes + ['Antibody Capture'] * n_proteins
        adata.obsm['spatial'] = np.random.rand(n_spots, 2) * 100

        # Initialize and split
        model = CitegeistModel(
            sample_name=sample_name,
            adata=adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        model.split_adata()

        # Verify data integrity
        assert model.gene_expression_adata.n_obs == n_spots
        assert model.antibody_capture_adata.n_obs == n_spots
        assert model.gene_expression_adata.n_vars == n_genes
        assert model.antibody_capture_adata.n_vars == n_proteins

        # Verify spot names are preserved
        assert list(model.gene_expression_adata.obs_names) == spot_names
        assert list(model.antibody_capture_adata.obs_names) == spot_names

    def test_spatial_coordinates_preserved(self, temp_output_dir, sample_name,
                                          mock_combined_adata):
        """Test that spatial coordinates are preserved through operations."""
        original_coords = mock_combined_adata.obsm['spatial'].copy()

        model = CitegeistModel(
            sample_name=sample_name,
            adata=mock_combined_adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        model.split_adata()

        # Verify spatial coordinates preserved
        assert np.array_equal(
            model.gene_expression_adata.obsm['spatial'],
            original_coords
        )
        assert np.array_equal(
            model.antibody_capture_adata.obsm['spatial'],
            original_coords
        )


# ==================== Utility Integration Tests ====================

@pytest.mark.integration
class TestUtilityIntegration:
    """Integration tests for utility functions in context."""

    def test_normalization_pipeline(self):
        """Test sequence of normalization operations."""
        from CITEgeist.model.citegeist_model import CitegeistModel

        # Create test data
        matrix = np.random.rand(50, 100) * 1000

        # Apply sequence of operations
        winsorized = CitegeistModel.winsorize(matrix, lower_percentile=5, upper_percentile=95)
        normalized = CitegeistModel.row_normalize(winsorized, target_sum=10000)
        clr_normalized = CitegeistModel.global_clr(normalized, epsilon=1e-6)

        # Verify final output is valid
        assert clr_normalized.shape == matrix.shape
        assert np.all(np.isfinite(clr_normalized))
        assert np.all(clr_normalized > 0)

    def test_benchmarking_pipeline(self, cell_type_names, n_spots):
        """Test benchmarking with generated predictions."""
        from CITEgeist.model.utils import benchmark_cell_proportions

        # Generate ground truth and predictions
        np.random.seed(42)
        true_props = np.random.dirichlet(np.ones(len(cell_type_names)), size=n_spots)

        # Generate predictions with some noise
        noise = np.random.randn(n_spots, len(cell_type_names)) * 0.05
        pred_props = true_props + noise

        # Ensure proportions are valid
        pred_props = np.maximum(pred_props, 0)
        pred_props = pred_props / pred_props.sum(axis=1, keepdims=True)

        # Run benchmarking
        results = benchmark_cell_proportions(true_props, pred_props, cell_type_names)

        # Verify results structure
        assert 'JSD' in results
        assert 'RMSE' in results
        assert 'MAE' in results
        assert 'Sum_RMSE' in results
        assert 'Sum_MAE' in results
        assert 'corr' in results

        # Verify reasonable metrics
        assert 0 <= results['JSD'] <= 1
        assert results['Sum_RMSE'] >= 0
        assert results['Sum_MAE'] >= 0


# ==================== Error Handling Integration Tests ====================

@pytest.mark.integration
class TestErrorHandling:
    """Integration tests for error handling across components."""

    def test_invalid_workflow_sequence(self, temp_output_dir, sample_name,
                                      mock_gex_adata, mock_protein_adata):
        """Test that invalid workflow sequences are caught."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        # Try to use model without loading cell profiles
        # (This would fail in actual optimization steps)
        assert model.cell_profile_dict is None

    def test_checkpoint_corruption_recovery(self, temp_output_dir, sample_name):
        """Test recovery from corrupted checkpoints."""
        checkpoint_dir = os.path.join(temp_output_dir, "checkpoints")
        manager = CheckpointManager(checkpoint_dir, sample_name)

        N, T, M = 50, 9, 100

        # Save valid checkpoint
        profiles = {i: np.random.rand(T, M) for i in range(25)}
        manager.save_checkpoint(set(profiles.keys()), profiles, N, T, M)

        # Corrupt the checkpoint file
        checkpoints = list(manager.output_dir.glob(f"{sample_name}_gene_expression_checkpoint_*.npz"))
        with open(checkpoints[0], 'w') as f:
            f.write("corrupted")

        # Should handle corruption gracefully
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert loaded_spots == set()
        assert loaded_profiles == {}

    def test_missing_spatial_coordinates(self, temp_output_dir, sample_name):
        """Test handling of data without spatial coordinates."""
        from anndata import AnnData

        # Create AnnData without spatial coordinates
        adata = AnnData(X=np.random.rand(10, 20))
        adata.var['feature_types'] = ['Gene Expression'] * 10 + ['Antibody Capture'] * 10

        model = CitegeistModel(
            sample_name=sample_name,
            adata=adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        model.split_adata()

        # Spatial operations would fail, but model should initialize
        assert model.gene_expression_adata is not None


# ==================== Performance Integration Tests ====================

@pytest.mark.integration
@pytest.mark.slow
class TestPerformance:
    """Integration tests for performance with larger datasets."""

    def test_large_dataset_handling(self, temp_output_dir, sample_name):
        """Test handling of larger dataset (but still manageable for tests)."""
        from anndata import AnnData
        import scipy.sparse as sp

        # Create moderately large dataset
        n_spots = 1000
        n_genes = 500
        n_proteins = 18

        gene_data = sp.random(n_spots, n_genes, density=0.3, format='csr')
        protein_data = np.random.rand(n_spots, n_proteins)

        gex_adata = AnnData(X=gene_data)
        gex_adata.var['feature_types'] = 'Gene Expression'
        gex_adata.obsm['spatial'] = np.random.rand(n_spots, 2) * 100

        protein_adata = AnnData(X=protein_data)
        protein_adata.var['feature_types'] = 'Antibody Capture'
        protein_adata.obsm['spatial'] = np.random.rand(n_spots, 2) * 100

        # Initialize model
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=gex_adata,
            antibody_capture_adata=protein_adata
        )

        # Should handle large dataset without errors
        assert model.gene_expression_adata.n_obs == n_spots
        assert model.gene_expression_adata.n_vars == n_genes

    def test_checkpoint_performance_many_saves(self, temp_output_dir, sample_name):
        """Test checkpoint system with many save operations."""
        checkpoint_dir = os.path.join(temp_output_dir, "checkpoints")
        manager = CheckpointManager(checkpoint_dir, sample_name)

        N, T, M = 200, 9, 100

        # Simulate incremental progress with many checkpoints
        for i in range(0, 100, 10):
            profiles = {j: np.random.rand(T, M) for j in range(i, i + 10)}
            manager.save_checkpoint(set(profiles.keys()), profiles, N, T, M)

        # Should maintain only latest checkpoint
        checkpoints = list(manager.output_dir.glob(f"{sample_name}_gene_expression_checkpoint_*.npz"))
        assert len(checkpoints) == 1

        # Should be able to load latest
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)
        assert len(loaded_spots) > 0
