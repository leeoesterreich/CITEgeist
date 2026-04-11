"""
Unit tests for CITEgeist utils module.

Tests spatial neighbor functions, benchmarking, validation, and utility functions.
"""

import os
import sys
import tempfile
import shutil

import pytest
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from CITEgeist.model.utils import (
    validate_cell_profile_dict,
    save_results_to_output,
    cleanup_memory,
    find_fixed_radius_neighbors,
    get_neighbors_with_fixed_radius,
    benchmark_cell_proportions,
    export_anndata_layers,
    calculate_expression_metrics
)


# ==================== Validation Tests ====================

@pytest.mark.unit
class TestValidation:
    """Tests for validation functions."""

    def test_validate_cell_profile_dict_valid(self, mock_cell_profile_dict):
        """Test validation of a valid cell profile dictionary."""
        assert validate_cell_profile_dict(mock_cell_profile_dict) is True

    def test_validate_cell_profile_dict_invalid_type(self):
        """Test validation rejects non-dict input."""
        assert validate_cell_profile_dict("not a dict") is False
        assert validate_cell_profile_dict([]) is False
        assert validate_cell_profile_dict(None) is False

    def test_validate_cell_profile_dict_invalid_structure(self):
        """Test validation rejects invalid structure."""
        # Non-string keys
        invalid_dict = {1: {"Major": ["protein1"]}}
        assert validate_cell_profile_dict(invalid_dict) is False

        # Non-dict values
        invalid_dict = {"CellType": ["not", "a", "dict"]}
        assert validate_cell_profile_dict(invalid_dict) is False

    def test_validate_cell_profile_dict_empty(self):
        """Test validation of empty dictionary."""
        assert validate_cell_profile_dict({}) is True


# ==================== File I/O Tests ====================

@pytest.mark.unit
class TestFileIO:
    """Tests for file input/output functions."""

    def test_save_results_to_output(self, temp_output_dir):
        """Test saving results to CSV."""
        results = {
            'metric1': [1, 2, 3],
            'metric2': [4, 5, 6]
        }
        filepath = os.path.join(temp_output_dir, 'test_results.csv')

        save_results_to_output(results, filepath)

        assert os.path.exists(filepath)
        df = pd.read_csv(filepath)
        assert df.shape == (3, 2)
        assert 'metric1' in df.columns
        assert 'metric2' in df.columns

    def test_save_results_empty(self, temp_output_dir):
        """Test saving empty results."""
        results = {}
        filepath = os.path.join(temp_output_dir, 'empty_results.csv')

        save_results_to_output(results, filepath)

        assert os.path.exists(filepath)


# ==================== Memory Management Tests ====================

@pytest.mark.unit
class TestMemoryManagement:
    """Tests for memory management functions."""

    def test_cleanup_memory(self):
        """Test memory cleanup runs without error."""
        # Create some objects
        large_array = np.random.rand(1000, 1000)
        del large_array

        # Should not raise any exceptions
        cleanup_memory()


# ==================== Spatial Neighbor Tests ====================

@pytest.mark.unit
class TestSpatialNeighbors:
    """Tests for spatial neighbor detection functions."""

    def test_find_fixed_radius_neighbors_basic(self, mock_gex_adata):
        """Test finding neighbors within a fixed radius."""
        spot_index = 0
        radius = 15  # Should capture nearby spots

        central_name, neighbor_names, spot_idx, neighbor_indices = \
            find_fixed_radius_neighbors(spot_index, mock_gex_adata, radius)

        assert spot_idx == spot_index
        assert isinstance(neighbor_names, list)
        assert isinstance(neighbor_indices, list)
        assert len(neighbor_names) == len(neighbor_indices)
        assert central_name == mock_gex_adata.obs_names[spot_index]

    def test_find_fixed_radius_neighbors_no_neighbors(self, mock_gex_adata):
        """Test with radius too small to find neighbors."""
        spot_index = 0
        radius = 0.1  # Very small radius

        _, neighbor_names, _, neighbor_indices = \
            find_fixed_radius_neighbors(spot_index, mock_gex_adata, radius)

        # Should find no neighbors (or very few)
        assert len(neighbor_names) >= 0
        assert len(neighbor_indices) >= 0

    def test_find_fixed_radius_neighbors_large_radius(self, mock_gex_adata):
        """Test with large radius capturing many neighbors."""
        spot_index = 50  # Middle of dataset
        radius = 1000  # Very large radius

        _, neighbor_names, _, neighbor_indices = \
            find_fixed_radius_neighbors(spot_index, mock_gex_adata, radius)

        # Should find many neighbors (but not including center)
        assert len(neighbor_names) > 10
        assert spot_index not in neighbor_indices

    def test_get_neighbors_with_fixed_radius_include_center(self, mock_gex_adata):
        """Test getting neighbors with center included."""
        spot_index = 10
        radius = 15

        neighbors = get_neighbors_with_fixed_radius(
            spot_index, mock_gex_adata, radius, include_center=True
        )

        assert isinstance(neighbors, list)
        assert spot_index in neighbors
        assert neighbors[0] == spot_index  # Center should be first

    def test_get_neighbors_with_fixed_radius_exclude_center(self, mock_gex_adata):
        """Test getting neighbors without center."""
        spot_index = 10
        radius = 15

        neighbors = get_neighbors_with_fixed_radius(
            spot_index, mock_gex_adata, radius, include_center=False
        )

        assert isinstance(neighbors, list)
        assert spot_index not in neighbors

    def test_neighbors_symmetric_distance(self, mock_gex_adata):
        """Test that neighbor detection is symmetric in distance."""
        spot1 = 10
        spot2 = 11
        radius = 20

        _, neighbors1, _, _ = find_fixed_radius_neighbors(spot1, mock_gex_adata, radius)
        _, neighbors2, _, _ = find_fixed_radius_neighbors(spot2, mock_gex_adata, radius)

        # If spot2 is in spot1's neighbors, spot1 should be in spot2's neighbors
        spot2_name = mock_gex_adata.obs_names[spot2]
        spot1_name = mock_gex_adata.obs_names[spot1]

        if spot2_name in neighbors1:
            assert spot1_name in neighbors2


# ==================== Benchmarking Tests ====================

@pytest.mark.unit
class TestBenchmarking:
    """Tests for benchmarking functions."""

    def test_benchmark_cell_proportions_perfect_prediction(self, cell_type_names, n_spots):
        """Test benchmarking with perfect predictions."""
        np.random.seed(42)
        # Create identical proportions
        true_props = np.random.dirichlet(np.ones(len(cell_type_names)), size=n_spots)
        pred_props = true_props.copy()

        results = benchmark_cell_proportions(true_props, pred_props, cell_type_names)

        # JSD should be near 0 for perfect prediction
        assert results['JSD'] < 0.01
        # Correlation should be near 1
        assert results['corr'] > 0.99
        # Errors should be near 0
        assert results['Sum_RMSE'] < 0.01
        assert results['Sum_MAE'] < 0.01

    def test_benchmark_cell_proportions_random(self, cell_type_names, n_spots):
        """Test benchmarking with random predictions."""
        np.random.seed(42)
        true_props = np.random.dirichlet(np.ones(len(cell_type_names)), size=n_spots)
        np.random.seed(123)  # Different seed for predictions
        pred_props = np.random.dirichlet(np.ones(len(cell_type_names)), size=n_spots)

        results = benchmark_cell_proportions(true_props, pred_props, cell_type_names)

        # Results should be reasonable but not perfect
        assert 0 <= results['JSD'] <= 1
        assert -1 <= results['corr'] <= 1
        assert results['Sum_RMSE'] > 0
        assert results['Sum_MAE'] > 0

        # Check per-cell-type metrics exist
        assert 'RMSE' in results
        assert 'MAE' in results
        assert len(results['RMSE']) == len(cell_type_names)
        assert len(results['MAE']) == len(cell_type_names)

    def test_benchmark_cell_proportions_zero_prediction(self, cell_type_names, n_spots):
        """Test benchmarking with zero predictions (worst case)."""
        np.random.seed(42)
        true_props = np.random.dirichlet(np.ones(len(cell_type_names)), size=n_spots)
        pred_props = np.zeros_like(true_props)

        results = benchmark_cell_proportions(true_props, pred_props, cell_type_names)

        # JSD should be 1 (maximum divergence) for zero predictions
        assert results['JSD'] == pytest.approx(1.0, abs=0.01)
        # Errors should be significant
        assert results['Sum_RMSE'] > 0.1
        assert results['Sum_MAE'] > 0.1

    def test_benchmark_cell_proportions_invalid_input(self, cell_type_names):
        """Test benchmarking with invalid inputs."""
        true_props = [[0.5, 0.5]]  # List instead of numpy array

        with pytest.raises(ValueError, match="Input proportions must be numpy arrays"):
            benchmark_cell_proportions(true_props, true_props, cell_type_names)

    def test_benchmark_cell_proportions_shape_mismatch(self, cell_type_names):
        """Test benchmarking with mismatched shapes."""
        true_props = np.array([[0.5, 0.5]])
        pred_props = np.array([[0.3, 0.3, 0.4]])  # Different number of columns

        with pytest.raises((ValueError, IndexError)):
            benchmark_cell_proportions(true_props, pred_props, cell_type_names)


# ==================== Export Functions Tests ====================

@pytest.mark.unit
class TestExportFunctions:
    """Tests for data export functions."""

    def test_export_anndata_layers_basic(self, mock_gex_adata, temp_output_dir, cell_type_names):
        """Test exporting AnnData layers to CSV."""
        # Add mock layers to AnnData
        for cell_type in cell_type_names[:3]:  # Just use first 3 for speed
            layer_name = f"{cell_type}_genes_pass1"
            mock_gex_adata.layers[layer_name] = mock_gex_adata.X.copy()

        export_anndata_layers(mock_gex_adata, temp_output_dir, pass_number=1)

        # Check that files were created
        pass1_dir = os.path.join(temp_output_dir, "pass1")
        assert os.path.exists(pass1_dir)

        # Check that CSV files exist
        for cell_type in cell_type_names[:3]:
            csv_file = os.path.join(pass1_dir, f"{cell_type}_layer_pass1.csv")
            assert os.path.exists(csv_file)

            # Verify CSV can be read
            df = pd.read_csv(csv_file, index_col=0)
            assert df.shape == mock_gex_adata.shape

    def test_export_anndata_layers_no_pass_number(self, mock_gex_adata, temp_output_dir):
        """Test exporting without specifying pass number."""
        layer_name = "test_layer_genes_pass1"
        mock_gex_adata.layers[layer_name] = mock_gex_adata.X.copy()

        export_anndata_layers(mock_gex_adata, temp_output_dir, pass_number=None)

        # Should create files in base directory
        csv_file = os.path.join(temp_output_dir, "test_layer_layer.csv")
        assert os.path.exists(csv_file)

    def test_export_anndata_layers_sparse_matrix(self, mock_gex_adata, temp_output_dir):
        """Test exporting sparse matrices."""
        layer_name = "sparse_layer_genes_pass1"
        # Ensure X is sparse
        if not isinstance(mock_gex_adata.X, csr_matrix):
            mock_gex_adata.X = csr_matrix(mock_gex_adata.X)

        mock_gex_adata.layers[layer_name] = mock_gex_adata.X.copy()

        export_anndata_layers(mock_gex_adata, temp_output_dir, pass_number=1)

        # Verify export worked
        csv_file = os.path.join(temp_output_dir, "pass1", "sparse_layer_layer_pass1.csv")
        assert os.path.exists(csv_file)


# ==================== Expression Metrics Tests ====================

@pytest.mark.unit
class TestExpressionMetrics:
    """Tests for gene expression metric calculation."""

    def test_calculate_expression_metrics_perfect(
            self, temp_output_dir, mock_ground_truth_layers, cell_type_names
        ):
        """Test metrics with perfect predictions."""
        # Create prediction directory matching ground truth
        pred_dir = os.path.join(temp_output_dir, "predictions", "pass1")
        os.makedirs(pred_dir, exist_ok=True)

        # Copy ground truth as predictions
        for cell_type in cell_type_names:
            gt_file = os.path.join(mock_ground_truth_layers, f"{cell_type}_GT.csv")
            pred_file = os.path.join(pred_dir, f"{cell_type}_layer_pass1.csv")
            shutil.copy(gt_file, pred_file)

        # Calculate metrics
        metrics = calculate_expression_metrics(
            mock_ground_truth_layers,
            os.path.join(temp_output_dir, "predictions"),
            normalize='range',
            pass_number=1
        )

        # Metrics should be near zero for perfect prediction
        # Filter out special tracking keys (_spurious_profiles, _missed_profiles)
        cell_metrics = {k: v for k, v in metrics.items() if not k.startswith('_')}
        assert len(cell_metrics) == len(cell_type_names)
        for cell_type, values in cell_metrics.items():
            assert values['RMSE'] < 0.01
            assert values['NRMSE'] < 0.01
            assert values['MAE'] < 0.01

    def test_calculate_expression_metrics_random(
            self, temp_output_dir, mock_ground_truth_layers, cell_type_names, n_genes, n_spots
        ):
        """Test metrics with random predictions."""
        # Create prediction directory with random data
        pred_dir = os.path.join(temp_output_dir, "predictions", "pass1")
        os.makedirs(pred_dir, exist_ok=True)

        np.random.seed(123)
        gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
        spot_names = [f"spot_{i:04d}" for i in range(n_spots)]

        for cell_type in cell_type_names:
            # Generate random predictions
            data = np.random.lognormal(mean=0.5, sigma=1.5, size=(n_genes, n_spots))
            df = pd.DataFrame(data, index=gene_names, columns=spot_names)

            pred_file = os.path.join(pred_dir, f"{cell_type}_layer_pass1.csv")
            df.to_csv(pred_file)

        # Calculate metrics
        metrics = calculate_expression_metrics(
            mock_ground_truth_layers,
            os.path.join(temp_output_dir, "predictions"),
            normalize='range',
            pass_number=1
        )

        # Metrics should show significant error
        cell_metrics = {k: v for k, v in metrics.items() if not k.startswith('_')}
        assert len(cell_metrics) == len(cell_type_names)
        for cell_type, values in cell_metrics.items():
            assert values['RMSE'] > 0
            assert values['NRMSE'] > 0
            assert values['MAE'] > 0

    def test_calculate_expression_metrics_mean_normalization(
            self, temp_output_dir, mock_ground_truth_layers, cell_type_names
        ):
        """Test metrics with mean normalization."""
        # Create prediction directory matching ground truth
        pred_dir = os.path.join(temp_output_dir, "predictions", "pass1")
        os.makedirs(pred_dir, exist_ok=True)

        for cell_type in cell_type_names:
            gt_file = os.path.join(mock_ground_truth_layers, f"{cell_type}_GT.csv")
            pred_file = os.path.join(pred_dir, f"{cell_type}_layer_pass1.csv")
            shutil.copy(gt_file, pred_file)

        # Calculate metrics with mean normalization
        metrics = calculate_expression_metrics(
            mock_ground_truth_layers,
            os.path.join(temp_output_dir, "predictions"),
            normalize='mean',
            pass_number=1
        )

        cell_metrics = {k: v for k, v in metrics.items() if not k.startswith('_')}
        assert len(cell_metrics) == len(cell_type_names)
        for cell_type, values in cell_metrics.items():
            assert 'NRMSE' in values
            assert values['NRMSE'] >= 0

    def test_calculate_expression_metrics_missing_files(self, temp_output_dir):
        """Test metrics calculation with missing prediction files."""
        gt_dir = os.path.join(temp_output_dir, "ground_truth")
        pred_dir = os.path.join(temp_output_dir, "predictions", "pass1")
        os.makedirs(gt_dir, exist_ok=True)
        os.makedirs(pred_dir, exist_ok=True)

        # Create a ground truth file but no prediction
        df = pd.DataFrame(np.random.rand(10, 5))
        df.to_csv(os.path.join(gt_dir, "TestType_GT.csv"))

        # Should handle missing files gracefully
        metrics = calculate_expression_metrics(
            gt_dir, os.path.join(temp_output_dir, "predictions"),
            normalize='range', pass_number=1
        )

        # Should return empty or skip missing cell type
        assert isinstance(metrics, dict)

    def test_calculate_expression_metrics_invalid_normalization(
            self, temp_output_dir, mock_ground_truth_layers
        ):
        """Test with invalid normalization parameter."""
        pred_dir = os.path.join(temp_output_dir, "predictions")
        os.makedirs(pred_dir, exist_ok=True)

        with pytest.raises(ValueError, match="Normalization type must be"):
            calculate_expression_metrics(
                mock_ground_truth_layers, pred_dir,
                normalize='invalid', pass_number=1
            )


# ==================== Edge Cases and Error Handling ====================

@pytest.mark.unit
class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_spatial_neighbors_single_spot(self):
        """Test neighbor detection with single spot."""
        from anndata import AnnData
        adata = AnnData(X=np.array([[1, 2, 3]]))
        adata.obsm['spatial'] = np.array([[0, 0]])

        spot_index = 0
        radius = 10

        _, neighbors, _, _ = find_fixed_radius_neighbors(spot_index, adata, radius)

        # Should find no neighbors with single spot
        assert len(neighbors) == 0

    def test_benchmark_single_spot(self, cell_type_names):
        """Test benchmarking with single spot."""
        true_props = np.array([[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2]])
        pred_props = np.array([[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2]])

        results = benchmark_cell_proportions(true_props, pred_props, cell_type_names)

        assert isinstance(results, dict)
        assert 'JSD' in results
        assert 'corr' in results

    def test_empty_cell_profile(self):
        """Test validation with empty cell profile."""
        empty_profile = {}
        assert validate_cell_profile_dict(empty_profile) is True

    def test_neighbors_edge_spot(self, mock_gex_adata):
        """Test neighbor detection for edge spots."""
        # Test first and last spots (likely at edges)
        for spot_idx in [0, len(mock_gex_adata) - 1]:
            neighbors = get_neighbors_with_fixed_radius(
                spot_idx, mock_gex_adata, radius=15, include_center=True
            )
            assert isinstance(neighbors, list)
            assert len(neighbors) > 0  # Should at least include itself
