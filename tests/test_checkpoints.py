"""
Unit tests for CITEgeist checkpoint management module.

Tests checkpoint saving, loading, and recovery functionality.
"""

import os
import sys
from pathlib import Path

import pytest
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from CITEgeist.model.checkpoints import CheckpointManager


# ==================== Fixtures ====================

@pytest.fixture
def checkpoint_params():
    """Standard parameters for checkpoint testing."""
    return {
        'N': 50,  # Number of spots
        'T': 9,   # Number of cell types
        'M': 100  # Number of genes
    }


@pytest.fixture
def sample_profiles(checkpoint_params):
    """Generate sample spotwise profiles."""
    np.random.seed(42)
    N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']
    profiles = {}
    for i in range(N // 2):  # Only half completed
        profiles[i] = np.random.rand(T, M) * 100
    return profiles


# ==================== CheckpointManager Tests ====================

@pytest.mark.unit
class TestCheckpointManagerInit:
    """Tests for CheckpointManager initialization."""

    def test_init_creates_directory(self, temp_output_dir, sample_name):
        """Test that initialization creates output directory."""
        checkpoint_dir = os.path.join(temp_output_dir, "checkpoints")
        manager = CheckpointManager(checkpoint_dir, sample_name)

        assert manager.output_dir.exists()
        assert manager.sample_name == sample_name

    def test_init_nested_directory(self, temp_output_dir, sample_name):
        """Test initialization with nested directory path."""
        checkpoint_dir = os.path.join(temp_output_dir, "a", "b", "c", "checkpoints")
        manager = CheckpointManager(checkpoint_dir, sample_name)

        assert manager.output_dir.exists()

    def test_init_existing_directory(self, temp_output_dir, sample_name):
        """Test initialization with existing directory."""
        checkpoint_dir = os.path.join(temp_output_dir, "checkpoints")
        os.makedirs(checkpoint_dir, exist_ok=True)

        # Should not raise error
        manager = CheckpointManager(checkpoint_dir, sample_name)
        assert manager.output_dir.exists()


@pytest.mark.unit
class TestCheckpointSaving:
    """Tests for checkpoint saving functionality."""

    def test_save_checkpoint_basic(self, temp_output_dir, sample_name, checkpoint_params, sample_profiles):
        """Test basic checkpoint saving."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        completed_spots = set(sample_profiles.keys())
        manager.save_checkpoint(completed_spots, sample_profiles, N, T, M)

        # Check that checkpoint file was created
        checkpoints = list(manager.output_dir.glob(f"{sample_name}_gene_expression_checkpoint_*.npz"))
        assert len(checkpoints) > 0

    def test_save_multiple_checkpoints(self, temp_output_dir, sample_name, checkpoint_params):
        """Test saving multiple checkpoints."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Save multiple checkpoints
        for i in range(3):
            profiles = {j: np.random.rand(T, M) for j in range(i * 10, (i + 1) * 10)}
            completed_spots = set(profiles.keys())
            manager.save_checkpoint(completed_spots, profiles, N, T, M)

        # Check that multiple checkpoint files exist
        checkpoints = list(manager.output_dir.glob(f"{sample_name}_gene_expression_checkpoint_*.npz"))
        assert len(checkpoints) == 3

    def test_save_checkpoint_empty(self, temp_output_dir, sample_name, checkpoint_params):
        """Test saving checkpoint with no completed spots."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        completed_spots = set()
        profiles = {}
        manager.save_checkpoint(completed_spots, profiles, N, T, M)

        # Should still create checkpoint
        checkpoints = list(manager.output_dir.glob(f"{sample_name}_gene_expression_checkpoint_*.npz"))
        assert len(checkpoints) > 0

    def test_save_checkpoint_data_integrity(self, temp_output_dir, sample_name, checkpoint_params, sample_profiles):
        """Test that saved data can be loaded correctly."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        completed_spots = set(sample_profiles.keys())
        manager.save_checkpoint(completed_spots, sample_profiles, N, T, M)

        # Load the checkpoint manually
        checkpoints = list(manager.output_dir.glob(f"{sample_name}_gene_expression_checkpoint_*.npz"))
        assert len(checkpoints) > 0

        data = np.load(checkpoints[0])
        assert 'profiles' in data
        assert 'completed_spots' in data
        assert data['profiles'].shape == (N, T, M)


@pytest.mark.unit
class TestCheckpointLoading:
    """Tests for checkpoint loading functionality."""

    def test_load_latest_checkpoint_basic(self, temp_output_dir, sample_name, checkpoint_params, sample_profiles):
        """Test loading the latest checkpoint."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Save checkpoint
        completed_spots = set(sample_profiles.keys())
        manager.save_checkpoint(completed_spots, sample_profiles, N, T, M)

        # Load checkpoint
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert len(loaded_spots) == len(completed_spots)
        assert len(loaded_profiles) == len(sample_profiles)

    def test_load_latest_of_multiple(self, temp_output_dir, sample_name, checkpoint_params):
        """Test that load_latest_checkpoint loads the most recent checkpoint."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Save multiple checkpoints
        for i in range(3):
            profiles = {j: np.random.rand(T, M) * (i + 1) for j in range(10)}
            completed_spots = set(profiles.keys())
            manager.save_checkpoint(completed_spots, profiles, N, T, M)

        # Load should get the latest (i=2)
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert len(loaded_spots) == 10
        # Latest profiles should have higher values due to multiplication
        mean_value = np.mean([prof.mean() for prof in loaded_profiles.values()])
        assert mean_value > 1.0  # Should be from latest iteration

    def test_load_no_checkpoints(self, temp_output_dir, sample_name, checkpoint_params):
        """Test loading when no checkpoints exist."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert loaded_spots == set()
        assert loaded_profiles == {}

    def test_load_checkpoint_filters_nan(self, temp_output_dir, sample_name, checkpoint_params):
        """Test that loading filters out profiles with NaN values."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Create profiles with some NaN values
        profiles = {}
        for i in range(20):
            prof = np.random.rand(T, M)
            if i % 5 == 0:  # Every 5th profile has NaN
                prof[0, 0] = np.nan
            profiles[i] = prof

        completed_spots = set(profiles.keys())
        manager.save_checkpoint(completed_spots, profiles, N, T, M)

        # Load checkpoint
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        # Should have filtered out NaN profiles
        assert len(loaded_profiles) < len(profiles)
        for prof in loaded_profiles.values():
            assert not np.any(np.isnan(prof))

    def test_load_checkpoint_wrong_dimensions(self, temp_output_dir, sample_name):
        """Test loading checkpoint with mismatched dimensions."""
        manager = CheckpointManager(temp_output_dir, sample_name)

        # Save with one set of dimensions
        profiles = {i: np.random.rand(9, 100) for i in range(10)}
        completed_spots = set(profiles.keys())
        manager.save_checkpoint(completed_spots, profiles, 50, 9, 100)

        # Try to load with different dimensions
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(50, 9, 200)  # Different M

        # Should return empty due to dimension mismatch
        assert loaded_spots == set()
        assert loaded_profiles == {}


@pytest.mark.unit
class TestCompleteRun:
    """Tests for complete run checking."""

    def test_save_and_check_complete_run(self, temp_output_dir, sample_name, checkpoint_params, sample_profiles):
        """Test saving and checking complete run."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Create complete profiles for all spots
        complete_profiles = {i: np.random.rand(T, M) for i in range(N)}
        completed_spots = set(complete_profiles.keys())

        manager.save_final_results(complete_profiles, completed_spots, N, T, M)

        # Check for complete run
        loaded_profiles = manager.check_complete_run(N, T, M)

        assert loaded_profiles is not None
        assert len(loaded_profiles) == N

    def test_check_complete_run_not_exists(self, temp_output_dir, sample_name, checkpoint_params):
        """Test checking for complete run when it doesn't exist."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        result = manager.check_complete_run(N, T, M)

        assert result is None

    def test_check_complete_run_wrong_dimensions(self, temp_output_dir, sample_name):
        """Test checking complete run with wrong dimensions."""
        manager = CheckpointManager(temp_output_dir, sample_name)

        # Save with one set of dimensions
        N, T, M = 50, 9, 100
        profiles = {i: np.random.rand(T, M) for i in range(N)}
        completed_spots = set(profiles.keys())
        manager.save_final_results(profiles, completed_spots, N, T, M)

        # Check with different dimensions
        result = manager.check_complete_run(50, 9, 200)  # Different M

        assert result is None


@pytest.mark.unit
class TestCleanup:
    """Tests for checkpoint cleanup functionality."""

    def test_cleanup_old_checkpoints(self, temp_output_dir, sample_name, checkpoint_params):
        """Test cleanup of old checkpoint files."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Create multiple checkpoints
        for i in range(5):
            profiles = {j: np.random.rand(T, M) for j in range((i * 10), (i * 10) + 10)}
            completed_spots = set(profiles.keys())
            manager.save_checkpoint(completed_spots, profiles, N, T, M)

        # The _cleanup_old_checkpoints is called automatically on each save
        # So we should only have one checkpoint (the latest)
        checkpoints = list(manager.output_dir.glob(f"{sample_name}_gene_expression_checkpoint_*.npz"))
        assert len(checkpoints) == 1


@pytest.mark.unit
class TestCheckpointRecovery:
    """Tests for checkpoint recovery and corruption handling."""

    def test_recovery_from_partial_checkpoint(self, temp_output_dir, sample_name, checkpoint_params):
        """Test recovery from partial checkpoint."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Save partial checkpoint
        profiles = {i: np.random.rand(T, M) for i in range(N // 2)}
        completed_spots = set(profiles.keys())
        manager.save_checkpoint(completed_spots, profiles, N, T, M)

        # Load and verify we can continue
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert len(loaded_spots) == N // 2
        assert len(loaded_profiles) == N // 2

        # Verify we can add more spots
        remaining_spots = set(range(N // 2, N))
        assert loaded_spots.isdisjoint(remaining_spots)

    def test_corrupted_checkpoint_handling(self, temp_output_dir, sample_name, checkpoint_params):
        """Test handling of corrupted checkpoint files."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = checkpoint_params['N'], checkpoint_params['T'], checkpoint_params['M']

        # Create a corrupted checkpoint file
        corrupted_file = manager.output_dir / f"{sample_name}_gene_expression_checkpoint_1.npz"
        with open(corrupted_file, 'w') as f:
            f.write("corrupted data")

        # Loading should handle corruption gracefully
        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert loaded_spots == set()
        assert loaded_profiles == {}


@pytest.mark.unit
class TestEdgeCases:
    """Tests for edge cases in checkpoint management."""

    def test_checkpoint_with_zero_spots(self, temp_output_dir, sample_name):
        """Test checkpoint with zero spots."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = 0, 9, 100

        profiles = {}
        completed_spots = set()
        manager.save_checkpoint(completed_spots, profiles, N, T, M)

        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert loaded_spots == set()
        assert loaded_profiles == {}

    def test_checkpoint_with_single_spot(self, temp_output_dir, sample_name):
        """Test checkpoint with single spot."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = 1, 9, 100

        profiles = {0: np.random.rand(T, M)}
        completed_spots = {0}
        manager.save_checkpoint(completed_spots, profiles, N, T, M)

        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert len(loaded_spots) == 1
        assert 0 in loaded_profiles

    def test_checkpoint_with_large_data(self, temp_output_dir, sample_name):
        """Test checkpoint with large data dimensions."""
        manager = CheckpointManager(temp_output_dir, sample_name)
        N, T, M = 5000, 20, 2000

        # Create a few profiles (not all to save memory)
        profiles = {i: np.random.rand(T, M) for i in range(10)}
        completed_spots = set(profiles.keys())

        # Should handle large dimensions
        manager.save_checkpoint(completed_spots, profiles, N, T, M)

        loaded_spots, loaded_profiles = manager.load_latest_checkpoint(N, T, M)

        assert len(loaded_spots) == 10

    def test_different_sample_names(self, temp_output_dir):
        """Test that different sample names create separate checkpoints."""
        manager1 = CheckpointManager(temp_output_dir, "sample1")
        manager2 = CheckpointManager(temp_output_dir, "sample2")

        N, T, M = 50, 9, 100

        # Save checkpoints for different samples
        profiles1 = {i: np.random.rand(T, M) * 1 for i in range(10)}
        profiles2 = {i: np.random.rand(T, M) * 2 for i in range(10)}

        manager1.save_checkpoint(set(profiles1.keys()), profiles1, N, T, M)
        manager2.save_checkpoint(set(profiles2.keys()), profiles2, N, T, M)

        # Load should get correct checkpoint for each sample
        loaded1_spots, loaded1_profiles = manager1.load_latest_checkpoint(N, T, M)
        loaded2_spots, loaded2_profiles = manager2.load_latest_checkpoint(N, T, M)

        assert len(loaded1_spots) == 10
        assert len(loaded2_spots) == 10

        # Verify values are different
        mean1 = np.mean([p.mean() for p in loaded1_profiles.values()])
        mean2 = np.mean([p.mean() for p in loaded2_profiles.values()])
        assert abs(mean1 - mean2) > 0.1  # Should be significantly different
