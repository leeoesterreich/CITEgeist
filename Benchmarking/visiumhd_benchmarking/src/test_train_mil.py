"""Tests for MIL training pipeline."""
import sys
from pathlib import Path
import importlib.util
import numpy as np
import torch
import pytest
import tempfile

# Add paths
_src_dir = Path(__file__).parent
_repo_root = _src_dir.parent.parent.parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from train_mil import (
    SpotDataset,
    train_epoch,
    evaluate,
    train,
    collate_spots,
)


# Load ProportionGuidedMIL via importlib to avoid __init__.py issues
_spec = importlib.util.spec_from_file_location(
    'proportion_mil',
    _repo_root / 'CITEgeist/model/proportion_mil.py'
)
_mil_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mil_module)
ProportionGuidedMIL = _mil_module.ProportionGuidedMIL


class TestSpotDataset:
    """Tests for SpotDataset."""

    def test_spot_dataset_basic(self):
        """Test spot dataset loading."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create mock spot files
            for i in range(3):
                spot_dir = tmpdir / f"spot_{i}"
                spot_dir.mkdir()

                # Mock embeddings and proportions
                np.save(
                    spot_dir / "embeddings.npy",
                    np.random.randn(10, 768).astype(np.float32)
                )
                np.save(
                    spot_dir / "proportions.npy",
                    np.random.dirichlet(np.ones(5)).astype(np.float32)
                )

            dataset = SpotDataset(tmpdir, n_cell_types=5)

            assert len(dataset) == 3

            embeddings, proportions = dataset[0]
            assert embeddings.shape[1] == 768
            assert proportions.shape == (5,)
            assert isinstance(embeddings, torch.Tensor)
            assert isinstance(proportions, torch.Tensor)

    def test_spot_dataset_min_nuclei_filter(self):
        """Test that spots with too few nuclei are filtered."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create spot with few nuclei (should be filtered)
            spot_dir = tmpdir / "spot_0"
            spot_dir.mkdir()
            np.save(
                spot_dir / "embeddings.npy",
                np.random.randn(3, 768).astype(np.float32)  # Only 3 nuclei
            )
            np.save(
                spot_dir / "proportions.npy",
                np.random.dirichlet(np.ones(5)).astype(np.float32)
            )

            # Create spot with enough nuclei
            spot_dir = tmpdir / "spot_1"
            spot_dir.mkdir()
            np.save(
                spot_dir / "embeddings.npy",
                np.random.randn(10, 768).astype(np.float32)  # 10 nuclei
            )
            np.save(
                spot_dir / "proportions.npy",
                np.random.dirichlet(np.ones(5)).astype(np.float32)
            )

            dataset = SpotDataset(tmpdir, n_cell_types=5, min_nuclei=5)

            # Only spot_1 should be included
            assert len(dataset) == 1

    def test_spot_dataset_missing_files(self):
        """Test that spots with missing files are skipped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Spot with only embeddings
            spot_dir = tmpdir / "spot_0"
            spot_dir.mkdir()
            np.save(
                spot_dir / "embeddings.npy",
                np.random.randn(10, 768).astype(np.float32)
            )

            # Spot with both files
            spot_dir = tmpdir / "spot_1"
            spot_dir.mkdir()
            np.save(
                spot_dir / "embeddings.npy",
                np.random.randn(10, 768).astype(np.float32)
            )
            np.save(
                spot_dir / "proportions.npy",
                np.random.dirichlet(np.ones(5)).astype(np.float32)
            )

            dataset = SpotDataset(tmpdir, n_cell_types=5)

            # Only spot_1 should be included
            assert len(dataset) == 1


class TestCollateSpots:
    """Tests for collate function."""

    def test_collate_returns_list(self):
        """Test that collate returns batch as list."""
        batch = [
            (torch.randn(10, 768), torch.randn(5)),
            (torch.randn(15, 768), torch.randn(5)),
        ]

        result = collate_spots(batch)

        assert isinstance(result, list)
        assert len(result) == 2
        assert result[0][0].shape == (10, 768)
        assert result[1][0].shape == (15, 768)


class TestTrainEpoch:
    """Tests for training epoch."""

    def test_train_epoch_basic(self):
        """Test single training epoch."""
        model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

        # Mock dataloader
        mock_data = [
            (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
            (torch.randn(15, 768), torch.tensor([0.1, 0.4, 0.2, 0.2, 0.1])),
        ]

        loss = train_epoch(model, mock_data, optimizer)

        assert loss > 0
        assert isinstance(loss, float)

    def test_train_epoch_updates_weights(self):
        """Test that training actually updates model weights."""
        model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)

        # Get initial weights
        initial_weights = model.attention_V[0].weight.clone()

        mock_data = [
            (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
        ]

        train_epoch(model, mock_data, optimizer)

        # Check weights changed
        assert not torch.allclose(initial_weights, model.attention_V[0].weight)


class TestEvaluate:
    """Tests for evaluation."""

    def test_evaluate_basic(self):
        """Test evaluation."""
        model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)
        model.eval()

        mock_data = [
            (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
        ]

        metrics = evaluate(model, mock_data)

        assert 'loss' in metrics
        assert 'pearson_r' in metrics
        assert 'pearson_p' in metrics
        assert isinstance(metrics['loss'], float)
        assert -1 <= metrics['pearson_r'] <= 1

    def test_evaluate_multiple_spots(self):
        """Test evaluation with multiple spots."""
        model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)
        model.eval()

        mock_data = [
            (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
            (torch.randn(15, 768), torch.tensor([0.1, 0.4, 0.2, 0.2, 0.1])),
            (torch.randn(8, 768), torch.tensor([0.3, 0.2, 0.1, 0.3, 0.1])),
        ]

        metrics = evaluate(model, mock_data)

        assert metrics['loss'] > 0


class TestTrain:
    """Tests for full training loop."""

    def test_train_reduces_loss(self):
        """Test that training reduces loss over epochs."""
        torch.manual_seed(42)

        model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)

        # Create consistent training data
        train_data = [
            (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
            (torch.randn(15, 768), torch.tensor([0.1, 0.4, 0.2, 0.2, 0.1])),
        ]
        val_data = [
            (torch.randn(10, 768), torch.tensor([0.15, 0.35, 0.2, 0.2, 0.1])),
        ]

        history = train(
            model,
            train_data,
            val_data,
            n_epochs=20,
            lr=1e-3,
        )

        assert 'train_loss' in history
        assert 'val_loss' in history
        assert 'val_r' in history

        # Check loss generally decreases
        assert len(history['train_loss']) == 20
        # Early loss should be higher than late loss (on average)
        early_loss = np.mean(history['train_loss'][:5])
        late_loss = np.mean(history['train_loss'][-5:])
        assert late_loss < early_loss

    def test_train_saves_best_model(self):
        """Test that training saves best model checkpoint."""
        with tempfile.TemporaryDirectory() as tmpdir:
            save_path = Path(tmpdir) / "best_model.pt"

            model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)

            train_data = [
                (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
            ]
            val_data = [
                (torch.randn(10, 768), torch.tensor([0.15, 0.35, 0.2, 0.2, 0.1])),
            ]

            train(
                model,
                train_data,
                val_data,
                n_epochs=5,
                save_path=save_path,
            )

            assert save_path.exists()

            # Load and verify it's valid
            state_dict = torch.load(save_path, weights_only=True)
            assert 'attention_V.0.weight' in state_dict


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
