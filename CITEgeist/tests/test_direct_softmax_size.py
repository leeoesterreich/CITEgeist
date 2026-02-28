"""Tests for DirectSoftmaxModel with size features."""
import torch
import pytest

from CITEgeist.model.vae import VAEEncoder
from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel


class TestDirectSoftmaxWithSize:
    """Test size feature integration."""

    @pytest.fixture
    def encoder(self):
        return VAEEncoder(in_channels=2, latent_dim=128)

    def test_model_accepts_size_features(self, encoder):
        """Test that model can be created with size features enabled."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=True,
        )
        assert model.use_size_features is True
        # Projection input should be 128 + 3 = 131
        assert model.projection[0].in_features == 131

    def test_forward_with_size_features(self, encoder):
        """Test forward pass with size features."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=True,
        )

        patches = torch.randn(10, 2, 96, 96)
        proportions = torch.ones(7) / 7
        size_features = torch.randn(10, 3)

        loss, assignments = model(patches, proportions, size_features=size_features)

        assert loss.shape == ()
        assert assignments.shape == (10, 7)

    def test_forward_without_size_features_raises(self, encoder):
        """Test that missing size features raises error when enabled."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=True,
        )

        patches = torch.randn(10, 2, 96, 96)
        proportions = torch.ones(7) / 7

        with pytest.raises(ValueError, match="size_features required"):
            model(patches, proportions, size_features=None)

    def test_model_without_size_features(self, encoder):
        """Test that model works without size features."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=False,
        )

        patches = torch.randn(10, 2, 96, 96)
        proportions = torch.ones(7) / 7

        loss, assignments = model(patches, proportions)

        assert loss.shape == ()
        assert assignments.shape == (10, 7)

    def test_assign_with_size_features(self, encoder):
        """Test assign method with size features."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=True,
        )

        patches = torch.randn(10, 2, 96, 96)
        proportions = torch.ones(7) / 7
        size_features = torch.randn(10, 3)

        assignments, confidence = model.assign(
            patches, proportions=proportions, size_features=size_features
        )

        assert assignments.shape == (10,)
        assert confidence.shape == (10,)
