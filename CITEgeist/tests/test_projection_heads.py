"""Tests for projection heads and prototypes.

This module tests the projection heads and prototypes used for cell type
assignment in the VAE + Sinkhorn system. The projection heads learn which
features in the VAE latent space are relevant for each cell type, and the
prototypes represent ideal representations of each cell type.
"""
import sys
import os
import importlib.util

import torch
import pytest


def _import_projection_heads():
    """Import projection_heads module directly without triggering model __init__.py."""
    module_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "model",
        "projection_heads.py"
    )
    spec = importlib.util.spec_from_file_location("projection_heads", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestProjectionHeads:
    """Tests for ProjectionHeads class."""

    def test_projection_heads_output_shape(self):
        """Each head should project to correct dimension."""
        proj_module = _import_projection_heads()
        ProjectionHeads = proj_module.ProjectionHeads

        heads = ProjectionHeads(
            input_dim=128,
            projection_dim=32,
            n_types=7
        )

        z = torch.randn(10, 128)  # 10 nuclei

        projected = heads(z)  # List of K tensors

        assert len(projected) == 7
        for p in projected:
            assert p.shape == (10, 32)

    def test_projection_heads_gradients(self):
        """Projection heads should have learnable parameters."""
        proj_module = _import_projection_heads()
        ProjectionHeads = proj_module.ProjectionHeads

        heads = ProjectionHeads(input_dim=128, projection_dim=32, n_types=5)

        z = torch.randn(10, 128, requires_grad=True)
        projected = heads(z)

        # Sum all outputs and backprop
        loss = sum(p.sum() for p in projected)
        loss.backward()

        # Check that gradients flow to input
        assert z.grad is not None
        assert z.grad.shape == (10, 128)

    def test_projection_heads_different_configs(self):
        """Should work with different input/output dimensions."""
        proj_module = _import_projection_heads()
        ProjectionHeads = proj_module.ProjectionHeads

        # Larger dimensions
        heads = ProjectionHeads(input_dim=256, projection_dim=64, n_types=10)
        z = torch.randn(5, 256)
        projected = heads(z)

        assert len(projected) == 10
        for p in projected:
            assert p.shape == (5, 64)


class TestPrototypes:
    """Tests for Prototypes class."""

    def test_prototypes_shape(self):
        """Prototypes should have correct shape."""
        proj_module = _import_projection_heads()
        Prototypes = proj_module.Prototypes

        prototypes = Prototypes(projection_dim=32, n_types=7)
        proto = prototypes()

        assert proto.shape == (7, 32)
        assert proto.requires_grad  # Should be learnable

    def test_prototypes_different_config(self):
        """Should work with different dimensions."""
        proj_module = _import_projection_heads()
        Prototypes = proj_module.Prototypes

        prototypes = Prototypes(projection_dim=64, n_types=10)
        proto = prototypes()

        assert proto.shape == (10, 64)

    def test_prototypes_init_from_kmeans(self):
        """K-means initialization should update prototype values."""
        # Skip if sklearn has library issues (e.g., GLIBCXX mismatch on HPC)
        try:
            from sklearn.cluster import KMeans
        except ImportError:
            pytest.skip("sklearn not available or has library issues")

        proj_module = _import_projection_heads()
        Prototypes = proj_module.Prototypes

        prototypes = Prototypes(projection_dim=32, n_types=3)

        # Create synthetic projected latents with clear clusters
        projected = [
            torch.randn(100, 32) + torch.tensor([10.0] + [0.0] * 31),  # cluster at x=10
            torch.randn(100, 32) + torch.tensor([0.0, 10.0] + [0.0] * 30),  # cluster at y=10
            torch.randn(100, 32) + torch.tensor([0.0, 0.0, 10.0] + [0.0] * 29),  # cluster at z=10
        ]

        proto_before = prototypes().clone().detach()
        prototypes.init_from_kmeans(projected)
        proto_after = prototypes()

        # Prototypes should have changed
        assert not torch.allclose(proto_before, proto_after)

        # Centroids should be near the expected values
        assert proto_after[0, 0].item() > 5  # First prototype near x=10
        assert proto_after[1, 1].item() > 5  # Second prototype near y=10
        assert proto_after[2, 2].item() > 5  # Third prototype near z=10


class TestComputeDistances:
    """Tests for compute_distances function."""

    def test_compute_distances_shape(self):
        """Distance computation should have correct shape."""
        proj_module = _import_projection_heads()
        ProjectionHeads = proj_module.ProjectionHeads
        Prototypes = proj_module.Prototypes
        compute_distances = proj_module.compute_distances

        heads = ProjectionHeads(input_dim=128, projection_dim=32, n_types=7)
        prototypes = Prototypes(projection_dim=32, n_types=7)

        z = torch.randn(10, 128)
        projected = heads(z)
        proto = prototypes()

        distances = compute_distances(projected, proto)

        assert distances.shape == (10, 7)
        assert (distances >= 0).all()

    def test_compute_distances_zero_for_identical(self):
        """Distance should be zero when projected equals prototype."""
        proj_module = _import_projection_heads()
        compute_distances = proj_module.compute_distances

        # Create projected latents identical to prototypes
        K, D = 3, 16
        proto = torch.randn(K, D)
        projected = [proto[k:k+1, :].repeat(5, 1) for k in range(K)]

        distances = compute_distances(projected, proto)

        # Diagonal should be zero (each projected[k] matches proto[k])
        for k in range(K):
            assert torch.allclose(distances[:, k], torch.zeros(5), atol=1e-6)

    def test_compute_distances_non_zero_for_different(self):
        """Distance should be non-zero for different vectors."""
        proj_module = _import_projection_heads()
        compute_distances = proj_module.compute_distances

        K, D, N = 3, 16, 10
        proto = torch.zeros(K, D)
        projected = [torch.ones(N, D) for _ in range(K)]  # All ones, proto all zeros

        distances = compute_distances(projected, proto)

        # Distance from ones to zeros should be sqrt(D)
        expected_dist = torch.sqrt(torch.tensor(D, dtype=torch.float32))
        assert torch.allclose(distances, expected_dist.expand(N, K), atol=1e-5)


class TestIntegration:
    """Integration tests for projection heads and prototypes."""

    def test_full_pipeline(self):
        """Test full projection -> distance -> softmax pipeline."""
        proj_module = _import_projection_heads()
        ProjectionHeads = proj_module.ProjectionHeads
        Prototypes = proj_module.Prototypes
        compute_distances = proj_module.compute_distances

        N, K = 100, 5
        input_dim, proj_dim = 128, 32

        heads = ProjectionHeads(input_dim=input_dim, projection_dim=proj_dim, n_types=K)
        prototypes = Prototypes(projection_dim=proj_dim, n_types=K)

        z = torch.randn(N, input_dim)
        projected = heads(z)
        proto = prototypes()
        distances = compute_distances(projected, proto)

        # Convert distances to soft assignments via softmax
        temperature = 0.1
        logits = -distances / temperature
        soft_assignments = torch.softmax(logits, dim=1)

        assert soft_assignments.shape == (N, K)
        assert torch.allclose(soft_assignments.sum(dim=1), torch.ones(N), atol=1e-5)
        assert (soft_assignments >= 0).all()
        assert (soft_assignments <= 1).all()

    def test_gradient_flow(self):
        """Gradients should flow through the entire pipeline."""
        proj_module = _import_projection_heads()
        ProjectionHeads = proj_module.ProjectionHeads
        Prototypes = proj_module.Prototypes
        compute_distances = proj_module.compute_distances

        N, K = 10, 3
        input_dim, proj_dim = 64, 16

        heads = ProjectionHeads(input_dim=input_dim, projection_dim=proj_dim, n_types=K)
        prototypes = Prototypes(projection_dim=proj_dim, n_types=K)

        z = torch.randn(N, input_dim, requires_grad=True)
        projected = heads(z)
        proto = prototypes()
        distances = compute_distances(projected, proto)

        # Simulate a loss
        loss = distances.mean()
        loss.backward()

        # Check gradient flow
        assert z.grad is not None

        # Check that prototype gradients exist
        assert prototypes.prototypes.grad is not None

        # Check that head parameters have gradients
        for head in heads.heads:
            for param in head.parameters():
                assert param.grad is not None
