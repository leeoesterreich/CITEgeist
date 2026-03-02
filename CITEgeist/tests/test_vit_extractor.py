"""Tests for ViT feature extraction.

This module tests the ViTFeatureExtractor class which wraps timm ViT models
for extracting features from 224x224 H&E patches. It supports multiple model
variants including UNI foundation model weights.

Unlike vit_encoder.py (custom 2-channel ViT for DAPI patches), this module
uses pretrained models from timm for 3-channel RGB histopathology images.

Note on running tests:
    On some HPC clusters, the CITEgeist.model package's __init__.py imports
    libraries (pyarrow, pandas) that require newer libstdc++ versions than
    available on login nodes. If tests fail with GLIBCXX errors:

    1. Run tests via SLURM compute nodes (recommended)
    2. Or clear LD_LIBRARY_PATH before running:
       $ unset LD_LIBRARY_PATH && pytest CITEgeist/tests/test_vit_extractor.py
    3. Or test via direct import (see standalone test at end of file)
"""
import numpy as np
import torch
import pytest


class TestViTFeatureExtractor:
    """Tests for ViT feature extraction from H&E patches."""

    def test_vit_extractor_shape(self):
        """Test ViT output shape.

        vit_small_patch16_224 produces 384-dimensional embeddings.
        """
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

        # Mock batch of 4 RGB images at 224x224
        x = torch.randn(4, 3, 224, 224)

        with torch.no_grad():
            features = extractor(x)

        assert features.shape == (4, 384), f"Expected (4, 384), got {features.shape}"

    def test_vit_extractor_numpy(self):
        """Test numpy input handling.

        extract_numpy() should accept numpy arrays and return numpy arrays.
        """
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

        # Numpy input in [0, 1] range
        x = np.random.rand(2, 3, 224, 224).astype(np.float32)

        features = extractor.extract_numpy(x)

        assert isinstance(features, np.ndarray), f"Expected np.ndarray, got {type(features)}"
        assert features.shape == (2, 384), f"Expected (2, 384), got {features.shape}"

    def test_vit_extractor_imagenet_normalization(self):
        """Test ImageNet normalization is applied.

        Input values in [0, 1] should be normalized using ImageNet statistics.
        """
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

        # Input in [0, 1]
        x = torch.ones(1, 3, 224, 224) * 0.5

        # After normalization, values should be different
        normalized = extractor.normalize(x)

        # ImageNet mean ~ 0.485, so (0.5 - 0.485) / 0.229 ~ 0.065 for R channel
        assert not torch.allclose(normalized, x), "Normalized output should differ from input"

        # Check that normalization is applied correctly for red channel
        # (0.5 - 0.485) / 0.229 = 0.0655...
        expected_red = (0.5 - 0.485) / 0.229
        assert torch.allclose(normalized[0, 0, 0, 0], torch.tensor(expected_red), atol=1e-4), \
            f"Expected red channel normalized value ~{expected_red}, got {normalized[0, 0, 0, 0]}"

    def test_vit_extractor_embed_dim_attribute(self):
        """Test that embed_dim attribute is correctly set."""
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        # vit_small has 384-dim output
        extractor_small = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)
        assert extractor_small.embed_dim == 384, f"Expected 384, got {extractor_small.embed_dim}"

        # vit_base has 768-dim output
        extractor_base = ViTFeatureExtractor(model_name='vit_base_patch16_224', pretrained=False)
        assert extractor_base.embed_dim == 768, f"Expected 768, got {extractor_base.embed_dim}"

    def test_vit_extractor_deterministic(self):
        """Test that extractor is deterministic in eval mode."""
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

        torch.manual_seed(42)
        x = torch.randn(2, 3, 224, 224)

        with torch.no_grad():
            output1 = extractor(x)
            output2 = extractor(x)

        assert torch.allclose(output1, output2, atol=1e-6), \
            "Extractor should be deterministic in eval mode"

    def test_vit_extractor_batched_numpy(self):
        """Test extract_numpy with batch processing."""
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

        # 10 patches, should be processed in batches
        x = np.random.rand(10, 3, 224, 224).astype(np.float32)

        features = extractor.extract_numpy(x, batch_size=4)

        assert features.shape == (10, 384), f"Expected (10, 384), got {features.shape}"


class TestViTExtractorModels:
    """Tests for different ViT model variants."""

    def test_vit_base_shape(self):
        """Test vit_base_patch16_224 produces 768-dim embeddings."""
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_base_patch16_224', pretrained=False)

        x = torch.randn(2, 3, 224, 224)

        with torch.no_grad():
            features = extractor(x)

        assert features.shape == (2, 768), f"Expected (2, 768), got {features.shape}"

    def test_vit_large_shape(self):
        """Test vit_large_patch16_224 produces 1024-dim embeddings."""
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_large_patch16_224', pretrained=False)

        x = torch.randn(2, 3, 224, 224)

        with torch.no_grad():
            features = extractor(x)

        assert features.shape == (2, 1024), f"Expected (2, 1024), got {features.shape}"


class TestLoadUNI:
    """Tests for UNI foundation model loading."""

    def test_load_uni_extractor_function_exists(self):
        """Test that load_uni_extractor function is importable."""
        from CITEgeist.model.vit_extractor import load_uni_extractor

        assert callable(load_uni_extractor), "load_uni_extractor should be callable"

    def test_load_uni_without_weights(self):
        """Test UNI model structure without actual weights.

        UNI is a ViT-L/16 so should have 1024-dim output.
        """
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        # Create ViT-L without pretrained weights (mimics UNI architecture)
        extractor = ViTFeatureExtractor(
            model_name='vit_large_patch16_224',
            pretrained=False,
            weights_path=None,
        )

        x = torch.randn(2, 3, 224, 224)

        with torch.no_grad():
            features = extractor(x)

        assert features.shape == (2, 1024), f"Expected (2, 1024), got {features.shape}"
        assert extractor.embed_dim == 1024, f"Expected embed_dim 1024, got {extractor.embed_dim}"


class TestViTExtractorEdgeCases:
    """Edge case tests for ViT extractor."""

    def test_single_image(self):
        """Test with batch size 1."""
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

        x = torch.randn(1, 3, 224, 224)

        with torch.no_grad():
            features = extractor(x)

        assert features.shape == (1, 384)

    def test_empty_batch_numpy(self):
        """Test extract_numpy handles empty input gracefully."""
        from CITEgeist.model.vit_extractor import ViTFeatureExtractor

        extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

        # Empty input
        x = np.empty((0, 3, 224, 224), dtype=np.float32)

        features = extractor.extract_numpy(x)

        assert features.shape == (0, 384), f"Expected (0, 384), got {features.shape}"


def run_standalone_tests():
    """Run tests with direct imports (bypasses CITEgeist.model.__init__.py).

    Use this when pytest fails due to GLIBCXX library version issues on HPC clusters.
    """
    import importlib.util
    import os

    # Load module directly
    module_path = os.path.join(
        os.path.dirname(__file__),
        "..", "model", "vit_extractor.py"
    )
    spec = importlib.util.spec_from_file_location("vit_extractor", module_path)
    vit_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vit_module)

    ViTFeatureExtractor = vit_module.ViTFeatureExtractor
    load_uni_extractor = vit_module.load_uni_extractor

    print("=" * 60)
    print("Running standalone tests (bypassing package __init__.py)")
    print("=" * 60)

    # Test 1: Shape
    print("\nTest 1: vit_small output shape...")
    extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)
    x = torch.randn(4, 3, 224, 224)
    with torch.no_grad():
        features = extractor(x)
    assert features.shape == (4, 384)
    print("  PASSED")

    # Test 2: NumPy
    print("\nTest 2: NumPy input...")
    x_np = np.random.rand(2, 3, 224, 224).astype(np.float32)
    features_np = extractor.extract_numpy(x_np)
    assert isinstance(features_np, np.ndarray)
    assert features_np.shape == (2, 384)
    print("  PASSED")

    # Test 3: Normalization
    print("\nTest 3: ImageNet normalization...")
    x_norm = torch.ones(1, 3, 224, 224) * 0.5
    normalized = extractor.normalize(x_norm)
    expected_red = (0.5 - 0.485) / 0.229
    assert torch.allclose(normalized[0, 0, 0, 0], torch.tensor(expected_red), atol=1e-4)
    print("  PASSED")

    # Test 4: embed_dim
    print("\nTest 4: embed_dim attribute...")
    assert extractor.embed_dim == 384
    extractor_base = ViTFeatureExtractor(model_name='vit_base_patch16_224', pretrained=False)
    assert extractor_base.embed_dim == 768
    print("  PASSED")

    # Test 5: vit_large
    print("\nTest 5: vit_large output shape...")
    extractor_large = ViTFeatureExtractor(model_name='vit_large_patch16_224', pretrained=False)
    x_large = torch.randn(2, 3, 224, 224)
    with torch.no_grad():
        features_large = extractor_large(x_large)
    assert features_large.shape == (2, 1024)
    assert extractor_large.embed_dim == 1024
    print("  PASSED")

    # Test 6: Empty batch
    print("\nTest 6: Empty batch...")
    x_empty = np.empty((0, 3, 224, 224), dtype=np.float32)
    features_empty = extractor.extract_numpy(x_empty)
    assert features_empty.shape == (0, 384)
    print("  PASSED")

    # Test 7: load_uni_extractor
    print("\nTest 7: load_uni_extractor callable...")
    assert callable(load_uni_extractor)
    print("  PASSED")

    print("\n" + "=" * 60)
    print("ALL STANDALONE TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    import sys

    if "--standalone" in sys.argv:
        run_standalone_tests()
    else:
        pytest.main([__file__, "-v"])
