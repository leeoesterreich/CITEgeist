"""Tests for morphology backbone abstraction."""
import numpy as np
import pytest
import torch

from CITEgeist.model.morphology.morphology_backbone import DAPIBackbone, HEBackbone


def test_dapi_backbone_output_shape():
    """DAPI backbone should output (N, 384) from (N, 2, 96, 96) input."""
    backbone = DAPIBackbone()
    patches = torch.randn(4, 2, 96, 96)
    with torch.no_grad():
        embeddings = backbone.extract(patches)
    assert embeddings.shape == (4, 384)


def test_he_backbone_output_shape():
    """H&E backbone should output (N, 384) from (N, 3, 224, 224) input."""
    backbone = HEBackbone(model_name='vit_small_patch16_224')
    patches = torch.randn(4, 3, 224, 224)
    with torch.no_grad():
        embeddings = backbone.extract(patches)
    assert embeddings.shape == (4, 384)


def test_dapi_backbone_extract_numpy():
    """extract_numpy should handle numpy input and batching."""
    backbone = DAPIBackbone()
    patches_np = np.random.randn(10, 2, 96, 96).astype(np.float32)
    embeddings = backbone.extract_numpy(patches_np, batch_size=4)
    assert embeddings.shape == (10, 384)
    assert embeddings.dtype == np.float32
