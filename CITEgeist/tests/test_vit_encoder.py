"""
Tests for ViT-Small encoder for nucleus morphology patches.

The ViT encoder processes 2-channel (DAPI + boundary) nucleus patches
and produces embeddings for downstream self-supervised learning (MAE, DINO).
"""

import pytest
import torch
import numpy as np


class TestViTEncoder:
    """Tests for the ViT-Small encoder architecture."""

    def test_vit_encoder_output_shape(self):
        """Test that encoder produces correct output shape.

        Input: (B, 2, 96, 96) - batch of 2-channel 96x96 patches
        Output: (B, 384) - [CLS] token embeddings
        """
        from CITEgeist.model.vit_encoder import ViTEncoder

        model = ViTEncoder(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
            depth=12,
            num_heads=6,
            mlp_ratio=4.0,
        )
        model.eval()

        # Batch of 4 images, 2 channels (DAPI + boundary), 96x96
        x = torch.randn(4, 2, 96, 96)

        with torch.no_grad():
            output = model(x)

        # Should return [CLS] token embedding of dimension 384
        assert output.shape == (4, 384), f"Expected (4, 384), got {output.shape}"

    def test_vit_encoder_patch_tokens(self):
        """Test that encoder can return patch tokens for MAE decoder.

        When return_patch_tokens=True:
        - Returns tuple (cls_token, patch_tokens)
        - cls_token: (B, 384)
        - patch_tokens: (B, num_patches, 384) where num_patches = (96/16)^2 = 36
        """
        from CITEgeist.model.vit_encoder import ViTEncoder

        model = ViTEncoder(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
            depth=12,
            num_heads=6,
            mlp_ratio=4.0,
        )
        model.eval()

        x = torch.randn(4, 2, 96, 96)

        with torch.no_grad():
            cls_token, patch_tokens = model(x, return_patch_tokens=True)

        # CLS token shape
        assert cls_token.shape == (4, 384), f"Expected cls_token (4, 384), got {cls_token.shape}"

        # Patch tokens: 96/16 = 6 patches per side, 6*6 = 36 patches total
        num_patches = (96 // 16) ** 2
        assert num_patches == 36, f"Expected 36 patches, calculated {num_patches}"
        assert patch_tokens.shape == (4, 36, 384), f"Expected patch_tokens (4, 36, 384), got {patch_tokens.shape}"

    def test_vit_encoder_deterministic(self):
        """Test that encoder is deterministic in eval mode.

        Same input should produce same output when model is in eval mode.
        """
        from CITEgeist.model.vit_encoder import ViTEncoder

        model = ViTEncoder(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
            depth=12,
            num_heads=6,
            mlp_ratio=4.0,
        )
        model.eval()

        # Fix seed for reproducible input
        torch.manual_seed(42)
        x = torch.randn(4, 2, 96, 96)

        with torch.no_grad():
            output1 = model(x)
            output2 = model(x)

        # Outputs should be identical
        assert torch.allclose(output1, output2, atol=1e-6), "Encoder should be deterministic in eval mode"

    def test_patch_embed_output_shape(self):
        """Test PatchEmbed produces correct output shape."""
        from CITEgeist.model.vit_encoder import PatchEmbed

        patch_embed = PatchEmbed(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
        )

        x = torch.randn(4, 2, 96, 96)
        patches = patch_embed(x)

        # Should produce (B, num_patches, embed_dim)
        # num_patches = (96/16)^2 = 36
        assert patches.shape == (4, 36, 384), f"Expected (4, 36, 384), got {patches.shape}"

    def test_attention_output_shape(self):
        """Test multi-head attention output shape."""
        from CITEgeist.model.vit_encoder import Attention

        attn = Attention(
            dim=384,
            num_heads=6,
        )

        # Sequence of 37 tokens (1 CLS + 36 patches)
        x = torch.randn(4, 37, 384)
        output = attn(x)

        assert output.shape == x.shape, f"Expected {x.shape}, got {output.shape}"

    def test_mlp_output_shape(self):
        """Test MLP block output shape."""
        from CITEgeist.model.vit_encoder import MLP

        mlp = MLP(
            in_features=384,
            hidden_features=384 * 4,  # mlp_ratio=4
        )

        x = torch.randn(4, 37, 384)
        output = mlp(x)

        assert output.shape == x.shape, f"Expected {x.shape}, got {output.shape}"

    def test_transformer_block_output_shape(self):
        """Test TransformerBlock output shape (pre-norm architecture)."""
        from CITEgeist.model.vit_encoder import TransformerBlock

        block = TransformerBlock(
            dim=384,
            num_heads=6,
            mlp_ratio=4.0,
        )

        x = torch.randn(4, 37, 384)
        output = block(x)

        assert output.shape == x.shape, f"Expected {x.shape}, got {output.shape}"

    def test_vit_encoder_gradient_flow(self):
        """Test that gradients flow through the encoder.

        We check that:
        1. All parameters receive gradient tensors (grad is not None)
        2. The output requires gradients (connected to computational graph)
        3. Backward pass completes without error
        """
        from CITEgeist.model.vit_encoder import ViTEncoder

        model = ViTEncoder(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
            depth=12,
            num_heads=6,
            mlp_ratio=4.0,
        )
        model.train()

        x = torch.randn(2, 2, 96, 96)
        output = model(x)

        # Check output requires grad (connected to graph)
        assert output.requires_grad, "Output should require gradients in train mode"

        # Compute a simple loss and backpropagate
        loss = output.sum()
        loss.backward()

        # Check that all model parameters have gradient tensors (not None)
        # Note: The actual gradient values may be very small or occasionally zero
        # due to initialization, but the gradient tensor should exist
        params_without_grad = []
        for name, param in model.named_parameters():
            if param.grad is None:
                params_without_grad.append(name)

        assert len(params_without_grad) == 0, f"Parameters without gradients: {params_without_grad}"

    def test_vit_encoder_parameter_count(self):
        """Test that model has expected number of parameters (ViT-Small scale)."""
        from CITEgeist.model.vit_encoder import ViTEncoder

        model = ViTEncoder(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
            depth=12,
            num_heads=6,
            mlp_ratio=4.0,
        )

        # Count parameters
        num_params = sum(p.numel() for p in model.parameters())

        # ViT-Small should have ~21-22M parameters
        # With 2-channel input instead of 3, slightly fewer in patch embed
        # Expected: ~21M params (give or take based on exact config)
        assert 20_000_000 < num_params < 25_000_000, \
            f"Expected ~21M parameters for ViT-Small, got {num_params:,}"


class TestViTEncoderEdgeCases:
    """Edge case tests for ViT encoder."""

    def test_single_batch(self):
        """Test with batch size 1."""
        from CITEgeist.model.vit_encoder import ViTEncoder

        model = ViTEncoder(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
            depth=12,
            num_heads=6,
            mlp_ratio=4.0,
        )
        model.eval()

        x = torch.randn(1, 2, 96, 96)

        with torch.no_grad():
            output = model(x)

        assert output.shape == (1, 384)

    def test_large_batch(self):
        """Test with larger batch size."""
        from CITEgeist.model.vit_encoder import ViTEncoder

        model = ViTEncoder(
            img_size=96,
            patch_size=16,
            in_chans=2,
            embed_dim=384,
            depth=12,
            num_heads=6,
            mlp_ratio=4.0,
        )
        model.eval()

        x = torch.randn(32, 2, 96, 96)

        with torch.no_grad():
            output = model(x)

        assert output.shape == (32, 384)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
