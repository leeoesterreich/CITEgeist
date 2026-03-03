# CITEgeist/tests/test_mae.py
"""Tests for Masked Autoencoder (MAE) module.

The MAE module implements self-supervised learning via masked image modeling
for 2-channel (DAPI + boundary) nucleus patches from spatial transcriptomics.

Architecture:
- Encoder: ViT-Small (384-dim, 12 layers, 6 heads)
- Decoder: Lightweight transformer (192-dim, 4 layers, 3 heads)
- Masking: 75% of patches masked during training
"""
import pytest
import torch
import numpy as np


class TestMAE:
    """Tests for the full Masked Autoencoder."""

    def test_mae_forward_returns_loss_and_reconstruction(self):
        """MAE forward pass should return loss, predictions, and mask.

        Forward pass:
        - Input: (B, 2, 96, 96) images
        - Returns: (loss, pred, mask)
          - loss: scalar MSE loss on masked patches
          - pred: (B, N, patch_pixels) where N=36, patch_pixels=16*16*2=512
          - mask: (B, N) binary mask (1=masked, 0=kept)
        """
        from CITEgeist.model.mae import MAE

        model = MAE(
            img_size=96,
            patch_size=16,
            in_chans=2,
            encoder_embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            decoder_embed_dim=192,
            decoder_depth=4,
            decoder_num_heads=3,
            mask_ratio=0.75,
        )
        model.train()

        # Batch of images
        x = torch.randn(4, 2, 96, 96)

        loss, pred, mask = model(x)

        # Loss should be a scalar
        assert loss.ndim == 0, f"Expected scalar loss, got shape {loss.shape}"
        assert loss.item() >= 0, "Loss should be non-negative"

        # Predictions should be (B, N, patch_pixels)
        num_patches = (96 // 16) ** 2  # 36
        patch_pixels = 16 * 16 * 2  # 512
        assert pred.shape == (4, num_patches, patch_pixels), \
            f"Expected pred shape (4, {num_patches}, {patch_pixels}), got {pred.shape}"

        # Mask should be (B, N)
        assert mask.shape == (4, num_patches), \
            f"Expected mask shape (4, {num_patches}), got {mask.shape}"

        # Mask should have ~75% ones (masked patches)
        mask_ratio_actual = mask.float().mean().item()
        assert 0.70 < mask_ratio_actual < 0.80, \
            f"Expected ~75% mask ratio, got {mask_ratio_actual:.2%}"

    def test_mae_encode_produces_embeddings(self):
        """MAE encode() should return CLS token embeddings for inference.

        encode() method:
        - Input: (B, 2, 96, 96) images
        - Output: (B, encoder_embed_dim) CLS token embeddings
        - No masking applied (full image encoded)
        """
        from CITEgeist.model.mae import MAE

        model = MAE(
            img_size=96,
            patch_size=16,
            in_chans=2,
            encoder_embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            decoder_embed_dim=192,
            decoder_depth=4,
            decoder_num_heads=3,
            mask_ratio=0.75,
        )
        model.eval()

        x = torch.randn(4, 2, 96, 96)

        with torch.no_grad():
            embeddings = model.encode(x)

        # Should return CLS token embeddings
        assert embeddings.shape == (4, 384), \
            f"Expected embeddings shape (4, 384), got {embeddings.shape}"

        # Embeddings should be deterministic in eval mode
        with torch.no_grad():
            embeddings2 = model.encode(x)

        assert torch.allclose(embeddings, embeddings2, atol=1e-6), \
            "encode() should be deterministic in eval mode"

    def test_mae_loss_decreases_with_training(self):
        """MAE loss should decrease over training steps.

        This validates that:
        1. Gradients flow correctly through encoder and decoder
        2. The model can learn to reconstruct masked patches
        3. The optimization objective is well-defined
        """
        from CITEgeist.model.mae import MAE

        model = MAE(
            img_size=96,
            patch_size=16,
            in_chans=2,
            encoder_embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            decoder_embed_dim=192,
            decoder_depth=4,
            decoder_num_heads=3,
            mask_ratio=0.75,
        )
        model.train()

        # Use fixed data for consistent training signal
        torch.manual_seed(42)
        x = torch.randn(8, 2, 96, 96)

        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)

        # Record losses over training steps
        losses = []
        for step in range(10):
            optimizer.zero_grad()
            loss, _, _ = model(x)
            loss.backward()
            optimizer.step()
            losses.append(loss.item())

        # Loss should decrease from start to end
        # Allow some noise but overall trend should be down
        initial_loss = np.mean(losses[:3])
        final_loss = np.mean(losses[-3:])

        assert final_loss < initial_loss, \
            f"Loss should decrease during training: initial={initial_loss:.4f}, final={final_loss:.4f}"


class TestMAEDecoder:
    """Tests for the MAE decoder architecture."""

    def test_decoder_output_shape(self):
        """Decoder should produce correct output shape.

        Input: (B, num_patches, decoder_embed_dim) + mask tokens
        Output: (B, num_patches, patch_pixels)
        """
        from CITEgeist.model.mae import MAEDecoder

        decoder = MAEDecoder(
            num_patches=36,
            patch_size=16,
            in_chans=2,
            embed_dim=192,
            depth=4,
            num_heads=3,
        )

        # Decoder input: encoded unmasked patches + mask tokens
        # For simplicity, test with full sequence
        x = torch.randn(4, 36, 192)
        ids_restore = torch.arange(36).unsqueeze(0).expand(4, -1)

        output = decoder(x, ids_restore)

        # Output should predict pixels for all patches
        patch_pixels = 16 * 16 * 2  # 512
        assert output.shape == (4, 36, patch_pixels), \
            f"Expected output shape (4, 36, {patch_pixels}), got {output.shape}"

    def test_decoder_handles_masked_input(self):
        """Decoder should handle partially masked input correctly.

        With 75% masking, only 25% of patches (9 out of 36) are visible.
        Decoder fills in mask tokens and predicts all patches.
        """
        from CITEgeist.model.mae import MAEDecoder

        decoder = MAEDecoder(
            num_patches=36,
            patch_size=16,
            in_chans=2,
            embed_dim=192,
            depth=4,
            num_heads=3,
        )

        # Only 9 visible patches (25% of 36)
        num_visible = 9
        x = torch.randn(4, num_visible, 192)

        # Random restore indices
        ids_restore = torch.stack([
            torch.randperm(36) for _ in range(4)
        ])

        output = decoder(x, ids_restore)

        # Should still output predictions for all 36 patches
        patch_pixels = 16 * 16 * 2
        assert output.shape == (4, 36, patch_pixels), \
            f"Expected output shape (4, 36, {patch_pixels}), got {output.shape}"


class TestMAEComponents:
    """Component tests for MAE internals."""

    def test_mask_token_is_learnable(self):
        """MAE decoder should have a learnable mask token."""
        from CITEgeist.model.mae import MAEDecoder

        decoder = MAEDecoder(
            num_patches=36,
            patch_size=16,
            in_chans=2,
            embed_dim=192,
            depth=4,
            num_heads=3,
        )

        # mask_token should be a learnable parameter
        assert hasattr(decoder, 'mask_token'), "Decoder should have mask_token"
        assert isinstance(decoder.mask_token, torch.nn.Parameter), \
            "mask_token should be a nn.Parameter"
        assert decoder.mask_token.shape == (1, 1, 192), \
            f"Expected mask_token shape (1, 1, 192), got {decoder.mask_token.shape}"

    def test_encoder_decoder_projection(self):
        """MAE should project encoder output to decoder dimension."""
        from CITEgeist.model.mae import MAE

        model = MAE(
            img_size=96,
            patch_size=16,
            in_chans=2,
            encoder_embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            decoder_embed_dim=192,
            decoder_depth=4,
            decoder_num_heads=3,
        )

        # Check encoder-to-decoder projection exists
        assert hasattr(model, 'enc_to_dec'), "MAE should have enc_to_dec projection"
        assert isinstance(model.enc_to_dec, torch.nn.Linear), \
            "enc_to_dec should be a Linear layer"
        assert model.enc_to_dec.in_features == 384, \
            f"Expected enc_to_dec input 384, got {model.enc_to_dec.in_features}"
        assert model.enc_to_dec.out_features == 192, \
            f"Expected enc_to_dec output 192, got {model.enc_to_dec.out_features}"

    def test_mae_gradient_flow(self):
        """Gradients should flow through both encoder and decoder."""
        from CITEgeist.model.mae import MAE

        model = MAE(
            img_size=96,
            patch_size=16,
            in_chans=2,
            encoder_embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            decoder_embed_dim=192,
            decoder_depth=4,
            decoder_num_heads=3,
        )
        model.train()

        x = torch.randn(2, 2, 96, 96)
        loss, _, _ = model(x)

        # Backward pass
        loss.backward()

        # Check encoder has gradients
        encoder_params_with_grad = sum(
            1 for p in model.encoder.parameters() if p.grad is not None
        )
        total_encoder_params = sum(1 for _ in model.encoder.parameters())
        assert encoder_params_with_grad == total_encoder_params, \
            f"Not all encoder params have gradients: {encoder_params_with_grad}/{total_encoder_params}"

        # Check decoder has gradients
        decoder_params_with_grad = sum(
            1 for p in model.decoder.parameters() if p.grad is not None
        )
        total_decoder_params = sum(1 for _ in model.decoder.parameters())
        assert decoder_params_with_grad == total_decoder_params, \
            f"Not all decoder params have gradients: {decoder_params_with_grad}/{total_decoder_params}"


class TestMAEEdgeCases:
    """Edge case tests for MAE."""

    def test_different_mask_ratios(self):
        """MAE should work with different mask ratios."""
        from CITEgeist.model.mae import MAE

        for mask_ratio in [0.5, 0.75, 0.9]:
            model = MAE(
                img_size=96,
                patch_size=16,
                in_chans=2,
                encoder_embed_dim=384,
                encoder_depth=12,
                encoder_num_heads=6,
                decoder_embed_dim=192,
                decoder_depth=4,
                decoder_num_heads=3,
                mask_ratio=mask_ratio,
            )

            x = torch.randn(2, 2, 96, 96)
            loss, pred, mask = model(x)

            # Check mask ratio is approximately correct
            actual_ratio = mask.float().mean().item()
            assert abs(actual_ratio - mask_ratio) < 0.1, \
                f"Expected mask ratio ~{mask_ratio}, got {actual_ratio}"

    def test_override_mask_ratio_in_forward(self):
        """Forward should allow overriding mask_ratio."""
        from CITEgeist.model.mae import MAE

        model = MAE(
            img_size=96,
            patch_size=16,
            in_chans=2,
            encoder_embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            decoder_embed_dim=192,
            decoder_depth=4,
            decoder_num_heads=3,
            mask_ratio=0.75,  # Default
        )

        x = torch.randn(2, 2, 96, 96)

        # Override to 50% masking
        loss, pred, mask = model(x, mask_ratio=0.5)

        actual_ratio = mask.float().mean().item()
        assert abs(actual_ratio - 0.5) < 0.1, \
            f"Expected mask ratio ~0.5 when overridden, got {actual_ratio}"

    def test_single_sample_batch(self):
        """MAE should work with batch size 1."""
        from CITEgeist.model.mae import MAE

        model = MAE(
            img_size=96,
            patch_size=16,
            in_chans=2,
            encoder_embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            decoder_embed_dim=192,
            decoder_depth=4,
            decoder_num_heads=3,
        )

        x = torch.randn(1, 2, 96, 96)
        loss, pred, mask = model(x)

        assert pred.shape == (1, 36, 512)
        assert mask.shape == (1, 36)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
