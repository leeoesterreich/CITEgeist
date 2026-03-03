"""
Tests for DINO self-supervised learning module.

DINO uses teacher-student self-distillation with:
- Multi-crop training (global + local crops)
- EMA teacher update
- Centering to prevent collapse
"""

import pytest
import torch
import numpy as np


class TestDINOHead:
    """Tests for the DINO projection head."""

    def test_dino_head_output_shape(self):
        """Test that DINOHead produces correct output shape.

        Input: (B, 384) - encoder embeddings
        Output: (B, 4096) - projection head outputs
        """
        from CITEgeist.model.dino import DINOHead

        head = DINOHead(
            in_dim=384,
            hidden_dim=2048,
            bottleneck_dim=256,
            out_dim=4096,
        )
        head.eval()

        # Batch of 8 embeddings from encoder
        x = torch.randn(8, 384)

        with torch.no_grad():
            output = head(x)

        assert output.shape == (8, 4096), f"Expected (8, 4096), got {output.shape}"

    def test_dino_head_custom_dims(self):
        """Test DINOHead with custom dimensions."""
        from CITEgeist.model.dino import DINOHead

        head = DINOHead(
            in_dim=256,
            hidden_dim=1024,
            bottleneck_dim=128,
            out_dim=2048,
        )

        x = torch.randn(4, 256)
        output = head(x)

        assert output.shape == (4, 2048), f"Expected (4, 2048), got {output.shape}"

    def test_dino_head_gradient_flow(self):
        """Test that gradients flow through DINOHead."""
        from CITEgeist.model.dino import DINOHead

        head = DINOHead(in_dim=384, hidden_dim=2048, bottleneck_dim=256, out_dim=4096)
        head.train()

        x = torch.randn(4, 384, requires_grad=True)
        output = head(x)

        assert output.requires_grad, "Output should require gradients"

        loss = output.sum()
        loss.backward()

        # Check all parameters have gradients
        for name, param in head.named_parameters():
            assert param.grad is not None, f"Parameter {name} has no gradient"


class TestDINO:
    """Tests for the full DINO model."""

    def test_dino_forward_returns_loss(self):
        """Test that DINO forward returns a scalar loss.

        DINO takes global crops (processed by both teacher and student)
        and optional local crops (student only), returning cross-entropy loss.
        """
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            out_dim=4096,
            hidden_dim=2048,
            bottleneck_dim=256,
        )
        model.train()

        # Create 2 global crops (required)
        batch_size = 4
        global_crops = [
            torch.randn(batch_size, 2, 96, 96),
            torch.randn(batch_size, 2, 96, 96),
        ]

        # Optional: 2 local crops (smaller, but we use same size for simplicity)
        local_crops = [
            torch.randn(batch_size, 2, 96, 96),
            torch.randn(batch_size, 2, 96, 96),
        ]

        loss = model(global_crops, local_crops)

        # Loss should be a scalar tensor
        assert loss.dim() == 0, f"Expected scalar loss, got shape {loss.shape}"
        assert loss.requires_grad, "Loss should require gradients for backprop"
        assert not torch.isnan(loss), "Loss should not be NaN"
        assert not torch.isinf(loss), "Loss should not be infinite"

    def test_dino_forward_global_only(self):
        """Test DINO forward with global crops only (no local crops)."""
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
        )
        model.train()

        batch_size = 4
        global_crops = [
            torch.randn(batch_size, 2, 96, 96),
            torch.randn(batch_size, 2, 96, 96),
        ]

        # No local crops
        loss = model(global_crops, local_crops=None)

        assert loss.dim() == 0, "Loss should be scalar"
        assert not torch.isnan(loss), "Loss should not be NaN"

    def test_dino_encode_produces_embeddings(self):
        """Test that encode() returns embeddings from student encoder.

        The encode method is used for inference to get [CLS] token embeddings.
        """
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
        )

        batch_size = 8
        imgs = torch.randn(batch_size, 2, 96, 96)

        embeddings = model.encode(imgs)

        # Should return (B, embed_dim) embeddings
        assert embeddings.shape == (batch_size, 384), \
            f"Expected ({batch_size}, 384), got {embeddings.shape}"
        assert not embeddings.requires_grad, "Encode should return detached embeddings"

    def test_dino_teacher_momentum_update(self):
        """Test that update_teacher() applies EMA correctly.

        After update: teacher = momentum * teacher + (1 - momentum) * student
        """
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
        )

        # Get initial teacher parameter (deep copy)
        teacher_param = next(model.teacher_encoder.parameters())
        initial_teacher = teacher_param.data.clone()

        # Get corresponding student parameter
        student_param = next(model.student_encoder.parameters())

        # Modify student parameter
        with torch.no_grad():
            student_param.data.fill_(1.0)

        momentum = 0.99

        # Expected: momentum * initial + (1 - momentum) * student
        expected = momentum * initial_teacher + (1 - momentum) * student_param.data

        # Apply EMA update
        model.update_teacher(momentum)

        # Check teacher was updated correctly
        updated_teacher = teacher_param.data

        assert torch.allclose(updated_teacher, expected, atol=1e-6), \
            "Teacher EMA update incorrect"

    def test_dino_teacher_frozen(self):
        """Test that teacher parameters do not require gradients."""
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
        )

        # Check all teacher encoder parameters
        for name, param in model.teacher_encoder.named_parameters():
            assert not param.requires_grad, f"Teacher encoder param {name} should be frozen"

        # Check all teacher head parameters
        for name, param in model.teacher_head.named_parameters():
            assert not param.requires_grad, f"Teacher head param {name} should be frozen"

    def test_dino_center_update(self):
        """Test that center is updated via EMA of teacher outputs."""
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
            out_dim=4096,
            center_momentum=0.9,
        )

        # Initial center should be zeros
        assert torch.allclose(model.center, torch.zeros(1, 4096)), \
            "Initial center should be zeros"

        # Create mock teacher output
        teacher_output = torch.randn(8, 4096)
        batch_mean = teacher_output.mean(dim=0, keepdim=True)

        # Expected: momentum * 0 + (1 - momentum) * batch_mean
        # With initial center = 0: expected = 0.1 * batch_mean
        expected = (1 - 0.9) * batch_mean

        model.update_center(teacher_output)

        assert torch.allclose(model.center, expected, atol=1e-6), \
            "Center EMA update incorrect"

    def test_dino_student_trainable(self):
        """Test that student parameters require gradients."""
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
        )

        # Check all student encoder parameters
        for name, param in model.student_encoder.named_parameters():
            assert param.requires_grad, f"Student encoder param {name} should be trainable"

        # Check all student head parameters
        for name, param in model.student_head.named_parameters():
            assert param.requires_grad, f"Student head param {name} should be trainable"

    def test_dino_loss_backprop(self):
        """Test that loss backward propagates to student parameters."""
        from CITEgeist.model.dino import DINO

        model = DINO(
            in_channels=2,
            embed_dim=384,
            encoder_depth=12,
            encoder_num_heads=6,
        )
        model.train()

        batch_size = 2
        global_crops = [
            torch.randn(batch_size, 2, 96, 96),
            torch.randn(batch_size, 2, 96, 96),
        ]

        loss = model(global_crops)
        loss.backward()

        # Check student encoder has gradients
        student_grads = 0
        for param in model.student_encoder.parameters():
            if param.grad is not None:
                student_grads += 1

        assert student_grads > 0, "Student encoder should have gradients after backward"

        # Check student head has gradients
        head_grads = 0
        for param in model.student_head.parameters():
            if param.grad is not None:
                head_grads += 1

        assert head_grads > 0, "Student head should have gradients after backward"


class TestDINOFactory:
    """Tests for DINO factory function."""

    def test_create_dino_model_defaults(self):
        """Test create_dino_model with default parameters."""
        from CITEgeist.model.dino import create_dino_model

        model = create_dino_model()

        # Check default configuration
        assert model.student_encoder.embed_dim == 384
        assert model.out_dim == 4096
        assert model.teacher_temp == 0.04
        assert model.student_temp == 0.1

    def test_create_dino_model_custom(self):
        """Test create_dino_model with custom parameters."""
        from CITEgeist.model.dino import create_dino_model

        model = create_dino_model(
            in_channels=2,
            embed_dim=256,
            encoder_depth=6,
            encoder_num_heads=4,
            out_dim=2048,
            teacher_temp=0.05,
        )

        assert model.student_encoder.embed_dim == 256
        assert model.out_dim == 2048
        assert model.teacher_temp == 0.05


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
