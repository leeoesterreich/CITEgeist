"""
DINO (Self-Distillation with No Labels) for nucleus morphology SSL.

This module implements DINO self-supervised learning for training the ViT encoder
on nucleus morphology patches. DINO uses a teacher-student framework where:
- Student: trained to match teacher outputs via cross-entropy
- Teacher: exponential moving average (EMA) of student
- Centering: prevents mode collapse by centering teacher outputs

Multi-crop training:
- Global crops (2x): Large crops (e.g., 0.4-1.0 of image) processed by both networks
- Local crops (6-8x): Small crops (e.g., 0.05-0.4) processed by student only
- Loss: Cross-entropy from all student outputs to all teacher (global) outputs

Reference:
    Caron et al. (2021) "Emerging Properties in Self-Supervised Vision Transformers"
    https://arxiv.org/abs/2104.14294
"""

import copy
from typing import List, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from .vit_encoder import ViTEncoder
except ImportError:
    from vit_encoder import ViTEncoder


class DINOHead(nn.Module):
    """DINO projection head with bottleneck architecture.

    The head projects encoder outputs to a high-dimensional space for
    self-distillation. Architecture:
        in_dim -> hidden_dim -> bottleneck_dim (L2-norm) -> out_dim

    Key design choices:
    - GELU activation (smoother than ReLU)
    - L2-normalization at bottleneck (stabilizes training)
    - Weight normalization on last layer (prevents collapse)

    Args:
        in_dim: Input dimension from encoder (e.g., 384 for ViT-Small).
        hidden_dim: Hidden layer dimension. Default: 2048.
        bottleneck_dim: Bottleneck dimension (L2-normalized). Default: 256.
        out_dim: Output dimension. Default: 4096.
    """

    def __init__(
        self,
        in_dim: int = 384,
        hidden_dim: int = 2048,
        bottleneck_dim: int = 256,
        out_dim: int = 4096,
    ):
        super().__init__()
        self.in_dim = in_dim
        self.hidden_dim = hidden_dim
        self.bottleneck_dim = bottleneck_dim
        self.out_dim = out_dim

        # First layer: in_dim -> hidden_dim
        self.fc1 = nn.Linear(in_dim, hidden_dim)
        self.act1 = nn.GELU()

        # Second layer: hidden_dim -> bottleneck_dim
        self.fc2 = nn.Linear(hidden_dim, bottleneck_dim)
        self.act2 = nn.GELU()

        # Final layer: bottleneck_dim -> out_dim (with weight normalization)
        self.fc3 = nn.Linear(bottleneck_dim, out_dim, bias=False)
        # Apply weight normalization to prevent collapse
        # Use parametrizations API for compatibility with deepcopy
        self.fc3 = nn.utils.parametrizations.weight_norm(self.fc3)

        self._init_weights()

    def _init_weights(self):
        """Initialize weights with truncated normal distribution."""
        for m in [self.fc1, self.fc2]:
            nn.init.trunc_normal_(m.weight, std=0.02)
            if m.bias is not None:
                nn.init.zeros_(m.bias)
        # fc3 has parametrized weight_norm, initialize the original weight
        nn.init.trunc_normal_(self.fc3.parametrizations.weight.original0, std=0.02)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass through projection head.

        Args:
            x: Input tensor of shape (B, in_dim).

        Returns:
            Output tensor of shape (B, out_dim).
        """
        # First layer
        x = self.fc1(x)
        x = self.act1(x)

        # Second layer with L2-normalization at bottleneck
        x = self.fc2(x)
        x = self.act2(x)
        x = F.normalize(x, dim=-1, p=2)

        # Final layer
        x = self.fc3(x)

        return x


class DINO(nn.Module):
    """DINO self-supervised learning framework for nucleus morphology.

    DINO trains a ViT encoder using self-distillation between a student and
    teacher network. The teacher is an exponential moving average of the student,
    and centering prevents mode collapse.

    Training procedure:
    1. Generate multi-crop augmentations (global + local crops)
    2. Pass global crops through both student and teacher
    3. Pass local crops through student only
    4. Compute cross-entropy loss between all student outputs and teacher outputs
    5. Update student via backprop, teacher via EMA

    Args:
        in_channels: Number of input channels. Default: 2 (DAPI + boundary).
        img_size: Input image size. Default: 96 (nucleus patches).
        patch_size: Size of each patch for ViT. Default: 16.
        embed_dim: Encoder embedding dimension. Default: 384.
        encoder_depth: Number of transformer blocks. Default: 12.
        encoder_num_heads: Number of attention heads. Default: 6.
        mlp_ratio: MLP hidden dim ratio in transformer blocks. Default: 4.0.
        out_dim: Projection head output dimension. Default: 4096.
        hidden_dim: Projection head hidden dimension. Default: 2048.
        bottleneck_dim: Projection head bottleneck dimension. Default: 256.
        teacher_temp: Temperature for teacher softmax. Default: 0.04.
        student_temp: Temperature for student softmax. Default: 0.1.
        center_momentum: Momentum for center EMA update. Default: 0.9.
    """

    def __init__(
        self,
        in_channels: int = 2,
        img_size: int = 96,
        patch_size: int = 16,
        embed_dim: int = 384,
        encoder_depth: int = 12,
        encoder_num_heads: int = 6,
        mlp_ratio: float = 4.0,
        out_dim: int = 4096,
        hidden_dim: int = 2048,
        bottleneck_dim: int = 256,
        teacher_temp: float = 0.04,
        student_temp: float = 0.1,
        center_momentum: float = 0.9,
    ):
        super().__init__()

        # Store hyperparameters
        self.out_dim = out_dim
        self.teacher_temp = teacher_temp
        self.student_temp = student_temp
        self.center_momentum = center_momentum

        # Create student encoder and head
        self.student_encoder = ViTEncoder(
            img_size=img_size,
            patch_size=patch_size,
            in_chans=in_channels,
            embed_dim=embed_dim,
            depth=encoder_depth,
            num_heads=encoder_num_heads,
            mlp_ratio=mlp_ratio,
        )
        self.student_head = DINOHead(
            in_dim=embed_dim,
            hidden_dim=hidden_dim,
            bottleneck_dim=bottleneck_dim,
            out_dim=out_dim,
        )

        # Create teacher encoder and head (copies of student, no gradients)
        self.teacher_encoder = copy.deepcopy(self.student_encoder)
        self.teacher_head = copy.deepcopy(self.student_head)

        # Freeze teacher parameters
        for param in self.teacher_encoder.parameters():
            param.requires_grad = False
        for param in self.teacher_head.parameters():
            param.requires_grad = False

        # Center for preventing collapse (registered as buffer for state_dict)
        self.register_buffer("center", torch.zeros(1, out_dim))

    @torch.no_grad()
    def update_teacher(self, momentum: float) -> None:
        """Update teacher parameters via exponential moving average.

        teacher = momentum * teacher + (1 - momentum) * student

        Args:
            momentum: EMA momentum (typically 0.996-0.9999, increases over training).
        """
        # Update encoder parameters
        for teacher_param, student_param in zip(
            self.teacher_encoder.parameters(),
            self.student_encoder.parameters(),
        ):
            teacher_param.data.mul_(momentum).add_(
                student_param.data, alpha=1 - momentum
            )

        # Update head parameters
        for teacher_param, student_param in zip(
            self.teacher_head.parameters(),
            self.student_head.parameters(),
        ):
            teacher_param.data.mul_(momentum).add_(
                student_param.data, alpha=1 - momentum
            )

    @torch.no_grad()
    def update_center(self, teacher_output: torch.Tensor) -> None:
        """Update center via exponential moving average of teacher outputs.

        center = center_momentum * center + (1 - center_momentum) * mean(teacher_output)

        Centering prevents mode collapse where all outputs converge to same point.

        Args:
            teacher_output: Teacher network outputs of shape (B, out_dim).
        """
        batch_center = teacher_output.mean(dim=0, keepdim=True)
        self.center = (
            self.center_momentum * self.center
            + (1 - self.center_momentum) * batch_center
        )

    def forward_student(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass through student network.

        Args:
            x: Input images of shape (B, C, H, W).

        Returns:
            Student outputs of shape (B, out_dim).
        """
        features = self.student_encoder(x)  # (B, embed_dim)
        output = self.student_head(features)  # (B, out_dim)
        return output

    @torch.no_grad()
    def forward_teacher(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass through teacher network (no gradients).

        Args:
            x: Input images of shape (B, C, H, W).

        Returns:
            Teacher outputs of shape (B, out_dim), centered and temperature-scaled.
        """
        features = self.teacher_encoder(x)  # (B, embed_dim)
        output = self.teacher_head(features)  # (B, out_dim)

        # Center and apply temperature
        output = output - self.center
        output = F.softmax(output / self.teacher_temp, dim=-1)

        return output

    def forward(
        self,
        global_crops: List[torch.Tensor],
        local_crops: Optional[List[torch.Tensor]] = None,
    ) -> torch.Tensor:
        """Compute DINO loss for multi-crop training.

        The loss is computed as cross-entropy between:
        - All student outputs (global + local crops)
        - All teacher outputs (global crops only)

        For each teacher output t and student output s (from different crops):
            loss += -sum(t * log(s))  # Cross-entropy

        Args:
            global_crops: List of 2 global crop tensors, each (B, C, H, W).
            local_crops: Optional list of local crop tensors, each (B, C, H, W).

        Returns:
            Scalar loss tensor.
        """
        if local_crops is None:
            local_crops = []

        # Process global crops through both networks
        teacher_outputs = []
        student_outputs = []

        for crop in global_crops:
            # Teacher forward (no gradients, updates center)
            t_out = self.forward_teacher(crop)
            teacher_outputs.append(t_out)

            # Student forward
            s_out = self.forward_student(crop)
            student_outputs.append(s_out)

        # Process local crops through student only
        for crop in local_crops:
            s_out = self.forward_student(crop)
            student_outputs.append(s_out)

        # Update center using teacher outputs
        all_teacher = torch.cat(teacher_outputs, dim=0)
        self.update_center(all_teacher)

        # Compute cross-entropy loss
        # For each teacher output, compute loss against all student outputs
        # (except the one from the same crop)
        loss = 0.0
        n_loss_terms = 0

        n_global = len(global_crops)

        for i, t_out in enumerate(teacher_outputs):
            for j, s_out in enumerate(student_outputs):
                # Skip same global crop (teacher i matches student i for global crops)
                if j < n_global and i == j:
                    continue

                # Apply temperature to student output and compute cross-entropy
                s_prob = F.log_softmax(s_out / self.student_temp, dim=-1)
                # Cross-entropy: -sum(t * log(s))
                loss_term = -torch.sum(t_out * s_prob, dim=-1).mean()
                loss += loss_term
                n_loss_terms += 1

        # Average over all loss terms
        if n_loss_terms > 0:
            loss = loss / n_loss_terms

        return loss

    def encode(self, imgs: torch.Tensor) -> torch.Tensor:
        """Encode images using the student encoder (for inference).

        Returns the [CLS] token embedding from the student encoder,
        suitable for downstream tasks.

        Args:
            imgs: Input images of shape (B, C, H, W).

        Returns:
            Embeddings of shape (B, embed_dim).
        """
        self.student_encoder.eval()
        with torch.no_grad():
            embeddings = self.student_encoder(imgs)
        return embeddings


def create_dino_model(
    in_channels: int = 2,
    img_size: int = 96,
    patch_size: int = 16,
    embed_dim: int = 384,
    encoder_depth: int = 12,
    encoder_num_heads: int = 6,
    mlp_ratio: float = 4.0,
    **kwargs,
) -> DINO:
    """Factory function to create a DINO model for nucleus morphology.

    This is the recommended way to instantiate DINO for training on
    2-channel (DAPI + boundary) nucleus patches.

    Args:
        in_channels: Number of input channels. Default: 2.
        img_size: Input image size. Default: 96 (nucleus patches).
        patch_size: Size of each patch for ViT. Default: 16.
        embed_dim: Encoder embedding dimension. Default: 384.
        encoder_depth: Number of transformer blocks. Default: 12.
        encoder_num_heads: Number of attention heads. Default: 6.
        mlp_ratio: MLP hidden dim ratio in transformer blocks. Default: 4.0.
        **kwargs: Additional arguments passed to DINO.

    Returns:
        DINO model configured for nucleus morphology SSL.
    """
    return DINO(
        in_channels=in_channels,
        img_size=img_size,
        patch_size=patch_size,
        embed_dim=embed_dim,
        encoder_depth=encoder_depth,
        encoder_num_heads=encoder_num_heads,
        mlp_ratio=mlp_ratio,
        **kwargs,
    )
