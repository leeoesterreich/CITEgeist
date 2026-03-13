"""Protein-Conditioned Multiple Instance Learning (PC-MIL) model.

Fuses spot-level protein proportions with per-nucleus image features
via per-class gated attention to produce single-cell type assignments.

Architecture:
    image_features (N, 384) -> projection -> (N, 64)
    protein_props (K,) -> encoder -> (32,) -> broadcast to N nuclei
    concat -> (N, 96) -> per-class gated attention -> (N, K) logits
    softmax(dim=1) -> attention -> mean(dim=0) -> proportions (K,)
    proportions @ profile_matrix -> reconstructed protein (M,)

See: docs/superpowers/specs/2026-03-13-protein-conditioned-mil-design.md
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Dict, List, Optional, Tuple


def build_profile_matrix(
    cell_profile_dict: Dict[str, List[str]],
    marker_names: List[str],
) -> np.ndarray:
    """Build (K, M) binary profile matrix from cell_profile_dict.

    Args:
        cell_profile_dict: {cell_type: [marker_name, ...]}
        marker_names: Ordered list of all marker names (defines column order)

    Returns:
        (K, M) numpy array where profile[k, m] = 1.0 if marker m belongs to type k
    """
    cell_types = list(cell_profile_dict.keys())
    K = len(cell_types)
    M = len(marker_names)
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    profile = np.zeros((K, M), dtype=np.float32)
    for k, ct in enumerate(cell_types):
        for marker in cell_profile_dict[ct]:
            if marker in marker_to_idx:
                profile[k, marker_to_idx[marker]] = 1.0

    return profile


class PCMILModel(nn.Module):
    """Protein-Conditioned MIL for single-cell type assignment.

    Args:
        image_dim: Input dimension from frozen ViT-S (384)
        n_types: Number of cell types K
        n_markers: Number of protein markers M
        image_proj_dim: Image projection output dimension (64)
        protein_context_dim: Protein encoder output dimension (32)
        hidden_dim: Hidden dimension in gated attention (128)
        dropout: Dropout rate in attention heads
        init_profile_matrix: Optional (K, M) tensor to initialize protein profiles
    """

    def __init__(
        self,
        image_dim: int = 384,
        n_types: int = 7,
        n_markers: int = 15,
        image_proj_dim: int = 64,
        protein_context_dim: int = 32,
        hidden_dim: int = 128,
        dropout: float = 0.1,
        init_profile_matrix: Optional[torch.Tensor] = None,
    ):
        super().__init__()
        self.n_types = n_types
        self.n_markers = n_markers

        # Image projection: 384 -> 64
        self.image_projection = nn.Sequential(
            nn.Linear(image_dim, image_proj_dim),
            nn.LayerNorm(image_proj_dim),
        )

        # Protein encoder: K -> 32
        self.protein_encoder = nn.Sequential(
            nn.Linear(n_types, protein_context_dim),
            nn.ReLU(),
            nn.Linear(protein_context_dim, protein_context_dim),
        )

        # Fused dimension
        fused_dim = image_proj_dim + protein_context_dim

        # Per-class gated attention: K separate gate+score heads
        self.W_gate = nn.ModuleList([
            nn.Linear(fused_dim, hidden_dim) for _ in range(n_types)
        ])
        self.W_score = nn.ModuleList([
            nn.Linear(fused_dim, hidden_dim) for _ in range(n_types)
        ])
        self.W_out = nn.ModuleList([
            nn.Sequential(nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
            for _ in range(n_types)
        ])

        # Protein reconstruction profile matrix: (K, M) learnable
        if init_profile_matrix is not None:
            self.protein_profiles = nn.Parameter(init_profile_matrix.clone())
        else:
            self.protein_profiles = nn.Parameter(torch.randn(n_types, n_markers) * 0.1)

        # Xavier init for projection and attention weights
        self._init_weights()

    def _init_weights(self):
        """Xavier initialization for linear layers."""
        for module in [self.image_projection, self.protein_encoder]:
            for m in module:
                if isinstance(m, nn.Linear):
                    nn.init.xavier_uniform_(m.weight)
                    nn.init.zeros_(m.bias)
        for k in range(self.n_types):
            nn.init.xavier_uniform_(self.W_gate[k].weight)
            nn.init.zeros_(self.W_gate[k].bias)
            nn.init.xavier_uniform_(self.W_score[k].weight)
            nn.init.zeros_(self.W_score[k].bias)
            for m in self.W_out[k]:
                if isinstance(m, nn.Linear):
                    nn.init.xavier_uniform_(m.weight)
                    nn.init.zeros_(m.bias)

    def forward_with_logits(
        self,
        image_features: torch.Tensor,
        protein_proportions: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass returning pre-softmax logits (for inference masking).

        Args:
            image_features: (N, image_dim) pre-extracted ViT features per nucleus
            protein_proportions: (K,) spot-level protein proportions

        Returns:
            logits: (N, K) pre-softmax logits
            attention: (N, K) per-nucleus type probability matrix
            proportions: (K,) predicted spot-level proportions
            reconstructed: (M,) reconstructed protein signal
        """
        N = image_features.shape[0]

        # Project image features
        image_emb = self.image_projection(image_features)  # (N, 64)

        # Encode protein context and broadcast to all nuclei
        protein_context = self.protein_encoder(protein_proportions)  # (32,)
        protein_context_broadcast = protein_context.unsqueeze(0).expand(N, -1)  # (N, 32)

        # Fuse: concat image + protein
        fused = torch.cat([image_emb, protein_context_broadcast], dim=1)  # (N, 96)

        # Per-class gated attention
        logits_list = []
        for k in range(self.n_types):
            gate_k = torch.sigmoid(self.W_gate[k](fused))    # (N, hidden)
            score_k = torch.tanh(self.W_score[k](fused))     # (N, hidden)
            gated = gate_k * score_k                          # (N, hidden)
            logit_k = self.W_out[k](gated)                    # (N, 1)
            logits_list.append(logit_k)

        logits = torch.cat(logits_list, dim=1)  # (N, K)

        # Softmax over types (dim=1) — CRITICAL: NOT dim=0
        attention = F.softmax(logits, dim=1)  # (N, K)

        # Spot-level proportions
        proportions = attention.mean(dim=0)  # (K,)

        # Protein reconstruction
        reconstructed = proportions @ self.protein_profiles  # (M,)

        return logits, attention, proportions, reconstructed

    def forward(
        self,
        image_features: torch.Tensor,
        protein_proportions: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass for a single spot.

        Args:
            image_features: (N, image_dim) pre-extracted ViT features per nucleus
            protein_proportions: (K,) spot-level protein proportions

        Returns:
            proportions: (K,) predicted spot-level proportions
            attention: (N, K) per-nucleus type probability matrix
            reconstructed: (M,) reconstructed protein signal
        """
        _, attention, proportions, reconstructed = self.forward_with_logits(
            image_features, protein_proportions,
        )
        return proportions, attention, reconstructed
