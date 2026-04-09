"""
Prototype-Contrastive LLP: Dual-head cell typing from spot-level proportion labels.

Implements a Learning from Label Proportions (LLP) model with:
- Prototype head: cosine-distance assignment to learnable prototypes
- Classifier head: MLP-based classification
- Variance-covariance regularization
- Two-stage training: proportion matching → sharpening
"""

import copy
import logging
import os
import sys
from typing import Dict, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# Lazy import to avoid circular deps; resolved at class instantiation time
def _get_vit_encoder():
    """Import ViTEncoder from the same model directory."""
    _model_dir = os.path.dirname(os.path.abspath(__file__))
    if _model_dir not in sys.path:
        sys.path.insert(0, _model_dir)
    import importlib
    vit_mod = importlib.import_module("vit_encoder")
    return vit_mod.ViTEncoder

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Loss functions
# ---------------------------------------------------------------------------


def proportion_kl_loss(
    oracle: torch.Tensor,
    predicted: torch.Tensor,
    eps: float = 1e-8,
    type_weights: Optional[torch.Tensor] = None,
) -> torch.Tensor:
    """KL divergence loss between oracle and predicted spot-level proportions.

    Computes KL(oracle || predicted) averaged over spots, with optional
    per-type weighting to upweight rare cell types.

    Args:
        oracle: Ground-truth proportion distributions, shape (S, T).
        predicted: Model-predicted proportion distributions, shape (S, T).
        eps: Small constant for numerical stability.
        type_weights: Optional (T,) tensor of per-type weights.
            If provided, each type's KL contribution is scaled by its weight.
            Weights are normalized to sum to T (so unweighted mean is default).

    Returns:
        Scalar tensor: mean KL divergence over spots.
    """
    oracle = oracle.clamp(min=eps)
    predicted = predicted.clamp(min=eps)
    # Normalise to valid distributions
    oracle = oracle / oracle.sum(dim=-1, keepdim=True)
    predicted = predicted / predicted.sum(dim=-1, keepdim=True)
    # KL(oracle || predicted) = Σ oracle * log(oracle / predicted)
    kl_per_type = oracle * (oracle.log() - predicted.log())  # (S, T)
    if type_weights is not None:
        # Normalize weights to sum to T (preserve scale)
        w = type_weights / type_weights.sum() * type_weights.shape[0]
        kl_per_type = kl_per_type * w.unsqueeze(0)
    kl = kl_per_type.sum(dim=-1)  # (S,)
    return kl.mean()


def consistency_loss(
    q_class: torch.Tensor,
    q_proto: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Symmetric KL consistency loss between classifier and prototype heads.

    Uses stop-gradient on the target distribution in each direction:
        (1 / (2N)) * [KL(q_class || sg(q_proto)) + KL(q_proto || sg(q_class))]

    Args:
        q_class: Classifier-head soft assignments, shape (N, T).
        q_proto: Prototype-head soft assignments, shape (N, T).
        eps: Small constant for numerical stability.

    Returns:
        Scalar tensor: symmetric KL consistency loss.
    """
    q_class = q_class.clamp(min=eps)
    q_proto = q_proto.clamp(min=eps)
    q_class = q_class / q_class.sum(dim=-1, keepdim=True)
    q_proto = q_proto / q_proto.sum(dim=-1, keepdim=True)

    # KL(q_class || sg(q_proto))
    kl_cp = (q_class * (q_class.log() - q_proto.detach().log())).sum(dim=-1)
    # KL(q_proto || sg(q_class))
    kl_pc = (q_proto * (q_proto.log() - q_class.detach().log())).sum(dim=-1)

    return (kl_cp + kl_pc).mean() / 2.0


def variance_covariance_loss(
    z: torch.Tensor,
    gamma: float = 1.0,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """VICReg-style variance and covariance regularisation on embeddings.

    Applied on UNNORMALIZED embeddings z before L2 normalisation.

    Variance term encourages each dimension to have std >= gamma:
        var_loss = (1/D) Σ_d max(0, γ − std(Z[:,d]))

    Covariance term penalises off-diagonal covariance:
        cov_loss = (1/D²) Σ_{d≠d'} cov(Z[:,d], Z[:,d'])²

    Args:
        z: Unnormalized embeddings, shape (N, D).
        gamma: Target standard deviation floor (default 1.0).

    Returns:
        Tuple of (var_loss, cov_loss) scalar tensors.
    """
    N, D = z.shape
    # Center embeddings
    z_centered = z - z.mean(dim=0, keepdim=True)

    # Variance loss
    std = z_centered.std(dim=0)  # (D,)
    var_loss = F.relu(gamma - std).mean()

    # Covariance loss
    cov = (z_centered.T @ z_centered) / (N - 1)  # (D, D)
    # Zero out diagonal, penalise squared off-diagonal entries
    mask = ~torch.eye(D, dtype=torch.bool, device=z.device)
    cov_loss = cov[mask].pow(2).sum() / (D * D)

    return var_loss, cov_loss


def separation_loss(
    prototypes: torch.Tensor,
    margin: float = 0.1,
) -> torch.Tensor:
    """Margin-based cosine hinge loss that pushes prototypes apart.

    Encourages pairwise cosine similarity between L2-normalised prototypes
    to stay below `margin`:
        L_sep = mean of max(0, cos(μ_t, μ_t') − margin)²
                over all upper-triangle pairs (t < t').

    Args:
        prototypes: L2-normalised prototype vectors, shape (T, D).
        margin: Maximum allowed cosine similarity (default 0.1).

    Returns:
        Scalar tensor: mean hinge loss over prototype pairs.
    """
    # Cosine similarity matrix: (T, T)
    sim = prototypes @ prototypes.T  # assumes already L2-normalised

    T = prototypes.shape[0]
    # Upper-triangle indices (i < j)
    i_idx, j_idx = torch.triu_indices(T, T, offset=1, device=prototypes.device)
    pairwise_sim = sim[i_idx, j_idx]  # (num_pairs,)

    loss = F.relu(pairwise_sim - margin).pow(2).mean()
    return loss


def sharpening_loss(
    q: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Mean entropy sharpening loss (minimise to sharpen assignments).

    A lower value indicates sharper (more confident) assignments.
        L_sharp = -(1/N) Σ_i Σ_t q_it * log(q_it)

    Args:
        q: Soft assignment distributions, shape (N, T).
        eps: Small constant for numerical stability.

    Returns:
        Scalar tensor: mean entropy over cells.
    """
    q = q.clamp(min=eps)
    q = q / q.sum(dim=-1, keepdim=True)
    entropy = -(q * q.log()).sum(dim=-1)  # (N,)
    return entropy.mean()


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class PrototypeContrastiveModel(nn.Module):
    """Dual-head prototype-contrastive model for cell-type proportion LLP.

    Architecture:
        ViT-S encoder (384-d) → Projector MLP (128-d) → z (unnormalized)
                                                            │
                                                   ẑ = z / ||z||  (L2-normalized)
                                                            │
                                         ┌──────────────────┼──────────────────┐
                                         │                                      │
                                 Prototype Head                        Classifier Head
                                 cos_sim(ẑ, μ̂) / τ                   MLP(z→64→T) + softmax
                                         │                                      │
                                    q_proto (T)                            q_class (T)

    When ``from_embeddings=True``, the ViT encoder is omitted and the model
    accepts precomputed (B, encoder_embed_dim) feature vectors directly.

    Args:
        num_types: Number of cell types. Default: 6.
        embed_dim: Projector output dimension. Default: 128.
        in_channels: Number of input image channels. Default: 3.
        temperature: Prototype softmax temperature (τ). Modifiable externally. Default: 0.1.
        encoder_depth: Number of ViT transformer blocks. Default: 12.
        encoder_embed_dim: ViT embedding dimension. Default: 384.
        encoder_num_heads: ViT attention heads. Default: 6.
        from_embeddings: If True, skip the ViT encoder and accept precomputed
            (B, encoder_embed_dim) embeddings as input. Default: False.
    """

    def __init__(
        self,
        num_types: int = 6,
        embed_dim: int = 128,
        in_channels: int = 3,
        temperature: float = 0.1,
        encoder_depth: int = 12,
        encoder_embed_dim: int = 384,
        encoder_num_heads: int = 6,
        encoder_img_size: int = 96,
        from_embeddings: bool = False,
    ):
        super().__init__()
        self.num_types = num_types
        self.embed_dim = embed_dim
        self.temperature = temperature
        self.from_embeddings = from_embeddings
        self._encoder_embed_dim = encoder_embed_dim

        if from_embeddings:
            # No ViT encoder — caller supplies precomputed feature vectors
            self.encoder = None
        else:
            # ViT-S encoder
            ViTEncoder = _get_vit_encoder()
            self.encoder = ViTEncoder(
                img_size=encoder_img_size,
                patch_size=16,
                in_chans=in_channels,
                embed_dim=encoder_embed_dim,
                depth=encoder_depth,
                num_heads=encoder_num_heads,
            )

        # Projector MLP: 384 → 128
        self.projector = nn.Sequential(
            nn.Linear(encoder_embed_dim, embed_dim),
            nn.LayerNorm(embed_dim),
            nn.ReLU(inplace=True),
            nn.Linear(embed_dim, embed_dim),
        )

        # Learnable prototype matrix (T, D) — kept unnormalized; normalized at use
        # Use a temporary tensor for xavier init to avoid in-place ops on a leaf Parameter
        _proto_init = torch.empty(1, num_types, embed_dim)
        nn.init.xavier_uniform_(_proto_init)
        self.prototypes = nn.Parameter(_proto_init.squeeze(0))

        # Classifier head: z (D) → 64 → T
        self.classifier_head = nn.Sequential(
            nn.Linear(embed_dim, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, num_types),
        )

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    def forward(self, patches: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Forward pass.

        Args:
            patches: Either (B, C, H, W) image patches when ``from_embeddings``
                is False, or (B, encoder_embed_dim) precomputed feature vectors
                when ``from_embeddings`` is True.

        Returns:
            Dictionary with keys:
                - ``q_class``: (B, T) softmax class probabilities from classifier head.
                - ``q_proto``: (B, T) softmax prototype assignment scores.
                - ``z``: (B, D) unnormalized projector embeddings.
                - ``z_norm``: (B, D) L2-normalized embeddings.
        """
        # Encoder → projector
        if self.from_embeddings:
            h = patches                     # (B, encoder_embed_dim) — already encoded
        else:
            h = self.encoder(patches)       # (B, encoder_embed_dim)
        z = self.projector(h)               # (B, embed_dim)
        z_norm = F.normalize(z, dim=-1)     # (B, embed_dim)

        # Prototype head: cos_sim(z_norm, proto_norm) / τ
        proto_norm = self.get_normalized_prototypes()   # (T, D)
        logits_proto = (z_norm @ proto_norm.T) / self.temperature  # (B, T)
        q_proto = F.softmax(logits_proto, dim=-1)       # (B, T)

        # Classifier head (operates on unnormalized z)
        logits_class = self.classifier_head(z)          # (B, T)
        q_class = F.softmax(logits_class, dim=-1)       # (B, T)

        return {
            "q_class": q_class,
            "q_proto": q_proto,
            "z": z,
            "z_norm": z_norm,
        }

    # ------------------------------------------------------------------
    # Prototype utilities
    # ------------------------------------------------------------------

    def get_normalized_prototypes(self) -> torch.Tensor:
        """Return L2-normalized prototype matrix.

        Returns:
            Tensor of shape (T, D) with unit-norm rows.
        """
        return F.normalize(self.prototypes, dim=-1)

    # ------------------------------------------------------------------
    # Encoder freeze / unfreeze helpers
    # ------------------------------------------------------------------

    def freeze_encoder(self) -> None:
        """Freeze all encoder parameters (set requires_grad=False).

        No-op when ``from_embeddings=True`` (no encoder present).
        """
        if self.encoder is None:
            return
        for param in self.encoder.parameters():
            param.requires_grad = False

    def unfreeze_last_n_blocks(self, n: int) -> None:
        """Unfreeze the last n transformer blocks and the final LayerNorm.

        Typically called after ``freeze_encoder()`` to enable partial fine-tuning.
        No-op when ``from_embeddings=True`` (no encoder present).

        Args:
            n: Number of trailing blocks to unfreeze.
        """
        if self.encoder is None:
            return
        blocks = self.encoder.blocks
        total = len(blocks)
        # Unfreeze last n blocks
        for block in blocks[max(0, total - n):]:
            for param in block.parameters():
                param.requires_grad = True
        # Unfreeze final LayerNorm
        for param in self.encoder.norm.parameters():
            param.requires_grad = True

    # ------------------------------------------------------------------
    # Checkpoint loading
    # ------------------------------------------------------------------

    def load_simclr_encoder(
        self,
        checkpoint_path: str,
        device: Optional[str] = None,
    ) -> None:
        """Load encoder weights from a SimCLR checkpoint.

        Handles two common checkpoint formats:
        1. Direct state dict with ``encoder.*`` prefix keys.
        2. Wrapped dict with a ``"model_state_dict"`` key that contains (1).

        No-op when ``from_embeddings=True`` (no encoder present); logs a warning.

        Args:
            checkpoint_path: Path to the ``.pt`` or ``.pth`` checkpoint file.
            device: Device string for loading (e.g. ``"cpu"``). Defaults to
                the current device of the first encoder parameter.
        """
        if self.encoder is None:
            logger.warning(
                "load_simclr_encoder: model was created with from_embeddings=True "
                "(no encoder). Checkpoint %s will not be loaded.",
                checkpoint_path,
            )
            return

        if device is None:
            try:
                device = next(self.encoder.parameters()).device
            except StopIteration:
                device = "cpu"

        ckpt = torch.load(checkpoint_path, map_location=device)

        # Unwrap nested checkpoint formats
        if isinstance(ckpt, dict) and "model_state_dict" in ckpt:
            ckpt = ckpt["model_state_dict"]

        # Extract keys with "encoder." prefix
        encoder_state = {}
        for k, v in ckpt.items():
            if k.startswith("encoder."):
                encoder_state[k[len("encoder."):]] = v

        if not encoder_state:
            logger.warning(
                "load_simclr_encoder: no 'encoder.*' keys found in checkpoint %s; "
                "attempting to load full state dict directly.",
                checkpoint_path,
            )
            encoder_state = ckpt

        missing, unexpected = self.encoder.load_state_dict(encoder_state, strict=False)
        if missing:
            logger.warning("load_simclr_encoder: missing keys in encoder: %s", missing)
        if unexpected:
            logger.warning("load_simclr_encoder: unexpected keys in checkpoint: %s", unexpected)
        logger.info("Loaded SimCLR encoder weights from %s", checkpoint_path)

    # ------------------------------------------------------------------
    # Prototype initialisation from K-Means
    # ------------------------------------------------------------------

    def init_prototypes_from_kmeans(
        self,
        embeddings: np.ndarray,
        num_types: Optional[int] = None,
    ) -> None:
        """Initialise prototypes using K-Means on a batch of embeddings.

        Runs sklearn KMeans, L2-normalises the centroids, and copies them
        into ``self.prototypes``.

        Args:
            embeddings: Array of shape (N, D) containing L2-normalised or
                raw embeddings from the projector.
            num_types: Number of clusters. Defaults to ``self.num_types``.
        """
        from sklearn.cluster import KMeans  # local import to keep torch-only tests fast

        k = num_types if num_types is not None else self.num_types
        km = KMeans(n_clusters=k, n_init=10, random_state=0)
        km.fit(embeddings)

        centroids = torch.tensor(km.cluster_centers_, dtype=self.prototypes.dtype)
        centroids = F.normalize(centroids, dim=-1)          # (k, D)

        with torch.no_grad():
            self.prototypes[:k].copy_(centroids)


# ---------------------------------------------------------------------------
# scatter_mean helper (works with torch_scatter OR torch_geometric)
# ---------------------------------------------------------------------------

try:
    from torch_scatter import scatter_mean as _scatter_mean_impl

    def _scatter_mean(src, index, dim=0, dim_size=None):
        return _scatter_mean_impl(src, index, dim=dim, dim_size=dim_size)

except ImportError:
    try:
        from torch_geometric.utils import scatter as _tg_scatter

        def _scatter_mean(src, index, dim=0, dim_size=None):
            return _tg_scatter(src, index, dim=dim, dim_size=dim_size, reduce="mean")

    except ImportError:
        # Pure-PyTorch fallback (no external dependency)
        def _scatter_mean(src, index, dim=0, dim_size=None):
            if dim_size is None:
                dim_size = int(index.max().item()) + 1
            shape = list(src.shape)
            shape[dim] = dim_size
            out = src.new_zeros(shape)
            cnt = src.new_zeros(dim_size)
            # Expand index to match src dims
            idx = index
            for _ in range(src.dim() - 1):
                idx = idx.unsqueeze(-1)
            idx = idx.expand_as(src)
            out.scatter_add_(dim, idx, src)
            cnt.scatter_add_(0, index, src.new_ones(src.shape[0]))
            cnt = cnt.clamp(min=1)
            # Divide out by count along the scatter dim
            cnt_shape = [1] * src.dim()
            cnt_shape[dim] = dim_size
            out = out / cnt.view(cnt_shape)
            return out


# ---------------------------------------------------------------------------
# Task 3: Data Augmentation
# ---------------------------------------------------------------------------

import random  # noqa: E402


class CellPatchAugmentation:
    """Stochastic augmentation pipeline for cell-patch tensors.

    Supports single patches (C, H, W) and batches (B, C, H, W).

    Args:
        continuous_rotation: If True, rotate by a uniform random angle in
            [0, 360°) using bilinear interpolation.  If False, rotate by a
            random multiple of 90°.
        intensity_scale: Half-range for multiplicative intensity jitter.
            Applied as ``x * U(1-s, 1+s)``.  Default: 0.2.
        noise_std: Standard deviation of additive Gaussian noise.  Default: 0.05.
    """

    def __init__(
        self,
        continuous_rotation: bool = True,
        intensity_scale: float = 0.2,
        noise_std: float = 0.05,
    ):
        self.continuous_rotation = continuous_rotation
        self.intensity_scale = intensity_scale
        self.noise_std = noise_std

    def _augment_single(self, x: torch.Tensor) -> torch.Tensor:
        """Augment a single (C, H, W) patch."""
        # Random horizontal flip
        if random.random() < 0.5:
            x = torch.flip(x, dims=[-1])
        # Random vertical flip
        if random.random() < 0.5:
            x = torch.flip(x, dims=[-2])
        # Rotation
        if self.continuous_rotation:
            try:
                import torchvision.transforms.functional as TF  # lazy import (avoids PIL link issues)
                angle = random.uniform(0.0, 360.0)
                x = TF.rotate(x, angle, interpolation=TF.InterpolationMode.BILINEAR)
            except ImportError:
                # Fall back to 90° multiples if torchvision/PIL unavailable
                k = random.randint(0, 3)
                x = torch.rot90(x, k=k, dims=[-2, -1])
        else:
            k = random.randint(0, 3)
            x = torch.rot90(x, k=k, dims=[-2, -1])
        # Intensity scaling
        if self.intensity_scale > 0:
            scale = 1.0 + random.uniform(-self.intensity_scale, self.intensity_scale)
            x = x * scale
        # Gaussian noise
        if self.noise_std > 0:
            x = x + torch.randn_like(x) * self.noise_std
        return x

    def __call__(self, x: torch.Tensor) -> torch.Tensor:
        """Apply augmentation.

        Args:
            x: (C, H, W) single patch or (B, C, H, W) batch.

        Returns:
            Augmented tensor with the same shape as input.
        """
        if x.dim() == 3:
            return self._augment_single(x)
        elif x.dim() == 4:
            return torch.stack([self._augment_single(xi) for xi in x], dim=0)
        else:
            raise ValueError(f"Expected 3-D or 4-D tensor, got {x.dim()}-D")


def tta_8x(patch: torch.Tensor) -> list:
    """Generate 8 deterministic test-time augmentation views of a patch.

    Produces 4 rotations (0°, 90°, 180°, 270°) × 2 flips (none, horizontal).

    Args:
        patch: Single patch tensor of shape (C, H, W).

    Returns:
        List of 8 tensors each of shape (C, H, W).
    """
    views = []
    for k in range(4):
        rotated = torch.rot90(patch, k=k, dims=[-2, -1])
        views.append(rotated)
        views.append(torch.flip(rotated, dims=[-1]))
    return views


# ---------------------------------------------------------------------------
# Task 4: Training loop helpers
# ---------------------------------------------------------------------------


def _spatial_train_val_split(
    spot_coords: Optional[np.ndarray],
    num_spots: int,
    val_frac: float = 0.2,
    buffer_um: float = 150.0,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """Split spots into train and validation sets.

    Args:
        spot_coords: (S, 2) array of spot coordinates in microns, or None
            for a random split.
        num_spots: Total number of spots.
        val_frac: Fraction of spots to hold out for validation.
        buffer_um: Exclusion buffer (microns) around the split boundary.
        seed: Random seed.

    Returns:
        Tuple (train_spot_indices, val_spot_indices) as numpy arrays.
    """
    rng = np.random.RandomState(seed)
    all_idx = np.arange(num_spots)

    if spot_coords is None or len(spot_coords) == 0:
        rng.shuffle(all_idx)
        n_val = max(1, int(num_spots * val_frac))
        val_idx = all_idx[:n_val]
        train_idx = all_idx[n_val:]
        return train_idx, val_idx

    # Spatial split along a random axis
    axis = rng.randint(0, 2)
    coords_1d = spot_coords[:, axis]
    threshold = np.quantile(coords_1d, 1.0 - val_frac)

    val_mask = coords_1d > threshold
    # Buffer: exclude spots within ±buffer_um of threshold
    in_buffer = np.abs(coords_1d - threshold) < buffer_um
    val_mask_buffered = val_mask & ~in_buffer
    train_mask = ~val_mask & ~in_buffer

    train_idx = all_idx[train_mask]
    val_idx = all_idx[val_mask_buffered]

    # Fallback: if either split is empty, do random split
    if len(train_idx) == 0 or len(val_idx) == 0:
        rng.shuffle(all_idx)
        n_val = max(1, int(num_spots * val_frac))
        val_idx = all_idx[:n_val]
        train_idx = all_idx[n_val:]

    return train_idx, val_idx


@torch.no_grad()
def _inference_all(
    model: "PrototypeContrastiveModel",
    patches: torch.Tensor,
    device: str,
    batch_size: int = 256,
    backbone: Optional[nn.Module] = None,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Run inference over all patches in mini-batches.

    Args:
        model: The PrototypeContrastiveModel (in eval mode).
        patches: (C, ch, 96, 96) all cell patches, or (C, 3, 224, 224) raw
            images when ``backbone`` is provided.
        device: Device string.
        batch_size: Mini-batch size.
        backbone: Optional external backbone (e.g. CTransPath). When provided,
            each mini-batch is passed through ``backbone.extract()`` to produce
            embeddings before being fed to ``model``.

    Returns:
        Tuple (q_class, q_proto, z) each on ``device``.
    """
    model.eval()
    q_class_list, q_proto_list, z_list = [], [], []
    N = patches.shape[0]
    for start in range(0, N, batch_size):
        batch = patches[start: start + batch_size].to(device)
        if backbone is not None:
            batch = backbone.extract(batch)  # no_grad already active from decorator
        out = model(batch)
        q_class_list.append(out["q_class"])
        q_proto_list.append(out["q_proto"])
        z_list.append(out["z"])
    return (
        torch.cat(q_class_list, dim=0),
        torch.cat(q_proto_list, dim=0),
        torch.cat(z_list, dim=0),
    )


@torch.no_grad()
def _validate(
    model: "PrototypeContrastiveModel",
    patches: torch.Tensor,
    c2s: torch.Tensor,
    oracle_props: torch.Tensor,
    val_spots: np.ndarray,
    device: str,
    batch_size: int = 256,
    backbone: Optional[nn.Module] = None,
) -> float:
    """Compute proportion KL loss on validation spots.

    Args:
        model: Model to evaluate.
        patches: (C, ch, 96, 96) all patches, or (C, 3, 224, 224) raw images
            when ``backbone`` is provided.
        c2s: (C,) cell-to-spot mapping.
        oracle_props: (S, T) ground-truth proportions.
        val_spots: Validation spot indices.
        device: Device string.
        batch_size: Mini-batch size for inference.
        backbone: Optional external backbone forwarded to ``_inference_all``.

    Returns:
        Mean proportion KL loss (float).
    """
    model.eval()
    num_spots = oracle_props.shape[0]
    T = oracle_props.shape[1]

    # Get all q_class predictions
    q_class, _, _ = _inference_all(model, patches, device, batch_size=batch_size, backbone=backbone)

    # Scatter mean to spots
    c2s_dev = c2s.to(device)
    pred_props = _scatter_mean(q_class, c2s_dev, dim=0, dim_size=num_spots)  # (S, T)

    # Compute KL only on val spots
    val_t = torch.tensor(val_spots, dtype=torch.long)
    kl = proportion_kl_loss(
        oracle_props[val_t].to(device),
        pred_props[val_t.to(device)],
    )
    return kl.item()


def _refresh_prototypes(
    model: "PrototypeContrastiveModel",
    patches: torch.Tensor,
    c2s: torch.Tensor,
    device: str,
    top_frac: float = 0.1,
    batch_size: int = 256,
    backbone: Optional[nn.Module] = None,
) -> None:
    """Refresh prototype vectors using the most-confident cells per type.

    For each type t, selects the top ``top_frac`` fraction of cells by
    ``q_proto[:, t]`` confidence and averages their L2-normalised embeddings
    to form the new prototype.

    Args:
        model: Model to update in-place.
        patches: (C, ch, 96, 96) all patches, or (C, 3, 224, 224) raw images
            when ``backbone`` is provided.
        c2s: (C,) cell-to-spot mapping (unused, kept for API consistency).
        device: Device string.
        top_frac: Fraction of cells to use per type. Default: 0.1.
        batch_size: Mini-batch size for inference.
        backbone: Optional external backbone forwarded to ``_inference_all``.
    """
    _, q_proto, z = _inference_all(model, patches, device, batch_size=batch_size, backbone=backbone)
    z_norm = F.normalize(z, dim=-1)  # (C, D)

    T = model.num_types
    new_protos = []
    for t in range(T):
        scores = q_proto[:, t]  # (C,)
        k = max(1, int(len(scores) * top_frac))
        top_idx = scores.topk(k).indices
        centroid = z_norm[top_idx].mean(dim=0)
        centroid = F.normalize(centroid, dim=-1)
        new_protos.append(centroid)

    new_proto_mat = torch.stack(new_protos, dim=0)  # (T, D)
    with torch.no_grad():
        model.prototypes.copy_(new_proto_mat)


def _train_one_epoch(
    model: "PrototypeContrastiveModel",
    patches: torch.Tensor,
    c2s: torch.Tensor,
    oracle_props: torch.Tensor,
    train_spots: np.ndarray,
    spots_per_batch: int,
    augmentation: "CellPatchAugmentation",
    optimizer: torch.optim.Optimizer,
    device: str,
    lambda_consist: float,
    lambda_var: float,
    lambda_cov: float,
    lambda_sep: float,
    sep_margin: float,
    lambda_sharp: float,
    type_weights: Optional[torch.Tensor] = None,
    backbone: Optional[nn.Module] = None,
) -> float:
    """Train for one epoch over all training spots.

    Args:
        model: PrototypeContrastiveModel.
        patches: (C, ch, 96, 96) all cell patches (on CPU).
        c2s: (C,) cell-to-spot mapping.
        oracle_props: (S, T) ground-truth proportions.
        train_spots: Array of training spot indices.
        spots_per_batch: Number of spots per mini-batch.
        augmentation: CellPatchAugmentation instance.
        optimizer: Torch optimiser.
        device: Device string.
        lambda_consist: Weight for consistency loss.
        lambda_var: Weight for variance loss.
        lambda_cov: Weight for covariance loss.
        lambda_sep: Weight for prototype separation loss.
        sep_margin: Cosine margin for separation loss.
        lambda_sharp: Weight for sharpening loss (0 in Stage 2A).
        backbone: Optional external encoder (e.g. CTransPathBackbone) to encode
            patches per-batch with gradients. Mutually exclusive with augmentation.

    Returns:
        Mean total loss over mini-batches.
    """
    if backbone is not None and augmentation is not None:
        raise ValueError(
            "backbone and augmentation are mutually exclusive in _train_one_epoch"
        )
    model.train()
    rng = np.random.RandomState()
    shuffled_spots = train_spots.copy()
    rng.shuffle(shuffled_spots)

    total_loss = 0.0
    num_batches = 0

    num_spots = oracle_props.shape[0]

    for batch_start in range(0, len(shuffled_spots), spots_per_batch):
        batch_spot_ids = shuffled_spots[batch_start: batch_start + spots_per_batch]
        batch_spot_set = set(batch_spot_ids.tolist())

        # Find cells belonging to these spots
        c2s_np = c2s.numpy() if c2s.is_cuda else c2s.cpu().numpy()
        cell_mask = np.isin(c2s_np, batch_spot_ids)
        cell_idx = np.where(cell_mask)[0]

        if len(cell_idx) == 0:
            continue

        # Remap global spot IDs to local batch indices for scatter_mean
        spot_id_to_local = {sid: li for li, sid in enumerate(batch_spot_ids)}
        local_c2s = np.array([spot_id_to_local[c2s_np[ci]] for ci in cell_idx])
        local_c2s_t = torch.tensor(local_c2s, dtype=torch.long, device=device)

        # Gather patches/embeddings for these cells and encode/augment
        cell_patches = patches[cell_idx]  # (n, 3, H, W) or (n, embed_dim)
        if backbone is not None:
            # Sub-batch extraction to avoid OOM: spots can have 700+ cells; running
            # all through Swin-Tiny with gradients at once exceeds 39 GB VRAM.
            # torch.cat preserves grad_fn so backward still flows to backbone weights.
            # For unfrozen CTransPath (Condition 5), gradient checkpointing + _SUB=32 needed.
            _SUB = 32
            emb_chunks = [
                backbone.extract(cell_patches[s:s + _SUB].to(device))
                for s in range(0, len(cell_patches), _SUB)
            ]
            aug_patches = torch.cat(emb_chunks, dim=0)  # (n, 768), grads flow
        elif augmentation is not None:
            aug_patches = augmentation(cell_patches).to(device)
        else:
            aug_patches = cell_patches.to(device)

        # Forward
        out = model(aug_patches)
        q_class = out["q_class"]      # (n, T)
        q_proto = out["q_proto"]      # (n, T)
        z = out["z"]                  # (n, D)

        # Scatter mean to local spots
        n_local = len(batch_spot_ids)
        pred_props = _scatter_mean(q_class, local_c2s_t, dim=0, dim_size=n_local)  # (n_local, T)

        # Oracle proportions for this batch (local order)
        oracle_batch = oracle_props[batch_spot_ids].to(device)  # (n_local, T)

        # Proportion KL loss (with optional per-type weighting)
        loss_prop = proportion_kl_loss(oracle_batch, pred_props, type_weights=type_weights)

        # Consistency
        loss_consist = consistency_loss(q_class, q_proto)

        # Variance-covariance
        var_loss, cov_loss = variance_covariance_loss(z)

        # Prototype separation
        proto_norm = model.get_normalized_prototypes()
        loss_sep = separation_loss(proto_norm, margin=sep_margin)

        # Sharpening (Stage 2B only; 0 in 2A via lambda_sharp=0)
        loss_sharp = sharpening_loss(q_class) + sharpening_loss(q_proto)

        loss = (
            loss_prop
            + lambda_consist * loss_consist
            + lambda_var * var_loss
            + lambda_cov * cov_loss
            + lambda_sep * loss_sep
            + lambda_sharp * loss_sharp
        )

        optimizer.zero_grad()
        loss.backward()
        all_params = list(model.parameters())
        if backbone is not None:
            all_params += list(backbone._model.parameters())
        torch.nn.utils.clip_grad_norm_(all_params, max_norm=5.0)
        optimizer.step()

        total_loss += loss.item()
        num_batches += 1

    return total_loss / max(num_batches, 1)


def train_prototype_contrastive(
    patches: torch.Tensor,
    cell_to_spot: torch.Tensor,
    oracle_props: torch.Tensor,
    num_types: int = 6,
    embed_dim: int = 128,
    in_channels: int = 3,
    encoder_depth: int = 12,
    n_epochs_2a: int = 100,
    n_epochs_2b: int = 50,
    encoder_warmup_epochs: int = 10,
    finetune_layers: int = 2,
    spots_per_batch: int = 200,
    lr_2a: float = 1e-4,
    lr_2b: float = 1e-5,
    weight_decay: float = 1e-4,
    tau_start: float = 0.1,
    tau_end: float = 0.05,
    lambda_consist: float = 1.0,
    lambda_var: float = 0.1,
    lambda_cov: float = 0.01,
    lambda_sep: float = 0.1,
    lambda_sharp: float = 0.5,
    sep_margin: float = 0.1,
    val_frac: float = 0.2,
    spot_coords: Optional[np.ndarray] = None,
    device: str = "cuda",
    simclr_checkpoint: Optional[str] = None,
    prototype_refresh: bool = True,
    use_type_weights: bool = False,
    oracle_label_smoothing: float = 0.0,
    encoder_img_size: int = 96,
    augmentation_continuous_rotation: bool = True,
    from_embeddings: bool = False,
    encoder_embed_dim: int = 384,
    external_backbone: Optional[nn.Module] = None,
    warmup_freeze_backbone: bool = True,
    early_stop_patience: int = 20,
    seed: int = 42,
) -> Dict:
    """Train the PrototypeContrastiveModel with two-stage LLP.

    Stage 2A: frozen encoder → warm up heads + prototypes.
    Partway through 2A, unfreeze last ``finetune_layers`` encoder blocks.
    Stage 2B: lower LR, linear temperature + sharpening ramp.

    Args:
        patches: (C, ch, 96, 96) cell patch tensor, or (C, encoder_embed_dim)
            precomputed embeddings when ``from_embeddings=True``.
        cell_to_spot: (C,) integer spot indices.
        oracle_props: (S, T) ground-truth proportion labels.
        num_types: Number of cell types.
        embed_dim: Projector output dimension.
        in_channels: Number of image channels (ignored when from_embeddings=True).
        encoder_depth: Number of ViT transformer blocks (ignored when from_embeddings=True).
        n_epochs_2a: Epochs for Stage 2A.
        n_epochs_2b: Epochs for Stage 2B.
        encoder_warmup_epochs: Epoch at which to unfreeze last N encoder blocks.
            Ignored when from_embeddings=True.
        finetune_layers: Number of encoder blocks to unfreeze.
            Ignored when from_embeddings=True.
        spots_per_batch: Spots per gradient step.
        lr_2a: Learning rate for Stage 2A.
        lr_2b: Learning rate for Stage 2B.
        weight_decay: AdamW weight decay.
        tau_start: Prototype temperature at start of Stage 2B.
        tau_end: Prototype temperature at end of Stage 2B.
        lambda_consist: Consistency loss weight.
        lambda_var: Variance loss weight.
        lambda_cov: Covariance loss weight.
        lambda_sep: Prototype separation loss weight.
        lambda_sharp: Maximum sharpening loss weight (ramped in 2B).
        sep_margin: Cosine margin for separation loss.
        val_frac: Fraction of spots for validation.
        spot_coords: (S, 2) optional spatial coordinates for spatial split.
        device: Torch device string.
        simclr_checkpoint: Path to SimCLR checkpoint to warm-start encoder.
            Ignored when from_embeddings=True.
        prototype_refresh: Whether to refresh prototypes between stages.
        use_type_weights: Whether to apply inverse-frequency type weighting.
        oracle_label_smoothing: If > 0, blend oracle proportions with a uniform
            prior: ``oracle = (1 - alpha) * oracle + alpha / num_types``.
            This creates a nonzero gradient floor for rare cell types whose
            oracle proportion is 0 in most spots. Recommended: 0.05–0.1.
        from_embeddings: If True, ``patches`` contains precomputed (C, encoder_embed_dim)
            feature vectors; the ViT encoder is bypassed and augmentation /
            encoder unfreezing are disabled.
        encoder_embed_dim: Encoder output dimension. Default: 384 (ViT-S).
            Set to 768 when using CTransPath precomputed embeddings.
            Ignored when ``from_embeddings=False`` (ViT encoder sets its own dim).
        external_backbone: Optional external encoder (e.g. CTransPathBackbone with
            frozen=False) to encode raw patches per-batch during training.
            When provided, patches must be (N, 3, H, W) image tensors;
            augmentation is disabled. Backbone is frozen during warmup epochs
            then unfrozen with lr = lr_2a / 10 if warmup_freeze_backbone=True.
        warmup_freeze_backbone: If True (default) and external_backbone is given,
            keep backbone frozen for the first encoder_warmup_epochs epochs, then
            unfreeze with a 10x lower learning rate than the head.
        early_stop_patience: Stop training a stage if val_kl does not improve for
            this many consecutive epochs. Set to 0 to disable. Default 20.
        seed: Random seed.

    Returns:
        Dict with keys:
            - ``model_state``: model state dict.
            - ``train_losses``: list of per-epoch training losses.
            - ``val_losses``: list of per-epoch validation losses.
            - ``spot_proportions``: (S, T) predicted proportions.
            - ``cell_assignments``: (C, T) soft cell assignments (q_class).
            - ``cell_proto_assignments``: (C, T) prototype soft assignments (q_proto).
            - ``embeddings``: (C, D) projector embeddings.
            - ``temperature_final``: final temperature value.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    num_spots = oracle_props.shape[0]

    # ------------------------------------------------------------------
    # Oracle label smoothing — prevents zero-gradient for rare types
    # ------------------------------------------------------------------
    if oracle_label_smoothing > 0.0:
        uniform = torch.ones_like(oracle_props) / num_types
        oracle_props = (1.0 - oracle_label_smoothing) * oracle_props + oracle_label_smoothing * uniform
        logger.info(
            "Oracle label smoothing alpha=%.2f applied; type floors: %s",
            oracle_label_smoothing,
            {f"type_{i}": f"{oracle_props[:, i].min().item():.4f}" for i in range(num_types)},
        )

    # ------------------------------------------------------------------
    # Build model
    # ------------------------------------------------------------------
    model = PrototypeContrastiveModel(
        num_types=num_types,
        embed_dim=embed_dim,
        in_channels=in_channels,
        temperature=tau_start,
        encoder_depth=encoder_depth,
        encoder_embed_dim=encoder_embed_dim,
        encoder_img_size=encoder_img_size,
        from_embeddings=from_embeddings,
    ).to(device)

    if simclr_checkpoint is not None:
        # load_simclr_encoder is a no-op when from_embeddings=True
        model.load_simclr_encoder(simclr_checkpoint, device=device)

    # Freeze encoder initially (no-op when from_embeddings=True)
    model.freeze_encoder()

    # ------------------------------------------------------------------
    # Train / val split
    # ------------------------------------------------------------------
    train_spots, val_spots = _spatial_train_val_split(
        spot_coords, num_spots, val_frac=val_frac, seed=seed
    )

    # Augmentation is spatial (rot/flip) — meaningless for 1D embeddings
    # augmentation_continuous_rotation=False: rot90 (~0s) vs TF.rotate bilinear (slow)
    # When from_embeddings=True augmentation is always None (no image tensors)
    if from_embeddings:
        augmentation = None
    elif augmentation_continuous_rotation:
        augmentation = CellPatchAugmentation(continuous_rotation=True)
    else:
        # Disable noise: torch.randn_like per-patch in Python loop is slow at N=3500+
        # Flip + rot90 only — both are view-ops, near-zero cost
        augmentation = CellPatchAugmentation(
            continuous_rotation=False, noise_std=0.0
        )

    # ------------------------------------------------------------------
    # Per-type weighting (inverse-frequency from oracle proportions)
    # ------------------------------------------------------------------
    type_weights_t = None
    if use_type_weights:
        mean_props = oracle_props.mean(dim=0).clamp(min=1e-4)  # (T,)
        # sqrt inverse-frequency: moderate upweighting without extreme ratios
        type_weights_t = (1.0 / mean_props).sqrt()
        type_weights_t = type_weights_t.to(device)
        logger.info(
            "Per-type weights (inv-freq): %s",
            {f"type_{i}": f"{w:.2f}" for i, w in enumerate(type_weights_t.tolist())},
        )

    # ------------------------------------------------------------------
    # Stage 2A optimizer
    # ------------------------------------------------------------------
    optimizer = torch.optim.AdamW(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=lr_2a,
        weight_decay=weight_decay,
    )

    train_losses: list = []
    val_losses: list = []

    # ------------------------------------------------------------------
    # Stage 2A
    # ------------------------------------------------------------------
    encoder_unfrozen = False
    _es_best_val = float("inf")
    _es_no_improve = 0
    _best_ckpt: Dict = {"val_kl": float("inf"), "model": None, "backbone": None}
    for epoch in range(n_epochs_2a):
        # Unfreeze encoder at warmup epoch (skipped when from_embeddings=True)
        if (
            not from_embeddings
            and epoch == encoder_warmup_epochs
            and not encoder_unfrozen
        ):
            model.unfreeze_last_n_blocks(finetune_layers)
            encoder_unfrozen = True
            # Enable gradient checkpointing for unfrozen blocks — avoids storing
            # all intermediate activations (saves ~24 GB per block on 224×224 H&E)
            if model.encoder is not None:
                model.encoder.use_gradient_checkpointing = True
            # Rebuild optimizer to include newly unfrozen params
            optimizer = torch.optim.AdamW(
                filter(lambda p: p.requires_grad, model.parameters()),
                lr=lr_2a,
                weight_decay=weight_decay,
            )
        elif (
            external_backbone is not None
            and warmup_freeze_backbone
            and epoch == encoder_warmup_epochs
            and not encoder_unfrozen
        ):
            external_backbone._model.requires_grad_(True)
            external_backbone._model.train()
            encoder_unfrozen = True
            optimizer = torch.optim.AdamW(
                [
                    {"params": list(filter(lambda p: p.requires_grad, model.parameters())), "lr": lr_2a},
                    {"params": list(external_backbone._model.parameters()), "lr": lr_2a / 10},
                ],
                weight_decay=weight_decay,
            )
            logger.info(
                "CTransPath backbone unfrozen at epoch %d (backbone lr=%.2e)",
                epoch,
                lr_2a / 10,
            )

        loss_val = _train_one_epoch(
            model=model,
            patches=patches,
            c2s=cell_to_spot,
            oracle_props=oracle_props,
            train_spots=train_spots,
            spots_per_batch=spots_per_batch,
            augmentation=augmentation,
            optimizer=optimizer,
            device=device,
            lambda_consist=lambda_consist,
            lambda_var=lambda_var,
            lambda_cov=lambda_cov,
            lambda_sep=lambda_sep,
            sep_margin=sep_margin,
            lambda_sharp=0.0,  # no sharpening in 2A
            type_weights=type_weights_t,
            backbone=external_backbone,
        )
        train_losses.append(loss_val)

        val_kl = _validate(
            model, patches, cell_to_spot, oracle_props, val_spots, device,
            backbone=external_backbone,
        )
        val_losses.append(val_kl)

        if val_kl < _es_best_val:
            _es_best_val = val_kl
            _es_no_improve = 0
            if val_kl < _best_ckpt["val_kl"]:
                _best_ckpt["val_kl"] = val_kl
                _best_ckpt["model"] = copy.deepcopy(model.state_dict())
                if external_backbone is not None:
                    _best_ckpt["backbone"] = copy.deepcopy(external_backbone._model.state_dict())
        else:
            _es_no_improve += 1
            if early_stop_patience > 0 and _es_no_improve >= early_stop_patience:
                logger.info(
                    "Stage 2A early stop at epoch %d/%d (no val_kl improvement for %d epochs)",
                    epoch + 1, n_epochs_2a, early_stop_patience,
                )
                break

        if (epoch + 1) % 10 == 0:
            logger.info(
                "Stage 2A epoch %d/%d | train_loss=%.4f | val_kl=%.4f",
                epoch + 1, n_epochs_2a, loss_val, val_kl,
            )

    # ------------------------------------------------------------------
    # Optional prototype refresh between stages
    # ------------------------------------------------------------------
    if prototype_refresh:
        _refresh_prototypes(model, patches, cell_to_spot, device, backbone=external_backbone)

    # ------------------------------------------------------------------
    # Stage 2B optimizer (lower LR)
    # ------------------------------------------------------------------
    optimizer_2b = torch.optim.AdamW(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=lr_2b,
        weight_decay=weight_decay,
    )

    # ------------------------------------------------------------------
    # Stage 2B
    # ------------------------------------------------------------------
    _es_best_val = float("inf")
    _es_no_improve = 0
    for epoch in range(n_epochs_2b):
        frac = epoch / max(n_epochs_2b - 1, 1)
        # Linear temperature anneal
        model.temperature = tau_start + frac * (tau_end - tau_start)
        # Linear sharpening ramp
        cur_lambda_sharp = frac * lambda_sharp

        loss_val = _train_one_epoch(
            model=model,
            patches=patches,
            c2s=cell_to_spot,
            oracle_props=oracle_props,
            train_spots=train_spots,
            spots_per_batch=spots_per_batch,
            augmentation=augmentation,
            optimizer=optimizer_2b,
            device=device,
            lambda_consist=lambda_consist,
            lambda_var=lambda_var,
            lambda_cov=lambda_cov,
            lambda_sep=lambda_sep,
            sep_margin=sep_margin,
            lambda_sharp=cur_lambda_sharp,
            type_weights=type_weights_t,
            backbone=external_backbone,
        )
        train_losses.append(loss_val)

        val_kl = _validate(
            model, patches, cell_to_spot, oracle_props, val_spots, device,
            backbone=external_backbone,
        )
        val_losses.append(val_kl)

        if val_kl < _es_best_val:
            _es_best_val = val_kl
            _es_no_improve = 0
            if val_kl < _best_ckpt["val_kl"]:
                _best_ckpt["val_kl"] = val_kl
                _best_ckpt["model"] = copy.deepcopy(model.state_dict())
                if external_backbone is not None:
                    _best_ckpt["backbone"] = copy.deepcopy(external_backbone._model.state_dict())
        else:
            _es_no_improve += 1
            if early_stop_patience > 0 and _es_no_improve >= early_stop_patience:
                logger.info(
                    "Stage 2B early stop at epoch %d/%d (no val_kl improvement for %d epochs)",
                    epoch + 1, n_epochs_2b, early_stop_patience,
                )
                break

        if (epoch + 1) % 10 == 0:
            logger.info(
                "Stage 2B epoch %d/%d | train_loss=%.4f | val_kl=%.4f | tau=%.4f",
                epoch + 1, n_epochs_2b, loss_val, val_kl, model.temperature,
            )

    # ------------------------------------------------------------------
    # Restore best checkpoint before inference
    # ------------------------------------------------------------------
    if _best_ckpt["model"] is not None:
        model.load_state_dict(_best_ckpt["model"])
        if external_backbone is not None and _best_ckpt["backbone"] is not None:
            external_backbone._model.load_state_dict(_best_ckpt["backbone"])
        logger.info("Restored best model checkpoint (val_kl=%.4f)", _best_ckpt["val_kl"])

    # ------------------------------------------------------------------
    # Final inference
    # ------------------------------------------------------------------
    q_class, q_proto, z = _inference_all(model, patches, device, backbone=external_backbone)

    c2s_dev = cell_to_spot.to(device)
    spot_proportions = _scatter_mean(q_class, c2s_dev, dim=0, dim_size=num_spots)

    return {
        "model_state": model.state_dict(),
        "train_losses": train_losses,
        "val_losses": val_losses,
        "spot_proportions": spot_proportions.cpu(),
        "cell_assignments": q_class.cpu(),
        "cell_proto_assignments": q_proto.cpu(),
        "embeddings": z.cpu(),
        "temperature_final": model.temperature,
    }


# ---------------------------------------------------------------------------
# Task 5: Inference with TTA
# ---------------------------------------------------------------------------


@torch.no_grad()
def run_inference_tta(
    model: "PrototypeContrastiveModel",
    patches: torch.Tensor,
    cell_to_spot: torch.Tensor,
    num_spots: int,
    device: str = "cuda",
    batch_size: int = 256,
) -> Dict:
    """Run inference with 8-fold test-time augmentation.

    Averages output probabilities (q_class, q_proto) across 8 TTA views:
    4 rotations × 2 flips.

    Args:
        model: Trained PrototypeContrastiveModel.
        patches: (C, ch, 96, 96) cell patch tensor.
        cell_to_spot: (C,) integer spot indices.
        num_spots: Total number of spots S.
        device: Device string.
        batch_size: Mini-batch size for inference.

    Returns:
        Dict with keys:
            - ``q_class``: (C, T) averaged classifier probabilities.
            - ``q_proto``: (C, T) averaged prototype probabilities.
            - ``cell_types_hard``: (C,) argmax cell type assignments.
            - ``spot_proportions``: (S, T) scatter-mean spot proportions.
    """
    model.eval()
    model.to(device)

    N = patches.shape[0]
    T = model.num_types

    if model.from_embeddings:
        # TTA (rot90/flip) is meaningless for 1D embeddings — single forward pass
        q_class_list, q_proto_list = [], []
        for start in range(0, N, batch_size):
            batch = patches[start: start + batch_size].to(device)
            out = model(batch)
            q_class_list.append(out["q_class"])
            q_proto_list.append(out["q_proto"])
        q_class_avg = torch.cat(q_class_list, dim=0)
        q_proto_avg = torch.cat(q_proto_list, dim=0)
    else:
        # 8 TTA views as index transforms: (k_rot, flip)
        tta_params = [(k, flip) for k in range(4) for flip in [False, True]]  # 8 views

        # Accumulate averaged probabilities
        q_class_avg = torch.zeros(N, T, device=device)
        q_proto_avg = torch.zeros(N, T, device=device)

        for k_rot, do_flip in tta_params:
            q_class_list, q_proto_list = [], []
            for start in range(0, N, batch_size):
                batch = patches[start: start + batch_size].to(device)
                # Apply deterministic transforms
                batch = torch.rot90(batch, k=k_rot, dims=[-2, -1])
                if do_flip:
                    batch = torch.flip(batch, dims=[-1])
                out = model(batch)
                q_class_list.append(out["q_class"])
                q_proto_list.append(out["q_proto"])

            q_class_avg += torch.cat(q_class_list, dim=0)
            q_proto_avg += torch.cat(q_proto_list, dim=0)

        q_class_avg /= len(tta_params)
        q_proto_avg /= len(tta_params)

    # Re-normalize (averaging softmaxes may not sum to exactly 1 due to fp)
    q_class_avg = q_class_avg / q_class_avg.sum(dim=-1, keepdim=True).clamp(min=1e-8)
    q_proto_avg = q_proto_avg / q_proto_avg.sum(dim=-1, keepdim=True).clamp(min=1e-8)

    cell_types_hard = q_class_avg.argmax(dim=-1)

    c2s_dev = cell_to_spot.to(device)
    spot_proportions = _scatter_mean(q_class_avg, c2s_dev, dim=0, dim_size=num_spots)
    spot_proportions = spot_proportions / spot_proportions.sum(dim=-1, keepdim=True).clamp(min=1e-8)

    return {
        "q_class": q_class_avg.cpu(),
        "q_proto": q_proto_avg.cpu(),
        "cell_types_hard": cell_types_hard.cpu(),
        "spot_proportions": spot_proportions.cpu(),
    }
