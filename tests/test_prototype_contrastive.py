"""Tests for Prototype-Contrastive LLP loss functions and model."""

import sys
import os
import importlib
import pytest
import torch
import numpy as np

# Import directly from the module file to avoid the heavy CITEgeist/__init__.py
# which pulls in pandas/gurobi/etc and fails in lightweight test environments.
_model_dir = os.path.join(os.path.dirname(__file__), "..", "CITEgeist", "model", "morphology")
sys.path.insert(0, os.path.abspath(_model_dir))
try:
    _pc = importlib.import_module("prototype_contrastive")
except ImportError:
    pytest.skip("CITEgeist/model/morphology/prototype_contrastive.py not available", allow_module_level=True)
proportion_kl_loss = _pc.proportion_kl_loss
consistency_loss = _pc.consistency_loss
variance_covariance_loss = _pc.variance_covariance_loss
separation_loss = _pc.separation_loss
sharpening_loss = _pc.sharpening_loss
PrototypeContrastiveModel = _pc.PrototypeContrastiveModel


# ---------------------------------------------------------------------------
# proportion_kl_loss
# ---------------------------------------------------------------------------

def test_proportion_kl_loss_perfect_match():
    """KL(oracle || predicted) = 0 when oracle == predicted."""
    torch.manual_seed(0)
    S, T = 8, 5
    q = torch.softmax(torch.randn(S, T), dim=-1)
    loss = proportion_kl_loss(q, q)
    assert loss.item() < 1e-5, f"Expected ~0, got {loss.item()}"


def test_proportion_kl_loss_uniform_vs_peaked():
    """KL(uniform || peaked) should be substantially > 0.5."""
    S, T = 16, 6
    uniform = torch.full((S, T), 1.0 / T)
    # peaked: all mass on class 0
    peaked = torch.zeros(S, T)
    peaked[:, 0] = 1.0 - 1e-6
    peaked[:, 1:] = 1e-6 / (T - 1)
    peaked = peaked / peaked.sum(dim=-1, keepdim=True)
    loss = proportion_kl_loss(uniform, peaked)
    assert loss.item() > 0.5, f"Expected loss > 0.5, got {loss.item()}"


# ---------------------------------------------------------------------------
# consistency_loss
# ---------------------------------------------------------------------------

def test_consistency_loss_identical():
    """Symmetric KL = 0 when both distributions are identical."""
    torch.manual_seed(1)
    N, T = 32, 7
    q = torch.softmax(torch.randn(N, T), dim=-1)
    loss = consistency_loss(q, q)
    assert loss.item() < 1e-5, f"Expected ~0, got {loss.item()}"


def test_consistency_loss_disagreement():
    """Two very different distributions should produce loss > 1.0."""
    N, T = 32, 7
    # q_class: all mass on class 0
    q_class = torch.zeros(N, T)
    q_class[:, 0] = 1.0 - 1e-6
    q_class[:, 1:] = 1e-6 / (T - 1)
    q_class = q_class / q_class.sum(dim=-1, keepdim=True)
    # q_proto: all mass on last class
    q_proto = torch.zeros(N, T)
    q_proto[:, -1] = 1.0 - 1e-6
    q_proto[:, :-1] = 1e-6 / (T - 1)
    q_proto = q_proto / q_proto.sum(dim=-1, keepdim=True)
    loss = consistency_loss(q_class, q_proto)
    assert loss.item() > 1.0, f"Expected loss > 1.0, got {loss.item()}"


# ---------------------------------------------------------------------------
# variance_covariance_loss
# ---------------------------------------------------------------------------

def test_variance_covariance_loss_collapsed():
    """All-identical embeddings → variance loss > 0.9."""
    N, D = 64, 32
    z = torch.ones(N, D)  # zero std → max variance penalty
    var_loss, cov_loss = variance_covariance_loss(z, gamma=1.0)
    assert var_loss.item() > 0.9, f"Expected var_loss > 0.9, got {var_loss.item()}"


def test_variance_covariance_loss_decorrelated():
    """Standard normal embeddings → low variance and covariance loss."""
    torch.manual_seed(42)
    N, D = 512, 64
    z = torch.randn(N, D)
    var_loss, cov_loss = variance_covariance_loss(z, gamma=1.0)
    assert var_loss.item() < 0.1, f"Expected var_loss < 0.1, got {var_loss.item()}"
    assert cov_loss.item() < 0.05, f"Expected cov_loss < 0.05, got {cov_loss.item()}"


# ---------------------------------------------------------------------------
# separation_loss
# ---------------------------------------------------------------------------

def test_separation_loss_identical_protos():
    """Identical prototype vectors → cosine sim = 1 >> margin → loss > 0.5."""
    T, D = 5, 16
    proto = torch.randn(1, D)
    proto = proto / proto.norm(dim=-1, keepdim=True)
    prototypes = proto.expand(T, D)  # all identical
    loss = separation_loss(prototypes, margin=0.1)
    assert loss.item() > 0.5, f"Expected loss > 0.5, got {loss.item()}"


def test_separation_loss_orthogonal_protos():
    """Orthogonal prototypes have cosine sim = 0 < margin=0.1 → loss ≈ 0."""
    # Build T orthogonal unit vectors (T <= D)
    T, D = 4, 16
    Q, _ = torch.linalg.qr(torch.randn(D, T))
    prototypes = Q.T  # (T, D), rows are orthonormal
    loss = separation_loss(prototypes, margin=0.1)
    assert loss.item() < 1e-5, f"Expected ~0, got {loss.item()}"


# ---------------------------------------------------------------------------
# sharpening_loss
# ---------------------------------------------------------------------------

def test_sharpening_loss():
    """Uniform distribution has higher entropy (sharpening loss) than peaked."""
    N, T = 32, 6
    uniform_q = torch.full((N, T), 1.0 / T)
    # peaked: near one-hot
    peaked_q = torch.zeros(N, T)
    peaked_q[:, 0] = 1.0 - 1e-6
    peaked_q[:, 1:] = 1e-6 / (T - 1)
    peaked_q = peaked_q / peaked_q.sum(dim=-1, keepdim=True)

    loss_uniform = sharpening_loss(uniform_q)
    loss_peaked = sharpening_loss(peaked_q)
    assert loss_uniform.item() > loss_peaked.item(), (
        f"Uniform loss ({loss_uniform.item():.4f}) should exceed "
        f"peaked loss ({loss_peaked.item():.6f})"
    )


# ---------------------------------------------------------------------------
# PrototypeContrastiveModel
# ---------------------------------------------------------------------------


def test_model_forward_shapes():
    """Forward pass produces correctly shaped outputs; distributions sum to 1; z_norm unit."""
    torch.manual_seed(0)
    B, T, D = 8, 6, 128
    model = PrototypeContrastiveModel(
        num_types=T,
        embed_dim=D,
        in_channels=3,
        temperature=0.1,
        encoder_depth=4,
        encoder_embed_dim=384,
        encoder_num_heads=6,
    )
    model.eval()
    patches = torch.randn(B, 3, 96, 96)
    with torch.no_grad():
        out = model(patches)

    # Shape checks
    assert out["q_class"].shape == (B, T), f"q_class shape {out['q_class'].shape}"
    assert out["q_proto"].shape == (B, T), f"q_proto shape {out['q_proto'].shape}"
    assert out["z"].shape == (B, D), f"z shape {out['z'].shape}"
    assert out["z_norm"].shape == (B, D), f"z_norm shape {out['z_norm'].shape}"

    # Distributions sum to 1
    assert torch.allclose(out["q_class"].sum(dim=-1), torch.ones(B), atol=1e-5), \
        "q_class rows do not sum to 1"
    assert torch.allclose(out["q_proto"].sum(dim=-1), torch.ones(B), atol=1e-5), \
        "q_proto rows do not sum to 1"

    # z_norm should be unit vectors
    norms = out["z_norm"].norm(dim=-1)
    assert torch.allclose(norms, torch.ones(B), atol=1e-5), \
        f"z_norm norms not all 1: {norms}"


def test_model_prototype_normalization():
    """get_normalized_prototypes() returns unit vectors after a forward pass."""
    torch.manual_seed(1)
    T, D = 6, 128
    model = PrototypeContrastiveModel(
        num_types=T,
        embed_dim=D,
        in_channels=3,
        encoder_depth=4,
    )
    model.eval()
    # Run a forward pass to ensure model is in a consistent state
    with torch.no_grad():
        _ = model(torch.randn(4, 3, 96, 96))
        proto_norm = model.get_normalized_prototypes()

    assert proto_norm.shape == (T, D), f"Expected ({T}, {D}), got {proto_norm.shape}"
    norms = proto_norm.norm(dim=-1)
    assert torch.allclose(norms, torch.ones(T), atol=1e-5), \
        f"Prototype norms not all 1: {norms}"


def test_model_freeze_unfreeze():
    """freeze_encoder() + unfreeze_last_n_blocks(2) freezes first blocks, unfreezes last 2."""
    torch.manual_seed(2)
    T, depth = 6, 4
    model = PrototypeContrastiveModel(
        num_types=T,
        encoder_depth=depth,
    )

    model.freeze_encoder()
    model.unfreeze_last_n_blocks(2)

    # blocks[0] and blocks[1] should be frozen
    for i in range(depth - 2):
        for p in model.encoder.blocks[i].parameters():
            assert not p.requires_grad, \
                f"encoder.blocks[{i}] param should be frozen after freeze+unfreeze(2)"

    # blocks[2] and blocks[3] (last 2) should be unfrozen
    for i in range(depth - 2, depth):
        for p in model.encoder.blocks[i].parameters():
            assert p.requires_grad, \
                f"encoder.blocks[{i}] param should be trainable after unfreeze_last_n_blocks(2)"

    # encoder.norm should be unfrozen (included in unfreeze_last_n_blocks)
    for p in model.encoder.norm.parameters():
        assert p.requires_grad, "encoder.norm should be trainable after unfreeze_last_n_blocks"

    # projector and heads should always be trainable
    for name, p in model.projector.named_parameters():
        assert p.requires_grad, f"projector param '{name}' should always be trainable"
    for name, p in model.classifier_head.named_parameters():
        assert p.requires_grad, f"classifier_head param '{name}' should always be trainable"


# ---------------------------------------------------------------------------
# Task 3: CellPatchAugmentation and tta_8x
# ---------------------------------------------------------------------------

CellPatchAugmentation = _pc.CellPatchAugmentation
tta_8x = _pc.tta_8x


def test_continuous_augmentation_shape():
    """Single (C, H, W) patch augmentation preserves shape."""
    torch.manual_seed(0)
    aug = CellPatchAugmentation(continuous_rotation=True, intensity_scale=0.2, noise_std=0.05)
    x = torch.randn(3, 96, 96)
    x_aug = aug(x)
    assert x_aug.shape == x.shape, f"Expected {x.shape}, got {x_aug.shape}"


def test_continuous_augmentation_batch():
    """Batch (B, C, H, W) augmentation preserves shape."""
    torch.manual_seed(1)
    aug = CellPatchAugmentation(continuous_rotation=True, intensity_scale=0.2, noise_std=0.05)
    x = torch.randn(8, 3, 96, 96)
    x_aug = aug(x)
    assert x_aug.shape == x.shape, f"Expected {x.shape}, got {x_aug.shape}"


def test_tta_8x():
    """tta_8x produces exactly 8 views, all with correct shape."""
    torch.manual_seed(2)
    patch = torch.randn(3, 96, 96)
    views = tta_8x(patch)
    assert len(views) == 8, f"Expected 8 views, got {len(views)}"
    for i, v in enumerate(views):
        assert v.shape == (3, 96, 96), f"View {i} has wrong shape: {v.shape}"


# ---------------------------------------------------------------------------
# Task 4: Training loop smoke tests
# ---------------------------------------------------------------------------

train_prototype_contrastive = _pc.train_prototype_contrastive


def _make_small_dataset(n_cells=50, n_spots=5, n_types=3, ch=3, seed=0):
    """Create minimal synthetic training data."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    patches = torch.randn(n_cells, ch, 96, 96)
    c2s = torch.randint(0, n_spots, (n_cells,))
    oracle = torch.softmax(torch.randn(n_spots, n_types), dim=-1)
    return patches, c2s, oracle


def test_train_one_epoch_smoke():
    """train_prototype_contrastive returns model_state and correct loss list length."""
    patches, c2s, oracle = _make_small_dataset(n_cells=50, n_spots=5, n_types=3)
    result = train_prototype_contrastive(
        patches=patches,
        cell_to_spot=c2s,
        oracle_props=oracle,
        num_types=3,
        embed_dim=32,
        in_channels=3,
        encoder_depth=2,
        n_epochs_2a=2,
        n_epochs_2b=1,
        encoder_warmup_epochs=10,   # won't trigger in 2 epochs
        finetune_layers=1,
        spots_per_batch=5,
        device="cpu",
        prototype_refresh=False,
        seed=0,
    )
    assert "model_state" in result, "Result must contain model_state"
    # total epochs = n_epochs_2a + n_epochs_2b = 3
    assert len(result["train_losses"]) == 3, (
        f"Expected 3 train loss entries, got {len(result['train_losses'])}"
    )
    assert len(result["val_losses"]) == 3, (
        f"Expected 3 val loss entries, got {len(result['val_losses'])}"
    )


def test_train_outputs_proportions():
    """spot_proportions has correct shape and rows sum to 1."""
    patches, c2s, oracle = _make_small_dataset(n_cells=50, n_spots=5, n_types=3)
    result = train_prototype_contrastive(
        patches=patches,
        cell_to_spot=c2s,
        oracle_props=oracle,
        num_types=3,
        embed_dim=32,
        in_channels=3,
        encoder_depth=2,
        n_epochs_2a=2,
        n_epochs_2b=1,
        encoder_warmup_epochs=10,
        finetune_layers=1,
        spots_per_batch=5,
        device="cpu",
        prototype_refresh=False,
        seed=1,
    )
    sp = result["spot_proportions"]
    assert sp.shape == (5, 3), f"Expected (5, 3), got {sp.shape}"
    row_sums = sp.sum(dim=-1)
    assert torch.allclose(row_sums, torch.ones(5), atol=1e-5), (
        f"spot_proportions rows do not sum to 1: {row_sums}"
    )


# ---------------------------------------------------------------------------
# Task 5: run_inference_tta
# ---------------------------------------------------------------------------

run_inference_tta = _pc.run_inference_tta


def test_inference_tta_shapes():
    """run_inference_tta returns correct shapes and valid probability ranges."""
    torch.manual_seed(3)
    n_cells, n_spots, n_types = 20, 5, 3
    model = PrototypeContrastiveModel(
        num_types=n_types,
        embed_dim=32,
        in_channels=3,
        encoder_depth=2,
    )
    patches = torch.randn(n_cells, 3, 96, 96)
    c2s = torch.randint(0, n_spots, (n_cells,))

    result = run_inference_tta(
        model=model,
        patches=patches,
        cell_to_spot=c2s,
        num_spots=n_spots,
        device="cpu",
        batch_size=8,
    )

    assert result["q_class"].shape == (n_cells, n_types), (
        f"q_class shape wrong: {result['q_class'].shape}"
    )
    assert result["q_proto"].shape == (n_cells, n_types), (
        f"q_proto shape wrong: {result['q_proto'].shape}"
    )
    assert result["cell_types_hard"].shape == (n_cells,), (
        f"cell_types_hard shape wrong: {result['cell_types_hard'].shape}"
    )
    assert result["spot_proportions"].shape == (n_spots, n_types), (
        f"spot_proportions shape wrong: {result['spot_proportions'].shape}"
    )
    # Values should be in [0, 1]
    assert result["q_class"].min() >= 0.0
    assert result["q_class"].max() <= 1.0 + 1e-5


def test_inference_tta_improves_over_single():
    """q_class rows sum to 1 (valid distributions after TTA averaging)."""
    torch.manual_seed(4)
    n_cells, n_spots, n_types = 20, 5, 3
    model = PrototypeContrastiveModel(
        num_types=n_types,
        embed_dim=32,
        in_channels=3,
        encoder_depth=2,
    )
    patches = torch.randn(n_cells, 3, 96, 96)
    c2s = torch.randint(0, n_spots, (n_cells,))

    result = run_inference_tta(
        model=model,
        patches=patches,
        cell_to_spot=c2s,
        num_spots=n_spots,
        device="cpu",
        batch_size=8,
    )

    row_sums = result["q_class"].sum(dim=-1)
    assert torch.allclose(row_sums, torch.ones(n_cells), atol=1e-4), (
        f"q_class rows do not sum to 1 after TTA: {row_sums}"
    )
