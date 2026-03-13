"""Tests for Protein-Conditioned MIL model."""
import torch
import torch.nn.functional as F
import numpy as np
import pytest


def test_pcmil_forward_shape():
    """PCMILModel forward pass produces correct output shapes."""
    from model.pc_mil import PCMILModel

    K = 7   # cell types
    M = 15  # protein markers
    N = 12  # nuclei in this spot

    model = PCMILModel(
        image_dim=384,
        n_types=K,
        n_markers=M,
        image_proj_dim=64,
        protein_context_dim=32,
        hidden_dim=128,
    )

    image_features = torch.randn(N, 384)
    protein_props = torch.rand(K)
    protein_props = protein_props / protein_props.sum()

    proportions, attention, reconstructed = model(image_features, protein_props)

    assert proportions.shape == (K,), f"Expected ({K},), got {proportions.shape}"
    assert attention.shape == (N, K), f"Expected ({N},{K}), got {attention.shape}"
    assert reconstructed.shape == (M,), f"Expected ({M},), got {reconstructed.shape}"

    # Proportions should sum to ~1
    assert abs(proportions.sum().item() - 1.0) < 1e-5

    # Each row of attention should sum to ~1 (softmax over types)
    row_sums = attention.sum(dim=1)
    assert torch.allclose(row_sums, torch.ones(N), atol=1e-5)


def test_pcmil_gradient_flow():
    """Gradients flow through all learnable components."""
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=7, n_markers=15)
    image_features = torch.randn(10, 384)
    protein_props = torch.rand(7)
    protein_props = protein_props / protein_props.sum()

    proportions, attention, reconstructed = model(image_features, protein_props)

    # Use combined loss (proportion + reconstruction) to reach all params
    target = torch.rand(7)
    target = target / target.sum()
    target_recon = torch.randn(15)
    loss = (torch.nn.functional.mse_loss(proportions, target)
            + torch.nn.functional.mse_loss(reconstructed, target_recon))
    loss.backward()

    # Check key parameter groups have gradients
    for name, param in model.named_parameters():
        if param.requires_grad:
            assert param.grad is not None, f"No gradient for {name}"
            assert not torch.all(param.grad == 0), f"Zero gradient for {name}"


def test_pcmil_profile_initialization():
    """Protein profile matrix can be initialized from cell_profile_dict."""
    from model.pc_mil import PCMILModel, build_profile_matrix

    cell_profile_dict = {
        "T_cells": ["CD3", "CD4"],
        "B_cells": ["CD19", "CD20"],
        "Macrophages": ["CD68"],
    }
    marker_names = ["CD3", "CD4", "CD19", "CD20", "CD68"]

    profile_matrix = build_profile_matrix(cell_profile_dict, marker_names)

    assert profile_matrix.shape == (3, 5)
    # T_cells should have 1.0 for CD3 (idx 0) and CD4 (idx 1)
    assert profile_matrix[0, 0] == 1.0
    assert profile_matrix[0, 1] == 1.0
    assert profile_matrix[0, 2] == 0.0
    # B_cells should have 1.0 for CD19 (idx 2) and CD20 (idx 3)
    assert profile_matrix[1, 2] == 1.0
    assert profile_matrix[1, 3] == 1.0

    # Use it to initialize model
    model = PCMILModel(
        image_dim=384,
        n_types=3,
        n_markers=5,
        init_profile_matrix=torch.tensor(profile_matrix, dtype=torch.float32),
    )
    assert torch.allclose(
        model.protein_profiles.data,
        torch.tensor(profile_matrix, dtype=torch.float32),
    )


def test_pcmil_variable_nuclei():
    """Model handles different numbers of nuclei per spot."""
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=7, n_markers=15)
    protein_props = torch.rand(7)
    protein_props = protein_props / protein_props.sum()

    for n in [1, 5, 20, 50]:
        feats = torch.randn(n, 384)
        props, attn, recon = model(feats, protein_props)
        assert props.shape == (7,)
        assert attn.shape == (n, 7)


def test_pcmil_softmax_dim():
    """CRITICAL: Verify softmax is over types (dim=1), not nuclei (dim=0)."""
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    feats = torch.randn(10, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])

    _, attention, _ = model(feats, protein_props)

    # Each nucleus row sums to 1 (softmax over types)
    row_sums = attention.sum(dim=1)
    assert torch.allclose(row_sums, torch.ones(10), atol=1e-5), \
        "Softmax must be over types (dim=1)! Row sums should be 1.0"

    # Column sums should NOT all be equal (would indicate dim=0 bug)
    col_sums = attention.sum(dim=0)
    assert not torch.allclose(col_sums, col_sums[0].expand(3), atol=0.1), \
        "Column sums are suspiciously uniform — possible dim=0 softmax bug"


# === Loss function tests ===

def test_pcmil_loss_components():
    """All 4 loss components compute without errors."""
    from model.pc_mil_training import compute_pc_mil_loss

    K, M, N = 7, 15, 10
    pred_props = torch.rand(K)
    pred_props = pred_props / pred_props.sum()
    true_props = torch.rand(K)
    true_props = true_props / true_props.sum()
    attention = F.softmax(torch.randn(N, K), dim=1)
    reconstructed = torch.randn(M)
    observed_protein = torch.randn(M)

    loss, components = compute_pc_mil_loss(
        pred_proportions=pred_props,
        true_proportions=true_props,
        attention=attention,
        reconstructed_protein=reconstructed,
        observed_protein=observed_protein,
    )

    assert loss.ndim == 0, "Loss should be scalar"
    assert not torch.isnan(loss), "Loss should not be NaN"
    assert "proportion" in components
    assert "reconstruction" in components
    assert "entropy" in components
    assert "diversity" in components


def test_diversity_loss_penalizes_absent_types():
    """Diversity loss increases when active types have near-zero prediction."""
    from model.pc_mil_training import compute_pc_mil_loss

    K, M, N = 3, 5, 10
    true_props = torch.tensor([0.4, 0.4, 0.2])  # all types present
    observed = torch.randn(M)

    # Good prediction: all types represented
    attn_good = F.softmax(torch.randn(N, K), dim=1)
    props_good = attn_good.mean(dim=0)
    recon_good = torch.randn(M)
    _, comp_good = compute_pc_mil_loss(props_good, true_props, attn_good, recon_good, observed)

    # Bad prediction: one type collapsed to 0
    attn_bad = torch.zeros(N, K)
    attn_bad[:, 0] = 0.7
    attn_bad[:, 1] = 0.3  # type 2 is absent
    props_bad = attn_bad.mean(dim=0)
    recon_bad = torch.randn(M)
    _, comp_bad = compute_pc_mil_loss(props_bad, true_props, attn_bad, recon_bad, observed)

    assert comp_bad["diversity"] > comp_good["diversity"], \
        "Diversity loss should be higher when an active type is predicted as absent"


def test_entropy_loss_penalizes_uniform():
    """Entropy loss increases when assignments are uniform (high entropy)."""
    from model.pc_mil_training import compute_pc_mil_loss

    K, M, N = 5, 10, 10
    true_props = torch.ones(K) / K
    observed = torch.randn(M)

    # Confident: each nucleus assigned to one type
    attn_confident = torch.zeros(N, K)
    for i in range(N):
        attn_confident[i, i % K] = 1.0
    # Add small epsilon for log stability
    attn_confident = attn_confident.clamp(min=1e-6)
    attn_confident = attn_confident / attn_confident.sum(dim=1, keepdim=True)

    # Uniform: all nuclei spread evenly
    attn_uniform = torch.ones(N, K) / K

    _, comp_confident = compute_pc_mil_loss(
        attn_confident.mean(0), true_props, attn_confident, torch.randn(M), observed
    )
    _, comp_uniform = compute_pc_mil_loss(
        attn_uniform.mean(0), true_props, attn_uniform, torch.randn(M), observed
    )

    assert comp_uniform["entropy"] > comp_confident["entropy"], \
        "Entropy loss should be higher for uniform assignments"


# === Inference tests ===

def test_pcmil_inference_basic():
    """Inference produces per-nucleus assignments with correct format."""
    from model.pc_mil_inference import pc_mil_infer_spot
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    model.eval()

    N = 8
    image_features = torch.randn(N, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])
    detected = np.array([True, True, True])

    result = pc_mil_infer_spot(
        model=model,
        image_features=image_features,
        protein_proportions=protein_props,
        detected_types=detected,
        cell_type_names=["TypeA", "TypeB", "TypeC"],
    )

    assert len(result) == N
    assert "cell_type" in result.columns
    assert "confidence" in result.columns
    assert set(result["cell_type"].unique()).issubset({"TypeA", "TypeB", "TypeC"})


def test_pcmil_inference_detection_mask():
    """Undetected types get zero assignments."""
    from model.pc_mil_inference import pc_mil_infer_spot
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    model.eval()

    N = 10
    image_features = torch.randn(N, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])
    # Only type 0 detected
    detected = np.array([True, False, False])

    result = pc_mil_infer_spot(
        model=model,
        image_features=image_features,
        protein_proportions=protein_props,
        detected_types=detected,
        cell_type_names=["TypeA", "TypeB", "TypeC"],
    )

    # All assignments should be TypeA (only detected type)
    assert (result["cell_type"] == "TypeA").all()


def test_pcmil_inference_all_false_detection():
    """All-false detection mask falls back to no masking (no NaN)."""
    from model.pc_mil_inference import pc_mil_infer_spot
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    model.eval()

    N = 5
    image_features = torch.randn(N, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])
    detected = np.array([False, False, False])  # nothing detected

    result = pc_mil_infer_spot(
        model=model,
        image_features=image_features,
        protein_proportions=protein_props,
        detected_types=detected,
        cell_type_names=["TypeA", "TypeB", "TypeC"],
    )

    assert len(result) == N
    assert not result["confidence"].isna().any(), "Should not have NaN confidences"
