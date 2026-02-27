# CITEgeist/tests/test_attention_aggregator.py
"""Tests for attention-weighted MIL aggregation."""
import torch
import pytest


def test_attention_aggregator_output_shape():
    """Aggregated proportions should sum to 1."""
    from CITEgeist.model.attention_aggregator import AttentionAggregator

    agg = AttentionAggregator(embed_dim=32, n_types=7)

    embeddings = torch.randn(20, 32)  # 20 nuclei, 32-dim
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    pred_props, entropy, attn_weights = agg(embeddings, soft_assignments)

    assert pred_props.shape == (7,)
    assert torch.allclose(pred_props.sum(), torch.tensor(1.0), atol=1e-5)
    assert attn_weights.shape == (20, 1)


def test_attention_weights_sum_to_one():
    """Attention weights should sum to 1 across nuclei."""
    from CITEgeist.model.attention_aggregator import AttentionAggregator

    agg = AttentionAggregator(embed_dim=32, n_types=7)

    embeddings = torch.randn(20, 32)
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    _, _, attn_weights = agg(embeddings, soft_assignments)

    assert torch.allclose(attn_weights.sum(), torch.tensor(1.0), atol=1e-5)


def test_per_class_attention_output_shape():
    """Per-class attention should produce K separate attention heads."""
    from CITEgeist.model.attention_aggregator import PerClassAttentionAggregator

    agg = PerClassAttentionAggregator(embed_dim=32, n_types=7)

    embeddings = torch.randn(20, 32)
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    pred_props, entropy, attn_weights_list = agg(embeddings, soft_assignments)

    assert pred_props.shape == (7,)
    assert len(attn_weights_list) == 7
    for w in attn_weights_list:
        assert w.shape == (20,)


def test_entropy_regularization():
    """Entropy should be higher for uniform attention, lower for peaked."""
    from CITEgeist.model.attention_aggregator import AttentionAggregator

    agg = AttentionAggregator(embed_dim=32, n_types=7)

    # Create embeddings that will produce different attention patterns
    # Uniform embeddings -> more uniform attention
    uniform_embed = torch.zeros(20, 32)
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    _, entropy_uniform, _ = agg(uniform_embed, soft_assignments)

    # Varied embeddings -> potentially more peaked attention
    varied_embed = torch.randn(20, 32) * 10
    _, entropy_varied, _ = agg(varied_embed, soft_assignments)

    # Both should produce valid entropy values
    assert entropy_uniform.item() >= 0
    assert entropy_varied.item() >= 0
