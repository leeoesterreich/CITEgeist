"""Tests for argmax global inference mode in PC-MIL."""
import numpy as np
import torch
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from unittest.mock import MagicMock
from model.pc_mil_inference import pc_mil_infer_spot


def _make_mock_model(n_nuclei, n_types, device="cpu"):
    model = MagicMock()
    logits = torch.zeros(n_nuclei, n_types)
    for i in range(n_nuclei):
        logits[i, i % n_types] = 5.0
    attention = torch.softmax(logits, dim=1)
    proportions = attention.mean(dim=0)
    reconstructed = torch.zeros(5)
    model.forward_with_logits.return_value = (logits, attention, proportions, reconstructed)
    model.eval = MagicMock(return_value=model)
    model.parameters.return_value = iter([torch.zeros(1)])
    return model


def test_argmax_global_no_constraints():
    n_nuclei, n_types = 6, 3
    model = _make_mock_model(n_nuclei, n_types)
    cell_type_names = ["TypeA", "TypeB", "TypeC"]
    detected = np.ones(n_types, dtype=bool)

    result = pc_mil_infer_spot(
        model=model,
        image_features=torch.randn(n_nuclei, 384),
        protein_proportions=torch.ones(n_types) / n_types,
        detected_types=detected,
        cell_type_names=cell_type_names,
        inference_mode="argmax_global",
    )

    assert len(result) == n_nuclei
    assert "cell_type" in result.columns
    assert "confidence" in result.columns
    assert "protein_score" in result.columns
    expected_types = [cell_type_names[i % n_types] for i in range(n_nuclei)]
    assert list(result["cell_type"]) == expected_types


def test_hungarian_is_default():
    """Default inference_mode should be hungarian_constrained."""
    n_nuclei, n_types = 4, 2
    model = _make_mock_model(n_nuclei, n_types)
    cell_type_names = ["TypeA", "TypeB"]
    detected = np.ones(n_types, dtype=bool)

    result = pc_mil_infer_spot(
        model=model,
        image_features=torch.randn(n_nuclei, 384),
        protein_proportions=torch.tensor([0.75, 0.25]),
        detected_types=detected,
        cell_type_names=cell_type_names,
    )

    assert len(result) == n_nuclei
    type_counts = result["cell_type"].value_counts()
    assert type_counts.get("TypeA", 0) == 3
    assert type_counts.get("TypeB", 0) == 1
