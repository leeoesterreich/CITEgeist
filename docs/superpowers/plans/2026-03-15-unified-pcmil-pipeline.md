# Unified PC-MIL Pipeline Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unify the PC-MIL single-cell assignment pipeline across Xenium (DAPI+boundary) and patient Visium (H&E) modalities, re-run Module 3 with 9 cell types, and validate assignments via marker gene expression.

**Architecture:** Standalone orchestrator scripts (not wired into core CitegeistModel API) that compose existing modules: Module 3 runner → modality-aware Cellpose → ViT-S feature extraction → unified PC-MIL → marker gene validation. Each step writes intermediate outputs and marker files for SLURM gating.

**Tech Stack:** Python 3.10, PyTorch, timm (ViT-S), Cellpose v3, Gurobi 12, SLURM array jobs

**Spec:** `docs/superpowers/specs/2026-03-14-unified-pcmil-pipeline-design.md`

**Step Numbering:** The spec defines 5 conceptual steps. This plan consolidates Steps 2+3 into one runner:

| Plan Script | Spec Steps | SLURM Script |
|------------|------------|-------------|
| `run_unified_step1_module3.py` | Step 1 (Module 3) | `sbatch_unified_step1.sh` |
| `run_unified_step2_features.py` | Steps 2+3 (Cellpose + ViT) | `sbatch_unified_step2.sh` |
| `run_unified_step3_pcmil.py` | Step 4 (PC-MIL) | `sbatch_unified_step3.sh` |
| `run_unified_step4_validate.py` | Step 5 (Validation) | `sbatch_unified_step4.sh` |

**Code Gap Resolution:**

| Gap | Severity | Resolution |
|-----|----------|-----------|
| #1 Profile dict schema | Critical | Task 1: `flatten_profile_dict()` adapter |
| #2 Inference hardcoded | Critical | Task 2: `inference_mode` parameter |
| #3 Cellpose API | Critical | **Bypassed**: Step 2 runner calls Cellpose directly (not via core API) |
| #4 Trainer val_dataset | Warning | Task 3: make optional |
| #5 Loss hyperparameters | Warning | Task 5: pinned in `unified_config.py` |
| #6 Feature output format | Warning | **Replaced**: Step 2 runner writes unified arrays directly |
| #7 Core API not wired | Warning | **Bypassed**: standalone orchestrator scripts |

**Xenium support:** SLURM scripts cover 12 patient samples only. Xenium is a follow-up task — runners support `--modality dapi` but Xenium-specific SLURM scripts and spatial coord loading are deferred.

---

## File Structure

### Files to Modify

| File | Responsibility | Changes |
|------|---------------|---------|
| `CITEgeist/model/pc_mil.py` | PC-MIL model + profile matrix | Add `flatten_profile_dict()` adapter |
| `CITEgeist/model/pc_mil_inference.py` | Per-spot inference | Add `inference_mode` param (`argmax_global` vs `hungarian_constrained`) |
| `CITEgeist/model/pc_mil_training.py` | Training loop | Make `val_dataset` optional; add loss-plateau early stopping |
| `CITEgeist/model/vit_extractor.py` | ViT feature extraction | Change default to `vit_small_patch16_224` |

### Files to Create

| File | Responsibility |
|------|---------------|
| `CITEgeist/model/unified_config.py` | Shared constants: 9-type profile dict, sample lists, paths, marker validation dict |
| `CITEgeist/model/marker_validation.py` | Single-cell marker gene scoring (Step 5) |
| `CITEgeist/examples/run_unified_step1_module3.py` | Step 1 runner: Module 3 with 9 types + Cellpose nuclei |
| `CITEgeist/examples/run_unified_step2_features.py` | Step 2+3 runner: Cellpose patches + ViT extraction |
| `CITEgeist/examples/run_unified_step3_pcmil.py` | Step 4 runner: PC-MIL train + infer per sample |
| `CITEgeist/examples/run_unified_step4_validate.py` | Step 5 runner: Marker gene validation |
| `CITEgeist/examples/sbatch_unified_step1.sh` | SLURM: Module 3 array job |
| `CITEgeist/examples/sbatch_unified_step2.sh` | SLURM: Feature extraction array job (GPU) |
| `CITEgeist/examples/sbatch_unified_step3.sh` | SLURM: PC-MIL training array job (GPU) |
| `CITEgeist/examples/sbatch_unified_step4.sh` | SLURM: Validation array job (CPU) |
| `CITEgeist/tests/test_profile_adapter.py` | Test profile dict flattening |
| `CITEgeist/tests/test_argmax_inference.py` | Test argmax global inference mode |
| `CITEgeist/tests/test_marker_validation.py` | Test marker gene scoring logic |

---

## Chunk 1: Core Code Fixes (Close Code Gaps)

### Task 1: Profile Dict Adapter

Closes **Code Gap #1**: `build_profile_matrix()` expects flat `{type: [markers]}` but Module 3 uses nested `{"Major": [...]}`.

**Files:**
- Modify: `CITEgeist/model/pc_mil.py:1-30`
- Test: `CITEgeist/tests/test_profile_adapter.py`

- [ ] **Step 1: Write the failing test**

```python
# CITEgeist/tests/test_profile_adapter.py
"""Tests for profile dict adapter between Module 3 and PC-MIL formats."""
import numpy as np
import pytest
from model.pc_mil import flatten_profile_dict, build_profile_matrix


def test_flatten_nested_profile_dict():
    """Nested {'Major': [...]} format should flatten to {type: [markers]}."""
    nested = {
        "Cancer": {"Major": ["EPCAM-1"]},
        "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    }
    flat = flatten_profile_dict(nested)
    assert flat == {
        "Cancer": ["EPCAM-1"],
        "Macrophages": ["CD68-1", "CD163-1"],
    }


def test_flatten_already_flat_profile_dict():
    """Already flat {type: [markers]} should pass through unchanged."""
    flat_input = {
        "Cancer": ["EPCAM-1"],
        "Macrophages": ["CD68-1", "CD163-1"],
    }
    flat = flatten_profile_dict(flat_input)
    assert flat == flat_input


def test_flatten_with_major_minor():
    """Nested with Major + Minor should merge both lists."""
    nested = {
        "Cancer": {"Major": ["EPCAM-1"], "Minor": ["SDC1-1"]},
    }
    flat = flatten_profile_dict(nested)
    assert set(flat["Cancer"]) == {"EPCAM-1", "SDC1-1"}


def test_build_profile_matrix_with_nested_input():
    """build_profile_matrix should accept nested dicts via auto-flattening."""
    nested = {
        "Cancer": {"Major": ["EPCAM-1"]},
        "Macrophages": {"Major": ["CD68-1"]},
    }
    all_markers = ["EPCAM-1", "CD68-1", "CD163-1"]
    matrix = build_profile_matrix(nested, all_markers)
    assert matrix.shape == (2, 3)
    # Cancer row: EPCAM-1=1, CD68-1=0, CD163-1=0
    np.testing.assert_array_equal(matrix[0], [1, 0, 0])
    # Macrophages row: EPCAM-1=0, CD68-1=1, CD163-1=0
    np.testing.assert_array_equal(matrix[1], [0, 1, 0])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_profile_adapter.py -v`
Expected: FAIL with `ImportError: cannot import name 'flatten_profile_dict'`

- [ ] **Step 3: Implement `flatten_profile_dict` and update `build_profile_matrix`**

In `CITEgeist/model/pc_mil.py`, add at the top (after imports):

```python
def flatten_profile_dict(profile_dict: dict) -> dict:
    """Flatten nested Module 3 profile dict to flat {type: [markers]} for PC-MIL.

    Accepts both formats:
      - Nested: {"Cancer": {"Major": ["EPCAM-1"], "Minor": ["SDC1-1"]}}
      - Flat:   {"Cancer": ["EPCAM-1"]}

    Returns flat format, merging Major + Minor if both present.
    """
    flat = {}
    for cell_type, value in profile_dict.items():
        if isinstance(value, dict):
            markers = []
            for category in ("Major", "Minor"):
                markers.extend(value.get(category, []))
            flat[cell_type] = markers
        elif isinstance(value, list):
            flat[cell_type] = value
        else:
            raise ValueError(f"Unexpected profile format for {cell_type}: {type(value)}")
    return flat
```

Update existing `build_profile_matrix()` to auto-flatten:

```python
def build_profile_matrix(cell_profile_dict, all_markers):
    """Build (K, M) binary profile matrix. Auto-flattens nested dicts."""
    cell_profile_dict = flatten_profile_dict(cell_profile_dict)
    # ... rest of existing implementation unchanged
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd CITEgeist && python -m pytest tests/test_profile_adapter.py -v`
Expected: All 4 tests PASS

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/pc_mil.py CITEgeist/tests/test_profile_adapter.py
git commit -m "feat: add profile dict adapter for nested Module 3 format"
```

---

### Task 2: Argmax Global Inference Mode

Closes **Code Gap #2**: `pc_mil_infer_spot()` is hardcoded to Hungarian constrained assignment.

**Files:**
- Modify: `CITEgeist/model/pc_mil_inference.py:40-120`
- Test: `CITEgeist/tests/test_argmax_inference.py`

- [ ] **Step 1: Write the failing test**

```python
# CITEgeist/tests/test_argmax_inference.py
"""Tests for argmax global inference mode in PC-MIL."""
import numpy as np
import torch
import pytest
from unittest.mock import MagicMock
from model.pc_mil_inference import pc_mil_infer_spot


def _make_mock_model(n_nuclei: int, n_types: int):
    """Create a mock PCMILModel that returns predictable logits."""
    model = MagicMock()
    # Each nucleus strongly prefers a different type (round-robin)
    logits = torch.zeros(n_nuclei, n_types)
    for i in range(n_nuclei):
        logits[i, i % n_types] = 5.0  # Strong preference
    attention = torch.softmax(logits, dim=1)
    proportions = attention.mean(dim=0)
    reconstructed = torch.zeros(5)
    model.forward_with_logits.return_value = (logits, attention, proportions, reconstructed)
    model.eval = MagicMock(return_value=model)
    return model


def test_argmax_global_no_constraints():
    """argmax_global mode should assign each nucleus to its highest-scoring type."""
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
    # Round-robin: nucleus 0→TypeA, 1→TypeB, 2→TypeC, 3→TypeA, ...
    expected_types = [cell_type_names[i % n_types] for i in range(n_nuclei)]
    assert list(result["cell_type"]) == expected_types


def test_hungarian_constrained_is_default():
    """Default inference_mode should be hungarian_constrained (backward compat)."""
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
        # No inference_mode specified → default hungarian_constrained
    )

    assert len(result) == n_nuclei
    # With 0.75/0.25 proportions and 4 nuclei: expect 3 TypeA, 1 TypeB
    type_counts = result["cell_type"].value_counts()
    assert type_counts.get("TypeA", 0) == 3
    assert type_counts.get("TypeB", 0) == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_argmax_inference.py -v`
Expected: FAIL with `TypeError: pc_mil_infer_spot() got an unexpected keyword argument 'inference_mode'`

- [ ] **Step 3: Add `inference_mode` parameter to `pc_mil_infer_spot()`**

In `CITEgeist/model/pc_mil_inference.py`, modify the function signature and add argmax branch:

```python
def pc_mil_infer_spot(
    model,
    image_features,
    protein_proportions,
    detected_types,
    cell_type_names,
    morph_features=None,
    nucleus_ids=None,
    barcode=None,
    inference_mode="hungarian_constrained",  # NEW: "argmax_global" or "hungarian_constrained"
):
```

After computing logits and attention (existing code), add before the Hungarian block:

```python
    if inference_mode == "argmax_global":
        # Simple argmax: each nucleus → highest-scoring detected type
        assignments = attention_np.argmax(axis=1)
        confidences = attention_np.max(axis=1)
        assigned_types = [cell_type_names[a] for a in assignments]

        records = []
        for i in range(n_nuclei):
            records.append({
                "nucleus_id": nucleus_ids[i] if nucleus_ids else f"nucleus_{i}",
                "barcode": barcode or "",
                "cell_type": assigned_types[i],
                "confidence": float(confidences[i]),
                "protein_score": float(confidences[i]),  # Schema compat with hungarian mode
            })
        return pd.DataFrame(records)

    # Existing Hungarian constrained assignment code follows...
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd CITEgeist && python -m pytest tests/test_argmax_inference.py -v`
Expected: Both tests PASS

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/pc_mil_inference.py CITEgeist/tests/test_argmax_inference.py
git commit -m "feat: add argmax_global inference mode to PC-MIL"
```

---

### Task 3: Make Trainer val_dataset Optional

Closes **Code Gap #4**: Trainer requires `val_dataset` but spec says no train/val split.

**Files:**
- Modify: `CITEgeist/model/pc_mil_training.py:270-506`

- [ ] **Step 1: Modify `train_pc_mil()` signature**

Change `val_dataset: SpotDataset` to `val_dataset: Optional[SpotDataset] = None`.

- [ ] **Step 2: Add loss-plateau early stopping when val_dataset is None**

Inside the training loop, after computing epoch loss, add:

```python
    # Early stopping logic
    if val_dataset is not None:
        # Existing val_r-based early stopping (unchanged)
        ...
    else:
        # Loss-plateau early stopping (relative threshold)
        if epoch > recon_warmup_epochs:
            recent_losses = history["train_loss"][-patience:]
            if len(recent_losses) >= patience:
                mean_loss = np.mean(recent_losses)
                loss_std = np.std(recent_losses)
                if mean_loss > 0 and loss_std / mean_loss < 1e-3:
                    logger.info(f"Early stop: training loss plateaued (rel_std={loss_std/mean_loss:.2e})")
                    break

        # Save best model on training loss (since no val_r to track)
        if save_path and epoch_loss < best_train_loss:
            best_train_loss = epoch_loss
            torch.save(model.state_dict(), save_path)
```

Initialize `best_train_loss = float('inf')` before the training loop.

- [ ] **Step 3: Guard all val_dataset usage with `if val_dataset is not None:`**

Wrap the existing validation block (val_loss computation, val_r, best_val_r tracking) with the guard.

- [ ] **Step 4: Test with synthetic data**

Run: `cd CITEgeist && python -c "
import numpy as np, torch
from model.pc_mil import PCMILModel
from model.pc_mil_training import train_pc_mil, SpotDataset
ds = SpotDataset(
    features_per_spot=[np.random.randn(3, 384).astype('f') for _ in range(5)],
    protein_props=np.random.dirichlet([1]*3, 5).astype('f'),
    protein_signals=np.random.randn(5, 4).astype('f'),
    true_props=np.random.dirichlet([1]*3, 5).astype('f'),
)
model = PCMILModel(image_dim=384, n_types=3, n_markers=4)
h = train_pc_mil(model, ds, val_dataset=None, n_epochs=5)
print(f'ok, {len(h[\"train_loss\"])} epochs')
"`
Expected: `ok, 5 epochs` (no crash)

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/pc_mil_training.py
git commit -m "feat: make val_dataset optional in PC-MIL trainer"
```

---

### Task 4: Change ViT Default to ViT-S

**Files:**
- Modify: `CITEgeist/model/vit_extractor.py:50`

- [ ] **Step 1: Change default model_name**

In `ViTFeatureExtractor.__init__()`, change:

```python
# Before
model_name: str = 'vit_large_patch16_224',
# After
model_name: str = 'vit_small_patch16_224',
```

- [ ] **Step 2: Verify no downstream breakage**

Run: `cd CITEgeist && python -c "from model.vit_extractor import ViTFeatureExtractor; print('ok')"`

- [ ] **Step 3: Commit**

```bash
git add CITEgeist/model/vit_extractor.py
git commit -m "chore: default ViT extractor to vit_small_patch16_224 (384-dim)"
```

---

## Chunk 2: Unified Config + Marker Validation Module

### Task 5: Unified Pipeline Config

**Files:**
- Create: `CITEgeist/model/unified_config.py`

- [ ] **Step 1: Create config file with all shared constants**

```python
# CITEgeist/model/unified_config.py
"""Shared configuration for the unified PC-MIL pipeline."""

from pathlib import Path

# === Cell Type Profile (9 types, K=9) ===
# Nested format for Module 3 (validated by validate_cell_profile_dict)
CELL_PROFILES_NESTED = {
    "Cancer": {"Major": ["EPCAM-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1", "HLA-DRA-1"]},
}

CELL_TYPE_NAMES = list(CELL_PROFILES_NESTED.keys())
K = len(CELL_TYPE_NAMES)  # 9

# === RNA Marker Validation Dictionary ===
RNA_MARKERS = {
    "Cancer": ["EPCAM"],
    "Macrophages": ["CD68", "CD163", "CSF1R"],
    "CD8_T_Cells": ["CD8A", "CD8B", "GZMB"],
    "CD4_T_Cells": ["CD4", "IL7R", "CCR7"],
    "B_Cells": ["CD19", "MS4A1", "CD79A"],
    "Fibroblasts": ["COL1A1", "DCN", "VIM"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Monocytes": ["CD14", "FCGR3A", "S100A8"],
    "Dendritic_Cells": ["ITGAX", "HLA-DRA", "CLEC10A"],
}

# === Patient Sample List (12 canonical) ===
PATIENT_SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

# === Paths ===
DATA_DIR = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
OUTPUT_BASE = Path("output/unified_pipeline")

# === Training Defaults ===
MAX_EPOCHS = 200
PATIENCE = 30
LAMBDA_RECON = 1.0
LAMBDA_ENTROPY = 0.1
LAMBDA_DIVERSITY = 0.5
LAMBDA_HUNGARIAN = 0.0  # Disabled for unified pipeline
RECON_WARMUP_EPOCHS = 20
PROTEIN_DROPOUT = 0.3

# === ViT Config ===
VIT_MODEL = "vit_small_patch16_224"
VIT_DIM = 384
PATCH_SIZE = 224
```

- [ ] **Step 2: Commit**

```bash
git add CITEgeist/model/unified_config.py
git commit -m "feat: add unified pipeline config with 9-type profile and constants"
```

---

### Task 6: Marker Gene Validation Module

**Files:**
- Create: `CITEgeist/model/marker_validation.py`
- Test: `CITEgeist/tests/test_marker_validation.py`

- [ ] **Step 1: Write the failing test**

```python
# CITEgeist/tests/test_marker_validation.py
"""Tests for single-cell marker gene validation."""
import numpy as np
import pandas as pd
import pytest
from model.marker_validation import compute_marker_scores, summarize_validation


def test_compute_marker_scores_positive_case():
    """Nuclei assigned to correct type should have positive marker scores."""
    # Simulated deconvolved GEX: 2 spots × 2 types × 3 genes
    # Genes: CD68, EPCAM, VIM
    gex_data = {
        "spot_A:::Macrophages": {"CD68": 5.0, "EPCAM": 0.1, "VIM": 0.2},
        "spot_A:::Cancer": {"CD68": 0.1, "EPCAM": 4.0, "VIM": 0.3},
    }
    gex_df = pd.DataFrame(gex_data).T

    assignments = pd.DataFrame({
        "nucleus_id": ["n1", "n2"],
        "barcode": ["spot_A", "spot_A"],
        "cell_type": ["Macrophages", "Cancer"],
    })

    rna_markers = {
        "Macrophages": ["CD68"],
        "Cancer": ["EPCAM"],
    }

    scores = compute_marker_scores(assignments, gex_df, rna_markers)

    assert len(scores) == 2
    assert scores.loc[scores["nucleus_id"] == "n1", "marker_score"].values[0] > 0
    assert scores.loc[scores["nucleus_id"] == "n2", "marker_score"].values[0] > 0


def test_compute_marker_scores_negative_case():
    """Nuclei assigned to wrong type should have negative marker scores."""
    gex_data = {
        "spot_A:::Macrophages": {"CD68": 5.0, "EPCAM": 0.1},
        "spot_A:::Cancer": {"CD68": 0.1, "EPCAM": 4.0},
    }
    gex_df = pd.DataFrame(gex_data).T

    # Deliberately wrong assignments
    assignments = pd.DataFrame({
        "nucleus_id": ["n1", "n2"],
        "barcode": ["spot_A", "spot_A"],
        "cell_type": ["Cancer", "Macrophages"],  # Swapped
    })

    rna_markers = {
        "Macrophages": ["CD68"],
        "Cancer": ["EPCAM"],
    }

    scores = compute_marker_scores(assignments, gex_df, rna_markers)
    # Wrong assignment: marker for assigned type should be LOW
    assert scores.loc[scores["nucleus_id"] == "n1", "marker_score"].values[0] < 0
    assert scores.loc[scores["nucleus_id"] == "n2", "marker_score"].values[0] < 0


def test_summarize_validation():
    """Summary should report per-type and overall metrics."""
    scores_df = pd.DataFrame({
        "nucleus_id": ["n1", "n2", "n3", "n4"],
        "barcode": ["s1", "s1", "s2", "s2"],
        "cell_type": ["Cancer", "Cancer", "Macrophages", "Macrophages"],
        "marker_score": [2.0, 1.5, 3.0, -0.5],
        "markers_above_others": [True, True, True, False],
    })

    summary = summarize_validation(scores_df)

    assert "per_type" in summary
    assert "overall" in summary
    assert summary["per_type"]["Cancer"]["fraction_correct"] == 1.0
    assert summary["per_type"]["Macrophages"]["fraction_correct"] == 0.5
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_marker_validation.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'model.marker_validation'`

- [ ] **Step 3: Implement marker validation module**

```python
# CITEgeist/model/marker_validation.py
"""Single-cell marker gene validation for PC-MIL assignments.

For each nucleus assigned a cell type, checks whether the deconvolved GEX
from its parent spot supports that assignment by comparing marker gene
expression for the assigned type vs other types.
"""
import logging
from typing import Dict, List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_marker_scores(
    assignments: pd.DataFrame,
    gex_df: pd.DataFrame,
    rna_markers: Dict[str, List[str]],
) -> pd.DataFrame:
    """Compute per-nucleus marker gene scores.

    Args:
        assignments: DataFrame with columns [nucleus_id, barcode, cell_type].
        gex_df: Deconvolved GEX DataFrame indexed as 'barcode:::cell_type',
                columns are gene names.
        rna_markers: Dict mapping cell type → list of expected RNA marker genes.

    Returns:
        DataFrame with columns [nucleus_id, barcode, cell_type,
        assigned_marker_mean, other_marker_mean, marker_score,
        markers_above_others].
    """
    records = []
    all_types = list(rna_markers.keys())
    # Collect all marker genes across all types
    all_marker_genes = set()
    for genes in rna_markers.values():
        all_marker_genes.update(genes)
    # Filter to genes present in GEX data
    available_genes = set(gex_df.columns) & all_marker_genes

    for _, row in assignments.iterrows():
        nid = row["nucleus_id"]
        barcode = row["barcode"]
        assigned_type = row["cell_type"]

        gex_key = f"{barcode}:::{assigned_type}"
        if gex_key not in gex_df.index:
            logger.debug(f"GEX key {gex_key} not found, skipping nucleus {nid}")
            continue

        gex_row = gex_df.loc[gex_key]

        # Marker expression for assigned type
        assigned_markers = [g for g in rna_markers.get(assigned_type, []) if g in available_genes]
        if not assigned_markers:
            continue
        assigned_mean = gex_row[assigned_markers].mean()

        # Marker expression for other types (negative control)
        other_means = []
        for other_type in all_types:
            if other_type == assigned_type:
                continue
            other_markers = [g for g in rna_markers.get(other_type, []) if g in available_genes]
            if other_markers:
                other_means.append(gex_row[other_markers].mean())

        other_mean = np.mean(other_means) if other_means else 0.0

        # Marker score: log fold-change (assigned vs others)
        marker_score = float(np.log2((assigned_mean + 1e-6) / (other_mean + 1e-6)))

        records.append({
            "nucleus_id": nid,
            "barcode": barcode,
            "cell_type": assigned_type,
            "assigned_marker_mean": float(assigned_mean),
            "other_marker_mean": float(other_mean),
            "marker_score": marker_score,
            "markers_above_others": assigned_mean > other_mean,
        })

    return pd.DataFrame(records)


def summarize_validation(scores_df: pd.DataFrame) -> dict:
    """Summarize marker validation results per type and overall.

    Args:
        scores_df: Output from compute_marker_scores().

    Returns:
        Dict with 'per_type' and 'overall' summary metrics.
    """
    per_type = {}
    for cell_type, group in scores_df.groupby("cell_type"):
        per_type[cell_type] = {
            "n_nuclei": len(group),
            "median_marker_score": float(group["marker_score"].median()),
            "mean_marker_score": float(group["marker_score"].mean()),
            "fraction_correct": float(group["markers_above_others"].mean()),
        }

    # Overall weighted average
    total_nuclei = len(scores_df)
    if total_nuclei > 0:
        overall_fraction = float(scores_df["markers_above_others"].mean())
        overall_median = float(scores_df["marker_score"].median())
    else:
        overall_fraction = 0.0
        overall_median = 0.0

    return {
        "per_type": per_type,
        "overall": {
            "n_nuclei": total_nuclei,
            "fraction_correct": overall_fraction,
            "median_marker_score": overall_median,
        },
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd CITEgeist && python -m pytest tests/test_marker_validation.py -v`
Expected: All 3 tests PASS

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/marker_validation.py CITEgeist/tests/test_marker_validation.py
git commit -m "feat: add single-cell marker gene validation module"
```

---

## Chunk 3: Pipeline Runner Scripts

### Task 7: Step 1 Runner — Module 3 Re-run

Re-runs Module 3 with 9-type profile, unknown disabled, Cellpose nuclei counts injected. Saves outputs under `output/unified_pipeline/{sample}/module3/`.

**Files:**
- Create: `CITEgeist/examples/run_unified_step1_module3.py`

- [ ] **Step 1: Create the runner script**

```python
# CITEgeist/examples/run_unified_step1_module3.py
"""Step 1: Re-run Module 3 with 9-type profile for unified pipeline.

Usage:
    python run_unified_step1_module3.py --sample HCC22-088-P1-S1
    python run_unified_step1_module3.py --sample HCC22-088-P1-S1 --modality he
    python run_unified_step1_module3.py --sample xenium_region_0 --modality dapi \
        --xenium-gex path/to/gex.h5ad --xenium-protein path/to/protein.h5ad
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import scanpy as sc
import squidpy as sq
from model import CitegeistModel
from model.unified_config import (
    CELL_PROFILES_NESTED, CELL_TYPE_NAMES, DATA_DIR, OUTPUT_BASE,
    PATIENT_SAMPLES,
)


def run_step1(sample_name: str, modality: str = "he",
              xenium_gex: str = None, xenium_protein: str = None):
    """Run Module 3 with unified 9-type profile."""
    output_dir = OUTPUT_BASE / sample_name / "module3"
    output_dir.mkdir(parents=True, exist_ok=True)

    marker_file = OUTPUT_BASE / sample_name / ".step1_complete"
    if marker_file.exists():
        logger.info(f"Step 1 already complete for {sample_name}, skipping")
        return

    if modality == "he":
        # Patient Visium: load combined, split
        sample_path = DATA_DIR / sample_name / "outs"
        logger.info(f"Loading patient data from {sample_path}")
        adata = sq.read.visium(
            str(sample_path),
            counts_file="filtered_feature_bc_matrix.h5",
            load_images=True,
            gex_only=False,
        )
        model = CitegeistModel(
            sample_name=sample_name,
            adata=adata,
            output_folder=str(output_dir),
        )
        model.split_adata()
    elif modality == "dapi":
        # Xenium pseudo-Visium: pre-split data
        logger.info(f"Loading Xenium data: {xenium_gex}, {xenium_protein}")
        adata_gex = sc.read_h5ad(xenium_gex)
        adata_cite = sc.read_h5ad(xenium_protein)
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=str(output_dir),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_cite,
        )
    else:
        raise ValueError(f"Unknown modality: {modality}")

    # Preprocess
    model.preprocess_gex()
    model.preprocess_antibody()

    # Load 9-type profile
    model.load_cell_profile_dict(CELL_PROFILES_NESTED)

    # Cellpose nuclei counts (cached masks saved for Step 2)
    if modality == "he":
        nuclei_counts = model.compute_spot_nuclei_counts_cellpose(
            resolution_mode="hires",
            use_gpu=False,
            save_masks=True,
        )
        logger.info(f"Cellpose: {nuclei_counts.sum():.0f} total nuclei across {len(nuclei_counts)} spots")

    # Run Pass 1: Cell proportions
    model.run_cell_proportion_model(
        validation_warn_only=True,
    )

    # Run Pass 2: Gene expression deconvolution
    model.run_cell_expression_pass1()

    # Write marker file
    marker_file.touch()
    logger.info(f"Step 1 complete for {sample_name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 1: Module 3")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--modality", default="he", choices=["he", "dapi"])
    parser.add_argument("--xenium-gex", help="Path to Xenium GEX h5ad (dapi mode)")
    parser.add_argument("--xenium-protein", help="Path to Xenium protein h5ad (dapi mode)")
    args = parser.parse_args()
    run_step1(args.sample, args.modality, args.xenium_gex, args.xenium_protein)
```

- [ ] **Step 2: Verify syntax**

Run: `python -c "import ast; ast.parse(open('CITEgeist/examples/run_unified_step1_module3.py').read()); print('ok')"`

- [ ] **Step 3: Commit**

```bash
git add CITEgeist/examples/run_unified_step1_module3.py
git commit -m "feat: add Step 1 runner for Module 3 with 9-type profile"
```

---

### Task 8: Step 2+3 Runner — Cellpose Patches + ViT Features

Reuses Cellpose masks from Step 1 (patient) or runs fresh (Xenium DAPI), extracts 224x224 patches, computes ViT-S features.

**Files:**
- Create: `CITEgeist/examples/run_unified_step2_features.py`

- [ ] **Step 1: Create the runner script**

```python
# CITEgeist/examples/run_unified_step2_features.py
"""Step 2+3: Cellpose patch extraction + ViT-S feature extraction.

Usage:
    python run_unified_step2_features.py --sample HCC22-088-P1-S1 --modality he
    python run_unified_step2_features.py --sample xenium_region_0 --modality dapi \
        --dapi-path /path/to/dapi.tif --boundary-path /path/to/boundary.tif
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.spatial import cKDTree

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from model.unified_config import OUTPUT_BASE, VIT_MODEL, VIT_DIM, PATCH_SIZE, DATA_DIR
from model.vit_extractor import ViTFeatureExtractor


def load_image_and_segment(sample_name, modality, dapi_path=None, boundary_path=None):
    """Load image + run/reuse Cellpose segmentation. Returns masks, centroids, image."""
    from cellpose import models

    cellpose_dir = OUTPUT_BASE / sample_name / "cellpose"
    cellpose_dir.mkdir(parents=True, exist_ok=True)

    # Cellpose v3/v4 compatibility (v4 removed models.Cellpose)
    def _get_cellpose_model(model_type, use_gpu):
        try:
            from cellpose import models
            return models.Cellpose(model_type=model_type, gpu=use_gpu)
        except AttributeError:
            from cellpose.models import CellposeModel
            return CellposeModel(model_type=model_type, gpu=use_gpu)

    if modality == "he":
        import squidpy as sq
        sample_path = DATA_DIR / sample_name / "outs"

        # Try to reuse masks from Step 1
        module3_dir = OUTPUT_BASE / sample_name / "module3"
        cached_masks = list(module3_dir.glob("*_cellpose_masks.npy"))

        adata = sq.read.visium(str(sample_path), counts_file="filtered_feature_bc_matrix.h5",
                               load_images=True, gex_only=True)
        spatial_key = list(adata.uns["spatial"].keys())[0]
        img = adata.uns["spatial"][spatial_key]["images"]["hires"]
        scale = adata.uns["spatial"][spatial_key]["scalefactors"]["tissue_hires_scalef"]
        spatial_coords = adata.obsm["spatial"] * scale
        barcodes = list(adata.obs_names)

        if cached_masks:
            logger.info(f"Reusing cached Cellpose masks from {cached_masks[0]}")
            masks = np.load(cached_masks[0])
        else:
            logger.info("No cached masks, running Cellpose on H&E")
            cp_model = _get_cellpose_model("cyto2", torch.cuda.is_available())
            masks, _, _, _ = cp_model.eval(img, channels=[0, 0], diameter=None)

        image = img  # (H, W, 3) for H&E

    elif modality == "dapi":
        import tifffile
        logger.info(f"Loading DAPI: {dapi_path}, Boundary: {boundary_path}")
        dapi = tifffile.imread(dapi_path).astype(np.float32)
        boundary = tifffile.imread(boundary_path).astype(np.float32)
        image = np.stack([dapi, boundary], axis=-1)  # (H, W, 2)

        cp_model = _get_cellpose_model("nuclei", torch.cuda.is_available())
        masks, _, _, _ = cp_model.eval(image, channels=[1, 2], diameter=None)

        # Load spatial coords from pseudo-Visium AnnData (passed via --xenium-coords)
        # Deferred to Xenium follow-up task
        spatial_coords = None
        barcodes = None
    else:
        raise ValueError(f"Unknown modality: {modality}")

    # Extract centroids from masks
    from scipy.ndimage import center_of_mass
    nucleus_ids = np.unique(masks)
    nucleus_ids = nucleus_ids[nucleus_ids > 0]  # Remove background
    centroids = center_of_mass(masks > 0, masks, nucleus_ids)
    centroids_df = pd.DataFrame(centroids, columns=["y_pixel", "x_pixel"])
    centroids_df["nucleus_id"] = nucleus_ids

    # Save masks and centroids
    np.save(cellpose_dir / "nuclei_masks.npy", masks)
    centroids_df.to_csv(cellpose_dir / "nuclei_centroids.csv", index=False)

    return masks, centroids_df, image, spatial_coords, barcodes


def assign_nuclei_to_spots(centroids_df, spatial_coords, barcodes, spot_radius=None):
    """Assign each nucleus to nearest Visium spot within radius."""
    if spatial_coords is None:
        logger.warning("No spatial coords available, skipping spot assignment")
        return centroids_df

    tree = cKDTree(spatial_coords)
    nucleus_coords = centroids_df[["y_pixel", "x_pixel"]].values
    # Swap to (x, y) for KDTree if spatial_coords are (x, y)
    distances, indices = tree.query(nucleus_coords)

    centroids_df["spot_barcode"] = [barcodes[i] for i in indices]
    centroids_df["distance_to_spot"] = distances

    if spot_radius is not None:
        centroids_df = centroids_df[centroids_df["distance_to_spot"] <= spot_radius]
        logger.info(f"Kept {len(centroids_df)} nuclei within spot radius {spot_radius}")

    return centroids_df


def extract_patches(image, centroids_df, modality, patch_size=224):
    """Extract patches centered on each nucleus centroid."""
    h, w = image.shape[:2]
    half = patch_size // 2
    patches = []
    valid_ids = []

    for _, row in centroids_df.iterrows():
        cy, cx = int(row["y_pixel"]), int(row["x_pixel"])
        y1, y2 = cy - half, cy + half
        x1, x2 = cx - half, cx + half

        # Skip if patch would be out of bounds
        if y1 < 0 or x1 < 0 or y2 > h or x2 > w:
            continue

        patch = image[y1:y2, x1:x2]

        if modality == "dapi":
            # 2-channel → 3-channel (zero-pad third channel)
            if patch.ndim == 2:
                patch = patch[:, :, np.newaxis]
            if patch.shape[2] == 2:
                pad = np.zeros((patch_size, patch_size, 1), dtype=patch.dtype)
                patch = np.concatenate([patch, pad], axis=2)

        patches.append(patch)
        valid_ids.append(row["nucleus_id"])

    patches_arr = np.array(patches)  # (N, 224, 224, 3)
    logger.info(f"Extracted {len(patches_arr)} patches ({len(centroids_df) - len(patches_arr)} skipped, out of bounds)")
    return patches_arr, valid_ids


def extract_vit_features(patches, device="cuda"):
    """Extract ViT-S features from patches."""
    extractor = ViTFeatureExtractor(model_name=VIT_MODEL, device=device)

    # Convert from (N, H, W, C) to (N, C, H, W) and normalize to [0, 1]
    patches_chw = np.transpose(patches, (0, 3, 1, 2)).astype(np.float32)
    if patches_chw.max() > 1.0:
        patches_chw = patches_chw / 255.0

    features = extractor.extract_numpy(patches_chw, batch_size=64)
    logger.info(f"Extracted features: {features.shape}")
    return features


def run_step2(sample_name, modality="he", dapi_path=None, boundary_path=None):
    """Run Cellpose + patch extraction + ViT feature extraction."""
    step1_marker = OUTPUT_BASE / sample_name / ".step1_complete"
    if not step1_marker.exists():
        logger.error(f"Step 1 not complete for {sample_name}")
        return

    step2_marker = OUTPUT_BASE / sample_name / ".step2_complete"
    if step2_marker.exists():
        logger.info(f"Step 2 already complete for {sample_name}, skipping")
        return

    masks, centroids_df, image, spatial_coords, barcodes = load_image_and_segment(
        sample_name, modality, dapi_path, boundary_path
    )

    centroids_df = assign_nuclei_to_spots(centroids_df, spatial_coords, barcodes)

    patches, valid_ids = extract_patches(image, centroids_df, modality, PATCH_SIZE)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    features = extract_vit_features(patches, device)

    # Save outputs
    feat_dir = OUTPUT_BASE / sample_name / "features"
    feat_dir.mkdir(parents=True, exist_ok=True)
    np.save(feat_dir / "vit_features.npy", features)
    np.save(feat_dir / "nucleus_ids.npy", np.array(valid_ids))
    np.save(feat_dir / "patches.npy", patches)

    # Save updated centroids (with spot assignments)
    cellpose_dir = OUTPUT_BASE / sample_name / "cellpose"
    centroids_df.to_csv(cellpose_dir / "nuclei_centroids.csv", index=False)

    # Nuclei per spot
    if "spot_barcode" in centroids_df.columns:
        nps = centroids_df.groupby("spot_barcode").size().reset_index(name="n_nuclei")
        nps.to_csv(cellpose_dir / "nuclei_per_spot.csv", index=False)

    step2_marker.touch()
    logger.info(f"Step 2 complete for {sample_name}: {len(features)} nuclei")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 2+3: Features")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--modality", default="he", choices=["he", "dapi"])
    parser.add_argument("--dapi-path", help="DAPI TIFF path (dapi mode)")
    parser.add_argument("--boundary-path", help="Boundary TIFF path (dapi mode)")
    args = parser.parse_args()
    run_step2(args.sample, args.modality, args.dapi_path, args.boundary_path)
```

- [ ] **Step 2: Verify syntax**

Run: `python -c "import ast; ast.parse(open('CITEgeist/examples/run_unified_step2_features.py').read()); print('ok')"`

- [ ] **Step 3: Commit**

```bash
git add CITEgeist/examples/run_unified_step2_features.py
git commit -m "feat: add Step 2+3 runner for Cellpose patches + ViT features"
```

---

### Task 9: Step 4 Runner — PC-MIL Train + Infer

**Files:**
- Create: `CITEgeist/examples/run_unified_step3_pcmil.py`

- [ ] **Step 1: Create the runner script**

```python
# CITEgeist/examples/run_unified_step3_pcmil.py
"""Step 4: PC-MIL training + inference per sample.

Usage:
    python run_unified_step3_pcmil.py --sample HCC22-088-P1-S1
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from model.unified_config import (
    CELL_PROFILES_NESTED, CELL_TYPE_NAMES, OUTPUT_BASE, K,
    MAX_EPOCHS, PATIENCE, LAMBDA_RECON, LAMBDA_ENTROPY,
    LAMBDA_DIVERSITY, LAMBDA_HUNGARIAN, RECON_WARMUP_EPOCHS,
    PROTEIN_DROPOUT,
)
from model.pc_mil import PCMILModel, flatten_profile_dict, build_profile_matrix
from model.pc_mil_training import train_pc_mil, SpotDataset
from model.pc_mil_inference import pc_mil_infer_spot


def load_step2_outputs(sample_name):
    """Load features, centroids, and Module 3 proportions."""
    base = OUTPUT_BASE / sample_name

    features = np.load(base / "features" / "vit_features.npy")
    nucleus_ids = np.load(base / "features" / "nucleus_ids.npy")
    centroids = pd.read_csv(base / "cellpose" / "nuclei_centroids.csv")

    # Module 3 proportions
    m3_dir = base / "module3"
    prop_files = list(m3_dir.glob("*finetuned*.csv"))
    if not prop_files:
        prop_files = list(m3_dir.glob("*global*.csv"))
    props_df = pd.read_csv(prop_files[0], index_col=0)

    # Module 3 h5ad for antibody signals
    h5ad_files = list(m3_dir.glob("*.h5ad"))
    import scanpy as sc
    adata = sc.read_h5ad(h5ad_files[0])

    return features, nucleus_ids, centroids, props_df, adata


def build_spot_datasets(features, nucleus_ids, centroids, props_df, adata):
    """Build SpotDataset from per-nucleus features and Module 3 outputs."""
    # Map nuclei to spots
    if "spot_barcode" not in centroids.columns:
        raise ValueError("Centroids missing spot_barcode. Run Step 2 first.")

    # Filter centroids to valid nucleus_ids
    valid_mask = centroids["nucleus_id"].isin(nucleus_ids)
    centroids = centroids[valid_mask].copy()

    # Build per-spot feature lists
    spots = props_df.index.tolist()
    features_per_spot = []
    protein_props_list = []
    protein_signals_list = []

    # Get antibody data from adata
    has_antibody = hasattr(adata, "obsm") and "antibody" in adata.obsm
    flat_profile = flatten_profile_dict(CELL_PROFILES_NESTED)
    all_markers = []
    for markers in flat_profile.values():
        all_markers.extend(markers)
    all_markers = sorted(set(all_markers))

    nid_to_idx = {nid: i for i, nid in enumerate(nucleus_ids)}

    for spot in spots:
        spot_nuclei = centroids[centroids["spot_barcode"] == spot]
        if len(spot_nuclei) == 0:
            continue

        # Get feature indices for this spot's nuclei
        spot_nids = spot_nuclei["nucleus_id"].values
        feat_indices = [nid_to_idx[nid] for nid in spot_nids if nid in nid_to_idx]
        if not feat_indices:
            continue

        spot_features = features[feat_indices]
        features_per_spot.append(spot_features)

        # Proportions for this spot
        if spot in props_df.index:
            spot_props = props_df.loc[spot, CELL_TYPE_NAMES].values.astype(np.float32)
        else:
            continue
        protein_props_list.append(spot_props)

        # Protein signal (CLR-normalized antibody)
        if spot in adata.obs_names:
            spot_idx = list(adata.obs_names).index(spot)
            if has_antibody:
                signal = adata.obsm["antibody"][spot_idx]
            else:
                signal = np.zeros(len(all_markers), dtype=np.float32)
        else:
            signal = np.zeros(len(all_markers), dtype=np.float32)
        protein_signals_list.append(signal.astype(np.float32))

    protein_props = np.array(protein_props_list)
    protein_signals = np.array(protein_signals_list)

    dataset = SpotDataset(
        features_per_spot=features_per_spot,
        protein_props=protein_props,
        protein_signals=protein_signals,
        true_props=protein_props,  # Self-supervised: Module 3 props = target
    )

    return dataset


def run_step3(sample_name):
    """Train PC-MIL and run inference."""
    step2_marker = OUTPUT_BASE / sample_name / ".step2_complete"
    if not step2_marker.exists():
        logger.error(f"Step 2 not complete for {sample_name}")
        return

    step3_marker = OUTPUT_BASE / sample_name / ".step3_complete"
    if step3_marker.exists():
        logger.info(f"Step 3 already complete for {sample_name}, skipping")
        return

    pcmil_dir = OUTPUT_BASE / sample_name / "pcmil"
    pcmil_dir.mkdir(parents=True, exist_ok=True)

    features, nucleus_ids, centroids, props_df, adata = load_step2_outputs(sample_name)
    dataset = build_spot_datasets(features, nucleus_ids, centroids, props_df, adata)

    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Get marker info
    flat_profile = flatten_profile_dict(CELL_PROFILES_NESTED)
    all_markers = sorted(set(m for ms in flat_profile.values() for m in ms))

    # Build initial profile matrix
    profile_matrix = build_profile_matrix(CELL_PROFILES_NESTED, all_markers)
    init_profile = torch.tensor(profile_matrix, dtype=torch.float32)

    model = PCMILModel(
        image_dim=384,
        n_types=K,
        n_markers=len(all_markers),
        image_proj_dim=64,
        protein_context_dim=32,
        hidden_dim=128,
        init_profile_matrix=init_profile,
    )

    # Train (no val_dataset — uses loss plateau stopping)
    history = train_pc_mil(
        model=model,
        train_dataset=dataset,
        val_dataset=None,
        n_epochs=MAX_EPOCHS,
        lr=1e-3,
        lambda_recon=LAMBDA_RECON,
        lambda_entropy=LAMBDA_ENTROPY,
        lambda_diversity=LAMBDA_DIVERSITY,
        lambda_hungarian=LAMBDA_HUNGARIAN,
        patience=PATIENCE,
        recon_warmup_epochs=RECON_WARMUP_EPOCHS,
        protein_dropout=PROTEIN_DROPOUT,
        device=device,
        save_path=str(pcmil_dir / "model_weights.pt"),
    )

    # Save training log
    with open(pcmil_dir / "training_log.json", "w") as f:
        json.dump({k: [float(v) for v in vs] for k, vs in history.items()}, f)

    # Inference: argmax global, per spot
    model.eval()
    all_assignments = []

    for i in range(len(dataset)):
        sample = dataset[i]
        img_feats = sample["image_features"].to(device)
        prot_props = sample["protein_props"].to(device)
        detected = np.ones(K, dtype=bool)  # All types detected (no GMM gating)

        result = pc_mil_infer_spot(
            model=model,
            image_features=img_feats,
            protein_proportions=prot_props,
            detected_types=detected,
            cell_type_names=CELL_TYPE_NAMES,
            inference_mode="argmax_global",
        )
        all_assignments.append(result)

    assignments_df = pd.concat(all_assignments, ignore_index=True)
    assignments_df.to_csv(pcmil_dir / "assignments.csv", index=False)

    step3_marker.touch()
    logger.info(f"Step 3 complete for {sample_name}: {len(assignments_df)} nuclei assigned")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 4: PC-MIL")
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    run_step3(args.sample)
```

- [ ] **Step 2: Verify syntax**

Run: `python -c "import ast; ast.parse(open('CITEgeist/examples/run_unified_step3_pcmil.py').read()); print('ok')"`

- [ ] **Step 3: Commit**

```bash
git add CITEgeist/examples/run_unified_step3_pcmil.py
git commit -m "feat: add Step 4 runner for PC-MIL train + infer"
```

---

### Task 10: Step 5 Runner — Marker Gene Validation

**Files:**
- Create: `CITEgeist/examples/run_unified_step4_validate.py`

- [ ] **Step 1: Create the runner script**

```python
# CITEgeist/examples/run_unified_step4_validate.py
"""Step 5: Marker gene validation of PC-MIL assignments.

Usage:
    python run_unified_step4_validate.py --sample HCC22-088-P1-S1
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from model.unified_config import OUTPUT_BASE, RNA_MARKERS
from model.marker_validation import compute_marker_scores, summarize_validation


def run_step4(sample_name):
    """Run marker gene validation."""
    step3_marker = OUTPUT_BASE / sample_name / ".step3_complete"
    if not step3_marker.exists():
        logger.error(f"Step 3 not complete for {sample_name}")
        return

    val_dir = OUTPUT_BASE / sample_name / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)

    # Load assignments
    assignments = pd.read_csv(OUTPUT_BASE / sample_name / "pcmil" / "assignments.csv")

    # Load deconvolved GEX
    m3_dir = OUTPUT_BASE / sample_name / "module3"
    parquet_files = list(m3_dir.glob("*gene_expression_pass1.parquet"))
    if not parquet_files:
        logger.error(f"No GEX parquet found in {m3_dir}")
        return
    gex_df = pd.read_parquet(parquet_files[0])

    # Compute scores
    scores = compute_marker_scores(assignments, gex_df, RNA_MARKERS)
    scores.to_csv(val_dir / "marker_gene_scores.csv", index=False)

    # Summarize
    summary = summarize_validation(scores)
    with open(val_dir / "validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Print summary
    logger.info(f"=== Validation Summary for {sample_name} ===")
    logger.info(f"Total nuclei scored: {summary['overall']['n_nuclei']}")
    logger.info(f"Overall fraction correct: {summary['overall']['fraction_correct']:.3f}")
    logger.info(f"Overall median marker score: {summary['overall']['median_marker_score']:.3f}")
    for ct, metrics in summary["per_type"].items():
        logger.info(f"  {ct}: {metrics['fraction_correct']:.3f} correct "
                     f"({metrics['n_nuclei']} nuclei, median score={metrics['median_marker_score']:.2f})")

    step4_marker = OUTPUT_BASE / sample_name / ".step4_complete"
    step4_marker.touch()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 5: Validation")
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    run_step4(args.sample)
```

- [ ] **Step 2: Verify syntax**

Run: `python -c "import ast; ast.parse(open('CITEgeist/examples/run_unified_step4_validate.py').read()); print('ok')"`

- [ ] **Step 3: Commit**

```bash
git add CITEgeist/examples/run_unified_step4_validate.py
git commit -m "feat: add Step 5 runner for marker gene validation"
```

---

## Chunk 4: SLURM Array Job Scripts

### Task 11: SLURM Scripts for All 4 Steps

**Files:**
- Create: `CITEgeist/examples/sbatch_unified_step1.sh`
- Create: `CITEgeist/examples/sbatch_unified_step2.sh`
- Create: `CITEgeist/examples/sbatch_unified_step3.sh`
- Create: `CITEgeist/examples/sbatch_unified_step4.sh`

- [ ] **Step 1: Create Step 1 SLURM script (Module 3, CPU)**

```bash
#!/bin/bash
#SBATCH --job-name=unified_step1
#SBATCH --array=0-11
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --output=output/unified_pipeline/logs/step1_%A_%a.out
#SBATCH --error=output/unified_pipeline/logs/step1_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Running Step 1 for ${SAMPLE}"

mkdir -p output/unified_pipeline/logs

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist
python examples/run_unified_step1_module3.py --sample ${SAMPLE} --modality he
```

- [ ] **Step 2: Create Step 2 SLURM script (Features, GPU)**

```bash
#!/bin/bash
#SBATCH --job-name=unified_step2
#SBATCH --array=0-11
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --output=output/unified_pipeline/logs/step2_%A_%a.out
#SBATCH --error=output/unified_pipeline/logs/step2_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Gate on Step 1 completion
MARKER="output/unified_pipeline/${SAMPLE}/.step1_complete"
if [ ! -f "${MARKER}" ]; then
    echo "Step 1 not complete for ${SAMPLE}, exiting"
    exit 1
fi

echo "Running Step 2 for ${SAMPLE}"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist
python examples/run_unified_step2_features.py --sample ${SAMPLE} --modality he
```

- [ ] **Step 3: Create Step 3 SLURM script (PC-MIL, GPU)**

```bash
#!/bin/bash
#SBATCH --job-name=unified_step3
#SBATCH --array=0-11
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --output=output/unified_pipeline/logs/step3_%A_%a.out
#SBATCH --error=output/unified_pipeline/logs/step3_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

MARKER="output/unified_pipeline/${SAMPLE}/.step2_complete"
if [ ! -f "${MARKER}" ]; then
    echo "Step 2 not complete for ${SAMPLE}, exiting"
    exit 1
fi

echo "Running Step 3 for ${SAMPLE}"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist
python examples/run_unified_step3_pcmil.py --sample ${SAMPLE}
```

- [ ] **Step 4: Create Step 4 SLURM script (Validation, CPU)**

```bash
#!/bin/bash
#SBATCH --job-name=unified_step4
#SBATCH --array=0-11
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --output=output/unified_pipeline/logs/step4_%A_%a.out
#SBATCH --error=output/unified_pipeline/logs/step4_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

MARKER="output/unified_pipeline/${SAMPLE}/.step3_complete"
if [ ! -f "${MARKER}" ]; then
    echo "Step 3 not complete for ${SAMPLE}, exiting"
    exit 1
fi

echo "Running Step 4 for ${SAMPLE}"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist
python examples/run_unified_step4_validate.py --sample ${SAMPLE}
```

- [ ] **Step 5: Create logs directory and commit all scripts**

```bash
mkdir -p output/unified_pipeline/logs
git add CITEgeist/examples/sbatch_unified_step{1,2,3,4}.sh
git commit -m "feat: add SLURM array job scripts for unified pipeline (4 steps)"
```

---

## Chunk 5: Execution Sequence

### Task 12: Submit and Monitor

This task is manual — not automated. Run after all code is committed and tests pass.

- [ ] **Step 1: Submit Step 1 (Module 3)**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
sbatch CITEgeist/examples/sbatch_unified_step1.sh
# Note the job ID, e.g., 1620000
```

- [ ] **Step 2: Monitor Step 1 completion**

Wait for all 12 array tasks. Check:
```bash
ls output/unified_pipeline/HCC22-088-*/. step1_complete | wc -l
# Should be 12
```

- [ ] **Step 3: Submit Step 2 (Features)**

```bash
sbatch CITEgeist/examples/sbatch_unified_step2.sh
```

- [ ] **Step 4: Monitor and submit Steps 3 and 4 sequentially**

Repeat the pattern: wait for marker files, submit next step.

- [ ] **Step 5: Review validation results**

```bash
# Quick summary across all samples
for d in output/unified_pipeline/HCC22-088-*/validation/validation_summary.json; do
    sample=$(echo $d | cut -d/ -f3)
    frac=$(python -c "import json; d=json.load(open('$d')); print(f\"{d['overall']['fraction_correct']:.3f}\")")
    echo "$sample: $frac"
done
```
