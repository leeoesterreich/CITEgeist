"""End-to-end test: Module 3 proportions → Module 3b MIL assignment."""
import numpy as np
import pandas as pd
import pytest
import torch

from CITEgeist.model.morphology_backbone import DAPIBackbone
from CITEgeist.model.single_cell_mil import SingleCellMIL, train_mil
from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_mil


@pytest.fixture
def synthetic_xenium_data():
    """Create synthetic data mimicking Xenium benchmark structure."""
    torch.manual_seed(42)
    np.random.seed(42)

    n_spots = 20
    n_types = 7
    cell_types = [
        "B cells", "CD4+ T cells", "CD8+ T cells",
        "Macrophages", "Endothelial", "Epithelial", "Fibroblasts",
    ]

    # Generate proportions (Module 3 output)
    raw_props = np.random.dirichlet(np.ones(n_types) * 2, size=n_spots)
    spot_ids = [f"spot_{i}" for i in range(n_spots)]
    proportions = pd.DataFrame(raw_props, columns=cell_types)
    proportions['spot_id'] = spot_ids

    # Generate fake patches and nuclei
    embeddings = {}
    nuclei_records = []
    nuclei_counts = {}
    gt_assignments = {}

    for i, sid in enumerate(spot_ids):
        n_nuclei = np.random.randint(3, 15)
        nuclei_counts[sid] = n_nuclei
        # Fake 384-dim embeddings (as if backbone already ran)
        embeddings[sid] = np.random.randn(n_nuclei, 384).astype(np.float32)
        for j in range(n_nuclei):
            nid = f"{sid}_{j}"
            nuclei_records.append({'nucleus_id': nid, 'spot_id': sid})
            # Random ground truth
            gt_assignments[nid] = cell_types[np.random.randint(n_types)]

    nuclei_spot_map = pd.DataFrame(nuclei_records)
    nuclei_counts_series = pd.Series(nuclei_counts)

    return {
        'proportions': proportions,
        'embeddings': embeddings,
        'nuclei_spot_map': nuclei_spot_map,
        'nuclei_counts': nuclei_counts_series,
        'cell_types': cell_types,
        'gt_assignments': gt_assignments,
    }


def test_e2e_mil_pipeline(synthetic_xenium_data):
    """Full pipeline: embeddings → MIL train → Hungarian → assignments."""
    data = synthetic_xenium_data

    result = run_nucleus_assignment_mil(
        embeddings=data['embeddings'],
        nuclei_spot_map=data['nuclei_spot_map'],
        proportions=data['proportions'],
        nuclei_counts=data['nuclei_counts'],
        cell_types=data['cell_types'],
        n_epochs=10,
        lambda_prior=1.0,
        device="cpu",
    )

    # Verify all nuclei assigned
    total_nuclei = len(data['nuclei_spot_map'])
    assert len(result.assignments) == total_nuclei
    assert result.method == "mil"

    # Verify all assignments are valid cell types
    for nid, ct in result.assignments.items():
        assert ct in data['cell_types']

    # Verify attention matrix is returned
    assert result.assignment_probs is not None
    assert len(result.assignment_probs) == total_nuclei


def test_e2e_proportion_consistency(synthetic_xenium_data):
    """Re-aggregated proportions should correlate with input proportions."""
    data = synthetic_xenium_data

    result = run_nucleus_assignment_mil(
        embeddings=data['embeddings'],
        nuclei_spot_map=data['nuclei_spot_map'],
        proportions=data['proportions'],
        nuclei_counts=data['nuclei_counts'],
        cell_types=data['cell_types'],
        n_epochs=30,
        lambda_prior=1.0,
        device="cpu",
    )

    # Re-aggregate to proportions
    cell_types = data['cell_types']
    input_props = data['proportions'].set_index('spot_id')[cell_types]

    for spot_id in data['proportions']['spot_id']:
        spot_nuclei = data['nuclei_spot_map'][
            data['nuclei_spot_map']['spot_id'] == spot_id
        ]
        n = len(spot_nuclei)
        if n == 0:
            continue
        counts = np.zeros(len(cell_types))
        for nid in spot_nuclei['nucleus_id']:
            ct = result.assignments.get(nid)
            if ct in cell_types:
                counts[cell_types.index(ct)] += 1
        reagg_props = counts / n

        # Re-aggregated should roughly match input (within discretization error)
        input_p = input_props.loc[spot_id].values
        # RMSE should be reasonable (< 0.5 for 7 types with random embeddings)
        rmse = np.sqrt(np.mean((reagg_props - input_p) ** 2))
        assert rmse < 0.5, f"Spot {spot_id}: RMSE={rmse:.3f} too high"
