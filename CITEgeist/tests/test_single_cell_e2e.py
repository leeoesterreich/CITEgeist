"""End-to-end test for single-cell resolution on synthetic data."""
import numpy as np
import pandas as pd
import pytest

from CITEgeist.model.morphology_features import extract_nucleus_features
from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment


@pytest.mark.slow
def test_morphology_extraction_on_real_mask():
    """Test morphology extraction on realistic synthetic mask."""
    np.random.seed(42)
    mask = np.zeros((500, 500), dtype=np.int32)

    nucleus_id = 1
    for i in range(10):
        for j in range(10):
            cx, cy = 25 + i * 45, 25 + j * 45
            if nucleus_id % 3 == 0:
                # Elongated
                mask[cy-5:cy+5, cx-15:cx+15] = nucleus_id
            elif nucleus_id % 3 == 1:
                # Small round
                y, x = np.ogrid[:500, :500]
                circle = ((x - cx)**2 + (y - cy)**2) <= 64
                mask[circle] = nucleus_id
            else:
                # Large round
                y, x = np.ogrid[:500, :500]
                circle = ((x - cx)**2 + (y - cy)**2) <= 144
                mask[circle] = nucleus_id
            nucleus_id += 1

    features = extract_nucleus_features(mask)

    assert len(features) == 100
    # Check feature variation
    assert features['area'].std() > 0
    assert features['circularity'].std() > 0

    # Elongated nuclei should have lower circularity
    elongated_ids = [i for i in range(1, 101) if i % 3 == 0]
    round_ids = [i for i in range(1, 101) if i % 3 != 0]

    elongated_circ = features[features['nucleus_id'].isin(elongated_ids)]['circularity'].mean()
    round_circ = features[features['nucleus_id'].isin(round_ids)]['circularity'].mean()

    assert elongated_circ < round_circ


@pytest.mark.slow
def test_full_assignment_pipeline():
    """Test full assignment pipeline with synthetic data."""
    np.random.seed(123)

    # Create mask with 50 nuclei across 10 spots
    mask = np.zeros((200, 200), dtype=np.int32)
    nucleus_id = 1
    nuclei_spot_data = []

    for spot_i in range(10):
        spot_x = 20 + (spot_i % 5) * 35
        spot_y = 20 + (spot_i // 5) * 80
        spot_id = f'spot_{spot_i}'

        # 5 nuclei per spot
        for n in range(5):
            nx = spot_x + (n % 3) * 10 - 10
            ny = spot_y + (n // 3) * 15

            # Create nucleus
            y, x = np.ogrid[:200, :200]
            circle = ((x - nx)**2 + (y - ny)**2) <= 25
            mask[circle] = nucleus_id

            nuclei_spot_data.append({
                'nucleus_id': nucleus_id,
                'spot_id': spot_id,
            })
            nucleus_id += 1

    nuclei_spot_map = pd.DataFrame(nuclei_spot_data)

    # Create proportions (biased toward different types per spot)
    proportions_data = []
    for i in range(10):
        if i < 5:
            props = [0.6, 0.3, 0.1]  # mostly type 0
        else:
            props = [0.1, 0.3, 0.6]  # mostly type 2
        proportions_data.append({
            'spot_id': f'spot_{i}',
            'TypeA': props[0],
            'TypeB': props[1],
            'TypeC': props[2],
        })
    proportions = pd.DataFrame(proportions_data)

    nuclei_counts = pd.Series([5] * 10, index=[f'spot_{i}' for i in range(10)])

    result = run_nucleus_assignment(
        mask=mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=['TypeA', 'TypeB', 'TypeC'],
    )

    # All nuclei should be assigned
    assert len(result.assignments) == 50

    # Check distribution roughly matches proportions
    # Spots 0-4 should have more TypeA
    spot_0_4_nuclei = [nid for nid, sid in zip(
        nuclei_spot_map['nucleus_id'], nuclei_spot_map['spot_id']
    ) if int(sid.split('_')[1]) < 5]

    type_a_count = sum(1 for nid in spot_0_4_nuclei if result.assignments[nid] == 'TypeA')
    assert type_a_count >= 10  # At least 10 of 25 should be TypeA (expected ~15)
