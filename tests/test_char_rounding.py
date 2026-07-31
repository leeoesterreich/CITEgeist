"""Layer 2 — characterization of the two shipped largest-remainder rounders.

Canonical: cellularity_utils.round_counts_largest_remainder (normalizes input to sum N,
           stable tie-break, handles over-allocation).
Duplicate: morphology_features.largest_remainder_discretize (NO normalization,
           non-stable tie-break, wraps deficit, no over-allocation handling).

Tier 5 E2b routes the duplicate through the canonical. They AGREE byte-for-byte on
the REAL caller contract — proportions summing to ~1.0 (any N>0) — gated below by
test_rounders_agree_on_probability_vectors (50 random simplex draws). OFF-CONTRACT
(input not summing to 1.0) the duplicate has a latent BUG: it never reduces
over-allocation, so e.g. dup([0.5,0.5,0.5], 6) -> [3,3,3] (sum 9 != 6) while canon
-> [2,2,2]. The reroute intentionally FIXES that, so it is documented here but
deliberately NOT pinned as must-preserve (Tier 5 Task 5.1 first proves every caller
is in-contract, then the reroute keeps the pinned in-contract cases byte-identical).
"""
import numpy as np
import pytest

from CITEgeist.model.assignment.cellularity_utils import round_counts_largest_remainder as canon
from CITEgeist.model.morphology.morphology_features import largest_remainder_discretize as dup

# (name, input, N, canon_expected, dup_expected)  [captured on HEAD]
# All rows are IN-CONTRACT (canon_exp == dup_exp) so they stay byte-identical after
# the Task 5.1 reroute. The off-contract divergence (props summing != 1.0) is a
# duplicate bug that the reroute fixes; it is intentionally NOT pinned here.
CASES = [
    ("even_props_N10", [0.25, 0.25, 0.25, 0.25], 10, [3, 3, 2, 2], [3, 3, 2, 2]),
    ("tie_remainders", [0.5, 0.5], 5, [3, 2], [3, 2]),
    ("single_type", [1.0], 7, [7], [7]),
    ("N0", [0.3, 0.7], 0, [0, 0], [0, 0]),
    ("all_zero_N4", [0.0, 0.0, 0.0], 4, [2, 1, 1], [2, 1, 1]),
    ("props_sum_lt1_N6", [0.2, 0.2, 0.2], 6, [2, 2, 2], [2, 2, 2]),
    ("realistic_N12", [0.4, 0.35, 0.15, 0.1], 12, [5, 4, 2, 1], [5, 4, 2, 1]),
    ("three_way_tie_N4", [1 / 3, 1 / 3, 1 / 3], 4, [2, 1, 1], [2, 1, 1]),
]


@pytest.mark.unit
@pytest.mark.parametrize("name,c,N,canon_exp,dup_exp", CASES)
def test_rounders_characterization(name, c, N, canon_exp, dup_exp):
    """Pin canonical and duplicate rounder outputs on shared in-contract fixtures."""
    arr = np.asarray(c, dtype=float)
    assert canon(arr, N).tolist() == canon_exp
    assert dup(arr, N).tolist() == dup_exp


@pytest.mark.unit
def test_rounders_agree_on_probability_vectors():
    """Prove canon and dup agree byte-for-byte on the real caller contract (sum==1.0)."""
    # The provably-identical domain Tier 5 E2b may merge: sum==1.0, N>0.
    rng = np.random.default_rng(0)
    for _ in range(50):
        p = rng.dirichlet(np.ones(rng.integers(2, 6)))
        N = int(rng.integers(1, 30))
        assert canon(p, N).tolist() == dup(p, N).tolist()


@pytest.mark.unit
def test_largest_remainder_discretize_matches_canonical():
    """E2b gate: exported duplicate reproduces the canonical rounder byte-for-byte.

    Covers representative + all-tied-remainder edge inputs, and keeps the goldens
    captured on HEAD 02bf76ad.
    """
    from CITEgeist.model.assignment.cellularity_utils import round_counts_largest_remainder
    from CITEgeist.model.morphology.morphology_features import largest_remainder_discretize

    golden = {
        ((0.5, 0.5, 0.0), 4): [2, 2, 0],
        ((0.4, 0.35, 0.25), 5): [2, 2, 1],
        ((0.4, 0.35, 0.25), 1): [1, 0, 0],
        ((0.5, 0.3, 0.2), 0): [0, 0, 0],
        ((0.25, 0.25, 0.25, 0.25), 6): [2, 2, 1, 1],  # all-tied remainders
        ((1 / 3, 1 / 3, 1 / 3), 2): [1, 1, 0],
        ((0.1,) * 10, 5): [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
    }
    for (props, n), exp in golden.items():
        arr = np.asarray(props)
        got = largest_remainder_discretize(arr, n)
        assert got.tolist() == exp  # fixed HEAD golden
        assert got.tolist() == round_counts_largest_remainder(arr, n).tolist()  # equivalence
        assert int(got.sum()) == n
