import numpy as np
import pytest

from CITEgeist.model.assignment.cellularity_utils import round_counts_largest_remainder


@pytest.mark.unit
def test_canonical_rounder_sums_to_N():
    out = round_counts_largest_remainder(np.array([0.4, 0.4, 0.2]), 10)
    assert out.sum() == 10
    assert out.tolist() == [4, 4, 2]


@pytest.mark.unit
def test_sace_per_gene_conservation():
    # SACE integerizes a (T,) allocation per gene against an integer spot count.
    from CITEgeist.model.gex.sace_gex import _largest_remainder_round

    alloc = np.array([1.2, 2.7, 0.1])
    total = 4
    out = _largest_remainder_round(alloc, total)
    assert int(out.sum()) == total  # per-gene mass conserved
