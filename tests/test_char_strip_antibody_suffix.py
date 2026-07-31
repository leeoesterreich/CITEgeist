"""Layer 2 — characterization of the centralized strip_antibody_suffix (marker_utils:10)."""
import pytest

from CITEgeist.model.marker_utils import strip_antibody_suffix as strip

# (name, active_markers) -> expected  [captured on HEAD]
CASES = [
    (("EPCAM-1", None), "EPCAM"),
    (("CD68-1", None), "CD68"),
    (("PD-L1-1", None), "PD-L1"),
    (("PD-1", None), "PD-1"),  # natural -1 preserved
    (("EPCAM", None), "EPCAM"),  # no suffix
    (("CD3", None), "CD3"),
    (("CD8-1", frozenset({"CD8"})), "CD8"),  # caller-vocab: base present -> strip
    (("CD8-1", frozenset({"CD4"})), "CD8-1"),  # caller-vocab: base absent -> keep
    (("PD-1", frozenset({"PD"})), "PD-1"),  # preserved even in vocab mode
]


@pytest.mark.unit
@pytest.mark.parametrize("args,expected", CASES)
def test_strip_antibody_suffix_characterization(args, expected):
    """Pin strip_antibody_suffix output across natural/vocab-mode suffix cases."""
    assert strip(*args) == expected
