# tests/test_marker_utils.py
import pytest

from CITEgeist.model.marker_utils import strip_antibody_suffix


@pytest.mark.unit
@pytest.mark.parametrize("name", ["EPCAM-1", "CD14-1", "MS4A1-1", "CD19-1", "PECAM1-1", "ACTA2-1", "CD8A-1"])
def test_spaceranger_suffix_is_stripped(name):
    assert strip_antibody_suffix(name) == name[:-2]


@pytest.mark.unit
def test_natural_dash_one_marker_preserved():
    assert strip_antibody_suffix("PD-1") == "PD-1"


@pytest.mark.unit
def test_non_suffixed_name_unchanged():
    assert strip_antibody_suffix("CD68") == "CD68"
    assert strip_antibody_suffix("PD-L1") == "PD-L1"


@pytest.mark.unit
def test_active_markers_scopes_stripping():
    # Conservative mode: strip only when the de-suffixed base is in the caller vocabulary.
    assert strip_antibody_suffix("CD14-1", active_markers={"CD14"}) == "CD14"
    assert strip_antibody_suffix("XYZ-1", active_markers={"CD14"}) == "XYZ-1"  # unknown base preserved
    assert strip_antibody_suffix("PD-1", active_markers={"PD"}) == "PD-1"  # natural name still preserved
