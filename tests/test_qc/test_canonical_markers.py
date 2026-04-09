"""Tests for canonical_markers.py."""

import pytest

from CITEgeist.model.qc import QCResult
from CITEgeist.model.qc.canonical_markers import CANONICAL_MARKERS, get_available_markers


@pytest.mark.unit
class TestCanonicalMarkers:
    """Tests for the CANONICAL_MARKERS table."""

    def test_all_expected_types_present(self, cell_types):
        """Every QC cell type (except Cancer) has markers."""
        for ct in cell_types:
            assert ct in CANONICAL_MARKERS, f"Missing markers for {ct}"

    def test_cancer_excluded(self):
        """Cancer is intentionally excluded (tissue-specific markers)."""
        assert "Cancer" not in CANONICAL_MARKERS

    def test_each_type_has_5_to_8_markers(self):
        """Each type has 5-8 markers per spec."""
        for ct, markers in CANONICAL_MARKERS.items():
            assert 5 <= len(markers) <= 8, f"{ct} has {len(markers)} markers"

    def test_no_duplicate_markers_within_type(self):
        """No duplicates within a single cell type."""
        for ct, markers in CANONICAL_MARKERS.items():
            assert len(markers) == len(set(markers)), f"Duplicates in {ct}"

    def test_get_available_markers_filters_missing(self):
        """get_available_markers only returns markers present in gene list."""
        available_genes = ["CD68", "CD163", "EPCAM", "FAKE_GENE"]
        result = get_available_markers("Macrophages", available_genes)
        assert "CD68" in result
        assert "CD163" in result
        assert "FAKE_GENE" not in result

    def test_get_available_markers_unknown_type(self):
        """Unknown cell type returns empty list."""
        result = get_available_markers("UnknownType", ["CD68"])
        assert result == []


@pytest.mark.unit
class TestQCResult:
    """Tests for the QCResult dataclass."""

    def test_qcresult_creation(self):
        """QCResult can be instantiated with required fields."""
        result = QCResult(metrics={"r": 0.9}, flags=["test flag"], figures={})
        assert result.metrics["r"] == 0.9
        assert len(result.flags) == 1
        assert isinstance(result.figures, dict)

    def test_qcresult_defaults(self):
        """QCResult default fields are independent across instances."""
        r1 = QCResult(metrics={})
        r2 = QCResult(metrics={})
        r1.flags.append("flag")
        assert len(r2.flags) == 0  # default_factory ensures independence
