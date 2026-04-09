"""Smoke tests validating canonical benchmark outputs exist and meet performance thresholds.

These tests read pre-computed result files from the Benchmarking/ directory.
They run fast (JSON reads) and confirm that prior benchmark runs produced
sensible results. They do NOT re-run any benchmarks.

Mark: requires_data — skipped automatically when files are absent.
"""
import json
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
XENIUM_QP_SUMMARY = (
    REPO_ROOT
    / "Benchmarking/xenium_benchmarking/CITEgeist/output_qp_singler_7type/canonical_summary_7type.json"
)
XENIUM_SACE_AGGREGATE = (
    REPO_ROOT
    / "Benchmarking/xenium_benchmarking/CITEgeist/results_sace/sace_aggregate.json"
)

# ── QP proportion results ──────────────────────────────────────────────────

@pytest.mark.requires_data
def test_xenium_qp_summary_exists():
    if not XENIUM_QP_SUMMARY.exists():
        pytest.skip(f"Xenium QP summary not found: {XENIUM_QP_SUMMARY}")
    data = json.loads(XENIUM_QP_SUMMARY.read_text())
    assert "overall_mean_r" in data, "Missing overall_mean_r key"
    assert "per_region_r" in data, "Missing per_region_r key"


@pytest.mark.requires_data
def test_xenium_qp_overall_r_above_threshold():
    """Production QP r must be > 0.60 (well below expected ~0.72)."""
    if not XENIUM_QP_SUMMARY.exists():
        pytest.skip(f"Xenium QP summary not found: {XENIUM_QP_SUMMARY}")
    data = json.loads(XENIUM_QP_SUMMARY.read_text())
    mean_r = data["overall_mean_r"]
    assert mean_r > 0.60, f"Xenium QP overall_mean_r={mean_r:.4f} below threshold 0.60"


@pytest.mark.requires_data
def test_xenium_qp_all_regions_present():
    if not XENIUM_QP_SUMMARY.exists():
        pytest.skip(f"Xenium QP summary not found: {XENIUM_QP_SUMMARY}")
    data = json.loads(XENIUM_QP_SUMMARY.read_text())
    assert data.get("n_regions", 0) == 5, f"Expected 5 regions, got {data.get('n_regions')}"
    assert len(data["per_region_r"]) == 5, "Expected 5 per-region r values"


# ── SACE GEX results ───────────────────────────────────────────────────────

@pytest.mark.requires_data
def test_xenium_sace_aggregate_exists():
    if not XENIUM_SACE_AGGREGATE.exists():
        pytest.skip(f"Xenium SACE aggregate not found: {XENIUM_SACE_AGGREGATE}")
    data = json.loads(XENIUM_SACE_AGGREGATE.read_text())
    assert "methods" in data, "Missing methods key"
    assert "sace_iter0" in data["methods"], "Missing sace_iter0 in methods"


@pytest.mark.requires_data
def test_xenium_sace_percell_r_above_threshold():
    """Per-cell GEX r must be > 0.25 (well below expected ~0.31)."""
    if not XENIUM_SACE_AGGREGATE.exists():
        pytest.skip(f"Xenium SACE aggregate not found: {XENIUM_SACE_AGGREGATE}")
    data = json.loads(XENIUM_SACE_AGGREGATE.read_text())
    mean_r = data["methods"]["sace_iter0"]["percell_mean_gene_r"]["mean"]
    assert mean_r > 0.25, f"Xenium SACE percell_mean_gene_r={mean_r:.4f} below threshold 0.25"


@pytest.mark.requires_data
def test_xenium_sace_all_five_regions_loaded():
    if not XENIUM_SACE_AGGREGATE.exists():
        pytest.skip(f"Xenium SACE aggregate not found: {XENIUM_SACE_AGGREGATE}")
    data = json.loads(XENIUM_SACE_AGGREGATE.read_text())
    assert data.get("regions_missing") == [], f"Missing regions: {data.get('regions_missing')}"
    assert data.get("n_regions_total") == 5
