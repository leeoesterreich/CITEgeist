"""
QC module for validating CITEgeist single-cell outputs.

Two modes:
- benchmark: Xenium with ground truth (proportion + GEX error characterization)
- self_consistency: Real patient data (marker enrichment + cross-patient consistency)
"""

from __future__ import annotations

from ._types import QCResult

__all__ = ["QCResult", "run_qc"]


def run_qc(*args, **kwargs):
    """Orchestrate all QC checks. See report.run_qc for full signature."""
    try:
        from .report import run_qc as _run_qc  # pylint: disable=import-outside-toplevel
    except ImportError as e:
        raise NotImplementedError(
            "report module is not yet available. " "run_qc will be fully implemented in Task 6."
        ) from e
    return _run_qc(*args, **kwargs)
