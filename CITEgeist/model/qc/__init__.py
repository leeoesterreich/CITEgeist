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
        raise ImportError(
            "Failed to import CITEgeist.model.qc.report. This usually means the CITEgeist "
            "installation is incomplete/corrupted, or one of its core dependencies "
            "(matplotlib, numpy, pandas, anndata) failed to import — try reinstalling with "
            "`pip install --force-reinstall citegeist`, or run "
            '`python -c "import CITEgeist.model.qc.report"` to see the underlying error.'
        ) from e
    return _run_qc(*args, **kwargs)
