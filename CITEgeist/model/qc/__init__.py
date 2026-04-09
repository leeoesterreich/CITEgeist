"""
QC module for validating CITEgeist single-cell outputs.

Two modes:
- benchmark: Xenium with ground truth (proportion + GEX error characterization)
- self_consistency: Real patient data (marker enrichment + cross-patient consistency)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from matplotlib.figure import Figure

__all__ = ["QCResult", "run_qc"]


@dataclass
class QCResult:
    """Result container for a single QC module.

    Args:
        metrics: Nested dict of computed metrics.
        flags: Human-readable warnings/notes.
        figures: panel_id → matplotlib Figure for composition.
    """

    metrics: dict[str, Any]
    flags: list[str] = field(default_factory=list)
    figures: dict[str, Figure] = field(default_factory=dict)


def run_qc(*args, **kwargs):
    """Orchestrate all QC checks. See report.run_qc for full signature."""
    try:
        from .report import run_qc as _run_qc
    except ImportError as e:
        raise NotImplementedError(
            "report module is not yet available. "
            "run_qc will be fully implemented in Task 6."
        ) from e
    return _run_qc(*args, **kwargs)
