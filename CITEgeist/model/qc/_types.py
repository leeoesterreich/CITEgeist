"""Shared type definitions for the QC subpackage.

Placed in a separate module to avoid cyclic imports between qc/__init__.py
and qc/report.py (which imports from the submodules that re-import __init__).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from matplotlib.figure import Figure

__all__ = ["QCResult"]


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
