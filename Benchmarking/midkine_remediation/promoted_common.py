"""Shared pieces of the two promoted v42 stages B17 (``control_pairs.py``) and B18 (``celltype_map.py``).

Both were exploratory checks (``exploratory_v42``); the analysis code is unchanged. Only the fail-closed numeric gate lives here, so the stage
modules do not import from the exploratory folder.
"""
from __future__ import annotations

import numpy as np


class GateFailure(RuntimeError):
    pass


def gate_equal(name: str, got: float, expected: float, rtol: float = 1e-9, atol: float = 1e-12) -> dict:
    """Fail-closed numeric gate; raises GateFailure on mismatch, returns the record on match."""
    ok = bool(np.isclose(got, expected, rtol=rtol, atol=atol))
    rec = {"name": name, "got": float(got), "expected": float(expected), "match": ok}
    if not ok:
        raise GateFailure(f"GATE FAILED {name}: got {got!r} expected {expected!r}")
    return rec
