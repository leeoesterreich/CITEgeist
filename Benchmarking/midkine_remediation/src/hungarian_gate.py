"""Independent validation gate for the v4 Hungarian per-cell rebuild (A12).

Gate 1 of spec section 9.3: every represented spot's assigned type-count vector
must EXACTLY equal an independently recomputed expectation.

The expectation reproduces the two-stage rounding that the production path
(CITEgeist/model/assignment/cell_assignment.py, deterministic Hungarian with no
morphology) applies. This module deliberately does NOT import any CITEgeist
code, so it cannot inherit a bug from the code it checks.

Production rule (read from source, 2026-10-14):

Stage 1 (``discretize_proportions``), per spot in ``cell_prop.index``:
    N = int(nuclei_count_raw[spot])          # truncation, not rounding
    N == 0  -> integer row = zeros
    else    -> largest_remainder(N * p_row, N)
Stage 2 (``assign_cells_to_types``), per spot with n_mapped >= 1 mapped cells:
    if stage-1 row sum != n_mapped:
        largest_remainder(row / max(row.sum(), 1.0) * n_mapped, n_mapped)
    else: stage-1 row is used as is.
Largest remainder (``round_counts_largest_remainder``):
    scaled = c * (N / c.sum())   (uniform N/T each if c.sum() == 0)
    floors = floor(scaled); the (N - sum floors) units go to the largest
    fractional remainders; TIES go to the LOWER type index (stable sort on
    -remainder over the column order of ``cell_prop``).
Labels: with no morphology scores, the spot's cells receive
    [type_0] * k_0 + [type_1] * k_1 + ... in column order, so the per-spot
    count vector is exactly the stage-2 integer row.

Float note: the scaling expression is kept in the same operation order as
production (``c * (N / c.sum())``), because near-ties are decided by float
rounding. A mathematically "cleaner" exact-rational version would disagree on
such spots and flag a correct file.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd


def largest_remainder(c: Sequence[float], n: int) -> np.ndarray:
    """Hamilton apportionment of ``n`` units over weights ``c`` (lower index wins ties)."""
    c = np.asarray(c, dtype=np.float64)
    t = len(c)
    n = int(n)
    if n <= 0:
        return np.zeros(t, dtype=np.int64)
    total = c.sum()
    scaled = c * (n / total) if total > 0 else np.full(t, n / t)
    floors = np.floor(scaled).astype(np.int64)
    rem = scaled - floors
    deficit = int(n - floors.sum())
    # Rank types by (-remainder, index): deterministic, lower index wins a tie.
    ranking = sorted(range(t), key=lambda j: (-rem[j], j))
    if deficit > 0:
        for j in ranking[:deficit]:
            floors[j] += 1
    elif deficit < 0:
        for j in sorted(range(t), key=lambda j: (rem[j], j))[:-deficit]:
            if floors[j] > 0:
                floors[j] -= 1
    return floors


def expected_counts(
    cell_prop: pd.DataFrame,
    nuclei_count_raw: pd.Series,
    n_mapped: pd.Series,
) -> pd.DataFrame:
    """Expected (spots x types) integer count matrix after both rounding stages.

    Args:
        cell_prop: stored QP proportions (spots x types), unrenormalised.
        nuclei_count_raw: per-spot StarDist nuclei count as produced by segmentation.
        n_mapped: per-spot number of cells that actually carry a label in the output.
            Spots absent from this series count as 0 mapped cells.
    """
    if len(cell_prop) == 0:
        raise ValueError("cell_prop is empty")
    missing = cell_prop.index.difference(nuclei_count_raw.index)
    if len(missing) > 0:
        raise KeyError(f"{len(missing)} spots in cell_prop lack nuclei_count_raw, e.g. {list(missing[:3])}")
    out = np.zeros(cell_prop.shape, dtype=np.int64)
    mapped = n_mapped.reindex(cell_prop.index).fillna(0).astype(int)
    for i, spot in enumerate(cell_prop.index):
        n_m = int(mapped.loc[spot])
        if n_m == 0:
            continue  # production skips spots with no mapped cells
        n_raw = int(nuclei_count_raw.loc[spot])  # int() truncates, as in production
        p = cell_prop.loc[spot].values.astype(np.float64)
        stage1 = largest_remainder(n_raw * p, n_raw) if n_raw > 0 else np.zeros(len(p), dtype=np.int64)
        if int(stage1.sum()) != n_m:
            row = stage1.astype(np.float64)
            stage1 = largest_remainder(row / max(row.sum(), 1.0) * n_m, n_m)
        out[i] = stage1
    return pd.DataFrame(out, index=cell_prop.index, columns=cell_prop.columns)


@dataclass
class GateResult:
    passed: bool
    n_spots_checked: int
    n_cells: int
    n_mismatch_spots: int
    reasons: List[str] = field(default_factory=list)
    mismatches: List[Dict] = field(default_factory=list)  # first few, for the report

    def to_dict(self) -> Dict:
        return {
            "passed": self.passed,
            "n_spots_checked": self.n_spots_checked,
            "n_cells": self.n_cells,
            "n_mismatch_spots": self.n_mismatch_spots,
            "reasons": self.reasons,
            "mismatches": self.mismatches,
        }


def verify_allocation(
    cell_prop: pd.DataFrame,
    nuclei_count_raw: pd.Series,
    n_cells_per_spot: pd.Series,
    cell_labels: pd.DataFrame,
    *,
    spot_col: str = "spot_barcode",
    type_col: str = "cell_type",
    max_report: int = 10,
) -> GateResult:
    """Compare per-spot label counts to the independently recomputed expectation.

    Fails closed: empty input, labels outside the type set, or spots absent from
    ``cell_prop`` are failures, not skips. Comparison is exact over EVERY spot of
    ``cell_prop`` (a spot with no labelled cells must have expected count 0).

    ``n_cells_per_spot`` is the number of cells each spot SHOULD contain, taken from a
    source independent of the labels (the segmentation nucleus->spot table). It is
    required: deriving it from the labels themselves would let a deleted or
    duplicated cell pass (found by mutation: dropping the last label row passed).
    The labelled total per spot must equal it exactly.
    """
    reasons: List[str] = []
    if cell_prop is None or len(cell_prop) == 0:
        return GateResult(False, 0, 0, 0, ["empty cell_prop"])
    if cell_labels is None or len(cell_labels) == 0:
        return GateResult(False, len(cell_prop), 0, 0, ["empty cell_labels (no cells to validate)"])

    types = list(cell_prop.columns)
    bad_types = sorted(set(cell_labels[type_col].astype(str)) - set(types))
    if bad_types:
        reasons.append(f"labels outside type set: {bad_types}")
    unknown_spots = set(cell_labels[spot_col].astype(str)) - set(map(str, cell_prop.index))
    if unknown_spots:
        reasons.append(f"{len(unknown_spots)} labelled spots absent from cell_prop, e.g. {sorted(unknown_spots)[:3]}")
    if reasons:
        return GateResult(False, len(cell_prop), len(cell_labels), 0, reasons)

    lab = cell_labels[[spot_col, type_col]].astype(str)
    observed = (
        lab.groupby([spot_col, type_col]).size().unstack(fill_value=0).reindex(index=cell_prop.index.astype(str), columns=types, fill_value=0)
    )
    n_mapped = n_cells_per_spot.reindex(cell_prop.index).fillna(0).astype(int)
    n_obs = pd.Series(observed.sum(axis=1).values, index=cell_prop.index)
    count_bad = n_obs.values != n_mapped.values
    if count_bad.any():
        reasons.append(
            f"{int(count_bad.sum())} spots have a labelled cell total != independent cell count, "
            f"e.g. {[str(x) for x in cell_prop.index[count_bad][:3]]}"
        )
    extra = n_cells_per_spot.index.difference(cell_prop.index)
    if len(extra) and int(n_cells_per_spot.loc[extra].sum()) > 0:
        reasons.append(f"{len(extra)} spots with cells are absent from cell_prop")
    try:
        expected = expected_counts(cell_prop, nuclei_count_raw, n_mapped)
    except (KeyError, ValueError) as exc:
        return GateResult(False, len(cell_prop), len(cell_labels), 0, [f"expected-count failure: {exc}"])

    diff = (observed.values != expected.values).any(axis=1)
    mism = []
    for i in np.where(diff)[0][:max_report]:
        mism.append(
            {
                "spot": str(cell_prop.index[i]),
                "nuclei_count_raw": int(nuclei_count_raw.loc[cell_prop.index[i]]),
                "n_mapped": int(n_mapped.iloc[i]),
                "observed": observed.iloc[i].tolist(),
                "expected": expected.iloc[i].tolist(),
            }
        )
    n_bad = int(diff.sum())
    if n_bad:
        reasons.append(f"{n_bad} spots differ from the independent expectation")
    return GateResult(not reasons, len(cell_prop), len(cell_labels), n_bad, reasons, mism)
