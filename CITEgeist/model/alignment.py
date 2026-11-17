"""Strict label alignment helpers for the proportion -> SACE boundaries."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp


def align_to(index: pd.Index, frame, *, what: str, columns: Optional[Sequence[str]] = None):
    """Reindex ``frame`` to ``index``; raise on any label mismatch.

    Raises ValueError on duplicate labels in either, missing labels, or extra
    labels. The two label sets must be equal; callers subset explicitly.
    If ``columns`` is given, the frame's columns must be the same set and are
    reordered to ``columns`` (ValueError names missing and extra types).
    """
    index = pd.Index(index)
    if index.has_duplicates:
        raise ValueError(f"{what}: target index has duplicate labels")
    if frame.index.has_duplicates:
        raise ValueError(f"{what}: has duplicate labels")
    missing = index.difference(frame.index)
    if len(missing):
        raise ValueError(f"{what}: {len(missing)} labels missing (e.g. {list(missing[:3])})")
    extra = frame.index.difference(index)
    if len(extra):
        raise ValueError(f"{what}: {len(extra)} extra labels (e.g. {list(extra[:3])}); subset explicitly")
    out = frame.loc[index]
    if columns is not None:
        columns = list(columns)
        have = set(out.columns)
        want = set(columns)
        if have != want or len(columns) != len(want) or out.columns.has_duplicates:
            raise ValueError(f"{what}: type columns differ; missing {sorted(want - have)}, extra {sorted(have - want)}")
        out = out[columns]
    return out


def positions_of(source_index: pd.Index, wanted, *, what: str):
    """Integer positions of each label in ``wanted`` within ``source_index``.

    Raises ValueError on duplicate labels in either or on labels missing from ``source_index``.
    Extra labels in ``source_index`` are allowed (they are simply not selected).
    """
    source_index = pd.Index(source_index)
    wanted = pd.Index(wanted)
    if source_index.has_duplicates or wanted.has_duplicates:
        raise ValueError(f"{what}: duplicate spot labels")
    pos = source_index.get_indexer(wanted)
    missing = wanted[pos < 0]
    if len(missing):
        raise ValueError(f"{what}: {len(missing)} labels missing (e.g. {list(missing[:3])})")
    return pos


def densify(m) -> np.ndarray:
    """Return a dense float ndarray from a dense or scipy-sparse matrix."""
    if sp.issparse(m):
        return np.asarray(m.toarray(), dtype=float)
    return np.asarray(m, dtype=float)
