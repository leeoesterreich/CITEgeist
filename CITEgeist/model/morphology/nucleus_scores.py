"""Non-DL morphology score producer for single-cell assignment.

Turns a labeled nuclei mask + per-nucleus spot membership + spot-level
proportions into a (C, n_types) score matrix via skimage regionprops shape
features and an LLP soft-label classifier. Zero torch.
"""

from typing import Optional

import numpy as np

from .morphology_features import extract_nucleus_features
from .soft_label_classifier import SoftLabelClassifier

# The eight real shape columns emitted by extract_nucleus_features
# (verified against live source — NOT extent/nc_ratio).
SHAPE_COLS = [
    "area",
    "perimeter",
    "circularity",
    "eccentricity",
    "solidity",
    "major_axis_length",
    "minor_axis_length",
    "aspect_ratio",
]


def _uniform(n_rows: int, n_types: int) -> np.ndarray:
    if n_rows == 0:
        return np.zeros((0, n_types), dtype=float)
    return np.full((n_rows, n_types), 1.0 / n_types, dtype=float)


def _is_degenerate(y_soft: np.ndarray) -> bool:
    """True when fewer than two cell types carry any mass (LogReg can't fit)."""
    mass = y_soft.sum(axis=0)
    return int((mass > 1e-10).sum()) < 2


def _reindex_to_n_types(p_seen: np.ndarray, seen_classes: np.ndarray, n_types: int) -> np.ndarray:
    """Scatter predict_proba columns (only observed classes) into full width."""
    full = np.zeros((p_seen.shape[0], n_types), dtype=float)
    full[:, seen_classes] = p_seen
    return full


def compute_morphology_scores(
    labeled_mask: np.ndarray,
    nuclei_spot_map: np.ndarray,
    spot_proportions: np.ndarray,
    cell_type_names: list,
    *,
    nucleus_ids: Optional[np.ndarray] = None,
    classifier: str = "logreg",
) -> np.ndarray:
    """Compute per-nucleus (C, n_types) morphology scores.

    Args:
        labeled_mask: (H, W) int mask; each non-zero label is a nucleus.
        nuclei_spot_map: (C,) int spot index per nucleus, aligned to
            ``nucleus_ids`` (or to ascending mask-label order if None).
        spot_proportions: (I, n_types) QP-predicted spot proportions
            (the LLP soft labels — NOT ground-truth-derived).
        cell_type_names: length-n_types ordered type names.
        nucleus_ids: optional (C,) mask labels selecting/ordering the nuclei
            to score. Defaults to all labels in ascending order.
        classifier: "logreg" (default) or "rf".

    Returns:
        (C, n_types) float array, rows sum to 1.
    """
    n_types = len(cell_type_names)
    feats = extract_nucleus_features(labeled_mask).set_index("nucleus_id")

    if nucleus_ids is None:
        nucleus_ids = feats.index.to_numpy()
    nucleus_ids = np.asarray(nucleus_ids)
    nuclei_spot_map = np.asarray(nuclei_spot_map, dtype=int)

    n_cells = len(nucleus_ids)
    if n_cells == 0:
        return np.zeros((0, n_types), dtype=float)

    # Fail clearly on misalignment instead of silently median-filling phantom
    # cells or emitting a wrong-length matrix (Codex MEDIUM-5).
    if len(nuclei_spot_map) != n_cells:
        raise ValueError(f"nuclei_spot_map length {len(nuclei_spot_map)} != nucleus_ids length {n_cells}")
    if nuclei_spot_map.min() < 0 or nuclei_spot_map.max() >= spot_proportions.shape[0]:
        raise ValueError(f"nuclei_spot_map spot indices out of range for {spot_proportions.shape[0]} spots")
    missing = np.setdiff1d(nucleus_ids, feats.index.to_numpy())
    if len(missing):
        raise ValueError(f"{len(missing)} requested nucleus_ids absent from the mask: {missing[:5]}")

    feats = feats.reindex(nucleus_ids)
    X = feats[SHAPE_COLS].to_numpy(dtype=float)
    # Fill any missing/NaN feature rows with per-column medians (robust to a
    # requested id absent from the mask); if a whole column is NaN, use 0.
    col_median = np.nanmedian(X, axis=0)
    col_median = np.where(np.isfinite(col_median), col_median, 0.0)
    inds = np.where(~np.isfinite(X))
    X[inds] = np.take(col_median, inds[1])

    y_soft = np.asarray(spot_proportions, dtype=float)[nuclei_spot_map]

    if _is_degenerate(y_soft):
        return _uniform(n_cells, n_types)

    if classifier == "rf":
        from sklearn.ensemble import RandomForestClassifier

        y_hard = y_soft.argmax(axis=1)
        if len(np.unique(y_hard)) < 2:
            return _uniform(n_cells, n_types)
        rf = RandomForestClassifier(n_estimators=200, random_state=42)
        rf.fit(X, y_hard)
        return _reindex_to_n_types(rf.predict_proba(X), rf.classes_, n_types)

    clf = SoftLabelClassifier(n_types).fit(X, y_soft)
    return _reindex_to_n_types(clf.predict_proba(X), clf._model.classes_, n_types)


def compute_morphology_scores_safe(*args, **kwargs) -> Optional[np.ndarray]:
    """Wrap :func:`compute_morphology_scores` so failures degrade gracefully.

    Returns the ``(C, n_types)`` score matrix on success. If scoring raises for
    any reason (degenerate features, classifier failure, mask/id misalignment),
    logs a warning and returns ``None`` so the caller can fall back to a
    count-constrained assignment instead of crashing the pipeline. See
    ``docs/superpowers/specs/2026-07-12-morph-bayesian-default-design.md``.
    """
    import logging

    try:
        return compute_morphology_scores(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 — intentional broad fallback
        logging.warning(
            "Morphology scoring failed (%s); falling back to count-constrained assignment.",
            exc,
        )
        return None
