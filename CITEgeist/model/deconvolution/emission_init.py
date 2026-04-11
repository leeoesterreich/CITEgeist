# CITEgeist/model/emission_init.py
"""
Emission initialization for per-type beta QP.

Provides marker-type association table (ported from archived NB model),
beta matrix initialization, and per-pair regularization sigma.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

# Marker-type associations with strong/soft strength labels.
# Ported from archived nb_initialization.py MARKER_TYPE_TABLE.
# E-Cadherin EXCLUDED: anti-correlates with Epithelial in RCC (r=-0.37).
MARKER_TYPE_TABLE: Dict[str, List[Tuple[str, str]]] = {
    "CD20": [("B cells", "strong")],
    "CD3E": [("CD4+ T cells", "strong"), ("CD8+ T cells", "strong")],
    "CD4": [("CD4+ T cells", "strong")],
    "CD8A": [("CD8+ T cells", "strong")],
    "CD68": [("Macrophages", "strong")],
    "CD163": [("Macrophages", "strong")],
    "CD16": [("Macrophages", "strong")],
    "HLA-DR": [("B cells", "soft"), ("Macrophages", "strong")],
    "CD11c": [("Macrophages", "soft")],
    "CD31": [("Endothelial", "strong")],
    "PanCK": [("Epithelial", "strong")],
    "alphaSMA": [("Fibroblasts", "strong")],
    "Vimentin": [("Macrophages", "soft"), ("Endothelial", "soft"), ("Fibroblasts", "strong")],
    "CD45": [("B cells", "strong"), ("CD4+ T cells", "strong"), ("CD8+ T cells", "strong"), ("Macrophages", "strong")],
    "CD45RA": [("B cells", "strong")],
    "CD45RO": [("CD4+ T cells", "strong"), ("CD8+ T cells", "strong")],
    "CD138": [("B cells", "soft")],
}

FUNCTIONAL_MARKERS = {
    "PD-1",
    "PD-L1",
    "LAG-3",
    "VISTA",
    "Ki-67",
    "PCNA",
    "GranzymeB",
    "Granzyme B",
    "Beta-catenin",
    "PTEN",
    "PD1",
    "PDL1",
    "LAG3",
}


def _strip_suffix(name: str) -> str:
    """Strip -1 suffix from marker names (e.g., CD68-1 -> CD68)."""
    if name.endswith("-1"):
        base = name[:-2]
        if base in MARKER_TYPE_TABLE or base in FUNCTIONAL_MARKERS:
            return base
    return name


def build_marker_config(
    available_markers: List[str],
    marker_subset: Optional[List[str]] = None,
) -> Tuple[List[str], np.ndarray, List[str]]:
    """Build marker panel config from available markers in the data.

    Args:
        available_markers: Marker names from adata.var_names.
        marker_subset: If provided, only include these canonical marker names.
            Use for Sprint 1 (7 Major only) vs Sprint 2 (all 17).

    Returns:
        markers: List of canonical marker names found in data.
        active_mask: (T, M) bool array, True where type-marker pair is active.
        type_names: List of cell type names.
    """
    type_names = list(CELL_TYPES)
    T = len(type_names)

    markers = []
    seen = set()
    for raw_name in available_markers:
        canonical = _strip_suffix(raw_name)
        if canonical in FUNCTIONAL_MARKERS:
            continue
        if canonical not in MARKER_TYPE_TABLE:
            continue
        if marker_subset is not None and canonical not in marker_subset:
            continue
        if canonical in seen:
            continue
        markers.append(canonical)
        seen.add(canonical)

    M = len(markers)
    active_mask = np.zeros((T, M), dtype=bool)

    for m_idx, marker in enumerate(markers):
        for type_name, _strength in MARKER_TYPE_TABLE.get(marker, []):
            if type_name in type_names:
                t_idx = type_names.index(type_name)
                active_mask[t_idx, m_idx] = True

    logger.info("build_marker_config: %d markers, %d types, %d active pairs", M, T, int(active_mask.sum()))
    return markers, active_mask, type_names


def initialize_beta_matrix(
    raw_counts: np.ndarray,
    markers: List[str],
    type_names: List[str],
    median_N: float,
    *,
    soft_scale: float = 0.1,
    inactive_val: float = 1e-3,
) -> np.ndarray:
    """Initialize per-type-per-marker beta matrix from raw counts.

    Strong types get full base rate (split if multiple strong),
    soft types get soft_scale * base rate, inactive get inactive_val.

    Args:
        raw_counts: (I, M) raw count matrix (column-max normalized).
        markers: List of M canonical marker names.
        type_names: List of T cell type names.
        median_N: Median cellularity across spots.
        soft_scale: Scale factor for soft assignments (default 0.1).
        inactive_val: Beta value for inactive pairs (default 1e-3).

    Returns:
        (M, T) beta initialization matrix.
    """
    M = len(markers)
    T = len(type_names)
    beta = np.full((M, T), inactive_val, dtype=np.float64)

    for m_idx, marker in enumerate(markers):
        assignments = MARKER_TYPE_TABLE.get(marker, [])
        col_median = max(float(np.median(raw_counts[:, m_idx])), 1e-6)
        lam_base = col_median / max(median_N, 1.0)

        strong_types = [(type_names.index(t), t) for t, s in assignments if s == "strong" and t in type_names]
        soft_types = [(type_names.index(t), t) for t, s in assignments if s == "soft" and t in type_names]

        total_weight = len(strong_types) + len(soft_types) * soft_scale
        if total_weight <= 0:
            continue

        for t_idx, _ in strong_types:
            beta[m_idx, t_idx] = lam_base / total_weight
        for t_idx, _ in soft_types:
            beta[m_idx, t_idx] = lam_base * soft_scale / total_weight

    return beta


def build_beta_prior_sigma(
    markers: List[str],
    type_names: List[str],
    *,
    sigma_exclusive: float = 5.0,
    sigma_shared: float = 2.0,
    sigma_inactive: float = 0.1,
    sigma_scale: float = 1.0,
) -> np.ndarray:
    """Build per-pair regularization sigma matrix.

    Exclusive strong markers get loose prior, shared/soft get moderate,
    inactive get tight. Multiplied by sigma_scale for calibration sweeps.

    Args:
        markers: List of M canonical marker names.
        type_names: List of T cell type names.
        sigma_exclusive: Prior sigma for exclusive strong markers.
        sigma_shared: Prior sigma for shared strong / soft markers.
        sigma_inactive: Prior sigma for inactive pairs.
        sigma_scale: Global multiplier for all sigma values.

    Returns:
        (M, T) sigma matrix.
    """
    M = len(markers)
    T = len(type_names)
    sigma = np.full((M, T), sigma_inactive * sigma_scale, dtype=np.float64)

    for m_idx, marker in enumerate(markers):
        assignments = MARKER_TYPE_TABLE.get(marker, [])
        strong_types = [t for t, s in assignments if s == "strong" and t in type_names]
        is_exclusive = len(strong_types) == 1

        for type_name, strength in assignments:
            if type_name not in type_names:
                continue
            t_idx = type_names.index(type_name)
            if strength == "strong" and is_exclusive:
                sigma[m_idx, t_idx] = sigma_exclusive * sigma_scale
            else:
                sigma[m_idx, t_idx] = sigma_shared * sigma_scale

    return sigma
