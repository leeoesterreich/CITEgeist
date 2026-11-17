"""GT-free antibody-brightness estimation: fixed RNA gene panels and RNA score helper.

Panels are fixed per dataset (spec sections 4.3 / 4.3a) and are never edited after a candidate is scored.
The estimation panels share no gene with ``EVALUATION_PANEL`` (used only by the breast safety check).
"""
import logging
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)

MIN_PANEL_GENES = 3  # a type with fewer present genes is not estimable

_SIM_AND_BREAST = {
    "B": ["BANK1", "CD22", "CD37", "TNFRSF13B"],
    "CD4 T": ["CD40LG", "ICOS", "CD28", "LTB"],
    "CD8 T": ["GZMK", "GZMA", "CCL5", "PRF1"],
    "Endothelial": ["CLDN5", "EGFL7", "RAMP2", "ESAM"],
    "Fibroblasts": ["COL3A1", "COL6A1", "FBLN1", "PDGFRA"],
}
_EPI = ["KRT7", "ELF3", "CLDN4", "MUC1"]

ESTIMATION_PANELS: Dict[str, Dict[str, List[str]]] = {
    "sim": {
        "Epithelial": list(_EPI),
        "B": list(_SIM_AND_BREAST["B"]),
        "Plasmablasts": ["MZB1", "JCHAIN", "SDC1", "XBP1"],
        "T": ["TRAC", "TRBC2", "CD3G", "LCK"],
        "Macrophages": ["CD14", "MS4A7", "CSF1R", "AIF1"],
        "Endothelial": list(_SIM_AND_BREAST["Endothelial"]),
        "Fibroblasts": list(_SIM_AND_BREAST["Fibroblasts"]),
        "PVL": ["RGS5", "ACTA2", "PDGFRB", "NOTCH3"],
    },
    "breast": {
        "Cancer": list(_EPI),
        "B": list(_SIM_AND_BREAST["B"]),
        "CD4 T": list(_SIM_AND_BREAST["CD4 T"]),
        "CD8 T": list(_SIM_AND_BREAST["CD8 T"]),
        "Macrophages": ["MS4A7", "CSF1R", "AIF1", "MS4A6A"],
        "Endothelial": list(_SIM_AND_BREAST["Endothelial"]),
        "Fibroblasts": list(_SIM_AND_BREAST["Fibroblasts"]),
    },
    "xenium7": {
        "Epithelial": ["KRT7", "MUC1", "CDH1", "EHF"],
        "B": ["BANK1", "TNFRSF13B", "SPIB", "TCL1A"],
        "CD4 T": ["CD28", "CTLA4", "FOXP3", "TIGIT"],
        "CD8 T": ["GZMK", "GZMA", "CCL5", "PRF1"],
        "Macrophages": ["CD14", "AIF1", "MS4A6A", "MPEG1"],
        "Endothelial": ["EGFL7", "RAMP2", "ERG", "ECSCR"],
        "Fibroblasts": ["FBLN1", "PDGFRA", "DPT", "FBN1"],
    },
    "xenium6": {
        "Epithelial": ["KRT7", "MUC1", "CDH1", "EHF"],
        "B": ["BANK1", "TNFRSF13B", "SPIB", "TCL1A"],
        "T": ["TRAC", "CD3G", "CD3D", "CD2"],
        "Macrophages": ["CD14", "AIF1", "MS4A6A", "MPEG1"],
        "Endothelial": ["EGFL7", "RAMP2", "ERG", "ECSCR"],
        "Fibroblasts": ["FBLN1", "PDGFRA", "DPT", "FBN1"],
    },
}

EVALUATION_PANEL: Dict[str, List[str]] = {
    "Cancer": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "CD4 T": ["CD4", "CD3E", "IL7R"],
    "CD8 T": ["CD8A", "CD8B", "CD3E"],
    "Macrophages": ["CD68", "CD163", "C1QA", "C1QB"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM"],
}

# Panel key -> profile (dictionary) type name(s) it scores, per dataset. Panel keys are not profile
# type names; a wrong entry silently mis-scores a type, so this is a fixed, tested constant.
# Sim: one Epithelial panel scores both Cancer and Normal Epithelial (identical panel -> one column,
# so downstream *_adj regressions see it once, not as two collinear columns).
PANEL_TO_PROFILE_TYPES: Dict[str, Dict[str, List[str]]] = {
    "sim": {
        "Epithelial": ["Cancer Epithelial", "Normal Epithelial"],
        "B": ["B-cells"],
        "Plasmablasts": ["Plasmablasts"],
        "T": ["T-cells"],
        "Macrophages": ["Myeloid"],
        "Endothelial": ["Endothelial"],
        "Fibroblasts": ["CAFs"],
        "PVL": ["PVL"],
    },
    "breast": {
        "Cancer": ["Cancer_Cells"],
        "B": ["B_Cells"],
        "CD4 T": ["CD4_T_Cells"],
        "CD8 T": ["CD8_T_Cells"],
        "Macrophages": ["Macrophages"],
        "Endothelial": ["Endothelial"],
        "Fibroblasts": ["Fibroblasts"],
    },
    "xenium7": {
        "Epithelial": ["Epithelial"],
        "B": ["B cells"],
        "CD4 T": ["CD4+ T cells"],
        "CD8 T": ["CD8+ T cells"],
        "Macrophages": ["Macrophages"],
        "Endothelial": ["Endothelial"],
        "Fibroblasts": ["Fibroblasts"],
    },
    "xenium6": {
        "Epithelial": ["Epithelial"],
        "B": ["B cells"],
        "T": ["T cells"],
        "Macrophages": ["Macrophages"],
        "Endothelial": ["Endothelial"],
        "Fibroblasts": ["Fibroblasts"],
    },
}


def profile_type_to_panel(panel_set: str) -> Dict[str, str]:
    """Invert PANEL_TO_PROFILE_TYPES: profile type name -> panel key (score column)."""
    return {t: k for k, ts in PANEL_TO_PROFILE_TYPES[panel_set].items() for t in ts}


def rna_scores(
    gex: ad.AnnData, panel: Mapping[str, Sequence[str]], min_genes: int = MIN_PANEL_GENES
) -> Tuple[pd.DataFrame, dict]:
    """Per-spot RNA score of each panel type: scanpy score_genes on normalize_total(1e4)+log1p raw counts.

    gene_pool is all genes (per call). Absent genes are dropped and reported; a type with fewer than
    ``min_genes`` (default MIN_PANEL_GENES) present genes is not estimable and gets an all-NaN column.
    The breast safety check passes ``min_genes=1`` (its stored baseline scores CD8 T from 2 genes).
    The input is not mutated.
    Returns (spots x types DataFrame, {type: {"used", "dropped", "estimable"}}).
    """
    a = gex.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    pool = list(a.var_names)
    present = set(pool)
    out = {}
    report = {}
    for t, genes in panel.items():
        used = [g for g in genes if g in present]
        dropped = [g for g in genes if g not in present]
        ok = len(used) >= min_genes
        report[t] = {"used": used, "dropped": dropped, "estimable": ok}
        if ok:
            sc.tl.score_genes(a, gene_list=used, gene_pool=pool, score_name="_s", random_state=0)
            out[t] = a.obs["_s"].to_numpy(dtype=float)
        else:
            out[t] = np.full(a.n_obs, np.nan)
    return pd.DataFrame(out, index=gex.obs_names), report


def combined_proxy(scores: pd.DataFrame, types: Sequence[str]) -> np.ndarray:
    """Sum of the listed types' RNA scores; non-estimable (all-NaN) types are skipped.

    Returns an all-NaN vector if none of the listed types is estimable.
    """
    cols = [t for t in types if not scores[t].isna().all()]
    if not cols:
        return np.full(len(scores), np.nan)
    return scores[cols].sum(axis=1).to_numpy(dtype=float)


def effective_brightness(
    Xc: np.ndarray,
    var_names: Sequence[str],
    profile: Mapping,
    prop: pd.DataFrame,
    min_slope: float = 0.001,
) -> Tuple[Dict[str, float], Dict[str, dict]]:
    """Effective brightness b*_m per Major marker (estimand of spec S4.1).

    b*_m is the OLS slope (with intercept, ``np.polyfit(ab, col, 1)[0]`` as in
    the oracle-brightness check of the development benchmark)
    of marker m's cosine-normalised spot value on the SUMMED proportion of every type listing m as Major.
    Spots with NaN in the listing types' ``prop`` columns are excluded. A slope <= ``min_slope`` is
    replaced by the median of the other markers' slopes that exceed it (``median_fallback`` True). A marker whose
    summed listing-type proportion has no variation across the used spots (e.g. all types absent) has an unidentifiable
    slope: ``raw_slope`` is NaN, ``degenerate_abundance`` True, and it takes the same neutral median.
    Xc rows must align with prop rows. Returns (b*, diag); diag has raw_slope, median_fallback,
    n_spots_used, types. Not geometric-mean normalised (see ``geomean_normalize``).
    """
    from CITEgeist.model.deconvolution.marker_share import _major_markers

    vn = list(var_names)
    Xc = np.asarray(Xc, dtype=float)
    slope: Dict[str, float] = {}
    diag: Dict[str, dict] = {}
    for m in sorted({m for v in profile.values() for m in _major_markers(v)}):
        types = [t for t, v in profile.items() if m in _major_markers(v)]
        sub = prop[types]
        keep = sub.notna().all(axis=1).to_numpy()
        ab = sub.to_numpy(dtype=float)[keep].sum(axis=1)
        col = Xc[keep, vn.index(m)]
        diag[m] = {"types": types, "n_spots_used": int(keep.sum())}
        if ab.size < 2 or np.ptp(ab) == 0:
            # Marker's listing types have no abundance variation (e.g. a type absent from the whole sample):
            # the slope is unidentifiable (np.polyfit scales by |ab|, so all-zero ab gives NaN and an SVD failure).
            # NaN is never > min_slope, so the marker takes the neutral median below.
            slope[m] = float("nan")
            diag[m]["degenerate_abundance"] = True
        else:
            slope[m] = float(np.polyfit(ab, col, 1)[0])
        diag[m]["raw_slope"] = slope[m]
    bstar: Dict[str, float] = {}
    for m, s in slope.items():
        if s > min_slope:
            bstar[m] = s
            diag[m]["median_fallback"] = False
    for m in slope:
        if m not in bstar:
            others = [s for k, s in slope.items() if k != m and s > min_slope]
            if not others:
                raise RuntimeError(f"no positive slopes to take a median from (marker {m})")
            bstar[m] = float(np.median(others))
            diag[m]["median_fallback"] = True
    return bstar, diag


def geomean_normalize(b: Mapping[str, float]) -> Dict[str, float]:
    """Rescale so the geometric mean over markers is 1. Raises ValueError on non-positive/non-finite values."""
    vals = np.array(list(b.values()), dtype=float)
    if vals.size == 0 or not np.all(np.isfinite(vals)) or np.any(vals <= 0):
        raise ValueError(f"brightness must be finite and positive, got {dict(b)}")
    gm = float(np.exp(np.mean(np.log(vals))))
    return {m: float(v) / gm for m, v in b.items()}


def neutral_fill(b: Mapping[str, Optional[float]], flagged: set) -> Dict[str, float]:
    """Fill flagged markers with the median of the OTHER markers' finite positive values (never 1.0)."""
    good = [float(v) for m, v in b.items() if m not in flagged and v is not None and np.isfinite(v) and v > 0]
    if not good:
        raise ValueError("no finite positive brightness values to take a neutral median from")
    fill = float(np.median(good))
    return {m: (fill if m in flagged else float(v)) for m, v in b.items()}


MIN_GATE_SPOTS = 20
GATE_QUANTILE = 0.9


def _is_estimable(entry) -> bool:
    return bool(entry["estimable"]) if isinstance(entry, Mapping) else bool(entry)


def _n_genes(entry) -> Optional[int]:
    return len(entry["used"]) if isinstance(entry, Mapping) and "used" in entry else None


def _residualize(proxy: np.ndarray, others: List[np.ndarray]) -> np.ndarray:
    """OLS residual of proxy on the other types' scores (with intercept), over spots finite in all."""
    out = np.full_like(proxy, np.nan)
    if not others:
        return proxy
    D = np.column_stack([np.ones(len(proxy))] + others)
    ok = np.isfinite(proxy) & np.all(np.isfinite(D), axis=1)
    coef, *_ = np.linalg.lstsq(D[ok], proxy[ok], rcond=None)
    out[ok] = proxy[ok] - D[ok] @ coef
    return out


def _other_columns(scores: pd.DataFrame, estimable: Mapping, own: Sequence[str]) -> List[np.ndarray]:
    """Estimable, non-constant score columns of the OTHER types, identical-valued columns kept once.

    A column identical to one of the marker's own columns (same panel scoring several types, e.g. Cancer
    and Normal Epithelial) is the marker's own signal and is excluded too.
    """
    seen = [scores[t].to_numpy(dtype=float) for t in own if t in scores]
    cols: List[np.ndarray] = []
    for t in scores.columns:
        if t in own or t not in estimable or not _is_estimable(estimable[t]):
            continue
        v = scores[t].to_numpy(dtype=float)
        if not np.all(np.isfinite(v)) or np.ptp(v) == 0:
            continue
        if any(np.array_equal(v, s) for s in seen) or any(np.array_equal(v, s) for s in cols):
            continue
        cols.append(v)
    return cols


def est_rna_gate(
    Xc: np.ndarray,
    var_names: Sequence[str],
    profile: Mapping,
    scores: pd.DataFrame,
    estimable: Mapping,
    *,
    adjusted: bool,
) -> Tuple[Dict[str, float], Dict[str, dict]]:
    """RNA-gated brightness: b_m = median of m's cosine-normalised value over the top-decile RNA-proxy spots.

    ``scores`` columns are PROFILE type names (a panel scoring several types is repeated under each name);
    ``estimable`` maps profile type -> bool or the ``rna_scores`` report entry ({"used", "estimable"}).
    The proxy of marker m is the summed score of the estimable types listing m as Major. Gate = proxy
    >= the 90th percentile. ``adjusted=True`` first replaces the proxy by its OLS residual on the scores
    of all other types. A marker whose proxy has no estimable type, is constant, has < 20 usable gate
    spots, or whose gate median is non-positive is flagged ``fallback`` and filled with the median of the
    other markers (never 1.0). The reported vector (fallbacks included) has geometric mean 1. If EVERY
    marker falls back there is no median to take: all b are equal (1.0 after normalisation, i.e. no
    correction, same as plain) and every diag entry has ``all_fallback`` True (no exception).
    Conventions: ``scores`` columns are profile type names (a panel scoring several types is repeated
    under each); ``estimable`` is Mapping[type, bool] (or the rna_scores report entry); n_spots_used
    counts gate spots with a finite marker value.
    diag per marker: used_brightness, method, proxy_types, n_genes_used, n_spots_used, fallback, reason,
    selected_spot_purity_proxy (mean standardised proxy in gate spots), selected_idx (gate row indices), all_fallback.
    """
    from CITEgeist.model.deconvolution.marker_share import _major_markers

    vn = list(var_names)
    Xc = np.asarray(Xc, dtype=float)
    method = "rna_gate_adj" if adjusted else "rna_gate"
    raw: Dict[str, Optional[float]] = {}
    diag: Dict[str, dict] = {}
    flagged = set()
    for m in sorted({m for v in profile.values() for m in _major_markers(v)}):
        listing = [t for t, v in profile.items() if m in _major_markers(v)]
        types = [t for t in listing if t in scores and t in estimable and _is_estimable(estimable[t])]
        d = {
            "used_brightness": None,
            "method": method,
            "proxy_types": types,
            "n_genes_used": None if not types else sum(_n_genes(estimable[t]) or 0 for t in types) or None,
            "n_spots_used": 0,
            "fallback": True,
            "reason": None,
            "selected_spot_purity_proxy": float("nan"),
            "selected_idx": np.array([], dtype=int),
        }
        diag[m] = d
        b = None
        if not types:
            d["reason"] = "no_estimable_type"
        else:
            proxy = combined_proxy(scores, types)
            if adjusted:
                proxy = _residualize(proxy, _other_columns(scores, estimable, listing))
            fin = np.isfinite(proxy)
            sd = float(np.std(proxy[fin])) if fin.any() else 0.0
            if not fin.any() or np.ptp(proxy[fin]) == 0 or not np.isfinite(sd) or sd == 0.0:
                d["reason"] = "proxy_constant"
            else:
                z = (proxy - np.mean(proxy[fin])) / sd
                thr = np.nanpercentile(proxy[fin], 100 * GATE_QUANTILE)
                sel = np.flatnonzero(fin & (proxy >= thr))
                col = Xc[sel, vn.index(m)]
                good = np.isfinite(col)
                d["selected_idx"] = sel
                d["n_spots_used"] = int(good.sum())
                d["selected_spot_purity_proxy"] = float(np.mean(z[sel]))
                if good.sum() < MIN_GATE_SPOTS:
                    d["reason"] = "n_spots<20"
                else:
                    med = float(np.median(col[good]))
                    if np.isfinite(med) and med > 0:
                        b = med
                        d["fallback"] = False
                    else:
                        d["reason"] = "nonpositive_b"
        raw[m] = b
        if b is None:
            flagged.add(m)
    all_fallback = len(flagged) == len(raw)
    if all_fallback:  # nothing to take a median from: uniform b, i.e. no correction (equals plain)
        out = {m: 1.0 for m in raw}
    else:
        out = geomean_normalize(neutral_fill(raw, flagged))
    for m in out:
        diag[m]["used_brightness"] = out[m]
        diag[m]["all_fallback"] = all_fallback
    return out, diag


def est_rna_slope(
    Xc: np.ndarray,
    var_names: Sequence[str],
    profile: Mapping,
    scores: pd.DataFrame,
    estimable: Mapping,
    *,
    adjusted: bool,
) -> Tuple[Dict[str, float], Dict[str, dict]]:
    """RNA-slope brightness: b_m = slope of m's cosine-normalised value on its RNA proxy.

    Same conventions as ``est_rna_gate`` (``scores`` columns are profile type names, ``estimable`` is
    Mapping[type, bool] or the ``rna_scores`` report entry). The proxy of marker m is the summed score
    of the estimable types listing m as Major. Unadjusted: Theil-Sen slope (``scipy.stats.theilslopes``).
    Adjusted: coefficient of the own proxy in one ``HuberRegressor`` fit of the value on the own proxy
    plus every other estimable type's score; constant/non-finite columns, duplicate columns and columns
    identical to the marker's own types' are dropped first. A non-positive/non-finite slope, a constant
    proxy, < 20 usable spots or no estimable type flags the marker ``fallback`` (reason
    ``non_positive_slope`` / ``proxy_constant`` / ``n_spots<20`` / ``no_estimable_type``); it is filled
    with the median of the other markers, then the whole vector is geometric-mean normalised. If every
    marker falls back: uniform b = 1.0, ``all_fallback`` True, no exception.
    diag per marker: as ``est_rna_gate`` (selected_idx = spots used in the fit;
    selected_spot_purity_proxy is None).
    """
    from scipy.stats import theilslopes
    from sklearn.linear_model import HuberRegressor

    from CITEgeist.model.deconvolution.marker_share import _major_markers

    vn = list(var_names)
    Xc = np.asarray(Xc, dtype=float)
    method = "rna_slope_adj" if adjusted else "rna_slope"
    raw: Dict[str, Optional[float]] = {}
    diag: Dict[str, dict] = {}
    flagged = set()
    for m in sorted({m for v in profile.values() for m in _major_markers(v)}):
        listing = [t for t, v in profile.items() if m in _major_markers(v)]
        types = [t for t in listing if t in scores and t in estimable and _is_estimable(estimable[t])]
        d = {
            "used_brightness": None,
            "method": method,
            "proxy_types": types,
            "n_genes_used": None if not types else sum(_n_genes(estimable[t]) or 0 for t in types) or None,
            "n_spots_used": 0,
            "fallback": True,
            "reason": None,
            "selected_spot_purity_proxy": None,
            "selected_idx": np.array([], dtype=int),
        }
        diag[m] = d
        b = None
        if not types:
            d["reason"] = "no_estimable_type"
        else:
            proxy = combined_proxy(scores, types)
            others = _other_columns(scores, estimable, listing) if adjusted else []
            y = Xc[:, vn.index(m)]
            ok = np.isfinite(proxy) & np.isfinite(y)
            for o in others:
                ok &= np.isfinite(o)
            sel = np.flatnonzero(ok)
            d["selected_idx"] = sel
            d["n_spots_used"] = int(sel.size)
            if sel.size < MIN_GATE_SPOTS:
                d["reason"] = "n_spots<20"
            elif np.ptp(proxy[sel]) == 0:
                d["reason"] = "proxy_constant"
            else:
                if adjusted:
                    D = np.column_stack([proxy[sel]] + [o[sel] for o in others])
                    mu, sd = D.mean(axis=0), D.std(axis=0)
                    sd[sd == 0] = 1.0  # a column constant over the used spots carries no slope
                    ys = y[sel].std() or 1.0
                    fit = HuberRegressor(alpha=0.0, max_iter=10000, tol=1e-12)
                    fit.fit((D - mu) / sd, y[sel] / ys)
                    s = float(fit.coef_[0] / sd[0] * ys)
                else:
                    s = float(theilslopes(y[sel], proxy[sel])[0])
                if np.isfinite(s) and s > 0:
                    b = s
                    d["fallback"] = False
                else:
                    d["reason"] = "non_positive_slope"
        raw[m] = b
        if b is None:
            flagged.add(m)
    all_fallback = len(flagged) == len(raw)
    if all_fallback:
        out = {m: 1.0 for m in raw}
    else:
        out = geomean_normalize(neutral_fill(raw, flagged))
    for m in out:
        diag[m]["used_brightness"] = out[m]
        diag[m]["all_fallback"] = all_fallback
    return out, diag


METHODS = ("plain", "gmm", "rna_gate", "rna_gate_adj", "rna_slope", "rna_slope_adj")


@dataclass
class BrightnessResult:
    """Brightness estimate: applied_b feeds ``marker_share_proportions(brightness=...)`` (for gmm: the raw GMM values, so
    predictions stay identical to the current estimator); reported_b is applied_b rescaled to geometric
    mean 1 and is what diagnostics, the harness and E_xen use. diag maps marker -> per-marker dict."""

    applied_b: Dict[str, float]
    reported_b: Dict[str, float]
    diag: Dict[str, dict]


def _markers(profile: Mapping) -> List[str]:
    from CITEgeist.model.deconvolution.marker_share import _major_markers

    return sorted({m for v in profile.values() for m in _major_markers(v)})


def _plain(profile: Mapping) -> BrightnessResult:
    ms = _markers(profile)
    b = {m: 1.0 for m in ms}
    diag = {m: {"used_brightness": 1.0, "method": "plain", "fallback": False, "reason": None} for m in ms}
    return BrightnessResult(dict(b), dict(b), diag)


def estimate_from_scores(
    method: str,
    X: np.ndarray,
    var_names: Sequence[str],
    profile: Mapping,
    scores: pd.DataFrame,
    estimable: Mapping[str, bool],
    *,
    seed: int = 1234,
) -> BrightnessResult:
    """The ONLY implementation of the five methods. ``scores`` columns are PROFILE type names.

    X is the raw antibody matrix (cosine-normalised here); rows of ``scores`` must align with X rows.
    """
    from CITEgeist.model.deconvolution.marker_share import cosine_normalize, estimate_marker_brightness_with_diagnostics

    if method not in METHODS:
        raise ValueError(f"unknown brightness method {method!r}; expected one of {METHODS}")
    if method == "plain":
        return _plain(profile)
    if method == "gmm":
        applied, diag = estimate_marker_brightness_with_diagnostics(
            X, var_names, _markers(profile), seed=seed, max_components=3
        )
        diag = {m: {**d, "method": "gmm"} for m, d in diag.items()}
        return BrightnessResult(applied, geomean_normalize(applied), diag)
    if len(scores) != np.asarray(X).shape[0]:
        raise ValueError(f"scores has {len(scores)} rows but X has {np.asarray(X).shape[0]}")
    Xc = cosine_normalize(X)
    fn = {
        "rna_gate": est_rna_gate,
        "rna_gate_adj": est_rna_gate,
        "rna_slope": est_rna_slope,
        "rna_slope_adj": est_rna_slope,
    }[method]
    b, diag = fn(Xc, var_names, profile, scores, estimable, adjusted=method.endswith("_adj"))
    return BrightnessResult(dict(b), dict(b), diag)


def estimate_brightness(
    method: str,
    X: np.ndarray,
    var_names: Sequence[str],
    profile: Mapping,
    gex: Optional[ad.AnnData] = None,
    rna_panel: Optional[Mapping[str, Sequence[str]]] = None,
    *,
    panel_set: Optional[str] = None,
    spot_names: Optional[Sequence[str]] = None,
    allow_plain_fallback: bool = False,
    seed: int = 1234,
) -> BrightnessResult:
    """GEX-level wrapper: validate, compute RNA scores, delegate to ``estimate_from_scores``.

    ``spot_names`` are the antibody matrix's spot names (X row order); RNA methods require them and
    ``gex.obs_names`` must equal them exactly (same order, unique), else ValueError (no silent re-order).
    If ``rna_panel`` is None it defaults to ``ESTIMATION_PANELS[panel_set]`` and ``panel_set`` is required.
    With ``panel_set`` set, panel-key score columns are mapped to profile type names via
    ``PANEL_TO_PROFILE_TYPES``; without it the panel keys must already be profile type names.
    """
    if method not in METHODS:
        raise ValueError(f"unknown brightness method {method!r}; expected one of {METHODS}")
    if method in ("plain", "gmm"):
        return estimate_from_scores(method, X, var_names, profile, pd.DataFrame(), {}, seed=seed)
    if gex is None:
        if not allow_plain_fallback:
            raise ValueError(f"method {method!r} needs gex; pass allow_plain_fallback=True to use plain")
        logger.warning("no gex given for %s: falling back to plain brightness", method)
        return _plain(profile)
    if rna_panel is None:
        if panel_set is None:
            raise ValueError("panel_set is required for RNA methods when rna_panel is not given")
        if panel_set not in ESTIMATION_PANELS:
            raise ValueError(f"unknown panel_set {panel_set!r}")
        rna_panel = ESTIMATION_PANELS[panel_set]
    if spot_names is None:
        raise ValueError("spot_names (antibody spot names, X row order) are required for RNA methods")
    spot_names = list(spot_names)
    gnames = list(gex.obs_names)
    if len(set(gnames)) != len(gnames) or len(set(spot_names)) != len(spot_names):
        raise ValueError("gex obs_names and spot_names must be unique (duplicates found)")
    if len(spot_names) != np.asarray(X).shape[0]:
        raise ValueError("spot_names length does not match X rows")
    if gnames != spot_names:
        raise ValueError("gex obs_names do not match the antibody spot names (same names, same order required)")
    panel_scores, report = rna_scores(gex, rna_panel)
    if panel_set is not None:
        mapping = PANEL_TO_PROFILE_TYPES[panel_set]
        cols, est = {}, {}
        for key in panel_scores.columns:
            for t in mapping.get(key, []):
                cols[t] = panel_scores[key].to_numpy()
                est[t] = report[key]
        scores = pd.DataFrame(cols, index=panel_scores.index)
    else:
        scores = panel_scores
        est = dict(report)
    return estimate_from_scores(method, X, var_names, profile, scores, est, seed=seed)
