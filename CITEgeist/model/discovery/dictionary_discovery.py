"""Data-driven hierarchical marker dictionary: parents from Module 2, child states inside parents.

Parents are the profiles Module 2 finds among the interesting markers. Inside each
parent's positive spots, markers that are themselves bimodal and that rise and fall
together are grouped into child "states". The hierarchical marker share then splits
each parent's share among its children. No cell-type labels are used anywhere.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

from CITEgeist.model.deconvolution.marker_share import _major_markers, cosine_normalize, marker_share_proportions
from CITEgeist.model.discovery.marker_interest import identify_interesting_markers
from CITEgeist.model.discovery.nb_mixture import fit_nb_mixture, size_factors_from_counts
from CITEgeist.model.discovery.spatial_colocalization import analyze_marker_colocalization, discover_profiles_continuous

ISOTYPES = frozenset({"mouse_IgG2a", "mouse_IgG1k", "mouse_IgG2bk", "rat_IgG2a"})
_COLOC_PERMUTATIONS = 199  # profile discovery uses the Moran's I value, not its p-value


@dataclass
class DiscoveredDictionary:
    """parents: name -> markers. children: parent -> {child_name: markers}.

    Child names look like ``S1_a1+a2``; the full child label is ``"{parent}/{child_name}"``.
    """

    parents: dict[str, list[str]]
    children: dict[str, dict[str, list[str]]]
    diagnostics: dict = field(default_factory=dict)


def _base(name: str) -> str:
    return name[:-2] if name.endswith("-1") else name


def _is_isotype(name: str) -> bool:
    return _base(name) in ISOTYPES


def _find_child_states(
    X,
    var_names,
    parent_markers,
    other_parent_markers,
    parent_score,
    *,
    seed,
    min_subset_frac,
    max_subset_frac,
    state_corr,
    min_parent_spots,
    size_factors=None,
    extra_reject=None,
    report_other_parent=False,
):
    """Returns (children, diagnostics) for one parent. X is the raw count matrix.

    extra_reject(j, pos, sf) -> reason or None is an optional last check on a candidate that
    passed the NB-mixture and subset-fraction checks (used by the guided mixing guard).
    report_other_parent records other-parent markers in the diagnostics as rejected.
    """
    var_names = list(var_names)
    X = np.asarray(X, dtype=float)
    sf = size_factors_from_counts(X) if size_factors is None else np.asarray(size_factors, dtype=float)
    diag: dict = {"n_positive_spots": 0, "candidates": {}}
    pidx = [var_names.index(m) for m in parent_markers]
    # anchor: parent marker most correlated (Spearman, all spots) with the parent score
    if len(pidx) == 1:
        anchor = parent_markers[0]
    else:
        rs = rankdata(parent_score)
        rhos = []
        for j in pidx:
            c = np.corrcoef(rankdata(X[:, j] / sf), rs)[0, 1]
            rhos.append(-np.inf if np.isnan(c) else c)
        anchor = parent_markers[int(np.argmax(rhos))]  # argmax: first on ties
    diag["anchor"] = anchor
    fit = fit_nb_mixture(X[:, pidx].sum(axis=1), sf, seed=seed)
    n_pos = int(fit.positive_mask.sum())
    diag["n_positive_spots"] = n_pos
    diag["parent_score_n_components"] = fit.n_components
    if fit.n_components < 2:
        diag["no_children_reason"] = "parent counts have no positive population"
        return {}, diag
    if n_pos < min_parent_spots:
        diag["no_children_reason"] = "too few positive spots"
        return {}, diag
    pos = fit.positive_mask
    state_idx = []
    for j, name in enumerate(var_names):
        if name in other_parent_markers and report_other_parent:
            diag["candidates"][name] = {"positive_fraction": None, "rejected": "other_major_marker"}
            continue
        if name in other_parent_markers or name == anchor or _is_isotype(name):
            continue
        mfit = fit_nb_mixture(X[pos, j], sf[pos], seed=seed)
        if mfit.n_components < 2:
            diag["candidates"][name] = {
                "positive_fraction": None,
                "rejected": "single component in parent-positive spots",
            }
            continue
        frac = float(mfit.positive_mask.mean())
        if not (min_subset_frac <= frac <= max_subset_frac):
            diag["candidates"][name] = {"positive_fraction": frac, "rejected": "positive fraction outside range"}
            continue
        reason = extra_reject(j, pos, sf) if extra_reject is not None else None
        diag["candidates"][name] = {"positive_fraction": frac, "rejected": reason}
        if reason is None:
            state_idx.append(j)
    if len(state_idx) < 2:
        diag["no_children_reason"] = "fewer than 2 state markers"
        return {}, diag
    sub = X[np.ix_(pos, state_idx)] / sf[pos][:, None]
    rho = np.corrcoef(rankdata(sub, axis=0), rowvar=False)  # Spearman; nan for a constant channel
    dist = 1.0 - np.nan_to_num(rho, nan=0.0)  # constant channel: no correlation
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)
    labels = fcluster(
        linkage(squareform(np.clip(dist, 0.0, 2.0), checks=False), method="average"),
        t=1.0 - state_corr,
        criterion="distance",
    )
    clusters = [sorted(np.array(state_idx)[labels == c]) for c in np.unique(labels)]
    if len(clusters) < 2:
        diag["no_children_reason"] = "fewer than 2 state clusters"
        return {}, diag
    clusters.sort(key=lambda c: c[0])
    children = {}
    for k, idxs in enumerate(clusters, start=1):
        mk = [var_names[i] for i in idxs]
        children[f"S{k}_{'+'.join(mk)}"] = mk
    return children, diag


def find_child_states(
    X: np.ndarray,
    var_names: Sequence[str],
    parent_markers: list[str],
    other_parent_markers: set[str],
    parent_score: np.ndarray,
    *,
    seed: int,
    min_subset_frac: float,
    max_subset_frac: float,
    state_corr: float,
    min_parent_spots: int,
) -> dict[str, list[str]]:
    """Find child states of one parent from RAW counts X (size factors from X).

    parent_score (summed parent counts / size factor, all spots) is used only to
    pick the anchor marker of multi-marker parents. Parent-positive spots come
    from an NB mixture on the summed raw parent counts. Returns {child_name: markers};
    empty if the parent has no positive population, too few positive spots, or
    fewer than two state clusters.
    """
    children, _ = _find_child_states(
        X,
        var_names,
        parent_markers,
        other_parent_markers,
        parent_score,
        seed=seed,
        min_subset_frac=min_subset_frac,
        max_subset_frac=max_subset_frac,
        state_corr=state_corr,
        min_parent_spots=min_parent_spots,
    )
    return children


def discover_dictionary(
    X: np.ndarray,
    coords: np.ndarray,
    var_names: Sequence[str],
    *,
    exclude: Iterable[str] = (),
    seed: int = 1234,
    min_subset_frac: float = 0.1,
    max_subset_frac: float = 0.9,
    state_corr: float = 0.3,
    min_parent_spots: int = 50,
) -> DiscoveredDictionary:
    """Discover parents (Module 2 profiles) and child states without labels.

    Isotype channels and ``exclude`` are dropped before anything else, so the
    cosine normalization here is over the kept channels only. Pass the same kept
    channels to ``hierarchical_marker_share``.
    """
    var_names = list(var_names)
    X = np.asarray(X, dtype=float)
    drop = set(exclude)
    keep = [j for j, n in enumerate(var_names) if n not in drop and not _is_isotype(n)]
    names = [var_names[j] for j in keep]
    Xk = X[:, keep]
    diagnostics: dict = {"marker_n_components": {}, "per_parent": {}}

    interest = identify_interesting_markers(Xk, coords, names, seed=seed, verbose=False, auto_components=True)
    for m in interest.markers:
        diagnostics["marker_n_components"][m.name] = m.n_components
    interesting = [n for n in names if n in set(interest.interesting_markers)]

    profiles: list[list[str]] = []
    if len(interesting) >= 2:
        coloc = analyze_marker_colocalization(
            Xk,
            coords,
            names,
            markers_to_analyze=interesting,
            n_permutations=_COLOC_PERMUTATIONS,
            seed=seed,
            verbose=False,
        )
        found = discover_profiles_continuous(coloc, seed=seed, verbose=False)
        profiles = [list(p) for p in found.profiles]
        profiles += [[s] for s in found.singletons if not any(s in p for p in profiles)]
    else:
        profiles = [[m] for m in interesting]
    order = {n: i for i, n in enumerate(names)}
    profiles = [sorted(p, key=order.__getitem__) for p in profiles]
    profiles.sort(key=lambda p: order[p[0]])
    parents = {f"P{i}_{'+'.join(p)}": p for i, p in enumerate(profiles, start=1)}

    sf = size_factors_from_counts(Xk)
    children: dict[str, dict[str, list[str]]] = {}
    for pname, pm in parents.items():
        others = {m for q, mk in parents.items() if q != pname for m in mk}
        score = Xk[:, [names.index(m) for m in pm]].sum(axis=1) / sf  # summed counts / size factor
        ch, d = _find_child_states(
            Xk,
            names,
            pm,
            others,
            score,
            seed=seed,
            min_subset_frac=min_subset_frac,
            max_subset_frac=max_subset_frac,
            state_corr=state_corr,
            min_parent_spots=min_parent_spots,
            size_factors=sf,
        )
        diagnostics["per_parent"][pname] = d
        if ch:
            children[pname] = ch
    logging.info("Discovered %d parents, %d with child states", len(parents), len(children))
    return DiscoveredDictionary(parents=parents, children=children, diagnostics=diagnostics)


def discover_states_guided(
    X: np.ndarray,
    var_names: Sequence[str],
    hand_profile: Mapping[str, Mapping | list],
    parent_types: Sequence[str],
    *,
    exclude: Iterable[str] = (),
    seed: int = 1234,
    min_subset_frac: float = 0.1,
    max_subset_frac: float = 0.9,
    state_corr: float = 0.3,
    min_parent_spots: int = 50,
) -> DiscoveredDictionary:
    """Find child states only inside chosen hand-dictionary types.

    Parents are all hand types (Major markers). A state candidate must not be a Major
    marker of any other hand type, must pass the same NB-mixture and subset checks as
    ``discover_dictionary`` inside the parent's positive spots, and must not track another
    hand type's marker-share fraction (Spearman >= ``state_corr`` within parent-positive
    spots), which would mean other cell types in the same spots are posing as a state.
    Isotype channels and ``exclude`` are dropped first; shares and the guard use the kept channels.
    """
    var_names = list(var_names)
    X = np.asarray(X, dtype=float)
    drop = set(exclude)
    keep = [j for j, n in enumerate(var_names) if n not in drop and not _is_isotype(n)]
    names = [var_names[j] for j in keep]
    Xk = X[:, keep]
    parents = {t: list(_major_markers(e)) for t, e in hand_profile.items()}
    unknown = [t for t in parent_types if t not in parents]
    if unknown:
        raise ValueError(f"parent_types not in hand_profile: {unknown}")
    frac, _ = marker_share_proportions(Xk, names, parents)  # raises on a missing marker
    sf = size_factors_from_counts(Xk)
    diagnostics: dict = {"per_parent": {}}
    children: dict[str, dict[str, list[str]]] = {}
    for pname in parent_types:
        pm = parents[pname]
        others = {m for q, mk in parents.items() if q != pname for m in mk}
        other_types = [q for q in parents if q != pname]
        score = Xk[:, [names.index(m) for m in pm]].sum(axis=1) / sf

        def mixing_guard(j, pos, sf_, other_types=other_types):
            cand = rankdata(Xk[pos, j] / sf_[pos])
            best_t, best_r = None, -np.inf
            for q in other_types:
                r = np.corrcoef(cand, rankdata(frac[q].to_numpy()[pos]))[0, 1]
                r = 0.0 if np.isnan(r) else float(r)
                if r > best_r:
                    best_t, best_r = q, r
            if best_t is not None and best_r >= state_corr:
                return f"mixing_corr={best_t}:{best_r:.2f}"
            return None

        ch, d = _find_child_states(
            Xk,
            names,
            pm,
            others,
            score,
            seed=seed,
            min_subset_frac=min_subset_frac,
            max_subset_frac=max_subset_frac,
            state_corr=state_corr,
            min_parent_spots=min_parent_spots,
            size_factors=sf,
            extra_reject=mixing_guard,
            report_other_parent=True,
        )
        diagnostics["per_parent"][pname] = d
        if ch:
            children[pname] = ch
    return DiscoveredDictionary(parents=parents, children=children, diagnostics=diagnostics)


def hierarchical_marker_share(
    X: np.ndarray,
    var_names: Sequence[str],
    dictionary: DiscoveredDictionary,
    brightness: Mapping[str, float] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parent fractions, plus child fractions that split each parent's share.

    Returns (parent fractions, child fractions with columns "{parent}/{child}").
    Child fraction = parent fraction * child score / sum of the parent's child
    scores; an even split is used where those scores are all zero. The count of
    even-split spots per parent is stored in ``dictionary.diagnostics["even_split_spots"]``.
    """
    var_names = list(var_names)
    parent_frac, n_zero = marker_share_proportions(X, var_names, dictionary.parents, brightness=brightness)
    Xc = cosine_normalize(X)
    cols: dict[str, np.ndarray] = {}
    even: dict[str, int] = {}
    for pname, kids in dictionary.children.items():
        scores = []
        for markers in kids.values():
            c = Xc[:, [var_names.index(m) for m in markers]]
            if brightness is not None:
                c = c / np.array([brightness.get(m, 1.0) for m in markers], dtype=float)
            scores.append(c.mean(axis=1))
        S = np.column_stack(scores)
        tot = S.sum(axis=1)
        zero = tot == 0
        share = np.divide(S, tot[:, None], out=np.zeros_like(S), where=tot[:, None] > 0)
        share[zero] = 1.0 / S.shape[1]
        even[pname] = int(zero.sum())
        pf = parent_frac[pname].to_numpy()
        for j, cname in enumerate(kids):
            cols[f"{pname}/{cname}"] = pf * share[:, j]
    dictionary.diagnostics["even_split_spots"] = even
    dictionary.diagnostics["n_all_zero_spots"] = n_zero
    return parent_frac, pd.DataFrame(cols, index=parent_frac.index)
