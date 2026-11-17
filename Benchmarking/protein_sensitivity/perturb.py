"""Pure protein-count perturbations for the v40 protein-input stress test.

All functions act on raw spot x marker count matrices BEFORE
preprocess_antibody(), so CLR (.X) and layers["raw_counts"] stay consistent.
"QP-active" = Major + Soft markers, matching qp_solver.map_antibodies_to_profiles_v2.
"""
from __future__ import annotations

import copy

import numpy as np


def as_integer_counts(X) -> np.ndarray:
    A = np.asarray(X)
    if not np.issubdtype(A.dtype, np.integer):
        if not np.allclose(A, np.round(A)):
            raise ValueError("protein counts must be integer-valued for binomial thinning")
        A = np.round(A)
    if (A < 0).any():
        raise ValueError("protein counts must be non-negative")
    return A.astype(np.int64)


def thin_counts(X, keep_frac: float, rng: np.random.Generator) -> np.ndarray:
    A = as_integer_counts(X)
    if keep_frac >= 1.0:
        return A.copy()
    return rng.binomial(A, keep_frac).astype(np.int64)


def dropout_channels(X, var_names, candidates, frac, rng):
    A = as_integer_counts(X).copy()
    names = list(var_names)
    pool = [m for m in candidates if m in names]
    k = int(round(frac * len(pool)))
    dropped = [str(m) for m in rng.choice(pool, size=k, replace=False)] if k else []
    for m in dropped:
        A[:, names.index(m)] = 0
    return A, dropped


def estimate_background(X) -> np.ndarray:
    A = as_integer_counts(X).astype(float)
    med = np.median(A, axis=0)
    return np.array([A[A[:, j] <= med[j], j].mean() for j in range(A.shape[1])])


def add_background_noise(X, multiple, background, rng):
    A = as_integer_counts(X)
    lam = np.broadcast_to(multiple * np.asarray(background, float), A.shape)
    return (A + rng.poisson(lam)).astype(np.int64)


def _soft_names(spec):
    return [s[0] for s in spec.get("Soft", [])]


def qp_active_markers(profile_dict) -> list[str]:
    out: list[str] = []
    for spec in profile_dict.values():
        for m in spec.get("Major", []):
            if m not in out:
                out.append(m)
    for spec in profile_dict.values():
        for m in _soft_names(spec):
            if m not in out:
                out.append(m)
    return out


def active_profile(profile_dict):
    out = {}
    for t, s in profile_dict.items():
        d = {"Major": list(s.get("Major", []))}
        if s.get("Soft"):
            d["Soft"] = [tuple(x) for x in s["Soft"]]
        out[t] = d
    return out


def _one_per_type(full):
    counts: dict[str, int] = {}
    for s in full.values():
        for m in s["Major"] + _soft_names(s):
            counts[m] = counts.get(m, 0) + 1
    one = {}
    for t, s in full.items():
        cands = s["Major"] or _soft_names(s)
        specific = [m for m in cands if counts[m] == 1]
        pick = (specific or cands)[0]
        one[t] = {"Major": [pick]} if pick in s["Major"] else {"Major": [], "Soft": [x for x in s["Soft"] if x[0] == pick]}
    return one


def ablation_panels(profile_dict):
    full = active_profile(profile_dict)
    return [("full", copy.deepcopy(full)), ("one_per_type", _one_per_type(full)), ("zero", {})]


def type_removal_series(profile_dict):
    one = _one_per_type(active_profile(profile_dict))
    series = [("one_per_type", copy.deepcopy(one))]
    cur, removed = copy.deepcopy(one), []
    for t in sorted(one):  # fixed alphabetical order, never re-ranked
        cur.pop(t)
        removed.append(t)
        series.append(("minus_" + "+".join(removed) if cur else "zero", copy.deepcopy(cur)))
    return series


def leave_one_type_out(profile_dict):
    full = active_profile(profile_dict)
    return [(f"loto_{t}", {u: s for u, s in full.items() if u != t}) for t in full]
