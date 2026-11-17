"""Explicit run variant for the v41p protein-sensitivity rerun (marker-share profile + brightness refit).

`None` means the legacy path (ACHIEVABLE profile, no brightness); nothing in the legacy code path reads this module.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np

PROFILES = ("achievable", "ms_v41")


@dataclass(frozen=True)
class RunVariant:
    profile: str  # "achievable" | "ms_v41"
    brightness: bool


def from_args(profile: str, brightness) -> RunVariant | None:
    """CLI args -> variant; ("achievable", 0) is the legacy path and returns None."""
    if profile not in PROFILES:
        raise ValueError(f"profile must be one of {PROFILES}, got {profile!r}")
    brightness = bool(int(brightness))
    if profile == "achievable" and not brightness:
        return None
    return RunVariant(profile=profile, brightness=brightness)


def profile_for(n_types: int, variant: RunVariant | None) -> dict:
    """Legacy -> baseline.PROFILE[n]; ms_v41 -> MS_V41_{n}_CELL_PROFILE_DICT."""
    from protein_sensitivity import baseline as B

    if variant is None or variant.profile == "achievable":
        return B.PROFILE[n_types]
    import benchmark_constants as BC

    return {6: BC.MS_V41_6_CELL_PROFILE_DICT, 7: BC.MS_V41_7_CELL_PROFILE_DICT}[n_types]


def meta(variant: RunVariant | None, n_types: int) -> dict:
    prof = profile_for(n_types, variant)
    return {
        "profile": "achievable" if variant is None else variant.profile,
        "brightness": False if variant is None else variant.brightness,
        "profile_sha256": hashlib.sha256(json.dumps(prof, sort_keys=True).encode()).hexdigest(),
        "n_types": n_types,
    }


def perturbed_sha256(adata_cite) -> str:
    X = adata_cite.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    h = hashlib.sha256()
    h.update(json.dumps(list(map(str, adata_cite.var_names))).encode())
    h.update(json.dumps(list(map(str, adata_cite.obs_names))).encode())
    h.update(np.ascontiguousarray(X, dtype=np.float64).tobytes())
    return h.hexdigest()
