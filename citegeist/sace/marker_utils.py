"""Antibody marker-name utilities."""
from __future__ import annotations

from typing import Collection, Optional

# Markers whose real name naturally ends in "-1" and must NOT be de-suffixed.
NATURAL_DASH_ONE_MARKERS = frozenset({"PD-1"})


def strip_antibody_suffix(
    name: str,
    active_markers: Optional[Collection[str]] = None,
    preserve: frozenset[str] = NATURAL_DASH_ONE_MARKERS,
) -> str:
    """Strip the SpaceRanger ``-1`` antibody suffix, preserving natural ``-1`` names.

    Default mode (``active_markers=None``) matches the pipeline's historical strip-all
    rule: strip a trailing ``-1`` (``EPCAM-1`` -> ``EPCAM``) EXCEPT for names in ``preserve``
    (the real checkpoint marker ``PD-1`` stays ``PD-1``).

    Caller-vocabulary mode (``active_markers`` given): strip ``-1`` only when the de-suffixed
    base is in ``active_markers`` (and not preserved); otherwise leave ``name`` unchanged.

    Args:
        name: antibody/marker name, possibly carrying the ``-1`` SpaceRanger suffix.
        active_markers: optional caller vocabulary of base marker names; when provided,
            only names whose base is in it are stripped.
        preserve: names to leave untouched even though they end in ``-1``.

    Returns:
        The de-suffixed name, or ``name`` unchanged if it does not end in ``-1``, is preserved,
        or (in caller-vocabulary mode) its base is not in ``active_markers``.
    """
    if not name.endswith("-1") or name in preserve:
        return name
    base = name[:-2]
    if active_markers is not None and base not in active_markers:
        return name
    return base
