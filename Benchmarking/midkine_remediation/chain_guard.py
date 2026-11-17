"""Optional v42 chain guard: what the shipped stage modules import instead of the private ``chain_preflight``.

``chain_preflight.py`` (v42 run ids, stamped-input checks, lineage stamps) is private lineage tooling and is NOT part of the public export
(``Benchmarking/public_exclude.txt``). The stage modules (B3-B12) ship, so they must import and run without it. This module re-exports the
four names they use:

  * private ``chain_preflight`` present (the private repo, every v2 run): the real functions, unchanged;
  * absent (a public copy): ``guard_entry`` and ``stage_stamp_begin`` do nothing, ``is_v42`` is False and ``legacy_or_required`` returns the
    value or the legacy default. A public copy has no v42 run ids, so it is a plain legacy run.

Only a missing ``chain_preflight`` is tolerated, and only when ``V2_NAMESPACE_ROOT`` is NOT set (every v2 chain run sets it: there a missing file raises). Any other import error (a broken dependency of the real module) is raised, so a private run
can never fall back to "no checks" by accident.
"""
from __future__ import annotations

import os

try:
    import chain_preflight as _real
except ModuleNotFoundError as _e:
    if _e.name != "chain_preflight":
        raise
    if os.environ.get("V2_NAMESPACE_ROOT"):  # a v2 chain run: the checks must exist; never fall back to "no checks"
        raise RuntimeError(
            "V2_NAMESPACE_ROOT is set (a v2 chain run) but chain_preflight.py is missing from this tree: refusing to run without the v42 guard") from _e
    _real = None

AVAILABLE = _real is not None


def guard_entry(stage, protocol_path, **kw):
    return _real.guard_entry(stage, protocol_path, **kw) if _real else None


def stage_stamp_begin(stage, protocol_path, artifact, producer, run_id, **kw):
    return _real.stage_stamp_begin(stage, protocol_path, artifact, producer, run_id, **kw) if _real else None


def is_v42(sha12, runs=None):
    return _real.is_v42(sha12, runs) if _real else False


def legacy_or_required(value, protocol_path, what, legacy):
    if _real:
        return _real.legacy_or_required(value, protocol_path, what, legacy)
    return value if value is not None else legacy
