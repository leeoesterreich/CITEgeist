"""Select the breast hand dictionary by protocol (single source: CITEgeist/examples/breast_dictionary.py).

``protocol["inputs"].get("dictionary_version", "v40")``: "v40" -> BREAST_PROFILE_V40 (the frozen v8 protocol has no such key,
so it keeps v40 behaviour), "v41" -> BREAST_PROFILE. Anything else raises.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Optional

_REL = Path("CITEgeist") / "examples" / "breast_dictionary.py"
# Walk up from this file (also works from a copied src/ tree if it sits under the repo); fall back to CWD.
_MOD_PATH = next((q / _REL for q in (*Path(__file__).resolve().parents, Path.cwd()) if (q / _REL).is_file()), None)
if _MOD_PATH is None:
    raise FileNotFoundError(f"cannot locate {_REL}")


def _load():
    spec = importlib.util.spec_from_file_location("_breast_dictionary_single_source", _MOD_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_BD = _load()
BREAST_PROFILE = _BD.BREAST_PROFILE
BREAST_PROFILE_V40 = _BD.BREAST_PROFILE_V40


def profile_for_protocol(protocol: Optional[dict]) -> dict:
    version = ((protocol or {}).get("inputs") or {}).get("dictionary_version", "v40")
    if version == "v40":
        return BREAST_PROFILE_V40
    if version == "v41":
        return BREAST_PROFILE
    raise ValueError(f"unknown inputs.dictionary_version {version!r} (expected 'v40' or 'v41')")


def types_for_protocol(protocol: Optional[dict]) -> list:
    return list(profile_for_protocol(protocol))
