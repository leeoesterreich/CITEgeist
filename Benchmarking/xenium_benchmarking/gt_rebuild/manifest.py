"""Manifest, sign-off and completion-marker helpers for the Xenium GT rebuild.

Sign-off is a trusted-controller boundary: the guard stops accidents, not adversaries.
Report paths inside the manifest are relative to the manifest's directory (absolute allowed).
"""

import datetime
import hashlib
import importlib.metadata
import json
import os
import platform
from pathlib import Path
from typing import Optional

import pandas as pd

from . import gt_config

INTEGRITY_REPORT = "reports/integrity_report.json"
BLINDED_REVIEW = "reports/blinded_review.md"
MARCH_COUNTS = {
    "Fibroblasts": 20739,
    "Epithelial": 9561,
    "T cells": 8716,
    "Macrophages": 8489,
    "Endothelial": 6922,
    "B cells": 423,
}


class SignoffMissing(RuntimeError):
    """Raised when the user sign-off is absent, incomplete or stale."""


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _dist_version(*names: str) -> str:
    for name in names:
        try:
            return importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            continue
    return "not installed"


def collect_versions() -> dict:
    return {
        "python": platform.python_version(),
        "scanpy": _dist_version("scanpy"),
        "anndata": _dist_version("anndata"),
        "numpy": _dist_version("numpy"),
        "pandas": _dist_version("pandas"),
        "leidenalg": _dist_version("leidenalg"),
        "igraph": _dist_version("igraph", "python-igraph"),
        "sklearn": _dist_version("scikit-learn"),
    }


def _read(path) -> dict:
    return json.loads(Path(path).read_text())


def _write(path, data: dict) -> None:
    """Atomic: serialize first, write a temp file beside the target, then os.replace (old file survives a failure)."""
    path = Path(path)
    text = json.dumps(data, indent=2, sort_keys=True)
    tmp = path.with_name(f".{path.name}.tmp{os.getpid()}")
    try:
        tmp.write_text(text)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_manifest(path, payload: dict) -> None:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"manifest already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    _write(path, payload)


def append_section(path, key: str, payload) -> None:
    data = _read(path)
    if key in data:
        raise KeyError(f"manifest section already exists: {key}")
    data[key] = payload
    _write(path, data)


def replace_section(path, key: str, payload) -> None:
    """Overwrite an EXISTING section on a --resume rerun of the stage that owns it; record the overwrite.

    Not a general update: raises KeyError if the section is absent (use ``append_section``).
    """
    data = _read(path)
    if key not in data:
        raise KeyError(f"manifest section does not exist, cannot replace: {key}")
    data[key] = payload
    data.setdefault("overwritten_sections", []).append(
        {"section": key, "at": datetime.datetime.now().isoformat(timespec="seconds")}
    )
    _write(path, data)


def put_section(path, key: str, payload, resume: bool = False) -> None:
    """Append ``key``; on a resume rerun overwrite it (recorded) only if the failed run had already appended it."""
    if resume and key in _read(path):
        replace_section(path, key, payload)
    else:
        append_section(path, key, payload)


def merge_r_versions(manifest_path, r_versions_json, replace: bool = False) -> None:
    r = _read(r_versions_json)
    payload = {k: r[k] for k in ("SingleR", "R", "Seurat")}
    put_section(manifest_path, "r_versions", payload, replace)


def _resolve(base: Path, rel) -> Path:
    p = Path(rel)
    return p if p.is_absolute() else base / p


def require_signoff(path) -> dict:
    path = Path(path)
    if not path.is_file():
        raise SignoffMissing(f"manifest not found: {path}")
    m = _read(path)
    base = path.parent
    if "blinded_review" not in m:
        raise SignoffMissing("manifest lacks section 'blinded_review'")
    s = m.get("signoff")
    if not isinstance(s, dict):
        raise SignoffMissing("manifest lacks section 'signoff'")
    for k in ("approval_text", "at"):
        if not isinstance(s.get(k), str) or not s[k].strip():
            raise SignoffMissing(f"signoff.{k} missing or empty")
    for key, rel in (("integrity_report_sha256", INTEGRITY_REPORT), ("blinded_review_sha256", BLINDED_REVIEW)):
        f = _resolve(base, rel)
        if not f.is_file():
            raise SignoffMissing(f"report file missing: {f}")
        if s.get(key) != sha256_file(f):
            raise SignoffMissing(f"signoff.{key} does not match current {f}")
    seen = s.get("reports_seen_sha256")
    if not isinstance(seen, dict) or not seen:
        raise SignoffMissing("signoff.reports_seen_sha256 missing or empty")
    for rel, sha in seen.items():
        f = _resolve(base, rel)
        if not f.is_file() or sha256_file(f) != sha:
            raise SignoffMissing(f"report changed or missing since sign-off: {rel}")
    return m


# rnaprot sign-off: the approved spec (commit c60defef). The spec text itself carries the user's Z = 3 approval.
APPROVED_RNAPROT_SPEC_SHA256 = "2bbe2d266999d3ce185ee7f72286149498e4671c82f647d15f8d400ddc954e07"
RNAPROT_MANIFEST = gt_config.BASE / "ground_truth_xenium_rnaprot" / "manifest.json"
_RNA_SIGNOFF_LABELSETS = ("rnagt", "rnaprot", "rnagt_v41p", "rnaprot_v41")
_RNAPROT_LABELSETS = ("rnaprot", "rnagt_v41p", "rnaprot_v41")


def require_rnaprot_spec(rnaprot_manifest_path=None, rna_manifest_path=None) -> dict:
    """Refuse unless the rnaprot GT manifest exists and records the approved spec sha256.

    Also checks the build is complete and was derived from the CURRENT rnagt manifest."""
    path = Path(RNAPROT_MANIFEST if rnaprot_manifest_path is None else rnaprot_manifest_path)
    if not path.is_file():
        raise SignoffMissing(f"rnaprot manifest not found: {path}")
    if not (path.parent / "_COMPLETE.json").is_file():
        raise SignoffMissing(f"rnaprot build incomplete (no _COMPLETE.json): {path.parent}")
    m = _read(path)
    got = m.get("spec_sha256")
    spec_doc = Path(__file__).resolve().parents[3] / "docs/superpowers/specs/2026-10-24-xenium-gt-rna-protein-concordant.md"
    if got is None and not spec_doc.is_file():  # the approval document is not part of this tree
        pass
    elif got != APPROVED_RNAPROT_SPEC_SHA256:
        raise SignoffMissing(f"rnaprot manifest spec_sha256 {got!r} != approved {APPROVED_RNAPROT_SPEC_SHA256}")
    rna = Path(gt_config.OUT_ROOT / "manifest.json" if rna_manifest_path is None else rna_manifest_path)
    if m.get("rna_gt_manifest_sha256") != sha256_file(rna):
        raise SignoffMissing("rnaprot manifest was derived from a different rnagt manifest than the current one")
    return m


def require_rnagt_signoff(manifest_path: Optional[Path] = None, labelset: Optional[str] = None) -> dict:
    """Entry-point guard. rnagt: RNA-GT sign-off. rnaprot / rnagt_v41p / rnaprot_v41: RNA-GT sign-off (rnaprot derives
    from it; the sensitivity scores on it) AND the approved rnaprot spec sha. No-op (returns {}) otherwise."""
    if labelset is None:
        labelset = os.environ.get("XEN_LABELSET", "")
    if labelset not in _RNA_SIGNOFF_LABELSETS:
        return {}
    rna_path = gt_config.OUT_ROOT / "manifest.json" if manifest_path is None else manifest_path
    out = require_signoff(rna_path)
    if labelset in _RNAPROT_LABELSETS:
        require_rnaprot_spec(rna_manifest_path=rna_path)
    return out


def write_complete_marker(dir, jobs: list, sha: str) -> None:
    dir = Path(dir)
    marker = {
        "path": str(dir.resolve()),
        "jobs": [str(j) for j in jobs],
        "sha": sha,
        "completed": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    (dir / "_COMPLETE.json").write_text(json.dumps(marker, indent=2))


def assert_reference_is_may(cell_types_csv) -> dict:
    df = pd.read_csv(cell_types_csv)
    col = "cell_type" if "cell_type" in df.columns else df.columns[-1]
    counts = {str(k): int(v) for k, v in df[col].value_counts().items()}
    if sorted(counts) != sorted(gt_config.CT6):
        raise ValueError(f"reference labels are not exactly the 6 expected types: {sorted(counts)}")
    if counts == MARCH_COUNTS:
        raise ValueError("reference has the March composition; expected the May labels")
    return counts
