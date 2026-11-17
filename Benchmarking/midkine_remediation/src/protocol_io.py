"""A12 protocol plumbing: hash-checked loader, typed statuses, output-namespace and v3 gates.

Gates covered here: G9 (stamps carry the protocol sha256; producers abort on hash
mismatch), G11 (refuse non-empty output target, atomic write, v3 hashes unchanged,
R00 has no ``obsp['spatial_distance']``), and the typed-status contract of
``gates_fail_closed`` (a missing/NaN value without a status fails the run).
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

VALID_STATUS = ("ok", "not_estimable", "failed")
_FROZEN_RE = re.compile(rb'("frozen_sha256"\s*:\s*)(null|"[0-9a-f]{64}")')
TO_FILL_TOKEN = "TO_FILL_PRE_FREEZE"


class ProtocolHashError(RuntimeError):
    """The protocol file does not hash to its recorded ``frozen_sha256`` (or is not frozen)."""


class OutputTargetError(RuntimeError):
    """Refused output target (non-empty, wrong namespace, or the v3 directory)."""


class GateFailure(RuntimeError):
    """A fail-closed gate tripped."""


# --------------------------------------------------------------------------- protocol hash
def protocol_hash(raw: bytes) -> str:
    """sha256 of the file bytes with the ``frozen_sha256`` value replaced by ``null``."""
    if len(_FROZEN_RE.findall(raw)) != 1:
        raise ProtocolHashError("protocol must contain exactly one frozen_sha256 field")
    nulled = _FROZEN_RE.sub(rb"\1null", raw, count=1)
    return hashlib.sha256(nulled).hexdigest()


def load_protocol(path, allow_draft: bool = False) -> dict:
    """Load the protocol and verify its hash. Aborts (raises) on any mismatch.

    ``status == FROZEN`` requires ``frozen_sha256`` to equal the recomputed hash.
    A DRAFT is refused unless ``allow_draft`` (pre-freeze tooling only).
    """
    raw = Path(path).read_bytes()
    proto = json.loads(raw)
    status = proto.get("status")
    if status == "FROZEN":
        recorded = proto.get("frozen_sha256")
        actual = protocol_hash(raw)
        if recorded != actual:
            raise ProtocolHashError(f"protocol hash mismatch: recorded={recorded} actual={actual}")
        proto["_protocol_sha256"] = actual
        proto["_protocol_sha12"] = actual[:12]
        return proto
    if not allow_draft:
        raise ProtocolHashError(f"protocol status is {status!r}, not FROZEN; producers must not run")
    actual = protocol_hash(raw)
    proto["_protocol_sha256"] = actual
    proto["_protocol_sha12"] = actual[:12]
    return proto


ACTIVE_PROTOCOL_FILE = Path(__file__).resolve().parents[1] / "ACTIVE_PROTOCOL"


def active_sha12(path=None) -> str:
    """The sha12 of the CURRENT protocol (one line in ``ACTIVE_PROTOCOL``); consumers that mean "current" read this.

    Paths that must stay pinned to an older protocol (v40 manuscript artefacts) keep their literal sha12 instead.
    """
    text = Path(path or ACTIVE_PROTOCOL_FILE).read_text().strip()
    if not re.fullmatch(r"[0-9a-f]{12}", text):
        raise ProtocolHashError(f"ACTIVE_PROTOCOL must hold one 12-hex sha12, got {text!r}")
    return text


def manuscript_sha12() -> str:
    """sha12 of the protocol run the manuscript data/figures are built from (``MANUSCRIPT_PROTOCOL``). Equal to
    ``active_sha12()`` once the manuscript data is rebuilt from the current run; v40 (0344d54f9d3d) until then."""
    return active_sha12(ACTIVE_PROTOCOL_FILE.with_name("MANUSCRIPT_PROTOCOL"))


def with_protocol_hash(payload: dict, proto: dict) -> dict:
    """G9: copy of ``payload`` carrying the protocol sha256 (a producer writes this inside its provenance stamp)."""
    if proto.get("status") != "FROZEN" or "_protocol_sha256" not in proto:
        raise ProtocolHashError("only a verified FROZEN protocol may be stamped into outputs")
    return {**payload, "protocol_sha256": proto["_protocol_sha256"], "protocol_sha12": proto["_protocol_sha12"]}


def freeze_text(raw: bytes) -> bytes:
    """Return ``raw`` with status FROZEN and the frozen_sha256 filled (used by the freeze step and tests)."""
    raw = re.sub(rb'("status"\s*:\s*)"DRAFT"', rb'\1"FROZEN"', raw, count=1)
    digest = protocol_hash(raw)
    return _FROZEN_RE.sub(rb'\1"' + digest.encode() + b'"', raw, count=1)


def find_to_fill(obj, path: str = "") -> list:
    """Paths of every string VALUE containing the TO_FILL token (dict keys are not scanned)."""
    hits = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            hits += find_to_fill(v, f"{path}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            hits += find_to_fill(v, f"{path}[{i}]")
    elif isinstance(obj, str) and TO_FILL_TOKEN in obj:
        hits.append(path)
    return hits


META_PATH_PREFIXES = ("/freeze_rule", "/pre_freeze_checklist")   # these two fields NAME the token while describing the rule


def unfilled(proto: dict) -> list:
    """Unfilled placeholder values, excluding the two meta fields that merely mention the token."""
    return [h for h in find_to_fill(proto) if not h.startswith(META_PATH_PREFIXES)]


# --------------------------------------------------------------------------- typed statuses
@dataclass(frozen=True)
class Typed:
    status: str
    reason: str = ""
    value: Optional[float] = None

    def __post_init__(self):
        if self.status not in VALID_STATUS:
            raise ValueError(f"status {self.status!r} not in {VALID_STATUS}")
        if self.status != "ok" and not self.reason:
            raise ValueError("a non-ok status must carry a reason")


def validate_typed(d: dict) -> None:
    """Fail closed: every outcome is {status, reason, value?}; ``ok`` needs a finite value."""
    if "status" not in d:
        raise GateFailure("outcome has no status")
    t = Typed(d["status"], d.get("reason", ""), d.get("value"))
    if t.status == "ok":
        v = t.value
        if v is None or (isinstance(v, float) and not math.isfinite(v)):
            raise GateFailure("status ok but value missing/NaN")


# --------------------------------------------------------------------------- G11
def check_output_target(target, protocol_sha12: str, run_id: str, allow_finalized: bool = False) -> Path:
    """Refuse a non-empty target, the v3 directory, and any path outside ``<sha12>/<run_id>``.

    ``allow_finalized=True`` bypasses ONLY the non-empty-target refusal, and only when the target is a directory holding ``specimen.json``
    (a finalized specimen that a ``--resume`` run will skip). Namespace and v3 checks always apply.
    """
    t = Path(target)
    parts = t.parts
    if "morphology_assignment_v3" in parts:
        raise OutputTargetError("never write into morphology_assignment_v3")
    if protocol_sha12 not in parts or run_id not in parts:
        raise OutputTargetError(f"target must be namespaced <{protocol_sha12}>/<{run_id}>: {t}")
    if allow_finalized and t.is_dir() and (t / "specimen.json").is_file():
        return t
    if t.exists() and (not t.is_dir() or any(t.iterdir())):
        raise OutputTargetError(f"target exists and is not an empty directory: {t}")
    return t


def atomic_write_json(path, payload) -> None:
    """Write via temp file in the same directory, then ``os.replace``."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=p.parent, prefix=".tmp_", suffix=".json")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
        os.replace(tmp, p)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def v3_hashes_ok(manifest_path, base_dir) -> Optional[list]:
    """Return a list of mismatching paths (empty = unchanged). Raises on an empty/absent manifest."""
    lines = [ln for ln in Path(manifest_path).read_text().splitlines() if ln.strip()]
    if not lines:
        raise GateFailure("v3 manifest is empty")
    bad = []
    for ln in lines:
        digest, rel = ln.split(None, 1)
        f = Path(base_dir) / rel.strip()
        h = hashlib.sha256()
        if not f.is_file():
            bad.append(str(f))
            continue
        with open(f, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != digest:
            bad.append(str(f))
    return bad


def assert_no_spatial_distance(adata_like, run_id: str) -> None:
    """R00 (and every run but R10) must not carry ``obsp['spatial_distance']``."""
    has = "spatial_distance" in getattr(adata_like, "obsp", {})
    if run_id == "R10":
        if not has:
            raise GateFailure("R10 must set obsp['spatial_distance']")
    elif has:
        raise GateFailure(f"{run_id} must not carry obsp['spatial_distance']")
