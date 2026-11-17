"""Provenance stamps for producer-written results files.

A stamp records who produced a file (`_produced_by`), at which code state
(`_code_commit`, `_code_dirty`), and a sha256 of the payload
(`_payload_sha256`). `verify` recomputes the hash, so a hand-edited value is
detected. JSON files carry the stamp inline; CSVs get a `<csv>.provenance.json`
sidecar holding a sha256 over the CSV bytes plus the three code-state fields.

Tamper-evident, not tamper-proof: there is no secret, so someone who recomputes
the hash can forge a stamp. It catches accidental and hand edits.
"""

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

STAMP_KEYS = ("_produced_by", "_code_commit", "_code_dirty", "_payload_sha256")
# Only these top-level keys are left out of the hash; every other key (even `_`-prefixed) is hashed.
META_KEYS = frozenset({"_payload_sha256", "_source", "_sources", "_source_md5", "_generated_at"})
CODE_KEYS = ("_produced_by", "_code_commit", "_code_dirty")
REPO_ROOT = Path(__file__).resolve().parent.parent
# Written at the repo root of the public export (no .git in a downloaded zip); never present in the dev checkout.
STAMP_FILE = "CODE_VERSION.json"


class CodeProvenanceError(RuntimeError):
    """The code version can be determined neither from git nor from CODE_VERSION.json."""


def _git_toplevel(repo_root):
    """Resolved git work-tree root containing repo_root, or None (no git binary / not in a work tree)."""
    try:
        r = subprocess.run(
            ["git", "--no-optional-locks", "rev-parse", "--show-toplevel"],
            cwd=repo_root,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, NotADirectoryError):
        return None
    return Path(r.stdout.strip()).resolve() if r.returncode == 0 and r.stdout.strip() else None


def in_git_checkout(repo_root=REPO_ROOT):
    """True only when repo_root IS a git work-tree root (an export unpacked inside an unrelated repo is not)."""
    return _git_toplevel(repo_root) == Path(repo_root).resolve()


def code_commit_info(repo_root=REPO_ROOT):
    """Return (commit_sha, source); source is "git" or "stamp". Never returns an empty sha.

    git HEAD when repo_root is a git checkout; otherwise `source_commit` from <repo_root>/CODE_VERSION.json
    (written by the public sync). Raises CodeProvenanceError if neither is usable.
    """
    repo_root = Path(repo_root)
    if in_git_checkout(repo_root):
        r = subprocess.run(
            ["git", "--no-optional-locks", "rev-parse", "HEAD"], cwd=repo_root, capture_output=True, text=True
        )
        sha = r.stdout.strip()
        if r.returncode != 0 or not sha:
            raise CodeProvenanceError(f"git rev-parse HEAD failed in {repo_root}: {r.stderr.strip() or 'empty output'}")
        return sha, "git"
    stamp = repo_root / STAMP_FILE
    if not stamp.is_file():
        raise CodeProvenanceError(
            f"no git checkout and no {STAMP_FILE} stamp in {repo_root}; cannot record code provenance"
        )
    try:
        sha = json.loads(stamp.read_text()).get("source_commit")
    except (ValueError, AttributeError) as e:
        raise CodeProvenanceError(f"malformed {stamp}: {e}") from e
    if not isinstance(sha, str) or not sha.strip():
        raise CodeProvenanceError(f"malformed {stamp}: missing non-empty string 'source_commit'")
    return sha.strip(), "stamp"


def code_commit(repo_root=REPO_ROOT):
    """Code version (git HEAD, else the CODE_VERSION.json source_commit). See code_commit_info."""
    return code_commit_info(repo_root)[0]


def _coerce(obj):
    """json `default=` hook: numpy scalars/arrays -> plain Python."""
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _git(*args):
    return subprocess.run(
        ["git", "--no-optional-locks", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout


def _canonical(payload):
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), default=_coerce)


def _hash_payload(payload):
    body = {k: v for k, v in payload.items() if k not in META_KEYS}
    return hashlib.sha256(_canonical(body).encode("utf-8")).hexdigest()


def _hash_csv(csv_path, fields):
    h = hashlib.sha256(Path(csv_path).read_bytes())
    h.update(_canonical({k: fields.get(k) for k in CODE_KEYS}).encode("utf-8"))
    return h.hexdigest()


def _producer_dirty(produced_by):
    """True unless `produced_by` is an existing, git-tracked file identical to HEAD.

    `-uno` (used for the repo-wide check) hides untracked files, so a never-committed
    producer would otherwise stamp clean. Unresolvable producers count as dirty.
    """
    path = Path(str(produced_by))
    path = path if path.is_absolute() else REPO_ROOT / path
    if not path.is_file():
        return True
    try:
        rel = str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return True
    if not _git("ls-files", "--", rel).strip():  # untracked (or ignored)
        return True
    return bool(_git("status", "--porcelain", "--", rel).strip())


def _code_fields(produced_by):
    commit, source = code_commit_info(REPO_ROOT)
    # Without git a clean tree cannot be verified: record "unknown" (None), never a faked False.
    dirty = (bool(_git("status", "--porcelain", "-uno").strip()) or _producer_dirty(produced_by)) if source == "git" else None
    return {
        "_produced_by": produced_by,
        "_code_commit": commit,
        "_code_commit_source": source,
        "_code_dirty": dirty,
    }


def stamp_payload(payload, produced_by):
    """Return payload (numpy coerced) plus the four stamp fields."""
    if not isinstance(payload, dict):
        raise TypeError(f"payload must be a dict, got {type(payload).__name__}")
    out = json.loads(_canonical(payload))  # coerces numpy; keeps NaN
    out.update(_code_fields(produced_by))
    out["_payload_sha256"] = _hash_payload(out)
    return out


def write_stamped_json(path, payload, produced_by):
    stamped = stamp_payload(payload, produced_by)
    Path(path).write_text(json.dumps(stamped, indent=2, sort_keys=True) + "\n")
    return stamped


def _sidecar_path(csv_path):
    csv_path = Path(csv_path)
    return csv_path.with_name(csv_path.name + ".provenance.json")


def write_csv_sidecar(csv_path, produced_by):
    sidecar = _code_fields(produced_by)
    sidecar["_payload_sha256"] = _hash_csv(csv_path, sidecar)
    _sidecar_path(csv_path).write_text(json.dumps(sidecar, indent=2, sort_keys=True) + "\n")
    return sidecar


def _check_fields(stamp, actual_hash):
    if not any(k in stamp for k in STAMP_KEYS):
        return False, "unstamped: no provenance fields"
    missing = [k for k in STAMP_KEYS if k not in stamp]
    if missing:
        return False, f"missing stamp fields: {missing}"
    if stamp["_payload_sha256"] != actual_hash:
        return False, "hash mismatch: payload was modified after stamping"
    return True, "ok"


def verify(path):
    """Return (ok, message) for a stamped JSON file or a CSV with sidecar."""
    path = Path(path)
    if path.suffix.lower() == ".csv":
        sidecar = _sidecar_path(path)
        if not sidecar.exists():
            return False, f"unstamped: no sidecar {sidecar.name}"
        try:
            stamp = json.loads(sidecar.read_text())
        except ValueError as e:
            return False, f"corrupt sidecar: {e}"
        if not isinstance(stamp, dict):
            return False, "sidecar is not a dict"
        return _check_fields(stamp, _hash_csv(path, stamp))
    try:
        data = json.loads(path.read_text())
    except ValueError as e:
        return False, f"corrupt JSON: {e}"
    if not isinstance(data, dict):
        return False, "top level is not a dict"
    return _check_fields(data, _hash_payload(data))
