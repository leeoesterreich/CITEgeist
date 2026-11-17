"""The one shared stamp writer for v2 (v42) producers.

A producer calls ``write_stamp(artifact, ...)`` once per artifact it hands downstream. The writer computes everything
itself at write time (a producer never types a hash): it writes the sidecar ``<artifact>.stamp.json`` that
``scripts/check_v42_lineage.py`` (``validate_stamp`` + ``stamp_closure_problems``) validates and re-derives.

What it records beyond the checker's core fields: ``producer_path`` + ``producer_blob_sha`` (the git blob of the producer script
at ``code_commit``), ``slurm_job_id``, ``host``, ``started_at``, ``ended_at`` (and an optional ``stage_role`` note and
``external_inputs_sha256`` for inputs outside the repo, e.g. the raw Visium h5, which the checker cannot pin by repo-relative path).

What it refuses (raises StampError, writes nothing):
  * a producer file, or any tracked file of the work tree, that differs from HEAD (stamp the code that actually ran);
  * a producer that is not tracked at HEAD;
  * not running inside a SLURM job (SLURM_JOB_ID unset): heavy producers are sbatch jobs, and the id is the recorded-job pointer;
  * an input path that is absolute/outside the repo, names a symlink component, or names a v41 / stand-in artifact
    (``check_v42_lineage.FORBIDDEN_TOKENS``, resolved path included);
  * an input or the artifact itself whose sha256 is in the v41 denylist, unless its path is a named reuse-allowlist entry.

Tamper-evident, not tamper-proof: someone who edits this writer can forge a stamp; the checker's job is consistency.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import socket
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_v42_lineage as L  # noqa: E402

REPO = Path(__file__).resolve().parents[1]

# The frozen v2 Phase 1 settings every stamp asserts (the checker enforces the same values).
FROZEN = {"phase1_method": "rna_refined", "c": 0.1, "n_passes": 3, "anchor_thr": 0}


NAMESPACE_ENV = "V2_NAMESPACE_ROOT"


class StampError(RuntimeError):
    pass


def utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _git(root: Path, *args) -> str:
    r = subprocess.run(["git", "--no-optional-locks", "-C", str(root), *args], capture_output=True, text=True)
    if r.returncode != 0:
        raise StampError(f"git {' '.join(args)} failed in {root}: {r.stderr.strip()[:200]}")
    return r.stdout.strip()


def _slurm_job_id() -> str:
    jid = os.environ.get("SLURM_JOB_ID")
    if not jid:
        raise StampError("SLURM_JOB_ID is not set: producers that stamp must run as sbatch jobs (the job id is the provenance pointer)")
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    if os.environ.get("SLURM_ARRAY_JOB_ID") and task:
        return f"{os.environ['SLURM_ARRAY_JOB_ID']}_{task}"
    return jid


def _rel(root: Path, p) -> str:
    """Repo-relative LEXICAL path (a symlink is not resolved here: the checker's path_problems looks at it)."""
    p = Path(p)
    p = p if p.is_absolute() else root / p
    try:
        return str(Path(os.path.normpath(p)).relative_to(os.path.normpath(root)))
    except ValueError:
        raise StampError(f"{p} is outside the repo {root}; list it in external_inputs, not inputs")


def _common_git_dir(root: Path) -> Path:
    return (Path(root) / _git(root, "rev-parse", "--git-common-dir")).resolve()


def is_linked_worktree(root: Path) -> bool:
    return (Path(root) / _git(root, "rev-parse", "--git-dir")).resolve() != _common_git_dir(root)


def namespace_root(code_root=REPO) -> Path:
    """The tree the stamped paths are relative to (what ``check_v42_lineage`` resolves ``input_sha256`` keys against): the MAIN checkout.

    Code may run from a detached worktree (``code_root``) while outputs live under the main checkout, so the namespace root is
    REQUIRED, not guessed: ``$V2_NAMESPACE_ROOT`` (the wrapper exports it). Unset is accepted only when ``code_root`` is itself a
    main (non-linked) checkout. It must be a directory of the same git repository as ``code_root``."""
    code_root = Path(code_root)
    env = os.environ.get(NAMESPACE_ENV)
    if not env:
        if is_linked_worktree(code_root):
            raise StampError(
                f"{NAMESPACE_ENV} is not set and the code at {code_root} is a linked git worktree: the stamped paths are relative to the "
                f"main checkout, which cannot be guessed from a worktree. Export {NAMESPACE_ENV}=<main checkout> (submit-v2-breast-chain.sh does)"
            )
        return code_root
    ns = Path(env)
    if not ns.is_dir():
        raise StampError(f"{NAMESPACE_ENV}={env} is not a directory")
    try:
        same = _common_git_dir(ns) == _common_git_dir(code_root)
    except StampError as e:
        raise StampError(f"{NAMESPACE_ENV}={env} is not a git work tree: {e}")
    if not same:
        raise StampError(f"{NAMESPACE_ENV}={env} belongs to a different git repository than the code at {code_root}")
    if is_linked_worktree(ns):
        raise StampError(f"{NAMESPACE_ENV}={env} is a linked worktree; it must be the main checkout the checker runs in")
    return ns


def resolve_roots(repo_root=None, git_root=None) -> tuple:
    """(namespace_root, git_root). Explicit ``repo_root`` alone means one tree plays both roles (a plain checkout / tests)."""
    if repo_root is None:
        git = Path(git_root) if git_root is not None else REPO
        return namespace_root(git), git
    return Path(repo_root), Path(git_root if git_root is not None else repo_root)


def preflight(producer, out_paths=(), *, repo_root=None, git_root=None, require_clean_tree: bool = True) -> tuple:
    """Everything ``write_stamp`` will demand, checked at JOB START (before any expensive compute): SLURM job, namespace root,
    clean code tree, tracked+clean producer, and every output path inside the namespace root with no path problem. Returns the roots."""
    ns, git = resolve_roots(repo_root, git_root)
    _slurm_job_id()
    _code_state(git, _rel(git, producer), require_clean_tree)
    problems = []
    for p in out_paths:
        try:
            rel = _rel(ns, p)
        except StampError as e:
            problems.append(str(e))
            continue
        problems += [f"output {rel}: {x}" for x in L.path_problems(ns, rel, check_tokens=True)]
    if problems:
        raise StampError("; ".join(problems))
    return ns, git


def _code_state(root: Path, producer_rel: str, require_clean_tree: bool = True) -> tuple:
    head = _git(root, "rev-parse", "HEAD")
    try:
        blob = _git(root, "rev-parse", "--verify", f"HEAD:{producer_rel}")
        if _git(root, "cat-file", "-t", blob) != "blob":
            raise StampError("not a blob")
    except StampError:
        raise StampError(f"producer {producer_rel} is not tracked at HEAD ({head[:12]}): commit it before stamping")
    here = _git(root, "hash-object", "--", producer_rel)
    if here != blob:
        raise StampError(f"producer {producer_rel} differs from HEAD (dirty): commit it before stamping")
    if require_clean_tree and _git(root, "status", "--porcelain", "--untracked-files=no"):
        raise StampError("work tree has uncommitted changes to tracked files; stamp from a clean (detached) worktree")
    return head, blob


def write_stamp(
    artifact,
    *,
    run_id: str,
    protocol_sha12: str,
    inputs,
    producer,
    started_at: str,
    repo_root=None,
    git_root=None,
    stage_role: str | None = None,
    external_inputs=(),
    reuse_paths=frozenset(),
    reuse_hashes=None,
    denylist_path=L.DEFAULT_DENYLIST,
    require_clean_tree: bool = True,
) -> dict:
    """Write ``<artifact>.stamp.json`` and return the stamp. See the module docstring for what is refused.

    ``reuse_hashes`` ({path: sha256}, optional) PINS the reuse: when given, every reuse path must have a pin and an input at that path
    must hash to it, else the stamp is refused (a named reuse is a named FILE, not a named path).

    ``repo_root`` = the namespace root (paths in the stamp are relative to it); ``git_root`` = the tree whose HEAD is the code that
    ran (``code_commit``, producer blob, clean-tree check). Both default via ``resolve_roots`` (``$V2_NAMESPACE_ROOT`` + this file's tree)."""
    root, git = resolve_roots(repo_root, git_root)
    art = Path(artifact)
    if not art.is_file() or art.is_symlink():
        raise StampError(f"artifact {art} must be a real (non-symlink) file")
    inputs = list(inputs)
    if not inputs:
        raise StampError("a stamp needs at least one input (input_sha256 must not be empty)")
    deny = L.load_denylist(denylist_path)
    reuse = frozenset(reuse_paths)
    pins = None if reuse_hashes is None else dict(reuse_hashes)
    jid = _slurm_job_id()
    producer_rel = _rel(git, producer)
    head, blob = _code_state(git, producer_rel, require_clean_tree)

    out_rel = _rel(root, art)
    problems = [f"artifact {x}" for x in L.path_problems(root, out_rel, check_tokens=True)]
    if pins is not None:
        problems += [f"reuse path {x} has no pinned sha256" for x in sorted(reuse - set(pins))]
    in_hashes = {}
    for p in inputs:
        rel = _rel(root, p)
        problems += [f"input {x}" for x in L.path_problems(root, rel, reuse_paths=reuse, check_tokens=True)]
        f = root / rel
        if not f.is_file():
            problems.append(f"input {rel} is not a file")
            continue
        h = L.sha256_file(f)
        if pins is not None and rel in reuse and rel in pins and h != pins[rel]:
            problems.append(f"input {rel}: reuse is pinned to sha256 {pins[rel][:12]} but the file is {h[:12]}")
        twin = git / rel  # a code file read from the namespace root must be the one the code tree ran with
        if git != root and twin.is_file() and L.sha256_file(twin) != h:
            problems.append(f"input {rel} differs between the namespace root {root} and the code tree {git}")
        if h in deny and rel not in reuse:
            problems.append(f"input {rel} is the v41 artifact {deny[h]} (denylist)")
        in_hashes[rel] = h
    out_h = L.sha256_file(art)
    if out_h in deny and out_rel not in reuse:
        problems.append(f"artifact {out_rel} has the hash of the v41 artifact {deny[out_h]} (denylist)")
    if problems:
        raise StampError("; ".join(problems))

    stamp = {
        **FROZEN,
        "run_id": run_id,
        "protocol_sha12": protocol_sha12,
        "code_commit": head,
        "output_sha256": out_h,
        "input_sha256": in_hashes,
        "producer_path": producer_rel,
        "producer_blob_sha": blob,
        "slurm_job_id": jid,
        "host": socket.gethostname(),
        "started_at": started_at,
        "ended_at": utc_now(),
    }
    if stage_role:
        stamp["stage_role"] = stage_role
    if external_inputs:
        stamp["external_inputs_sha256"] = {str(Path(p)): L.sha256_file(Path(p)) for p in external_inputs}
    side = Path(str(art) + ".stamp.json")
    tmp = side.with_name(side.name + f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(stamp, indent=2, sort_keys=True))
    os.replace(tmp, side)
    return stamp


def verify_stamp(artifact, *, run_id: str, protocol_sha12: str, repo_root=None, reuse_paths=frozenset(), denylist_path=L.DEFAULT_DENYLIST) -> dict:
    """Consumer-side gate: the artifact is a real file with a valid sidecar for ``run_id`` whose hashes still hold.

    Cheap (no ancestry): ``check_v42_lineage --stage complete`` does the full closure. Raises StampError."""
    root, art = resolve_roots(repo_root)[0], Path(artifact)
    rel = _rel(root, art)
    probs = L.path_problems(root, rel, reuse_paths=frozenset(reuse_paths))
    if probs:
        raise StampError("; ".join(probs))
    stamp = L.load_stamp(art) if (art.is_file() and L._has_sidecar(art)) else None
    probs = L.validate_stamp(stamp, {run_id}, protocol_sha12)
    if not probs:
        if stamp["output_sha256"] != L.sha256_file(art):
            probs.append("output_sha256 does not match the artifact")
        deny = L.load_denylist(denylist_path)
        if L.sha256_file(art) in deny and rel not in reuse_paths:
            probs.append(f"artifact is the v41 artifact {deny[L.sha256_file(art)]} (denylist)")
        for irel, want in stamp["input_sha256"].items():
            f = root / irel
            if not f.is_file() or L.sha256_file(f) != want:
                probs.append(f"input {irel} is missing or changed since stamping")
    if probs:
        raise StampError(f"{rel}: " + "; ".join(probs))
    return stamp
