#!/usr/bin/env python
"""Fail-closed lineage contracts for the manuscript v42 (frozen v2 Phase 1) data build.

Goal: a v42 build that reads ANY v41 number must fail, never render.

Usage
-----
    python scripts/check_v42_lineage.py --data-dir <built data dir> --stage scaffold|complete

  scaffold  structure only. Every data key is classified; every P/G/BOTH row has the expected fields;
            run ids exist; every reuse-allowlist row has a hash; the manifest covers the data keys on disk
            and the template placeholders, and nothing else.
  complete  scaffold + closure: every P/G/BOTH source exists on disk, is a real file (not a symlink, not
            under a v41 path), matches its sha256 and carries a valid stamp; the built data is labelset
            v42; no template cites a retired key or figure. Real-data closure is exercised in Task 11.

Files
-----
  v42_lineage_manifest.yaml (CITEGEIST_LINEAGE_MANIFEST)  one row per top-level key of every citable data JSON
  Benchmarking/midkine_remediation/v42_run_ids.yaml    {protocol_key: {B1..B12, B17..B26, B221, X1..X11, X51: run id}}

Classes: P proportions, G Xenium GEX scoring, BOTH, U unaffected.

Manifest row
------------
  expression         placeholder expressions in the v42 templates that read this key (derived, informational)
  data_key           "<json path under data dir without .json>::<top-level key>"
  class              P | G | BOTH | U
  expected_method    P/G/BOTH: "rna_refined"
  expected_run       P/G/BOTH: a stage id (e.g. "X2") or a list of stage ids
  expected_protocol  P/G/BOTH: a protocol key of v42_run_ids.yaml (breast stages B*) or "none" (Xenium stages X*)
  reuse_category     U only: gt | competitor | wet_lab | public_data | gene_removal (must match a reuse_allowlist entry by source).
                     The allowlist category `design` is NOT a U-row category: it names a hash-pinned experiment DESIGN file (the list of
                     perturbation conditions of the X7 protein stress sweep, byte-identical to the v41p design by intent, holding no result)
                     that a v2 stamp may take as an input; the stamp writer and the closure both compare the file with the pinned sha256.
  u_reason           U only, required: why the key is unaffected by the v2 Phase 1 swap. One of U_REASONS; "retired" iff
                     the row is retired. A P key relabelled U without a reason fails both stages; Task 11 also calls
                     check_u_unchanged(), which fails any U key whose built value differs from v41-final.
  source             repo-relative path of the artifact the key is built from (P/G/BOTH: a NEW v2 artifact; may contain
                     {protocol} and {B1}..{X9}, formatted from the run ids by --stage complete)
  baseline           U only, optional: "absent" = v42 ADDS this key (it is not in v41-final). `check_u_unchanged` then requires it in the new build and
                     fails if v41-final already has it, instead of reporting "not in the baseline".
  v42_expected       U only, optional: the exact value v42 builds for a key that CHANGED vs v41-final by design (figure panel relabel). `check_u_unchanged`
                     compares the new build to this pin instead of to v41-final, and fails if the pin equals the v41-final value.
  source_sha256      sha256 of `source` (null until the producing task fills it; required by --stage complete)
  retired            optional true: v41-lineage number; a v42 template must not cite it (checked by complete)

Stamp schema (a source's "stamp"; later tasks write it)
-------------------------------------------------------
A sidecar `<source>.stamp.json`, or an embedded top-level `_stamp` object in a JSON source:
  phase1_method  must equal "rna_refined"
  c              must equal 0.1
  n_passes       must equal 3
  anchor_thr     must equal 0
  run_id         must be one of the run ids of the row's expected_run stages (v42_run_ids.yaml)
  protocol_sha12 breast stages: the 12-hex sha the row's expected_protocol resolves to; Xenium stages: "none"
  code_commit    40-hex git commit of the producer; must exist in THIS repo (git cat-file -e <sha>^{commit})
  output_sha256  sha256 of the stamped artifact (sidecar stamp: the file bytes; embedded _stamp: canonical JSON of the
                 file without `_stamp`). `complete` recomputes it.
  input_sha256   {repo-relative path: sha256} of the producer's inputs, non-empty. `complete` recomputes every hash, rejects
                 symlinks / forbidden tokens, and (for every stage except the roots B2 and X1) requires at least one input
                 under another v2 run id, so a stage cannot be stamped without naming its v2 upstream.
  producer_path  repo-relative path of the producer script; producer_blob_sha is its git blob sha AT code_commit, which `complete`
                 recomputes with `git rev-parse <code_commit>:<producer_path>` (a stamp from a dirty or different producer fails)
  slurm_job_id   "<id>" or "<id>_<array task>" of the job that ran the producer
  host, started_at, ended_at   where and when (ISO 8601; started_at <= ended_at)
code_commit must be a descendant-or-equal of `min_code_commit` (v42_run_ids.yaml top-level key: the commit that holds the v2
producers) and must NOT be an ancestor-or-equal of the `v41-final` tag. Any stamped artifact, or any stamp input, whose sha256 is in
the committed v41 denylist (Benchmarking/midkine_remediation/v42_v41_denylist.json, built by scripts/build_v42_denylist.py) fails
unless its path is a named reuse-allowlist entry. Stamps are written by scripts/v42_stamp.py (the one shared writer).
Types are strict: c is a float, n_passes an int, anchor_thr an int/float; bool never passes.
A stamp that says marker_share (a v41 artifact under a v2-looking name) fails.
RESIDUAL (cannot be proved by a checker): someone who recomputes output_sha256 and input_sha256 over v41 bytes and
names real v2 upstream files still passes. Closing that needs the producer to run in a recorded job whose log/DVC pointer
the checker can read (B-stage closure in Task 3 / Task 8); this file proves consistency, not provenance.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = Path(os.environ.get("CITEGEIST_LINEAGE_MANIFEST", REPO / "v42_lineage_manifest.yaml"))
DEFAULT_RUN_IDS = REPO / "Benchmarking" / "midkine_remediation" / "v42_run_ids.yaml"
DEFAULT_MS_DIR = Path(os.environ.get("CITEGEIST_MS_DIR", REPO))
DEFAULT_DENYLIST = REPO / "Benchmarking" / "midkine_remediation" / "v42_v41_denylist.json"
BASELINE_REF = "v41-final"
SLURM_ID_RE = re.compile(r"^\d+(_\d+)?$")

PENDING = "PENDING_TASK3"
# B13-B16 are named in the scope map (trajectory under B9, B14/B16 unaffected by input hash, B15 Xenium annotation sensitivity) and have no id of their own here.
# B17 matched gene controls and B18 MDK-NCL cell-type map are exploratory checks promoted to reported stages (Section 2.4); their ids are v2b171 and v2b18.
# B19 ORA Test B on the rebuilt spots, B20 ORA on the SACE layers, B21 keratin presence and estrogen covariation (checks 7, 8) and B22 bivariate Moran's I
# robustness (checks 9, 10, 11, 13, 14, 15) and B23 single-covariate keratin association are post-hoc descriptive checks promoted the same way; their ids are v2b191, v2b201, v2b21, v2b22 and v2b23. They are descriptive only: a verdict label
# they print is reported, never a gate. Note: v2b21 and v2b22 are also the spelling of a first and second rerun of B2 (<stage id> + digit); B2 reruns use a generation timestamp, not a new id.
# X10 is the withheld-gene cross-fit (heldout_genes.json) and X11 the leave-one-region-out tuning audit (loro_audit.json); they are listed with them
# because they also reach the manuscript through ``stage_sources``.
PROMOTED_STAGES = ("B17", "B18", "B19", "B20", "B21", "B22", "B23", "X10", "X11")
# A promoted stage may be replaced in ``stage_sources`` by its registered rerun (the X5 -> X51 pattern): the row for B22 may be named B221 (its source then holds ``{B221}``), once the manifest is
# repointed. The superseded stage's stamp is then no longer walked by the closure (check_complete walks manifest rows only), so v2b22 (whose stamp pins the superseded v2b21 / v2b121) stays registered
# without being superseded. Exactly one row per promoted stage, under its own name or its rerun's, never both.
PROMOTED_RERUN = {"B22": "B221"}
HELD_OUT_STAGE = "X10"
HELD_OUT_OK = "HELD_OUT_HOLDS"
STAGES = tuple([f"B{i}" for i in range(1, 13)] + ["B17", "B18", "B19", "B20", "B21", "B22", "B23"] + [f"X{i}" for i in range(1, 10)] + ["X10", "X11", "X12", "X51", "B221", "B24", "B25", "B26"])
RUN_NAME_RE = re.compile(r"^v2[bx]\d+$")
SHA12_RE = re.compile(r"^[0-9a-f]{12}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
CLASSES = ("P", "G", "BOTH", "U")
AFFECTED = ("P", "G", "BOTH")
# `design` = a hash-pinned experiment design (condition list), reused on purpose by an X stage; see the header (reuse_category).
REUSE_CATEGORIES = ("gt", "competitor", "wet_lab", "public_data", "gene_removal", "design")
GENE_REMOVAL_METHODS = ("RCTD", "Cell2Location")
U_REASONS = ("competitor", "gt_only", "wet_lab_manual", "public_data", "cohort_constant", "dictionary_only", "protein_only", "layout_only", "retired")
ROOT_STAGES = ("B2", "X1", "B24")  # stages with no v2 upstream (B24: synthetic worlds; its stamp pins the protocol and the rna_total deviation file)
# v41 / v40 protocol hashes and run names: never valid for v42.
V41_SHA12 = frozenset({"ec0cf5b80715", "0344d54f9d3d"})
OLD_RUN_NAMES = frozenset(
    {"h1", "h2", "h3", "h3b", "h3c", "h3d", "ks1", "m1", "m2", "n1", "pc1", "pr2", "pr4", "pr5", "r3", "r4", "r5"}
)
# Path fragments that mark a v41 / P138 / stand-in artifact. Matched case-insensitively (see forbidden_hit).
FORBIDDEN_TOKENS = (
    "morphology_assignment_v3",
    "ec0cf5b80715",
    "0344d54f9d3d",
    "p138",
    "rnaprot_v41",
    "sensitivity_v41p",
    "marker_share",
    "output_sc_assignment_p138",
    "results_v41",
    "olddict",
)
# Built data JSONs that legitimately carry no `_generated_at`/`_sources`. Anything else lacking them fails `complete`.
UNSTAMPED_DATA_FILES = {
    "benchmarks/marker_inventory": "hand-authored panel-size audit (input to marker_inventory_resolved); counts antibody panels, no Phase 1 output",
    "refs": "figure-panel letter map extracted from the figure manifests; carries _sources but no numbers",
    "geomx_benchmark/pseudo_geomx_gex": "retired GeoMx stub holding only _source/_source_md5; no citable keys",
}
# Data files that load_data() puts in the render namespace but that hold no citable numbers.
NON_CITABLE_FILE_RE = re.compile(
    r"^(v3\d_.*|references|zotero_uri_map|MANIFEST)$|^_archive"
)
# Retired figures: not in the v42 supplementary list (see supplementary_order_for_v42).
STAMP_FIELDS = (
    "phase1_method", "c", "n_passes", "anchor_thr", "run_id", "protocol_sha12", "code_commit", "output_sha256", "input_sha256",
    "producer_path", "producer_blob_sha", "slurm_job_id", "host", "started_at", "ended_at",
)


class LineageError(RuntimeError):
    """A lineage contract failed. `.problems` lists every individual failure."""

    def __init__(self, problems):
        self.problems = list(problems)
        super().__init__("; ".join(self.problems))


class PendingProtocol(LineageError):
    pass


class UnknownProtocol(LineageError):
    pass


class MissingStage(LineageError):
    pass


class NotRun(LineageError):
    pass


class OldRunName(LineageError):
    pass


class MissingInput(LineageError):
    pass


def _err(cls, msg):
    return cls([msg])


# ---------------------------------------------------------------- run ids
def load_run_ids(path: Path = DEFAULT_RUN_IDS) -> dict:
    return yaml.safe_load(Path(path).read_text())["runs"]


V2_TS_RE = re.compile(r"^\d{8}T\d{6}Z$")


def load_generation(path: Path = DEFAULT_RUN_IDS) -> dict:
    """The one record of which B1/B2 input generation is current: ``{"current_v2_ts", "superseded_v2_ts", "superseded_run_ids"}`` (all from the
    run-ids file). Inputs stamped under a superseded V2_TS (the B1 refit / B2 nuclei dirs) or a superseded run id (B3-B5) must never be mixed in."""
    d = yaml.safe_load(Path(path).read_text())
    sup = d.get("superseded") or {}
    return {"current_v2_ts": d.get("current_v2_ts"), "superseded_v2_ts": [str(x) for x in sup.get("v2_ts") or []],
            "superseded_run_ids": [str(x) for x in sup.get("run_ids") or []]}


def superseded_path_parts(rel: str, gen: dict) -> list:
    """Components of ``rel`` that are a superseded run id, or end in a superseded V2_TS (``module3_rna_refined_7type_<ts>``)."""
    out = []
    for part in Path(rel).parts:
        if part in gen["superseded_run_ids"]:
            out.append(f"run id {part}")
        for ts in gen["superseded_v2_ts"]:
            if part == ts or part.endswith("_" + ts):
                out.append(f"V2_TS {ts}")
    return out


def generation_problems(gen: dict, runs: dict) -> list:
    p = []
    cur = gen.get("current_v2_ts")
    if cur is not None and not V2_TS_RE.match(str(cur)):
        p.append(f"current_v2_ts {cur!r} is not a UTC timestamp like 20261102T200000Z")
    if cur is not None and str(cur) in gen["superseded_v2_ts"]:
        p.append(f"current_v2_ts {cur} is also listed as superseded")
    for t in gen["superseded_v2_ts"]:
        if not V2_TS_RE.match(t):
            p.append(f"superseded v2_ts {t!r} is not a UTC timestamp")
    for key, stages in runs.items():
        owner = {}
        for st, n in (stages or {}).items():
            if n in owner:     # e.g. a B1 rerun named v2b11 would equal the live B11 id: the "<stage>1" scheme cannot tell them apart
                p.append(f"run id {n!r} is assigned to both {owner[n]} and {st} under protocol {key} (ids must be unique per protocol)")
            owner[n] = st
    live = {n for stages in runs.values() for n in (stages or {}).values()}
    for n in gen["superseded_run_ids"]:
        if n in live:
            p.append(f"run id {n!r} is superseded but is still assigned in runs")
    return p


def load_min_code_commit(path: Path = DEFAULT_RUN_IDS):
    """The commit that holds the v2 producers (top-level `min_code_commit` of the run-ids file), or None."""
    return yaml.safe_load(Path(path).read_text()).get("min_code_commit")


def load_stage_min_code_commit(path: Path = DEFAULT_RUN_IDS) -> dict:
    """{stage: commit} from the optional top-level `stage_min_code_commit` of the run-ids file: a stage's stamp must be built from code that
    descends from this commit (a producer fix that lives outside the stamped producer file, e.g. B9/B12 and src/morans.py)."""
    return {str(k): str(v) for k, v in (yaml.safe_load(Path(path).read_text()).get("stage_min_code_commit") or {}).items()}


def stage_min_problems(stage: str, code_commit: str, stage_min: dict | None, git_root: Path) -> list:
    """The stage's code_commit must descend from its `stage_min_code_commit` (no entry for the stage = no extra rule)."""
    sm = (stage_min or {}).get(stage)
    if not sm:
        return []
    if not COMMIT_RE.match(sm) or not commit_exists(sm, git_root):
        return [f"stage_min_code_commit[{stage}] {sm!r} (v42_run_ids.yaml) is not a commit in this repo"]
    if is_ancestor(sm, code_commit, git_root) is not True:
        return [f"code_commit {code_commit} is not a descendant of stage_min_code_commit[{stage}] {sm} (the stage's fix is not in the code that made it)"]
    return []


def load_denylist(path: Path = DEFAULT_DENYLIST) -> dict:
    """{sha256: repo-relative path} of the v41 artifacts that feed P, G or BOTH keys. Missing or empty fails closed."""
    p = Path(path)
    if not p.is_file():
        raise _err(LineageError, f"v41 denylist {p} is missing (build it with scripts/build_v42_denylist.py)")
    d = json.loads(p.read_text()).get("sha256")
    if not isinstance(d, dict) or not d or not all(SHA256_RE.match(str(k)) for k in d):
        raise _err(LineageError, f"v41 denylist {p} is empty or malformed (a vacuous denylist is not accepted)")
    return d


def resolve_protocol(sha12: str, runs: dict) -> dict:
    """Return the stage->run id map for protocol `sha12`; raise rather than fall back."""
    if sha12 == PENDING:
        raise _err(PendingProtocol, f"protocol is still {PENDING}: Task 3 must replace it with the real sha12")
    if sha12 in V41_SHA12:
        raise _err(UnknownProtocol, f"protocol {sha12!r} is a v40/v41 protocol, not valid for v42")
    if not SHA12_RE.match(str(sha12)) or sha12 not in runs:
        raise _err(UnknownProtocol, f"unknown protocol sha12 {sha12!r}; known: {sorted(runs)}")
    return runs[sha12]


def resolve_run(sha12: str, stage: str, runs: dict) -> str:
    """Run id for (protocol, stage). Unknown sha, missing stage, NOT_RUN and old run names all raise."""
    stages = resolve_protocol(sha12, runs)
    if stage not in stages:
        raise _err(MissingStage, f"protocol {sha12} has no run id for stage {stage}")
    name = stages[stage]
    if name is None or str(name) == "NOT_RUN":
        raise _err(NotRun, f"stage {stage} of protocol {sha12} is NOT_RUN")
    if name in OLD_RUN_NAMES or not RUN_NAME_RE.match(str(name)):
        raise _err(OldRunName, f"stage {stage} of protocol {sha12}: {name!r} is not a v2 run name (^v2[bx]\\d+$)")
    return name


def require_protocol_ready(runs: dict, sha12: str | None = None) -> str:
    """Gate for a v42 data build: the active protocol must be a real sha with every stage resolvable.

    With no `sha12`, requires exactly one protocol key in the run-ids file and returns it."""
    if sha12 is None:
        if len(runs) != 1:
            raise _err(UnknownProtocol, f"expected exactly one protocol key in v42_run_ids.yaml, found {sorted(runs)}")
        (sha12,) = runs
    for stage in STAGES:
        resolve_run(sha12, stage, runs)
    return sha12


def check_run_ids(runs: dict) -> list:
    problems = []
    if not runs:
        problems.append("v42_run_ids.yaml has no protocol key")
    seen = {}
    for proto, stages in runs.items():
        if proto != PENDING and not SHA12_RE.match(str(proto)):
            problems.append(f"protocol key {proto!r} is neither {PENDING} nor 12 hex")
        if proto in V41_SHA12:
            problems.append(f"protocol key {proto!r} is a v40/v41 protocol")
        missing = [s for s in STAGES if s not in (stages or {})]
        if missing:
            problems.append(f"protocol {proto}: missing stages {missing}")
        for s, name in (stages or {}).items():
            if s not in STAGES:
                problems.append(f"protocol {proto}: unknown stage {s}")
            if name == "NOT_RUN":
                continue
            if name in OLD_RUN_NAMES or not RUN_NAME_RE.match(str(name)):
                problems.append(f"protocol {proto} stage {s}: {name!r} is not a v2 run name")
            if name in seen:
                problems.append(f"run id {name!r} used by both {seen[name]} and {proto}/{s}")
            seen[name] = f"{proto}/{s}"
    return problems


def require_inputs(paths: dict) -> None:
    """Every non-None v42 source path must exist. A missing one raises; the caller never substitutes another path."""
    missing = sorted(f"{k}: {p}" for k, p in paths.items() if p is not None and not Path(p).exists())
    if missing:
        raise MissingInput([f"missing v2 input {m}" for m in missing])


# ---------------------------------------------------------------- data keys
def _file_id(data_dir: Path, path: Path) -> str:
    return str(path.relative_to(data_dir).with_suffix(""))


def data_keys(data_dir: Path) -> set:
    """Citable top-level keys on disk: {"<file id>::<key>"} (keys starting with "_" are never rendered)."""
    data_dir = Path(data_dir)
    keys = set()
    for p in sorted(data_dir.rglob("*.json")):
        rel = p.relative_to(data_dir)
        if "figures" in rel.parts or NON_CITABLE_FILE_RE.match(rel.stem):
            continue
        d = json.loads(p.read_text())
        if isinstance(d, dict):
            keys.update(f"{_file_id(data_dir, p)}::{k}" for k in d if not str(k).startswith("_"))
    return keys


def _chains(ast) -> list:
    """Maximal dotted/subscript chains rooted at a Name, as tuples of strings (Jinja AST, so {{ }} and {% %} alike)."""
    from jinja2 import nodes

    inner = set()
    for n in ast.find_all((nodes.Getattr, nodes.Getitem)):
        inner.add(id(n.node))
    out = []
    for n in ast.find_all((nodes.Getattr, nodes.Getitem, nodes.Name)):
        if id(n) in inner:
            continue
        parts, cur = [], n
        while True:
            if isinstance(cur, nodes.Getattr):
                parts.append(cur.attr)
                cur = cur.node
            elif isinstance(cur, nodes.Getitem) and isinstance(cur.arg, nodes.Const) and isinstance(cur.arg.value, str):
                parts.append(cur.arg.value)
                cur = cur.node
            elif isinstance(cur, nodes.Getitem):  # dynamic subscript: the chain below it is what is read
                parts = []
                cur = cur.node
            elif isinstance(cur, nodes.Name):
                parts.append(cur.name)
                break
            else:
                parts = []
                break
        if parts:
            out.append(tuple(reversed(parts)))
    return out


def template_refs(template_text: str, file_ids: set, keys_by_file: dict | None = None) -> tuple:
    """-> (n_placeholders, {data_key: [expressions]}). `file_ids` are "dir/stem" ids of the data files.

    Scans the Jinja2 AST, so `{{ }}` expressions, `{% if/for/set %}` blocks, filters and subscripts are all covered.
    n_placeholders counts "{{" and "{%" independently of the parser, so a broken scan cannot pass vacuously.
    A reference to a whole data file (e.g. `{% for k in clinical.cohort %}`) expands to every key in `keys_by_file[file]`
    (or to "<file>::*" when no key map is given)."""
    from jinja2 import Environment, TemplateSyntaxError

    n_open = template_text.count("{{") + template_text.count("{%")
    try:
        ast = Environment().parse(template_text)
    except TemplateSyntaxError as e:
        raise _err(LineageError, f"template does not parse (line {e.lineno}): {e.message}")
    by_parts = {tuple(f.split("/")): f for f in file_ids}
    refs: dict = {}
    for parts in _chains(ast):
        label = ".".join(parts)
        for n in range(len(parts) - 1, 0, -1):
            fid = by_parts.get(parts[:n])
            if fid:
                refs.setdefault(f"{fid}::{parts[n]}", []).append(label)
                break
        else:
            fid = by_parts.get(parts)
            if fid:
                for k in sorted((keys_by_file or {}).get(fid, {f"{fid}::*"})):
                    refs.setdefault(k, []).append(label)
    return n_open, refs


def templates_of(ms_dir: Path) -> list:
    ms_dir = Path(ms_dir)
    # v43 = the author's hand edits ported onto the same v42 data; both documents read the v42 data keys
    docs = sorted(ms_dir.glob("CITEgeist_CRM_v42*.md")) + sorted(ms_dir.glob("CITEgeist_CRM_v43*.md")) + sorted(ms_dir.glob("CITEgeist_CRM_v44*.md"))
    return docs + sorted((ms_dir / "response").glob("v42_response_to_reviewers.md"))


# ---------------------------------------------------------------- manifest
def load_manifest(path: Path = DEFAULT_MANIFEST) -> dict:
    return yaml.safe_load(Path(path).read_text())


def forbidden_hit(text: str):
    """The first FORBIDDEN_TOKEN in `text` (case-insensitive), or None."""
    low = str(text).lower()
    return next((t for t in FORBIDDEN_TOKENS if t in low), None)


def _stages_of(row) -> list:
    e = row.get("expected_run")
    return [e] if isinstance(e, str) else list(e or [])


def check_manifest_rows(manifest: dict, runs: dict, provenance_allowlist_paths=frozenset()) -> list:
    problems = []
    rows = manifest.get("rows") or []
    allow = manifest.get("reuse_allowlist") or []
    allow_by_path = {}
    for a in allow:
        if a.get("category") not in REUSE_CATEGORIES:
            problems.append(f"reuse_allowlist {a.get('path')!r}: category {a.get('category')!r} not in {REUSE_CATEGORIES}")
        if not SHA256_RE.match(str(a.get("sha256"))):
            problems.append(f"reuse_allowlist {a.get('path')!r}: missing/invalid sha256")
        if a.get("category") == "gene_removal" and a.get("method") not in GENE_REMOVAL_METHODS:
            problems.append(f"reuse_allowlist {a.get('path')!r}: gene_removal reuse is limited to RCTD and Cell2Location (method={a.get('method')!r})")
        allow_by_path[a.get("path")] = a
    seen = set()
    for r in rows:
        k = r.get("data_key")
        if k in seen:
            problems.append(f"duplicate data_key {k}")
        seen.add(k)
        cls = r.get("class")
        if cls not in CLASSES:
            problems.append(f"{k}: unclassified (class={cls!r})")
            continue
        src = r.get("source")
        if cls in AFFECTED:
            if r.get("expected_method") != "rna_refined":
                problems.append(f"{k}: expected_method must be rna_refined")
            stages = _stages_of(r)
            if not stages or any(s not in STAGES or s == "X9" for s in stages):
                problems.append(f"{k}: expected_run {r.get('expected_run')!r} must be stage ids from {[x for x in STAGES if x != 'X9']}")
            else:
                proto = r.get("expected_protocol")
                is_b = any(s.startswith("B") for s in stages)
                if is_b and not all(s.startswith("B") for s in stages):
                    problems.append(f"{k}: expected_run {stages} mixes breast (B) and Xenium (X) stages; split the row")
                if is_b and proto not in runs:
                    problems.append(f"{k}: expected_protocol {proto!r} is not a protocol key of the run ids")
                if not is_b and proto != "none":
                    problems.append(f"{k}: Xenium-only row must have expected_protocol 'none'")
            if r.get("reuse_category"):
                problems.append(f"{k}: class {cls} must not use reuse_category {r['reuse_category']!r} (cached/reused source)")
            if not src:
                problems.append(f"{k}: class {cls} needs a source path")
            elif src in provenance_allowlist_paths:
                problems.append(f"{k}: class {cls} source {src} is in provenance_allowlist.yaml")
            if src and forbidden_hit(src):
                problems.append(f"{k}: class {cls} source {src} names a v41/stand-in artifact")
            if r.get("retired"):
                problems.append(f"{k}: a retired row must be class U")
            if r.get("u_reason") is not None:
                problems.append(f"{k}: class {cls} must not carry a u_reason")
            if "v42_expected" in r:
                problems.append(f"{k}: class {cls} must not carry v42_expected (only a U row can be pinned against v41-final)")
            if r.get("baseline") is not None:
                problems.append(f"{k}: class {cls} must not carry a baseline marker (only a U row can be absent from v41-final)")
        else:
            reason = r.get("u_reason")
            if not reason:
                problems.append(f"{k}: U row needs a u_reason from {U_REASONS}")
            elif reason not in U_REASONS:
                problems.append(f"{k}: u_reason {reason!r} not in {U_REASONS}")
            elif r.get("retired") and reason != "retired":
                problems.append(f"{k}: a retired row must have u_reason 'retired' (got {reason!r})")
            elif reason == "retired" and not r.get("retired"):
                problems.append(f"{k}: u_reason 'retired' on a row that is not retired")
            if "v42_expected" in r and r.get("baseline") is not None:
                problems.append(f"{k}: v42_expected and baseline are exclusive")
            if r.get("baseline") not in (None, "absent"):
                problems.append(f"{k}: baseline {r.get('baseline')!r} must be 'absent' or omitted")
            cat = r.get("reuse_category")
            if cat is not None and cat not in REUSE_CATEGORIES:
                problems.append(f"{k}: reuse_category {cat!r} outside the allowlist categories")
            if src:
                a = allow_by_path.get(src)
                if cat is None or a is None:
                    problems.append(f"{k}: reuses {src} outside the reuse allowlist")
                elif a.get("category") != cat:
                    problems.append(f"{k}: reuse_category {cat!r} != allowlist category {a.get('category')!r}")
                elif r.get("source_sha256") != a.get("sha256"):
                    problems.append(f"{k}: source_sha256 does not match the reuse allowlist hash")
            elif cat is not None:
                problems.append(f"{k}: reuse_category without a source")
    return problems


def check_stage_sources(manifest: dict, runs: dict) -> list:
    """Structure of the manifest's top-level ``stage_sources`` list: one row per promoted stage (B17, B18), each naming the stamped artifact the
    stage hands to a consumer. These stages feed manuscript keys that Task 11 builds, so they have no ``data_key`` row yet (``check_coverage``
    would call it a key missing from the data). A row: ``stage``, ``expected_protocol`` (a protocol key of the run ids), ``source`` (contains
    ``{protocol}`` and the stage's ``{B17}``/``{B18}`` placeholder, no v41 token), ``source_sha256`` (null until the stage ran). ``--stage complete``
    runs the same stamp closure on them as on a P row, and a stage that is missing from the run ids or NOT_RUN raises there (resolve_run)."""
    problems = []
    rows = manifest.get("stage_sources") or []
    seen = [r.get("stage") for r in rows]
    for st in PROMOTED_STAGES:
        n_rows = seen.count(st) + (seen.count(PROMOTED_RERUN[st]) if st in PROMOTED_RERUN else 0)
        if n_rows != 1:
            problems.append(f"stage_sources: stage {st} must have exactly one row (found {n_rows})")
    for r in rows:
        st, src = r.get("stage"), r.get("source")
        if st not in PROMOTED_STAGES and st not in PROMOTED_RERUN.values():
            problems.append(f"stage_sources: {st!r} is not a promoted stage {PROMOTED_STAGES}")
            continue
        is_x = st.startswith("X")  # X stamps carry protocol "none" (vo.stamp_artifact), like the X rows
        if is_x:
            if r.get("expected_protocol") != "none":
                problems.append(f"stage_sources {st}: an X stage must have expected_protocol 'none' (got {r.get('expected_protocol')!r})")
        elif r.get("expected_protocol") not in runs:
            problems.append(f"stage_sources {st}: expected_protocol {r.get('expected_protocol')!r} is not a protocol key of the run ids")
        need = ["{" + st + "}"] + ([] if is_x else ["{protocol}"])
        if not src or any(n not in src for n in need):
            problems.append(f"stage_sources {st}: source {src!r} must contain {' and '.join(need)}")
        elif forbidden_hit(src):
            problems.append(f"stage_sources {st}: source {src} names a v41/stand-in artifact")
        sha = r.get("source_sha256")
        if sha is not None and not SHA256_RE.match(str(sha)):
            problems.append(f"stage_sources {st}: source_sha256 is neither null nor a sha256")
    return problems


def held_out_verdict_problems(src: Path) -> list:
    """``heldout_genes.json`` may reach the manuscript only with verdict HELD_OUT_HOLDS. HELD_OUT_WORSE or NO_VERDICT is a real, stamped record of
    the pre-registered test, so it is not a lineage defect: the build stops and the author decides what the manuscript says."""
    try:
        v = (json.loads(Path(src).read_text()).get("verdict") or {}).get("verdict")
    except (OSError, ValueError, AttributeError) as e:
        return [f"{src}: cannot read the verdict ({e}); needs an author decision"]
    if v != HELD_OUT_OK:
        return [f"{src}: withheld-gene verdict is {v!r}, not {HELD_OUT_OK}; this needs an author decision (report the numbers, do not tune or rephrase)"]
    return []


def stage_source_rows(manifest: dict) -> list:
    """The ``stage_sources`` rows in the shape ``check_complete`` closes over (class P, one stage, no data key on disk)."""
    return [{"data_key": f"stage_sources::{r['stage']}", "class": "P", "expected_run": r["stage"], "expected_protocol": r["expected_protocol"],
             "source": r["source"], "source_sha256": r.get("source_sha256")} for r in manifest.get("stage_sources") or []]


def check_coverage(manifest: dict, data_dir: Path, templates: list) -> list:
    problems = []
    rows = {r.get("data_key") for r in manifest.get("rows") or []}
    on_disk = data_keys(data_dir)
    if not on_disk:
        problems.append(f"no citable data keys found under {data_dir} (vacuous coverage check)")
    if not rows:
        problems.append("manifest has no rows")
    problems += [f"data key missing from manifest: {k}" for k in sorted(on_disk - rows)]
    problems += [f"manifest key missing from data: {k}" for k in sorted(rows - on_disk)]
    file_ids = {k.split("::")[0] for k in on_disk | rows}
    if not templates:
        problems.append("no v42 templates found")
    keys_by_file: dict = {}
    for k in on_disk | rows:
        keys_by_file.setdefault(k.split("::")[0], set()).add(k)
    for t in templates:
        text = Path(t).read_text()
        n_open, refs = template_refs(text, file_ids, keys_by_file)
        if n_open == 0 or not refs:
            problems.append(f"{Path(t).name}: no data placeholders found (n={n_open})")
        problems += [f"{Path(t).name}: placeholder reads {k} not in manifest" for k in sorted(set(refs) - rows)]
        problems += [f"{Path(t).name}: placeholder reads {k} not on disk" for k in sorted(set(refs) - on_disk - rows)]
    return problems


def check_scaffold(data_dir: Path, manifest_path=DEFAULT_MANIFEST, run_ids_path=DEFAULT_RUN_IDS, ms_dir=DEFAULT_MS_DIR) -> None:
    manifest = load_manifest(manifest_path)
    runs = load_run_ids(run_ids_path)
    allow_paths = _provenance_paths(Path(data_dir))
    problems = check_run_ids(runs)
    problems += generation_problems(load_generation(run_ids_path), runs)
    problems += check_manifest_rows(manifest, runs, allow_paths)
    problems += check_stage_sources(manifest, runs)
    problems += check_coverage(manifest, Path(data_dir), templates_of(ms_dir))
    if problems:
        raise LineageError(problems)


def _provenance_paths(data_dir: Path) -> frozenset:
    p = Path(data_dir) / "provenance_allowlist.yaml"
    if not p.exists():
        return frozenset()
    return frozenset(e.get("path") for e in (yaml.safe_load(p.read_text()) or []))


# ---------------------------------------------------------------- stamps and closure
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_sha256(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _has_sidecar(source: Path) -> bool:
    return Path(str(source) + ".stamp.json").exists()


def load_stamp(source: Path):
    side = Path(str(source) + ".stamp.json")
    if side.exists():
        return json.loads(side.read_text())
    if source.suffix == ".json":
        d = json.loads(source.read_text())
        if isinstance(d, dict) and isinstance(d.get("_stamp"), dict):
            return d["_stamp"]
    return None


def artifact_sha256(source: Path) -> str:
    """Hash a stamp covers: the file bytes (sidecar stamp), or canonical JSON without `_stamp` (embedded stamp)."""
    if not _has_sidecar(source) and source.suffix == ".json":
        d = json.loads(source.read_text())
        if isinstance(d, dict) and "_stamp" in d:
            return canonical_sha256({k: v for k, v in d.items() if k != "_stamp"})
    return sha256_file(source)


def _is_num(x, typ) -> bool:
    return type(x) in typ


def validate_stamp(stamp, run_names: set, protocol_sha12: str) -> list:
    """Schema and value checks on a stamp (types are strict: n_passes 3.0, anchor_thr False, c True all fail)."""
    if not isinstance(stamp, dict):
        return ["no stamp (sidecar .stamp.json or embedded _stamp)"]
    p = [f"stamp missing {f}" for f in STAMP_FIELDS if f not in stamp]
    if p:
        return p
    if stamp["phase1_method"] != "rna_refined":
        p.append(f"phase1_method {stamp['phase1_method']!r} != rna_refined")
    if not (_is_num(stamp["c"], (float,)) and stamp["c"] == 0.1):
        p.append(f"c {stamp['c']!r} != 0.1 (float)")
    if not (_is_num(stamp["n_passes"], (int,)) and stamp["n_passes"] == 3):
        p.append(f"n_passes {stamp['n_passes']!r} != 3 (int)")
    if not (_is_num(stamp["anchor_thr"], (int, float)) and stamp["anchor_thr"] == 0):
        p.append(f"anchor_thr {stamp['anchor_thr']!r} != 0 (int or float, not bool)")
    if stamp["run_id"] not in run_names:
        p.append(f"run_id {stamp['run_id']!r} not in {sorted(run_names)}")
    if stamp["protocol_sha12"] != protocol_sha12:
        p.append(f"protocol_sha12 {stamp['protocol_sha12']!r} != {protocol_sha12!r}")
    if not COMMIT_RE.match(str(stamp["code_commit"])):
        p.append("code_commit is not a 40-hex commit")
    if not SHA256_RE.match(str(stamp["output_sha256"])):
        p.append("output_sha256 is not a 64-hex sha256")
    pp = stamp["producer_path"]
    if not isinstance(pp, str) or not pp or Path(pp).is_absolute() or ".." in Path(pp).parts:
        p.append("producer_path must be a repo-relative path")
    if not COMMIT_RE.match(str(stamp["producer_blob_sha"])):
        p.append("producer_blob_sha is not a 40-hex git blob sha")
    if not SLURM_ID_RE.match(str(stamp["slurm_job_id"])):
        p.append(f"slurm_job_id {stamp['slurm_job_id']!r} is not a SLURM job id (<id> or <id>_<task>)")
    if not isinstance(stamp["host"], str) or not stamp["host"]:
        p.append("host must be a non-empty string")
    try:
        t0, t1 = _dt.datetime.fromisoformat(str(stamp["started_at"])), _dt.datetime.fromisoformat(str(stamp["ended_at"]))
        if t0.tzinfo != t1.tzinfo and (t0.tzinfo is None or t1.tzinfo is None):
            p.append("started_at/ended_at must both carry (or both lack) a timezone")
        elif t0 > t1:
            p.append("started_at is after ended_at")
    except ValueError:
        p.append("started_at/ended_at are not ISO 8601 timestamps")
    ins = stamp["input_sha256"]
    if not isinstance(ins, dict) or not ins:
        p.append("input_sha256 must be a non-empty {path: sha256} map")
    elif not all(isinstance(k, str) and SHA256_RE.match(str(v)) for k, v in ins.items()):
        p.append("input_sha256 has an entry that is not {path: 64-hex sha256}")
    return p


_COMMIT_CACHE: dict = {}


def commit_exists(sha: str, git_root: Path) -> bool:
    key = (str(git_root), sha)
    if key not in _COMMIT_CACHE:
        r = subprocess.run(["git", "-C", str(git_root), "cat-file", "-e", f"{sha}^{{commit}}"], capture_output=True)
        _COMMIT_CACHE[key] = r.returncode == 0
    return _COMMIT_CACHE[key]


def _git_rc(git_root: Path, *args):
    r = subprocess.run(["git", "--no-optional-locks", "-C", str(git_root), *args], capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def is_ancestor(a: str, b: str, git_root: Path):
    """True/False for `a` is an ancestor-or-equal of `b`; None when git cannot decide (unknown ref)."""
    rc, _ = _git_rc(git_root, "merge-base", "--is-ancestor", a, b)
    return True if rc == 0 else False if rc == 1 else None


def code_provenance_problems(stamp: dict, git_root: Path, min_code_commit, baseline_ref: str) -> list:
    """producer blob sha at code_commit, descendant of min_code_commit, not an ancestor of the v41 baseline tag."""
    p = []
    cc = stamp["code_commit"]
    rc, blob = _git_rc(git_root, "rev-parse", "--verify", f"{cc}:{stamp['producer_path']}")
    if rc == 0 and _git_rc(git_root, "cat-file", "-t", blob)[1] != "blob":
        rc = 1
    if rc != 0:
        p.append(f"producer {stamp['producer_path']} does not exist at code_commit {cc}")
    elif blob != stamp["producer_blob_sha"]:
        p.append(f"producer_blob_sha {stamp['producer_blob_sha']} != {blob} (git blob of {stamp['producer_path']} at code_commit)")
    if not min_code_commit or not COMMIT_RE.match(str(min_code_commit)) or not commit_exists(min_code_commit, git_root):
        p.append(f"min_code_commit {min_code_commit!r} (v42_run_ids.yaml) is not a commit in this repo")
    elif is_ancestor(min_code_commit, cc, git_root) is not True:
        p.append(f"code_commit {cc} is not a descendant of min_code_commit {min_code_commit} (producer predates the v2 code)")
    base = is_ancestor(cc, baseline_ref, git_root)
    if base is None:
        p.append(f"cannot compare code_commit with {baseline_ref} (ref missing); failing closed")
    elif base:
        p.append(f"code_commit {cc} is an ancestor of {baseline_ref}: the producer is v41 code")
    return p


def path_problems(repo_root: Path, rel: str, reuse_paths=frozenset(), check_tokens: bool = True) -> list:
    """Problems with a repo-relative path: absolute / `..`, any symlink component up to the repo root (not only the leaf),
    a resolved path that differs from the lexical one, a forbidden token in the lexical path, and a forbidden token in the
    fully RESOLVED path. The checks do not short-circuit each other, so each one is independently load-bearing."""
    probs = []
    root = Path(repo_root)
    rp = Path(rel)
    if rp.is_absolute() or ".." in rp.parts:
        return [f"path {rel} is absolute or contains '..'"]
    cur = root
    for part in rp.parts:
        cur = cur / part
        if cur.is_symlink():
            probs.append(f"path {rel}: component {cur.relative_to(root)} is a symlink")
    if check_tokens and rel not in reuse_paths and forbidden_hit(rel):
        probs.append(f"path {rel} names a v41/stand-in artifact ({forbidden_hit(rel)})")
    real, real_root = (root / rel).resolve(), root.resolve()
    try:
        real_rel = str(real.relative_to(real_root))
    except ValueError:
        real_rel = str(real)
        probs.append(f"path {rel} resolves outside the repo to {real}")
    if real != real_root / rel:
        probs.append(f"path {rel} resolves to {real_rel}, not to itself")
    hit = forbidden_hit(real_rel)
    if check_tokens and rel not in reuse_paths and hit:
        probs.append(f"path {rel} resolves to a v41/stand-in path {real} ({hit})")
    return probs


def stamp_closure_problems(
    stamp: dict,
    src: Path,
    repo_root: Path,
    git_root: Path,
    reuse_paths,
    run_names: set,
    own_stage: str,
    min_code_commit=None,
    baseline_ref: str = BASELINE_REF,
    denylist: dict | None = None,
    generation: dict | None = None,
    stage_min_code_commit: dict | None = None,
    reuse_hashes: dict | None = None,
) -> list:
    """Recompute everything a stamp asserts: output hash, commit (+ blob, ancestry), input hashes, denylist, v2 upstream, and that no path
    (the artifact or any input) sits under a SUPERSEDED run id or V2_TS of the generation record (``load_generation``, v42_run_ids.yaml).

    By design only this closure (``--stage complete``) applies the generation rule: ``--stage scaffold`` promises the structure of the run-ids
    file and the manifest, and the committed old-generation stamps are history that no manifest source reaches, so they must not fail it."""
    p = []
    gen = generation if generation is not None else load_generation()
    try:
        art_rel = str(Path(src).resolve().relative_to(Path(repo_root).resolve()))
    except ValueError:
        art_rel = str(src)
    for what, rel in [("artifact", art_rel)] + [("input", r) for r in stamp["input_sha256"]]:
        bad = superseded_path_parts(rel, gen)
        if bad:
            p.append(f"{what} {rel} is under the superseded {bad} (generation record in v42_run_ids.yaml)")
    if stamp["output_sha256"] != artifact_sha256(src):
        p.append("stamp output_sha256 does not match the artifact (v41 bytes under a copied stamp?)")
    if denylist is not None:
        for h in {stamp["output_sha256"], sha256_file(src)}:
            if h in denylist:
                p.append(f"artifact hash {h[:12]} is the v41 artifact {denylist[h]} (denylist)")
    if not commit_exists(stamp["code_commit"], git_root):
        p.append(f"code_commit {stamp['code_commit']} does not resolve to a commit in this repo")
    else:
        p += code_provenance_problems(stamp, git_root, min_code_commit, baseline_ref)
        p += stage_min_problems(own_stage, stamp["code_commit"], stage_min_code_commit, git_root)
    upstream = False
    for rel, want in stamp["input_sha256"].items():
        pp = path_problems(repo_root, rel, reuse_paths)
        p += [f"stamp input {x}" for x in pp]
        f = Path(repo_root) / rel
        if not f.is_file():
            p.append(f"stamp input {rel} not on disk")
        elif not pp and sha256_file(f) != want:
            p.append(f"stamp input {rel}: input_sha256 does not match the file")
        if reuse_hashes is not None and rel in reuse_hashes and want != reuse_hashes[rel]:
            p.append(f"stamp input {rel}: reused input is pinned to sha256 {reuse_hashes[rel][:12]} in the reuse allowlist but the stamp recorded {want[:12]}")
        if denylist is not None and rel not in reuse_paths and want in denylist:
            p.append(f"stamp input {rel}: hash is the v41 artifact {denylist[want]} (denylist)")
        if any(RUN_NAME_RE.match(x) and x != stamp["run_id"] for x in Path(rel).parts) and rel not in reuse_paths:
            upstream = True
    if own_stage not in ROOT_STAGES and not upstream:
        p.append(f"stamp for stage {own_stage} has no v2 upstream input (an input path under another v2 run id)")
    return p


def check_complete(
    data_dir: Path,
    manifest_path=DEFAULT_MANIFEST,
    run_ids_path=DEFAULT_RUN_IDS,
    ms_dir=DEFAULT_MS_DIR,
    repo_root: Path = REPO,
    git_root: Path = REPO,
    denylist_path: Path = DEFAULT_DENYLIST,
    baseline_ref: str = BASELINE_REF,
    u_values: bool = True,
) -> None:
    """All lineage checks on built data. `u_values` compares every non-retired U key's built VALUE with the baseline tag (fail-closed:
    an unreadable tag or file is a problem). Only synthetic-world unit tests that have no baseline tag turn it off."""
    check_scaffold(data_dir, manifest_path, run_ids_path, ms_dir)
    manifest = load_manifest(manifest_path)
    runs = load_run_ids(run_ids_path)
    min_code_commit = load_min_code_commit(run_ids_path)
    stage_min = load_stage_min_code_commit(run_ids_path)
    denylist = load_denylist(denylist_path)
    reuse_paths = frozenset(a["path"] for a in manifest.get("reuse_allowlist") or [])
    reuse_hashes = {a["path"]: a["sha256"] for a in manifest.get("reuse_allowlist") or []}
    problems = []
    for r in list(manifest["rows"]) + stage_source_rows(manifest):
        k = r["data_key"]
        if r["class"] in AFFECTED:
            stages = _stages_of(r)
            proto = r["expected_protocol"]
            if not all(s.startswith("B") for s in stages) and not all(s.startswith("X") for s in stages):
                problems.append(f"{k}: expected_run {stages} mixes breast (B) and Xenium (X) stages")
                continue
            try:
                names = set()
                for s in stages:
                    names.add(resolve_run(proto if s.startswith("B") else _x_proto(runs), s, runs))
            except LineageError as e:
                problems += [f"{k}: {x}" for x in e.problems]
                continue
            key = proto if proto != "none" else _x_proto(runs)
            rel = r["source"].format(protocol=key, **runs[key])
            src = Path(repo_root) / rel
            pp = path_problems(repo_root, rel, check_tokens=True)
            problems += [f"{k}: source {x}" for x in pp]
            if not src.is_file():
                problems.append(f"{k}: source {rel} not on disk")
                continue
            if any("is a symlink" in x for x in pp):
                continue
            if r.get("source_sha256") != sha256_file(src):
                problems.append(f"{k}: source_sha256 does not match {rel}")
            if stages == [HELD_OUT_STAGE]:
                problems += [f"{k}: {x}" for x in held_out_verdict_problems(src)]
            stamp = load_stamp(src)
            sp = validate_stamp(stamp, names, proto if proto != "none" else "none")
            problems += [f"{k}: {x}" for x in sp]
            if not sp:
                own_stage = next(st for st in stages if runs[key][st] == stamp["run_id"])
                problems += [f"{k}: {x}" for x in stamp_closure_problems(stamp, src, repo_root, git_root, reuse_paths, set(runs[key].values()), own_stage, min_code_commit, baseline_ref, denylist, stage_min_code_commit=stage_min, reuse_hashes=reuse_hashes)]
        elif r.get("source"):
            src = Path(repo_root) / r["source"]
            if not src.is_file() or sha256_file(src) != r.get("source_sha256"):
                problems.append(f"{k}: reused source {r['source']} missing or hash changed")
    problems += _check_built_labelset(Path(data_dir))
    problems += _check_manual_keys(manifest, Path(data_dir))
    problems += _check_retired(manifest, templates_of(ms_dir))
    if u_values:  # source-less U rows (refs::*, derived constants) carry no hash to check; their VALUES must still equal the baseline
        problems += u_value_problems(Path(data_dir), manifest, baseline_ref, git_root)
    if problems:
        raise LineageError(problems)


# ---------------------------------------------------------------- is this one source covered by a verified v42 stamp?
# A stamped MANIFEST that lists the sha256 of the files it covers. File name -> the key holding {repo-relative path: sha256}.
COVERING_MANIFESTS = {"layers_manifest.json": "layers"}
# A stamped result that lists the sha256 of SIBLING files in its own directory, keyed by bare file name (e.g. the ORA Test B
# results, which pin variants.csv). It covers only files in the same directory as itself.
SIBLING_MANIFESTS = {"ora_test_b_results.json": "sibling_files_sha256"}
_STAMP_VERDICTS: dict = {}


def verify_source_stamp(
    src: Path,
    repo_root: Path = REPO,
    git_root: Path | None = None,
    run_ids_path: Path = DEFAULT_RUN_IDS,
    manifest_path: Path = DEFAULT_MANIFEST,
    denylist_path: Path = DEFAULT_DENYLIST,
    baseline_ref: str = BASELINE_REF,
) -> tuple:
    """(ok, message): does `src` carry a v42 stamp that passes the SAME checks ``check_complete`` runs on a manifest source?

    The stamp must exist (sidecar ``<src>.stamp.json`` or embedded ``_stamp``), its ``run_id`` must be a current run id of the run-ids
    file (a superseded or unknown id is rejected), and ``validate_stamp`` + ``stamp_closure_problems`` must find nothing: output hash
    recomputed, code commit / producer blob / ancestry, every input hash recomputed, denylist, v2 upstream, no superseded path.
    A file that merely sits next to a ``.stamp.json`` of another file does not pass: the sidecar is looked up by exact name."""
    src = Path(src)
    root = Path(repo_root)
    git_root = Path(git_root) if git_root is not None else root
    if not src.is_file():
        return False, "not a file"
    try:
        stat = (src.stat().st_mtime_ns, src.stat().st_size, Path(str(src) + ".stamp.json").exists() and Path(str(src) + ".stamp.json").stat().st_mtime_ns)
    except OSError as e:
        return False, f"cannot stat: {e}"
    ck = (str(src), str(root), str(git_root), str(run_ids_path), str(manifest_path), stat)
    if ck in _STAMP_VERDICTS:
        return _STAMP_VERDICTS[ck]
    verdict = _verify_source_stamp_uncached(src, root, git_root, run_ids_path, manifest_path, denylist_path, baseline_ref)
    _STAMP_VERDICTS[ck] = verdict
    return verdict


def _verify_source_stamp_uncached(src, root, git_root, run_ids_path, manifest_path, denylist_path, baseline_ref) -> tuple:
    try:
        stamp = load_stamp(src)
    except (ValueError, OSError) as e:
        return False, f"unreadable stamp: {e}"
    if stamp is None:
        return False, "no v42 stamp (no sidecar .stamp.json, no embedded _stamp)"
    if not isinstance(stamp, dict):
        return False, "stamp is not an object"
    runs = load_run_ids(run_ids_path)
    rid = stamp.get("run_id")
    hits = [(k, st) for k, m in runs.items() for st, v in m.items() if v == rid]
    if len(hits) != 1:
        return False, f"run_id {rid!r} is not a current run id of {Path(run_ids_path).name} (superseded or unknown)"
    key, stage = hits[0]
    proto = key if stage.startswith("B") else "none"
    sp = validate_stamp(stamp, {rid}, proto)
    if sp:
        return False, "; ".join(sp)
    try:
        rel = str(Path(os.path.abspath(src)).relative_to(Path(os.path.abspath(root))))
    except ValueError:
        rel = None
    problems = []
    if rel is not None:
        problems += [f"source {x}" for x in path_problems(root, rel)]
    man = load_manifest(manifest_path)
    problems += stamp_closure_problems(
        stamp, src, root, git_root, frozenset(a["path"] for a in man.get("reuse_allowlist") or []), set(runs[key].values()), stage,
        load_min_code_commit(run_ids_path), baseline_ref, load_denylist(denylist_path),
        stage_min_code_commit=load_stage_min_code_commit(run_ids_path),
        reuse_hashes={a["path"]: a["sha256"] for a in man.get("reuse_allowlist") or []},
    )
    return (False, "; ".join(problems)) if problems else (True, f"stamped by run {rid} ({stage})")


def v42_stamp_covers(path: Path, repo_root: Path = REPO, verify=None) -> tuple:
    """(ok, message): is `path` covered by a v42 stamp that actually verifies?

    1. The file has its own stamp and ``verify_source_stamp`` accepts it, or
    2. a stamped manifest in one of its parent directories (COVERING_MANIFESTS, e.g. the X3 ``layers_manifest.json``) lists the file's
       repo-relative path with the file's CURRENT sha256, the manifest names the same run id as its stamp, and the manifest's own stamp
       verifies (case 1). A pinned summary does not cover files it does not list.
    Symlinks are never covered. ``verify`` is injectable for tests; the default is ``verify_source_stamp`` on ``repo_root``."""
    path = Path(path)
    root = Path(repo_root)
    verify = verify or (lambda p: verify_source_stamp(p, repo_root=root))
    if path.is_symlink():
        return False, "symlink"
    if not path.is_file():
        return False, "missing"
    own = verify(path)
    if own[0]:
        return own
    reasons = [f"own stamp: {own[1]}"]
    try:
        rel = str(Path(os.path.abspath(path)).relative_to(Path(os.path.abspath(root))))
    except ValueError:
        return False, "; ".join(reasons) + "; outside the repo root"
    here = Path(os.path.abspath(path)).parent
    for parent in Path(os.path.abspath(path)).parents:
        if parent == Path(os.path.abspath(root)).parent:
            break
        candidates = [(name, field, rel) for name, field in COVERING_MANIFESTS.items()]
        if parent == here:
            candidates += [(name, field, path.name) for name, field in SIBLING_MANIFESTS.items()]
        for name, field, key in candidates:
            m = parent / name
            if not m.is_file() or m == Path(os.path.abspath(path)):
                continue
            try:
                d = json.loads(m.read_text())
            except ValueError:
                reasons.append(f"{name} is not JSON")
                continue
            listed = d.get(field) if isinstance(d, dict) else None
            if not isinstance(listed, dict) or key not in listed:
                continue
            if listed[key] != sha256_file(path):
                return False, "; ".join(reasons) + f"; sha256 differs from the one {name} lists"
            ok, msg = verify(m)
            if not ok:
                return False, "; ".join(reasons) + f"; {name} stamp: {msg}"
            st = load_stamp(m)
            if not isinstance(d.get("run_id"), str) or d.get("run_id") != st.get("run_id"):
                return False, "; ".join(reasons) + f"; {name} run_id {d.get('run_id')!r} != its stamp's run_id {st.get('run_id')!r}"
            return True, f"listed in {name} (sha256 matches), whose stamp verifies"
    return False, "; ".join(reasons) + "; no stamped manifest lists it"


def _x_proto(runs: dict) -> str:
    """The protocol key under which the X stages are keyed (the same single key as the B stages)."""
    return require_single_key(runs)


def require_single_key(runs: dict) -> str:
    if len(runs) != 1:
        raise _err(UnknownProtocol, f"expected exactly one protocol key, found {sorted(runs)}")
    return next(iter(runs))


def _built_json_files(data_dir: Path):
    for p in sorted(data_dir.rglob("*.json")):
        rel = p.relative_to(data_dir)
        if "figures" in rel.parts or NON_CITABLE_FILE_RE.match(rel.stem):
            continue
        yield _file_id(data_dir, p), p


def _check_built_labelset(data_dir: Path) -> list:
    found, problems = {}, []
    for fid, p in _built_json_files(data_dir):
        d = json.loads(p.read_text())
        if isinstance(d, dict) and "_generated_at" in d and "_sources" in d:
            found.setdefault(d.get("_labelset"), []).append(p.name)
        elif fid not in UNSTAMPED_DATA_FILES:
            problems.append(f"{fid}.json lacks _generated_at/_sources and is not in UNSTAMPED_DATA_FILES")
    if set(found) != {"v42"}:
        problems.append(f"built data must record _labelset v42 only, found { {k: len(v) for k, v in found.items()} }")
    return problems


def _check_manual_keys(manifest: dict, data_dir: Path) -> list:
    """A P, G or BOTH key may not be fed by a hand-entered constant (`_manual_keys` of a built file).

    Entries look like "<group>.<key>[.<sub>]"; the top-level data key is the second component."""
    cls = {r["data_key"]: r["class"] for r in manifest["rows"]}
    out = []
    for fid, p in _built_json_files(data_dir):
        d = json.loads(p.read_text())
        for entry in (d.get("_manual_keys") or []) if isinstance(d, dict) else []:
            parts = str(entry).split(".")
            top = parts[1] if len(parts) > 1 else parts[0]
            if cls.get(f"{fid}::{top}") in AFFECTED:
                out.append(f"{fid}::{top}: class {cls[f'{fid}::{top}']} key is built from a hand-entered constant (_manual_keys {entry!r}); compute it from v2 output")
    return out


def _check_retired(manifest: dict, templates: list) -> list:
    retired = {r["data_key"] for r in manifest["rows"] if r.get("retired")}
    ids = {k.split("::")[0] for k in retired}
    keys_by_file: dict = {}
    for r in manifest["rows"]:
        keys_by_file.setdefault(r["data_key"].split("::")[0], set()).add(r["data_key"])
    out = []
    for t in templates:
        text = Path(t).read_text()
        _, refs = template_refs(text, ids, keys_by_file)
        out += [f"{Path(t).name}: cites retired key {k}" for k in sorted(set(refs) & retired)]
        for fig in manifest.get("retired_figures") or []:
            if fig["file"] in text:
                out.append(f"{Path(t).name}: still embeds retired figure {fig['file']}")
    return out


# ---------------------------------------------------------------- U keys must not move
def u_value_drift(manifest: dict, new: dict, base: dict) -> list:
    """Every non-retired U key must have the same value in `new` as in `base` ({file id: parsed JSON dict}).

    A U key that moves on a v2 rebuild was misclassified (it depends on Phase 1). Task 11 calls this with the v42 build
    against v41-final; see check_u_unchanged for the loader."""
    rows = [r for r in manifest.get("rows") or [] if r["class"] == "U" and not r.get("retired")]
    if not rows:
        return ["no U rows to compare (vacuous u-drift check)"]
    out = []
    for r in rows:
        fid, key = r["data_key"].split("::")
        b, n = base.get(fid), new.get(fid)
        if "v42_expected" in r:  # a U key that v42 CHANGED by design (e.g. a figure panel relabel): pinned to the explicit new value, so further drift still fails
            if not isinstance(n, dict) or key not in n:
                out.append(f"{r['data_key']}: U key missing from the new build")
            elif n[key] != r["v42_expected"]:
                out.append(f"{r['data_key']}: U key is {n[key]!r:.80}, the manifest pins {r['v42_expected']!r:.80}")
            elif isinstance(b, dict) and key in b and b[key] == r["v42_expected"]:
                out.append(f"{r['data_key']}: v42_expected equals the baseline value (drop the pin so the key is compared)")
            continue
        if r.get("baseline") == "absent":  # a U key that v42 ADDS (public data, no Phase 1 input): it must be built, and v41-final must not have it
            if isinstance(b, dict) and key in b:
                out.append(f"{r['data_key']}: marked baseline absent but the baseline has it (drop the marker so it is compared)")
            elif not isinstance(n, dict) or key not in n:
                out.append(f"{r['data_key']}: U key missing from the new build")
            continue
        if not isinstance(b, dict) or key not in b:
            out.append(f"{r['data_key']}: not in the baseline (cannot confirm it is unchanged)")
        elif not isinstance(n, dict) or key not in n:
            out.append(f"{r['data_key']}: U key missing from the new build")
        elif n[key] != b[key]:
            out.append(f"{r['data_key']}: U key changed vs baseline ({b[key]!r:.60} -> {n[key]!r:.60}); it depends on v2 output, reclassify")
    return out


def u_value_problems(data_dir: Path, manifest: dict, baseline_ref: str = "v41-final", repo_root: Path = REPO) -> list:
    """Load the built data under `data_dir` and the manuscript data dir at git ref `baseline_ref`; return the problems for any U key that moved.

    Fail-closed: a baseline tag or file that cannot be read from git is a problem, never a skip."""
    fids = sorted({r["data_key"].split("::")[0] for r in manifest["rows"] if r["class"] == "U" and not r.get("retired")})
    new, base, problems = {}, {}, []
    for fid in fids:
        f = Path(data_dir) / f"{fid}.json"
        new[fid] = json.loads(f.read_text()) if f.exists() else None
        r = subprocess.run(["git", "-C", str(repo_root), "show", f"{baseline_ref}:{os.environ.get('CITEGEIST_MS_DATA_REL', 'data')}/{fid}.json"], capture_output=True, text=True)
        if r.returncode != 0:
            problems.append(f"{fid}: not in {baseline_ref} ({r.stderr.strip()[:80]})")
        else:
            base[fid] = json.loads(r.stdout)
    problems += u_value_drift(manifest, new, base)
    return problems


def check_u_unchanged(data_dir: Path, manifest_path=DEFAULT_MANIFEST, baseline_ref: str = "v41-final", repo_root: Path = REPO) -> None:
    """Raise if any non-retired U key of the build under `data_dir` differs from `baseline_ref`."""
    problems = u_value_problems(Path(data_dir), load_manifest(manifest_path), baseline_ref, repo_root)
    if problems:
        raise LineageError(problems)


# ---------------------------------------------------------------- supplementary list for v42
def supplementary_order_for_v42(manifest_json: dict, lineage: dict) -> list:
    """MANIFEST.json supplementary_order with the lineage manifest's retired figures dropped and ids renumbered S1..Sn.

    Applies to labelset v42 only; v40/v41 keep reading MANIFEST.json unchanged."""
    gone = {f["file"] for f in lineage.get("retired_figures") or []}
    kept = [e for e in manifest_json["supplementary_order"] if e["file"] not in gone]
    return [{**e, "id": f"S{i}"} for i, e in enumerate(kept, start=1)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--stage", choices=("scaffold", "complete"), required=True)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--run-ids", type=Path, default=DEFAULT_RUN_IDS)
    ap.add_argument("--ms-dir", type=Path, default=DEFAULT_MS_DIR)
    a = ap.parse_args(argv)
    fn = check_scaffold if a.stage == "scaffold" else check_complete
    try:
        fn(a.data_dir, a.manifest, a.run_ids, a.ms_dir)
    except LineageError as e:
        for p in e.problems:
            print(f"FAIL: {p}", file=sys.stderr)
        print(f"{a.stage}: {len(e.problems)} problem(s)", file=sys.stderr)
        return 1
    print(f"{a.stage}: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
