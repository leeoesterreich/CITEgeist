"""Fail-closed input resolution and checks for the v2 breast chain stages B3-B12 (v42 rollout, Task 5).

Why: the A12 wrappers and entry points carried v41 defaults (``c1``, ``m3``, ``r1``, ``n1``, ``ec0cf5b80715``,
``output/module3_7type``). A v42 run that forgot one variable would silently read v41. This module is the one place that

  * resolves every run id from ``v42_run_ids.yaml`` (never from a default) and refuses an explicit id that differs;
  * checks that each input directory is exactly ``<family>/<protocol_sha12>/<run id>`` (or the stamped B1/B2 layout), has no
    v41 / stand-in token or symlink component, and is complete (stamped B1 CSVs and B2 nuclei, conserved SACE layers, ...);
  * refuses an incomplete COMMOT set before the compartment hurdle can read it as complete;
  * records what was resolved in a receipt next to the submission.

Two callers:
  * shell wrappers: ``chain_preflight.py env --stage B5 --run-id v2b5`` prints ``export`` lines (all checks run, also under DRYRUN);
    ``chain_preflight.py receipt ...`` writes the receipt at submit time.
  * Python entry points: ``guard_entry(stage, protocol_path, **paths)`` at the top of ``main`` (job-start check; defence in depth
    for a job launched without the wrapper). A protocol that is not in the run-ids file (a v40/v41 protocol file passed
    explicitly, e.g. by legacy tests) is not checked: ``guard_entry`` returns None.

Stage map (scope map B1-B12): B1 module3 refit, B2 nuclei, B3 Hungarian rebuild, B4 SACE layers, B5 compartments, B6 COMMOT,
B7 positive control, B8 hurdle, B9 Moran's I (+trajectory), B10 NMF, B11 program recovery, B12 keratin-secretory.
B17 matched gene controls and B18 MDK-NCL cell-type map are the promoted exploratory checks 2 and 3 (they read B5+B9 and B6+B8; control_pairs.py, celltype_map.py).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE / "src"))

import check_v42_lineage as L  # noqa: E402
import nuclei_artifact as NA  # noqa: E402
import protocol_io as PIO  # noqa: E402
import v42_stamp as W  # noqa: E402

DEFAULT_PROTOCOL = HERE / "protocol.json"
RESULTS = "Benchmarking/midkine_remediation/results"
HEADLINE = "HCC22-088-P4-S2_1i_rep"
DEVIATION_FILE = "PROTOCOL_DEVIATION_2026-11-06_rna_total.md"   # = hurdle.DEVIATION_FILE (a literal here so this module does not import the analysis code; tested equal)
MODULE3_PARENT_RE = re.compile(r"^module3_rna_refined_7type_(\d{8}T\d{6}Z)$")
NUCLEI_PARENT_RE = re.compile(r"^nuclei_rebuilt_7type_(\d{8}T\d{6}Z)$")
V41_DEFAULT_IDS = frozenset({"c1", "m3", "r1", "n1", "pc1", "h3d", "m1", "ks1", "pr4", "r3", "r4", "r5"})

# stage -> (output family relative to the namespace root, upstream stages whose outputs it reads)
FAMILY = {
    "B3": "output/morphology_assignment_v4_hungarian",
    "B4": "output/sace_layers_module3_7type",
    "B5": "output/compartments_module3_7type",
    "B6": "output/commot_runs_module3_7type",
    "B7": "output/positive_control_module3_7type",
    "B8": RESULTS,
    "B9": RESULTS,
    "B10": "output/nmf_counts_module3_7type",
    "B11": RESULTS,
    "B12": RESULTS,
    "B17": RESULTS,
    "B18": RESULTS,
    "B19": "Benchmarking/ora_check/runs",   # not Benchmarking/ora_check/results: that folder is DVC-protected
    "B20": "Benchmarking/ora_check/runs",
    "B21": RESULTS,
    "B22": RESULTS,
    "B221": RESULTS,   # first rerun of B22 (adds check 17)
    "B23": RESULTS,
    "B24": RESULTS,   # coverage of the nominal GMR interval under rna="total" (prod_coverage.py; a stamped input of B8)
}
# input key -> (producing stage, env var that may carry an explicit id for it)
INPUT_STAGE = {"layers_run": "B4", "comp_run": "B5", "commot_root": "B6", "v4_run": "B3", "nmf_run": "B10", "pc_run": "B7",
               "morans_run": "B9", "hurdle_run": "B8", "controls_run": "B17", "keratin_run": "B12", "e2_run": "B21", "coverage_run": "B24"}
ID_ENV = {"B24": "COVERAGE_RUN_ID", "B17": "CONTROLS_RUN_ID", "B12": "KERATIN_RUN_ID", "B21": "E2_RUN_ID", "B8": "HURDLE_RUN_ID", "B9": "MORANS_RUN_ID", "B3": "V4_RUN_ID", "B4": "LAYERS_RUN_ID", "B5": "COMP_RUN_ID", "B6": "COMMOT_RUN_ID", "B7": "POSITIVE_CONTROL_RUN_ID", "B10": "NMF_RUN_ID"}
REQUIRED = {
    "B3": ("module3_dir", "nuclei_dir"),
    "B4": ("module3_dir", "nuclei_dir"),
    "B5": ("module3_dir", "layers_run"),
    "B6": ("module3_dir", "layers_run", "comp_run"),
    "B7": ("module3_dir", "layers_run", "comp_run"),
    "B8": ("comp_run", "commot_root", "positive_control", "coverage_run"),
    "B9": ("comp_run", "v4_run"),
    "B10": ("module3_dir", "layers_run"),
    "B11": ("module3_dir", "layers_run", "nmf_run"),
    "B12": ("comp_run",),
    "B17": ("comp_run", "morans_run"),
    "B18": ("commot_root", "hurdle_run"),
    "B19": ("module3_dir", "nuclei_dir"),
    "B20": ("module3_dir", "nuclei_dir"),
    "B21": ("comp_run", "commot_root", "hurdle_run"),
    "B22": ("comp_run", "morans_run", "controls_run", "keratin_run", "e2_run"),
    "B221": ("comp_run", "morans_run", "controls_run", "keratin_run", "e2_run"),
    "B23": ("comp_run", "keratin_run"),
    "B24": (),
}


ENV_NAME = {"module3_dir": "MODULE3_DIR", "nuclei_dir": "NUCLEI_DIR", "layers_run": "LAYERS_ROOT", "comp_run": "COMP_RUN / COMPARTMENTS_RUN",
            "commot_root": "COMMOT_ROOT", "v4_run": "V4_ROOT", "nmf_run": "NMF_RUN", "positive_control": "POSITIVE_CONTROL",
            "morans_run": "MORANS_RUN", "hurdle_run": "HURDLE_RUN", "controls_run": "CONTROLS_RUN", "keratin_run": "KERATIN_RUN", "e2_run": "E2_RUN", "coverage_run": "COVERAGE_RUN"}


class ChainPreflightError(RuntimeError):
    """A v42 stage input is missing, defaulted, from the wrong run, or incomplete."""


class CommotIncomplete(ChainPreflightError):
    """The COMMOT set has a missing, staged-only or infrastructure-failed specimen."""


# ------------------------------------------------------------------ ids and protocol
def protocol_sha12(protocol_path=None) -> str:
    return PIO.load_protocol(protocol_path or DEFAULT_PROTOCOL, allow_draft=True)["_protocol_sha12"]


def is_v42(sha12: str, runs=None) -> bool:
    return sha12 in (runs if runs is not None else L.load_run_ids())


def run_id_for(stage: str, sha12: str, runs=None) -> str:
    try:
        return L.resolve_run(sha12, stage, runs if runs is not None else L.load_run_ids())
    except L.LineageError as e:
        raise ChainPreflightError(f"cannot resolve the {stage} run id for protocol {sha12}: {e}") from e


def check_explicit_id(stage: str, given, sha12: str, runs=None) -> str:
    """The yaml id for ``stage``; an explicit ``given`` must equal it (a v41 default such as c1/m3/r1/n1 is named in the message)."""
    want = run_id_for(stage, sha12, runs)
    if given in (None, ""):
        return want
    if given != want:
        hint = " (a v41 default name)" if given in V41_DEFAULT_IDS else ""
        if given in L.load_generation()["superseded_run_ids"]:
            hint = " (a SUPERSEDED v42 run: built from the counts-only B2 / the in-job re-segmentation; its folder stays on disk but is not an input)"
        raise ChainPreflightError(f"{stage}: run id {given!r}{hint} is not the v42 id {want!r} in v42_run_ids.yaml for protocol {sha12}")
    return want


def _ns(repo_root=None) -> Path:
    if repo_root is not None:
        return Path(repo_root)
    try:
        return W.resolve_roots(None)[0]
    except W.StampError as e:
        raise ChainPreflightError(str(e)) from e


def _rel(ns: Path, p) -> str:
    n = Path(os.path.normpath(ns))
    q = Path(os.path.normpath(p if Path(p).is_absolute() else n / p))
    try:
        return str(q.relative_to(n))
    except ValueError:
        pass
    try:
        return str(q.resolve().relative_to(n.resolve()))
    except ValueError:
        raise ChainPreflightError(f"{p} is outside the namespace root {ns}")


def _abs(ns: Path, p) -> Path:
    return Path(os.path.normpath(p if Path(p).is_absolute() else Path(ns) / p))


def _expect_path(ns: Path, what: str, given, expected_rel: str) -> Path:
    exp = Path(os.path.normpath(Path(ns) / expected_rel))
    got = _abs(ns, given)
    if got != exp:
        raise ChainPreflightError(f"{what}: {got} is not the v42 location {exp}")
    probs = L.path_problems(Path(ns), expected_rel)
    if probs:
        raise ChainPreflightError(f"{what}: " + "; ".join(probs))
    if not exp.is_dir():
        raise ChainPreflightError(f"{what}: {exp} does not exist (has the upstream stage finished?)")
    return exp


def stage_dir_rel(stage: str, sha12: str, run_id: str) -> str:
    return f"{FAMILY[stage]}/{sha12}/{run_id}"


def _specimens(protocol) -> list:
    return list(protocol["specimens"]["ordered"])


def _sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------ B1 / B2 inputs
def check_generation_ts(ts: str, what: str, gen=None) -> None:
    """The V2_TS of a B1/B2 dir must be the current one in v42_run_ids.yaml; a superseded one is named as such."""
    gen = gen if gen is not None else L.load_generation()
    if ts in gen["superseded_v2_ts"]:
        raise ChainPreflightError(f"{what}: V2_TS {ts} is SUPERSEDED (v42_run_ids.yaml): its B2 is counts-only and B3-B5 built from it were superseded; use V2_TS {gen['current_v2_ts']}")
    if gen["current_v2_ts"] is not None and ts != str(gen["current_v2_ts"]):
        raise ChainPreflightError(f"{what}: V2_TS {ts} is not the current V2_TS {gen['current_v2_ts']} (v42_run_ids.yaml)")


def check_module3_dir(module3_dir, sha12: str, specimens, ns=None, runs=None) -> dict:
    """B1 output dir: ``output/module3_rna_refined_7type_<ts>/<B1 id>``, one stamped csv per specimen, read by exact name (never a glob)."""
    if module3_dir in (None, ""):
        raise ChainPreflightError("MODULE3_DIR is required for a v42 run (no default; the old default output/module3_7type is v40)")
    ns = _ns(ns)
    b1 = run_id_for("B1", sha12, runs)
    d = _abs(ns, module3_dir)
    if d.name != b1 or not MODULE3_PARENT_RE.match(d.parent.name):
        raise ChainPreflightError(f"MODULE3_DIR {d} must be output/module3_rna_refined_7type_<UTC ts>/{b1} (the stamped B1 dir); a v41 or p0 dir is refused")
    check_generation_ts(MODULE3_PARENT_RE.match(d.parent.name).group(1), "MODULE3_DIR")
    rel = _rel(ns, d)
    if str(Path(rel).parent.parent) != "output":
        raise ChainPreflightError(f"MODULE3_DIR {d} must sit directly under {ns}/output")
    probs = L.path_problems(ns, rel)
    if probs:
        raise ChainPreflightError("MODULE3_DIR: " + "; ".join(probs))
    if not d.is_dir():
        raise ChainPreflightError(f"MODULE3_DIR {d} does not exist")
    out = {}
    for s in specimens:
        csv = d / f"{s}_cell_prop_finetuned_results.csv"
        if not csv.is_file() or csv.is_symlink():
            raise ChainPreflightError(f"MODULE3_DIR: {csv.name} is missing or a symlink in {d}")
        try:
            W.verify_stamp(csv, run_id=b1, protocol_sha12=sha12, repo_root=ns)
        except W.StampError as e:
            raise ChainPreflightError(f"MODULE3_DIR: {csv.name} is not a valid stamped B1 output: {e}") from e
        out[s] = {"csv_sha256": _sha256(csv), "stamp_sha256": _sha256(str(csv) + ".stamp.json")}
    return out


def check_nuclei_dir(nuclei_dir, sha12: str, specimens, module3_dir=None, ns=None, runs=None) -> dict:
    """B2 output dir: ``output/nuclei_rebuilt_7type_<ts>/<B2 id>`` (same ts as MODULE3_DIR), stamped, never the v3 stand-in."""
    if nuclei_dir in (None, ""):
        raise ChainPreflightError("NUCLEI_DIR is required for a v42 run (no default; output/morphology_assignment_v3 is a v4 stand-in)")
    ns = _ns(ns)
    b2 = run_id_for("B2", sha12, runs)
    d = _abs(ns, nuclei_dir)
    try:
        NA.refuse_standin(d)
    except NA.StandInRefused as e:
        raise ChainPreflightError(str(e)) from e
    m = NUCLEI_PARENT_RE.match(d.parent.name)
    if d.name != b2 or not m:
        raise ChainPreflightError(f"NUCLEI_DIR {d} must be output/nuclei_rebuilt_7type_<UTC ts>/{b2} (the stamped B2 dir)")
    check_generation_ts(m.group(1), "NUCLEI_DIR")
    if module3_dir not in (None, ""):
        mm = MODULE3_PARENT_RE.match(_abs(ns, module3_dir).parent.name)
        if mm and mm.group(1) != m.group(1):
            raise ChainPreflightError(f"NUCLEI_DIR timestamp {m.group(1)} differs from MODULE3_DIR timestamp {mm.group(1)}: two different V2_TS runs")
    try:
        NA.refuse_standin(d)
    except NA.StandInRefused as e:
        raise ChainPreflightError(str(e)) from e
    probs = L.path_problems(ns, _rel(ns, d))
    if probs:
        raise ChainPreflightError("NUCLEI_DIR: " + "; ".join(probs))
    out = {}
    import pandas as pd

    for s in specimens:
        f = NA.nuclei_csv(d, s)
        if not f.is_file():
            raise ChainPreflightError(f"NUCLEI_DIR: {f.name} missing in {d}")
        p0_index = None
        if module3_dir not in (None, ""):
            csv = _abs(ns, module3_dir) / f"{s}_cell_prop_finetuned_results.csv"
            if csv.is_file():
                p0_index = pd.read_csv(csv, index_col=0).index
        if p0_index is None:
            p0_index = NA.read_counts(f).index
        try:  # counts + positions only together, both stamped for this run, counts == the table's row counts, no stand-in
            NA.load_nuclei_artifact(d, s, p0_index, protocol_sha12=sha12, run_id=b2, repo_root=ns, runs=runs)
        except (NA.NucleiCoverageError, NA.StandInRefused, ValueError) as e:
            raise ChainPreflightError(f"NUCLEI_DIR: {s}: {e}") from e
        t = NA.table_csv(d, s)
        out[s] = {"csv_sha256": _sha256(f), "stamp_sha256": _sha256(str(f) + ".stamp.json"),
                  "table_sha256": _sha256(t), "table_stamp_sha256": _sha256(str(t) + ".stamp.json")}
    return out


# ------------------------------------------------------------------ upstream run dirs
def check_layers_run(path, sha12, specimens, ns, runs=None) -> Path:
    d = _expect_path(ns, "layers_run (B4)", path, stage_dir_rel("B4", sha12, run_id_for("B4", sha12, runs)))
    for s in specimens:
        mf = d / s / "manifest.json"
        if not mf.is_file():
            raise ChainPreflightError(f"layers_run: {mf} missing (B4 incomplete for {s})")
        if not json.loads(mf.read_text()).get("conservation", {}).get("passed"):
            raise ChainPreflightError(f"layers_run: {s} did not pass the SACE conservation gate")
    return d


def check_comp_run(path, sha12, specimens, ns, runs=None) -> Path:
    d = _expect_path(ns, "comp_run (B5)", path, stage_dir_rel("B5", sha12, run_id_for("B5", sha12, runs)))
    for name in ["M.json", "a12_7_summary.json"] + [f"{s}.h5ad" for s in specimens]:
        if not (d / name).is_file():
            raise ChainPreflightError(f"comp_run: {d / name} missing (B5 incomplete)")
    return d


def check_v4_run(path, sha12, ns, protocol, runs=None) -> Path:
    d = _expect_path(ns, "v4_run (B3)", path, stage_dir_rel("B3", sha12, run_id_for("B3", sha12, runs)))
    gated = set(protocol["rebuild"]["no_rebuild_below_gate"]["specimens"])
    for s in _specimens(protocol):
        if s in gated:
            continue  # recorded not_estimable by Definition A; the below-gate route (BELOW_GATE=1) fills them for figures
        if not (d / s / f"{s}_single_cell.h5ad").is_file():
            raise ChainPreflightError(f"v4_run: {s}_single_cell.h5ad missing in {d} (B3 incomplete; the normal route covers the 10 in-gate specimens)")
    return d


def check_nmf_run(path, sha12, ns, runs=None) -> Path:
    d = _expect_path(ns, "nmf_run (B10)", path, stage_dir_rel("B10", sha12, run_id_for("B10", sha12, runs)))
    if not (d / "nmf_counts.json").is_file():
        raise ChainPreflightError(f"nmf_run: {d / 'nmf_counts.json'} missing (B10 incomplete)")
    return d


def check_positive_control(path, sha12, ns, runs=None) -> Path:
    exp = f"{stage_dir_rel('B7', sha12, run_id_for('B7', sha12, runs))}/{HEADLINE}/positive_control.json"
    got = _abs(ns, path)
    if got != Path(os.path.normpath(Path(ns) / exp)):
        raise ChainPreflightError(f"positive_control: {got} is not the B7 file {Path(ns) / exp}")
    probs = L.path_problems(Path(ns), exp)
    if probs:
        raise ChainPreflightError("positive_control: " + "; ".join(probs))
    if not got.is_file():
        raise ChainPreflightError(f"positive_control: {got} does not exist (B7 incomplete)")
    return got


def check_results_run(path, sha12, ns, stage: str, fname: str, runs=None) -> Path:
    """An upstream results/<sha12>/<id> dir (B8, B9) that a promoted stage reads a reference value from: the v42 location, the stamped file and its stamp present."""
    what = f"{stage} results"
    d = _expect_path(ns, what, path, stage_dir_rel(stage, sha12, run_id_for(stage, sha12, runs)))
    for name in (fname, fname + ".stamp.json"):
        if not (d / name).is_file():
            raise ChainPreflightError(f"{what}: {d / name} missing ({stage} incomplete or never stamped)")
    return d


def require_commot_complete(commot_root, specimens) -> dict:
    """Every specimen finalized (specimen.json), no staged ``.tmp_<spec>_*`` left, no run lost to the infrastructure
    (reason ``child-timeout`` / ``child-aborted-without-run.json`` / ``child-run.json-corrupt``), every listed run has its run.json.

    A typed scientific failure (e.g. ``commot-ot-nonconvergence``) is NOT incomplete: it is a result the hurdle reports as typed."""
    root = Path(commot_root)
    problems = []
    for s in specimens:
        d = root / s
        staged = sorted(p.name for p in root.glob(f".tmp_{s}_*") if re.match(rf"^\.tmp_{re.escape(s)}_\d+$", p.name)) if root.is_dir() else []
        if staged:
            problems.append(f"{s}: staged {staged} (not finalized; resume with RESUME=1)")
        sj = d / "specimen.json"
        if not sj.is_file():
            problems.append(f"{s}: no specimen.json (COMMOT not finalized)")
            continue
        runs = json.loads(sj.read_text()).get("runs") or []
        if not runs:
            problems.append(f"{s}: specimen.json lists no runs")
        for r in runs:
            rid = r.get("run_id")
            if not (d / str(rid) / "run.json").is_file():
                problems.append(f"{s}/{rid}: run.json missing")
            if str(r.get("reason") or "").startswith("child-"):
                problems.append(f"{s}/{rid}: lost to the infrastructure ({r.get('reason')}); archive the specimen dir and rerun")
    if problems:
        raise CommotIncomplete("COMMOT set is incomplete: " + "; ".join(problems))
    return {"commot_root": str(root), "specimens": list(specimens)}


# ------------------------------------------------------------------ the one check function
def check_inputs(stage: str, sha12: str, ns=None, protocol=None, runs=None, **paths) -> dict:
    """Check every input of ``stage``. ``paths`` keys: module3_dir, nuclei_dir, layers_run, comp_run, commot_root, v4_run, nmf_run,
    positive_control, v3_root (must be empty). Missing required keys, a wrong location or an incomplete upstream raises
    ChainPreflightError. Returns {"inputs": {...hashes...}, "resolved": {key: path}}."""
    if stage not in REQUIRED:
        raise ChainPreflightError(f"unknown stage {stage!r}; known: {sorted(REQUIRED)}")
    runs = runs if runs is not None else L.load_run_ids()
    ns = _ns(ns)
    if protocol is None:
        protocol = PIO.load_protocol(DEFAULT_PROTOCOL)
    if protocol["_protocol_sha12"] != sha12:
        raise ChainPreflightError(f"protocol sha12 {protocol['_protocol_sha12']} != {sha12}")
    specimens = _specimens(protocol)
    for k in REQUIRED[stage]:
        if paths.get(k) in (None, ""):
            raise ChainPreflightError(f"{stage}: {k} ({ENV_NAME[k]}) is required for a v42 run (no default; the old defaults are v40/v41)")
    for k, v in paths.items():  # no v2 stage reads output/morphology_assignment_v3 (the v4 stand-in), whichever input it arrives as
        if v not in (None, "") and k != "v3_root":
            try:
                NA.refuse_standin(_abs(ns, v))
            except NA.StandInRefused as e:
                raise ChainPreflightError(f"{k} ({ENV_NAME.get(k, k)}): {e}") from e
    if paths.get("v3_root") not in (None, ""):
        raise ChainPreflightError(f"{stage}: v3_root must stay unset (output/morphology_assignment_v3 is a v4 stand-in)")
    inputs, resolved = {}, {}
    if paths.get("module3_dir") not in (None, ""):
        inputs["module3"] = check_module3_dir(paths["module3_dir"], sha12, specimens, ns, runs)
        resolved["module3_dir"] = str(_abs(ns, paths["module3_dir"]))
    if paths.get("nuclei_dir") not in (None, ""):
        inputs["nuclei"] = check_nuclei_dir(paths["nuclei_dir"], sha12, specimens, paths.get("module3_dir"), ns, runs)
        resolved["nuclei_dir"] = str(_abs(ns, paths["nuclei_dir"]))
    if paths.get("layers_run") not in (None, ""):
        resolved["layers_run"] = str(check_layers_run(paths["layers_run"], sha12, specimens, ns, runs))
    if paths.get("comp_run") not in (None, ""):
        resolved["comp_run"] = str(check_comp_run(paths["comp_run"], sha12, specimens, ns, runs))
    if paths.get("v4_run") not in (None, ""):
        resolved["v4_run"] = str(check_v4_run(paths["v4_run"], sha12, ns, protocol, runs))
    if paths.get("nmf_run") not in (None, ""):
        resolved["nmf_run"] = str(check_nmf_run(paths["nmf_run"], sha12, ns, runs))
    if paths.get("positive_control") not in (None, ""):
        resolved["positive_control"] = str(check_positive_control(paths["positive_control"], sha12, ns, runs))
    if paths.get("morans_run") not in (None, ""):
        resolved["morans_run"] = str(check_results_run(paths["morans_run"], sha12, ns, "B9", STAMPED_ARTIFACT["B9"], runs))
    if paths.get("hurdle_run") not in (None, ""):
        resolved["hurdle_run"] = str(check_results_run(paths["hurdle_run"], sha12, ns, "B8", STAMPED_ARTIFACT["B8"], runs))
    for key, up in (("controls_run", "B17"), ("keratin_run", "B12"), ("e2_run", "B21"), ("coverage_run", "B24")):
        if paths.get(key) not in (None, ""):
            resolved[key] = str(check_results_run(paths[key], sha12, ns, up, STAMPED_ARTIFACT[up], runs))
    if paths.get("commot_root") not in (None, ""):
        d = _expect_path(ns, "commot_root (B6)", paths["commot_root"], stage_dir_rel("B6", sha12, run_id_for("B6", sha12, runs)))
        resolved["commot_root"] = str(d)
        if stage in ("B8", "B18", "B21"):
            require_commot_complete(d, protocol["specimens"]["hurdle_specimens_all"])
    return {"inputs": inputs, "resolved": resolved}


def _v42_protocol(protocol_path):
    """(protocol, sha12, runs) when ``protocol_path`` is a v42 protocol, else None (legacy: a v40/v41 protocol or a toy protocol of a unit test).

    In a CHAIN RUN (``V2_NAMESPACE_ROOT`` set: every v2 wrapper exports it) None is not allowed: a protocol that is not hash-stamped, or whose
    sha12 is not a key of ``v42_run_ids.yaml``, raises. Without this a mis-set ``--protocol`` would run with no check at all and the legacy
    (v40/v41) defaults."""
    chain = bool(os.environ.get(W.NAMESPACE_ENV))
    try:
        raw = json.loads(Path(protocol_path).read_text())
    except (OSError, ValueError):
        raw = None
    if not isinstance(raw, dict) or "frozen_sha256" not in raw:
        if chain:
            raise ChainPreflightError(f"{W.NAMESPACE_ENV} is set (a v2 chain run) but {protocol_path} is not a hash-stamped protocol file")
        return None  # not a hash-stamped protocol (a toy/synthetic protocol in a unit test): it cannot be a v42 protocol, so there is nothing to check
    proto = PIO.load_protocol(protocol_path, allow_draft=True)  # only the sha12 and the specimen lists are read here; the producer verifies the frozen hash itself
    sha12 = proto["_protocol_sha12"]
    runs = L.load_run_ids()
    if sha12 not in runs:
        if chain:
            raise ChainPreflightError(
                f"{W.NAMESPACE_ENV} is set (a v2 chain run) but protocol {sha12} ({protocol_path}) is not in v42_run_ids.yaml: refusing to run it unchecked")
        return None
    return proto, sha12, runs


def guard_entry(stage: str, protocol_path, *, run_id=None, out_root=None, repo_root=None, **paths):
    """Job-start check for a Python entry point. Returns None for a protocol that is not a v42 protocol (legacy behaviour unchanged) unless this
    is a chain run (``V2_NAMESPACE_ROOT`` set), where an unknown protocol raises; otherwise checks ``run_id``/``out_root`` against the yaml and
    every input path, raising ChainPreflightError."""
    got = _v42_protocol(protocol_path)
    if got is None:
        return None
    proto, sha12, runs = got
    ns = _ns(repo_root)
    if run_id is not None:
        check_explicit_id(stage, run_id, sha12, runs)
    if out_root is not None:
        _expect_path_parent(ns, stage, out_root, sha12, runs)
    return check_inputs(stage, sha12, ns=ns, protocol=proto, runs=runs, **paths)


def _expect_path_parent(ns, stage, out_root, sha12, runs):
    exp = Path(os.path.normpath(Path(ns) / stage_dir_rel(stage, sha12, run_id_for(stage, sha12, runs))))
    got = _abs(ns, out_root)
    if got != exp:
        raise ChainPreflightError(f"{stage}: out-root {got} is not the v42 location {exp}")


def legacy_or_required(value, protocol_path, what: str, legacy):
    """Argparse default helper: ``value`` if given; else ``legacy`` for a non-v42 protocol; else raise (no v42 default)."""
    if value is not None:
        return value
    if is_v42(protocol_sha12(protocol_path)):
        raise ChainPreflightError(f"{what} is required for a v42 run (no default; the legacy default is a v40/v41 path)")
    return legacy


# ------------------------------------------------------------------ stamps for the lineage sources (B8-B12)
# The one source file each of these stages hands to the v42 lineage manifest (rows with ``expected_run`` B8..B12).
STAMPED_ARTIFACT = {
    "B8": "commot_compartment_hurdle.json", "B9": "morans_results.json", "B10": "nmf_counts.json",
    "B11": "program_recovery.json", "B12": "keratin_secretory_summary.json",
    "B17": "control_pairs_results.json", "B18": "celltype_map.json",
    "B19": "ora_test_b_results.json", "B20": "ora_sace_results.json",
    "B21": "e2_covariation_results.json", "B22": "moran_robustness_results.json", "B221": "moran_robustness_results.json",
    "B23": "keratin_covariate_results.json", "B24": "prod_coverage_total.json",
}
# input keys each stage reads (same names as ``check_inputs``); every file read from them is hashed into the stamp
STAMP_INPUT_KEYS = {
    "B8": ("comp_run", "commot_root", "positive_control", "coverage_run"), "B9": ("comp_run", "v4_run"), "B10": ("module3_dir", "layers_run"),
    "B11": ("module3_dir", "layers_run", "nmf_run"), "B12": ("comp_run",),
    "B17": ("comp_run", "morans_run"), "B18": ("commot_root", "hurdle_run"),
    "B19": ("module3_dir", "nuclei_dir"), "B20": ("module3_dir", "nuclei_dir"),
    "B21": ("comp_run", "commot_root", "hurdle_run"), "B22": ("comp_run", "morans_run", "controls_run", "keratin_run", "e2_run"),
    "B221": ("comp_run", "morans_run", "controls_run", "keratin_run", "e2_run"),
    "B23": ("comp_run", "keratin_run"),
    "B24": (),   # no upstream stage: the inputs of the stamp are the protocol and the deviation file (below)
}


def stamp_input_files(stage: str, protocol, ns, **dirs) -> list:
    """The files a stage's stamp pins as inputs (namespace-root paths, all must exist; never a glob): per stage and key
    B1 csv / B4 ``manifest.json`` / B5 ``M.json``, ``a12_7_summary.json`` and ``<spec>.h5ad`` / B6 ``specimen.json`` + each ``run.json`` /
    B7 ``positive_control.json`` / B3 ``<spec>_single_cell.h5ad`` of the in-gate specimens / B10 ``nmf_counts.json``."""
    if stage not in STAMP_INPUT_KEYS:
        raise ChainPreflightError(f"stage {stage!r} writes no lineage stamp; stamped stages: {sorted(STAMP_INPUT_KEYS)}")
    specs = _specimens(protocol)
    gated = set(protocol["rebuild"]["no_rebuild_below_gate"]["specimens"])
    files = []
    for key in STAMP_INPUT_KEYS[stage]:
        v = dirs.get(key)
        if v in (None, ""):
            raise ChainPreflightError(f"{stage}: stamp input {key} was not given")
        d = _abs(ns, v)
        if key == "module3_dir":
            files += [d / f"{s}_cell_prop_finetuned_results.csv" for s in specs]
            if stage in ("B19", "B20"):  # the ORA stages also read the v2b1 p0 (the marker-share arm) of every specimen
                files += [d / "p0" / s / f"{s}_cell_prop_global_results.csv" for s in specs]
        elif key == "layers_run":
            files += [d / s / "manifest.json" for s in specs]
        elif key == "comp_run":
            files += [d / "M.json", d / "a12_7_summary.json"] + [d / f"{s}.h5ad" for s in specs]
        elif key == "commot_root":
            for s in protocol["specimens"]["hurdle_specimens_all"]:
                sj = d / s / "specimen.json"
                files.append(sj)
                if stage in ("B18", "B21"):  # the cell-type map (B18) and the keratin/E2 stage (B21) reads the R00 compartment rows and the sender scores themselves
                    files += [d / s / "rows_R00.csv.gz", d / s / "R00" / "sender.csv.gz"]
                if sj.is_file():
                    files += [d / s / str(r.get("run_id")) / "run.json" for r in json.loads(sj.read_text()).get("runs") or []]
        elif key == "positive_control":
            files.append(d)
        elif key == "v4_run":
            files += [d / s / f"{s}_single_cell.h5ad" for s in specs if s not in gated]
        elif key == "nmf_run":
            files.append(d / "nmf_counts.json")
        elif key == "morans_run":
            files.append(d / STAMPED_ARTIFACT["B9"])
        elif key == "hurdle_run":
            files.append(d / STAMPED_ARTIFACT["B8"])
            if stage == "B21":  # B21 gates the keratin-present rule against the B8 headline sweep file (A07 label group sizes), which the stamped B8 json does not hold
                files.append(d / HEADLINE / "sweep" / f"sweep_{HEADLINE}.json")
        elif key == "controls_run":
            files.append(d / STAMPED_ARTIFACT["B17"])
        elif key == "coverage_run":
            files.append(d / STAMPED_ARTIFACT["B24"])
        elif key == "keratin_run":
            files.append(d / STAMPED_ARTIFACT["B12"])
        elif key == "e2_run":
            files.append(d / STAMPED_ARTIFACT["B21"])
        elif key == "nuclei_dir":
            files += [d / s / f"{s}_nuclei_counts_fullres_patchwise.csv" for s in specs]
    if stage == "B24":  # synthetic worlds: no upstream data. The stamp pins the protocol and the deviation file the covariate mode is licensed by (B24 is a root stage of the lineage)
        files += [Path(ns) / "Benchmarking" / "midkine_remediation" / "protocol.json", Path(ns) / "Benchmarking" / "midkine_remediation" / DEVIATION_FILE]
    missing = [str(f) for f in files if not Path(f).is_file()]
    if missing:
        raise ChainPreflightError(f"{stage}: stamp inputs missing: {missing[:5]}{' ...' if len(missing) > 5 else ''}")
    return files


class StageStamp:
    """Job-start preflight + end-of-job stamp for a B8-B12 lineage source. ``stage_stamp_begin`` creates it (None for a legacy protocol)."""

    def __init__(self, stage, protocol, sha12, artifact, producer, run_id, ns, git, started):
        self.stage, self.protocol, self.sha12, self.artifact = stage, protocol, sha12, Path(artifact)
        self.producer, self.run_id, self.ns, self.git, self.started = producer, run_id, ns, git, started

    def finish(self, **dirs) -> dict:
        """Write ``<artifact>.stamp.json`` (inputs = every file ``stamp_input_files`` names). Raises StampError / ChainPreflightError."""
        files = stamp_input_files(self.stage, self.protocol, self.ns, **dirs)
        return W.write_stamp(
            self.artifact, run_id=self.run_id, protocol_sha12=self.sha12, inputs=files, producer=self.producer, started_at=self.started,
            repo_root=self.ns, git_root=self.git, stage_role=f"{self.stage}_lineage_source")


def stage_stamp_begin(stage: str, protocol_path, artifact, producer, run_id, *, repo_root=None, git_root=None):
    """Job-start gate, BEFORE any compute: SLURM job, namespace root, clean code tree, tracked+clean producer, artifact path inside the namespace
    root with no path problem (``v42_stamp.preflight``). Returns a StageStamp, or None for a legacy (non-v42) protocol."""
    if stage not in STAMPED_ARTIFACT:
        raise ChainPreflightError(f"stage {stage!r} writes no lineage stamp; stamped stages: {sorted(STAMPED_ARTIFACT)}")
    got = _v42_protocol(protocol_path)
    if got is None:
        return None
    proto, sha12, runs = got
    if Path(artifact).name != STAMPED_ARTIFACT[stage]:
        raise ChainPreflightError(f"{stage}: the lineage source is {STAMPED_ARTIFACT[stage]}, not {Path(artifact).name}")
    check_explicit_id(stage, run_id, sha12, runs)
    try:
        ns, git = W.preflight(producer, [artifact], repo_root=repo_root, git_root=git_root)
    except W.StampError as e:
        raise ChainPreflightError(f"{stage}: cannot stamp this run: {e}") from e
    return StageStamp(stage, proto, sha12, artifact, producer, run_id, ns, git, W.utc_now())


# ------------------------------------------------------------------ one code commit per stage run (I1)
CODE_SHA_RE = re.compile(r"^[0-9a-f]{40}$")
RECEIPTS_ENV = "V2_RECEIPTS_ROOT"  # test seam only: where pins are READ from (default <ns>/output/v2_receipts); tests point it at an empty tmp dir
PRODUCER_PATHS = ("Benchmarking/midkine_remediation", "scripts", "CITEgeist", ".commands")  # the same list submit-v2-breast-chain.sh checks


def pin_code_commit(stage: str, run_id: str, ns, head: str, *, dryrun: bool = False, protocol_path=None):
    """The one code commit every submission of (stage, run id) runs: a first submission pins ``head``; later ones reuse the pin even when ``head``
    moved (a first/rest array split, a resume). Returns (sha, how) with how in {"pinned", "new", "from-receipt"}.

    Pin file: ``<ns>/output/v2_receipts/<sha12>/<stage>_<run_id>.code_sha.txt`` (exclusive create, never overwritten). A run submitted before
    pins existed (v2b3-v2b5) has none: its EARLIEST receipt's ``code_commit`` is the pin, so a late submission into it cannot silently use a new commit."""
    sha12 = protocol_sha12(protocol_path)
    d = (Path(os.environ[RECEIPTS_ENV]) if os.environ.get(RECEIPTS_ENV) else Path(ns) / "output" / "v2_receipts") / sha12
    pin = d / f"{stage}_{run_id}.code_sha.txt"
    how = None
    if pin.is_file():
        sha, how = pin.read_text().strip(), "pinned"
    else:
        early = []
        if d.is_dir():
            early = sorted(f for f in (p.name for p in d.iterdir()) if re.match(rf"^{re.escape(stage)}_{re.escape(run_id)}_\d{{8}}T\d+Z\.json$", f))
        sha = None
        for f in early:
            try:
                c = json.loads((d / f).read_text()).get("code_commit")
            except (OSError, ValueError):
                c = None
            if c:
                sha, how = c, "from-receipt"
                break
        if sha is None:
            sha, how = head, "new"
    if not CODE_SHA_RE.match(sha or ""):
        raise ChainPreflightError(f"{stage} {run_id}: pinned code commit {sha!r} is not a 40-hex commit id ({pin})")
    r = subprocess.run(["git", "--no-optional-locks", "-C", str(ns), "cat-file", "-e", f"{sha}^{{commit}}"], capture_output=True, text=True)
    if r.returncode != 0:
        raise ChainPreflightError(f"{stage} {run_id}: pinned code commit {sha} is not a commit of this repository ({pin})")
    if sha != head:  # HEAD moved since the pin: harmless only if no producer path changed; otherwise a rerun would silently run the OLD code
        d_ = subprocess.run(["git", "--no-optional-locks", "-C", str(ns), "diff", "--quiet", sha, head, "--", *PRODUCER_PATHS], capture_output=True, text=True)
        if d_.returncode == 1:
            raise ChainPreflightError(
                f"{stage} {run_id}: the pinned code commit {sha[:12]} ({how}) differs from HEAD {head[:12]} in producer code "
                f"({', '.join(PRODUCER_PATHS)}). Resubmitting this run id would run the OLD code. Use a new run id for this stage "
                f"(add it to v42_run_ids.yaml) or, for the p0/nuclei/refit chain, a new V2_TS. The pin file {pin} is never edited.")
        if d_.returncode not in (0, 1):
            raise ChainPreflightError(f"{stage} {run_id}: cannot compare the pinned commit with HEAD: {d_.stderr.strip()[:200]}")
    if how != "pinned" and not dryrun:
        d.mkdir(parents=True, exist_ok=True)
        try:
            with open(pin, "x") as fh:
                fh.write(sha + "\n")
        except FileExistsError:  # a concurrent first submission pinned it between our read and our create: the pin wins
            return pin_code_commit(stage, run_id, ns, head, dryrun=dryrun, protocol_path=protocol_path)
    return sha, how


# ------------------------------------------------------------------ B3 route
def check_b3_route(specimens, below_gate: bool, protocol) -> None:
    """B3 normal route takes only in-gate specimens; BELOW_GATE=1 takes only the below-gate ones (never the bypass for a normal specimen)."""
    gated = set(protocol["rebuild"]["no_rebuild_below_gate"]["specimens"])
    specimens = list(specimens)
    if not specimens:
        raise ChainPreflightError("B3: no specimens given")
    unknown = [s for s in specimens if s not in _specimens(protocol)]
    if unknown:
        raise ChainPreflightError(f"B3: unknown specimens {unknown}")
    wrong = [s for s in specimens if (s in gated) != below_gate]
    if wrong:
        route = "BELOW_GATE=1 (only the below-gate specimens)" if not below_gate else "the normal route (BELOW_GATE unset)"
        raise ChainPreflightError(f"B3: {wrong} belong to the other route; this submission is {'BELOW_GATE=1' if below_gate else 'normal'}, use {route} for them")


# ------------------------------------------------------------------ env / receipt (shell side)
UPSTREAM_KEYS = {
    "B3": {}, "B4": {},
    "B5": {"layers_run": "B4"},
    "B6": {"layers_run": "B4", "comp_run": "B5"},
    "B7": {"layers_run": "B4", "comp_run": "B5"},
    "B8": {"comp_run": "B5", "commot_root": "B6", "pc": "B7", "coverage_run": "B24"},
    "B9": {"comp_run": "B5", "v4_run": "B3"},
    "B10": {"layers_run": "B4"},
    "B11": {"layers_run": "B4", "nmf_run": "B10"},
    "B12": {"comp_run": "B5"},
    "B17": {"comp_run": "B5", "morans_run": "B9"},
    "B18": {"commot_root": "B6", "hurdle_run": "B8"},
    "B19": {}, "B20": {},
    "B21": {"comp_run": "B5", "commot_root": "B6", "hurdle_run": "B8"},
    "B22": {"comp_run": "B5", "morans_run": "B9", "controls_run": "B17", "keratin_run": "B12", "e2_run": "B21"},
    "B221": {"comp_run": "B5", "morans_run": "B9", "controls_run": "B17", "keratin_run": "B12", "e2_run": "B21"},
    "B23": {"comp_run": "B5", "keratin_run": "B12"},
    "B24": {},
}


def resolve_env(stage: str, run_id, environ=None, protocol_path=None, ns=None) -> dict:
    """Resolve and check everything a wrapper needs for ``stage``. Returns {"exports": {VAR: value}, "report": {...}}."""
    env = os.environ if environ is None else environ
    proto = PIO.load_protocol(protocol_path or DEFAULT_PROTOCOL)
    sha12 = proto["_protocol_sha12"]
    runs = L.load_run_ids()
    if sha12 not in runs:
        raise ChainPreflightError(f"protocol {sha12} is not in v42_run_ids.yaml: the v2 wrappers run only for a v42 protocol")
    ns = _ns(ns)
    own = check_explicit_id(stage, run_id, sha12, runs)
    if run_id in (None, ""):
        raise ChainPreflightError(f"{stage}: RUN_ID must be given explicitly (expected {own})")
    ids = {stage: own}
    paths = {"module3_dir": env.get("MODULE3_DIR"), "nuclei_dir": env.get("NUCLEI_DIR")}
    for key, up in UPSTREAM_KEYS[stage].items():
        k = "positive_control" if key == "pc" else key
        envvar = ID_ENV.get(up)
        ids[up] = check_explicit_id(up, env.get(envvar) if envvar else None, sha12, runs)
        if key == "pc":
            paths["positive_control"] = str(Path(ns) / stage_dir_rel("B7", sha12, ids[up]) / HEADLINE / "positive_control.json")
            pc_env = env.get("POSITIVE_CONTROL")
            if pc_env and _abs(ns, pc_env) != _abs(ns, paths["positive_control"]):
                raise ChainPreflightError(f"POSITIVE_CONTROL {pc_env} is not the B7 file {paths['positive_control']}")
        else:
            paths[k] = str(Path(ns) / stage_dir_rel(up, sha12, ids[up]))
    for need_up in ("B1", "B2"):
        ids[need_up] = run_id_for(need_up, sha12, runs)
    if env.get("V3_ROOT"):
        raise ChainPreflightError(f"V3_ROOT is set ({env['V3_ROOT']}): it must stay unset in a v42 run")
    if stage not in ("B3", "B4", "B5", "B6", "B7", "B10", "B11", "B19", "B20"):
        paths.pop("module3_dir")
    if stage not in ("B3", "B4", "B19", "B20"):
        paths.pop("nuclei_dir")
    rep = check_inputs(stage, sha12, ns=ns, protocol=proto, runs=runs, **{k: v for k, v in paths.items() if v})
    out_root = str(Path(ns) / stage_dir_rel(stage, sha12, own))
    ex = {"SHA12": sha12, "RUN_ID": own, "OUT_ROOT": out_root, "V2_NAMESPACE_ROOT": str(ns)}
    for k, v in rep["resolved"].items():
        ex[{"module3_dir": "MODULE3_DIR", "nuclei_dir": "NUCLEI_DIR", "layers_run": "LAYERS_ROOT", "comp_run": "COMP_RUN",
            "commot_root": "COMMOT_ROOT", "v4_run": "V4_ROOT", "nmf_run": "NMF_RUN", "positive_control": "POSITIVE_CONTROL",
            "morans_run": "MORANS_RUN", "hurdle_run": "HURDLE_RUN", "controls_run": "CONTROLS_RUN", "keratin_run": "KERATIN_RUN", "e2_run": "E2_RUN", "coverage_run": "COVERAGE_RUN"}[k]] = v
    ex["COMPARTMENTS_RUN"] = ex.get("COMP_RUN", "")
    if "LAYERS_ROOT" in ex:
        ex["LAYERS_RUN_ID"] = ids["B4"]
    if "NMF_RUN" in ex:
        ex["NMF_RUN_ID"] = ids["B10"]
    ex = {k: v for k, v in ex.items() if v != ""}
    return {"exports": ex, "report": {"stage": stage, "protocol_sha12": sha12, "ids": ids, **rep}}


def write_receipt(report: dict, *, ns, array="", cmd="", code_dir=None, dryrun=False, now=None) -> Path:
    """One JSON per submission under ``<ns>/output/v2_receipts/<sha12>/`` (exclusive create; never overwritten)."""
    ts = (now or _dt.datetime.now(_dt.timezone.utc)).strftime("%Y%m%dT%H%M%S%fZ")
    code_commit = None
    if code_dir:
        r = subprocess.run(["git", "--no-optional-locks", "-C", str(code_dir), "rev-parse", "HEAD"], capture_output=True, text=True)
        code_commit = r.stdout.strip() or None
    rec = {**report, "namespace_root": str(ns), "code_dir": str(code_dir) if code_dir else None, "code_commit": code_commit,
           "array": array, "submit_cmd": cmd, "written_at": ts, "dryrun": bool(dryrun)}
    d = Path(ns) / "output" / "v2_receipts" / report["protocol_sha12"]
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{report['stage']}_{report['ids'][report['stage']]}_{ts}.json"
    with open(f, "x") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)
    return f


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("env", help="check a stage and print export lines for the wrapper")
    r = sub.add_parser("receipt", help="check a stage and write the submission receipt")
    for s in (e, r):
        s.add_argument("--stage", required=True, choices=sorted(REQUIRED))
        s.add_argument("--run-id", required=True)
        s.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    r.add_argument("--array", default="")
    r.add_argument("--cmd", default="")
    r.add_argument("--code-dir", default=None)
    pn = sub.add_parser("pin", help="print the one code commit for (stage, run id): pins HEAD on the first submission, reuses the pin afterwards")
    pn.add_argument("--stage", required=True, choices=sorted(REQUIRED))
    pn.add_argument("--run-id", required=True)
    pn.add_argument("--head", required=True)
    pn.add_argument("--dryrun", action="store_true")
    pn.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    gn = sub.add_parser("generation", help="print the current V2_TS (v42_run_ids.yaml); with --v2-ts, refuse any other (and any superseded) one")
    gn.add_argument("--v2-ts", default=None)
    b3 = sub.add_parser("b3-route", help="check the B3 specimen list against the normal / BELOW_GATE route")
    b3.add_argument("--specimens", required=True)
    b3.add_argument("--below-gate", action="store_true")
    b3.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "b3-route":
            check_b3_route(a.specimens.split(), a.below_gate, PIO.load_protocol(a.protocol))
            return 0
        if a.cmd == "generation":
            gen = L.load_generation()
            if gen["current_v2_ts"] is None:
                raise ChainPreflightError("v42_run_ids.yaml has no current_v2_ts")
            if a.v2_ts:
                check_generation_ts(a.v2_ts, "V2_TS", gen)
            print(gen["current_v2_ts"])
            return 0
        if a.cmd == "pin":
            ns = _ns(None)
            sha, how = pin_code_commit(a.stage, a.run_id, ns, a.head, dryrun=a.dryrun, protocol_path=a.protocol)
            if sha != a.head:
                print(f"chain_preflight: {a.stage} {a.run_id}: using the {how} code commit {sha[:12]}, not HEAD {a.head[:12]}", file=sys.stderr)
            print(sha)
            return 0
        res = resolve_env(a.stage, a.run_id, protocol_path=a.protocol)
        if a.cmd == "env":
            for k, v in sorted(res["exports"].items()):
                print(f"export {k}={shlex.quote(v)}")
            return 0
        f = write_receipt(res["report"], ns=res["exports"]["V2_NAMESPACE_ROOT"], array=a.array, cmd=a.cmd, code_dir=a.code_dir)
        print(f)
        return 0
    except ChainPreflightError as ex:
        print(f"chain_preflight: REFUSED: {ex}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
