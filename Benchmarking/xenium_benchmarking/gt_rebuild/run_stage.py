"""Per-stage orchestration for the Xenium GT rebuild (called by .commands/gt-rebuild-run.sh and the submit script).

Subcommands (``python -m gt_rebuild.run_stage <cmd> ...``, run from Benchmarking/xenium_benchmarking):
  root                      print the GT root (gt_config.OUT_ROOT; override with env GT_REBUILD_OUT_ROOT)
  check <stage|next> [--resume]   submit-time refuse / resume logic (uses DEPS); prints the stage, exit 1 to refuse
  plan [--resume] [--running a,b]   one 'stage<TAB>deps<TAB>done|running|todo' line per stage, in topological order
  time <stage>              print the Slurm walltime for the stage (STAGE_TIME table)
  variants                  print one 'name resolution seed' line per stability variant (RECIPE)
  name [--clusters CSV --out-dir DIR]   marker-gene cluster naming (marker_names); defaults to the real work/ paths
  tsplit | assemble-labels | gates | stability-variant <name> | stability-report | rnalimits | packet   per-stage work
  rerun-stages <stage> [--also a,b]    stages a `rerun-from <stage>` supersedes (stage + dependents + --also, in order)
  rerun-plan <stage> [--also a,b]      'stage<TAB>relpath' per existing output / done-marker to move (never deletes)
  rerun-manifest <stage> --ts T --reason R [--also a,b]   append superseded_<T> as 'pending' (before any move)
  rerun-complete --ts T     after the moves: write _superseded_<T>/MOVED_OK, flip the section to 'done' (moved list)
  rerun-relocate-sections <ts>   repair: move stage-owned sections (gt_built) into superseded_<ts>['previous_sections']
  rerun-pending             exit 1 (with a recovery hint) while any superseded_* section is still 'pending'
  finish <stage>            post-stage validation / manifest sections, then write work/.done_<stage>.json
Module interfaces elsewhere in gt_rebuild are unchanged; heavy imports are lazy so ``check`` stays instant.
"""

import argparse
import datetime
import json
import os
import subprocess
import sys
from pathlib import Path

from . import gt_config

BASE = gt_config.BASE
REPO = BASE.parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import provenance  # noqa: E402
ROOT = gt_config.OUT_ROOT
WORK = ROOT / "work"
LABELS_DIR = ROOT / "labels"
REPORTS = ROOT / "reports"
COUNTS_CSV = gt_config.OLD_LABELS_DIR / "xenium_cell_gex.csv"
SPIKE_N_CLUSTERS = 15
# stability variants come from RECIPE (single definition): name -> {"resolution", "seed"}
VARIANTS = gt_config.RECIPE["stability_variants"]
# per-stage Slurm walltime (the only place); each stability variant re-clusters, re-names, t-splits and aggregates
DEFAULT_TIME = "06:00:00"
VARIANT_TIME = "08:00:00"
STAGE_TIME = {"stability": "02:00:00", "stability_variant": VARIANT_TIME}

# Ordered stages -> outputs (paths relative to ROOT; a directory means every file under it).
STAGES = {
    "manifest": ["work/manifest_initial.json"],
    "cluster": ["work/level1_clusters.csv", "work/level1_qc.json"],
    "name": ["work/cluster_names/cluster_names.csv", "work/cluster_names/marker_gene_z.csv"],
    "tsplit": [
        "work/tcell_labels.csv",
        "work/tcell_subclusters.csv",
        "work/tcell_gene_contributions.csv",
        "work/tcell_absent_genes.json",
    ],
    "assemble": [
        "labels",
        "ground_truth_6type/ground_truth",
        "ground_truth_7type/ground_truth",
        "ground_truth_6type/excluded_cells_per_spot.csv",
        "ground_truth_7type/excluded_cells_per_spot.csv",
    ],
    "gtgex": ["ground_truth_7type/ground_truth_gex"],
    "gates": ["reports/integrity_report.json", "reports/integrity_info.json"],
    "quality": ["reports/quality"],
    **{
        f"stability_{v}": [
            f"work/stability/{v}/level1_clusters.csv",
            f"work/stability/{v}/cluster_names/cluster_names.csv",
            f"work/stability/{v}/labels",
            f"work/stability/{v}/gt7/ground_truth",
            f"work/stability/{v}/gt7/excluded_cells_per_spot.csv",
        ]
        for v in VARIANTS
    },
    "stability": ["reports/stability.json"],
    "rnalimits": ["reports/rna_limits.json", "reports/tcell_dotplot.png"],
    "packet": ["reports/review_packet"],
}
ORDER = list(STAGES)  # also a topological order of DEPS
# explicit dependency map: a stage may run once every stage in DEPS[stage] has a valid done-marker
DEPS = {
    "manifest": [],
    "cluster": ["manifest"],
    "name": ["cluster"],
    "tsplit": ["name"],
    "assemble": ["tsplit"],
    "gtgex": ["assemble"],
    "gates": ["gtgex"],
    "quality": ["assemble"],
    **{f"stability_{v}": ["manifest"] for v in VARIANTS},
    "stability": ["assemble", *[f"stability_{v}" for v in VARIANTS]],
    "rnalimits": ["assemble"],
    "packet": ["name"],
}
# top-level manifest sections written by a stage's finish(); a rerun that supersedes the stage must relocate them
# (into superseded_<ts>["previous_sections"]) or the fresh run's append_section raises "already exists".
# cluster_cells.py's sections belong to the 'cluster' stage, which rerun-from never supersedes.
STAGE_SECTIONS = {"gtgex": ["gt_built"]}
assert list(DEPS) == ORDER and all(ORDER.index(d) < ORDER.index(s) for s, ds in DEPS.items() for d in ds)


def _sha(path):
    from . import manifest

    return manifest.sha256_file(path)


def _files(rel):
    p = ROOT / rel
    if p.is_dir():
        return sorted(f for f in p.rglob("*") if f.is_file())
    return [p] if p.exists() else []


def marker_path(stage):
    return WORK / f".done_{stage}.json"


def _outputs_exist(stage):
    return any(_files(rel) for rel in STAGES[stage])


def _marker_ok(stage):
    """Return None if the marker is present and every recorded hash still matches, else the reason."""
    m = marker_path(stage)
    try:
        rec = json.loads(m.read_text())
        outs = rec["outputs"]
        if not isinstance(outs, dict) or not outs:
            return f"{stage}: marker has no outputs"
    except Exception as exc:
        return f"{stage}: unreadable marker ({type(exc).__name__})"
    for rel, sha in outs.items():
        f = ROOT / rel
        if not f.is_file():
            return f"{stage}: recorded output missing: {rel}"
        if _sha(f) != sha:
            return f"{stage}: hash changed: {rel}"
    return None


# ----------------------------------------------------------------------------- check
def _orphan(done):
    """Reason string if a done stage has a dependency without a marker, else None."""
    for s in done:
        for d in DEPS[s]:
            if d not in done:
                return f"stage '{d}' has no marker but dependent stage '{s}' does"
    return None


def check(stage_arg, resume):
    """Return (stage, None) to submit, or (None, reason) to refuse."""
    if stage_arg != "next" and stage_arg not in STAGES:
        return None, f"unknown stage '{stage_arg}'; one of {ORDER} or 'next'"
    done = [s for s in ORDER if marker_path(s).exists()]
    if not resume:
        if stage_arg == "next":
            return None, "'next' needs --resume"
        stage = stage_arg
        if stage == "manifest" and (ROOT / "manifest.json").exists():
            return None, f"{ROOT / 'manifest.json'} already exists (use --resume to continue)"
        if marker_path(stage).exists() or _outputs_exist(stage):
            return None, f"output of stage '{stage}' already exists (use --resume only after checking)"
        for d in DEPS[stage]:
            if d not in done:
                return None, f"prerequisite stage '{d}' has no done-marker"
        return stage, None
    # --resume: verify every completed stage, require no gaps, continue with a stage whose deps are all done
    for s in done:
        bad = _marker_ok(s)
        if bad:
            return None, f"resume refused: {bad}"
    orphan = _orphan(done)
    if orphan:
        return None, f"resume refused: {orphan}"
    ready = [s for s in ORDER if s not in done and all(d in done for d in DEPS[s])]
    if not ready:
        return None, "all stages already complete"
    if stage_arg == "next":
        stage = ready[0]
    elif stage_arg in done:
        return None, f"resume refused: stage '{stage_arg}' is already complete"
    elif stage_arg not in ready:
        missing = [d for d in DEPS[stage_arg] if d not in done]
        return None, f"resume refused: stage '{stage_arg}' is not ready; prerequisite(s) {missing} not done"
    else:
        stage = stage_arg
    if stage == "manifest" and (ROOT / "manifest.json").exists():
        return None, "manifest.json exists without a manifest marker; inspect it by hand"
    return stage, None


def plan(resume, running=(), assume_todo=()):
    """Return (rows, None) with (stage, deps, status) in topological order, or (None, reason) to refuse (dag mode).

    ``assume_todo``: stages treated as not done although they have a marker (a `rerun-from` dry run, where nothing was
    moved yet); their markers are not verified.
    """
    bad_names = [r for r in list(running) + list(assume_todo) if r not in STAGES]
    if bad_names:
        return None, f"unknown running stage(s) {bad_names}"
    done = [s for s in ORDER if marker_path(s).exists() and s not in assume_todo]
    for s in done:
        bad = _marker_ok(s)
        if bad:
            return None, f"refused: {bad}"
    orphan = _orphan(done)
    if orphan:
        return None, f"refused: {orphan}"
    if set(running) & set(done):
        return None, f"stage(s) {sorted(set(running) & set(done))} are already done, not running"
    rows = []
    for s in ORDER:
        if s in done:
            status = "done"
        elif s in running:
            status = "running"
        else:
            status = "todo"
            if not resume:
                if s == "manifest" and (ROOT / "manifest.json").exists():
                    return None, f"{ROOT / 'manifest.json'} already exists (use --resume to continue)"
                if _outputs_exist(s):
                    return None, f"output of stage '{s}' already exists (use --resume only after checking)"
        rows.append((s, DEPS[s], status))
    return rows, None


# ----------------------------------------------------------------------------- rerun-from
def rerun_stages(stage, also=()):
    """Stages superseded by ``rerun-from stage``: the stage, ``also``, and every dependent of either, in ORDER.

    ``also`` is explicit on purpose (stability variants and the packet depend on manifest/name, not on assemble, but
    their contents depend on the exclusion rule). ``manifest`` is never re-run. Raises ValueError on an unknown stage.
    """
    bad = [x for x in [stage, *also] if x not in STAGES]
    if bad:
        raise ValueError(f"unknown stage(s) {bad}; one of {ORDER}")
    todo = {stage, *also}
    for s in ORDER:  # ORDER is topological, so one pass reaches every dependent
        if any(d in todo for d in DEPS[s]):
            todo.add(s)
    kept = sorted(todo & {"manifest", "cluster", "name"})
    if kept:
        raise ValueError(
            f"{kept} cannot be re-run by rerun-from (manifest holds the run history; cluster/name are fixed)"
        )
    return [s for s in ORDER if s in todo]


def rerun_paths(stages):
    """(stage, path relative to ROOT) for every existing output and done-marker of ``stages``.

    A stability variant moves its whole ``work/stability/<variant>`` directory (leftovers there would trip a rerun).
    """
    out = []
    for s in stages:
        v = s[len("stability_") :] if s.startswith("stability_") else None
        rels = [f"work/stability/{v}"] if v in VARIANTS else list(STAGES[s])
        for rel in [*rels, f"work/.done_{s}.json"]:
            if (ROOT / rel).exists():
                out.append((s, rel))
    return out


def rerun_pending():
    """Names of manifest sections ``superseded_*`` whose moves were never completed (status 'pending')."""
    from . import manifest

    m = manifest._read(ROOT / "manifest.json")
    return [
        k for k, v in m.items() if k.startswith("superseded_") and isinstance(v, dict) and v.get("status") == "pending"
    ]


def _refuse_if_pending():
    pend = rerun_pending()
    if pend:
        ts = pend[0][len("superseded_") :]
        raise ValueError(
            f"manifest section(s) {pend} are still 'pending': a previous rerun-from stopped part-way. Recovery: look "
            f"in {ROOT / ('_superseded_' + ts)}, finish moving the paths in "
            f"manifest.json[{pend[0]!r}]['moved_planned'] "
            f"there (mv only), then run: python -m gt_rebuild.run_stage rerun-complete --ts {ts}"
        )


def _stage_sections(stages):
    return [k for s in stages for k in STAGE_SECTIONS.get(s, [])]


def rerun_manifest(stage, also, ts, reason):
    """Append manifest section ``superseded_<ts>`` (status 'pending', planned moves) before any move; returns it."""
    from . import manifest

    if not reason.strip():
        raise ValueError("--reason is required")
    _refuse_if_pending()
    stages = rerun_stages(stage, also)
    review = ROOT / manifest.BLINDED_REVIEW
    planned = [rel for _, rel in rerun_paths(stages)]
    if not planned:
        raise ValueError("nothing to move: no output or done-marker of the superseded stages exists")
    payload = {
        "status": "pending",
        "at": _now(),
        "reason": reason,
        "from_stage": stage,
        "also": list(also),
        "stages": stages,
        "superseded_dir": f"_superseded_{ts}",
        "moved_planned": planned,
        "moved": [],
        "kept_stages": [s for s in ORDER if s not in stages],
        "blinded_review_superseded": {
            "path": manifest.BLINDED_REVIEW,
            "sha256": manifest.sha256_file(review) if review.is_file() else None,
            "note": "reviewed the pre-exclusion packet; the manifest section 'blinded_review' is kept as history. "
            "The review must be redone on the new reports/review_packet; this file is NOT moved by rerun-from.",
        },
    }
    mpath = ROOT / "manifest.json"
    data = manifest._read(mpath)
    if f"superseded_{ts}" in data:
        raise KeyError(f"manifest section already exists: superseded_{ts}")
    # relocate (never delete) stage-owned sections into the record, in the same write that creates it
    payload["previous_sections"] = {k: data.pop(k) for k in _stage_sections(stages) if k in data}
    data[f"superseded_{ts}"] = payload
    manifest._write(mpath, data)
    return payload


def rerun_relocate_sections(ts):
    """Repair a root whose superseded_<ts> section predates section relocation: move the stage-owned top-level
    sections of its superseded stages into ``previous_sections``. Idempotent (already relocated -> no-op).
    Raises ValueError if the section is missing, nothing is to relocate and nothing was relocated before, or a
    section exists both at top level and in ``previous_sections`` (a fresh one; refuse rather than clobber)."""
    from . import manifest

    mpath = ROOT / "manifest.json"
    data = manifest._read(mpath)
    key = f"superseded_{ts}"
    sec = data.get(key)
    if not isinstance(sec, dict):
        raise ValueError(f"no manifest section {key}")
    if not isinstance(sec.get("stages"), list) or not isinstance(sec.get("at"), str):
        raise ValueError(f"manifest section {key} is malformed (needs 'stages' list and 'at' string)")
    prev = sec.get("previous_sections", {})
    names = _stage_sections(sec["stages"])
    clash = [k for k in names if k in data and k in prev]
    if clash:
        raise ValueError(f"{clash} exist both at top level and in {key}['previous_sections']; not overwriting")
    todo = [k for k in names if k in data]
    if not todo:
        if prev:
            print(f"already relocated: {sorted(prev)}")
            return []
        raise ValueError(f"nothing to relocate: none of {names} is a top-level manifest section")
    for k in todo:  # only superseded (older) sections may be buried; a fresh one from the rerun must stay
        stamp = data[k].get(f"{k}_at") if isinstance(data[k], dict) else None
        if not isinstance(stamp, str) or stamp >= sec["at"]:
            raise ValueError(
                f"top-level section {k!r} (timestamp {stamp!r}) is not older than {key}['at'] ({sec['at']}); "
                "it looks like the fresh rerun's section, refusing to relocate it"
            )
    n = 0
    while True:
        suffix = _now().replace(":", "").replace("-", "") + (f"-{n}" if n else "")
        backup = mpath.with_name("manifest.json.bak-" + suffix)
        if not backup.exists():
            break
        n += 1
    backup.write_text(mpath.read_text())
    for k in todo:
        prev[k] = data.pop(k)
    data[key] = {**sec, "previous_sections": prev}
    manifest._write(mpath, data)
    print(f"relocated {todo} into {key}['previous_sections'] (backup: {backup.name})")
    return todo


def rerun_complete(ts):
    """After the moves: verify every planned path is at the destination and gone from the source, write MOVED_OK and
    flip the section to status 'done' with the actually-moved list. Raises ValueError (section stays 'pending')."""
    from . import manifest

    key = f"superseded_{ts}"
    sec = manifest._read(ROOT / "manifest.json").get(key)
    if not isinstance(sec, dict) or sec.get("status") != "pending":
        raise ValueError(f"no pending manifest section {key}")
    dest = ROOT / sec["superseded_dir"]
    moved = [r for r in sec["moved_planned"] if (dest / r).exists() and not (ROOT / r).exists()]
    left = [r for r in sec["moved_planned"] if r not in moved]
    if left:
        raise ValueError(f"{len(left)} planned path(s) not moved yet, e.g. {left[:3]}; section {key} stays 'pending'")
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "MOVED_OK").write_text(json.dumps({"at": _now(), "moved": moved}, indent=2))
    manifest.replace_section(
        ROOT / "manifest.json", key, {**sec, "status": "done", "moved": moved, "completed_at": _now()}
    )


# ----------------------------------------------------------------------------- helpers
def _json_default(o):
    return o.item() if hasattr(o, "item") else str(o)


def _dump(obj, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default))


def _resume():
    """True only on a --resume submission (submit script exports GT_RESUME=1)."""
    return os.environ.get("GT_RESUME") == "1"


def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def _read_level1(path=None):
    import pandas as pd

    df = pd.read_csv(path or WORK / "level1_clusters.csv", dtype=str, index_col="cell_id")
    return df["cluster"]


def _read_names(path=None):
    from .cluster_names import read_cluster_names

    return read_cluster_names(path or WORK / "cluster_names" / "cluster_names.csv")


def _load_adata():
    from .cluster_cells import load_xenium_rna

    xdir = Path(gt_config.XENIUM_DIR)
    return load_xenium_rna(xdir / "cell_feature_matrix.h5", xdir / "cells.parquet")


def _log_adata(adata):
    import scanpy as sc

    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=gt_config.RECIPE["normalize_total_target_sum"])
    sc.pp.log1p(a)
    return a


# ----------------------------------------------------------------------------- stages
def stage_manifest():
    from . import manifest

    if (ROOT / "manifest.json").exists():
        raise FileExistsError("manifest.json already exists")
    counts = manifest.assert_reference_is_may(gt_config.REF_LABELS_6)
    xdir = Path(gt_config.XENIUM_DIR)
    inputs = {
        "ref_labels_6": gt_config.REF_LABELS_6,
        "reference_h5ad": gt_config.REF_LABELS_6.parents[1] / "GSE156632_combined.h5ad",
        "xenium_cell_feature_matrix_h5": xdir / "cell_feature_matrix.h5",
        "mapping_csv": gt_config.MAPPING_CSV,
    }
    for f in sorted(Path(gt_config.REF_MTX_DIR).iterdir()):
        if f.is_file():
            inputs[f"ref_mtx/{f.name}"] = f
    payload = {
        "created_at": _now(),
        "git_head": provenance.code_commit(REPO),
        "recipe": gt_config.RECIPE,
        "seeds": {"seed": gt_config.RECIPE["seed"], "seed_stability": gt_config.RECIPE["seed_stability"]},
        "versions": manifest.collect_versions(),
        "reference_label_counts": counts,
        "input_paths": {k: str(v) for k, v in inputs.items()},
        "inputs_sha256": {k: manifest.sha256_file(v) for k, v in inputs.items()},
    }
    manifest.write_manifest(ROOT / "manifest.json", payload)
    # manifest.json grows (later stages append sections), so the marker hashes this frozen copy instead
    (WORK / "manifest_initial.json").write_text((ROOT / "manifest.json").read_text())


def stage_name(clusters_csv=None, out_dir=None):
    """Name clusters from the fixed RNA marker rules (marker_names); the variant stages pass their own paths."""
    import pandas as pd

    from .marker_names import score_clusters, write_outputs

    adata, _ = _load_adata()
    a = _log_adata(adata)
    X = a.X.toarray() if hasattr(a.X, "toarray") else a.X
    logx = pd.DataFrame(X, index=a.obs_names, columns=a.var_names)
    level1 = _read_level1(clusters_csv)
    res = score_clusters(logx.loc[level1.index], level1, gt_config.MARKER_SETS, gt_config.MARKER_AMBIGUOUS_MARGIN)
    write_outputs(res, out_dir or WORK / "cluster_names")
    print(res.to_string(index=False))


def stage_tsplit():
    from .tcell_split import run_split_full

    adata, _ = _load_adata()
    labels, table, contrib, absent = run_split_full(adata, _read_level1(), _read_names(), gt_config.RECIPE["seed"])
    labels.to_csv(WORK / "tcell_labels.csv", index=False)
    table.to_csv(WORK / "tcell_subclusters.csv", index=False)
    contrib.to_csv(WORK / "tcell_gene_contributions.csv", index=False)
    _dump({"absent_genes": absent}, WORK / "tcell_absent_genes.json")


def _assign(level1, names, tsplit_labels, out_dir, spot_dirs=()):
    """Label the kept cells (excluded clusters' cells are dropped); write the exclusion files (assemble_labels)."""
    import pandas as pd

    from .assemble_labels import assign_cell_labels, write_exclusion_files, write_label_files

    df = assign_cell_labels(level1, names, tsplit_labels)
    write_label_files(df, out_dir)
    df.to_csv(Path(out_dir) / "cell_labels_full.csv", index=False)
    mapping = pd.read_csv(gt_config.MAPPING_CSV, index_col=0)
    summary = write_exclusion_files(level1, names, df, mapping, out_dir, spot_dirs)
    print(
        f"excluded clusters {summary['clusters']}: {summary['n_excluded']} of {summary['n_total']} cells "
        f"({summary['excluded_fraction']:.4f}); spots with every cell excluded: {summary['n_spots_all_excluded']}"
    )
    return df


def stage_assemble_labels():
    import pandas as pd

    _assign(
        _read_level1(),
        _read_names(),
        pd.read_csv(WORK / "tcell_labels.csv", dtype={"cell_id": str}),
        LABELS_DIR,
        [ROOT / "ground_truth_6type", ROOT / "ground_truth_7type"],
    )


def stage_gates():
    from . import integrity_gates as ig

    names = _read_names()
    qc = json.loads((WORK / "level1_qc.json").read_text())
    excluded_csv = LABELS_DIR / "excluded_cells.csv"
    results = ig.run_gates(
        ROOT, LABELS_DIR / "cell_labels_7type.csv", gt_config.MAPPING_CSV, names, qc["dropped_ids"], excluded_csv
    )
    ig.write_report(results, REPORTS / "integrity_report.json")
    _dump(ig.report_only(names, ROOT, qc, excluded_csv), REPORTS / "integrity_info.json")
    bad = [g.name for g in results if not g.ok]
    print(f"gates: all_ok={not bad} failed={bad}")  # a failed gate is reported, never a stage failure


def stage_stability_variant(name):
    """Per-variant work after cluster + name (done in the shell): t-split, assign labels, aggregate to spots."""
    from .tcell_split import run_split_full

    seed = VARIANTS[name]["seed"]
    vdir = WORK / "stability" / name
    adata, _ = _load_adata()
    level1 = _read_level1(vdir / "level1_clusters.csv")
    names = _read_names(vdir / "cluster_names" / "cluster_names.csv")
    tl, _, _, _ = run_split_full(adata, level1, names, seed)
    _assign(level1, names, tl, vdir / "labels", [vdir / "gt7"])
    agg = BASE / "ground_truth_singler" / "src" / "aggregate_to_spots.py"
    subprocess.run(
        [sys.executable, str(agg), "--labels-csv", str(vdir / "labels" / "cell_labels_7type.csv"),
         "--mapping-csv", str(gt_config.MAPPING_CSV), "--output-dir", str(vdir / "gt7"), "--n-types", "7",
         "--prefix", "Xenium"],
        check=True,
    )  # fmt: skip


def cell_overlap(base_ids, variant_ids) -> dict:
    """Cells kept by both runs: counts, and ``overlap_fraction`` = both / (kept by either run)."""
    b, v = set(base_ids), set(variant_ids)
    both, either = len(b & v), len(b | v)
    return {
        "n_compared": both,
        "n_only_base": len(b) - both,
        "n_only_variant": len(v) - both,
        "overlap_fraction": both / either if either else 0.0,
    }


def stage_stability_report():
    import pandas as pd

    from .stability import spot_prop_change, stability_summary

    # variant-count check: every RECIPE variant must have produced its outputs (variant stages are separate jobs)
    present = [v for v in VARIANTS if (WORK / "stability" / v / "labels" / "cell_labels_full.csv").is_file()]
    if len(VARIANTS) < 1 or len(present) != len(VARIANTS):
        raise RuntimeError(f"stability: found {len(present)} of {len(VARIANTS)} variants: {present}")
    base7 = pd.read_csv(LABELS_DIR / "cell_labels_7type.csv", dtype={"cell_id": str}).set_index("cell_id").singler_label
    base6 = pd.read_csv(LABELS_DIR / "cell_labels_6type.csv", dtype={"cell_id": str}).set_index("cell_id").singler_label
    var7, var6, change, n_clusters, overlap = {}, {}, {}, {}, {}
    for name in VARIANTS:
        vdir = WORK / "stability" / name
        df = pd.read_csv(vdir / "labels" / "cell_labels_full.csv", dtype={"cell_id": str})
        # each run drops its own excluded clusters' cells: agreement is over the cells kept by both runs
        common = base7.index.intersection(df["cell_id"])
        overlap[name] = cell_overlap(base7.index, df["cell_id"])
        var7[name] = df.set_index("cell_id").label7.loc[common]
        var6[name] = df.set_index("cell_id").label6.loc[common]
        n_clusters[name] = int(_read_level1(vdir / "level1_clusters.csv").nunique())
        change[name] = {}
        for r in range(5):
            a = pd.read_csv(ROOT / "ground_truth_7type" / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_col=0)
            b = pd.read_csv(vdir / "gt7" / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_col=0)
            change[name][str(r)] = spot_prop_change(a, b).to_dict("index")
    _dump(
        {
            "base_n_clusters": int(_read_level1().nunique()),
            "variant_n_clusters": n_clusters,
            "cell_overlap": overlap,
            "label7_agreement": {n: stability_summary({n: v}, base7.loc[v.index])[n] for n, v in var7.items()},
            "label6_agreement": {n: stability_summary({n: v}, base6.loc[v.index])[n] for n, v in var6.items()},
            "spot_prop_change_7type": change,
        },
        REPORTS / "stability.json",
    )


def stage_rnalimits():
    import pandas as pd

    from .rna_limits import dotplot, rna_limits

    counts = pd.read_csv(COUNTS_CSV, index_col=0)
    labels7 = (
        pd.read_csv(LABELS_DIR / "cell_labels_7type.csv", dtype={"cell_id": str}).set_index("cell_id").singler_label
    )
    table = pd.read_csv(WORK / "tcell_subclusters.csv", dtype={"subcluster": str})
    out = rna_limits(counts, labels7, _read_level1(), _read_names(), table)
    _dump(out, REPORTS / "rna_limits.json")
    tl = pd.read_csv(WORK / "tcell_labels.csv", dtype=str).set_index("cell_id").subcluster
    adata, _ = _load_adata()
    dotplot(_log_adata(adata), tl, REPORTS / "tcell_dotplot.png")


def stage_packet():
    from .review_packet import build_packet

    adata, _ = _load_adata()
    import pandas as pd

    z = pd.read_csv(WORK / "cluster_names" / "marker_gene_z.csv", index_col="cluster", dtype={"cluster": str})
    build_packet(_log_adata(adata), _read_level1(), _read_names(), REPORTS / "review_packet", z)


def finish(stage):
    """Post-stage validation and manifest sections, then the done-marker (hashes of every output file)."""
    from . import manifest

    extra = {}
    if stage == "cluster":
        n = int(_read_level1().nunique())
        extra = {"n_clusters": n, "spike_n_clusters": SPIKE_N_CLUSTERS, "equals_spike": n == SPIKE_N_CLUSTERS}
        print(f"cluster: n clusters = {n} (spike had {SPIKE_N_CLUSTERS}; report only)")
    elif stage == "name":
        names = _read_names()
        print(f"name: read_cluster_names OK, {len(names)} clusters")
    elif stage == "assemble":
        for n in (6, 7):
            got = sorted((ROOT / f"ground_truth_{n}type" / "ground_truth").glob("Xenium_region_*_prop.csv"))
            if len(got) != 5:
                raise RuntimeError(f"{n}-type: expected 5 prop CSVs, found {len(got)}")
        print("assemble: 10 prop CSVs present")
    elif stage == "gtgex":
        manifest.put_section(ROOT / "manifest.json", "gt_built", {"gt_built_at": _now()}, _resume())
    outs = {}
    for rel in STAGES[stage]:
        files = _files(rel)
        if not files:
            raise RuntimeError(f"stage '{stage}': expected output missing: {rel}")
        for f in files:
            if f.stat().st_size == 0:
                raise RuntimeError(f"stage '{stage}': empty output: {f}")
            outs[str(f.relative_to(ROOT))] = _sha(f)
    _dump({"stage": stage, "completed": _now(), "outputs": outs, **extra}, marker_path(stage))
    print(f"{stage}: marker written ({len(outs)} files hashed)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd")
    ap.add_argument("arg", nargs="?")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--clusters", help="name: level1_clusters.csv (default: work/level1_clusters.csv)")
    ap.add_argument("--out-dir", help="name: output directory (default: work/cluster_names)")
    ap.add_argument("--running", default="", help="plan: comma-separated stages already in flight")
    ap.add_argument("--assume-todo", default="", help="plan: comma-separated stages to treat as not done (dry run)")
    ap.add_argument("--also", default="", help="rerun-*: extra comma-separated stages to supersede")
    ap.add_argument("--ts", help="rerun-manifest: UTC timestamp of the _superseded_<ts> directory")
    ap.add_argument("--reason", default="", help="rerun-manifest: why the stages are superseded")
    a = ap.parse_args(argv)
    simple = {
        "manifest": stage_manifest,
        "tsplit": stage_tsplit,
        "assemble-labels": stage_assemble_labels,
        "gates": stage_gates,
        "stability-report": stage_stability_report,
        "rnalimits": stage_rnalimits,
        "packet": stage_packet,
    }
    if a.cmd == "root":
        print(ROOT)
    elif a.cmd == "check":
        stage, why = check(a.arg or "", a.resume)
        if why:
            print(f"REFUSED: {why}", file=sys.stderr)
            sys.exit(1)
        print(stage)
    elif a.cmd == "plan":
        rows, why = plan(a.resume, [r for r in a.running.split(",") if r], [r for r in a.assume_todo.split(",") if r])
        if why:
            print(f"REFUSED: {why}", file=sys.stderr)
            sys.exit(1)
        for stage, deps, status in rows:
            print(stage, ",".join(deps) or "-", status, sep="\t")
    elif a.cmd in ("rerun-stages", "rerun-plan", "rerun-manifest", "rerun-complete", "rerun-pending",
                   "rerun-relocate-sections"):
        also = [x for x in a.also.split(",") if x]
        try:
            if a.cmd == "rerun-pending":
                _refuse_if_pending()
                return
            if a.cmd == "rerun-relocate-sections":
                rerun_relocate_sections(a.arg or a.ts or "")
                return
            if a.cmd == "rerun-complete":
                if not a.ts:
                    raise ValueError("--ts is required")
                rerun_complete(a.ts)
                return
            stages = rerun_stages(a.arg or "", also)
            if a.cmd == "rerun-stages":
                print(",".join(stages))
            elif a.cmd == "rerun-plan":
                for s_, rel in rerun_paths(stages):
                    print(s_, rel, sep="\t")
            else:
                if not a.ts:
                    raise ValueError("--ts is required")
                rerun_manifest(a.arg, also, a.ts, a.reason)
        except (ValueError, FileNotFoundError, KeyError) as exc:
            print(f"REFUSED: {exc}", file=sys.stderr)
            sys.exit(1)
    elif a.cmd == "time":
        print(STAGE_TIME.get(a.arg, VARIANT_TIME if str(a.arg).startswith("stability_") else DEFAULT_TIME))
    elif a.cmd == "name":
        WORK.mkdir(parents=True, exist_ok=True)
        stage_name(a.clusters, a.out_dir)
    elif a.cmd == "stability-variant":
        if a.arg not in VARIANTS:
            ap.error(f"unknown stability variant {a.arg}")
        WORK.mkdir(parents=True, exist_ok=True)
        stage_stability_variant(a.arg)
    elif a.cmd == "variants":  # one "name resolution seed" line per stability variant (read by the shell)
        for n, v in VARIANTS.items():
            print(n, v["resolution"], v["seed"])
    elif a.cmd == "finish":
        finish(a.arg)
    elif a.cmd in simple:
        WORK.mkdir(parents=True, exist_ok=True)
        simple[a.cmd]()
    else:
        ap.error(f"unknown command {a.cmd}")


if __name__ == "__main__":
    main()
