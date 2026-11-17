"""Engine v2 breast refit (stage B1): RNA refinement of a FRESH protein start p0, on rebuilt nuclei, for ALL p0 spots.

Subcommands
  refit  One specimen. p0 = the fresh protein-start csv of this run (``v2_inputs.py p0``; stamped), nuclei = the B2 artifact
         (``v2_inputs.py nuclei``; stamped, real file, exactly the p0 barcodes), GEX = the raw filtered h5 of exactly the p0 spots.
         Calls ``run_rna_refinement`` with the frozen settings (multinomial, c=0.1, n_passes=3, anchor_thr=0, n_hvg=2000), the same
         call the validated engine-v2 experiment used, and promotes the result to the module3 layout
         ``<module3_dir>/<specimen>_cell_prop_finetuned_results.csv`` + stamp sidecar (``promote_b1``).
  gate   Aggregates the engine-v2 experiment candidates into breast_gate.json (historical; reads experiment run dirs only).

Nothing here reads ``output/morphology_assignment_v3`` or any v41 proportions: the earlier experiment read the v41 plain
proportions and the v3 stand-in nuclei and was therefore limited to the stand-in spot subset (8 of 12 specimens had fewer spots).

The breast check of the experiment is the v41 one: ``Benchmarking/brightness_estimator/breast_safety.py`` (Spearman of each type's
spot proportion with its EVALUATION-panel RNA score, median over specimens, pooled mean share) and its gate arithmetic
``selection.g3``. Both are imported, not copied.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "Benchmarking" / "brightness_estimator"))
sys.path.insert(0, str(HERE))

import check_v42_lineage as L  # noqa: E402
import nuclei_artifact as NA  # noqa: E402
import v42_stamp as W  # noqa: E402

PRODUCER = "Benchmarking/midkine_remediation/src/rna_refined_props.py"
P0_NAME = "{specimen}_cell_prop_global_results.csv"  # the fresh protein start, written by v2_inputs.py p0 under <p0_dir>/<specimen>/
OUT_NAME = "{specimen}_cell_prop_finetuned_results.csv"  # the module3 layout every downstream stage reads

# candidate id -> run_rna_refinement keywords. a = frozen default; b = V1 (one pass); hvg* = sensitivity of the default.
CANDIDATES = {
    "a_V3_multinomial": {"c": 0.1, "n_passes": 3, "n_hvg": 2000},
    "b_V1": {"c": 0.1, "n_passes": 1, "n_hvg": 2000},
    "a_V3_multinomial_hvg1000": {"c": 0.1, "n_passes": 3, "n_hvg": 1000},
    "a_V3_multinomial_hvg4000": {"c": 0.1, "n_passes": 3, "n_hvg": 4000},
}
FROZEN_CANDIDATE = "a_V3_multinomial"
SENSITIVITY = ("a_V3_multinomial_hvg1000", "a_V3_multinomial_hvg4000")  # reported, never eligible
MAX_WALL_S = 30 * 60
MAX_RSS_GB = 32.0
ROW_SUM_ATOL = 1e-9


def specimens() -> list:
    import run_module3_marker_share as base

    return sorted(base.SPECIMENS)


# ---------------------------------------------------------------------------------------------- validation (pure)
def validate_p0(p0: pd.DataFrame, types: list, specimen: str, expected_n: int | None = None) -> None:
    """The full p0 table: unique barcodes, the protocol's types in order, finite, non-negative, rows sum to 1, expected size."""
    if p0.index.has_duplicates:
        raise ValueError(f"{specimen}: p0 has duplicate barcodes")
    if list(p0.columns) != list(types):
        raise ValueError(f"{specimen}: p0 columns {list(p0.columns)} != protocol types {list(types)}")
    a = p0.to_numpy(dtype=float)
    if a.size == 0 or not np.isfinite(a).all() or (a < 0).any():
        raise ValueError(f"{specimen}: p0 is empty, non-finite or negative")
    if not np.allclose(a.sum(axis=1), 1.0, atol=ROW_SUM_ATOL):
        raise ValueError(f"{specimen}: p0 rows do not sum to 1 (a protein start is a share)")
    if expected_n is not None and len(p0) != expected_n:
        raise ValueError(f"{specimen}: p0 has {len(p0)} spots, protocol expects {expected_n} (specimens.expected_n_spots_qp_file)")


def validate_props_full_p0(props: pd.DataFrame, p0: pd.DataFrame, specimen: str = "?") -> None:
    """The refit covers EVERY p0 spot, in p0's order and columns, with finite shares that sum to 1 (no subsetting)."""
    if not props.index.equals(p0.index):
        gone, extra = p0.index.difference(props.index), props.index.difference(p0.index)
        raise ValueError(f"{specimen}: refit spots != p0 spots ({len(gone)} p0 spots dropped, {len(extra)} spots not in p0)")
    if list(props.columns) != list(p0.columns):
        raise ValueError(f"{specimen}: refit columns {list(props.columns)} != p0 columns {list(p0.columns)}")
    a = props.to_numpy(dtype=float)
    if not np.isfinite(a).all() or (a < -1e-12).any():
        raise ValueError(f"{specimen}: refit has non-finite or negative shares")
    if not np.allclose(a.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError(f"{specimen}: refit rows do not sum to 1")


def refine(gex, p0: pd.DataFrame, nuclei: pd.Series, **override):
    """The frozen v2 call (the same keywords as the validated experiment). Returns (p, info); asserts the settings that ran."""
    from CITEgeist.model.deconvolution.rna_refinement_adata import run_rna_refinement

    kw = {**CANDIDATES[FROZEN_CANDIDATE], **override}
    p, info = run_rna_refinement(gex, p0, nuclei, **kw)
    prm = info["params"]
    if (prm["c"], prm["n_passes"], prm["model"], prm["anchor_thr"]) != (0.1, 3, "multinomial", 0.0) and not override:
        raise RuntimeError(f"refinement ran with {prm}, not the frozen v2 settings")
    return p, info


# ---------------------------------------------------------------------------------------------- B1
NOOP_ATOL = 1e-12


def refuse_noop_refit(refined: pd.DataFrame, p0_csv, specimen: str = "?") -> None:
    """A refit equal to p0 (same spots and shares within 1e-12, at least one spot) did nothing; refuse to promote it."""
    p0 = pd.read_csv(p0_csv, index_col=0)
    if len(refined) > 0 and refined.index.equals(p0.index) and list(refined.columns) == list(p0.columns):
        if np.allclose(refined.to_numpy(dtype=float), p0.to_numpy(dtype=float), rtol=0, atol=NOOP_ATOL):
            raise ValueError(f"{specimen}: the refit equals p0 within {NOOP_ATOL} over {len(refined)} spots (a no-op refinement): not promoting")


def promote_b1(
    refined: pd.DataFrame,
    *,
    specimen: str,
    p0_csv,
    nuclei_csv,
    module3_dir,
    protocol_sha12: str,
    run_id: str,
    started_at: str,
    repo_root=None,
    git_root=None,
    producer: str = PRODUCER,
    external_inputs=(),
    reuse_paths=frozenset(),
    denylist_path=L.DEFAULT_DENYLIST,
) -> dict:
    """Write ``<module3_dir>/<specimen>_cell_prop_finetuned_results.csv`` and its stamp (the module3 layout). Never overwrites."""
    out = Path(module3_dir) / OUT_NAME.format(specimen=specimen)
    if out.exists() or Path(str(out) + ".stamp.json").exists():
        raise FileExistsError(f"{out} exists; a rerun writes a NEW run dir")
    refuse_noop_refit(refined, p0_csv, specimen)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    refined.to_csv(tmp)
    tmp.replace(out)
    try:
        return W.write_stamp(
            out, run_id=run_id, protocol_sha12=protocol_sha12, inputs=[p0_csv, nuclei_csv], producer=producer,
            started_at=started_at, repo_root=repo_root, git_root=git_root, stage_role="b1_rna_compartment_proportions",
            external_inputs=external_inputs, reuse_paths=reuse_paths, denylist_path=denylist_path,
        )
    except BaseException:
        out.unlink()  # an unstamped promoted file must not be left behind for a downstream stage to read
        raise


def refit_specimen(
    sample: str,
    *,
    p0_dir,
    nuclei_dir,
    module3_dir,
    protocol_sha12: str,
    types: list,
    expected_n: int | None = None,
    gex=None,
    data_root=None,
    repo_root=None,
    git_root=None,
    runs: dict | None = None,
    started_at: str | None = None,
    denylist_path=L.DEFAULT_DENYLIST,
) -> pd.DataFrame:
    """One specimen of stage B1. Refuses the v3 stand-in first; validates p0, nuclei coverage and the refit; promotes + stamps."""
    NA.refuse_standin(nuclei_dir)
    # job-start gate: namespace root, clean code tree, producer, output paths -- all before the (long) refinement runs
    ns, git = W.preflight(PRODUCER, [Path(module3_dir) / OUT_NAME.format(specimen=sample)], repo_root=repo_root, git_root=git_root)
    started_at = started_at or W.utc_now()
    runs = runs if runs is not None else L.load_run_ids()
    run_id = L.resolve_run(protocol_sha12, "B1", runs)
    p0_csv = Path(p0_dir) / sample / P0_NAME.format(specimen=sample)
    t0 = time.time()
    p0 = pd.read_csv(p0_csv, index_col=0)
    validate_p0(p0, types, sample, expected_n)
    W.verify_stamp(p0_csv, run_id=run_id, protocol_sha12=protocol_sha12, repo_root=ns, denylist_path=denylist_path)
    nuc = NA.load_nuclei(nuclei_dir, sample, p0.index, protocol_sha12=protocol_sha12, repo_root=ns, runs=runs)
    ext = []
    if gex is None:
        import breast_predict as bp

        data_root = Path(data_root or bp.DATA_ROOT)
        ext = [data_root / sample / "outs" / "filtered_feature_bc_matrix.h5"]
        gex = bp.load_raw_gex(sample, list(p0.index), data_root)
        gex.layers["raw_counts"] = gex.X.copy()
    if list(gex.obs_names) != list(p0.index):
        gex = gex[list(p0.index)].copy()
    p, info = refine(gex, p0, nuc)
    validate_props_full_p0(p, p0, sample)
    stamp = promote_b1(
        p, specimen=sample, p0_csv=p0_csv, nuclei_csv=NA.nuclei_csv(nuclei_dir, sample), module3_dir=module3_dir,
        protocol_sha12=protocol_sha12, run_id=run_id, started_at=started_at, repo_root=ns, git_root=git, external_inputs=ext,
        denylist_path=denylist_path,
    )
    work = Path(module3_dir) / "refit" / sample
    work.mkdir(parents=True)
    from CITEgeist.model.deconvolution.rna_refinement_adata import save_refinement_artifacts

    save_refinement_artifacts({"rna_refinement": info, "cell_prop_base": p0, "cell_prop": p}, work)
    (work / "refit_summary.json").write_text(
        json.dumps(
            {
                "specimen": sample,
                "n_spots": int(len(p)),
                "n_zero_nuclei_spots": int((nuc == 0).sum()),
                "fraction_zero_nuclei_spots": float((nuc == 0).mean()),
                "fallback": info["fallback"],
                "em_converged": [bool(x) for x in np.ravel(info["em_converged"])],
                "n_selected_genes": len(info["selected_genes"]),
                "params": info["params"],
                "code_commit": stamp["code_commit"],
                "seconds": round(time.time() - t0, 1),
            },
            indent=2,
            default=str,
        )
    )
    return p


def parse_time_v(path: Path) -> dict:
    """Wall seconds and peak RSS (kB) from a ``/usr/bin/time -v`` file."""
    txt = Path(path).read_text()
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", txt)
    parts = [float(x) for x in m.group(1).split(":")]
    wall = sum(v * 60**i for i, v in enumerate(reversed(parts)))
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt).group(1))
    return {"wall_s": wall, "max_rss_kb": rss}


def _git(*a) -> str:
    return subprocess.run(["git", "--no-optional-locks", "-C", str(REPO), *a], capture_output=True, text=True).stdout.strip()


def gate(run_dir: Path, stored_rho_json: Path = None) -> dict:
    import breast_safety as bs
    import selection as sel

    run_dir = Path(run_dir)
    specs = specimens()
    stored = json.loads(Path(stored_rho_json or bs.STORED_RHO).read_text())["rho"]["plain"]
    res = {"candidates": {}}
    for cand in CANDIDATES:
        per = {}
        for s in specs:
            f = run_dir / cand / s / "breast_specimen.json"
            per[s] = json.loads(f.read_text()) if f.is_file() else None
        done = [s for s in specs if per[s] is not None]
        entry = {"params": CANDIDATES[cand], "n_specimens_run": len(done), "missing": [s for s in specs if per[s] is None]}
        perf = {}
        for s in done:
            tv = run_dir / cand / s / "time_v.txt"
            if tv.is_file():
                perf[s] = parse_time_v(tv)
        entry["perf"] = {
            "max_wall_s": max((v["wall_s"] for v in perf.values()), default=None),
            "max_rss_gb": (max(v["max_rss_kb"] for v in perf.values()) / 1024 / 1024) if perf else None,
            "n_with_perf": len(perf),
            "per_specimen": perf,
        }
        entry["fallback"] = {s: per[s]["fallback"] for s in done if per[s]["fallback"]}
        entry["n_fallback"] = len(entry["fallback"])
        entry["em_converged"] = {
            "converged": int(sum(sum(per[s]["em_converged"]) for s in done)),
            "total": int(sum(len(per[s]["em_converged"]) for s in done)),
        }
        entry["n_spots_without_nuclei"] = {s: per[s]["n_spots_without_nuclei"] for s in done}
        if done:
            types = list(per[done[0]]["rho_refined"])
            rc = {s: per[s]["rho_refined"] for s in done}
            rp = {s: per[s]["rho_plain_same_spots"] for s in done}
            n_tot = sum(per[s]["n_spots_refined"] for s in done)
            sh_c = {t: sum(per[s]["share_sum_refined"][t] for s in done) / n_tot for t in types}
            sh_p = {t: sum(per[s]["share_sum_plain_same_spots"][t] for s in done) / n_tot for t in types}
            med_c = bs.median_over_specimens(rc)
            med_p = bs.median_over_specimens(rp)
            g3in = {
                "median_spearman": {"candidate": med_c, "plain": med_p},
                "pooled_mean_share": {"candidate": sh_c, "plain": sh_p},
            }
            entry["median_spearman_refined"] = med_c
            entry["median_spearman_v41_plain_same_spots"] = med_p
            entry["median_delta"] = {t: med_c[t] - med_p[t] for t in types}
            entry["pooled_mean_share_refined"] = sh_c
            entry["pooled_mean_share_v41_plain"] = sh_p
            entry["mean_share_shift"] = {t: sh_c[t] - sh_p[t] for t in types}
            entry["n_types_not_worse_by_more_than_0.05"] = int(
                sum(1 for t in types if med_p[t] - med_c[t] <= sel.G3_MAX_WORSE + sel._EPS)
            )
            entry["g3_pass"] = bool(sel.g3(g3in))
            if cand == "a_V3_multinomial":  # the evaluation path must reproduce the stored v41 plain numbers
                pairs = [(per[s]["rho_plain_full"][t], stored[t][s]) for s in done for t in types]
                diffs = [abs(a - b) for a, b in pairs if np.isfinite(a) and np.isfinite(b)]
                nan_mismatch = sum(bool(np.isfinite(a)) != bool(np.isfinite(b)) for a, b in pairs)
                res["stored_v41_plain_reproduction"] = {
                    "max_abs_diff": max(diffs),
                    "nan_mismatches": int(nan_mismatch),
                    "n_compared": len(diffs),
                    "ok": bool(max(diffs) <= bs.RHO_ATOL and nan_mismatch == 0),
                }
        mw, mr = entry["perf"]["max_wall_s"], entry["perf"]["max_rss_gb"]
        entry["all_12_run"] = len(done) == len(specs)
        entry["perf_ok"] = bool(
            entry["perf"]["n_with_perf"] == len(done) and mw is not None and mw <= MAX_WALL_S and mr <= MAX_RSS_GB
        )
        entry["gate_pass"] = bool(entry["all_12_run"] and entry["perf_ok"] and entry.get("g3_pass", False))
        res["candidates"][cand] = entry
    res["sensitivity_candidates"] = list(SENSITIVITY)
    res["formula"] = {
        "check": "Benchmarking/brightness_estimator/breast_safety.py (Spearman vs EVALUATION-panel RNA score, median over 12 specimens, pooled mean share)",
        "gate_arithmetic": "Benchmarking/brightness_estimator/selection.py::g3 (producer of the manuscript's phase1_choices.json)",
        "selection_py_last_commit": _git("log", "-1", "--format=%H", "--", "Benchmarking/brightness_estimator/selection.py"),
        "breast_safety_py_last_commit": _git("log", "-1", "--format=%H", "--", "Benchmarking/brightness_estimator/breast_safety.py"),
        "thresholds": {
            "max_worse": sel.G3_MAX_WORSE,
            "min_types_ok": sel.G3_MIN_TYPES_OK,
            "max_share_shift": sel.G3_MAX_SHARE_SHIFT,
            "min_improvement": sel.G3_MIN_IMPROVEMENT,
        },
        "scoring_note": "refined vs plain compared on the same spots (p0 spots that have nuclei counts); RNA scores computed on the full p0 spot set",
    }
    res["code_git_head"] = _git("rev-parse", "HEAD")
    res["run_dir"] = str(run_dir)
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("refit")
    r.add_argument("--specimen", required=True)
    r.add_argument("--p0-dir", type=Path, required=True)
    r.add_argument("--nuclei-dir", type=Path, required=True, help="the B2 artifact dir (<root>/v2b2); never output/morphology_assignment_v3")
    r.add_argument("--module3-dir", type=Path, required=True, help="NEW dir; receives <specimen>_cell_prop_finetuned_results.csv")
    r.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    r.add_argument("--data-root", type=Path, default=None)
    g = sub.add_parser("gate")
    g.add_argument("--run-dir", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "refit":
        import protocol_io as PIO

        proto = PIO.load_protocol(a.protocol)
        refit_specimen(
            a.specimen, p0_dir=a.p0_dir, nuclei_dir=a.nuclei_dir, module3_dir=a.module3_dir, protocol_sha12=proto["_protocol_sha12"],
            types=proto["inputs"]["types"], expected_n=proto["specimens"]["expected_n_spots_qp_file"][a.specimen], data_root=a.data_root,
        )
    else:
        res = gate(a.run_dir)
        (a.run_dir / "breast_gate.json").write_text(json.dumps(res, indent=2, default=str))
        keys = ("all_12_run", "perf_ok", "g3_pass", "gate_pass")
        print(json.dumps({c: {k: e.get(k) for k in keys} for c, e in res["candidates"].items()}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
