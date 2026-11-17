"""A12.9 real-data driver: glue the TDD'd producers into one frozen run.

Runs, at the FROZEN settings (``descriptive.run_h1_specimen``), the hurdle H1 for every
specimen of ``protocol.specimens.hurdle_specimens_all`` (headline P4-S2 gets the 5-run
claim tree; the others are descriptive R00 context), the A01-A07 mechanics sweep for the
headline specimen (``sweep.run_sweep``), H2 (``descriptive.rank_gmr``), H3 (side-by-side
P4-S1 / P4-S2 rows, no contrast), H4 (``descriptive.h4_all``), the negative control
(``controls.negative_control``, P4-S2, R00), and the assembled stamped
``commot_compartment_hurdle.json`` carrying the protocol sha256, every pre-registered
outcome's typed status and the deviations (incl. the 2026-10-15 user override).

The COMMOT-injection positive control runs in ``envs/COMMOT`` (``src/positive_control.py``
via ``slurm/a12_pc_real.sbatch``); this driver is CITEgeist_env and records it as
``pending`` unless ``--positive-control <path>`` points at the written json.

``--finalize`` re-assembles ONLY ``commot_compartment_hurdle.json`` from the already-written
per-specimen h1 / sweep / H2 / H3 / H4 / negative-control files (plus the positive-control
json), so the story can be re-stamped once the COMMOT job finishes, without recomputing.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import a12_9_inputs as A  # noqa: E402
import chain_guard as PF  # noqa: E402
import compartments as C  # noqa: E402
import controls as K  # noqa: E402
import descriptive as D  # noqa: E402
import hurdle as H  # noqa: E402
import protocol_io as PIO  # noqa: E402
import sweep as SW  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/src/a12_9_real.py"
POS_CONTROL_RUN_ID = "pc1"
DEVIATIONS = None  # set from protocol.hurdle.pre_freeze_deviations in main
RNA_COVARIATE = "total"   # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md); must equal descriptive.RNA_COVARIATE
COVERAGE_ARTIFACT = "prod_coverage_total.json"   # the stamped B24 artifact (prod_coverage.py)


def load_coverage(coverage_run) -> dict:
    """The measured synthetic GMR coverage of the stamped B24 run ``coverage_run`` as ``descriptive.coverage_record`` (a stamped INPUT of B8 in total mode).
    Fails closed: the artifact and its stamp must exist, the artifact must be the B24 stage under the same covariate mode, and carry four case values."""
    d = Path(coverage_run)
    f = d / COVERAGE_ARTIFACT
    for p in (f, Path(str(f) + ".stamp.json")):
        if not p.is_file():
            raise PIO.GateFailure(f"B24 coverage input missing: {p} (has B24 run and been stamped?)")
    j = json.loads(f.read_text())
    if j.get("stage") != "B24" or j.get("rna_covariate") != RNA_COVARIATE:
        raise PIO.GateFailure(f"{f}: not a B24 artifact for rna_covariate {RNA_COVARIATE!r} (stage {j.get('stage')!r}, rna_covariate {j.get('rna_covariate')!r})")
    cases = list(j["coverage_by_case"])
    import hashlib

    return D.coverage_record([j["coverage_by_case"][c] for c in cases],
                             {"stage": "B24", "run_id": j["run_id"], "cases": cases, "sha256": hashlib.sha256(f.read_bytes()).hexdigest(), "file": f.name})


# --------------------------------------------------------------------------- one specimen
def _specimen_runs(specimen: str) -> tuple:
    return ("R00", "R01", "R02", "R10", "R11") if specimen == D.HEADLINE_SPECIMEN else ("R00",)


def run_h1(specimen, comp_dir, commot_root, out_base, run_id, protocol_path, coverage=None) -> dict:
    """Compute and write one specimen's H1 (write_h1-gated: frozen settings, verified protocol, empty namespaced target)."""
    obs_all, xy_all = A.load_compartments(comp_dir / f"{specimen}.h5ad")
    el = C.eligible_cancer(obs_all).to_numpy()
    obs, xy = obs_all[el], xy_all[el]
    scores = {r: A.load_score(commot_root / specimen, r, list(obs.index)) for r in _specimen_runs(specimen)}
    res = D.run_h1_specimen(specimen, obs, xy, scores, coverage)
    D.write_h1(out_base / specimen / "h1", res, PRODUCED_BY, protocol_path, run_id)
    return res


# --------------------------------------------------------------------------- H4 N-P MDK values
def eligible_np_mdk(comp_dir: Path, specimen: str) -> np.ndarray:
    """N-P MDK of the specimen's ELIGIBLE cancer compartments (h5ad ``X`` = N_P; ``eligible_cancer`` col recomputed)."""
    import anndata as ad

    a = ad.read_h5ad(str(comp_dir / f"{specimen}.h5ad"), backed="r")
    el = C.eligible_cancer(a.obs).to_numpy()
    midx = int(list(a.var_names).index(C.MDK))
    v = np.asarray(a.X[el][:, midx].todense()).ravel().astype(float)
    a.file.close()
    return v


# --------------------------------------------------------------------------- primary components for the headline control
def primary_components(comp_dir, commot_root, specimen) -> dict:
    """obs (eligible), xy, R00 score, label high, covariate matrix, blocks, pitch — the primary run's objects."""
    obs_all, xy_all = A.load_compartments(comp_dir / f"{specimen}.h5ad")
    el = C.eligible_cancer(obs_all).to_numpy()
    obs, xy = obs_all[el], xy_all[el]
    sc = A.load_score(commot_root / specimen, "R00", list(obs.index))
    if isinstance(sc, dict):
        return {"status": "failed", "reason": f"R00: {sc.get('reason', '')}", "obs": obs, "xy": xy}
    lab = C.select_label(obs["k_share"].to_numpy())
    pitch = H.median_pitch(xy)
    return {"status": "ok", "obs": obs, "xy": xy, "score": np.asarray(sc, float), "label": lab, "pitch": pitch,
            "blocks": H.make_blocks(xy, pitch, H.BLOCK_SIDE_PITCHES, H.BLOCK_MIN_N), "cov": H.covariate_matrix(obs, rna=RNA_COVARIATE)}


# --------------------------------------------------------------------------- negative control (P4-S2, R00)
def run_negative_control(comp_dir, commot_root, specimen) -> dict:
    prim = primary_components(comp_dir, commot_root, specimen)
    if prim["status"] != "ok":
        return D._nonest("commot_run_failed", f"negative control needs an estimable R00 run: {prim['reason']}")
    lab = prim["label"]
    if lab["status"] != "ok":
        return D._nonest(lab["reason_code"], f"negative control needs the primary label: {lab['reason']}")
    return K.negative_control(prim["score"], lab["high"].astype(int), prim["cov"], prim["blocks"], prim["xy"], prim["pitch"])


# --------------------------------------------------------------------------- story assembly
def assemble_story(proto, out_base: Path, run_id: str, pc_path, coverage_source=None) -> dict:
    sha12 = proto["_protocol_sha12"]
    specs = proto["specimens"]["hurdle_specimens_all"]
    headline = D.HEADLINE_SPECIMEN
    h1 = {}
    for s in specs:
        f = out_base / s / "h1" / f"h1_{s}.json"
        h1[s] = json.loads(f.read_text()) if f.is_file() else {"status": "missing", "reason": f"no h1 file at {f}"}
    sweep_f = out_base / headline / "sweep" / f"sweep_{headline}.json"
    sweep = json.loads(sweep_f.read_text()) if sweep_f.is_file() else {"status": "missing", "reason": "no sweep file"}
    misc = {k: (json.loads((out_base / k).read_text()) if (out_base / k).is_file() else {"status": "missing", "reason": f"no {k}"})
            for k in ("h2_rank.json", "h3_side_by_side.json", "h4_all.json", "negative_control.json")}
    pc = (json.loads(Path(pc_path).read_text()) if pc_path and Path(pc_path).is_file()
          else {"status": "pending", "reason": "COMMOT injection runs in envs/COMMOT via slurm/a12_pc_real.sbatch (not a CITEgeist_env output)"})
    head = h1[headline]
    head_mag = (head.get("runs", {}).get("R00", {}).get("magnitude") or {})
    story = {
        "stage": "A12.9",
        "protocol_id": proto["protocol_id"],
        "headline_specimen": headline,
        "claim_tree": head.get("claim_tree"),
        "headline_h1": {"status": head_mag.get("status", head.get("status")), "runs": head.get("runs")},
        "per_specimen": {s: {"status": (h1[s].get("runs", {}).get("R00", {}).get("magnitude", {}) or {}).get("status", h1[s].get("status")),
                            "n_eligible": h1[s].get("n_eligible"), "reason": h1[s].get("reason")} for s in specs},
        "sweep": {"status": sweep.get("status"), "variants": {v: (r.get("magnitude", r) if isinstance(r, dict) else {}) for v, r in sweep.get("variants", {}).items()}},
        "H2": misc["h2_rank.json"], "H3": misc["h3_side_by_side.json"], "H4": misc["h4_all.json"],
        "controls": {"negative_control": misc["negative_control.json"], "positive_control_commot": pc},
        "deviations": DEVIATIONS if DEVIATIONS is not None else proto.get("hurdle", {}).get("pre_freeze_deviations", []),
        "run_id": run_id,
        "rna_covariate": RNA_COVARIATE,
        "protocol_deviation": H.deviation_record(),
        "covariates": list(H.primary_covariates(RNA_COVARIATE)),
        **({"coverage_source": coverage_source} if coverage_source else {}),
        "scope": "headline claim-tree decision + descriptive context; one branch per part for P4-S2; no p-value/q-value/permutation",
    }
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
    import provenance

    stamped = provenance.stamp_payload(PIO.with_protocol_hash(story, proto), PRODUCED_BY)
    target = out_base / "commot_compartment_hurdle.json"
    PIO.atomic_write_json(target, stamped)
    return stamped


# --------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compartments-run", type=Path, required=True, help="A12.7 output dir (holds <SPECIMEN>.h5ad + M.json)")
    ap.add_argument("--commot-root", type=Path, required=True, help="output/commot_runs_module3_7type/<sha12>/<run_id>")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<run_id> (already namespaced)")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--specimen", default=None, help="one specimen (default: all of specimens.hurdle_specimens_all)")
    ap.add_argument("--positive-control", default=None, help="path to the envs/COMMOT positive_control.json (optional)")
    ap.add_argument("--coverage-run", type=Path, default=None, help="the stamped B24 run dir (measured synthetic coverage under rna='total'); REQUIRED for a v42 protocol")
    ap.add_argument("--finalize", action="store_true", help="re-assemble only commot_compartment_hurdle.json from written files")
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    a = ap.parse_args(argv)
    global DEVIATIONS
    proto = PIO.load_protocol(a.protocol)
    # v42: COMMOT must be complete (a lost/staged specimen would otherwise be read as a typed failure and look like a finished hurdle)
    PF.guard_entry("B8", a.protocol, run_id=a.run_id, out_root=a.out_root, comp_run=a.compartments_run, commot_root=a.commot_root,
                   positive_control=a.positive_control, coverage_run=a.coverage_run)
    if a.coverage_run is None:
        raise SystemExit("--coverage-run (the stamped B24 run) is required: the total-covariate hurdle records the B24 coverage, not the per-cell constant")
    DEVIATIONS = proto.get("hurdle", {}).get("pre_freeze_deviations", [])
    out_base = a.out_root
    # v42: job-start stamp gate (before any compute). Only a complete hurdle (all specimens, or --finalize) is a lineage source; a --specimen run is not stamped.
    stamp = None if (a.specimen and not a.finalize) else PF.stage_stamp_begin("B8", a.protocol, out_base / "commot_compartment_hurdle.json", Path(__file__).resolve(), a.run_id)
    coverage = load_coverage(a.coverage_run)         # after the job-start stamp gate (no compute or input read before it), before any specimen is computed
    stamp_dirs = dict(comp_run=a.compartments_run, commot_root=a.commot_root, positive_control=a.positive_control, coverage_run=a.coverage_run)
    if a.finalize:
        assemble_story(proto, out_base, a.run_id, a.positive_control, coverage["source"])
        if stamp:
            stamp.finish(**stamp_dirs)
        print("finalized", out_base / "commot_compartment_hurdle.json", flush=True)
        return 0
    specs = [a.specimen] if a.specimen else proto["specimens"]["hurdle_specimens_all"]
    wt_ids = [s for s in proto["specimens"]["hurdle_specimens_all"] if proto["specimens"]["genotype"][s] != "ESR1_D538G"]
    surgical_wt = list(proto["specimens"]["surgical_wt_specimens"])
    h1 = {}
    for s in specs:
        h1[s] = run_h1(s, a.compartments_run, a.commot_root, out_base, a.run_id, a.protocol, coverage)
        print(f"h1 {s}: {h1[s]['runs']['R00']['magnitude']['status']} n={h1[s]['n_eligible']}", flush=True)

    headline = D.HEADLINE_SPECIMEN
    if headline in specs:
        # A01-A07 sweep (headline only): R00 score for ALL compartments, aligned to all rows
        obs_all, xy_all = A.load_compartments(a.compartments_run / f"{headline}.h5ad")
        score_all = A.load_score(a.commot_root / headline, "R00", list(obs_all.index))
        if isinstance(score_all, dict):
            print(f"sweep: R00 failed ({score_all.get('reason')})", flush=True)
            sweep_res = {"status": "failed", "reason": f"R00: {score_all.get('reason')}"}
        else:
            sweep_res = SW.run_sweep(headline, obs_all, xy_all, score_all, coverage)
            SW.write_sweep(out_base / headline / "sweep", sweep_res, PRODUCED_BY, a.protocol, a.run_id)
        # H2 (GMR rank, absolute-rule specimens), H3 (side-by-side), H4 (N-P MDK ranks)
        entries, stats = {}, {}
        for s, r in h1.items():
            mag = r["runs"]["R00"]["magnitude"]
            entries[s] = {"status": mag["status"], "value": mag.get("value"), "rule": r["label"].get("rule")}
            stats[s] = D.h4_statistics(eligible_np_mdk(a.compartments_run, s))
        H2 = D.rank_gmr(entries, headline)
        H3 = {"specimens": [p for p in ("HCC22-088-P4-S1", headline) if p in h1],
              "rows": {p: {"R00": h1[p]["runs"]["R00"], "R11": h1[p].get("runs", {}).get("R11")} for p in ("HCC22-088-P4-S1", headline) if p in h1},
              "note": "side-by-side descriptive rows only; no contrast, no p-value"}
        H4 = D.h4_all(stats, headline, wt_ids, surgical_wt)
        neg = run_negative_control(a.compartments_run, a.commot_root, headline)
        for name, payload in (("h2_rank.json", H2), ("h3_side_by_side.json", H3), ("h4_all.json", H4), ("negative_control.json", neg)):
            PIO.atomic_write_json(out_base / name, payload)
            print(f"{name} written", flush=True)

    assemble_story(proto, out_base, a.run_id, a.positive_control, coverage["source"])
    if stamp:
        stamp.finish(**stamp_dirs)
    print("done A12.9 out_root=", out_base, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
