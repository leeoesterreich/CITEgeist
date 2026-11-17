"""B21: keratin presence, MDK and the estrogen programme (check 7), and ESR1 regulon activity vs MDK (check 8), on eligible cancer compartments. Descriptive only.

Promoted from the exploratory checks 7 and 8 (``exploratory_v42/check7_keratin_presence_e2.py``, ``check8_er_regulon.py``, git history). The analysis code lives in
``keratin_e2.py`` (check 7) and ``er_regulon.py`` (check 8), unchanged; this file is the stage entry that runs check 7 and then check 8 in ONE job and writes ONE
artifact, ``e2_covariation_results.json``: per check the summary fields of the exploratory ``check{7,8}_summary.json`` (gates, cohort lines, per-specimen records), the
per-specimen tables of the exploratory CSVs as records, the check 8 regulon table and meta, the sha256 of the saved regulon files, and a ``config`` block.
The same tables are also written beside the artifact as csv (their sha256 is in the artifact).

Gates (fail closed; ``gate_failure.json``, no artifact, no stamp, exit 2)
  * check 7, headline specimen P4-S2: the keratin-present rule reproduces the B8 A07 group sizes (n_high, n_low) and active-sender counts.
    Reference: ``<hurdle-run>/<headline>/sweep/sweep_<headline>.json`` ``variants.A07`` (the B8 stage's own sweep file; pinned in the B21 stamp).
  * check 8 vs check 7: per specimen the eligible count and the MDK>0 count equal the values check 7 just computed in this run (the exploratory gate read the
    check 7 summary file; same tolerance, 0).
Seeds ``SeedSequence([20261107, crc32(label)])`` per test label, B = 9999 surrogates, 200 matched controls, k = 6: unchanged. No p-value is called a finding.

Inputs come from the CLI (set by ``.commands/submit-v2-e2.sh`` through ``chain_preflight.py``): the B5 run (``--compartments-run``), the B6 COMMOT run
(``--commot-root``), the B8 run (``--hurdle-run``) and ``--out-root``. Code dependency (not an input): ``_e2_scoring.py`` (this folder)
(MSigDB Li 2023 E2 gene sets, hard-coded in that module) and the committed ``regulon_inputs/`` tables.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import chain_guard as PF  # noqa: E402
import er_regulon as C8  # noqa: E402
import hurdle as H  # noqa: E402
import keratin_e2 as C7  # noqa: E402
import keratin_secretory_association as K  # noqa: E402
import promoted_common as CM  # noqa: E402
import protocol_io as PIO  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/e2_covariation.py"
STAGE = "B21"
ARTIFACT = "e2_covariation_results.json"
HEADLINE = C7.HEADLINE
READING = {"check7": "Descriptive only. Per-specimen numbers are not corrected for multiple looks except the named BH families. Nothing here is a claim about ESR1-D538G or about cause. "
                     "A null control percentile near 50 means the association is no stronger than depth-matched random genes; a percentile above 95 is reported as such and is not a finding.",
           "check8": "Descriptive only. Per-specimen numbers are not corrected for multiple looks except the named BH families. Nothing here is a claim about ESR1-D538G or about cause. "
                     "A control percentile near 50 means the association is no stronger than depth-matched random genes; above 95 is reported as such and is not a finding. "
                     "If validity checks (a) or (b) fail, the activity measure is reported as not readable in this data, and the main-question numbers are still reported but flagged."}
ESR1_SIGNATURE = "not used: media-3.xlsx absent from main tree and no column identified as the ESR1 signature (see rule file)"
RULE_FILES = {"check7": "Benchmarking/midkine_remediation/exploratory_v42/check7_DECISION_RULE.md", "check8": "Benchmarking/midkine_remediation/exploratory_v42/check8_DECISION_RULE.md"}


def config() -> dict:
    """Every seed, draw count, threshold, label and reading rule of the stage (written into the artifact)."""
    return {"seed": C7.SEED, "seed_rule": "SeedSequence([seed, crc32(label)]) per test label; control draws: control_pairs.make_draws with seed string '<specimen>|check7' / '<specimen>|check8'",
            "rna_covariate": C8.RNA_COVARIATE, "covariates": list(H.primary_covariates(C8.RNA_COVARIATE)), "protocol_deviation": H.deviation_record(),
            "n_surrogates": C7.B, "n_controls": C7.N_CTRL, "k": K.K, "min_eligible": C7.MIN_EL, "min_group": C7.MIN_GROUP, "min_mdk_nonzero": C7.MIN_NZ,
            "phases": list(C7.PHASES), "keratin_presence_threshold_main": C7.THR_MAIN, "keratin_presence_threshold_sens": C7.THR_SENS,
            "check8": {"min_targets": C8.MIN_TARGETS, "min_b_group": C8.MIN_B_GROUP, "b_min_prop": C8.B_MIN_PROP, "b_min_counts": C8.B_MIN_COUNTS, "tf": C8.TF,
                       "variants": {v: {"network_file": f, "excluded_from_features": sorted(e)} for v, (f, e) in C8.VARIANTS.items()}},
            "reading_rule": READING, "esr1_signature": ESR1_SIGNATURE, "rule_files": RULE_FILES, "labels": {"verdict_labels": "none; values only", "a07_gate_names": ["A07 n_high", "A07 n_low", "A07 n_active_high", "A07 n_active_low"]},
            "descriptive_only": True}


def sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_csv(path: Path, rows: list, first: tuple) -> None:
    names = sorted({k for r in rows for k in r}, key=lambda k: (k not in first, k))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=names)
        w.writeheader()
        w.writerows(rows)


def specimens_of(protocol: dict) -> list:
    return sorted(protocol["specimens"]["ordered"], key=lambda s: (s != HEADLINE, s))


def compute(comp_run, commot_root, hurdle_run, protocol: dict, specimens=None, B: int = C7.B, log=print) -> dict:
    """Run check 7 then check 8 on every specimen. Returns the artifact body (no I/O besides reading the inputs); raises ``CM.GateFailure`` on a gate failure."""
    specs = list(specimens) if specimens else specimens_of(protocol)
    genotype, pairs = protocol["specimens"]["genotype"], protocol["specimens"]["pairs_biopsy_surgery"]
    commot_root = Path(commot_root)
    a07 = C7.a07_reference(hurdle_run)
    per7, gates7 = {}, []
    for s in specs:
        t = time.time()
        per7[s] = C7.analyze(s, commot_root, B=B, comp_run=comp_run)
        log(f"check7 {s}: n={per7[s]['n_eligible']} B={per7[s]['partB'].get('status')} C={per7[s]['partC'].get('status')} {time.time()-t:.0f}s")
        if s == HEADLINE:
            gates7 = C7.gate_a07(per7[s], a07)
            log("GATE PASSED (P4-S2 reproduces A07 group sizes and active counts)")
    coh7 = C7.cohort(per7)
    rows7 = C7.flat_rows(per7, genotype)  # after cohort(): q values are attached
    nets = {v: pd.read_csv(C8.REG_DIR / f, sep="\t") for v, (f, _) in C8.VARIANTS.items()}
    per8 = {}
    for s in specs:
        t = time.time()
        per8[s] = C8.analyze8(s, nets, commot_root, B=B, comp_run=comp_run)
        v0 = per8[s].get("variants", {}).get("collectri", {})
        log(f"check8 {s}: n={per8[s]['n_eligible']} status={per8[s]['status']} collectri={v0.get('status')} measured={v0.get('regulon', {}).get('n_measured_on_panel')} {time.time()-t:.0f}s")
    gates8 = C8.gate_vs_check7(per8, per7)
    log(f"GATE PASSED ({len(gates8)} comparisons with check 7)")
    coh8 = C8.cohort8(per8, pairs)
    panel_reg = {}
    for s, r in per8.items():
        for v, blk in r.get("variants", {}).items():
            panel_reg.setdefault(v, blk.pop("_panel_regulon"))
    for r in per8.values():
        for blk in r.get("variants", {}).values():
            blk.pop("_panel_regulon", None)
    rows8 = C8.flat_rows8(per8, genotype)
    reg_rows = [{"variant": v, **row} for v, rows in panel_reg.items() for row in rows]
    return {"check7": {"gates": gates7, "a07_as_stored": a07, "cohort": coh7, "per_specimen": per7, "per_specimen_table": rows7},
            "check8": {"gates": gates8, "cohort": coh8, "per_specimen": per8, "per_specimen_table": rows8, "esr1_regulon_used": reg_rows}}


def run(comp_run, commot_root, hurdle_run, out_root, protocol_path, run_id: str = "", specimens=None, B: int = C7.B, meta_extra=None) -> int:
    """Compute, write ``e2_covariation_results.json`` + csv tables + ``_COMPLETE.json``. 0 ok, 2 gate failure (``gate_failure.json``, nothing else)."""
    protocol = PIO.load_protocol(protocol_path)
    out = Path(out_root)
    out.mkdir(parents=True, exist_ok=True)
    code = HERE.parents[1]
    git = lambda *a: subprocess.run(["git", "-C", str(code), *a], capture_output=True, text=True).stdout.strip()  # noqa: E731
    reg_files = {p.name: sha256(p) for p in sorted(C8.REG_DIR.iterdir()) if p.is_file()}
    meta = {"stage": STAGE, "run_id": run_id, "produced_by": PRODUCED_BY, "descriptive_only": True, "config": config(),
            "compartments_run": str(comp_run), "commot_run": str(commot_root), "hurdle_run": str(hurdle_run),
            "regulon_meta": json.loads((C8.REG_DIR / "meta.json").read_text()), "regulon_inputs_sha256": reg_files,
            "blas_threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, "git_commit": git("rev-parse", "HEAD"), **(meta_extra or {})}
    try:
        body = compute(comp_run, commot_root, hurdle_run, protocol, specimens, B, log=lambda m: print(m, flush=True))
    except CM.GateFailure as e:
        PIO.atomic_write_json(out / "gate_failure.json", {**meta, "gate_failure": str(e)})
        print(str(e), file=sys.stderr, flush=True)
        return 2
    write_csv(out / "check7_per_specimen.csv", body["check7"]["per_specimen_table"], ("specimen", "genotype", "n_eligible", "part", "test", "form"))
    write_csv(out / "check8_per_specimen.csv", body["check8"]["per_specimen_table"], ("specimen", "genotype", "n_eligible", "variant", "test", "form"))
    pd.DataFrame(body["check8"]["esr1_regulon_used"]).to_csv(out / "esr1_regulon_used.csv", index=False)
    meta["sidecar_files_sha256"] = {n: sha256(out / n) for n in ("check7_per_specimen.csv", "check8_per_specimen.csv", "esr1_regulon_used.csv")}
    body["check7"]["esr1_signature"] = ESR1_SIGNATURE
    PIO.atomic_write_json(out / ARTIFACT, json.loads(json.dumps({**meta, **body}, default=_json_default)))
    PIO.atomic_write_json(out / "_COMPLETE.json", {"path": str(out), "jobs": [os.environ.get("SLURM_JOB_ID") or "no-slurm"], "sha": meta["git_commit"],
                                                  "completed": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "produced_by": PRODUCED_BY})
    return 0


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError(f"not serialisable: {type(o)}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compartments-run", type=Path, required=True, help="the B5 run dir (id from v42_run_ids.yaml)")
    ap.add_argument("--commot-root", type=Path, required=True, help="the B6 COMMOT run dir")
    ap.add_argument("--hurdle-run", type=Path, required=True, help="the B8 run dir; its headline sweep file (A07) is the check 7 gate")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<B21 id>; never overwritten")
    ap.add_argument("--run-id", required=True, help="the B21 id")
    ap.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry(STAGE, a.protocol, run_id=a.run_id, out_root=a.out_root, comp_run=a.compartments_run, commot_root=a.commot_root, hurdle_run=a.hurdle_run)
    if a.out_root.exists():
        raise SystemExit(f"refusing: {a.out_root} exists (a run dir is never overwritten)")
    stamp = PF.stage_stamp_begin(STAGE, a.protocol, a.out_root / ARTIFACT, Path(__file__).resolve(), a.run_id)  # job-start stamp gate, before any compute
    rc = run(a.compartments_run, a.commot_root, a.hurdle_run, a.out_root, a.protocol, a.run_id)
    if rc == 0 and stamp:
        stamp.finish(comp_run=a.compartments_run, commot_root=a.commot_root, hurdle_run=a.hurdle_run)
    print("done" if rc == 0 else "gate failed", a.run_id, a.out_root, flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
