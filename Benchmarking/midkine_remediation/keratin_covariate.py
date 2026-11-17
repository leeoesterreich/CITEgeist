"""B23: which single covariate absorbs the basal-keratin vs secretory association (P4-S2 headline). Descriptive only.

Promoted from the exploratory check 1 (``exploratory_v42/check1_covariate_absorption.py``, git history). It re-runs the B12 primary test (Spearman of the basal-keratin
share ``k_share`` against the 11-gene secretory score, Moran spectral randomization p, B = 9999, seed 20261028, KNN k = 6) and its covariate-adjusted version (both
variables residualised on the covariates, then the same Spearman + surrogates) for EVERY subset of the three B12 covariates (log_spot_umi, log_rna_total (RNA content, protocol
deviation 2026-11-06; ``hurdle.PRIMARY_COVARIATES_TOTAL``), qp_prop): none (= B12 ``primary``), each alone, each pair, all three (= B12 ``partial_primary``), on the eligible cancer compartments of the B5 run, for every
specimen. The statistic, null, covariate matrix and seeds are those of ``keratin_secretory_association.py`` (``msr_spearman``, ``residualize``, ``spearman``) and
``hurdle.covariate_matrix``, called, not copied. ``adj_seed``, ``label`` and ``subsets`` are the exploratory functions unchanged (tests check AST identity).
Changed against the exploratory file: inputs come from the CLI (the B5 run and the stamped B12 run, not constants), the per-specimen work is split so a toy test can
drive it (``analyze_arrays``), a covariate-domain error is recorded as not_estimable (as B12 does) instead of raising, and the gate covers every specimen.

Seeds (unchanged): unadjusted and all-three use the B12 seed ``20261028`` (so they reproduce B12 value by value); every other subset uses
``SeedSequence([20261028, crc32(label)])``. B = 9999 surrogates, k = 6. ``canonical_eigh`` (inside ``morans.SpectralNull``) makes the eigenbasis independent of the BLAS
thread count; the sbatch pins 16 threads, as the B12 run did.

GATE (fail closed: ``gate_failure.json``, no artifact, no stamp, exit 2; stops at the first disagreement)
  Reference: the stamped B12 summary ``<keratin-run>/keratin_secretory_summary.json``.
  * config: seed, n_surrogates, k, covariates and the B5 compartments run equal the B12 config.
  * every specimen: the status (estimable or not) of the unadjusted test equals B12 ``primary``; if estimable, n, rho and p equal B12 ``primary`` and the all-three
    test equals B12 ``partial_primary`` (status, n, rho, p). rho: rtol 1e-9, atol 1e-12 (``promoted_common.gate_equal`` default, as B21/B22); p, n and n_surrogates exactly
    (rtol 0, atol 1e-12 for p). P4-S2 is analysed first, so a disagreement there stops the job before the cohort.
COMPARISON (not a gate): with ``--exploratory-results`` (check1_results.json) every (specimen, adjustment) rho and p is compared with the exploratory file; differences
  are recorded in the artifact (``exploratory_comparison``) and never stop the run. The exploratory run was gated against the SUPERSEDED B12 run v2b12 (before
  ``canonical_eigh``); its rho values are unaffected by the eigenbasis, its surrogate p-values can differ in the last digit. The stamped numbers are the ones in this artifact.

HOW THE RESULTS ARE READ (written before anything runs; from the exploratory check)
  * The single covariate that, ALONE, moves the P4-S2 rho most toward zero (largest shift of |rho| toward 0 relative to the unadjusted rho) is the leading explanation
    for the shrinkage. Nothing here is called a cause.
  * P4-S2 is compared descriptively with the min/median/max of the ESR1_WT specimens for each adjustment. One mutant specimen: no test, no p-value for that contrast.
  * Per-specimen p-values are single tests, not BH-adjusted.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import a12_10_real as R  # noqa: E402
import chain_guard as PF  # noqa: E402
import hurdle as H  # noqa: E402
import keratin_secretory_association as K  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
import protocol_io as PIO  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/keratin_covariate.py"
STAGE = "B23"
ARTIFACT = "keratin_covariate_results.json"
SIDECAR = "keratin_covariate_per_specimen.csv"
HEADLINE = K.HEADLINE
SEED = K.SEED      # 20261028, as B12
B = K.B            # 9999
RNA_COVARIATE = "total"   # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md): result keys and labels are built from these names
COVS = H.PRIMARY_COVARIATES_TOTAL
FULL = "+".join(COVS)
RHO_RTOL, RHO_ATOL, P_RTOL, P_ATOL = 1e-9, 1e-12, 0.0, 1e-12   # the tolerance of every promoted stage gate (promoted_common.gate_equal)
READING = ("Descriptive only. The leading explanation for the shrinkage is the single covariate that moves the P4-S2 rho most toward zero; this is not a causal claim. "
           "P4-S2 is one ESR1-D538G specimen from one patient: its comparison with the ESR1_WT range carries no test and no p-value. Per-specimen p-values are single tests, not BH-adjusted.")


# --------------------------------------------------------------------------- exploratory functions, unchanged
def adj_seed(lab: str) -> int:
    """Order-independent stream per adjustment. ``unadjusted`` and the all-three set use the v2b12 seed (so they reproduce
    v2b12, which calls msr_spearman with seed=SEED for both); every other adjustment gets a seed derived from (SEED, label)."""
    if lab in ("unadjusted", FULL):
        return SEED
    return int(np.random.SeedSequence([SEED, zlib.crc32(lab.encode())]).generate_state(1)[0])


def label(subset) -> str:
    return "unadjusted" if not subset else "+".join(subset)


def subsets():
    out = [()]
    for r in (1, 2, 3):
        out += list(itertools.combinations(COVS, r))
    return out


# --------------------------------------------------------------------------- one specimen
def analyze_arrays(spec: str, kshare, sec, xy, cov_matrix, B: int = B, seed: int = SEED) -> dict:
    """All adjustments for one specimen. ``cov_matrix(kept)`` returns the standardised covariate columns for the tuple ``kept`` (in COVS order);
    it may raise ``H.CovariateDomainError`` (recorded as not_estimable for every adjusted test, as B12 does for its partial tests)."""
    kshare, sec = np.asarray(kshare, float), np.asarray(sec, float)
    n = len(kshare)
    if n < K.MIN_EL or n > M.DENSE_EIGH_MAX_N or np.std(kshare) == 0 or np.std(sec) == 0:
        return {"specimen": spec, "status": "not_estimable", "n": int(n)}
    sn = M.SpectralNull(knn_weights(np.asarray(xy, float), K.K))
    try:
        full, cov_err = cov_matrix(tuple(COVS)), ""
    except H.CovariateDomainError as e:
        full, cov_err = None, str(e)
    corr = None
    if full is not None:
        corr = {c: {"rho_with_kshare": K.spearman(full[:, i], kshare), "rho_with_secretory": K.spearman(full[:, i], sec)} for i, c in enumerate(COVS)}
    res = {}
    for sub in subsets():
        lab = label(sub)
        if not sub:
            r = K.msr_spearman(kshare, sec, sn, B=B, seed=seed)      # seed == adj_seed('unadjusted')
        elif full is None:
            res[lab] = {"status": "not_estimable", "reason": f"covariate_domain: {cov_err}"}
            continue
        else:
            Z = cov_matrix(tuple(sub))
            rk, rs = K.residualize(kshare, Z), K.residualize(sec, Z)
            if min(rk.std(), rs.std()) < 1e-12:
                res[lab] = {"status": "not_estimable", "reason": "constant residual"}
                continue
            r = K.msr_spearman(rk, rs, sn, B=B, seed=seed if lab == FULL else adj_seed(lab))
        res[lab] = {"status": "ok", "rho": r["rho"], "p_two_sided": r["p_two_sided"], "n_surrogates": r["n_surrogates"]}
    return {"specimen": spec, "status": "ok", "n": int(n), "adjustments": res, "covariate_spearman": corr}


def analyze(spec: str, comp_run, B: int = B, seed: int = SEED) -> dict:
    """The real-data specimen: eligible cancer compartments of the B5 run (the loader B12 uses), ``k_share`` and the 11-gene secretory score."""
    if not (Path(comp_run) / f"{spec}.h5ad").is_file():
        return {"specimen": spec, "status": "not_estimable", "n": 0, "reason": "missing_source"}
    d = R.load_cgeist_eligible(Path(comp_run), spec)
    obs, genes = d["obs"], list(d["genes"])
    sec, _ = R._mean_z(d["X"], genes, R.SECRETORY)
    return analyze_arrays(spec, obs["k_share"].to_numpy(float), sec, d["xy"],
                          lambda kept: H.covariate_matrix(obs, drop=tuple(c for c in COVS if c not in kept), rna=RNA_COVARIATE), B=B, seed=seed)


# --------------------------------------------------------------------------- gate against the stamped B12 summary
def _flag(name: str, got, expected) -> dict:
    ok = got == expected
    rec = {"name": name, "got": got, "expected": expected, "match": bool(ok)}
    if not ok:
        raise CM.GateFailure(f"GATE FAILED {name}: got {got!r} expected {expected!r}")
    return rec


def _gate_test(spec: str, what: str, rec: dict, st: dict) -> list:
    """One of our tests (``rec`` = our adjustment record) against one stored B12 test (``st``): status, then rho and p."""
    out = [_flag(f"{spec} {what} status", rec.get("status"), st.get("status"))]
    if st.get("status") == "ok":
        out += [CM.gate_equal(f"{spec} {what} rho", rec["rho"], st["rho"], rtol=RHO_RTOL, atol=RHO_ATOL),
                CM.gate_equal(f"{spec} {what} p", rec["p_two_sided"], st["p_two_sided"], rtol=P_RTOL, atol=P_ATOL),
                _flag(f"{spec} {what} n_surrogates", rec["n_surrogates"], st["n_surrogates"])]
    return out


def gate_specimen(rec: dict, stored: dict) -> list:
    """``unadjusted`` = B12 ``primary`` and the all-three set = B12 ``partial_primary`` for this specimen; raises CM.GateFailure on any disagreement."""
    spec = rec["specimen"]
    p, pp = stored["per_specimen"]["primary"].get(spec), stored["per_specimen"]["partial_primary"].get(spec)
    if p is None or pp is None:
        raise CM.GateFailure(f"GATE FAILED {spec}: the stamped B12 summary has no primary / partial_primary row for this specimen")
    if p.get("status") != "ok":
        return [_flag(f"{spec} estimability (B12 primary not estimable)", rec["status"], "not_estimable")]
    out = [_flag(f"{spec} estimability", rec["status"], "ok"), _flag(f"{spec} n", rec["n"], p["n"])]
    out += _gate_test(spec, "unadjusted", rec["adjustments"]["unadjusted"], p)
    out += _gate_test(spec, "all-three", rec["adjustments"][FULL], pp)
    return out


def gate_config(stored: dict, comp_run, B_run: int, seed: int) -> list:
    c = stored["config"]
    return [_flag("config seed", seed, c["seed"]), _flag("config n_surrogates", B_run, c["n_surrogates"]), _flag("config k", K.K, c["k"]),
            _flag("config covariates", list(COVS), list(c["covariates_partial"])),
            _flag("config compartments_run", str(Path(comp_run).resolve()), str(Path(c["compartments_run"]).resolve()))]


# --------------------------------------------------------------------------- summary and comparison
def summarise(rows: dict, genotype: dict, headline: str = HEADLINE) -> dict:
    """P4-S2 vs ESR1_WT range per adjustment, and the shift toward zero for each adjustment (the exploratory ``summarise``)."""
    wt = [s for s, r in rows.items() if r["status"] == "ok" and genotype.get(s) == "ESR1_WT"]
    adjs = list(rows[headline]["adjustments"])
    h = rows[headline]["adjustments"]
    unadj = h["unadjusted"]["rho"]
    table = {}
    for a in adjs:
        v = [rows[s]["adjustments"][a]["rho"] for s in wt if rows[s]["adjustments"][a]["status"] == "ok"]
        rho = h[a]["rho"]
        table[a] = {"P4_S2_rho": rho, "P4_S2_p": h[a]["p_two_sided"],
                    "shift_toward_zero_fraction": float((abs(unadj) - abs(rho)) / abs(unadj)),
                    "wt_n": len(v), "wt_min": float(np.min(v)), "wt_median": float(np.median(v)), "wt_max": float(np.max(v)),
                    "P4_S2_inside_wt_range": bool(np.min(v) <= rho <= np.max(v))}
    singles = [a for a in adjs if "+" not in a and a != "unadjusted"]
    lead = max(singles, key=lambda a: table[a]["shift_toward_zero_fraction"])
    return {"p4s2_vs_wt_by_adjustment": table,
            "single_covariate_shift_ranking": sorted(singles, key=lambda a: -table[a]["shift_toward_zero_fraction"]), "leading_single_covariate": lead,
            "interpretation_rule": "leading explanation = the single covariate that moves P4-S2 rho most toward zero; not a causal claim"}


def compare_exploratory(rows: dict, ref: dict) -> dict:
    """Not a gate. Every (specimen, adjustment) present in both: rho at rtol 1e-9 / atol 1e-12, p exactly. Returns the counts and the list of differences."""
    diffs, n = [], 0
    for s, r in rows.items():
        er = ref.get("per_specimen", {}).get(s)
        if r["status"] != "ok" or not er or er.get("status") != "ok":
            continue
        for a, mine in r["adjustments"].items():
            theirs = er["adjustments"].get(a)
            if mine["status"] != "ok" or not theirs or theirs.get("status") != "ok":
                continue
            n += 1
            for f, tol in (("rho", dict(rtol=RHO_RTOL, atol=RHO_ATOL)), ("p_two_sided", dict(rtol=P_RTOL, atol=P_ATOL))):
                if not np.isclose(mine[f], theirs[f], **tol):
                    diffs.append({"specimen": s, "adjustment": a, "field": f, "stage": mine[f], "exploratory": theirs[f], "abs_diff": float(abs(mine[f] - theirs[f]))})
    lead_ref = (ref.get("summary") or {}).get("leading_single_covariate")
    return {"non_gating": True, "n_compared": n, "n_rho_differ": sum(d["field"] == "rho" for d in diffs), "n_p_differ": sum(d["field"] == "p_two_sided" for d in diffs),
            "max_abs_rho_diff": max([d["abs_diff"] for d in diffs if d["field"] == "rho"], default=0.0),
            "max_abs_p_diff": max([d["abs_diff"] for d in diffs if d["field"] == "p_two_sided"], default=0.0),
            "exploratory_reference_run": ref.get("reference"), "exploratory_git_commit": ref.get("git_commit"), "exploratory_leading_single_covariate": lead_ref,
            "tolerance": {"rho_rtol": RHO_RTOL, "rho_atol": RHO_ATOL, "p_rtol": P_RTOL, "p_atol": P_ATOL},
            "note": "The exploratory run was gated against the superseded B12 run v2b12 (before canonical_eigh): rho does not depend on the eigenbasis, surrogate p can differ. The stamped values are the ones in this artifact.",
            "differences": diffs}


def flat_rows(rows: dict, genotype: dict) -> list:
    out = []
    for s, r in rows.items():
        if r["status"] != "ok":
            continue
        for a, v in r["adjustments"].items():
            out.append({"specimen": s, "genotype": genotype.get(s, ""), "n_eligible": r["n"], "adjustment": a, "status": v["status"],
                        "rho": v.get("rho", ""), "p_two_sided": v.get("p_two_sided", ""), "n_surrogates": v.get("n_surrogates", ""), "reason": v.get("reason", "")})
    return out


def specimens_of(protocol: dict) -> list:
    return sorted(protocol["specimens"]["ordered"], key=lambda s: (s != HEADLINE, s))


def compute(comp_run, stored: dict, protocol: dict, specimens=None, B_run: int = B, seed: int = SEED, analyze_fn=None, log=print) -> dict:
    """Analyse every specimen (P4-S2 first) and gate each against the stored B12 summary right after it is computed. Raises ``CM.GateFailure``; the exception's
    ``record`` attribute holds the failing specimen. No I/O besides reading the inputs."""
    analyze_fn = analyze_fn or analyze
    specs = list(specimens) if specimens else specimens_of(protocol)
    genotype = protocol["specimens"]["genotype"]
    gates = gate_config(stored, comp_run, B_run, seed)
    rows = {}
    for s in specs:
        t = time.time()
        rows[s] = analyze_fn(s, comp_run, B=B_run, seed=seed)
        try:
            gates += gate_specimen(rows[s], stored)
        except CM.GateFailure as e:
            e.record = rows[s]
            raise
        log(f"{s}: {rows[s]['status']} n={rows[s].get('n')} {time.time()-t:.0f}s  gate passed")
    body = {"gates": gates, "per_specimen": rows, "per_specimen_table": flat_rows(rows, genotype)}
    if HEADLINE in rows and rows[HEADLINE]["status"] == "ok":
        body["summary"] = summarise(rows, genotype)
    return body


# --------------------------------------------------------------------------- stage entry
def sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def config() -> dict:
    return {"seed": SEED, "seed_rule": "unadjusted and all-three: the B12 seed; every other subset SeedSequence([seed, crc32(label)])", "n_surrogates": B, "k": K.K,
            "covariates": list(COVS), "rna_covariate": RNA_COVARIATE, "protocol_deviation": H.deviation_record(), "adjustments": [label(s) for s in subsets()], "min_eligible": K.MIN_EL,
            "gate_tolerance": {"rho_rtol": RHO_RTOL, "rho_atol": RHO_ATOL, "p_rtol": P_RTOL, "p_atol": P_ATOL}, "reading_rule": READING, "descriptive_only": True}


def run(comp_run, keratin_run, out_root, protocol_path, run_id: str = "", specimens=None, B_run: int = B, exploratory_results=None, meta_extra=None) -> int:
    """Compute, write ``keratin_covariate_results.json`` + csv + ``_COMPLETE.json``. 0 ok, 2 gate failure (``gate_failure.json``, nothing else)."""
    protocol = PIO.load_protocol(protocol_path)
    out = Path(out_root)
    out.mkdir(parents=True, exist_ok=True)
    code = HERE.parents[1]
    git = lambda *a: subprocess.run(["git", "-C", str(code), *a], capture_output=True, text=True).stdout.strip()  # noqa: E731
    ref_file = Path(keratin_run) / "keratin_secretory_summary.json"
    stored = json.loads(ref_file.read_text())
    meta = {"stage": STAGE, "run_id": run_id, "produced_by": PRODUCED_BY, "descriptive_only": True, "config": config(), "compartments_run": str(comp_run),
            "keratin_run": str(keratin_run), "reference_b12_sha256": sha256(ref_file),
            "blas_threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, "git_commit": git("rev-parse", "HEAD"), **(meta_extra or {})}
    try:
        body = compute(comp_run, stored, protocol, specimens, B_run, log=lambda m: print(m, flush=True))
    except CM.GateFailure as e:
        PIO.atomic_write_json(out / "gate_failure.json", {**meta, "gate_failure": str(e), "failing_specimen": getattr(e, "record", None)})
        print(str(e), file=sys.stderr, flush=True)
        return 2
    if exploratory_results:
        ep = Path(exploratory_results)
        ref = json.loads(ep.read_text())
        body["exploratory_comparison"] = {"file": str(ep), "sha256": sha256(ep), **compare_exploratory(body["per_specimen"], ref)}
    else:
        body["exploratory_comparison"] = {"non_gating": True, "skipped": "no --exploratory-results given"}
    with open(out / SIDECAR, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["specimen", "genotype", "n_eligible", "adjustment", "status", "rho", "p_two_sided", "n_surrogates", "reason"])
        w.writeheader()
        w.writerows(body["per_specimen_table"])
    meta["sidecar_files_sha256"] = {SIDECAR: sha256(out / SIDECAR)}
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
    ap.add_argument("--keratin-run", type=Path, required=True, help="the stamped B12 run dir; its keratin_secretory_summary.json is the gate reference")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<B23 id>; never overwritten")
    ap.add_argument("--run-id", required=True, help="the B23 id")
    ap.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    ap.add_argument("--exploratory-results", type=Path, default=None, help="check1_results.json of the exploratory check 1; compared, never gated (path and sha256 go into the artifact)")
    a = ap.parse_args(argv)
    PF.guard_entry(STAGE, a.protocol, run_id=a.run_id, out_root=a.out_root, comp_run=a.compartments_run, keratin_run=a.keratin_run)
    if a.out_root.exists():
        raise SystemExit(f"refusing: {a.out_root} exists (a run dir is never overwritten)")
    if a.exploratory_results is not None and not a.exploratory_results.is_file():
        raise SystemExit(f"--exploratory-results {a.exploratory_results} does not exist")
    stamp = PF.stage_stamp_begin(STAGE, a.protocol, a.out_root / ARTIFACT, Path(__file__).resolve(), a.run_id)  # job-start stamp gate, before any compute
    rc = run(a.compartments_run, a.keratin_run, a.out_root, a.protocol, a.run_id, exploratory_results=a.exploratory_results)
    if rc == 0 and stamp:
        stamp.finish(comp_run=a.compartments_run, keratin_run=a.keratin_run)
    print("done" if rc == 0 else "gate failed", a.run_id, a.out_root, flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
