"""B1 (rev 7): coverage of the nominal 95% GMR interval through the EXACT PRODUCTION path ``descriptive.gmr_interval``.

B24 (protocol deviation 2026-11-06, ``PROTOCOL_DEVIATION_2026-11-06_rna_total.md``): this module is also the stamped stage B24 (run ``v2b24``, ``main``). It reruns the same four
cases with the production covariates under ``rna="total"`` (``RNA_COVARIATE``): in the synthetic worlds ``rna_total`` is the simulated other-gene allocated total ``E_o`` (obs column
``rna_total``; ``hurdle.rna_total(obs)`` is gated against it for every dataset). The stamped artifact ``prod_coverage_total.json`` is a stamped INPUT of B8 (the coverage the hurdle
records). The pre-stated criterion below is reported per case against the binding historical values, never a gate of the stage; the per-cell series (0.884 ...) stays on disk
as the per-cell estimator's record.

The stored numbers 0.874/0.896/0.898/0.900 (synthetic/g6_fix/boot_*_s6.json) came from the historical ``block_bootstrap`` (ridge,
> 50% finite rule, 500 resamples). This module reruns the same four cases (same lattice, same dataset seeds 1000..1499, same
frozen k_share threshold, same 6-pitch tiles, 500 datasets) with the production estimator at the FROZEN production settings:
2000 resamples, the production bootstrap seed ``descriptive.BOOT_SEED`` for every dataset (production uses one fixed seed), 0.95 valid
fraction, production covariates (``hurdle.covariate_matrix``, strict-positive domain gate, 3 primary covariates).

PRE-STATED REPRODUCTION CRITERION (fixed before the run): the stored values reproduce iff, for ALL FOUR cases,
  (i)   |new coverage - stored coverage| <= 0.03   (about one binomial half-width at n = 500, p = 0.88: 1.96*sqrt(.88*.12/500) = 0.029),
  (ii)  the new coverage lies inside the stored Wilson 95% interval, and
  (iii) at most 25 of the 500 datasets (5%) are typed not_estimable.
New coverage = covered / estimable datasets (the same population definition as the stored n_ok); the conservative figure that
counts a not_estimable dataset as a miss is reported beside it.

HELD-OUT PROSPECTIVE RUN (rev 8, user decision): the first production-path run (array 11956368, seeds 1000..1499) STARTED BEFORE the
commit that fixed the criterion above (b3ab5a97), so it is exploratory only. The held-out run uses the SAME four cases, 500 datasets,
production estimator and settings and the SAME criterion, but dataset seeds HELDOUT_SEED0..HELDOUT_SEED0+499 = 70000000..70000499,
disjoint from every seed used earlier (synthetic lattice datasets 1000..1999 and their +7 offsets, injection 5000..5049 (+7), k0
reference draws 900000..900019, G7 seeds below 7000 + 100000*(n_lattices*9) + 1000, shift pools 20261015/20261016, bootstrap seed
20261014). The criterion is compared with the BINDING values 0.874/0.896/0.898/0.900 (``stored_record``), not with the exploratory
production-path values. This file, run_prod_coverage.py --heldout and slurm/prod_coverage_heldout.sbatch are committed BEFORE the job
is submitted; every task start time reported by SLURM must be later than that commit's committer time. If any case fails the criterion
the binding numbers are untouched and the decision returns to the user.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import descriptive as D
import hurdle as H
import synthetic_g6 as G
import synthetic_g6fix as F

CASES = ("null_autocorrelated", "null_autocorrelated_stress", "null_depth_confounded", "null_noisy_nst_kshare")
N_DATASETS = 500
SEED0 = 1000
HELDOUT_SEED0 = 70_000_000
TOL = 0.03
MAX_NOT_ESTIMABLE = 25
RNA_COVARIATE = "total"        # protocol deviation 2026-11-06: the covariate the production hurdle uses
ARTIFACT = "prod_coverage_total.json"
PRODUCED_BY = "Benchmarking/midkine_remediation/src/prod_coverage.py"
_C = {}


def _init(lat):
    _C["lat"] = lat
    G.K0_SHARE = F.K0_SHARE_FROZEN


def _one(args):
    case, seed = args
    lat = _C["lat"]
    sim = G.simulate(case, lat, np.random.default_rng(seed), inject_fold=1.0)
    obs = sim["obs"]
    if not np.allclose(H.rna_total(obs), obs["rna_total"].to_numpy(float), rtol=1e-9, atol=0.0):      # the simulator's E_o IS the total: fail closed if the two definitions drift
        raise RuntimeError(f"{case} seed {seed}: hurdle.rna_total(obs) differs from the simulated other-gene total E_o")
    try:
        cov = H.covariate_matrix(obs, rna=RNA_COVARIATE)
    except H.CovariateDomainError as e:
        return {"status": "not_estimable", "reason_code": "covariate_domain", "reason": str(e)}
    r = D.gmr_interval(sim["scores"][0], sim["group"], cov, lat.blocks)       # production defaults (2000 resamples, production seed, 0.95 valid fraction)
    return {k: r.get(k) for k in ("status", "reason_code", "lo", "hi", "valid_fraction", "n_valid", "n_requested", "seed")}


def summarise(res, truth: float = 1.0) -> dict:
    ok = [r for r in res if r["status"] == "ok"]
    k = sum(1 for r in ok if r["lo"] <= truth <= r["hi"])
    n_ok = len(ok)
    vf = [r["valid_fraction"] for r in res if r.get("valid_fraction") is not None]
    reasons = {}
    for r in res:
        if r["status"] != "ok":
            reasons[r["reason_code"]] = reasons.get(r["reason_code"], 0) + 1
    return {
        "n_datasets": len(res), "n_estimable": n_ok, "n_not_estimable": len(res) - n_ok, "not_estimable_reasons": reasons,
        "coverage": k / n_ok if n_ok else float("nan"), "wilson95": list(G.wilson(k, n_ok)) if n_ok else None,
        "coverage_not_estimable_as_miss": k / len(res), "wilson95_not_estimable_as_miss": list(G.wilson(k, len(res))),
        "valid_fraction": {"min": float(min(vf)), "median": float(np.median(vf)), "mean": float(np.mean(vf)), "n_below_0.95": int(sum(v < 0.95 for v in vf)), "n": len(vf)} if vf else None,
    }


def reproduces(new: dict, stored: dict) -> dict:
    """Apply the pre-stated criterion to one case. ``new`` = summarise(); ``stored`` = {coverage, ci95_wilson}."""
    d = abs(new["coverage"] - stored["coverage"])
    lo, hi = stored["ci95_wilson"]
    c = {"abs_diff": d, "i_within_tolerance": bool(d <= TOL), "ii_inside_stored_wilson": bool(lo <= new["coverage"] <= hi),
         "iii_not_estimable_ok": bool(new["n_not_estimable"] <= MAX_NOT_ESTIMABLE)}
    c["reproduces"] = all(c[k] for k in ("i_within_tolerance", "ii_inside_stored_wilson", "iii_not_estimable_ok"))
    return c


def run_case(case, n=N_DATASETS, procs=1, lat=None, seed0=SEED0):
    lat = lat or G.make_lattice()
    G.K0_SHARE = F.K0_SHARE_FROZEN
    jobs = [(case, seed0 + i) for i in range(n)]
    t0 = time.time()
    if procs > 1:
        with Pool(procs, initializer=_init, initargs=(lat,)) as p:
            res = p.map(_one, jobs, chunksize=4)
    else:
        _init(lat)
        res = [_one(j) for j in jobs]
    out = summarise(res)
    out.update({"case": case, "tile_side_pitches": 6.0, "fold": 1.0, "seconds": time.time() - t0,
                "settings": {"B": D.BOOT_N, "seed_rule": "production fixed seed for every dataset", "seed": D.BOOT_SEED, "min_valid_fraction": D.MIN_VALID_FRACTION,
                             "dataset_seeds": f"{seed0}..{seed0 + n - 1}", "covariates": list(H.primary_covariates(RNA_COVARIATE)), "rna_covariate": RNA_COVARIATE,
                             "k0_share": F.K0_SHARE_FROZEN}})
    return out


def stored_record(case, root: Path) -> dict:
    j = json.loads((Path(root) / f"boot_{case}_s6.json").read_text())["result"]["gmr"]
    return {"coverage": j["coverage"], "ci95_wilson": j["ci95_wilson"], "n_ok": j["n_ok"]}


# --------------------------------------------------------------------------- B24: the stamped stage
LATTICE_SPECIMEN = "HCC22-088-P4-S2_1i_rep"


def lattice_inputs(namespace_root=None) -> dict:
    """The untracked files the synthetic lattice is built from, resolved to the NAMESPACE ROOT (``V2_NAMESPACE_ROOT``, the main checkout) and not to the clean code worktree
    (same mechanism as B19 / B20, dff892913): ``qp_dir`` (the protocol's inputs.synthetic_lattice_dir, ``<root>/output/module3_marker_share_7type``), the QP csv of the lattice specimen and the Visium
    tissue_positions.csv. With no namespace root (unit tests, legacy) the module defaults are used. Raises FileNotFoundError if a file is missing."""
    import os

    ns = namespace_root or os.environ.get("V2_NAMESPACE_ROOT")
    if ns:
        rel = json.loads((Path(__file__).resolve().parents[1] / "protocol.json").read_text())["inputs"]["synthetic_lattice_dir"]
        qp_dir = Path(ns) / rel
    else:
        import spatial_utils as SU

        qp_dir = SU.QP_DIR
    import spatial_utils as SU

    files = {"qp_csv": qp_dir / f"{LATTICE_SPECIMEN}_cell_prop_finetuned_results.csv",
             "tissue_positions": SU.DATA_ROOT / LATTICE_SPECIMEN / "outs" / "spatial" / "tissue_positions.csv"}
    for k, f in files.items():
        if not f.is_file():
            raise FileNotFoundError(f"B24 lattice input {k} not found: {f}")
    import hashlib

    return {"qp_dir": qp_dir, "sha256": {k: hashlib.sha256(f.read_bytes()).hexdigest() for k, f in files.items()}, "files": {k: str(f) for k, f in files.items()}}


def stage_artifact(per_case: dict, run_id: str, fix_root: Path, n: int, seed0: int, lattice: dict = None) -> dict:
    """The B24 artifact: per case the ``run_case`` summary, the criterion against the binding historical values (reported, not a gate), the coverage per case in ``CASES`` order
    (``coverage_by_case``: what B8 reads), and the per-cell production-path series beside it."""
    cases = list(per_case)
    if cases != list(CASES):
        raise ValueError(f"B24 needs exactly the cases {list(CASES)} in order, got {cases}")
    for c, r in per_case.items():
        if not (r["n_estimable"] > 0 and np.isfinite(r["coverage"])):
            raise RuntimeError(f"B24 {c}: no estimable dataset (coverage {r['coverage']})")
    stored = {c: stored_record(c, fix_root) for c in cases}
    return {"stage": "B24", "run_id": run_id, "rna_covariate": RNA_COVARIATE,
            "unstamped_external_inputs_sha256": (lattice or {}).get("sha256"), "unstamped_external_inputs_files": (lattice or {}).get("files"),
            "unstamped_external_inputs_note": "the QP csv (untracked, path token marker_share: not a stamp input) and the Visium tissue_positions.csv the synthetic lattice is built from; read from the namespace root", "protocol_deviation": H.deviation_record(),
            "covariates": list(H.primary_covariates(RNA_COVARIATE)), "n_datasets": int(n), "dataset_seed0": int(seed0), "cases": cases,
            "coverage_by_case": {c: per_case[c]["coverage"] for c in cases},
            "per_case": per_case, "criterion_vs_binding_historical": {c: reproduces(per_case[c], stored[c]) for c in cases}, "stored_historical": stored,
            "per_cell_production_path_coverage": dict(D.PRODUCTION_PATH_COVERAGE),
            "criterion_text": __doc__.split("PRE-STATED", 1)[1].split("HELD-OUT", 1)[0].strip(),
            "note": "Reported as it comes out; the stage gates nothing on the coverage values. The B8 hurdle records coverage_by_case as its measured synthetic coverage."}


def main(argv=None) -> int:
    """B24 stage entry (sbatch): all four cases at the production estimator with rna='total', one stamped artifact. Creates a NEW run dir (never overwrites)."""
    HERE = Path(__file__).resolve().parent
    sys.path.insert(0, str(HERE.parent))
    sys.path.insert(0, str(HERE.parents[2] / "scripts"))
    import chain_guard as PF
    import protocol_io as PIO

    ap = argparse.ArgumentParser(description=main.__doc__)
    ap.add_argument("--run-id", required=True, help="the B24 id in v42_run_ids.yaml")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<B24 id>")
    ap.add_argument("--n", type=int, default=N_DATASETS)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry("B24", a.protocol, run_id=a.run_id, out_root=a.out_root)
    if a.out_root.exists():
        raise SystemExit(f"refusing: {a.out_root} exists (a run dir is never overwritten)")
    stamp = PF.stage_stamp_begin("B24", a.protocol, a.out_root / ARTIFACT, Path(__file__).resolve(), a.run_id)    # job-start stamp gate, before any compute
    lattice = lattice_inputs()
    lat = G.make_lattice(LATTICE_SPECIMEN, qp_dir=lattice["qp_dir"])
    per_case = {c: run_case(c, a.n, a.procs, lat=lat) for c in CASES}
    for c, r in per_case.items():
        print(f"{c}: coverage {r['coverage']:.3f} estimable {r['n_estimable']}/{r['n_datasets']} {r['seconds']:.0f}s", flush=True)
    art = stage_artifact(per_case, a.run_id, HERE.parent / "synthetic" / "g6_fix", a.n, SEED0, lattice)
    a.out_root.mkdir(parents=True)
    PIO.atomic_write_json(a.out_root / ARTIFACT, art)
    if stamp:
        stamp.finish()
    print("done B24", a.run_id, a.out_root, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
