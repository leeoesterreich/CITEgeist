"""Basal-keratin share vs secretory score / MDK: per-specimen Spearman with a Moran spectral randomization null.

(Rules fixed before any real result was seen).

Reuses the v41 pipeline pieces unchanged: eligible cancer compartments (``compartments.eligible_cancer``), the 11-gene
secretory score (``a12_10_real._mean_z``), ``obs.k_share`` (the h3 keratin variable), KNN k=6 weights and the
``morans.SpectralNull`` eigen-decomposition. The null surrogates the rank-transformed outcome field with sign-flipped
spectral components and recomputes Spearman rho on ranks; the other variable is held fixed.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import a12_10_real as R  # noqa: E402
import chain_guard as PF  # noqa: E402
import hurdle as H  # noqa: E402
import morans as M  # noqa: E402
import protocol_io as PIO  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/keratin_secretory_association.py"
ANALYSIS_ID = "ks1"
SEED = 20261028
B = 9999
K = R.K
MIN_EL = R.MIN_EL
MIN_NZ = R.MIN_NZ
ALPHA = R.ALPHA
SECRETORY = R.SECRETORY
HEADLINE = "HCC22-088-P4-S2_1i_rep"
TESTS = ("primary", "secondary", "partial_primary", "partial_secondary")
DEFAULT_COMP_RUN = "output/compartments_module3_7type/ec0cf5b80715/c1"
RNA_COVARIATE = "total"        # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md)


# --------------------------------------------------------------------------- statistic and null
def spearman(x, y) -> float:
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    return float(rx @ ry / np.sqrt((rx @ rx) * (ry @ ry)))


def msr_spearman(fixed, field, sn: M.SpectralNull, B: int = B, seed: int = SEED, chunk: int = 500) -> dict:
    """Spearman rho(fixed, field) with a two-sided MSR p-value; ``field`` is surrogated (rank-transformed), ``fixed`` is not."""
    n = sn.n
    rx = rankdata(fixed)
    zx = (rx - rx.mean()) / rx.std()
    ry = rankdata(field)
    yc = ry - ry.mean()
    rho = float(zx @ yc / (n * yc.std()))
    c = sn.V.T @ yc
    rng = np.random.default_rng(seed)
    hits = 0
    done = 0
    while done < B:
        m = min(chunk, B - done)
        S = rng.integers(0, 2, size=(n, m), dtype=np.int8) * 2.0 - 1.0
        Z = sn.V @ (S * c[:, None])
        Rk = rankdata(Z, axis=0)
        Rc = Rk - Rk.mean(axis=0)
        rstar = (zx @ Rc) / (n * Rk.std(axis=0))
        hits += int(np.sum(np.abs(rstar) >= abs(rho) - 1e-12))
        done += m
    return {"rho": rho, "p_two_sided": float((1 + hits) / (B + 1)), "n_surrogates": int(B)}


def residualize(v, Z) -> np.ndarray:
    """OLS residuals of v on [1, Z]."""
    X = np.column_stack([np.ones(len(v)), Z])
    beta, *_ = np.linalg.lstsq(X, np.asarray(v, float), rcond=None)
    return np.asarray(v, float) - X @ beta


def _bh(pvals):
    return R._bh(pvals) if len(pvals) else []


def _nonest(code: str, reason: str, **extra) -> dict:
    return {"status": "not_estimable", "reason_code": code, "reason": reason, **extra}


# --------------------------------------------------------------------------- one specimen
def analyze_specimen(spec: str, kshare, X, genes, xy, covariates=None, cov_error: str = "", B: int = B, seed: int = SEED) -> dict:
    """All four tests for one specimen's eligible compartments. ``X`` = N_P (n, genes); ``covariates`` = standardised (n, p) or None."""
    n = len(kshare)
    out = {"specimen": spec, "n_eligible": int(n)}
    gate = None
    if n < MIN_EL:
        gate = _nonest("eligibility_lt30", f"{n} eligible compartments < {MIN_EL}", n_eligible=n)
    elif n > M.DENSE_EIGH_MAX_N:
        gate = _nonest("exceeds_eigh", f"{n} eligible compartments > dense-eigh ceiling {M.DENSE_EIGH_MAX_N}", n_eligible=n)
    elif np.std(kshare) == 0:
        gate = _nonest("constant_keratin", "keratin share has zero variance", n_eligible=n)
    if gate is None:
        sec, nsec = R._mean_z(X, list(genes), SECRETORY)
        if np.std(sec) == 0:
            gate = _nonest("constant_secretory", "secretory score has zero variance", n_eligible=n)
    if gate is not None:
        return {**out, **{t: dict(gate) for t in TESTS}}
    mdk = R._column(X, list(genes), "MDK")
    nz = int((mdk > 0).sum())
    sn = M.SpectralNull(knn_weights(np.asarray(xy, float), K))

    def run(fixed, field):
        r = msr_spearman(fixed, field, sn, B=B, seed=seed)
        return {"status": "ok", "reason": "", "reason_code": "", "value": r["rho"], **r, "n": int(n)}

    out["n_secretory_genes"] = nsec
    out["mdk_nonzero"] = nz
    out["primary"] = run(kshare, sec)
    sec_gate = None
    if nz < MIN_NZ:
        sec_gate = _nonest("mdk_nonzero_lt10", f"{nz} compartments with MDK>0 < {MIN_NZ}", mdk_nonzero=nz)
    out["secondary"] = sec_gate or run(mdk, kshare)
    if covariates is None:
        pg = _nonest("covariate_domain", cov_error or "covariates unavailable")
        out["partial_primary"], out["partial_secondary"] = dict(pg), dict(pg)
    else:
        rk, rs = residualize(kshare, covariates), residualize(sec, covariates)
        out["partial_primary"] = _nonest("constant_residual", "residualized variable is constant") if min(rk.std(), rs.std()) < 1e-12 else run(rk, rs)
        if sec_gate:
            out["partial_secondary"] = dict(sec_gate)
        else:
            rm = residualize(mdk, covariates)
            out["partial_secondary"] = _nonest("constant_residual", "residualized variable is constant") if min(rm.std(), rk.std()) < 1e-12 else run(rm, rk)
    return out


# --------------------------------------------------------------------------- cohort assembly
def _patient(spec: str) -> str:
    m = re.search(r"-(P\d+)-", spec)
    return m.group(1) if m else ""


def assemble(per_spec: dict, genotype: dict, headline: str = HEADLINE) -> dict:
    """BH within each test family over estimable specimens, cohort counts, focal block, descriptive mutant-vs-WT."""
    tests = {}
    cohort = {}
    for t in TESTS:
        ok = [s for s, r in per_spec.items() if r[t]["status"] == "ok"]
        q = _bh([per_spec[s][t]["p_two_sided"] for s in ok])
        for s, qi in zip(ok, q):
            per_spec[s][t]["q"] = float(qi)
        rows = {s: per_spec[s][t] for s in per_spec}
        rhos = [per_spec[s][t]["rho"] for s in ok]
        sig = [s for s in ok if per_spec[s][t]["q"] < ALPHA]
        pos = [s for s in ok if per_spec[s][t]["rho"] > 0]
        cohort[t] = {"n_specimens": len(per_spec), "n_estimable": len(ok), "n_not_estimable": len(per_spec) - len(ok),
                     "n_positive": len(pos), "n_significant_q_lt_0.05": len(sig), "n_positive_and_significant": len([s for s in sig if s in pos]),
                     "n_negative_and_significant": len([s for s in sig if s not in pos]), "significant_specimens": sig,
                     "mean_rho": float(np.mean(rhos)) if rhos else None, "median_rho": float(np.median(rhos)) if rhos else None}
        tests[t] = rows
    wt = [s for s in per_spec if genotype.get(s) == "ESR1_WT"]
    desc = {"note": "DESCRIPTIVE only: one ESR1-D538G specimen from one patient; no test, no p-value for this contrast.", "mutant_specimen": headline}
    for t in TESTS:
        wt_rho = [per_spec[s][t]["rho"] for s in wt if per_spec[s][t]["status"] == "ok"]
        h = per_spec.get(headline, {}).get(t, {})
        desc[t] = {"mutant_rho": h.get("rho") if h.get("status") == "ok" else None, "mutant_status": h.get("status"),
                   "n_wt_estimable": len(wt_rho), "wt_median_rho": float(np.median(wt_rho)) if wt_rho else None,
                   "wt_min_rho": float(np.min(wt_rho)) if wt_rho else None, "wt_max_rho": float(np.max(wt_rho)) if wt_rho else None}
    focal = {t: {"specimen": headline, **tests[t].get(headline, {"status": "missing"})} for t in TESTS}
    return {"focal_P4_S2": focal, "cohort": cohort, "descriptive_mutant_vs_wt": desc, "per_specimen": tests}


def to_rows(per_spec: dict, genotype: dict) -> list:
    rows = []
    for s, r in per_spec.items():
        for t in TESTS:
            x = r[t]
            rows.append({"specimen": s, "patient": _patient(s), "genotype": genotype.get(s, ""), "test": t, "n_eligible": r["n_eligible"],
                         "status": x["status"], "reason_code": x.get("reason_code", ""), "rho": x.get("rho", ""),
                         "p_two_sided": x.get("p_two_sided", ""), "q_bh": x.get("q", ""), "n_surrogates": x.get("n_surrogates", "")})
    return rows


# --------------------------------------------------------------------------- real-data driver
def load_specimen(comp_dir: Path, spec: str):
    f = Path(comp_dir) / f"{spec}.h5ad"
    if not f.is_file():
        return None
    d = R.load_cgeist_eligible(Path(comp_dir), spec)
    try:
        cov, err = H.covariate_matrix(d["obs"], rna=RNA_COVARIATE), ""
    except H.CovariateDomainError as e:
        cov, err = None, str(e)
    return d, cov, err


def run(comp_dir: Path, out_root: Path, protocol_path: Path, B: int = B, seed: int = SEED, analysis_id: str = ANALYSIS_ID) -> dict:
    proto = PIO.load_protocol(protocol_path)
    sha12 = proto["_protocol_sha12"]
    PIO.check_output_target(out_root, sha12, analysis_id)
    if Path(out_root).exists() and any(Path(out_root).iterdir()):
        raise PIO.OutputTargetError(f"{out_root} already exists and is not empty; a run dir is never overwritten")
    specs = proto["specimens"]["ordered"]
    genotype = proto["specimens"]["genotype"]
    per_spec = {}
    for s in specs:
        got = load_specimen(comp_dir, s)
        if got is None:
            g = _nonest("missing_source", f"no {s}.h5ad in {comp_dir}")
            per_spec[s] = {"specimen": s, "n_eligible": 0, **{t: dict(g) for t in TESTS}}
            continue
        d, cov, err = got
        per_spec[s] = analyze_specimen(s, d["obs"]["k_share"].to_numpy(float), d["X"], d["genes"], d["xy"], cov, err, B=B, seed=seed)
        print(f"{s}: n={d['n']} primary={per_spec[s]['primary']['status']} rho={per_spec[s]['primary'].get('rho')}", flush=True)
    result = assemble(per_spec, genotype)
    story = {"stage": "ks1", "protocol_id": proto["protocol_id"], "analysis_id": analysis_id,
             "spec": "docs/superpowers/specs/2026-10-28-midkine-keratin-secretory-association.md",
             "config": {"k": K, "seed": seed, "n_surrogates": B, "secretory_genes": SECRETORY, "min_eligible": MIN_EL, "min_mdk_nonzero": MIN_NZ,
                        "alpha": ALPHA, "compartments_run": str(comp_dir), "statistic": "spearman",
                        "covariates_partial": list(H.primary_covariates(RNA_COVARIATE)), "rna_covariate": RNA_COVARIATE, "protocol_deviation": H.deviation_record(),
                        "headline_specimen": HEADLINE},
             **result}
    sys.path.insert(0, str(REPO / "scripts"))
    import provenance

    stamped = provenance.stamp_payload(PIO.with_protocol_hash(story, proto), PRODUCED_BY)
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    PIO.atomic_write_json(out_root / "keratin_secretory_summary.json", stamped)
    rows = to_rows(per_spec, genotype)
    with open(out_root / "keratin_secretory_per_specimen.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    sha = provenance.code_commit(REPO)
    rel = os.path.relpath(out_root.resolve(), REPO.resolve()) if str(out_root.resolve()).startswith(str(REPO.resolve())) else str(out_root)
    marker = {"path": rel, "jobs": [os.environ.get("SLURM_JOB_ID") or "no-slurm"], "sha": sha,
              "completed": datetime.datetime.now().astimezone().isoformat(timespec="seconds"), "produced_by": PRODUCED_BY,
              "protocol_sha12": sha12}
    PIO.atomic_write_json(out_root / "_COMPLETE.json", marker)
    return stamped


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compartments-run", type=Path, default=None, help="the B5 run dir; REQUIRED for a v42 protocol, legacy default is the v41 c1 run")
    ap.add_argument("--out-root", type=Path, default=None, help="results/<sha12>/<run id>; REQUIRED for a v42 protocol, legacy default is the v41 ks1 dir")
    ap.add_argument("--run-id", default=None, help="the B12 id; REQUIRED for a v42 protocol, legacy default ks1")
    ap.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry("B12", a.protocol, run_id=a.run_id, out_root=a.out_root, comp_run=a.compartments_run)
    comp = PF.legacy_or_required(a.compartments_run, a.protocol, "--compartments-run", REPO / DEFAULT_COMP_RUN)
    out = PF.legacy_or_required(a.out_root, a.protocol, "--out-root", HERE / "results" / "ec0cf5b80715" / ANALYSIS_ID)
    rid = PF.legacy_or_required(a.run_id, a.protocol, "--run-id", ANALYSIS_ID)
    stamp = PF.stage_stamp_begin("B12", a.protocol, Path(out) / "keratin_secretory_summary.json", Path(__file__).resolve(), rid)  # v42: job-start stamp gate, before any compute
    run(comp, out, a.protocol, analysis_id=rid)
    if stamp:
        stamp.finish(comp_run=comp)
    print("done", rid, out, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
