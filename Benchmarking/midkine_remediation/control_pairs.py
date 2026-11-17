"""B17: matched random gene controls for the CITEgeist MDK-secretory bivariate Moran's I (Definition C), 200 per specimen. Reported in Section 2.4 (Fig 5C).

Promoted from the exploratory check 2 (``exploratory_v42/check2_control_pairs.py``, git history). The analysis code below is that code, unchanged:
``deciles``, ``make_draws``, ``pctile_block``, ``method_block``, ``prepare_pool`` and ``citegeist_alone`` are byte-for-byte the exploratory
functions. What was removed: the Cell2location and Tangram blocks (the cross-method comparison is no longer reported, so no competitor file is an input).
The seeds, draw counts and matching rule are unchanged and are written into the output (``config``).

Question: is the observed MDK-secretory bivariate Moran's I of one cancer compartment larger than the same statistic for gene pairs that are
unrelated to the hypothesis, under the same graph and the same pipeline? Descriptive: no p-values, nothing here is called a cause.

Matching rule (also in the output)
  * Genes: all genes of the CITEgeist layer (``N_P``, eligible Cancer_Cells compartments of the B5 run), own kNN(k=6) graph (Definition C).
  * Pool = genes minus MDK, the 11 secretory genes, the 6 basal keratins, and any constant gene.
  * Expression decile = decile 0..9 of the gene's mean over the eligible spots (average ranks over all genes; ties share a decile).
  * Control (a), 200 draws: MDK replaced by a random pool gene from MDK's decile; y = the real secretory score.
  * Control (b), 200 draws: MDK kept; the 11 secretory genes replaced by 11 distinct random pool genes, each from its secretory gene's decile;
    y = mean of per-gene z-scores (as the real score).
  * Univariate reference: median of morans_bv(g, g, W) over a fixed random 500 pool genes.
  * Seed: ``np.random.SeedSequence([20261103, crc32(specimen + "|cg_alone")])`` spawned into three streams (a, b, univariate).

Specimens: those with a Definition C value in the B9 run (all but P4-S1, which is not estimable); P4-S2 first.
GATE (fail-closed): the recomputed real I must equal the B9 stored ``per_specimen`` Definition C ``value``, 1e-9 relative. Any specimen failing stops the
job: ``gate_failure.json`` is written, no result file is written, no stamp exists, exit code 2.

Inputs come from the CLI (set by ``.commands/submit-v2-controls.sh`` through ``chain_preflight.py``): the B5 run (``--compartments-run``), the B9 run
(``--morans-run``) and ``--out-root``; nothing is read from the checkout the code runs in.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import zlib
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import a12_10_real as R  # noqa: E402
import chain_guard as PF  # noqa: E402
import compartments as C  # noqa: E402
import morans as M  # noqa: E402
import promoted_common as CM  # noqa: E402
import protocol_io as PIO  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/control_pairs.py"
STAGE = "B17"
ARTIFACT = "control_pairs_results.json"
HEADLINE = "HCC22-088-P4-S2_1i_rep"
SEED = 20261103
N_CTRL = 200
N_UNI = 500
EXCLUDE = set(R.SECRETORY) | {C.MDK} | set(C.KERATINS)
MATCHING_RULE = (
    "Genes: all genes of the CITEgeist layer (N_P, eligible Cancer_Cells compartments of the B5 run) on its own kNN(k=6) graph (Definition C). "
    "Pool = genes minus MDK, the 11 secretory genes, the 6 basal keratins and any constant gene. Expression decile = decile 0..9 of the gene's mean "
    "over the eligible spots (average ranks over all genes; ties share a decile). Control (a), 200 draws: MDK replaced by a random pool gene from "
    "MDK's decile, y = the real secretory score. Control (b), 200 draws: MDK kept, the 11 secretory genes replaced by 11 distinct random pool genes, "
    "each from its secretory gene's decile, y = mean of per-gene z-scores. Univariate reference: median morans_bv(g, g) over a fixed random 500 pool genes. "
    "Seed: SeedSequence([20261103, crc32(specimen + '|cg_alone')]) spawned into three streams (a, b, univariate). Observed I only; no p-values."
)


# --------------------------------------------------------------------------- draws (unchanged from the exploratory check 2)
def deciles(means) -> np.ndarray:
    """Decile 0..9 of each gene by the rank of its mean (average ranks, ties share a decile)."""
    r = rankdata(means, method="average")
    return np.minimum(((r - 1) / len(means) * 10).astype(int), 9)


def make_draws(genes, dec, pool_idx, mdk_i, sec_idx, specimen: str, n_ctrl: int = N_CTRL, n_uni: int = N_UNI, seed: int = SEED) -> dict:
    """Fixed-seed gene-index draws (indices into ``genes``); identical for every method of a specimen."""
    ss = np.random.SeedSequence([seed, zlib.crc32(specimen.encode())])
    ra, rb, ru = (np.random.default_rng(s) for s in ss.spawn(3))
    pool_idx = np.asarray(pool_idx)
    by_dec = {k: pool_idx[dec[pool_idx] == k] for k in range(10)}
    cand_a = by_dec[dec[mdk_i]]
    if len(cand_a) == 0:
        raise ValueError(f"{specimen}: no pool gene in MDK's decile")
    a = ra.choice(cand_a, size=n_ctrl, replace=True)
    b = np.empty((n_ctrl, len(sec_idx)), dtype=int)
    for i in range(n_ctrl):
        used = set()
        for j, s in enumerate(sec_idx):
            c = [g for g in by_dec[dec[s]] if g not in used]
            if not c:
                raise ValueError(f"{specimen}: decile pool exhausted for secretory gene {genes[s]}")
            g = int(rb.choice(c))
            used.add(g)
            b[i, j] = g
    uni = ru.choice(pool_idx, size=min(n_uni, len(pool_idx)), replace=False)
    return {"a": a, "b": b, "uni": uni}


def pctile_block(v) -> dict:
    v = np.asarray(v, float)
    return {"mean": float(v.mean()), "median": float(np.median(v)), "p5": float(np.percentile(v, 5)), "p95": float(np.percentile(v, 95)),
            "n": int(len(v))}


def method_block(T, genes, W, mdk_i, sec_idx, draws, real_I=None) -> dict:
    """Real I and control distributions for one method matrix T (n x len(genes)) on graph W."""
    x_mdk = T[:, mdk_i]
    y_real = R._mean_z_vecs([T[:, j] for j in sec_idx])
    real = M.morans_bv(x_mdk, y_real, W)
    Ia = np.array([M.morans_bv(T[:, g], y_real, W) for g in draws["a"]])
    Ib = np.array([M.morans_bv(x_mdk, R._mean_z_vecs([T[:, g] for g in row]), W) for row in draws["b"]])
    uni = np.array([M.morans_bv(T[:, g], T[:, g], W) for g in draws["uni"]])
    out = {"real_I": float(real), "univariate_median_I": float(np.median(uni)), "n_univariate_genes": int(len(uni))}
    for tag, I in (("control_a_random_mdk", Ia), ("control_b_random_secretory", Ib)):
        blk = pctile_block(I)
        blk["real_percentile_in_controls"] = float(100.0 * np.mean(I < real))
        blk["real_gt_p95"] = bool(real > blk["p95"])
        blk["real_within_p5_p95"] = bool(blk["p5"] <= real <= blk["p95"])
        out[tag] = blk
    return out


def prepare_pool(T_list, genes, cg_T, mdk_i, sec_idx) -> tuple:
    """Pool indices: not excluded and non-constant in every matrix; deciles from the CITEgeist matrix."""
    dec = deciles(cg_T.mean(axis=0))
    ok = np.ones(len(genes), bool)
    for T in T_list:
        ok &= T.std(axis=0) > 0
    pool = np.array([i for i, g in enumerate(genes) if g not in EXCLUDE and ok[i]])
    return dec, pool


def citegeist_alone(spec: str, d: dict, stored: dict, ref_id: str) -> dict:
    """``ref_id`` is the id of the B9 run the gate reads (its directory name); it names the gate, so a label can never say a run that was not used."""
    genes, X, xy = list(d["genes"]), np.asarray(d["X"], float), d["xy"]
    n = d["n"]
    mdk = R._column(X, genes, C.MDK)
    if n < R.MIN_EL or n > M.DENSE_EIGH_MAX_N or int((mdk > 0).sum()) < R.MIN_NZ:
        return {"status": "not_estimable", "n_eligible": int(n)}
    W = knn_weights(np.asarray(xy, float), R.K)
    mdk_i = genes.index(C.MDK)
    sec_idx = [genes.index(g) for g in R.SECRETORY]
    real = M.morans_bv(X[:, mdk_i], R._mean_z_vecs([X[:, j] for j in sec_idx]), W)
    gate = CM.gate_equal(f"{spec} CITEgeist-alone Definition C I vs {ref_id} per_specimen", real, stored["per_specimen"][spec]["value"])
    dec, pool = prepare_pool([X], genes, X, mdk_i, sec_idx)
    dr = make_draws(genes, dec, pool, mdk_i, sec_idx, spec + "|cg_alone")
    return {"status": "ok", "n_eligible": int(n), "pool_size": int(len(pool)), "gate": gate, **method_block(X, genes, W, mdk_i, sec_idx, dr)}


def cohort_alone(rows: dict) -> dict:
    """Cohort summary of the CITEgeist-alone Definition C block (the exploratory ``cohort`` code for that block, unchanged)."""
    out = {}
    ok = {s: r for s, r in rows.items() if r["status"] == "ok"}
    if ok:
        cg = {"n_specimens": len(ok), "mean_real_I": float(np.mean([v["real_I"] for v in ok.values()])),
              "mean_univariate_median_I": float(np.mean([v["univariate_median_I"] for v in ok.values()]))}
        for c in ("control_a_random_mdk", "control_b_random_secretory"):
            cg[c] = {"mean_of_control_means": float(np.mean([v[c]["mean"] for v in ok.values()])), "mean_p95": float(np.mean([v[c]["p95"] for v in ok.values()])),
                     "n_real_gt_p95": int(sum(v[c]["real_gt_p95"] for v in ok.values())), "n_real_within_p5_p95": int(sum(v[c]["real_within_p5_p95"] for v in ok.values()))}
        out["citegeist_alone_definition_C"] = cg
    return out


def specimens_of(stored: dict) -> list:
    """Specimens with a Definition C value in the B9 run, P4-S2 first."""
    est = [s for s, r in stored["per_specimen"].items() if r.get("status") != "not_estimable"]
    return sorted(est, key=lambda s: (s != HEADLINE, s))


def run(comp_run, morans_run, out_root, run_id: str = "", specimens=None, meta_extra=None) -> int:
    """Compute the controls for every specimen, gate each against B9, write ``control_pairs_results.json`` + ``_COMPLETE.json``. 0 ok, 2 gate failure."""
    stored = json.loads((Path(morans_run) / "morans_results.json").read_text())
    specs = list(specimens) if specimens else specimens_of(stored)
    out = Path(out_root)
    out.mkdir(parents=True, exist_ok=True)
    code = HERE.parents[1]
    git = lambda *a: subprocess.run(["git", "-C", str(code), *a], capture_output=True, text=True).stdout.strip()  # noqa: E731
    meta = {"stage": STAGE, "run_id": run_id, "produced_by": PRODUCED_BY, "descriptive_only": True, "seed": SEED, "seed_rule": "SeedSequence([seed, crc32(specimen + '|cg_alone')]).spawn(3)",
            "n_controls_per_type": N_CTRL, "n_univariate_genes": N_UNI, "k": R.K, "matching_rule": MATCHING_RULE, "compartments_run": str(comp_run),
            "morans_reference": str(Path(morans_run) / "morans_results.json"), "excluded_from_pool": sorted(EXCLUDE), "git_commit": git("rev-parse", "HEAD"),
            **(meta_extra or {})}
    rows = {}
    for s in specs:
        t = time.time()
        try:
            d = R.load_cgeist_eligible(Path(comp_run), s)
            rows[s] = {"specimen": s, **citegeist_alone(s, d, stored, Path(morans_run).name)}
        except CM.GateFailure as e:
            PIO.atomic_write_json(out / "gate_failure.json", {**meta, "gate_failure": str(e), "specimens_done": list(rows)})
            print(str(e), file=sys.stderr, flush=True)
            return 2
        print(f"{s}: {rows[s]['status']} {time.time()-t:.0f}s gate passed", flush=True)
    PIO.atomic_write_json(out / ARTIFACT, {**meta, "per_specimen": rows, "cohort": cohort_alone(rows)})
    PIO.atomic_write_json(out / "_COMPLETE.json", {"path": str(out), "jobs": [os.environ.get("SLURM_JOB_ID") or "no-slurm"], "sha": meta["git_commit"],
                                                  "completed": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "produced_by": PRODUCED_BY})
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compartments-run", type=Path, required=True, help="the B5 run dir (id from v42_run_ids.yaml)")
    ap.add_argument("--morans-run", type=Path, required=True, help="the B9 run dir (id from v42_run_ids.yaml); its Definition C values are the gate")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<B17 id>; never overwritten")
    ap.add_argument("--run-id", required=True, help="the B17 id")
    ap.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry(STAGE, a.protocol, run_id=a.run_id, out_root=a.out_root, comp_run=a.compartments_run, morans_run=a.morans_run)
    if a.out_root.exists():
        raise SystemExit(f"refusing: {a.out_root} exists (a run dir is never overwritten)")
    stamp = PF.stage_stamp_begin(STAGE, a.protocol, a.out_root / ARTIFACT, Path(__file__).resolve(), a.run_id)  # job-start stamp gate, before any compute
    rc = run(a.compartments_run, a.morans_run, a.out_root, a.run_id)
    if rc == 0 and stamp:
        stamp.finish(comp_run=a.compartments_run, morans_run=a.morans_run)
    print("done" if rc == 0 else "gate failed", a.run_id, a.out_root, flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
