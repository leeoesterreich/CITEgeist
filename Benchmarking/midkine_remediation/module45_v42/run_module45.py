"""Module 4 (+4b) per specimen and Module 5 across specimens, on the v42 SACE layers. Package defaults only; nothing tuned.

    python -u run_module45.py --out-root <NEW dir>                 (unstamped; exploratory)
    python -u run_module45.py --stamp                              (stage B25: stamped; serial; see below)
    python -u run_module45.py --stamp --specimen <name|index>      (stage B25 array task: Module 4 + 4b for one specimen, stamped)
    python -u run_module45.py --stamp --finalize                   (stage B25 dependent job: Module 5 + readouts + final stamps + _COMPLETE.json)

Inputs (fixed, v42 / release v3):
  layers      output/sace_layers_module3_7type/fc43dbea67bf/v2b41/<specimen>/   (B4: SACE layers on RNA-refined proportions)
  proportions output/module3_rna_refined_7type_20261102T200000Z/v2b1/<specimen>_cell_prop_finetuned_results.csv   (B1)
The proportion file md5 must equal the layers manifest's inputs_md5.qp_csv (asserted).
Module 4: discover_programs_from_layers defaults (K=5, lambda_sparsity 0.01, top_n_genes 50, random_state 42), np.random.seed(42) per specimen.
Module 4b: analyze_program_relationships defaults. Module 5: integrate_samples defaults.
Module 5 groups programs within each cell type (package default stratify_by_cell_type=True).
Stage B25 (--stamp, via .commands/submit-m45-stage.sh from a clean detached worktree): out dir = results/<B25 run id of v42_run_ids.yaml>
(must not exist); refuses a dirty tree or a non-SLURM run before computing; writes a .stamp.json sidecar (code_commit, run_id, input sha256s of the
SACE layers and proportions read, slurm job id) for every per_specimen file, module5/module5.pkl and mdk_readout.json; runs mdk_readout.py and
analyze_mdk_secretory.py at the end of the same job. Data are read from $V2_NAMESPACE_ROOT (main checkout), code from the worktree.
Saved per specimen: pickle of the results objects, W/H arrays, programs table, 4b pairs table.
"""
import argparse, hashlib, json, os, pickle, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from CITEgeist.model.programs.anchored_program_discovery import discover_programs_from_layers, analyze_program_relationships
from CITEgeist.model.programs.cross_sample_integration import integrate_samples, save_integration_results

SPECIMENS = ["HCC22-088-P1-S1", "HCC22-088-P1-S2", "HCC22-088-P2-S1", "HCC22-088-P2-S2",
             "HCC22-088-P3-S1_A", "HCC22-088-P3-S2", "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
             "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep", "HCC22-088-P6-S1", "HCC22-088-P6-S2_D"]
# data root: the main checkout (a stamped run executes from a clean worktree that has no output/ tree)
NS = Path(os.environ.get("V2_NAMESPACE_ROOT", REPO))
LAYERS = NS / "output/sace_layers_module3_7type/fc43dbea67bf/v2b41"
PROPS = NS / "output/module3_rna_refined_7type_20261102T200000Z/v2b1"
SEED = 42
STAGE = "B25"
RUN_IDS = REPO / "Benchmarking/midkine_remediation/v42_run_ids.yaml"
REL_RESULTS = "Benchmarking/midkine_remediation/module45_v42/results"
HERE = Path(__file__).resolve().parent


def stage_run_id():
    import yaml

    runs = (yaml.safe_load(RUN_IDS.read_text()) or {}).get("runs") or {}
    assert len(runs) == 1, sorted(runs)
    rid = next(iter(runs.values())).get(STAGE)
    assert rid and rid != "NOT_RUN", f"stage {STAGE} has no run id in {RUN_IDS}"
    return str(rid), str(next(iter(runs)))


def specimen_inputs(spec):
    d = LAYERS / spec
    man = json.loads((d / "manifest.json").read_text())
    names = ["manifest.json", "genes.txt", "barcodes.txt", "coords.csv"] + [f"layer_{t}.npz" for t in man["types"]]
    return [d / n for n in names] + [PROPS / f"{spec}_cell_prop_finetuned_results.csv"]


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def build_adata(spec):
    d = LAYERS / spec
    man = json.loads((d / "manifest.json").read_text())
    assert man["conservation"]["passed"], spec
    qp = PROPS / f"{spec}_cell_prop_finetuned_results.csv"
    assert md5(qp) == man["inputs_md5"]["qp_csv"], f"{spec}: proportions md5 differs from the layers' input"
    genes = [l for l in (d / "genes.txt").read_text().splitlines() if l.strip()]
    bcs = [l for l in (d / "barcodes.txt").read_text().splitlines() if l.strip()]
    xy = pd.read_csv(d / "coords.csv").set_index("barcode").reindex(bcs)[["x", "y"]].to_numpy(float)
    types = man["types"]
    prop = pd.read_csv(qp, index_col=0)
    assert list(prop.columns) == types
    prop = prop.loc[bcs]
    a = ad.AnnData(X=sp.csr_matrix((len(bcs), len(genes)), dtype=np.float64),
                   obs=pd.DataFrame(index=pd.Index(bcs, name="barcodes")), var=pd.DataFrame(index=pd.Index(genes, name="genes")))
    a.obsm["spatial"] = xy
    for t in types:
        L = sp.load_npz(d / f"layer_{t}.npz").tocsr()
        assert L.shape == (len(bcs), len(genes)), (spec, t, L.shape)
        a.layers[f"{t}_genes_pass1"] = L
    return a, prop, types


def run_one(spec, out):
    pk = out / "per_specimen" / f"{spec}.pkl"
    if pk.exists():
        return pickle.loads(pk.read_bytes())
    t0 = time.time()
    a, prop, types = build_adata(spec)
    np.random.seed(SEED)
    m4 = discover_programs_from_layers(a, cell_type_proportions=prop)  # package defaults
    m4b = analyze_program_relationships(m4, a)                          # package defaults
    rec = {"m4": m4, "m4b": m4b, "n_spots": a.n_obs, "n_genes": a.n_vars, "types": types, "seconds": time.time() - t0}
    (out / "per_specimen").mkdir(parents=True, exist_ok=True)
    tmp = pk.with_suffix(".tmp")
    tmp.write_bytes(pickle.dumps(rec))
    tmp.rename(pk)
    # flat tables
    rows = []
    for t, r in m4.results_by_anchor.items():
        for k, p in enumerate(r.programs):
            rows.append({"specimen": spec, "cell_type": t, "program_id": p.program_id, "moran_i": p.spatial_moran_i,
                         "moran_p": p.spatial_moran_pvalue, "variance_explained": p.variance_explained,
                         "mean_activity": p.mean_activity, "n_spots_used": r.n_spots_used, "top_genes": ",".join(p.top_genes)})
    pd.DataFrame(rows).to_csv(out / "per_specimen" / f"{spec}_programs.csv", index=False)
    pd.DataFrame([vars(p) | {"top_genes_overlap": ",".join(p.top_genes_overlap)} for p in m4b.all_pairs]).to_csv(
        out / "per_specimen" / f"{spec}_module4b_pairs.csv", index=False)
    print(f"{spec}: spots={a.n_obs} genes={a.n_vars} programs={m4.total_programs} pairs={m4b.n_pairs_tested} sig={m4b.n_significant} secs={time.time()-t0:.0f}", flush=True)
    return rec


def _stamp_kw(rid, sha12, started):
    import v42_stamp as W  # noqa: F401

    prod = str(Path(__file__).resolve().relative_to(REPO))
    return dict(run_id=rid, protocol_sha12=sha12, producer=prod, started_at=started, repo_root=NS, git_root=REPO, stage_role="module45")


def stamp_specimen(out, s, rid, sha12, started):
    """Sidecars for one specimen's three files; each pins that specimen's layers + proportions.
    The clean-tree check (git status over the whole repo, ~1-2 min on NFS) runs at job start (sbatch + preflight) and again for the LAST stamp;
    the earlier stamps skip it (the producer blob is still checked against HEAD each time)."""
    import v42_stamp as W

    early = dict(_stamp_kw(rid, sha12, started), require_clean_tree=False)
    for n in (f"{s}.pkl", f"{s}_programs.csv", f"{s}_module4b_pairs.csv"):
        W.write_stamp(out / "per_specimen" / n, inputs=specimen_inputs(s), **early)


def stamp_final(out, rid, sha12, started):
    """module5.pkl pins the per-specimen pickles and each specimen's layers manifest + proportions; mdk_readout.json pins module5.pkl,
    the pickles and the cancer layers it reads."""
    import v42_stamp as W

    kw = _stamp_kw(rid, sha12, started)
    early = dict(kw, require_clean_tree=False)
    pk = [out / "per_specimen" / f"{s}.pkl" for s in SPECIMENS]
    light = [p for s in SPECIMENS for p in (LAYERS / s / "manifest.json", PROPS / f"{s}_cell_prop_finetuned_results.csv")]
    W.write_stamp(out / "module5" / "module5.pkl", inputs=pk + light, **early)
    cancer = [p for s in SPECIMENS for p in (LAYERS / s / "genes.txt", LAYERS / s / "layer_Cancer_Cells.npz")]
    W.write_stamp(out / "mdk_readout.json", inputs=[out / "module5" / "module5.pkl"] + pk + light + cancer, **kw)


def stamp_all(out, rid, sha12, started):
    for s in SPECIMENS:
        stamp_specimen(out, s, rid, sha12, started)
    stamp_final(out, rid, sha12, started)


def specimen_outputs(out, s):
    return [out / "per_specimen" / n for n in (f"{s}.pkl", f"{s}_programs.csv", f"{s}_module4b_pairs.csv")]


def run_specimen_task(spec, rid, sha12):
    """Array task: Module 4 + 4b for ONE specimen, stamped. Refuses to overwrite any of that specimen's outputs (or stamps)."""
    import v42_stamp as W

    assert spec in SPECIMENS, f"unknown specimen {spec!r}"
    out = NS / REL_RESULTS / rid
    arts = specimen_outputs(out, spec)
    existing = [p for p in arts if p.exists() or Path(str(p) + ".stamp.json").exists()]
    if existing:
        raise FileExistsError(f"{spec}: outputs exist, a rerun needs a NEW {STAGE} run id: {[str(p) for p in existing]}")
    W.preflight(str(Path(__file__).resolve().relative_to(REPO)), out_paths=arts, repo_root=NS, git_root=REPO)
    started = W.utc_now()
    run_one(spec, out)
    stamp_specimen(out, spec, rid, sha12, started)
    print(f"{STAGE} specimen {spec} stamped", flush=True)


def finalize(rid, sha12):
    """Dependent job: assert all 12 per-specimen outputs exist with valid stamps, then Module 5, readouts, final stamps, _COMPLETE.json."""
    import v42_stamp as W

    out = NS / REL_RESULTS / rid
    final = [out / "mdk_readout.json", out / "module5" / "module5.pkl", out / "_COMPLETE.json"]
    existing = [p for p in final if p.exists()]
    if existing:
        raise FileExistsError(f"{[str(p) for p in existing]} exist; a rerun needs a NEW {STAGE} run id")
    problems = []
    for s in SPECIMENS:
        for p in specimen_outputs(out, s):
            if not p.is_file():
                problems.append(f"missing {p}")
                continue
            try:
                W.verify_stamp(p, run_id=rid, protocol_sha12=sha12, repo_root=NS)
            except Exception as e:
                problems.append(f"{p.name}: {e}")
    if problems:
        raise FileNotFoundError(f"finalize refused ({len(problems)} problems): " + " | ".join(problems))
    W.preflight(str(Path(__file__).resolve().relative_to(REPO)), out_paths=final, repo_root=NS, git_root=REPO)
    started = W.utc_now()
    recs = {s: run_one(s, out) for s in SPECIMENS}  # loads the pickles (all exist)
    finish(out, recs)
    readouts(out)
    stamp_final(out, rid, sha12, started)
    write_complete(out, rid)


def finish(out, recs):
    m5 = integrate_samples({s: r["m4"] for s, r in recs.items()}, {s: r["m4b"] for s, r in recs.items()})  # package defaults
    save_integration_results(m5, out / "module5", prefix="module5_v42")
    (out / "module5" / "module5_summary.txt").write_text(m5.summary() + "\n")
    (out / "module5" / "module5.pkl").write_bytes(pickle.dumps(m5))
    meta = {"layers": str(LAYERS), "proportions": str(PROPS), "specimens": SPECIMENS, "harmony_converged": bool(m5.harmony_converged),
            "n_iterations": int(m5.n_iterations), "n_aligned": len(m5.aligned_programs), "params": m5.parameters,
            "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / "run_meta.json").write_text(json.dumps(meta, indent=1, default=str))
    print(m5.summary(), flush=True)


def write_complete(out, rid):
    (out / "_COMPLETE.json").write_text(json.dumps({"status": "COMPLETE", "stage": STAGE, "run_id": rid, "slurm_job_id": os.environ.get("SLURM_JOB_ID")}, indent=2))
    print(f"{STAGE} stamped {out}", flush=True)


def readouts(out):
    for script in ("mdk_readout.py", "analyze_mdk_secretory.py"):
        subprocess.run([sys.executable, "-u", str(HERE / script), str(out)], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", type=Path, default=None)
    ap.add_argument("--stamp", action="store_true", help="stage B25: out dir from v42_run_ids.yaml, clean tree, stamps, readouts")
    ap.add_argument("--specimen", default=None, help="with --stamp: array task, Module 4 + 4b for this one specimen only (name, or 0-based index)")
    ap.add_argument("--finalize", action="store_true", help="with --stamp: after all 12 specimen tasks, Module 5 + readouts + final stamps + _COMPLETE.json")
    a = ap.parse_args()
    if (a.specimen is not None or a.finalize) and not a.stamp:
        ap.error("--specimen/--finalize require --stamp")
    if a.specimen is not None and a.finalize:
        ap.error("--specimen and --finalize are separate jobs")
    if a.stamp:
        sys.path.insert(0, str(REPO / "scripts"))
        import v42_stamp as W

        rid, sha12 = stage_run_id()
        out = NS / REL_RESULTS / rid
        assert a.out_root in (None, out), f"--stamp writes to {out}, not {a.out_root}"
        if a.specimen is not None:
            spec = SPECIMENS[int(a.specimen)] if a.specimen.isdigit() else a.specimen
            return run_specimen_task(spec, rid, sha12)
        if a.finalize:
            return finalize(rid, sha12)
        if out.exists():
            raise FileExistsError(f"{out} exists; a rerun needs a NEW {STAGE} run id")
        W.preflight(str(Path(__file__).resolve().relative_to(REPO)), out_paths=[out / "mdk_readout.json"], repo_root=NS, git_root=REPO)
        started = W.utc_now()
        out.mkdir(parents=True)
    else:
        assert a.out_root, "--out-root is required without --stamp"
        out = a.out_root
        out.mkdir(parents=True, exist_ok=True)
    recs = {s: run_one(s, out) for s in SPECIMENS}
    finish(out, recs)
    if a.stamp:
        readouts(out)
        stamp_all(out, rid, sha12, started)
        write_complete(out, rid)


if __name__ == "__main__":
    main()
