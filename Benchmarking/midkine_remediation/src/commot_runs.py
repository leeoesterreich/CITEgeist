"""A12.8: COMMOT run matrix for one specimen (protocol ``commot.run_matrix``). Run with the COMMOT env python (commot 0.0.3, numpy 1.26).

Runs (all on the compartments of the A12.6 layers; one process, ONE ``import commot``):
    R00 N_P, p_construct 0.05, n_floor 0.5, dis_thr 500        R08 n_floor 0.25     R09 n_floor 1.0
    R01 N_A                                                     R10 true between-spot (same-spot pairs and the diagonal pushed beyond dis_thr)
    R11 N_P_nf (nuclei-free; same compartments/labels as R00)
    P4-S2 only: R02 N_S, R03/R04/R05 p_construct 0.02/0.10/0.20, R06/R07 dis_thr 1.5/3.0 median spot pitches.
LR list: ``filter_lr_database`` (CellChat, human, heteromeric, min_cell_pct 0.01) on the R00 N_P object, frozen and reused by every run.
Every run uses the FULL frozen LR list. (An MDK-only run was tried and is NOT equivalent: COMMOT divides all masses by max(S.sum(), D.sum()) over the
whole run, so absolute sender scores depend on the other pairs; ``--verify-subset`` records the full-vs-MDK-only difference as a diagnostic.)
Per run: sender/receiver tables (all compartments), the five obsp MDK-pair matrices (sparse), the G5 key-mapping check, and the exact
within-spot / between-spot decomposition (auto + within + between == total, asserted).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(HERE.parent))
import chain_guard as PF  # noqa: E402
import provenance  # noqa: E402
import commot_gates as G  # noqa: E402
import compartments as C  # noqa: E402
import protocol_io as PIO  # noqa: E402

MDK_PAIRS = [("MDK", "NCL"), ("MDK", "SDC2"), ("MDK", "SDC4"), ("MDK", "LRP1"), ("MDK", "ITGA6_ITGB1")]
DIS_THR = 500.0
HEADLINE = "HCC22-088-P4-S2_1i_rep"
R10_FAR = 1e6
# Per-run wall-clock budget for the forked child (3 h = ~1.77x the largest measured healthy full-LR run; a non-convergent run dies on its own
# as a nonzero exit well before this, so no legitimate run is ever timeout-killed). Overridable via the run_child_timeout_s argument or the
# RUN_CHILD_TIMEOUT_S env var (the hang regression test sets it small).
RUN_CHILD_TIMEOUT_S = 10800


def run_list(specimen: str) -> list:
    base = dict(input="N_P", p_construct=0.05, n_floor=0.5, dis_thr=DIS_THR, true_between=False)
    runs = [{"run_id": "R00", **base}, {"run_id": "R01", **base, "input": "N_A"}, {"run_id": "R08", **base, "n_floor": 0.25},
            {"run_id": "R09", **base, "n_floor": 1.0}, {"run_id": "R10", **base, "true_between": True}, {"run_id": "R11", **base, "input": "N_P_nf"}]
    if specimen == HEADLINE:
        runs += [{"run_id": "R02", **base, "input": "N_S"}, {"run_id": "R03", **base, "p_construct": 0.02}, {"run_id": "R04", **base, "p_construct": 0.10},
                 {"run_id": "R05", **base, "p_construct": 0.20}, {"run_id": "R06", **base, "dis_thr": "1.5pitch"}, {"run_id": "R07", **base, "dis_thr": "3.0pitch"}]
    return runs


# --------------------------------------------------------------------------- pure helpers (unit-tested without commot)
def spot_pitch(coords_by_row, spot_id) -> float:
    """Median nearest-neighbour distance between distinct spot centres (coordinate units)."""
    from scipy.spatial import cKDTree

    u = pd.DataFrame(np.asarray(coords_by_row, float)).assign(s=np.asarray(spot_id)).drop_duplicates("s").iloc[:, :2].to_numpy()
    d, _ = cKDTree(u).query(u, k=2)
    return float(np.median(d[:, 1]))


def derive_np(E_lr, n_st, depth_ratio, floor):
    """N_P at another floor: log1p(E / max(n_st, floor) / depth_ratio) (same formula as ``compartments.build_compartments``)."""
    return np.log1p(E_lr / np.maximum(n_st, floor)[:, None] / np.asarray(depth_ratio)[:, None])


def decomposition(P: sp.spmatrix, spot_id) -> dict:
    """Exact split of the transported mass (rows = sender compartments, cols = receivers): autocrine (diagonal), within-spot (other compartment
    in the same spot) and between-spot. ``auto + within + between == total`` per row and overall (asserted)."""
    P = sp.coo_matrix(P)
    spot = np.asarray(spot_id)
    same_spot = spot[P.row] == spot[P.col]
    diag = P.row == P.col
    n = P.shape[0]
    auto = np.bincount(P.row[diag], weights=P.data[diag], minlength=n)
    within = np.bincount(P.row[same_spot & ~diag], weights=P.data[same_spot & ~diag], minlength=n)
    between = np.bincount(P.row[~same_spot], weights=P.data[~same_spot], minlength=n)
    total_row = np.asarray(P.tocsr().sum(axis=1)).ravel()
    err = float(np.max(np.abs(auto + within + between - total_row))) if n else 0.0
    tot = float(total_row.sum())
    if err > 1e-9 * max(1.0, float(np.abs(total_row).max() if n else 1.0)):
        raise AssertionError(f"decomposition does not sum to the row totals (max err {err})")
    out = {"total": tot, "autocrine": float(auto.sum()), "within_spot": float(within.sum()), "between_spot": float(between.sum()),
           "max_row_error": err}
    for k in ("autocrine", "within_spot", "between_spot"):
        out[f"share_{k}"] = out[k] / tot if tot > 0 else None
    out["_rows"] = {"auto": auto, "within": within, "between": between}
    return out


def select_rows(obs: pd.DataFrame, p_construct: float) -> np.ndarray:
    """Compartments kept at a higher construction cutoff are an exact row subset of the base set (the stored proportions and row sums do not change)."""
    return (obs["qp_prop"].to_numpy() >= p_construct)


# --------------------------------------------------------------------------- COMMOT (needs the COMMOT env)
def lr_gene_set(db) -> set:
    g = set()
    for col in (0, 1):
        for x in db.iloc[:, col]:
            g.update(str(x).split("_"))
    return g


def _adata(X, spatial, genes, names):
    import anndata as ad

    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=list(names)), var=pd.DataFrame(index=list(genes)))
    a.obsm["spatial"] = np.asarray(spatial, float)
    return a


def freeze_lr(Xlr, spatial, genes, names, min_cell_pct=0.01):
    import commot as ct

    db = ct.pp.ligand_receptor_database(species="human", database="CellChat")
    return ct.pp.filter_lr_database(db, _adata(Xlr, spatial, genes, names), heteromeric=True, min_cell_pct=min_cell_pct).reset_index(drop=True)


def mdk_subset(lr: pd.DataFrame) -> pd.DataFrame:
    keep = [(str(r[0]), str(r[1])) in set(MDK_PAIRS) for r in lr.values]
    return lr[np.array(keep, bool)].reset_index(drop=True)


def commot_run(Xlr, spatial, genes, names, lr, dis_thr, dmat=None):
    """One COMMOT call -> dict(sender, receiver, obsp{pair: csr}, has_spatial_distance). Raises on a COMMOT error (callers type it)."""
    import commot as ct

    a = _adata(Xlr, spatial, genes, names)
    if dmat is not None:
        a.obsp["spatial_distance"] = dmat
    ct.tl.spatial_communication(a, database_name="cellchat", df_ligrec=lr, dis_thr=dis_thr, heteromeric=True, pathway_sum=True)
    keys = [k for k in a.obsp.keys() if k.startswith("commot-cellchat-") and k.count("-") >= 3 and not k.endswith("total-total")]
    return {"sender": a.obsm["commot-cellchat-sum-sender"].copy(), "receiver": a.obsm["commot-cellchat-sum-receiver"].copy(),
            "obsp": {k: sp.csr_matrix(a.obsp[k]) for k in keys}, "has_spatial_distance": "spatial_distance" in a.obsp}


def run_one(spec: dict, Xlr, spatial, genes, names, spot_id, lr_all, lr_run, dis_thr, run_dir: Path) -> dict:
    """Execute one run id; always writes run.json (typed failure on a COMMOT error); returns the record."""
    t0 = time.time()
    rec = {k: (v if not isinstance(v, np.generic) else v.item()) for k, v in spec.items()}
    rec.update({"dis_thr_value": float(dis_thr), "n_rows": int(Xlr.shape[0]), "n_lr_pairs_run": int(len(lr_run)), "n_lr_pairs_frozen": int(len(lr_all))})
    run_dir.mkdir(parents=True, exist_ok=True)
    try:
        dmat = None
        if spec["true_between"]:
            from scipy.spatial import distance_matrix

            d = distance_matrix(spatial, spatial)
            same = np.asarray(spot_id)[:, None] == np.asarray(spot_id)[None, :]
            dmat = np.where(same, R10_FAR, d)
            del d, same
            pre = G.r10_precheck(dmat, dis_thr)
            if pre["status"] != "ok":
                rec.update(status="not_estimable", reason=pre["reason"], seconds=round(time.time() - t0, 1))
                (run_dir / "run.json").write_text(json.dumps(rec, indent=2))
                return rec
            rec["r10_pairs_inside_cutoff"] = int(pre["pairs_inside_cutoff"])
        res = commot_run(Xlr, spatial, genes, names, lr_run, dis_thr, dmat)
        del dmat
        PIO.assert_no_spatial_distance(types_simple(res["has_spatial_distance"]), spec["run_id"] if spec["true_between"] else "R00")
        sender = {c: res["sender"][c].to_numpy() for c in res["sender"].columns}
        rs = {k: np.asarray(P.sum(axis=1)).ravel() for k, P in res["obsp"].items()}
        frozen = [tuple(map(str, r[:2])) for r in lr_all.values]
        g5 = G.key_mapping_check(rs, sender, frozen)
        rec["G5"] = {"status": g5["status"], "pairs": g5["pairs"]}
        keep_cols = [c for c in res["sender"].columns if c in G.PAIR_MAP]
        tab = pd.DataFrame({c: res["sender"][c].to_numpy() for c in keep_cols}, index=list(names))
        for c in [c for c in res["receiver"].columns if c.replace("r-", "s-", 1) in G.PAIR_MAP]:
            tab[c] = res["receiver"][c].to_numpy()
        dec = {}
        for scol, okey in G.PAIR_MAP.items():
            if okey in res["obsp"]:
                P = res["obsp"][okey]
                sp.save_npz(run_dir / f"obsp_{scol}.npz", P, compressed=True)
                d = decomposition(P, spot_id)
                if scol == "s-MDK-NCL":
                    for k, v in d["_rows"].items():
                        tab[f"ncl_{k}"] = v
                dec[scol] = {k: v for k, v in d.items() if k != "_rows"}
        tab.to_csv(run_dir / "sender.csv.gz")
        rec["decomposition"] = dec
        rec["status"] = "ok" if g5["pairs"].get("s-MDK-NCL", {}).get("ok") else "failed"
        rec["reason"] = "" if rec["status"] == "ok" else f"G5 key-mapping failed for s-MDK-NCL: {g5['pairs'].get('s-MDK-NCL')}"
    except Exception as e:  # noqa: BLE001
        why = "no pair inside the cutoff (COMMOT zero-size reduction)" if "zero-size" in str(e) else f"{type(e).__name__}: {e}"
        rec.update(status="failed", reason=why)
    rec["seconds"] = round(time.time() - t0, 1)
    (run_dir / "run.json").write_text(json.dumps(rec, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    return rec


class types_simple:
    def __init__(self, has):
        self.obsp = {"spatial_distance": None} if has else {}


# --------------------------------------------------------------------------- specimen driver
def _run_one_child(spec: dict, Xrun, sub_xy, lr_genes, sub_names, sub_spot, lr_all, lr_run, dis_thr, run_dir: Path) -> None:
    """Forked-child entrypoint: execute one run_id via ``run_one``, then detach. The parent recovers the result not from a return value
    (a forked child's return does not reach the parent) but by re-reading ``run_dir/run.json``, which ``run_one`` writes in every case
    (ok or typed failure)."""
    run_one(spec, Xrun, sub_xy, lr_genes, sub_names, sub_spot, lr_all, lr_run, dis_thr, run_dir)


def _synthesize_rec(spec: dict, Xrun, dis_thr: float, status: str, reason: str, seconds: float) -> dict:
    """Minimal record for a run whose child did not produce run.json (died / timed out / aborted), so the caller's ``records.append``,
    ``print(...)`` and ``meta["runs"]`` keep working non-None."""
    rec = {k: (v if not isinstance(v, np.generic) else v.item()) for k, v in spec.items()}
    rec.update({"dis_thr_value": float(dis_thr), "n_rows": int(Xrun.shape[0]), "status": status, "reason": reason, "seconds": seconds})
    return rec


def _write_run_json(run_dir: Path, rec: dict) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "run.json").write_text(json.dumps(rec, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


def _run_one_isolated(spec: dict, Xrun, sub_xy, lr_genes, sub_names, sub_spot, lr_all, lr_run, dis_thr, run_dir: Path,
                      run_child_timeout_s: int | None = None) -> dict:
    """Run one run_id in a forked child with a wall-clock budget, so a COMMOT non-convergence that surfaces only as a silent child death
    (no exception reaches the parent) is typed as a failed run.json instead of killing the batch. The parent owns run.json on every non-ok
    path (the child's return value never reaches us; we re-read run.json instead)."""
    timeout = run_child_timeout_s if run_child_timeout_s is not None else int(os.environ.get("RUN_CHILD_TIMEOUT_S", RUN_CHILD_TIMEOUT_S))
    t_run = time.time()
    ctx = multiprocessing.get_context("fork")
    p = ctx.Process(target=_run_one_child, args=(spec, Xrun, sub_xy, lr_genes, sub_names, sub_spot, lr_all, lr_run, dis_thr, run_dir))
    p.start()
    p.join(timeout)
    f = run_dir / "run.json"
    if p.exitcode == 0:
        if f.exists():
            try:
                rec = json.loads(f.read_text())
                if not isinstance(rec, dict):
                    raise ValueError("run.json is not a JSON object")
            except Exception:  # noqa: BLE001 — a corrupt child run.json must not propagate out of the batch loop
                rec = _synthesize_rec(spec, Xrun, dis_thr, "failed", "child-run.json-corrupt", round(time.time() - t_run, 1))
                _write_run_json(run_dir, rec)
        else:
            rec = _synthesize_rec(spec, Xrun, dis_thr, "failed", "child-aborted-without-run.json", round(time.time() - t_run, 1))
            _write_run_json(run_dir, rec)
    elif p.exitcode is not None:  # child died on its own (nonzero) — the silent non-convergence class
        rec = _synthesize_rec(spec, Xrun, dis_thr, "failed", f"commot-ot-nonconvergence-or-child-exit (exit {p.exitcode})",
                              round(time.time() - t_run, 1))
        _write_run_json(run_dir, rec)
    else:  # still alive after the budget -> timeout; terminate, and kill as the mandatory fallback for a SIGTERM-ignoring child
        p.terminate()
        p.join(10)
        if p.is_alive():
            p.kill()
        p.join()
        rec = _synthesize_rec(spec, Xrun, dis_thr, "failed", "child-timeout", round(time.time() - t_run, 1))
        _write_run_json(run_dir, rec)
    return rec


def run_matrix(specimen, protocol_sha256, layers, genes, qp, nuc, umi, coords, M, tmp: Path, verify_subset=False, dis_thr_primary=DIS_THR,
               run_child_timeout_s: int | None = None, resume: bool = False):
    """Run the full matrix for one specimen into ``tmp`` (inputs already in memory: this is what the toy test drives with ``dis_thr_primary`` 3)."""
    import commot as ct

    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1)  # no NEW multi-threaded BLAS pool is created after the guard (fork-safe)
    t0 = time.time()
    base = C.build_compartments(layers, genes, qp, nuc, umi, coords, M)
    el = C.eligible_cancer(base.obs).to_numpy()
    spot_id = base.obs["spot_id"].to_numpy()
    names = list(base.obs.index)
    pitch = spot_pitch(base.spatial, spot_id)
    db = ct.pp.ligand_receptor_database(species="human", database="CellChat")
    need = lr_gene_set(db)
    genes = list(genes)
    gi = np.array([i for i, g in enumerate(genes) if g in need])
    lr_genes = [genes[i] for i in gi]
    Elr = base.E[:, gi]
    lr_all = freeze_lr(sp.csr_matrix(base.X["N_P"][:, gi]), base.spatial, lr_genes, names)
    lr_mdk = mdk_subset(lr_all)
    tmp = Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    lr_all.to_csv(tmp / "lr_frozen.csv", index=False)
    pd.DataFrame({"id": names, "eligible_cancer": el, "k_share": base.obs["k_share"].to_numpy(), "spot_id": spot_id}).to_csv(tmp / "rows_R00.csv.gz", index=False)
    meta = {"specimen": specimen, "median_spot_pitch": pitch, "dis_thr_primary": dis_thr_primary, "dis_thr_in_pitch_units": dis_thr_primary / pitch, "M": M,
            "n_rows_R00": len(names), "n_eligible_R00": int(el.sum()), "n_lr_frozen": int(len(lr_all)), "n_lr_mdk_in_frozen": int(len(lr_mdk)),
            "mdk_pairs_run": [list(map(str, r[:2])) for r in lr_mdk.values], "protocol_sha256": protocol_sha256,
            "commot_version": getattr(ct, "__version__", "0.0.3"), "numpy": np.__version__}
    if verify_subset:
        Xn = sp.csr_matrix(base.X["N_P"][:, gi])
        full = commot_run(Xn, base.spatial, lr_genes, names, lr_all, dis_thr_primary)
        sub = commot_run(Xn, base.spatial, lr_genes, names, lr_mdk, dis_thr_primary)
        diag = {}
        for c in sub["sender"].columns:
            if c in G.PAIR_MAP and c in full["sender"].columns:
                a_, b_ = full["sender"][c].to_numpy(), sub["sender"][c].to_numpy()
                nz = (a_ > 0) & (b_ > 0)
                ratio = b_[nz] / a_[nz]
                diag[c] = {"max_abs_diff": float(np.max(np.abs(a_ - b_))), "n_both_active": int(nz.sum()),
                           "ratio_sub_over_full_cv": float(ratio.std() / ratio.mean()) if nz.sum() > 1 else None,
                           "ratio_sub_over_full_mean": float(ratio.mean()) if nz.sum() else None}
        meta["subset_diagnostic"] = diag
    records = []
    for spec in run_list(specimen):
        thr = spec["dis_thr"]
        thr = float(thr.replace("pitch", "")) * pitch if isinstance(thr, str) else float(dis_thr_primary)
        rid = spec["run_id"]
        if resume:
            _jf = tmp / rid / "run.json"
            if _jf.exists():
                try:
                    _prev = json.loads(_jf.read_text())
                    # Resume reuse requires the run's key artifacts to be present (rows.csv.gz is written by the parent AFTER the
                    # forked child returns, so run.json==ok alone does not guarantee it; sender.csv.gz is written by the child on
                    # every success). If any is missing, fall through and re-run rather than reusing a partially-materialised run.
                    if _prev.get("status") == "ok" and (tmp / rid / "rows.csv.gz").exists() and (tmp / rid / "sender.csv.gz").exists():
                        records.append(_prev)
                        print(specimen, rid, "resumed_ok", _prev.get("seconds"), flush=True)
                        continue
                except Exception:
                    pass
        if spec["input"] == "N_P" and spec["p_construct"] == 0.05:
            Xrun = derive_np(Elr, base.obs["n_st"].to_numpy(), base.obs["depth_ratio"].to_numpy(), spec["n_floor"]) if spec["n_floor"] != 0.5 else base.X["N_P"][:, gi]
            sub_names, sub_xy, sub_spot = names, base.spatial, spot_id
        elif spec["input"] == "N_P" and spec["p_construct"] > 0.05:
            rows = select_rows(base.obs, spec["p_construct"])
            Xrun = base.X["N_P"][rows][:, gi]
            sub_names, sub_xy, sub_spot = [n for n, k in zip(names, rows) if k], base.spatial[rows], spot_id[rows]
        elif spec["input"] == "N_P":   # p_construct 0.02: rebuild from the layers
            c2 = C.build_compartments(layers, genes, qp, nuc, umi, coords, M, p_construct=spec["p_construct"])
            Xrun, sub_names, sub_xy, sub_spot = c2.X["N_P"][:, gi], list(c2.obs.index), c2.spatial, c2.obs["spot_id"].to_numpy()
            del c2
        else:
            Xrun = base.X[spec["input"]][:, gi]
            sub_names, sub_xy, sub_spot = names, base.spatial, spot_id
        rec = _run_one_isolated(spec | {"dis_thr": thr}, sp.csr_matrix(Xrun), sub_xy, lr_genes, sub_names, sub_spot, lr_all, lr_all, thr, tmp / rid,
                                run_child_timeout_s)
        pd.DataFrame({"id": sub_names}).to_csv(tmp / rid / "rows.csv.gz", index=False)
        records.append(rec)
        print(specimen, rid, rec["status"], rec.get("reason", "")[:120], rec["seconds"], flush=True)
    meta["runs"] = [{k: r.get(k) for k in ("run_id", "status", "reason", "seconds", "n_rows", "dis_thr_value")} for r in records]
    meta["seconds"] = round(time.time() - t0, 1)
    meta["git_head"] = provenance.code_commit(REPO_ROOT)
    (tmp / "specimen.json").write_text(json.dumps(meta, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    return meta


def _stale_tmp_dirs(out_root, specimen) -> list:
    """Staged ``.tmp_<specimen>_<pid>`` dirs only (a specimen whose name prefixes another's, e.g. X vs X_rep, must not match)."""
    pat = re.compile(rf"^\.tmp_{re.escape(specimen)}_\d+$")
    return sorted(p for p in Path(out_root).glob(f".tmp_{specimen}_*") if pat.match(p.name))


def process(specimen, proto, layers_root, comp_run_dir, out_root, data_root, qp_dir, verify_subset, run_child_timeout_s=None,
            resume: bool = False):  # pragma: no cover - real data
    import build_compartment_adata as B
    import rebuild_single_cell as RB

    M = json.loads((Path(comp_run_dir) / "M.json").read_text())["M"]
    layers, genes, qp_full, qp, nuc, umi, coords = B.load_inputs(Path(layers_root) / specimen, Path(qp_dir) / RB.SPECIMENS[specimen]["qp_csv"])
    target = Path(out_root) / specimen
    if resume:
        # W1: a specimen with a finalized dir (no .tmp_* remains, e.g. a whole-array RESUME resubmit) must not be recomputed; load and
        # return the already-finalized meta. main() lets such a target through check_output_target(allow_finalized=True).
        if (target / "specimen.json").is_file():
            meta = json.loads((target / "specimen.json").read_text())
            print(specimen, "already finalized — skipping resume", flush=True)
            return meta
        matches = _stale_tmp_dirs(out_root, specimen)
        if len(matches) == 1:
            tmp = matches[0]
        elif len(matches) == 0:
            tmp = Path(out_root) / f".tmp_{specimen}_{os.getpid()}"   # nothing to reuse -> full run
        else:
            raise RuntimeError(f"resume ambiguity: {len(matches)} .tmp dirs for {specimen}: {matches}")
    else:
        # W3: a crash before finalize leaves .tmp_<spec>_<pid> holding already-computed ok runs. Never delete them; fail loudly.
        _stale = _stale_tmp_dirs(out_root, specimen)
        if _stale:
            raise RuntimeError(f"stale staged dir(s) for {specimen}: {[str(d) for d in _stale]}; re-run with --resume to reuse them, "
                               "or move them aside manually")
        tmp = Path(out_root) / f".tmp_{specimen}_{os.getpid()}"
    meta = run_matrix(specimen, proto["_protocol_sha256"], layers, genes, qp, nuc, umi, coords, M, tmp, verify_subset, run_child_timeout_s=run_child_timeout_s,
                      resume=resume)
    if target.exists():
        target.rmdir()
    os.rename(tmp, target)
    return meta


def toy_main(out: Path, dis_thr: float = 3.0, run_child_timeout_s: int | None = None) -> int:
    """Toy specimen through the same matrix code (nx 12 x ny 10 grid; used by the unit test in envs/COMMOT)."""
    import positive_control as PC

    layers, genes, qp, nuc, umi, coords, M = PC.toy_inputs(nx=12, ny=10)
    meta = run_matrix("TOY", "toy", layers, genes, qp, nuc, umi, coords, M, Path(out), verify_subset=True, dis_thr_primary=dis_thr,
                      run_child_timeout_s=run_child_timeout_s)
    print("TOY_DONE", json.dumps({r["run_id"]: r["status"] for r in meta["runs"]}))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--toy", type=Path, default=None, help="run the synthetic toy specimen into this directory (test only)")
    ap.add_argument("--specimen")
    ap.add_argument("--layers-root", type=Path)
    ap.add_argument("--compartments-run", type=Path, help="A12.7 output dir (holds M.json)")
    ap.add_argument("--out-root", type=Path, help="output/commot_runs_module3_7type/<sha12>/<run_id>")
    ap.add_argument("--run-id")
    ap.add_argument("--data-root", type=Path)
    ap.add_argument("--qp-dir", type=Path, default=None, help="the stamped B1 dir; REQUIRED for a v42 protocol, legacy default output/module3_7type otherwise")
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    ap.add_argument("--verify-subset", action="store_true")
    ap.add_argument("--resume", action="store_true", help="resume an existing staged .tmp matrix, reusing already-ok run_ids")
    ap.add_argument("--run-child-timeout-s", type=int, default=None, help="per-run child wall-clock budget (default RUN_CHILD_TIMEOUT_S / env)")
    a = ap.parse_args(argv)
    if a.toy:
        return toy_main(a.toy, run_child_timeout_s=a.run_child_timeout_s)
    for req in ("specimen", "layers_root", "compartments_run", "out_root", "run_id", "data_root"):
        if getattr(a, req) is None:
            ap.error(f"--{req.replace('_', '-')} is required")
    proto = PIO.load_protocol(a.protocol)
    PF.guard_entry("B6", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.qp_dir, layers_run=a.layers_root, comp_run=a.compartments_run)
    a.qp_dir = PF.legacy_or_required(a.qp_dir, a.protocol, "--qp-dir", REPO_ROOT / "output" / "module3_7type")
    if tuple(a.out_root.parts[-2:]) != (proto["_protocol_sha12"], a.run_id):
        raise PIO.OutputTargetError(f"out-root must end in <{proto['_protocol_sha12']}>/<{a.run_id}>")
    PIO.check_output_target(a.out_root / a.specimen, proto["_protocol_sha12"], a.run_id, allow_finalized=a.resume)
    a.out_root.mkdir(parents=True, exist_ok=True)
    process(a.specimen, proto, a.layers_root, a.compartments_run, a.out_root, a.data_root, a.qp_dir, a.verify_subset, a.run_child_timeout_s,
            resume=a.resume)
    return 0


if __name__ == "__main__":
    sys.exit(main())
