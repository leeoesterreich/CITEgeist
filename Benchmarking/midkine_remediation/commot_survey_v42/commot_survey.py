"""Stage B26: all-pair COMMOT survey at two levels, and the check 19 ranking of cancer-sent pathways. Package settings only; nothing tuned.

    python -u commot_survey.py --level l1 --specimen <name|index>   (array task: level 1, one specimen, stamped)
    python -u commot_survey.py --level l2 --specimen <name|index>   (array task: level 2, one specimen, stamped)
    python -u commot_survey.py --finalize                            (dependent job: check 19 ranking, stamped, _COMPLETE.json)

Level 1 (CITEgeist): the stage B6 R00 setting rerun with EVERY sender/receiver column kept. Input: the B5 compartment h5ad
(output/compartments_module3_7type/<protocol>/<B5 id>/<specimen>.h5ad), the frozen LR list (commot_runs.freeze_lr, CellChat human,
heteromeric, min_cell_pct 0.01), commot_runs.commot_run with dis_thr 500. Rows are spot|cell-type compartments.
Level 2 (raw Visium spots, no deconvolution): all in-tissue spots of the GEO processed files, normalize_total + log1p, the same settings as
check 18 run_commot (copied from exploratory_v42/check18_spot_commot_keratin.py, which copied them from commot_all_samples.py).
Check 19 (rule: exploratory_v42/check19_cancer_sent_pathways_DECISION_RULE.md): cancer-sent pathway signal = sum of each pathway's sender score over
the cancer rows (level 1: compartments with cell_type Cancer_Cells; level 2: spots with Cancer_Cells share >= 0.2 in the B1 proportions);
rank pathways within the specimen (1 = largest). Descriptive only.

The level-2 artifacts pin the B1 proportions as their stamp input: the stage's cancer-spot rule (finalize) is defined on the B1 barcodes, and the
raw Visium files (outside the repo) are hashed under external_inputs_sha256. Stage output dir = results/<B26 run id of v42_run_ids.yaml>
(never reused: a rerun needs a NEW id). Run from a clean detached worktree by .commands/submit-commot-survey-stage.sh; data are read from
$V2_NAMESPACE_ROOT (main checkout). Heavy/optional modules (commot, h5py, anndata, scanpy, v42_stamp) are imported inside the functions that need them.
"""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "Benchmarking/midkine_remediation/src"))
sys.path.insert(0, str(REPO / "scripts"))

SPECIMENS = ["HCC22-088-P1-S1", "HCC22-088-P1-S2", "HCC22-088-P2-S1", "HCC22-088-P2-S2",
             "HCC22-088-P3-S1_A", "HCC22-088-P3-S2", "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
             "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep", "HCC22-088-P6-S1", "HCC22-088-P6-S2_D"]
LEVELS = ("l1", "l2")
NS = Path(os.environ.get("V2_NAMESPACE_ROOT", REPO))
STAGE = "B26"
RUN_IDS = REPO / "Benchmarking/midkine_remediation/v42_run_ids.yaml"
REL_RESULTS = "Benchmarking/midkine_remediation/commot_survey_v42/results"
DATA_ROOT = Path(str(_cg.root("DATA") / "GSE289326" / "processed_files"))
DIS_THR = 500.0        # B6 R00 / commot_all_samples.py
MIN_CELL_PCT = 0.01
CANCER = "Cancer_Cells"
CANCER_MIN = 0.2       # level-2 cancer-spot rule (check 18 / check 19)
MDK_PAIRS = ["MDK-LRP1", "MDK-SDC4", "MDK-NCL", "MDK-SDC2", "MDK-ITGA6_ITGB1"]
PATHWAY = "MK"


# ----------------------------------------------------------------------------------------------- registry / paths
def load_registry():
    """(B26 run id, protocol sha12, {stage: run id}, current V2_TS) from v42_run_ids.yaml (one protocol key)."""
    import yaml

    d = yaml.safe_load(RUN_IDS.read_text()) or {}
    runs = d.get("runs") or {}
    assert len(runs) == 1, sorted(runs)
    sha12, ids = next(iter(runs.items()))
    rid = ids.get(STAGE)
    assert rid and rid != "NOT_RUN", f"stage {STAGE} has no run id in {RUN_IDS}"
    return str(rid), str(sha12), {k: str(v) for k, v in ids.items()}, str(d["current_v2_ts"])


def props_path(spec, ids, ts):
    return NS / f"output/module3_rna_refined_7type_{ts}/{ids['B1']}/{spec}_cell_prop_finetuned_results.csv"


def h5ad_path(spec, sha12, ids):
    return NS / f"output/compartments_module3_7type/{sha12}/{ids['B5']}/{spec}.h5ad"


def task_files(level, spec):
    names = ["sender.csv.gz", "receiver.csv.gz", "lr.csv"] + (["obs.csv.gz", "validation.json"] if level == "l1" else ["xy.csv.gz"])
    return [f"{spec}_{n}" for n in names]


def resolve_specimen(s):
    spec = SPECIMENS[int(s)] if str(s).isdigit() else s
    assert spec in SPECIMENS, f"unknown specimen {spec!r}"
    return spec


# ----------------------------------------------------------------------------------------------- ranking (pure)
def lr_names(lr):
    """(pair names 'LIGAND-RECEPTOR', pathway names) of a frozen LR table (columns: ligand, receptor, pathway, ...)."""
    pairs = [f"{r.iloc[0]}-{r.iloc[1]}" for _, r in lr.iterrows()]
    return pairs, sorted(set(lr.iloc[:, 2].astype(str)))


def sent_totals(sender, names, mask=None):
    """Sum of each ``s-<name>`` sender column over the masked rows (all rows if mask is None), largest first. Names absent from the table are dropped."""
    cols = [f"s-{n}" for n in dict.fromkeys(names) if f"s-{n}" in sender.columns]
    sub = sender if mask is None else sender.loc[np.asarray(mask, bool)]
    x = sub[cols].sum()
    x.index = [c[2:] for c in cols]
    return x.sort_values(ascending=False, kind="stable")


def _rank(x, name):
    return int(x.rank(ascending=False, method="min")[name]) if name in x.index else None


def rank_summary(sender, pair_names, pathway_names, cancer_mask, top_n=5):
    """The check 19 numbers for one specimen and level. ``cancer_mask``: boolean per sender row (cancer compartments / cancer spots)."""
    cancer_mask = np.asarray(cancer_mask, bool)
    pw = sent_totals(sender, pathway_names, cancer_mask)
    tot = float(pw.sum())
    share = pw / tot if tot > 0 else pw * np.nan
    pr = sent_totals(sender, pair_names, cancer_mask)
    pr_all = sent_totals(sender, pair_names)
    pw_all = sent_totals(sender, pathway_names)
    return {
        "n_rows": int(len(sender)),
        "n_cancer_rows": int(cancer_mask.sum()),
        "MK_rank": _rank(pw, PATHWAY),
        "MK_share": float(share[PATHWAY]) if PATHWAY in share.index else None,
        "top5": [[k, float(v)] for k, v in share.head(top_n).items()],
        "MK_rank_all_celltypes": _rank(pw_all, PATHWAY),
        "MDK_pair_rank_cancer": {p: _rank(pr, p) for p in MDK_PAIRS},
        "MDK_pair_rank_all_celltypes": {p: _rank(pr_all, p) for p in MDK_PAIRS},
        "n_pairs_in_lr": int(len(dict.fromkeys(pair_names))),
        "n_pathways_in_lr": int(len(pathway_names)),
        "n_pairs_ranked": int(len(pr)),
        "n_pathways_ranked": int(len(pw)),
    }


def top5_text(top5):
    return " | ".join(f"{k} {v:.2f}" for k, v in top5)


def short(spec):
    return spec.replace("HCC22-088-", "")


# ----------------------------------------------------------------------------------------------- stamping
def stamp_kw(rid, sha12, started):
    prod = str(Path(__file__).resolve().relative_to(REPO))
    return dict(run_id=rid, protocol_sha12=sha12, producer=prod, started_at=started, repo_root=NS, git_root=REPO, stage_role="commot_survey")


def stamp_files(paths, inputs, rid, sha12, started, external=(), last_clean=True):
    """Sidecar for every path. The whole-repo clean-tree check (1-2 min on NFS) runs at job start (preflight) and for the LAST stamp only."""
    import v42_stamp as W

    kw = stamp_kw(rid, sha12, started)
    for i, p in enumerate(paths):
        last = i == len(paths) - 1
        W.write_stamp(p, inputs=inputs, external_inputs=external, **dict(kw, require_clean_tree=bool(last and last_clean)))


def preflight(out_paths):
    import v42_stamp as W

    W.preflight(str(Path(__file__).resolve().relative_to(REPO)), out_paths=out_paths, repo_root=NS, git_root=REPO)


# ----------------------------------------------------------------------------------------------- level 1
def run_l1(spec, out, sha12, ids):
    import anndata as ad
    import commot as ct
    import commot_runs as CR
    import scipy.sparse as sp

    t0 = time.time()
    a = ad.read_h5ad(h5ad_path(spec, sha12, ids))
    db = ct.pp.ligand_receptor_database(species="human", database="CellChat")
    need = CR.lr_gene_set(db)
    genes = list(a.var_names)
    gi = np.array([i for i, g in enumerate(genes) if g in need])
    lg = [genes[i] for i in gi]
    X = sp.csr_matrix(a.X[:, gi])
    names = list(a.obs_names)
    xy = np.asarray(a.obsm["spatial"], float)
    lr = CR.freeze_lr(X, xy, lg, names, MIN_CELL_PCT)
    print(spec, "rows", len(names), "LR", len(lr), flush=True)
    res = CR.commot_run(X, xy, lg, names, lr, DIS_THR)
    res["sender"].index = names
    res["receiver"].index = names
    f = {n: out / f"{spec}_{n}" for n in ("sender.csv.gz", "receiver.csv.gz", "lr.csv", "obs.csv.gz", "validation.json")}
    res["sender"].to_csv(f["sender.csv.gz"])
    res["receiver"].to_csv(f["receiver.csv.gz"])
    lr.to_csv(f["lr.csv"], index=False)
    a.obs[["spot_id", "barcode", "cell_type", "qp_prop", "n_nuclei_spot", "eligible_cancer"]].to_csv(f["obs.csv.gz"])
    # comparison with the stored B6 R00 MDK sender columns (recorded, not a gate: the stored file is float32-rounded and unstamped)
    old_p = NS / f"output/commot_runs_module3_7type/{sha12}/{ids['B6']}/{spec}/R00/sender.csv.gz"
    diffs = {}
    if old_p.is_file():
        old = pd.read_csv(old_p, index_col=0)
        for c in [c for c in old.columns if c.startswith("s-MDK")]:
            if c in res["sender"].columns:
                diffs[c] = float(np.max(np.abs(old[c].to_numpy() - res["sender"][c].reindex(old.index).to_numpy())))
    f["validation.json"].write_text(json.dumps({"max_abs_diff_vs_B6_R00": diffs, "b6_sender": str(old_p), "seconds": time.time() - t0,
                                                "n_lr": int(len(lr)), "n_rows": len(names), "dis_thr": DIS_THR, "min_cell_pct": MIN_CELL_PCT}, indent=1))
    print("done", diffs, time.time() - t0, flush=True)
    return list(f.values()), [h5ad_path(spec, sha12, ids)], []


# ----------------------------------------------------------------------------------------------- level 2
def read_raw(spec):
    """All in-tissue spots with finite coords: (X spots x genes csr, barcodes, unique gene names, xy=(pxl_col, pxl_row)). Copy of check 18 read_raw."""
    import h5py
    import scipy.sparse as sp

    outs = DATA_ROOT / spec / "outs"
    with h5py.File(outs / "filtered_feature_bc_matrix.h5", "r") as f:
        m = f["matrix"]
        bc = np.array([b.decode() for b in m["barcodes"][:]])
        ft = m["features"]["feature_type"][:]
        names = np.array([b.decode() for b in m["features"]["name"][:]])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])).tocsr()   # genes x spots
    ge = np.where(ft == b"Gene Expression")[0]
    X, names = X[ge], names[ge]
    seen, uniq = {}, []
    for n in names:                              # same name-1, name-2 rule as the model's var_names_make_unique
        k = seen.get(n, 0)
        uniq.append(n if k == 0 else f"{n}-{k}")
        seen[n] = k + 1
    tp = pd.read_csv(outs / "spatial" / "tissue_positions.csv", index_col=0)
    tp = tp.loc[tp["in_tissue"] == 1]
    xy = np.column_stack([tp["pxl_col_in_fullres"], tp["pxl_row_in_fullres"]]).astype(float)
    tp = tp.loc[np.isfinite(xy).all(axis=1)]
    xy = xy[np.isfinite(xy).all(axis=1)]
    bi = pd.Index(bc).get_indexer(tp.index)
    assert (bi >= 0).all(), "in-tissue barcodes absent from the h5"
    return X[:, bi].T.tocsr().astype(np.float32), np.asarray(tp.index), np.array(uniq), xy


def run_l2(spec, out, sha12, ids, ts):
    import anndata as ad
    import commot as ct
    import scanpy as sc

    t0 = time.time()
    X, bc, genes, xy = read_raw(spec)
    a = ad.AnnData(X=X.copy(), var=pd.DataFrame(index=genes))     # identical to check 18 run_commot, but keeping all sender AND receiver columns
    a.obs_names = [f"s{i}" for i in range(a.n_obs)]
    a.obsm["spatial"] = xy
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    db = ct.pp.ligand_receptor_database(species="human", database="CellChat")
    lr = ct.pp.filter_lr_database(db, a, min_cell_pct=MIN_CELL_PCT)
    print(spec, "spots", a.n_obs, "LR", len(lr), flush=True)
    ct.tl.spatial_communication(a, database_name="cellchat", df_ligrec=lr, dis_thr=DIS_THR, heteromeric=True, pathway_sum=True)
    s = a.obsm["commot-cellchat-sum-sender"].copy()
    r = a.obsm["commot-cellchat-sum-receiver"].copy()
    s.index = bc
    r.index = bc
    f = {n: out / f"{spec}_{n}" for n in ("sender.csv.gz", "receiver.csv.gz", "lr.csv", "xy.csv.gz")}
    s.to_csv(f["sender.csv.gz"])
    r.to_csv(f["receiver.csv.gz"])
    lr.to_csv(f["lr.csv"], index=False)
    pd.DataFrame(xy, index=bc, columns=["x", "y"]).to_csv(f["xy.csv.gz"])
    pb = pd.read_csv(props_path(spec, ids, ts), index_col=0).index
    print(f"B1 proportion rows {len(pb)}, in the raw in-tissue set {int(pb.isin(bc).sum())}; done {time.time() - t0:.0f}s", flush=True)
    ext = [DATA_ROOT / spec / "outs/filtered_feature_bc_matrix.h5", DATA_ROOT / spec / "outs/spatial/tissue_positions.csv"]
    return list(f.values()), [props_path(spec, ids, ts)], ext


# ----------------------------------------------------------------------------------------------- array task
def run_task(level, spec, rid, sha12, ids, ts):
    assert level in LEVELS, level
    out = NS / REL_RESULTS / rid / level
    arts = [out / n for n in task_files(level, spec)]
    existing = [p for p in arts if p.exists() or Path(str(p) + ".stamp.json").exists()]
    if existing:
        raise FileExistsError(f"{level} {spec}: outputs exist, a rerun needs a NEW {STAGE} run id: {[str(p) for p in existing]}")
    import v42_stamp as W

    preflight(arts)
    started = W.utc_now()
    out.mkdir(parents=True, exist_ok=True)
    written, inputs, ext = run_l1(spec, out, sha12, ids) if level == "l1" else run_l2(spec, out, sha12, ids, ts)
    assert sorted(map(str, written)) == sorted(map(str, arts)), (written, arts)
    stamp_files(written, inputs, rid, sha12, started, ext)
    print(f"{STAGE} {level} {spec} stamped", flush=True)


# ----------------------------------------------------------------------------------------------- finalize (check 19)
def cancer_mask(level, spec, sender, out, ids, ts):
    if level == "l1":
        obs = pd.read_csv(out / "l1" / f"{spec}_obs.csv.gz", index_col=0)
        m = (obs["cell_type"].reindex(sender.index) == CANCER).to_numpy()
        assert (sender.index.str.endswith(f"|{CANCER}") == m).all(), "obs cell_type and the row-name suffix disagree"
        return m
    share = pd.read_csv(props_path(spec, ids, ts), index_col=0)[CANCER].reindex(sender.index)
    return (share >= CANCER_MIN).to_numpy()


def compute_check19(out, ids, ts):
    summary, rows = {}, []
    for level in LEVELS:
        summary[level] = {}
        for spec in SPECIMENS:
            sender = pd.read_csv(out / level / f"{spec}_sender.csv.gz", index_col=0)
            pairs, pathways = lr_names(pd.read_csv(out / level / f"{spec}_lr.csv"))
            d = rank_summary(sender, pairs, pathways, cancer_mask(level, spec, sender, out, ids, ts))
            summary[level][spec] = d
            rows.append({"level": level, "spec": short(spec), "n_cancer_rows": d["n_cancer_rows"], "MK_rank": d["MK_rank"],
                         "MK_share": d["MK_share"], "top5": top5_text(d["top5"])})
            print(level, short(spec), d["n_cancer_rows"], d["MK_rank"], d["MK_share"], flush=True)
    return summary, pd.DataFrame(rows)


def finalize(rid, sha12, ids, ts):
    import v42_stamp as W

    out = NS / REL_RESULTS / rid
    final = [out / "check19_cancer_sent.csv", out / "check19_summary.json", out / "_COMPLETE.json"]
    existing = [p for p in final if p.exists()]
    if existing:
        raise FileExistsError(f"{[str(p) for p in existing]} exist; a rerun needs a NEW {STAGE} run id")
    problems = []
    for level in LEVELS:
        for spec in SPECIMENS:
            for n in task_files(level, spec):
                p = out / level / n
                if not p.is_file():
                    problems.append(f"missing {p}")
                    continue
                try:
                    W.verify_stamp(p, run_id=rid, protocol_sha12=sha12, repo_root=NS)
                except Exception as e:
                    problems.append(f"{p.name}: {e}")
    if problems:
        raise FileNotFoundError(f"finalize refused ({len(problems)} problems): " + " | ".join(problems))
    preflight(final)
    started = W.utc_now()
    summary, table = compute_check19(out, ids, ts)
    table.to_csv(final[0], index=False)
    final[1].write_text(json.dumps({"stage": STAGE, "run_id": rid, "rule": "exploratory_v42/check19_cancer_sent_pathways_DECISION_RULE.md",
                                    "cancer_min_share_l2": CANCER_MIN, "levels": summary}, indent=1))
    inputs = [out / lv / f"{s}_{n}" for lv in LEVELS for s in SPECIMENS for n in ("sender.csv.gz", "lr.csv")]
    inputs += [out / "l1" / f"{s}_obs.csv.gz" for s in SPECIMENS] + [props_path(s, ids, ts) for s in SPECIMENS]
    stamp_files(final[:2], inputs, rid, sha12, started)
    final[2].write_text(json.dumps({"status": "COMPLETE", "stage": STAGE, "run_id": rid, "slurm_job_id": os.environ.get("SLURM_JOB_ID")}, indent=2))
    print(f"{STAGE} stamped {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", choices=LEVELS)
    ap.add_argument("--specimen", help="array task: specimen name or 0-based index (with --level)")
    ap.add_argument("--finalize", action="store_true", help="after all 24 tasks: check 19 ranking, stamps, _COMPLETE.json")
    a = ap.parse_args()
    if a.finalize == (a.specimen is not None):
        ap.error("give --level and --specimen (an array task) OR --finalize")
    if a.specimen is not None and not a.level:
        ap.error("--specimen needs --level")
    rid, sha12, ids, ts = load_registry()
    if a.finalize:
        return finalize(rid, sha12, ids, ts)
    return run_task(a.level, resolve_specimen(a.specimen), rid, sha12, ids, ts)


if __name__ == "__main__":
    main()
