"""Mechanical R1a/R1b/R2/R3 readouts defined in MDK_PROGRAM_DECISION_RULE.md (no definition changed).
    python mdk_readout.py <run_dir>      -> <run_dir>/mdk_readout.json and mdk_readout.md
Implementation choices where the rule is silent (all stated in the output):
  - rank = ordinal rank by loading, 1 = largest, ties broken by gene order (stable sort); same rule for MDK and null genes.
  - gene universe for ranks = the program's own gene_names; E = genes with nonzero total in the specimen's cancer layer,
    intersected with those gene_names.
  - matching bins = 10 equal-frequency bins (pandas qcut on rank of mean expression) of mean cancer-layer expression over E;
    null = 1000 genes drawn without replacement from MDK's bin (excluding MDK), numpy default_rng(0), one fresh rng per specimen.
  - MDK program = cancer program with the largest MDK loading; if that loading is 0 the specimen is reported missing.
  - R1b universe = E; top-200 = top 200 genes of the MDK program by loading among E; hypergeom sf(k-1, |E|, |S in E|, 200).
  - R2 group = aligned group whose per-sample signature for that specimen equals the MDK program's aligned loadings;
    n specimens = len(samples_present); MDK consensus rank among all union genes; null = same matched genes (in the union).
  - R3 = pairs in m4b.significant_pairs (package FDR 0.05) involving the MDK program (anchor Cancer_Cells, program id).
"""
import json, os, pickle, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import hypergeom

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
# data live in the main checkout; a stamped run executes from a clean worktree (V2_NAMESPACE_ROOT = main checkout)
LAYERS = Path(os.environ.get("V2_NAMESPACE_ROOT", REPO)) / "output/sace_layers_module3_7type/fc43dbea67bf/v2b41"
SEC = ["HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4", "PDIA6", "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1"]
CT = "Cancer_Cells"
run = Path(sys.argv[1])
SPECS = sorted(p.stem for p in (run / "per_specimen").glob("*.pkl"))
m5 = pickle.loads((run / "module5" / "module5.pkl").read_bytes())


def ordinal_rank(v):
    order = np.argsort(-v, kind="stable")
    r = np.empty(len(v), int)
    r[order] = np.arange(1, len(v) + 1)
    return r


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(len(p), np.nan)
    pv = p[ok]
    n = len(pv)
    o = np.argsort(pv)
    adj = pv[o] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    qq = np.empty(n)
    qq[o] = np.minimum(adj, 1)
    q[ok] = qq
    return q


rows = []
for s in SPECS:
    rec = pickle.loads((run / "per_specimen" / f"{s}.pkl").read_bytes())
    r = rec["m4"].results_by_anchor[CT]
    genes = list(r.gene_names)
    gi = {g: i for i, g in enumerate(genes)}
    W = np.asarray(r.W)
    row = {"specimen": s, "n_programs_cancer": W.shape[1], "n_allzero_cancer": int((np.abs(W).sum(0) == 0).sum())}
    if "MDK" not in gi or W[gi["MDK"]].max() <= 0:
        row["missing"] = "MDK absent or zero in every cancer program"
        rows.append(row)
        continue
    k = int(np.argmax(W[gi["MDK"]]))  # argmax -> first (lowest index) on ties
    row["mdk_program"] = k
    col = W[:, k]
    # expression of cancer layer
    d = LAYERS / s
    lg = [l for l in (d / "genes.txt").read_text().splitlines() if l.strip()]
    L = sp.load_npz(d / f"layer_{CT}.npz").tocsr()
    tot = np.asarray(L.sum(0)).ravel()
    mean = tot / L.shape[0]
    lmap = {g: i for i, g in enumerate(lg)}
    Eidx = np.array([gi[g] for g in genes if g in lmap and tot[lmap[g]] > 0 and g in gi])  # indices into program genes
    mexp = np.array([mean[lmap[genes[i]]] for i in Eidx])
    row["n_expressed"] = int(len(Eidx))
    mdk_i = gi["MDK"]
    assert mdk_i in set(Eidx.tolist()), "MDK not expressed in cancer layer"
    bins = pd.qcut(pd.Series(mexp).rank(method="first"), 10, labels=False).to_numpy()
    mb = bins[list(Eidx).index(mdk_i)]
    pool = np.array([i for i, b in zip(Eidx, bins) if b == mb and i != mdk_i])
    rng = np.random.default_rng(0)
    null = rng.choice(pool, size=1000, replace=False)
    rk = ordinal_rank(col)
    mrank = int(rk[mdk_i])
    row.update(mdk_rank=mrank, mdk_loading=float(col[mdk_i]), mdk_rank_pct=mrank / len(genes), bin_pool=int(len(pool)))
    row["R1a_p"] = (1 + int((rk[null] <= mrank).sum())) / 1001
    # R1b
    Eg = [genes[i] for i in Eidx]
    Eset = set(Eg)
    sec_in_E = [g for g in SEC if g in Eset]
    Erk = Eidx[np.argsort(-col[Eidx], kind="stable")][:200]
    top = [genes[i] for i in Erk]
    hit = [g for g in SEC if g in top]
    row.update(n_sec_in_E=len(sec_in_E), n_sec_top200=len(hit), sec_top200=",".join(hit))
    row["R1b_p"] = float(hypergeom.sf(len(hit) - 1, len(Eidx), len(sec_in_E), 200))
    row["sec_ranks"] = {g: (int(rk[gi[g]]) if g in gi else None) for g in SEC}
    # R2
    allg = sorted({g for t in range(1) for g in m5.aligned_programs[0].consensus_signature})
    aidx = {g: i for i, g in enumerate(allg)}
    found = None
    for a in m5.aligned_programs:
        sig = a.program_signatures.get((s, CT, k))  # keyed by (sample, cell_type, program): ids restart per cell type
        if sig is None:
            continue
        sg = np.array([sig[aidx[g]] for g in genes])
        if np.allclose(sg, col, atol=0, rtol=0):
            found = a
            break
    if found is None:
        row["R2_missing"] = "no aligned group holds this program (sample, cell_type, program) key"
    else:
        cons = np.array([found.consensus_signature[g] for g in allg])
        crk = ordinal_rank(cons)
        cm = int(crk[aidx["MDK"]])
        cn = np.array([crk[aidx[genes[i]]] for i in null])
        row.update(R2_group=found.program_id, R2_group_cell_type=found.cell_type, R2_n_specimens=len(found.samples_present),
                   R2_conservation=found.conservation_score, R2_mdk_consensus_rank=cm,
                   R2_p=(1 + int((cn <= cm).sum())) / 1001, R2_top_genes=",".join(found.top_genes[:10]),
                   # descriptive only, not used by any rule: programs per cell type in the R2 group
                   R2_group_cell_types={ct_: int(n_) for ct_, n_ in sorted(pd.Series([kk[1] for kk in found.program_keys]).value_counts().items())})
    # R3
    prs = []
    for q in rec["m4b"].significant_pairs:
        if (q.anchor1 == CT and q.program1_id == k) or (q.anchor2 == CT and q.program2_id == k):
            other = (q.anchor2, q.program2_id) if (q.anchor1 == CT and q.program1_id == k) else (q.anchor1, q.program1_id)
            prs.append({"partner_type": other[0], "partner_program": int(other[1]), "I": float(q.bivariate_morans_i), "p": float(q.bivariate_pvalue)})
    row["R3_n_fdr_pairs"] = len(prs)
    row["R3_pairs"] = prs
    row["R3_n_fdr_pairs_all_in_specimen"] = int(rec["m4b"].n_significant)
    rows.append(row)

q = bh([r.get("R1b_p", np.nan) for r in rows])
for r, qq in zip(rows, q):
    if "R1b_p" in r:
        r["R1b_q"] = float(qq)
        r["R1_both"] = bool(r["R1a_p"] < 0.05 and qq < 0.05)
n_both = sum(1 for r in rows if r.get("R1_both"))
spans = [r for r in rows if r.get("R2_n_specimens") is not None]
summary = {
    "n_specimens": len(rows), "n_specimens_computable": sum(1 for r in rows if "R1b_p" in r),
    "R1a_p_lt_0.05": sum(1 for r in rows if r.get("R1a_p", 1) < 0.05),
    "R1b_q_lt_0.05": sum(1 for r in rows if r.get("R1b_q", 1) < 0.05),
    "R1_both_count": n_both,
    "condition_MDK_secretory_load_together_met": bool(n_both >= 6),
    "R2_groups_spanning_ge6": sum(1 for r in spans if r["R2_n_specimens"] >= 6),
    "condition_conserved_met_per_specimen_group": {r["specimen"]: bool(r["R2_n_specimens"] >= 6) for r in spans},
}
(run / "mdk_readout.json").write_text(json.dumps({"summary": summary, "rows": rows}, indent=1, default=str))
cols = ["specimen", "mdk_program", "mdk_rank", "mdk_loading", "R1a_p", "n_sec_top200", "R1b_p", "R1b_q", "R1_both",
        "R2_n_specimens", "R2_conservation", "R2_mdk_consensus_rank", "R2_p", "R3_n_fdr_pairs"]
df = pd.DataFrame(rows)
for c in cols:
    if c not in df:
        df[c] = np.nan
md = ["# MDK program readout (rule: MDK_PROGRAM_DECISION_RULE.md)", "", df[cols].to_markdown(index=False, floatfmt=".4g"), "",
      "## Summary", "", "```", json.dumps(summary, indent=1), "```", "", "## R3 pairs", ""]
for r in rows:
    md.append(f"- {r['specimen']}: " + (json.dumps(r.get("R3_pairs")) if "R3_pairs" in r else r.get("missing", "")))
md += ["", "## R2 group cell-type composition (descriptive only; count of programs per cell type in the group; not used by any rule)", ""]
for r in rows:
    md.append(f"- {r['specimen']}: " + (json.dumps(r.get("R2_group_cell_types")) if "R2_group_cell_types" in r else r.get("R2_missing", r.get("missing", ""))))
md += ["", "## Implementation notes", "", __doc__]
(run / "mdk_readout.md").write_text("\n".join(md))
print("\n".join(md[:40]))
