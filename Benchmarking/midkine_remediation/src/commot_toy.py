"""Tiny synthetic COMMOT runs for gates G4/G5/G11. Run with the COMMOT env python (commot 0.0.3 needs numpy < 2).

The protocol names CITEgeist_env for commot, but ``import commot`` raises ``AttributeError: np.Inf was removed in
NumPy 2.0`` there (numpy 2.1.3); the working install is envs/COMMOT (numpy 1.26.4, commot 0.0.3).
Prints one JSON document on stdout.
"""
import json
import sys
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")
import anndata as ad  # noqa: E402
import commot as ct  # noqa: E402
from scipy.spatial import distance_matrix  # noqa: E402

LIG, RECS = "MDK", ["NCL", "SDC2", "SDC4", "LRP1", "ITGA6", "ITGB1"]
GENES = [LIG] + RECS + ["OTHER"]
PAIRS = [("MDK", "NCL"), ("MDK", "SDC2"), ("MDK", "SDC4"), ("MDK", "LRP1"), ("MDK", "ITGA6_ITGB1")]


def toy(rows, coords):
    X = np.zeros((len(rows), len(GENES)))
    for i, (lig, rec) in enumerate(rows):
        X[i, GENES.index(LIG)] = lig
        for r in RECS:
            X[i, GENES.index(r)] = rec
        X[i, -1] = 1.0
    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=[f"r{i}" for i in range(len(rows))]), var=pd.DataFrame(index=GENES))
    a.obsm["spatial"] = np.asarray(coords, float)
    return a


def lr_frame():
    db = ct.pp.ligand_receptor_database(species="human", database="CellChat")
    keep = [(l, r) for l, r in PAIRS]
    m = np.array([(row.iloc[0], row.iloc[1]) in keep for _, row in db.iterrows()])
    return db[m].reset_index(drop=True)


def run(a, df, dis_thr, dmat=None):
    b = a.copy()
    if dmat is not None:
        b.obsp["spatial_distance"] = dmat
    ct.tl.spatial_communication(b, database_name="cellchat", df_ligrec=df, dis_thr=dis_thr, heteromeric=True, pathway_sum=True)
    return b


def summarize(b):
    keys = [k for k in b.obsp.keys() if k.startswith("commot-cellchat-") and k.count("-") >= 3 and not k.endswith("total-total")]
    sender = b.obsm["commot-cellchat-sum-sender"]
    rs = {k: np.asarray(b.obsp[k].sum(axis=1)).ravel().tolist() for k in keys}
    S = {k: b.obsp[k].toarray().tolist() for k in keys if k.endswith("MDK-NCL")}
    return {"obsp_keys": sorted(b.obsp.keys()), "row_sums": rs, "sender": {c: sender[c].tolist() for c in sender.columns}, "S_MDK_NCL": S}


def compute(thr):
    df = lr_frame()
    # three rows: A ligand-only and B receptor-only CO-LOCATED at (0,0); C (both) one pitch away at (1,0)
    a3 = toy([(5.0, 0.0), (0.0, 5.0), (3.0, 3.0)], [[0, 0], [0, 0], [1, 0]])
    b3 = run(a3, df, dis_thr=thr)
    out = {"commot_version": ct.__version__ if hasattr(ct, "__version__") else "0.0.3", "numpy": np.__version__, "lr_pairs": [list(map(str, r[:2])) for r in df.values]}
    out["three_row"] = summarize(b3)
    out["three_row"]["has_spatial_distance"] = "spatial_distance" in b3.obsp
    # R10: same-spot pairs and the diagonal pushed beyond dis_thr
    d = distance_matrix(a3.obsm["spatial"], a3.obsm["spatial"])
    spot = np.array([0, 0, 1])
    same = spot[:, None] == spot[None, :]
    d10 = np.where(same, 1e6, d)
    try:
        b10 = run(a3, df, dis_thr=thr, dmat=d10)
        out["r10_three_row"] = summarize(b10)
        out["r10_three_row"]["has_spatial_distance"] = "spatial_distance" in b10.obsp
    except Exception as e:  # noqa: BLE001  (COMMOT raises when no pair lies inside the cutoff)
        out["r10_three_row"] = {"error": repr(e)[:200]}
    # degenerate two-row toy at identical coordinates (documented as unusable)
    a2 = toy([(5.0, 0.0), (0.0, 5.0)], [[0, 0], [0, 0]])
    try:
        b2 = run(a2, df, dis_thr=thr)
        out["degenerate_two_row"] = {"ran": True, "S_MDK_NCL": summarize(b2)["S_MDK_NCL"]}
    except Exception as e:  # noqa: BLE001
        out["degenerate_two_row"] = {"ran": False, "error": repr(e)[:200]}
    return out


def main():
    """One threshold -> that result; several -> {str(thr): result}. Several thresholds share ONE ``import commot``
    (about 100-220 s of mostly NFS-bound import: scanpy -> anndata -> torch, ot, networkx)."""
    thrs = [float(x) for x in sys.argv[1:]] or [2.0]
    out = compute(thrs[0]) if len(thrs) == 1 else {str(t): compute(t) for t in thrs}
    print("JSON_BEGIN" + json.dumps(out) + "JSON_END")


if __name__ == "__main__":
    main()
