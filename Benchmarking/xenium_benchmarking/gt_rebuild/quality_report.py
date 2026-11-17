"""Descriptive GT-quality report (protein lineage agreement, PanCK vs epithelial, RNA own-lineage gap).

Descriptive only, never a gate: nothing here returns pass/fail and ``integrity_gates`` must not import this module.
Logic is ported from ``reference_audit/gt_arbitration/arbitrate_gt.py`` (a top-level script that cannot be imported).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

from . import gt_config

LIN = {
    "Epithelial": ["PanCK", "E-Cadherin"],
    "Macrophages": ["CD68", "CD163"],
    "T cells": ["CD3E"],
    "B cells": ["CD20"],
    "Endothelial": ["CD31"],
    "Fibroblasts": ["alphaSMA"],
}
MK = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT7", "CDH1", "MUC1", "CA9", "CDH16", "AQP2", "UMOD"],
    "Endothelial": ["PECAM1", "VWF", "EGFL7", "RAMP2", "CLEC14A", "ERG", "SOX17", "CD34", "ACKR1"],
    "Macrophages": ["CD68", "CD163", "MARCO", "MRC1", "MS4A6A", "C1QA", "AIF1", "CD14", "MPEG1", "TREM2"],
    "Fibroblasts": ["DPT", "FBLN1", "FBN1", "PDGFRA", "ASPN", "OGN", "SFRP2", "COL1A1", "DCN", "COL5A2", "MFAP5"],
    "B cells": ["MS4A1", "CD79A", "CD19", "BANK1", "TCL1A", "CD79B"],
    "T cells": ["CD3E", "CD3D", "CD3G", "TRAC", "CD2", "CD247", "IL7R", "CD28"],
}
SCHEMES = ("raw", "log1p")
_T_SPLIT = {"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"}


def _transform(scheme):
    if scheme not in SCHEMES:
        raise ValueError(f"scheme must be one of {SCHEMES}, got {scheme!r}")
    return (lambda x: x) if scheme == "raw" else np.log1p


def _collapse_t(labels: pd.Series) -> pd.Series:
    return labels.replace(_T_SPLIT)


def gmm_thr(values: np.ndarray) -> float:
    """3-component GMM on non-zero values; midpoint of the two lowest means (random_state 42, 200k subsample)."""
    v = np.asarray(values, dtype=float)
    nz = v[v > 0]
    if len(nz) < 50:
        return float(np.percentile(nz, 50)) if len(nz) else 0.0
    rng = np.random.RandomState(0)
    if len(nz) > 200000:
        nz = rng.choice(nz, 200000, replace=False)
    m = np.sort(GaussianMixture(3, random_state=42).fit(nz.reshape(-1, 1)).means_.ravel())
    return float((m[0] + m[1]) / 2)


def _lineage_calls(protein: pd.DataFrame, scheme: str):
    f = _transform(scheme)
    chans = sorted({c for v in LIN.values() for c in v})
    thr = {c: gmm_thr(f(protein[c].values.astype(float))) for c in chans}
    ratio = {c: f(protein[c].values.astype(float)) / max(thr[c], 1e-9) for c in chans}
    score = pd.DataFrame({k: np.max([ratio[c] for c in ch], axis=0) for k, ch in LIN.items()}, index=protein.index)
    pos = score > 1
    call = score.where(pos).idxmax(axis=1).fillna("none")
    return call, pos, thr


def protein_lineage_agreement(labels: pd.Series, protein: pd.DataFrame, scheme: str) -> pd.DataFrame:
    """One-row table: fraction of cells with a protein lineage call, agreement with labels, own-lineage positivity.

    CD4+/CD8+ T labels are collapsed to ``T cells`` (the protein lineage channels do not split T cells).
    """
    protein = protein.loc[labels.index]
    lab = _collapse_t(labels)
    call, pos, _ = _lineage_calls(protein, scheme)
    called = call != "none"
    row = {
        "scheme": scheme,
        "n_cells": int(len(lab)),
        "frac_cells_with_protein_call": float(called.mean()),
        "agree_among_called": float((lab[called] == call[called]).mean()) if called.any() else float("nan"),
        "agree_all_cells": float((lab == call).mean()),
    }
    for lt in LIN:
        m = (lab == lt).values
        row[f"own_pos_{lt}"] = float(pos[lt].values[m].mean()) if m.any() else float("nan")
    return pd.DataFrame([row])


def panck_vs_epithelial_by_region(labels: pd.Series, protein: pd.DataFrame, mapping: pd.DataFrame, scheme="log1p"):
    """Per region: protein PanCK/E-Cadherin positive fractions next to the label Epithelial fraction."""
    f = _transform(scheme)
    protein = protein.loc[labels.index]
    thr = {c: gmm_thr(f(protein[c].values.astype(float))) for c in ("PanCK", "E-Cadherin")}
    panck = f(protein["PanCK"]) > thr["PanCK"]
    epi = panck | (f(protein["E-Cadherin"]) > thr["E-Cadherin"])
    mp = mapping[mapping["spot_idx"] >= 0]
    rows = []
    for r in sorted(mp["region_id"].dropna().unique()):
        ids = mp.index[mp["region_id"] == r].intersection(labels.index)
        rows.append(
            {
                "scheme": scheme,
                "region": int(r),
                "n_cells": int(len(ids)),
                "protein_epi_pos_PanCK_or_Ecad": float(epi[ids].mean()),
                "protein_PanCK_pos": float(panck[ids].mean()),
                "label_epi_frac": float((labels[ids] == "Epithelial").mean()),
            }
        )
    return pd.DataFrame(rows)


def rna_own_lineage_gap(labels: pd.Series, counts: pd.DataFrame) -> pd.DataFrame:
    """Detection of each label's own markers in vs out of the label (gap = in - out).

    Internal coherence only: the labels were derived from this same RNA, so this is NOT independent validation.
    """
    counts = counts.loc[labels.index]
    lab = _collapse_t(labels)
    rows = []
    for lt, genes in MK.items():
        g = [x for x in genes if x in counts.columns]
        m = (lab == lt).values
        if not g or not m.any():
            continue
        det = counts[g].values > 0
        d_in, d_out = det[m].mean(axis=0), det[~m].mean(axis=0)
        rows.append(
            {
                "label": lt,
                "n_markers": len(g),
                "mean_det_in": float(d_in.mean()),
                "mean_det_out": float(d_out.mean()),
                "gap": float((d_in - d_out).mean()),
                "n_cells": int(m.sum()),
            }
        )
    return pd.DataFrame(rows)


def cluster_score_table(names: pd.DataFrame) -> pd.DataFrame:
    """Per-cluster SingleR scores with the argmax score label next to the assigned (fine-tuned) label."""
    out = names.copy()
    sc = [c for c in out.columns if c.startswith("score_")]
    out["argmax_score_label"] = out[sc].idxmax(axis=1).str.replace("score_", "", regex=False)
    return out


def excluded_clusters_table(names: pd.DataFrame) -> pd.DataFrame:
    """Clusters left out of the GT (top type score <= the RECIPE threshold): cluster, label, top_score, n_cells."""
    from .cluster_names import excluded_mask

    t = names.loc[excluded_mask(names)]
    cols = [c for c in ("cluster", "label", "top_score", "n_cells") if c in t.columns]
    return t[cols].reset_index(drop=True)


def write_report(labels, protein, mapping, counts, out_dir, names=None):
    """Write all report CSVs (both schemes, always) under out_dir; return the list of paths."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    files = []
    for s in SCHEMES:
        for name, df in (
            (f"protein_agreement_{s}.csv", protein_lineage_agreement(labels, protein, s)),
            (f"panck_vs_epithelial_{s}.csv", panck_vs_epithelial_by_region(labels, protein, mapping, s)),
        ):
            df.to_csv(out / name, index=False)
            files.append(out / name)
    rna_own_lineage_gap(labels, counts).to_csv(out / "rna_own_lineage_gap.csv", index=False)
    files.append(out / "rna_own_lineage_gap.csv")
    if names is not None:
        cluster_score_table(names).to_csv(out / "cluster_scores.csv", index=False)
        excluded_clusters_table(names).to_csv(out / "excluded_clusters.csv", index=False)
        files += [out / "cluster_scores.csv", out / "excluded_clusters.csv"]
    return files


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labels", required=True, help="label CSV with cell_id,singler_label")
    ap.add_argument("--cluster-names", help="cluster_names.csv (optional)")
    ap.add_argument("--out-dir", default=str(gt_config.OUT_ROOT / "reports"))
    a = ap.parse_args(argv)
    labels = pd.read_csv(a.labels).set_index("cell_id").singler_label
    protein = pd.read_csv(gt_config.OLD_LABELS_DIR / "xenium_cell_proteins.csv", index_col=0)
    counts = pd.read_csv(gt_config.OLD_LABELS_DIR / "xenium_cell_gex.csv", index_col=0)
    mapping = pd.read_csv(gt_config.MAPPING_CSV, index_col=0)
    names = None
    if a.cluster_names:
        from .cluster_names import read_cluster_names

        names = read_cluster_names(a.cluster_names)
    for f in write_report(labels, protein, mapping, counts, a.out_dir, names):
        print(f)


if __name__ == "__main__":
    main()
