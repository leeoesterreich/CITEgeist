#!/usr/bin/env python
"""
Hurdle model for zero-inflated COMMOT sender scores.

Standard Mann-Whitney on all cells misses signal because ~64% of cancer cells
have zero sender scores. This script implements a two-part hurdle model:

  Part 1 (participation): Fisher exact test on proportion of active senders
          (nonzero COMMOT score) between D538G vs WT cancer cells.
  Part 2 (magnitude):     Mann-Whitney U on nonzero sender scores only,
          comparing D538G vs WT among active senders.

Loads the per-cell COMMOT sender scores from the single-cell h5ad produced by
analyze_p4s2_commot.py, applies UMI >= 20 filter, and tests all pathways.

Outputs:
  mdk_analysis/outputs/commot_p4s2/hurdle_results.json        — full results
  mdk_analysis/outputs/commot_p4s2/hurdle_results_summary.csv  — per-pathway table

Usage:
    python -u mdk_analysis/scripts/commot_hurdle_test.py
"""

import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
SAMPLE = "HCC22-088-P4-S2_1i_rep"
SC_H5AD = PROJECT / "output" / "morphology_assignment_v3" / SAMPLE / f"{SAMPLE}_single_cell.h5ad"
OUTPUT_DIR = PROJECT / "mdk_analysis" / "outputs" / "commot_p4s2"

MIN_UMI = 20
BASAL_KERATINS = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT16", "KRT17"]
MDK_PATHWAYS = ["s-MDK-NCL", "s-MDK-SDC2", "s-MDK-SDC4", "s-MDK-LRP1", "s-MDK-ITGA6_ITGB1"]


def load_and_classify(h5ad_path: Path) -> ad.AnnData:
    """Load single-cell h5ad and classify D538G via RNA basal cytokeratin > 0.

    Uses GEX-based classification (any KRT5/14/17 > 0 in cancer cells) to match
    the original hurdle model analysis. The M3.5 protein gate is more conservative
    and yields too few D538G cells for adequate statistical power.
    """
    log.info("Loading %s", h5ad_path)
    adata = sc.read_h5ad(str(h5ad_path))
    adata.var_names_make_unique()

    if "cell_type" in adata.obs.columns:
        adata.obs["celltype"] = adata.obs["cell_type"].values
    if "spot_barcode" in adata.obs.columns:
        adata.obs["barcode"] = adata.obs["spot_barcode"].values

    cancer_mask = adata.obs["celltype"].str.contains("Cancer", case=False)

    available_krt = [g for g in BASAL_KERATINS if g in adata.var_names]
    if not available_krt:
        raise ValueError(f"No basal keratins found: {BASAL_KERATINS}")
    X = adata[:, available_krt].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    krt_sum = np.asarray(X).sum(axis=1)
    adata.obs["D538G"] = (krt_sum > 0) & cancer_mask
    log.info(
        "D538G via RNA basal cytokeratin: %d / %d cancer (keratins: %s)",
        adata.obs["D538G"].sum(),
        cancer_mask.sum(),
        available_krt,
    )
    return adata


def hurdle_test(sender_df: pd.DataFrame, mut_mask: np.ndarray) -> list[dict]:
    """Run hurdle model on each pathway in sender_df."""
    rows = []
    for pathway in sender_df.columns:
        mut_vals = sender_df.loc[mut_mask, pathway].values.astype(float)
        wt_vals = sender_df.loc[~mut_mask, pathway].values.astype(float)

        # --- Part 1: Participation (Fisher exact) ---
        mut_active = int((mut_vals > 0).sum())
        mut_zero = int((mut_vals == 0).sum())
        wt_active = int((wt_vals > 0).sum())
        wt_zero = int((wt_vals == 0).sum())

        table = [[mut_active, mut_zero], [wt_active, wt_zero]]
        odds_ratio, fisher_p = fisher_exact(table, alternative="two-sided")

        mut_pct = mut_active / len(mut_vals) * 100 if len(mut_vals) > 0 else 0
        wt_pct = wt_active / len(wt_vals) * 100 if len(wt_vals) > 0 else 0

        # --- Part 2: Magnitude among active senders (MW on nonzero) ---
        mut_nonzero = mut_vals[mut_vals > 0]
        wt_nonzero = wt_vals[wt_vals > 0]

        if len(mut_nonzero) >= 2 and len(wt_nonzero) >= 2:
            mw_stat, mw_p = mannwhitneyu(mut_nonzero, wt_nonzero, alternative="two-sided")
            sender_fc = float(mut_nonzero.mean() / wt_nonzero.mean()) if wt_nonzero.mean() > 0 else float("inf")
        else:
            mw_stat, mw_p, sender_fc = float("nan"), float("nan"), float("nan")

        # --- Also compute standard MW on ALL cells for comparison ---
        if len(mut_vals) >= 2 and len(wt_vals) >= 2:
            std_stat, std_p = mannwhitneyu(mut_vals, wt_vals, alternative="two-sided")
            std_fc = float(mut_vals.mean() / wt_vals.mean()) if wt_vals.mean() > 0 else float("inf")
        else:
            std_stat, std_p, std_fc = float("nan"), float("nan"), float("nan")

        rows.append(
            {
                "pathway": pathway,
                "n_mut": len(mut_vals),
                "n_wt": len(wt_vals),
                "mut_active_pct": round(mut_pct, 1),
                "wt_active_pct": round(wt_pct, 1),
                "fisher_OR": round(float(odds_ratio), 2),
                "fisher_p": float(fisher_p),
                "n_mut_nonzero": len(mut_nonzero),
                "n_wt_nonzero": len(wt_nonzero),
                "sender_fc_nonzero": round(sender_fc, 2) if not np.isnan(sender_fc) else None,
                "mw_p_nonzero": float(mw_p) if not np.isnan(mw_p) else None,
                "std_fc_all": round(std_fc, 2) if not np.isnan(std_fc) else None,
                "std_p_all": float(std_p) if not np.isnan(std_p) else None,
            }
        )

    df = pd.DataFrame(rows)
    if len(df) > 0:
        valid_fisher = df["fisher_p"].dropna()
        if len(valid_fisher) > 0:
            _, fisher_fdr = fdrcorrection(valid_fisher)
            df.loc[valid_fisher.index, "fisher_fdr"] = fisher_fdr

        valid_mw = df["mw_p_nonzero"].dropna()
        if len(valid_mw) > 0:
            _, mw_fdr = fdrcorrection(valid_mw)
            df.loc[valid_mw.index, "mw_fdr_nonzero"] = mw_fdr

    return df.sort_values("mw_p_nonzero").to_dict("records")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = load_and_classify(SC_H5AD)

    cancer = adata[adata.obs["celltype"].str.contains("Cancer", case=False)].copy()

    total_umi = np.array(cancer.X.sum(axis=1)).ravel()
    umi_mask = total_umi >= MIN_UMI
    log.info(
        "UMI filter (>=%d): %d pass, %d filtered of %d cancer", MIN_UMI, umi_mask.sum(), (~umi_mask).sum(), len(cancer)
    )
    cancer = cancer[umi_mask].copy()

    # Load COMMOT sender scores from per-cell CSV (exported by analyze_p4s2_commot.py).
    # The CSV has one row per cancer cell with barcode as index (duplicates exist
    # since multiple cells share a spot). It was exported from ALL cancer cells
    # before UMI filtering, so we need to re-apply the same filter.
    sender_csv = OUTPUT_DIR / f"{SAMPLE}_commot_mdk_percell.csv"
    if not sender_csv.exists():
        raise FileNotFoundError(
            f"Missing COMMOT per-cell CSV: {sender_csv}\n" "Run analyze_p4s2_commot.py first to generate sender scores."
        )
    sender_all = pd.read_csv(sender_csv)
    log.info("Loaded %d per-cell sender scores from %s", len(sender_all), sender_csv.name)

    # The CSV was exported from all cancer cells in the same h5ad. Rebuild the
    # same cancer subset and apply UMI filter to get aligned indices.
    all_cancer = adata[adata.obs["celltype"].str.contains("Cancer", case=False)].copy()
    all_cancer_umi = np.array(all_cancer.X.sum(axis=1)).ravel()
    all_cancer_umi_mask = all_cancer_umi >= MIN_UMI

    if len(sender_all) != len(all_cancer):
        log.warning(
            "Sender CSV rows (%d) != cancer cells in h5ad (%d); using barcode join", len(sender_all), len(all_cancer)
        )
        # Fallback: broadcast spot-level scores to cells via barcode
        sender_all = sender_all.rename(columns={"barcode": "spot_barcode"})
        bc_col = "barcode" if "barcode" in cancer.obs.columns else "spot_barcode"
        cancer_bc = cancer.obs[bc_col].astype(str).values
        sender_all_indexed = (
            sender_all.set_index("spot_barcode") if "spot_barcode" in sender_all.columns else sender_all
        )
        pathway_cols = [c for c in sender_all_indexed.columns if c.startswith("s-")]
        sender_df = sender_all_indexed[pathway_cols].reindex(cancer_bc).fillna(0.0)
        sender_df.index = range(len(sender_df))
    else:
        # Positional alignment: same order as h5ad cancer subset
        pathway_cols = [c for c in sender_all.columns if c.startswith("s-")]
        sender_filtered = sender_all.loc[all_cancer_umi_mask, pathway_cols].reset_index(drop=True)
        sender_df = sender_filtered

    mut_mask = cancer.obs["D538G"].values.astype(bool)

    log.info("Testing %d pathways: %d D538G, %d WT", len(sender_df.columns), mut_mask.sum(), (~mut_mask).sum())

    results = hurdle_test(sender_df, mut_mask)

    mdk_results = {r["pathway"]: r for r in results if r["pathway"] in MDK_PATHWAYS}
    total_row = next((r for r in results if r["pathway"] == "s-total-total"), None)

    manuscript_numbers = {
        "sample": SAMPLE,
        "umi_filter": MIN_UMI,
        "n_cancer_umi_filtered": int(len(cancer)),
        "n_d538g": int(mut_mask.sum()),
        "n_wt": int((~mut_mask).sum()),
        "zero_inflation_pct": round(float((sender_df.values == 0).mean()) * 100, 1),
        "mdk_pathways": mdk_results,
        "total_sender_participation": {
            "fisher_OR": total_row["fisher_OR"] if total_row else None,
            "fisher_p": total_row["fisher_p"] if total_row else None,
            "mut_active_pct": total_row["mut_active_pct"] if total_row else None,
            "wt_active_pct": total_row["wt_active_pct"] if total_row else None,
        },
    }

    out_json = OUTPUT_DIR / "hurdle_results.json"
    with open(out_json, "w") as f:
        json.dump(manuscript_numbers, f, indent=2, default=str)
    log.info("Wrote %s", out_json)

    out_csv = OUTPUT_DIR / "hurdle_results_summary.csv"
    pd.DataFrame(results).to_csv(out_csv, index=False)
    log.info("Wrote %s", out_csv)

    print("\n" + "=" * 70)
    print("HURDLE MODEL RESULTS")
    print("=" * 70)
    print(f"Cancer cells (UMI>={MIN_UMI}): {len(cancer)} ({mut_mask.sum()} D538G, {(~mut_mask).sum()} WT)")
    if total_row:
        print(f"\nTotal sender participation:")
        print(f"  D538G: {total_row['mut_active_pct']}% active | WT: {total_row['wt_active_pct']}% active")
        print(f"  Fisher exact OR={total_row['fisher_OR']}, p={total_row['fisher_p']:.2e}")
    for pw in MDK_PATHWAYS:
        r = mdk_results.get(pw)
        if r:
            print(f"\n{pw}:")
            print(f"  Standard MW (all cells): FC={r['std_fc_all']}, p={r['std_p_all']:.2e}" if r["std_p_all"] else "")
            print(
                f"  Hurdle MW (senders only): FC={r['sender_fc_nonzero']}, p={r['mw_p_nonzero']:.2e}"
                if r["mw_p_nonzero"]
                else f"  Hurdle MW: insufficient nonzero (n_mut={r['n_mut_nonzero']}, n_wt={r['n_wt_nonzero']})"
            )
    print("=" * 70)


if __name__ == "__main__":
    main()
