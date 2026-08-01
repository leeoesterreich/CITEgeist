#!/usr/bin/env python
"""
Reproducible P4-S2 D538G mutant identification + COMMOT signaling analysis.

Consolidates vignette 2 analysis as a headless script with all results
persisted to mdk_analysis/outputs/commot_p4s2/.

Loads the single-cell h5ad from morphology_assignment_v3 (current CITEgeist
pipeline output), classifies D538G via basal cytokeratins in the cancer cell
layer, runs COMMOT, computes differential signaling, DEG, Enrichr, and ESR1
signature validation.

Usage:
    python -u mdk_analysis/scripts/analyze_p4s2_commot.py
"""

import json
import logging
import os
import sys
from pathlib import Path

import anndata as ad

try:
    import commot as ct
except ImportError:
    print("ERROR: commot not available. Run in COMMOT conda environment.")
    sys.exit(1)

import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, spearmanr, ttest_ind
from statsmodels.stats.multitest import fdrcorrection

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
SAMPLE = "HCC22-088-P4-S2_1i_rep"
SC_H5AD = PROJECT / "output" / "morphology_assignment_v3" / SAMPLE / f"{SAMPLE}_single_cell.h5ad"
PROTEIN_GATES_CSV = PROJECT / "output" / "sace_protein_12patient" / SAMPLE / f"{SAMPLE}_protein_gates.csv"
OUTPUT_DIR = PROJECT / "mdk_analysis" / "outputs" / "commot_p4s2"

BASAL_KERATINS = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT16", "KRT17"]

HPA_CELL_MAPPING = {
    "Cancer_Cells": "Breast glandular cells",
    "CD8_T_Cells": "T-cells",
    "CD4_T_Cells": "T-cells",
    "Macrophages": "Macrophages",
    "B_Cells": "B-cells",
    "Endothelial": "Endothelial cells",
    "Fibroblasts": "Fibroblasts",
}


def load_single_cell_adata() -> ad.AnnData:
    """Load single-cell AnnData from morphology_assignment_v3."""
    log.info("Loading single-cell h5ad from %s", SC_H5AD)
    adata = sc.read_h5ad(str(SC_H5AD))
    adata.var_names_make_unique()

    # Standardize column names for COMMOT compatibility
    # (commot_all_samples.py uses 'celltype' and 'barcode')
    if "cell_type" in adata.obs.columns:
        adata.obs["celltype"] = adata.obs["cell_type"].values
    if "spot_barcode" in adata.obs.columns:
        adata.obs["barcode"] = adata.obs["spot_barcode"].values

    log.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)
    log.info("Cell types: %s", adata.obs["celltype"].value_counts().to_dict())
    return adata


def classify_d538g(adata: ad.AnnData) -> ad.AnnData:
    """Classify cancer cells as D538G-mutant via Module 3.5 KRT5 protein gate.

    Uses SACE protein allocation + GMM gating (func_KRT5_Epithelial_gate)
    from the M3.5 pipeline. Falls back to raw GEX KRT > 0 if gates unavailable.
    """
    cancer_mask = adata.obs["celltype"].str.contains("Cancer", case=False)

    if PROTEIN_GATES_CSV.exists():
        gates = pd.read_csv(PROTEIN_GATES_CSV, index_col=0)
        gates.index = gates.index.astype(str)
        gate_col = "func_KRT5_Epithelial_gate"
        if gate_col in gates.columns:
            krt5_gate = gates[gate_col].reindex(adata.obs.index.astype(str)).fillna(0).astype(bool)
            adata.obs["D538G Mutation"] = krt5_gate.values & cancer_mask.values
            n_mut = adata.obs["D538G Mutation"].sum()
            log.info(
                "D538G classification (M3.5 KRT5 protein gate): %d mutant / %d cancer cells", n_mut, cancer_mask.sum()
            )
            return adata
        else:
            log.warning("KRT5 gate column not found in %s, falling back to GEX", PROTEIN_GATES_CSV)
    else:
        log.warning("Protein gates not found at %s, falling back to GEX", PROTEIN_GATES_CSV)

    # Fallback: raw GEX basal cytokeratin > 0
    available_krt = [g for g in BASAL_KERATINS if g in adata.var_names]
    if not available_krt:
        raise ValueError(f"No basal keratins found. Expected: {BASAL_KERATINS}")
    X = adata[:, available_krt].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    krt_sum = np.asarray(X).sum(axis=1)
    adata.obs["D538G Mutation"] = (krt_sum > 0) & cancer_mask
    n_mut = adata.obs["D538G Mutation"].sum()
    log.info(
        "D538G classification (GEX fallback): %d mutant / %d cancer cells (keratins: %s)",
        n_mut,
        cancer_mask.sum(),
        available_krt,
    )
    return adata


def run_commot(adata: ad.AnnData) -> ad.AnnData:
    """Run COMMOT spatial communication on single-cell AnnData."""
    df_ligrec = ct.pp.ligand_receptor_database(database="CellChat", species="human")
    df_filtered = ct.pp.filter_lr_database(df_ligrec, adata, min_cell_pct=0.01)
    log.info("Filtered L-R pairs: %d", len(df_filtered))
    for _, row in df_filtered.iterrows():
        log.info("  %s -> %s (%s)", row.iloc[0], row.iloc[1], row.iloc[3])

    ct.tl.spatial_communication(
        adata,
        database_name="cellchat",
        df_ligrec=df_filtered,
        dis_thr=500,
        heteromeric=True,
        pathway_sum=True,
    )
    log.info("COMMOT complete")
    return adata


def compute_commot_differential(adata: ad.AnnData) -> pd.DataFrame:
    """Compute differential signaling between D538G mutant and WT cancer cells."""
    cancer = adata[adata.obs["celltype"].str.contains("Cancer", case=False)].copy()

    # Filter out low-UMI cells (empty/dead cells with insufficient RNA)
    MIN_UMI = 20
    total_umi = np.array(cancer.X.sum(axis=1)).ravel()
    umi_mask = total_umi >= MIN_UMI
    n_filtered = (~umi_mask).sum()
    log.info(
        "UMI filter (>=%d): %d cells pass, %d filtered out of %d cancer cells",
        MIN_UMI,
        umi_mask.sum(),
        n_filtered,
        len(cancer),
    )
    cancer = cancer[umi_mask].copy()

    sender_key = "commot-cellchat-sum-sender"
    if sender_key not in cancer.obsm:
        log.error("COMMOT sender signals not found in obsm. Keys: %s", list(cancer.obsm.keys()))
        raise RuntimeError(f"Missing {sender_key} — COMMOT may have failed silently")
    sender = cancer.obsm[sender_key]
    mut_mask = cancer.obs["D538G Mutation"].values.astype(bool)

    rows = []
    for pathway in sender.columns:
        mut_vals = sender.loc[mut_mask, pathway].dropna()
        wt_vals = sender.loc[~mut_mask, pathway].dropna()
        if len(mut_vals) == 0 or len(wt_vals) == 0:
            continue
        stat, pval = mannwhitneyu(mut_vals, wt_vals, alternative="two-sided")
        mut_mean = float(mut_vals.mean())
        wt_mean = float(wt_vals.mean())
        fc = mut_mean / wt_mean if wt_mean > 0 else float("inf")
        rows.append(
            {
                "pathway": pathway,
                "mut_mean": mut_mean,
                "wt_mean": wt_mean,
                "fold_change": fc,
                "mean_diff": mut_mean - wt_mean,
                "p_value": pval,
                "statistic": stat,
                "n_mut": len(mut_vals),
                "n_wt": len(wt_vals),
            }
        )

    df = pd.DataFrame(rows)
    if len(df) > 0:
        df = df.sort_values("p_value")
        _, fdr = fdrcorrection(df["p_value"])
        df["fdr"] = fdr
    else:
        df["fdr"] = pd.Series(dtype=float)
    log.info("COMMOT differential:\n%s", df.to_string(index=False))
    return df


def run_deg_enrichr(adata: ad.AnnData) -> dict:
    """Run Wilcoxon DEG on cancer cells, then Enrichr on upregulated genes."""
    cancer = adata[adata.obs["celltype"].str.contains("Cancer", case=False)].copy()
    cancer.obs["D538G Mutation"] = cancer.obs["D538G Mutation"].astype("category")
    sc.pp.log1p(cancer)
    sc.tl.rank_genes_groups(cancer, "D538G Mutation", method="wilcoxon")

    results = cancer.uns["rank_genes_groups"]
    groups = list(results["names"].dtype.names)
    group_key = [g for g in groups if str(g) == "True" or g is True][0]
    df = pd.DataFrame(
        {
            "gene": results["names"][group_key],
            "score": results["scores"][group_key],
            "logfc": results["logfoldchanges"][group_key],
            "pval": results["pvals"][group_key],
            "pval_adj": results["pvals_adj"][group_key],
        }
    )
    up = df[(df["score"] > 0) & (df["pval_adj"] < 0.05) & (df["logfc"] > 0)]
    gene_list = up["gene"].tolist()
    log.info("D538G upregulated genes: %d (of %d tested)", len(gene_list), len(df))

    enrichr_results = {}
    for gs in ["MSigDB_Hallmark_2020", "KEGG_2021_Human"]:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gs,
                organism="Human",
                outdir=str(OUTPUT_DIR / "enrichr"),
            )
            enrichr_results[gs] = (
                enr.results[["Term", "Adjusted P-value", "Overlap", "Genes"]].head(15).to_dict("records")
            )
        except Exception as e:
            log.warning("Enrichr failed for %s: %s", gs, e)

    return {
        "deg_table": up.head(30).to_dict("records"),
        "n_upregulated": len(gene_list),
        "gene_list": gene_list,
        "enrichr": enrichr_results,
    }


def validate_esr1_signature(adata: ad.AnnData) -> dict:
    """Validate D538G classification with published ESR1 mutation signature."""
    cancer = adata[adata.obs["celltype"].str.contains("Cancer", case=False)].copy()
    sc.pp.log1p(cancer)

    excel_path = PROJECT / "examples" / "media-3.xlsx"
    if not excel_path.exists():
        excel_path = PROJECT / "CITEgeist" / "examples" / "media-3.xlsx"
    if not excel_path.exists():
        log.warning("ESR1 signature file not found")
        return {"status": "skipped", "reason": "media-3.xlsx not found"}

    df = pd.read_excel(excel_path, sheet_name=0)
    up_col = [c for c in df.columns if "Up genes" in c][0]
    down_col = [c for c in df.columns if "Down genes" in c][0]
    up_genes = [g for g in df[up_col].dropna() if g in cancer.var_names]
    down_genes = [g for g in df[down_col].dropna() if g in cancer.var_names]

    if not up_genes or not down_genes:
        log.warning("No ESR1 signature genes found in adata var_names")
        return {"status": "skipped", "reason": "no overlapping signature genes"}

    X_up = cancer[:, up_genes].X
    X_down = cancer[:, down_genes].X
    if hasattr(X_up, "toarray"):
        X_up, X_down = X_up.toarray(), X_down.toarray()
    score = np.asarray(X_up).mean(axis=1) - np.asarray(X_down).mean(axis=1)

    mut_mask = cancer.obs["D538G Mutation"].values.astype(bool)
    t_stat, p_val = ttest_ind(score[mut_mask], score[~mut_mask], equal_var=False)

    result = {
        "t_statistic": float(t_stat),
        "p_value": float(p_val),
        "mut_mean": float(score[mut_mask].mean()),
        "wt_mean": float(score[~mut_mask].mean()),
        "n_up_genes": len(up_genes),
        "n_down_genes": len(down_genes),
        "n_mut": int(mut_mask.sum()),
        "n_wt": int((~mut_mask).sum()),
    }
    log.info("ESR1 signature: t=%.3f, p=%.2e", t_stat, p_val)
    return result


def compute_hpa_correlation(adata: ad.AnnData) -> dict:
    """Correlate COMMOT sender signals with HPA cell-type expression."""
    hpa_path = PROJECT / "mdk_analysis" / "data" / "rna_single_cell_type.tsv"
    if not hpa_path.exists():
        log.warning("HPA data not found: %s", hpa_path)
        return {"status": "skipped", "reason": "rna_single_cell_type.tsv not found"}

    hpa_df = pd.read_csv(hpa_path, sep="\t")
    sender_df = adata.obsm.get("commot-cellchat-sum-sender")
    if sender_df is None:
        return {"status": "skipped", "reason": "no sender signals"}

    celltypes = adata.obs["celltype"].values

    global_sender, global_hpa = [], []
    for pathway in sender_df.columns:
        if not pathway.startswith("s-"):
            continue
        path = pathway[2:]  # strip "s-" prefix
        if path == "total-total":
            continue
        ligand = path.split("-")[0]
        for ct_name, hpa_name in HPA_CELL_MAPPING.items():
            mask = celltypes == ct_name
            if mask.sum() == 0:
                continue
            sender_mean = float(sender_df.loc[mask, pathway].mean())
            hpa_match = hpa_df[(hpa_df["Gene name"] == ligand) & (hpa_df["Cell type"].str.lower() == hpa_name.lower())]
            if hpa_match.empty:
                continue
            global_sender.append(sender_mean)
            global_hpa.append(float(hpa_match["nTPM"].mean()))

    if len(global_sender) < 2:
        return {"status": "insufficient_data", "n_points": len(global_sender)}

    rho, p = spearmanr(global_sender, global_hpa)
    result = {
        "spearman_rho": round(float(rho), 2),
        "p_value": float(p),
        "n_data_points": len(global_sender),
    }
    log.info("HPA correlation: rho=%.2f, p=%.2e (n=%d)", rho, p, len(global_sender))
    return result


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = load_single_cell_adata()
    adata = classify_d538g(adata)
    adata = run_commot(adata)

    commot_df = compute_commot_differential(adata)
    commot_df.to_csv(OUTPUT_DIR / "commot_differential.csv", index=False)

    deg_results = run_deg_enrichr(adata)
    esr1_results = validate_esr1_signature(adata)
    hpa_results = compute_hpa_correlation(adata)

    # Save cancer-cell sender scores
    cancer = adata[adata.obs["celltype"].str.contains("Cancer", case=False)]
    sender_out = cancer.obsm["commot-cellchat-sum-sender"].copy()
    if "barcode" in cancer.obs.columns:
        sender_out.index = cancer.obs["barcode"].values
    sender_out.index.name = "barcode"
    sender_out.to_csv(OUTPUT_DIR / f"{SAMPLE}_commot_mdk_percell.csv")

    # Build manuscript numbers
    mdk_row = commot_df[commot_df["pathway"] == "s-MDK-SDC4"]
    if len(mdk_row) == 0:
        log.error("MDK-SDC4 pathway not found in results!")
        mdk_numbers = {"error": "pathway not found"}
    else:
        mdk_row = mdk_row.iloc[0]
        mdk_numbers = {
            "fold_change": round(float(mdk_row["fold_change"]), 2),
            "fdr": float(mdk_row["fdr"]),
            "mut_mean": round(float(mdk_row["mut_mean"]), 3),
            "wt_mean": round(float(mdk_row["wt_mean"]), 3),
        }

    manuscript_numbers = {
        "sample": SAMPLE,
        "classification": "m35_krt5_protein_gate" if PROTEIN_GATES_CSV.exists() else "basal_cytokeratin_gt_zero",
        "n_cells_total": int(adata.n_obs),
        "n_cancer_cells": int((adata.obs["celltype"].str.contains("Cancer", case=False)).sum()),
        "n_d538g_mutant": int(adata.obs["D538G Mutation"].sum()),
        "n_d538g_wt": int(
            ((adata.obs["celltype"].str.contains("Cancer", case=False)) & ~adata.obs["D538G Mutation"]).sum()
        ),
        "commot": {
            "mdk_sdc4": mdk_numbers,
            "n_significant_pathways": int((commot_df["fdr"] < 0.05).sum()),
            "lr_pairs_tested": commot_df["pathway"].tolist(),
        },
        "deg": {"n_upregulated": deg_results["n_upregulated"]},
        "esr1_signature": esr1_results,
        "hpa_correlation": hpa_results,
    }
    with open(OUTPUT_DIR / "manuscript_numbers.json", "w") as f:
        json.dump(manuscript_numbers, f, indent=2)

    # Save DEG gene list
    pd.DataFrame(deg_results["deg_table"]).to_csv(OUTPUT_DIR / "d538g_upregulated_genes.csv", index=False)

    log.info("All results saved to %s", OUTPUT_DIR)
    print("\n" + "=" * 60)
    print("MANUSCRIPT NUMBERS")
    print("=" * 60)
    print(f"MDK-SDC4 fold change: {mdk_numbers.get('fold_change', 'N/A')}x")
    print(f"MDK-SDC4 FDR: {mdk_numbers.get('fdr', 'N/A')}")
    print(f"Significant pathways (FDR<0.05): {(commot_df['fdr'] < 0.05).sum()}")
    print(f"D538G upregulated genes: {deg_results['n_upregulated']}")
    print(f"ESR1 signature p: {esr1_results.get('p_value', 'N/A')}")
    print(f"HPA Spearman rho: {hpa_results.get('spearman_rho', 'N/A')}")
    print("=" * 60)


if __name__ == "__main__":
    main()
