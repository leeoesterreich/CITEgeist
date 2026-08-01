#!/usr/bin/env python
"""
Task 6: Bivariate Moran's I for competitor GEX validation.

Computes bivariate spatial autocorrelation for candidate genes on
pseudo-single-cell data from each method/sample.

For MDK specifically, we compute the SAME metric as CITEgeist:
  Moran_BV(MDK, secretory_composite, w)
where secretory_composite = mean of z-scored SECRETORY_GENES.
This ensures a like-for-like comparison on the S12 bar chart.

For other candidate genes, standard univariate Moran's I is computed.

Samples where MDK has <MDK_NONZERO_THRESHOLD nonzero cells are skipped
(saved as NaN), matching CITEgeist's sample-inclusion rule
(mdk_analysis/scripts/spatial_morans.py:293). This ensures equal n.

Usage:
    python -u run_bivariate_morans.py --method cell2location --sample HCC22-088-P1-S1
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import false_discovery_control, zscore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import (
    CANDIDATE_GENES,
    MDK_NONZERO_THRESHOLD,
    MORANS_K_NEIGH,
    MORANS_N_PERMS,
    OUTPUT_ROOT,
    SECRETORY_GENES,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def compute_morans(gene_expr: np.ndarray, second_var, w, permutations: int = 999):
    """Compute bivariate Moran's I between two variables.

    Parameters
    ----------
    gene_expr : array
        First variable (e.g. MDK expression).
    second_var : array or None
        Second variable (e.g. secretory composite). If None, computes
        univariate self-autocorrelation (gene_expr vs gene_expr).
    w : libpysal spatial weights
    permutations : int
        Number of random permutations for pseudo p-value.

    Returns
    -------
    I, p_sim : float, float
    """
    from esda.moran import Moran_BV

    if second_var is None:
        second_var = gene_expr
    result = Moran_BV(gene_expr, second_var, w, permutations=permutations)
    return result.I, result.p_sim


def build_secretory_composite(adata):
    """Construct per-cell secretory composite score: mean of z-scored
    SECRETORY_GENES present in adata.var_names.

    Returns None if no secretory genes are found.
    Matches CITEgeist's construction in mdk_analysis/scripts/spatial_morans.py:301-311.
    """
    import scipy.sparse

    present = [g for g in SECRETORY_GENES if g in adata.var_names]
    if not present:
        return None

    sec_x = adata[:, present].X
    if scipy.sparse.issparse(sec_x):
        sec_mat = sec_x.toarray().astype(float)
    else:
        sec_mat = np.asarray(sec_x, dtype=float)

    sec_z = zscore(sec_mat, axis=0, nan_policy="omit")
    sec_z = np.nan_to_num(sec_z, 0)
    sec_score = sec_z.mean(axis=1)
    return sec_score


def main():
    parser = argparse.ArgumentParser(description="Bivariate Moran's I analysis")
    parser.add_argument("--method", required=True, choices=["cell2location", "tangram"])
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    method = args.method
    sample = args.sample

    # Paths
    sc_path = OUTPUT_ROOT / method / sample / f"{sample}_single_cell.h5ad"
    out_dir = OUTPUT_ROOT / method / "morans"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{sample}_morans.csv"

    if not sc_path.exists():
        logger.error("Single-cell h5ad not found: %s", sc_path)
        sys.exit(1)

    logger.info("=== Moran's I: method=%s sample=%s ===", method, sample)

    # Load pseudo-single-cell data
    adata = sc.read_h5ad(sc_path)
    adata.var_names_make_unique()
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Filter out low-UMI cells (empty/near-empty pseudo-cells confuse spatial stats)
    import scipy.sparse

    total_umi = np.array(adata.X.sum(axis=1)).flatten()
    umi_mask = total_umi > 20
    n_removed = (~umi_mask).sum()
    if n_removed > 0:
        logger.info("Filtering %d/%d cells with <=20 UMI", n_removed, adata.n_obs)
        adata = adata[umi_mask].copy()
        logger.info("After filter: %d cells", adata.n_obs)

    # Get expression matrix (dense)

    if scipy.sparse.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)

    # Build spatial weights
    from libpysal.weights import KNN

    coords = np.array(adata.obsm["spatial"])
    logger.info("Building KNN spatial weights (k=%d)...", MORANS_K_NEIGH)
    w = KNN.from_array(coords, k=MORANS_K_NEIGH)
    w.transform = "R"

    gene_names = list(adata.var_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Pre-build secretory composite for MDK bivariate analysis
    sec_score = build_secretory_composite(adata)
    if sec_score is not None:
        n_sec = sum(1 for g in SECRETORY_GENES if g in adata.var_names)
        logger.info("Built secretory composite from %d/%d genes", n_sec, len(SECRETORY_GENES))
    else:
        logger.warning("No secretory genes found — MDK will use self-autocorrelation fallback")

    results = []

    # --- Candidate genes first ---
    logger.info("Computing Moran's I for %d candidate genes...", len(CANDIDATE_GENES))
    for gene in CANDIDATE_GENES:
        if gene not in gene_to_idx:
            logger.warning("Candidate gene %s not found in var_names, skipping", gene)
            results.append(
                {
                    "gene": gene,
                    "I": np.nan,
                    "p_value": np.nan,
                    "q_value": np.nan,
                    "method": method,
                    "sample": sample,
                    "is_candidate": True,
                    "note": "gene not in var_names",
                }
            )
            continue

        idx = gene_to_idx[gene]
        expr = X[:, idx].astype(np.float64)

        # Bug 2 fix: MDK nonzero threshold — skip samples where MDK has too
        # few nonzero cells, matching CITEgeist's inclusion rule.
        if gene == "MDK":
            mdk_nonzero = int(np.sum(expr > 0))
            if mdk_nonzero < MDK_NONZERO_THRESHOLD:
                logger.info(
                    "  MDK: SKIP sample %s — only %d nonzero cells (threshold=%d)",
                    sample,
                    mdk_nonzero,
                    MDK_NONZERO_THRESHOLD,
                )
                results.append(
                    {
                        "gene": gene,
                        "I": np.nan,
                        "p_value": np.nan,
                        "q_value": np.nan,
                        "method": method,
                        "sample": sample,
                        "is_candidate": True,
                        "note": f"MDK below detection: {mdk_nonzero} nonzero < {MDK_NONZERO_THRESHOLD}",
                    }
                )
                continue

        # For MDK: bivariate Moran's I against secretory composite
        # (matches CITEgeist's metric). For other genes: self-autocorrelation.
        if gene == "MDK" and sec_score is not None:
            second = sec_score
            logger.info("  MDK: computing bivariate Moran's I (MDK vs secretory composite)")
        else:
            second = None

        I_val, p_val = compute_morans(expr, second, w, permutations=MORANS_N_PERMS)
        logger.info("  %s: I=%.4f, p=%.4f", gene, I_val, p_val)
        results.append(
            {
                "gene": gene,
                "I": I_val,
                "p_value": p_val,
                "q_value": np.nan,  # filled after FDR on full transcriptome
                "method": method,
                "sample": sample,
                "is_candidate": True,
            }
        )

    # --- Full transcriptome ---
    logger.info("Computing Moran's I for full transcriptome (%d genes)...", len(gene_names))
    transcriptome_results = []
    for i, gene in enumerate(gene_names):
        if i > 0 and i % 500 == 0:
            logger.info("  Progress: %d / %d genes", i, len(gene_names))
        expr = X[:, i].astype(np.float64)
        # Skip zero-variance genes
        if np.std(expr) < 1e-10:
            transcriptome_results.append(
                {
                    "gene": gene,
                    "I": np.nan,
                    "p_value": np.nan,
                    "method": method,
                    "sample": sample,
                    "is_candidate": gene in CANDIDATE_GENES,
                }
            )
            continue
        # Transcriptome-wide uses self-autocorrelation (univariate)
        I_val, p_val = compute_morans(expr, None, w, permutations=MORANS_N_PERMS)
        transcriptome_results.append(
            {
                "gene": gene,
                "I": I_val,
                "p_value": p_val,
                "method": method,
                "sample": sample,
                "is_candidate": gene in CANDIDATE_GENES,
            }
        )

    # BH FDR correction on non-NaN transcriptome p-values
    tx_df = pd.DataFrame(transcriptome_results)
    valid_mask = tx_df["p_value"].notna()
    if valid_mask.sum() > 0:
        pvals = tx_df.loc[valid_mask, "p_value"].values
        qvals = false_discovery_control(pvals, method="bh")
        tx_df.loc[valid_mask, "q_value"] = qvals
    else:
        tx_df["q_value"] = np.nan

    # Merge: candidate results use the transcriptome q-values where available
    candidate_df = pd.DataFrame(results)
    for i, row in candidate_df.iterrows():
        gene = row["gene"]
        match = tx_df.loc[tx_df["gene"] == gene]
        if len(match) > 0:
            candidate_df.at[i, "q_value"] = match.iloc[0]["q_value"]

    # Remove candidate genes from transcriptome df to avoid duplicates
    candidate_set = set(CANDIDATE_GENES)
    tx_df_no_cand = tx_df[~tx_df["gene"].isin(candidate_set)]

    # Combine: candidates first, then rest of transcriptome
    final_df = pd.concat([candidate_df, tx_df_no_cand], ignore_index=True)
    final_df.to_csv(out_path, index=False)

    # Summary
    n_sig = (final_df["p_value"] <= 0.001).sum()
    n_fdr = (final_df["q_value"].dropna() <= 0.05).sum()
    logger.info("Saved %d genes to %s", len(final_df), out_path)
    logger.info("  Significant (p<=0.001): %d, FDR<0.05: %d", n_sig, n_fdr)

    # Print candidate summary
    for _, row in candidate_df.iterrows():
        sig = "***" if row["p_value"] <= 0.001 else ""
        logger.info(
            "  CANDIDATE %s: I=%.4f, p=%.4f, q=%.4f %s",
            row["gene"],
            row["I"] if not np.isnan(row["I"]) else 0.0,
            row["p_value"] if not np.isnan(row["p_value"]) else 1.0,
            row["q_value"] if not np.isnan(row["q_value"]) else 1.0,
            sig,
        )


if __name__ == "__main__":
    main()
