#!/usr/bin/env python
"""
MDK v3 Single-Cell Reanalysis

Re-runs Module 4 MDK secretory program discovery on v3 morphology pipeline
per-cell GEX outputs. Compares findings to original spot-level analysis.

Input: output/morphology_assignment_v3/HCC22-088-P4-S2_1i_rep/_single_cell.h5ad
Output: output/mdk_v3_reanalysis/ (CSVs + comparison report)
"""

import logging
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.stats import mannwhitneyu, pearsonr
from sklearn.decomposition import NMF

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================
SAMPLE_NAME = "HCC22-088-P4-S2_1i_rep"
V3_DIR = Path("output/morphology_assignment_v3") / SAMPLE_NAME
SC_H5AD = V3_DIR / f"{SAMPLE_NAME}_single_cell.h5ad"
OUTPUT_DIR = Path("output/mdk_v3_reanalysis")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Both QP subtypes are pooled as a single Epithelial population for analysis.
# Within-Epithelial stratification uses the KRT5 functional gate (Module 3.5)
# when available; falls back to QP cell type label for backward compatibility.
CANCER_TYPES = ['Cancer_Basal', 'Cancer_Luminal']
MIN_DETECTION_FRAC = 0.05  # Keep genes detected in >=5% of cancer cells
NMF_SEEDS = [42, 123, 456, 789, 1024]
NMF_K_VALUES = [3, 4, 5]
PRIMARY_K = 5
PRIMARY_SEED = 42

# Known gene sets for flagging (NOT force-included in NMF)
SECRETORY_GENES = ["HSP90B1", "XBP1", "CALR", "CANX", "PDIA4", "PDIA6", "ERO1A", "SEC61A1"]
MDK_RECEPTORS = ["SDC4", "NCL", "PTPRZ1", "LRP1", "GPC1"]
ER_RELATED = ["ESR1", "PGR", "GREB1", "TFF1", "FOXA1", "GATA3"]


def load_and_preprocess():
    """Load v3 single-cell h5ad, subset to cancer cells, normalize for NMF."""
    logger.info(f"Loading {SC_H5AD}")
    adata = sc.read_h5ad(SC_H5AD)
    logger.info(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Subset to cancer cells
    cancer_mask = adata.obs['cell_type'].isin(CANCER_TYPES)
    adata_cancer = adata[cancer_mask].copy()
    logger.info(f"Cancer cells: {adata_cancer.n_obs} "
                f"({dict(adata_cancer.obs['cell_type'].value_counts())})")

    # Store raw SACE counts before normalization
    adata_cancer.raw = adata_cancer.copy()

    # Size-normalize (NO log1p — NMF needs additive non-negative data)
    sc.pp.normalize_total(adata_cancer, target_sum=1e4)

    # Gene detection filter: keep genes in >=5% of cancer cells
    min_cells = int(MIN_DETECTION_FRAC * adata_cancer.n_obs)
    if hasattr(adata_cancer.X, 'toarray'):
        detection_counts = np.asarray((adata_cancer.X > 0).sum(axis=0)).ravel()
    else:
        detection_counts = (adata_cancer.X > 0).sum(axis=0)
        if hasattr(detection_counts, 'A1'):
            detection_counts = detection_counts.A1
        else:
            detection_counts = np.asarray(detection_counts).ravel()

    gene_mask = detection_counts >= min_cells
    n_genes_before = adata_cancer.n_vars
    adata_nmf = adata_cancer[:, gene_mask].copy()
    logger.info(f"Gene filter: {n_genes_before} -> {adata_nmf.n_vars} genes "
                f"(>={min_cells} cells, {MIN_DETECTION_FRAC:.0%} detection)")

    # Report which target genes survived
    for gene_set_name, gene_list in [
        ("MDK", ["MDK"]),
        ("Secretory", SECRETORY_GENES),
        ("MDK receptors", MDK_RECEPTORS),
        ("ER-related", ER_RELATED),
    ]:
        present = [g for g in gene_list if g in adata_nmf.var_names]
        missing = [g for g in gene_list if g not in adata_nmf.var_names]
        if missing:
            logger.warning(f"{gene_set_name}: present={present}, MISSING={missing}")
        else:
            logger.info(f"{gene_set_name}: all present ({present})")

    # Subtype annotation: use KRT5 protein gate (Module 3.5) if available,
    # otherwise fall back to QP cell type label as a proxy.
    krt5_cols = [c for c in adata_cancer.obs.columns
                 if 'KRT5' in c and c.endswith('_gate')]
    if krt5_cols:
        krt5_pos = adata_cancer.obs[krt5_cols].max(axis=1).fillna(0).astype(bool)
        adata_cancer.obs['cancer_subtype'] = np.where(krt5_pos, 'KRT5-high', 'KRT5-low')
        logger.info(f"Subtype by KRT5 gate: {krt5_pos.sum()} KRT5-high, "
                    f"{(~krt5_pos).sum()} KRT5-low")
    else:
        adata_cancer.obs['cancer_subtype'] = adata_cancer.obs['cell_type'].map({
            'Cancer_Basal': 'KRT5-high', 'Cancer_Luminal': 'KRT5-low',
        })
        logger.info("KRT5 gate not in obs; using QP cell type as subtype proxy")

    return adata_cancer, adata_nmf


def run_nmf_stability(adata_nmf):
    """Run NMF across seeds and K values. Return primary results + stability metrics."""
    X = np.asarray(adata_nmf.X.copy() if not hasattr(adata_nmf.X, 'toarray')
                   else adata_nmf.X.toarray())
    X = np.maximum(X, 0)  # safety clamp
    gene_names = list(adata_nmf.var_names)

    # --- Multi-seed runs at primary K ---
    logger.info(f"Running NMF K={PRIMARY_K} across {len(NMF_SEEDS)} seeds...")
    seed_results = {}
    for seed in NMF_SEEDS:
        model = NMF(n_components=PRIMARY_K, init='nndsvda',
                    random_state=seed, max_iter=500)
        H = model.fit_transform(X)       # (cells x programs)
        W = model.components_.T           # (genes x programs)
        recon_err = model.reconstruction_err_
        seed_results[seed] = {'W': W, 'H': H, 'recon_err': recon_err}
        logger.info(f"  Seed {seed}: recon_error={recon_err:.4f}")

    # Primary result
    W = seed_results[PRIMARY_SEED]['W']
    H = seed_results[PRIMARY_SEED]['H']

    # --- Multi-K sensitivity at primary seed ---
    logger.info(f"Running NMF K sensitivity (K={NMF_K_VALUES})...")
    k_results = {}
    for k in NMF_K_VALUES:
        model = NMF(n_components=k, init='nndsvda',
                    random_state=PRIMARY_SEED, max_iter=500)
        H_k = model.fit_transform(X)
        W_k = model.components_.T
        k_results[k] = {
            'W': W_k, 'H': H_k,
            'recon_err': model.reconstruction_err_,
        }
        # Check if MDK lands in a program
        if "MDK" in gene_names:
            mdk_idx = gene_names.index("MDK")
            mdk_loadings = W_k[mdk_idx, :]
            best_prog = int(np.argmax(mdk_loadings))
            logger.info(f"  K={k}: recon_error={model.reconstruction_err_:.4f}, "
                        f"MDK best program={best_prog} (loading={mdk_loadings[best_prog]:.4f})")
        else:
            logger.info(f"  K={k}: recon_error={model.reconstruction_err_:.4f}, MDK not in genes")

    # --- Stability metrics ---
    stability = compute_nmf_stability(seed_results, gene_names)

    return W, H, gene_names, seed_results, k_results, stability


def compute_nmf_stability(seed_results, gene_names):
    """Assess NMF stability across seeds for the MDK-containing program."""
    if "MDK" not in gene_names:
        logger.warning("MDK not in gene list — skipping stability analysis")
        return {'mdk_present': False}

    mdk_idx = gene_names.index("MDK")
    seeds = sorted(seed_results.keys())

    # For each seed, find which program has highest MDK loading
    mdk_programs = {}
    mdk_w_vectors = {}
    mdk_top50 = {}

    for seed in seeds:
        W = seed_results[seed]['W']
        mdk_loadings = W[mdk_idx, :]
        best_prog = int(np.argmax(mdk_loadings))
        mdk_programs[seed] = best_prog

        # W vector for the MDK program
        mdk_w_vectors[seed] = W[:, best_prog]

        # Top 50 genes in MDK program
        top50_idx = np.argsort(W[:, best_prog])[::-1][:50]
        mdk_top50[seed] = set(gene_names[i] for i in top50_idx)

    # Cross-seed Pearson correlation of MDK program W vectors
    ref_seed = seeds[0]
    w_correlations = {}
    for seed in seeds[1:]:
        r, _ = pearsonr(mdk_w_vectors[ref_seed], mdk_w_vectors[seed])
        w_correlations[f"{ref_seed}_vs_{seed}"] = r

    # Cross-seed Jaccard overlap of top 50 genes
    jaccard_overlaps = {}
    for seed in seeds[1:]:
        intersection = len(mdk_top50[ref_seed] & mdk_top50[seed])
        union = len(mdk_top50[ref_seed] | mdk_top50[seed])
        jaccard_overlaps[f"{ref_seed}_vs_{seed}"] = intersection / union if union > 0 else 0

    # Consensus: do all seeds agree on the same program rank?
    all_same_rank = len(set(mdk_programs.values())) == 1

    stability = {
        'mdk_present': True,
        'mdk_program_per_seed': mdk_programs,
        'all_seeds_same_rank': all_same_rank,
        'w_correlations': w_correlations,
        'jaccard_top50': jaccard_overlaps,
        'mean_w_correlation': np.mean(list(w_correlations.values())),
        'mean_jaccard': np.mean(list(jaccard_overlaps.values())),
    }

    logger.info(f"NMF stability: same_rank={all_same_rank}, "
                f"mean_W_corr={stability['mean_w_correlation']:.3f}, "
                f"mean_Jaccard={stability['mean_jaccard']:.3f}")

    return stability


def run_region_enrichment(H, adata_cancer):
    """Test each NMF program for KRT5-high vs KRT5-low enrichment within the Epithelial layer."""
    logger.info("Running region enrichment (KRT5-high vs KRT5-low)...")

    subtypes = adata_cancer.obs['cancer_subtype'].values
    basal_mask = subtypes == 'KRT5-high'
    luminal_mask = subtypes == 'KRT5-low'

    n_basal = basal_mask.sum()
    n_luminal = luminal_mask.sum()
    logger.info(f"  KRT5-high: {n_basal}, KRT5-low: {n_luminal}")

    if n_basal == 0 or n_luminal == 0:
        logger.error(f"Cannot run enrichment: Basal={n_basal}, Luminal={n_luminal}")
        return pd.DataFrame()

    n_programs = H.shape[1]
    rows = []

    for k in range(n_programs):
        h_k = H[:, k]
        basal_vals = h_k[basal_mask]
        luminal_vals = h_k[luminal_mask]

        mean_basal = float(np.mean(basal_vals))
        mean_luminal = float(np.mean(luminal_vals))

        # Mann-Whitney U
        stat, pval = mannwhitneyu(basal_vals, luminal_vals, alternative='two-sided')

        # Bonferroni correction
        pval_corrected = min(pval * n_programs, 1.0)

        # Fold change
        fc = (mean_basal + 1e-10) / (mean_luminal + 1e-10)

        # Specificity
        denom = mean_basal + mean_luminal
        specificity = abs(mean_basal - mean_luminal) / denom if denom > 0 else 0

        enriched_in = "KRT5-high" if mean_basal > mean_luminal else "KRT5-low"

        rows.append({
            'program': k,
            'mean_basal': mean_basal,
            'mean_luminal': mean_luminal,
            'fold_change': fc,
            'specificity': specificity,
            'pvalue_raw': pval,
            'pvalue_bonferroni': pval_corrected,
            'enriched_in': enriched_in,
            'significant': pval_corrected < 0.05 and max(fc, 1.0 / fc) > 1.5,
        })

        logger.info(f"  Program {k}: FC={fc:.2f} ({enriched_in}), "
                    f"p_bonf={pval_corrected:.2e}, sig={rows[-1]['significant']}")

    enrichment_df = pd.DataFrame(rows)
    return enrichment_df


def identify_mdk_program(W, gene_names, enrichment_df):
    """Find the MDK-containing program and extract context genes."""
    if "MDK" not in gene_names:
        logger.error("MDK not in gene list after filtering — cannot proceed")
        return None, None, pd.DataFrame()

    mdk_idx = gene_names.index("MDK")
    n_programs = W.shape[1]

    # MDK loading per program
    mdk_loadings = []
    for k in range(n_programs):
        mdk_loadings.append({
            'program': k,
            'mdk_loading': float(W[mdk_idx, k]),
            'basal_enriched': bool(enrichment_df.loc[
                enrichment_df['program'] == k, 'significant'
            ].values[0]) if len(enrichment_df) > 0 else False,
            'top_5_genes': ', '.join(
                [gene_names[i] for i in np.argsort(W[:, k])[::-1][:5]]
            ),
        })

    mdk_df = pd.DataFrame(mdk_loadings).sort_values('mdk_loading', ascending=False)
    logger.info("\nMDK loading by program:")
    for _, row in mdk_df.iterrows():
        logger.info(f"  Program {row['program']}: MDK={row['mdk_loading']:.4f}, "
                    f"Basal-enriched={row['basal_enriched']}, "
                    f"top5=[{row['top_5_genes']}]")

    # Pick MDK program: highest MDK loading
    mdk_program = int(mdk_df.iloc[0]['program'])
    mdk_loading = float(mdk_df.iloc[0]['mdk_loading'])
    logger.info(f"\nMDK program: {mdk_program} (loading={mdk_loading:.4f})")

    # Extract top 100 context genes (excluding MDK itself)
    program_loadings = W[:, mdk_program]
    sorted_idx = np.argsort(program_loadings)[::-1]

    context_genes = []
    for idx in sorted_idx:
        gene = gene_names[idx]
        if gene == "MDK":
            continue
        context_genes.append((gene, float(program_loadings[idx])))
        if len(context_genes) >= 100:
            break

    context_df = pd.DataFrame(context_genes, columns=['gene', 'loading'])

    # Flag known gene sets
    for gene_set_name, gene_list in [
        ("Secretory", SECRETORY_GENES),
        ("MDK_receptor", MDK_RECEPTORS),
        ("ER_related", ER_RELATED),
    ]:
        context_df[gene_set_name] = context_df['gene'].isin(gene_list)

    n_secretory = context_df['Secretory'].sum()
    n_receptor = context_df['MDK_receptor'].sum()
    n_er = context_df['ER_related'].sum()
    logger.info(f"Context genes (top 100): {n_secretory} secretory, "
                f"{n_receptor} MDK receptors, {n_er} ER-related")

    return mdk_program, mdk_df, context_df


def compute_spatial_autocorrelation(H, adata_cancer):
    """Compute Moran's I for each program using cell-level spatial coordinates."""
    logger.info("Computing spatial autocorrelation (Moran's I)...")

    if 'spatial' not in adata_cancer.obsm:
        logger.error("No spatial coordinates in adata_cancer.obsm — skipping Moran's I")
        return pd.DataFrame()

    # squidpy spatial_autocorr expects genes in .X (var_names), not .obs
    # Create a temporary AnnData with programs as "genes" and spatial coords
    n_programs = H.shape[1]
    program_names = [f'program_{k}' for k in range(n_programs)]

    import anndata
    adata_prog = anndata.AnnData(
        X=H,
        obs=adata_cancer.obs.copy(),
        var=pd.DataFrame(index=program_names),
        obsm={'spatial': adata_cancer.obsm['spatial'].copy()},
    )

    # Build spatial graph directly on the temp AnnData
    sq.gr.spatial_neighbors(adata_prog, coord_type='generic', n_neighs=15)
    sq.gr.spatial_autocorr(adata_prog, mode='moran')

    moran_df = adata_prog.uns['moranI'].copy()

    for k in range(n_programs):
        row = moran_df.loc[f'program_{k}']
        logger.info(f"  Program {k}: Moran's I={row['I']:.4f}, p={row['pval_norm']:.2e}")

    return moran_df


def write_outputs(W, H, gene_names, adata_cancer, enrichment_df, mdk_program,
                  mdk_df, context_df, moran_df, stability, k_results):
    """Write all CSV outputs and comparison report."""
    logger.info(f"Writing outputs to {OUTPUT_DIR}/")

    # Program loadings (genes x programs)
    loadings_df = pd.DataFrame(W, index=gene_names,
                               columns=[f'program_{k}' for k in range(W.shape[1])])
    loadings_df.to_csv(OUTPUT_DIR / "program_loadings.csv")

    # Program activities (cells x programs)
    cell_ids = list(adata_cancer.obs.index)
    activities_df = pd.DataFrame(H, index=cell_ids,
                                 columns=[f'program_{k}' for k in range(H.shape[1])])
    activities_df.to_csv(OUTPUT_DIR / "program_activities.csv")

    # Region enrichment
    enrichment_df.to_csv(OUTPUT_DIR / "region_enrichment_summary.csv", index=False)

    # MDK context genes
    context_df.to_csv(OUTPUT_DIR / "mdk_context_genes.csv", index=False)

    # Spatial autocorrelation
    moran_df.to_csv(OUTPUT_DIR / "spatial_autocorrelation.csv")

    # NMF stability
    stability_rows = []
    if stability.get('mdk_present'):
        for key, val in stability['w_correlations'].items():
            jaccard = stability['jaccard_top50'].get(key, np.nan)
            stability_rows.append({
                'comparison': key,
                'w_correlation': val,
                'jaccard_top50': jaccard,
            })
    stability_out = pd.DataFrame(stability_rows)
    stability_out.to_csv(OUTPUT_DIR / "nmf_stability.csv", index=False)

    # MDK program summary
    if mdk_df is not None:
        mdk_df.to_csv(OUTPUT_DIR / "mdk_program_loadings.csv", index=False)

    # K sensitivity
    k_rows = []
    for k, res in k_results.items():
        k_rows.append({'K': k, 'recon_error': res['recon_err']})
    pd.DataFrame(k_rows).to_csv(OUTPUT_DIR / "k_sensitivity.csv", index=False)

    # --- Comparison report ---
    write_comparison_report(
        W, gene_names, enrichment_df, mdk_program, context_df,
        moran_df, stability, adata_cancer,
    )

    logger.info(f"All outputs written to {OUTPUT_DIR}/")


def write_comparison_report(W, gene_names, enrichment_df, mdk_program,
                            context_df, moran_df, stability, adata_cancer):
    """Write human-readable comparison report."""
    lines = []
    lines.append("=" * 70)
    lines.append("MDK v3 SINGLE-CELL REANALYSIS — COMPARISON REPORT")
    lines.append("=" * 70)
    lines.append("")

    # Data summary
    n_basal = (adata_cancer.obs['cancer_subtype'] == 'KRT5-high').sum()
    n_luminal = (adata_cancer.obs['cancer_subtype'] == 'KRT5-low').sum()
    lines.append(f"Sample: {SAMPLE_NAME}")
    lines.append(f"Cancer cells: {adata_cancer.n_obs} ({n_basal} KRT5-high, {n_luminal} KRT5-low)")
    lines.append(f"Genes after filter: {len(gene_names)}")
    lines.append("")

    # MDK program
    if mdk_program is not None and "MDK" in gene_names:
        mdk_idx = gene_names.index("MDK")
        mdk_loading = W[mdk_idx, mdk_program]
        enrich_row = enrichment_df[enrichment_df['program'] == mdk_program].iloc[0]
        moran_row = moran_df.loc[f'program_{mdk_program}']

        lines.append("MDK PROGRAM RESULTS")
        lines.append("-" * 40)
        lines.append(f"MDK program: {mdk_program}")
        lines.append(f"MDK loading: {mdk_loading:.4f}")
        lines.append(f"KRT5-high enrichment: FC={enrich_row['fold_change']:.2f}, "
                     f"p_bonf={enrich_row['pvalue_bonferroni']:.2e}")
        lines.append(f"Enriched in: {enrich_row['enriched_in']}")
        lines.append(f"Significant (p<0.05, |FC|>1.5): {enrich_row['significant']}")
        lines.append(f"Moran's I: {moran_row['I']:.4f} (p={moran_row['pval_norm']:.2e})")
        lines.append("")

        # Context genes
        lines.append("TOP 20 CONTEXT GENES (co-loaded with MDK)")
        lines.append("-" * 40)
        for i, (_, row) in enumerate(context_df.head(20).iterrows()):
            flags = []
            if row.get('Secretory', False):
                flags.append('SECRETORY')
            if row.get('MDK_receptor', False):
                flags.append('RECEPTOR')
            if row.get('ER_related', False):
                flags.append('ER')
            flag_str = f"  [{', '.join(flags)}]" if flags else ""
            lines.append(f"  {i+1:2d}. {row['gene']:<12s} (loading: {row['loading']:.4f}){flag_str}")
        lines.append("")

        # Secretory gene summary
        secretory_in_top100 = context_df[context_df['Secretory']]['gene'].tolist()
        receptor_in_top100 = context_df[context_df['MDK_receptor']]['gene'].tolist()
        lines.append(f"Secretory genes in top 100: {secretory_in_top100 or 'none'}")
        lines.append(f"MDK receptors in top 100: {receptor_in_top100 or 'none'}")
        lines.append("")
    else:
        lines.append("MDK NOT FOUND IN FILTERED GENE SET")
        lines.append("")

    # Stability
    lines.append("NMF STABILITY")
    lines.append("-" * 40)
    if stability.get('mdk_present'):
        lines.append(f"All seeds same MDK program rank: {stability['all_seeds_same_rank']}")
        lines.append(f"Per-seed MDK programs: {stability['mdk_program_per_seed']}")
        lines.append(f"Mean W correlation: {stability['mean_w_correlation']:.3f}")
        lines.append(f"Mean Jaccard (top 50): {stability['mean_jaccard']:.3f}")
    else:
        lines.append("MDK not in gene list — stability not assessed")
    lines.append("")

    # All programs enrichment table
    lines.append("ALL PROGRAMS — REGION ENRICHMENT")
    lines.append("-" * 40)
    lines.append(enrichment_df.to_string(index=False))
    lines.append("")

    # Caveat
    lines.append("CAVEATS")
    lines.append("-" * 40)
    lines.append("- P-values are exploratory (single patient, spatially correlated cells)")
    lines.append("- KRT5-high/KRT5-low subtyping uses Module 3.5 protein gate when available, QP cell type as fallback")
    lines.append("- Three confounded changes vs original: unit, features, region definition")
    lines.append("- Compare direction and magnitude, not exact p-values, to original")
    lines.append("")

    # Comparison table
    lines.append("SIDE-BY-SIDE COMPARISON")
    lines.append("-" * 40)
    lines.append(f"{'Metric':<40s} {'Original (spot)':<20s} {'v3 (single-cell)':<20s}")
    lines.append("-" * 80)
    lines.append(f"{'Unit of analysis':<40s} {'~3000 spots':<20s} {f'{adata_cancer.n_obs} cells':<20s}")
    lines.append(f"{'Region definition':<40s} {'Keratin threshold':<20s} {'Cell type label':<20s}")
    lines.append(f"{'NMF input':<40s} {'Deconv layer (raw)':<20s} {'Size-norm SACE':<20s}")

    if mdk_program is not None and "MDK" in gene_names:
        mdk_idx = gene_names.index("MDK")
        enrich_row = enrichment_df[enrichment_df['program'] == mdk_program].iloc[0]
        moran_row = moran_df.loc[f'program_{mdk_program}']
        p_bonf = enrich_row['pvalue_bonferroni']
        fc_val = enrich_row['fold_change']
        moran_i = moran_row['I']
        stable = stability.get('all_seeds_same_rank', '?')
        mdk_load = W[mdk_idx, mdk_program]
        lines.append(f"{'MDK program rank':<40s} {'? (run vignette4)':<20s} {str(mdk_program):<20s}")
        lines.append(f"{'MDK loading':<40s} {'?':<20s} {mdk_load:<20.4f}")
        lines.append(f"{'Region enrichment p (Bonf)':<40s} {'?':<20s} {p_bonf:<20.2e}")
        lines.append(f"{'Region fold change':<40s} {'?':<20s} {fc_val:<20.2f}")
        lines.append(f"{'Morans I':<40s} {'?':<20s} {moran_i:<20.4f}")
        lines.append(f"{'Stable across seeds?':<40s} {'N/A':<20s} {str(stable):<20s}")

    lines.append("")

    report_text = "\n".join(lines)
    report_path = OUTPUT_DIR / "comparison_report.txt"
    report_path.write_text(report_text)
    logger.info(f"Comparison report: {report_path}")

    # Also print to stdout
    print("\n" + report_text)


if __name__ == "__main__":
    logger.info("=" * 60)
    logger.info("MDK v3 Single-Cell Reanalysis")
    logger.info("=" * 60)

    # Step 1-2: Load and preprocess
    adata_cancer, adata_nmf = load_and_preprocess()

    # Step 3: NMF with stability
    W, H, gene_names, seed_results, k_results, stability = run_nmf_stability(adata_nmf)

    # Step 4-5: Region enrichment
    enrichment_df = run_region_enrichment(H, adata_cancer)

    # Step 6: MDK program identification
    mdk_program, mdk_df, context_df = identify_mdk_program(W, gene_names, enrichment_df)

    # Step 7: Spatial autocorrelation
    moran_df = compute_spatial_autocorrelation(H, adata_cancer)

    # Step 8-9: Write all outputs + comparison report
    write_outputs(W, H, gene_names, adata_cancer, enrichment_df, mdk_program,
                  mdk_df, context_df, moran_df, stability, k_results)

    logger.info("MDK v3 reanalysis complete.")
