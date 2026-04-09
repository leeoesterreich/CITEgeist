#!/usr/bin/env python
"""
SACE Re-filter + MDK Cross-Patient Reanalysis

Re-runs SACE GEX allocation with relaxed gene filtering (mean>0.5 instead of
mean>1.1) to get ~7k genes, then runs MDK program discovery on each sample.
Produces per-sample reports and a cross-patient summary.

Usage: python run_sace_refilter_mdk.py <sample_index>  (0-11 for array jobs)
       python run_sace_refilter_mdk.py --summarize      (after all samples done)
"""

import argparse
import json
import logging
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

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================
PATIENT_DATA_ROOT = Path(
    "/ix1/alee/LO_LAB/General/Lab_Data/"
    "20250210_CITEGeistPublicData_GEO_Alex/processed_files"
)
V3_DIR = Path("output/morphology_assignment_v3")
OUTPUT_BASE = Path("output/mdk_v3_reanalysis_expanded")
OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]

MODEL_PROFILE_DICT = {
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "B_Cells": {"Major": ["MS4A1-1", "CD19-1"]},
    "Macrophages": {"Major": ["CD68-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "CD4_T_Cells": {"Major": ["CD3E-1", "CD4-1"]},
    # NOTE: Future re-runs should collapse to a single Epithelial type (EPCAM-1 only)
    # and use KRT5 as a Module 3.5 functional gate for basal-state identification.
    # Existing outputs validated equivalent to merged approach (mean |delta| I = 0.014).
    "Cancer_Luminal": {"Major": ["EPCAM-1"]},
    "Cancer_Basal": {"Major": ["EPCAM-1", "KRT5-1"]},
    "Dendritic_Cells": {"Major": ["HLA-DRA-1"]},
}

# Both subtypes are pooled as a single Epithelial population.
# Within-Epithelial stratification uses the KRT5 protein gate (Module 3.5).
CANCER_TYPES = ['Cancer_Basal', 'Cancer_Luminal']

# Relaxed gene filter: Option A (min_counts=50, 1% nonzero, mean>0.5)
FILTER_MIN_COUNTS = 50
FILTER_NONZERO_PCT = 0.01
FILTER_MEAN_THRESH = 0.5

# NMF settings
NMF_SEEDS = [42, 123, 456, 789, 1024]
PRIMARY_K = 5
PRIMARY_SEED = 42

# Gene sets for flagging
SECRETORY_GENES = ["HSP90B1", "XBP1", "CALR", "CANX", "PDIA4", "PDIA6", "ERO1A", "SEC61A1"]
MDK_RECEPTORS = ["SDC4", "NCL", "PTPRZ1", "LRP1", "GPC1"]
ER_RELATED = ["ESR1", "PGR", "GREB1", "TFF1", "FOXA1", "GATA3"]


# =============================================================================
# Stage 1: Re-run SACE with relaxed gene filter
# =============================================================================
def rerun_sace(sample_name):
    """Re-run SACE allocation with relaxed gene filtering."""
    logger.info(f"=== SACE re-run: {sample_name} ===")

    v3_sample_dir = V3_DIR / sample_name
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"

    # Load raw Visium data
    adata = sq.read.visium(
        str(patient_dir),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False,
    )
    logger.info(f"Loaded raw: {adata.shape[0]} spots x {adata.shape[1]} features")

    # Initialize CITEgeist model
    repo_root = str(Path(__file__).resolve().parent.parent.parent)
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from CITEgeist.model import CitegeistModel

    output_dir = OUTPUT_BASE / sample_name
    output_dir.mkdir(parents=True, exist_ok=True)

    sace_model = CitegeistModel(
        sample_name=sample_name,
        adata=adata,
        output_folder=str(output_dir),
    )
    sace_model.split_adata()

    # Relaxed gene filter
    sace_model.filter_gex(
        min_counts=FILTER_MIN_COUNTS,
        nonzero_percentage=FILTER_NONZERO_PCT,
        mean_expression_threshold=FILTER_MEAN_THRESH,
    )
    n_genes = sace_model.gene_expression_adata.n_vars
    n_spots = sace_model.gene_expression_adata.n_obs
    logger.info(f"After relaxed filter: {n_genes} genes, {n_spots} spots")

    sace_model.preprocess_gex()
    sace_model.preprocess_antibody()
    sace_model.load_cell_profile_dict(MODEL_PROFILE_DICT)

    # Load v3 ensemble proportions and cell assignments
    ensemble_df = pd.read_csv(v3_sample_dir / "ensemble_proportions.csv", index_col=0)
    assignments_df = pd.read_csv(v3_sample_dir / "cell_assignments.csv")

    # Subset to common spots
    common_spots = ensemble_df.index.intersection(
        sace_model.gene_expression_adata.obs_names
    )
    if len(common_spots) < len(ensemble_df):
        logger.warning(f"  {len(common_spots)}/{len(ensemble_df)} ensemble spots in GEX")
    sace_model.gene_expression_adata = sace_model.gene_expression_adata[common_spots].copy()
    sace_model.antibody_capture_adata = sace_model.antibody_capture_adata[common_spots].copy()
    ensemble_df = ensemble_df.loc[common_spots]

    # Inject proportions
    sace_model.results["cell_prop"] = ensemble_df

    # Build cell assignments and spot map
    cell_assignments = dict(zip(
        assignments_df["nucleus_id"].astype(str),
        assignments_df["assigned_type"],
    ))

    nucleus_spot_map = pd.read_csv(v3_sample_dir / "nucleus_spot_mapping.csv")
    cell_spot_map = nucleus_spot_map.rename(columns={
        "nucleus_id": "cell_id",
        "centroid_x": "x",
        "centroid_y": "y",
    })
    assigned_ids = set(assignments_df["nucleus_id"].values)
    cell_spot_map = cell_spot_map[cell_spot_map["cell_id"].isin(assigned_ids)].copy()
    cell_spot_map["cell_id"] = cell_spot_map["cell_id"].astype(str)

    # Run SACE
    spot_type_gex, cell_adata, diagnostics = sace_model.run_sace_allocation(
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        max_iter=1,
    )

    # Save expanded single-cell h5ad
    sc_path = output_dir / f"{sample_name}_single_cell.h5ad"
    cell_adata.write_h5ad(str(sc_path))
    logger.info(f"Saved: {cell_adata.shape[0]} cells x {cell_adata.shape[1]} genes -> {sc_path}")

    return cell_adata


# =============================================================================
# Stage 2: MDK NMF analysis (same as run_mdk_v3_reanalysis.py but parameterized)
# =============================================================================
def run_mdk_analysis(sample_name, cell_adata=None):
    """Run MDK NMF program discovery on one sample."""
    logger.info(f"=== MDK analysis: {sample_name} ===")

    output_dir = OUTPUT_BASE / sample_name

    if cell_adata is None:
        sc_path = output_dir / f"{sample_name}_single_cell.h5ad"
        cell_adata = sc.read_h5ad(sc_path)
        logger.info(f"Loaded: {cell_adata.shape}")

    # Subset to cancer cells
    cancer_mask = cell_adata.obs['cell_type'].isin(CANCER_TYPES)
    n_cancer = cancer_mask.sum()
    if n_cancer < 20:
        logger.warning(f"Only {n_cancer} cancer cells — skipping MDK analysis")
        return None

    adata_cancer = cell_adata[cancer_mask].copy()
    krt5_cols = [c for c in adata_cancer.obs.columns
                 if 'KRT5' in c and c.endswith('_gate')]
    if krt5_cols:
        krt5_pos = adata_cancer.obs[krt5_cols].max(axis=1).fillna(0).astype(bool)
        n_basal = int(krt5_pos.sum())
        n_luminal = int((~krt5_pos).sum())
        logger.info(f"Cancer cells: {adata_cancer.n_obs} ({n_basal} KRT5-high, {n_luminal} KRT5-low)")
    else:
        n_basal = int((adata_cancer.obs['cell_type'] == 'Cancer_Basal').sum())
        n_luminal = int((adata_cancer.obs['cell_type'] == 'Cancer_Luminal').sum())
        logger.info(f"Cancer cells: {adata_cancer.n_obs} ({n_basal} KRT5-high proxy, {n_luminal} KRT5-low proxy)")

    # Store raw, normalize
    adata_cancer.raw = adata_cancer.copy()
    sc.pp.normalize_total(adata_cancer, target_sum=1e4)

    # Gene detection filter (5% of cancer cells)
    min_cells = int(0.05 * adata_cancer.n_obs)
    X_dense = adata_cancer.X.toarray() if hasattr(adata_cancer.X, 'toarray') else np.asarray(adata_cancer.X)
    detection_counts = (X_dense > 0).sum(axis=0)
    gene_mask = detection_counts >= min_cells
    adata_nmf = adata_cancer[:, gene_mask].copy()
    logger.info(f"Gene filter: {adata_cancer.n_vars} -> {adata_nmf.n_vars} genes (>={min_cells} cells)")

    # Check target genes
    gene_names = list(adata_nmf.var_names)
    mdk_present = "MDK" in gene_names
    secretory_present = [g for g in SECRETORY_GENES if g in gene_names]
    secretory_missing = [g for g in SECRETORY_GENES if g not in gene_names]
    logger.info(f"MDK: {'present' if mdk_present else 'MISSING'}")
    logger.info(f"Secretory: present={secretory_present}, missing={secretory_missing}")

    if not mdk_present:
        logger.warning("MDK not in filtered genes — writing null result")
        _write_null_result(output_dir, sample_name, adata_cancer, gene_names,
                           secretory_present, secretory_missing)
        return None

    # NMF (multi-seed at K=5)
    X = np.maximum(np.asarray(adata_nmf.X.toarray() if hasattr(adata_nmf.X, 'toarray')
                              else adata_nmf.X), 0)

    seed_results = {}
    for seed in NMF_SEEDS:
        model = NMF(n_components=PRIMARY_K, init='nndsvda', random_state=seed, max_iter=500)
        H = model.fit_transform(X)
        W = model.components_.T
        seed_results[seed] = {'W': W, 'H': H, 'recon_err': model.reconstruction_err_}

    W = seed_results[PRIMARY_SEED]['W']
    H = seed_results[PRIMARY_SEED]['H']

    # Stability
    mdk_idx = gene_names.index("MDK")
    mdk_programs_per_seed = {}
    mdk_w_vectors = {}
    mdk_top50 = {}
    for seed in sorted(seed_results.keys()):
        Ws = seed_results[seed]['W']
        best = int(np.argmax(Ws[mdk_idx, :]))
        mdk_programs_per_seed[seed] = best
        mdk_w_vectors[seed] = Ws[:, best]
        top50_idx = np.argsort(Ws[:, best])[::-1][:50]
        mdk_top50[seed] = set(gene_names[i] for i in top50_idx)

    ref = sorted(seed_results.keys())[0]
    w_corrs = [pearsonr(mdk_w_vectors[ref], mdk_w_vectors[s])[0]
               for s in sorted(seed_results.keys())[1:]]
    jaccards = [len(mdk_top50[ref] & mdk_top50[s]) / len(mdk_top50[ref] | mdk_top50[s])
                for s in sorted(seed_results.keys())[1:]]
    all_same_rank = len(set(mdk_programs_per_seed.values())) == 1

    # Region annotation: use KRT5 protein gate if available, fall back to QP type
    krt5_cols = [c for c in adata_cancer.obs.columns
                 if 'KRT5' in c and c.endswith('_gate')]
    if krt5_cols:
        krt5_pos = adata_cancer.obs[krt5_cols].max(axis=1).fillna(0).astype(bool)
        adata_cancer.obs['cancer_subtype'] = np.where(krt5_pos, 'KRT5-high', 'KRT5-low')
    else:
        adata_cancer.obs['cancer_subtype'] = adata_cancer.obs['cell_type'].map({
            'Cancer_Basal': 'KRT5-high', 'Cancer_Luminal': 'KRT5-low',
        })

    # Region enrichment
    subtypes = adata_cancer.obs['cancer_subtype'].values
    basal_mask = subtypes == 'KRT5-high'
    luminal_mask = subtypes == 'KRT5-low'

    enrichment_rows = []
    for k in range(PRIMARY_K):
        h_k = H[:, k]
        if basal_mask.sum() == 0 or luminal_mask.sum() == 0:
            enrichment_rows.append({
                'program': k, 'mean_basal': 0, 'mean_luminal': 0,
                'fold_change': 1, 'pvalue_bonferroni': 1, 'enriched_in': 'N/A',
                'significant': False,
            })
            continue
        mb = float(np.mean(h_k[basal_mask]))
        ml = float(np.mean(h_k[luminal_mask]))
        _, pval = mannwhitneyu(h_k[basal_mask], h_k[luminal_mask], alternative='two-sided')
        pval_corr = min(pval * PRIMARY_K, 1.0)
        fc = (mb + 1e-10) / (ml + 1e-10)
        enrichment_rows.append({
            'program': k, 'mean_basal': mb, 'mean_luminal': ml,
            'fold_change': fc, 'pvalue_bonferroni': pval_corr,
            'enriched_in': 'KRT5-high' if mb > ml else 'KRT5-low',
            'significant': pval_corr < 0.05 and max(fc, 1.0 / fc) > 1.5,
        })
    enrichment_df = pd.DataFrame(enrichment_rows)

    # MDK program
    mdk_program = int(np.argmax(W[mdk_idx, :]))
    mdk_loading = float(W[mdk_idx, mdk_program])
    mdk_row = enrichment_df[enrichment_df['program'] == mdk_program].iloc[0]

    # Context genes
    prog_loadings = W[:, mdk_program]
    sorted_idx = np.argsort(prog_loadings)[::-1]
    context_genes = []
    for idx in sorted_idx:
        if gene_names[idx] == "MDK":
            continue
        context_genes.append((gene_names[idx], float(prog_loadings[idx])))
        if len(context_genes) >= 100:
            break
    context_df = pd.DataFrame(context_genes, columns=['gene', 'loading'])
    for name, gl in [("Secretory", SECRETORY_GENES), ("MDK_receptor", MDK_RECEPTORS),
                     ("ER_related", ER_RELATED)]:
        context_df[name] = context_df['gene'].isin(gl)

    # Moran's I
    import anndata
    moran_i_val = np.nan
    moran_p_val = np.nan
    if 'spatial' in adata_cancer.obsm:
        program_names = [f'program_{k}' for k in range(PRIMARY_K)]
        adata_prog = anndata.AnnData(
            X=H, obs=adata_cancer.obs.copy(),
            var=pd.DataFrame(index=program_names),
            obsm={'spatial': adata_cancer.obsm['spatial'].copy()},
        )
        sq.gr.spatial_neighbors(adata_prog, coord_type='generic', n_neighs=15)
        sq.gr.spatial_autocorr(adata_prog, mode='moran')
        moran_df = adata_prog.uns['moranI']
        moran_i_val = float(moran_df.loc[f'program_{mdk_program}', 'I'])
        moran_p_val = float(moran_df.loc[f'program_{mdk_program}', 'pval_norm'])

    # Top 5 genes in MDK program
    top5_idx = np.argsort(W[:, mdk_program])[::-1][:5]
    top5_genes = [gene_names[i] for i in top5_idx]

    # Build result dict
    result = {
        'sample': sample_name,
        'n_cancer_cells': int(adata_cancer.n_obs),
        'n_basal': int(n_basal),
        'n_luminal': int(n_luminal),
        'n_genes_sace': int(cell_adata.n_vars),
        'n_genes_nmf': len(gene_names),
        'mdk_present': True,
        'mdk_program': mdk_program,
        'mdk_loading': mdk_loading,
        'mdk_fc': float(mdk_row['fold_change']),
        'mdk_pval_bonf': float(mdk_row['pvalue_bonferroni']),
        'mdk_enriched_in': mdk_row['enriched_in'],
        'mdk_significant': bool(mdk_row['significant']),
        'moran_i': moran_i_val,
        'moran_p': moran_p_val,
        'stable_across_seeds': all_same_rank,
        'mean_w_corr': float(np.mean(w_corrs)) if w_corrs else np.nan,
        'mean_jaccard': float(np.mean(jaccards)) if jaccards else np.nan,
        'top5_genes': top5_genes,
        'secretory_in_top100': context_df[context_df['Secretory']]['gene'].tolist(),
        'receptor_in_top100': context_df[context_df['MDK_receptor']]['gene'].tolist(),
        'er_in_top100': context_df[context_df['ER_related']]['gene'].tolist(),
        'secretory_in_sace': [g for g in SECRETORY_GENES if g in cell_adata.var_names],
        'secretory_missing': [g for g in SECRETORY_GENES if g not in cell_adata.var_names],
    }

    # Save per-sample outputs
    enrichment_df.to_csv(output_dir / "region_enrichment_summary.csv", index=False)
    context_df.to_csv(output_dir / "mdk_context_genes.csv", index=False)
    with open(output_dir / "mdk_result.json", 'w') as f:
        json.dump(result, f, indent=2, default=str)

    logger.info(f"MDK program {mdk_program}: loading={mdk_loading:.2f}, "
                f"FC={result['mdk_fc']:.2f}, p_bonf={result['mdk_pval_bonf']:.2e}, "
                f"Moran's I={moran_i_val:.4f}, top5={top5_genes}")
    logger.info(f"Secretory in SACE: {result['secretory_in_sace']}")
    logger.info(f"Secretory in top 100 context: {result['secretory_in_top100']}")

    return result


def _write_null_result(output_dir, sample_name, adata_cancer, gene_names,
                       secretory_present, secretory_missing):
    """Write a null result JSON when MDK is not in the filtered gene set."""
    result = {
        'sample': sample_name,
        'n_cancer_cells': int(adata_cancer.n_obs),
        'n_genes_nmf': len(gene_names),
        'mdk_present': False,
        'secretory_in_nmf': secretory_present,
        'secretory_missing': secretory_missing,
    }
    with open(output_dir / "mdk_result.json", 'w') as f:
        json.dump(result, f, indent=2)


# =============================================================================
# Stage 3: Cross-patient summary
# =============================================================================
def summarize_cross_patient():
    """Aggregate MDK results across all samples."""
    logger.info("=== Cross-patient MDK summary ===")

    results = []
    for sample in SAMPLES:
        result_path = OUTPUT_BASE / sample / "mdk_result.json"
        if not result_path.exists():
            logger.warning(f"  {sample}: no result file")
            continue
        with open(result_path) as f:
            results.append(json.load(f))

    if not results:
        logger.error("No results found")
        return

    # Summary table
    lines = []
    lines.append("=" * 90)
    lines.append("MDK CROSS-PATIENT SUMMARY — v3 Single-Cell with Expanded Gene Filter")
    lines.append(f"Gene filter: min_counts={FILTER_MIN_COUNTS}, "
                 f"nonzero>={FILTER_NONZERO_PCT:.0%}, mean>{FILTER_MEAN_THRESH}")
    lines.append("=" * 90)
    lines.append("")

    # Per-sample table
    header = (f"{'Sample':<28s} {'Cells':>5s} {'Genes':>5s} {'MDK?':>4s} "
              f"{'Prog':>4s} {'FC':>6s} {'p_bonf':>10s} {'Enrich':>12s} "
              f"{'Moran':>6s} {'Stable':>6s} {'Top gene':>10s}")
    lines.append(header)
    lines.append("-" * len(header))

    n_mdk_found = 0
    n_basal_enriched = 0
    n_significant = 0
    n_stable = 0
    all_secretory_in_context = []
    all_top5 = []

    for r in results:
        sample_short = r['sample'].replace('HCC22-088-', '')
        if not r.get('mdk_present', False):
            lines.append(f"{sample_short:<28s} {r.get('n_cancer_cells', '?'):>5} "
                         f"{r.get('n_genes_nmf', '?'):>5} {'NO':>4s}")
            continue

        n_mdk_found += 1
        if r['mdk_enriched_in'] == 'KRT5-high':
            n_basal_enriched += 1
        if r['mdk_significant']:
            n_significant += 1
        if r['stable_across_seeds']:
            n_stable += 1
        all_secretory_in_context.extend(r.get('secretory_in_top100', []))
        all_top5.append(r.get('top5_genes', []))

        top1 = r['top5_genes'][0] if r.get('top5_genes') else '?'
        lines.append(
            f"{sample_short:<28s} {r['n_cancer_cells']:>5d} {r['n_genes_nmf']:>5d} "
            f"{'YES':>4s} {r['mdk_program']:>4d} {r['mdk_fc']:>6.2f} "
            f"{r['mdk_pval_bonf']:>10.2e} {r['mdk_enriched_in']:>12s} "
            f"{r['moran_i']:>6.3f} {str(r['stable_across_seeds']):>6s} {top1:>10s}"
        )

    lines.append("")
    lines.append("CROSS-PATIENT STATISTICS")
    lines.append("-" * 40)
    lines.append(f"Samples with MDK in NMF: {n_mdk_found}/{len(results)}")
    lines.append(f"MDK program Basal-enriched: {n_basal_enriched}/{n_mdk_found}")
    lines.append(f"Significant (p<0.05, |FC|>1.5): {n_significant}/{n_mdk_found}")
    lines.append(f"Stable across seeds: {n_stable}/{n_mdk_found}")
    lines.append("")

    # Secretory gene recovery
    lines.append("SECRETORY GENE RECOVERY (expanded filter)")
    lines.append("-" * 40)
    for r in results:
        if not r.get('mdk_present'):
            continue
        sample_short = r['sample'].replace('HCC22-088-', '')
        in_sace = r.get('secretory_in_sace', [])
        in_ctx = r.get('secretory_in_top100', [])
        lines.append(f"  {sample_short}: in SACE={in_sace}, in MDK context={in_ctx}")
    lines.append("")

    # Consensus top genes
    if all_top5:
        from collections import Counter
        gene_counts = Counter()
        for t5 in all_top5:
            gene_counts.update(t5)
        lines.append("CONSENSUS TOP GENES (across samples)")
        lines.append("-" * 40)
        for gene, count in gene_counts.most_common(15):
            lines.append(f"  {gene:<15s} appears in {count}/{n_mdk_found} samples")
        lines.append("")

    lines.append("CONCLUSION")
    lines.append("-" * 40)
    if n_significant > n_mdk_found / 2:
        lines.append("MDK secretory program is CONSISTENT across patients at single-cell level.")
    elif n_basal_enriched > n_mdk_found / 2:
        lines.append("MDK is Basal-enriched in most patients but not always significant.")
    else:
        lines.append("MDK finding is NOT consistent across patients.")

    report = "\n".join(lines)
    report_path = OUTPUT_BASE / "cross_patient_summary.txt"
    report_path.write_text(report)
    print("\n" + report)
    logger.info(f"Summary written to {report_path}")

    # Also save as CSV for easy parsing
    summary_df = pd.DataFrame(results)
    summary_df.to_csv(OUTPUT_BASE / "cross_patient_summary.csv", index=False)


# =============================================================================
# Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_index', nargs='?', type=int, default=None,
                        help='Sample index (0-11) for array jobs')
    parser.add_argument('--summarize', action='store_true',
                        help='Generate cross-patient summary from existing results')
    args = parser.parse_args()

    if args.summarize:
        summarize_cross_patient()
        return

    if args.sample_index is None:
        # Default: use SLURM_ARRAY_TASK_ID
        import os
        task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
        if task_id is not None:
            args.sample_index = int(task_id)
        else:
            logger.error("No sample_index and no SLURM_ARRAY_TASK_ID")
            sys.exit(1)

    if args.sample_index < 0 or args.sample_index >= len(SAMPLES):
        logger.error(f"Invalid sample index {args.sample_index} (0-{len(SAMPLES)-1})")
        sys.exit(1)

    sample_name = SAMPLES[args.sample_index]
    logger.info(f"Processing sample {args.sample_index}: {sample_name}")

    # Stage 1: SACE with expanded genes
    cell_adata = rerun_sace(sample_name)

    # Stage 2: MDK analysis
    run_mdk_analysis(sample_name, cell_adata)

    logger.info(f"=== {sample_name} complete ===")


if __name__ == "__main__":
    main()
