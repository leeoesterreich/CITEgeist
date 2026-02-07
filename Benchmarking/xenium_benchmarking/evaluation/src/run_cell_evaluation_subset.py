#!/usr/bin/env python3
"""
Run three-tier cell classification evaluation on the intersection of
Module 3 gating output (50K cells from 5 regions) and ScType annotations (465K cells).

This script handles the index alignment issue: Module 3 has 50K cells, ScType has 465K,
and Xenium GEX has 465K. We subset GEX to the intersection so Tier 2 positional
indexing works correctly.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / 'Benchmarking/xenium_pseudovisium/src'))
sys.path.insert(0, str(REPO_ROOT / 'Benchmarking/xenium_benchmarking/evaluation/src'))

from load_xenium import load_xenium_data, split_gex_protein
from evaluate_cell_classification import (
    MODULE3_LINEAGE_MAP, map_to_lineage,
    evaluate_concordance, evaluate_refinement, evaluate_discordance,
    plot_confusion_matrix, plot_refinement_summary,
)

logger = logging.getLogger(__name__)

XENIUM_DATA_DIR = Path(
    '/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/'
    'Xenium_RNA_Proteomic_RenalCellCarcinoma'
)


def concatenate_region_assignments(gating_output_dir):
    gating_dir = Path(gating_output_dir)
    dfs = []
    for region_dir in sorted(gating_dir.glob('region_*')):
        csv_path = region_dir / 'assignments_no_neg.csv'
        if csv_path.exists():
            df = pd.read_csv(csv_path, index_col=0)
            logger.info(f'Loaded {len(df)} cells from {region_dir.name}')
            dfs.append(df)
    if not dfs:
        raise FileNotFoundError(f'No files found in {gating_output_dir}/region_*/')
    combined = pd.concat(dfs)
    n_dup = combined.index.duplicated().sum()
    if n_dup > 0:
        logger.warning(f'{n_dup} duplicate cell IDs, keeping first')
        combined = combined[~combined.index.duplicated(keep='first')]
    logger.info(f'Combined Module 3 assignments: {len(combined)} cells')
    return combined


def main(gating_output_dir, sctype_dir, output_dir, seurat_dir=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(exist_ok=True)
    results = {'metadata': {}}

    # Step 1: Concatenate Module 3
    logger.info('Step 1: Loading Module 3 gating assignments')
    m3_df = concatenate_region_assignments(gating_output_dir)
    m3_df.to_csv(output_dir / 'combined_module3_assignments.csv')
    m3_assignments = m3_df['assigned_type']
    m3_lineage = map_to_lineage(m3_assignments, MODULE3_LINEAGE_MAP)

    # Step 2: Load ScType
    logger.info('Step 2: Loading ScType annotations')
    sctype_df = pd.read_csv(Path(sctype_dir) / 'sctype_annotations.csv', index_col=0)
    sctype_assignments = sctype_df['sctype_annotation']
    sctype_lineage = sctype_df['sctype_lineage']
    logger.info(f'ScType: {len(sctype_df)} cells')

    # Step 3: Find intersection
    common_cells = m3_assignments.index.intersection(sctype_assignments.index)
    logger.info(f'Module 3: {len(m3_assignments)}, ScType: {len(sctype_assignments)}, Intersection: {len(common_cells)}')
    results['metadata'].update({
        'n_module3_cells': int(len(m3_assignments)),
        'n_sctype_cells': int(len(sctype_assignments)),
        'n_intersection_cells': int(len(common_cells)),
    })
    if len(common_cells) == 0:
        logger.error('No overlapping cells!')
        return

    # Step 4: Load Xenium GEX, subset to intersection
    logger.info('Step 4: Loading Xenium GEX data')
    adata = load_xenium_data(XENIUM_DATA_DIR)
    adata_gex, _ = split_gex_protein(adata)
    xenium_cells = set(adata_gex.obs_names)
    final_common = sorted(set(common_cells) & xenium_cells)
    logger.info(f'Xenium: {len(xenium_cells)} cells, Final eval set: {len(final_common)}')
    results['metadata']['n_xenium_cells'] = int(len(xenium_cells))
    results['metadata']['n_final_evaluation_cells'] = int(len(final_common))
    if len(final_common) == 0:
        logger.error('No cells in 3-way intersection!')
        return

    # Subset GEX so positional indices align with DataFrames
    adata_gex_sub = adata_gex[final_common, :].copy()
    sc.pp.normalize_total(adata_gex_sub, target_sum=1e4)
    sc.pp.log1p(adata_gex_sub)
    gex_data = adata_gex_sub.X
    if hasattr(gex_data, 'toarray'):
        gex_data = gex_data.toarray()
    gex_data = np.asarray(gex_data, dtype=np.float64)
    gene_names = list(adata_gex_sub.var_names)

    m3_asgn = m3_assignments.loc[final_common]
    m3_lin = m3_lineage.loc[final_common]
    sc_asgn = sctype_assignments.loc[final_common]
    sc_lin = sctype_lineage.loc[final_common]

    logger.info('Module 3 distribution:')
    for ct, n in m3_asgn.value_counts().items():
        logger.info(f'  {ct}: {n} ({100*n/len(m3_asgn):.1f}%)')
    logger.info('ScType distribution (top 10):')
    for ct, n in sc_asgn.value_counts().head(10).items():
        logger.info(f'  {ct}: {n} ({100*n/len(sc_asgn):.1f}%)')

    # Tier 1: Concordance
    logger.info('Tier 1: Concordance')
    t1 = evaluate_concordance(m3_lin, sc_lin, 'ScType')
    results['tier1_concordance_sctype'] = t1
    logger.info(f'Agreement: {100*t1["agreement_rate"]:.1f}%')
    plot_confusion_matrix(t1, str(figures_dir / 'concordance_sctype.png'),
                          'Module 3 vs ScType: Lineage Concordance')

    # Tier 2: Refinement
    logger.info('Tier 2: Refinement')
    t2 = evaluate_refinement(m3_asgn, sc_asgn, gex_data, gene_names, 'ScType')
    results['tier2_refinement_sctype'] = t2
    logger.info(f'Validated: {t2["n_subtypes_validated"]}/{t2["n_subtypes_tested"]}')
    plot_refinement_summary(t2, str(figures_dir / 'refinement_sctype.png'))

    # Tier 3: Discordance
    logger.info('Tier 3: Discordance')
    t3 = evaluate_discordance(m3_lin, sc_lin, m3_asgn, sc_asgn, 'ScType')
    results['tier3_discordance_sctype'] = t3
    logger.info(f'Discordance: {100*t3["discordance_rate"]:.1f}%')

    # Optional: Seurat
    if seurat_dir:
        sp = Path(seurat_dir) / 'seurat_label_transfer.csv'
        if sp.exists():
            logger.info('Evaluating against Seurat')
            sdf = pd.read_csv(sp, index_col=0)
            sc2 = sorted(set(final_common) & set(sdf.index))
            logger.info(f'Seurat overlap: {len(sc2)}')
            if len(sc2) > 100:
                sa = sdf.loc[sc2, 'predicted_type']
                sl = sdf.loc[sc2, 'predicted_lineage']
                ml = m3_lineage.loc[sc2]
                ma = m3_assignments.loc[sc2]
                c = evaluate_concordance(ml, sl, 'Seurat')
                results['tier1_concordance_seurat'] = c
                plot_confusion_matrix(c, str(figures_dir / 'concordance_seurat.png'),
                                      'Module 3 vs Seurat: Lineage Concordance')
                ag = adata_gex[sc2, :].copy()
                sc.pp.normalize_total(ag, target_sum=1e4)
                sc.pp.log1p(ag)
                gd = ag.X
                if hasattr(gd, 'toarray'):
                    gd = gd.toarray()
                gd = np.asarray(gd, dtype=np.float64)
                r = evaluate_refinement(ma, sa.loc[sc2], gd, gene_names, 'Seurat')
                results['tier2_refinement_seurat'] = r
                plot_refinement_summary(r, str(figures_dir / 'refinement_seurat.png'))
                d = evaluate_discordance(ml, sl, ma, sa, 'Seurat')
                results['tier3_discordance_seurat'] = d

    # Save
    rp = output_dir / 'cell_classification_evaluation.json'
    with open(rp, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    # Print summary
    print('\n' + '=' * 70)
    print('THREE-TIER CELL CLASSIFICATION EVALUATION SUMMARY')
    print('=' * 70)
    print(f'\nEvaluation set: {len(final_common)} cells')

    for tk, td in results.items():
        if tk == 'metadata':
            continue
        if tk.startswith('tier1'):
            m = td['method']
            print(f'\n--- Tier 1: Concordance ({m}) ---')
            print(f'  Agreement: {100*td["agreement_rate"]:.1f}%')
            print(f'  Valid cells: {td["n_valid"]:,}, M3 Unassigned: {td["n_module3_unassigned"]:,}')
            for lin, met in td.get('per_lineage', {}).items():
                print(f'  {lin:12s}: P={met["precision"]:.3f} R={met["recall"]:.3f} F1={met["f1"]:.3f}')
        elif tk.startswith('tier2'):
            m = td['method']
            print(f'\n--- Tier 2: Refinement ({m}) ---')
            print(f'  Validated: {td["n_subtypes_validated"]}/{td["n_subtypes_tested"]}')
            for sub, info in td.get('per_subtype', {}).items():
                st = info['validation_status']
                n = info['n_refined_cells']
                print(f'  {sub:22s}: {st} ({n} cells)')
                for mk, ev in info.get('marker_evidence', {}).items():
                    d2 = 'HIGH' if ev.get('validated') else 'LOW'
                    print(f'    {mk}: {d2} (mean={ev["refined_mean"]:.2f}, pct={ev["pct_above_median"]:.2f})')
        elif tk.startswith('tier3'):
            m = td['method']
            print(f'\n--- Tier 3: Discordance ({m}) ---')
            print(f'  Rate: {100*td["discordance_rate"]:.1f}%')
            print(f'  Discordant: {td["n_discordant"]:,}/{td["n_valid"]:,}')
            for p, c in list(td.get('top_discordance_pairs', {}).items())[:5]:
                print(f'    {p}: {c}')

    print('\n' + '=' * 70)
    print(f'Results: {rp}')
    print(f'Figures: {figures_dir}')
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Three-tier cell classification eval (subset)')
    parser.add_argument('--gating_output_dir', type=str, required=True)
    parser.add_argument('--sctype_dir', type=str, required=True)
    parser.add_argument('--seurat_dir', type=str, default=None)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )
    main(args.gating_output_dir, args.sctype_dir, args.output_dir, args.seurat_dir)
