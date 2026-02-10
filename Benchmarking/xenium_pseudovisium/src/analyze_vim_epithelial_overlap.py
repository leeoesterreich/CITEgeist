"""
Analyze vimentin/epithelial overlap and fibroblast gate breakdown.

Questions:
A. How many Vimentin+ cells are also PanCK+ or E-Cadherin+?
B. What fraction of "Fibroblasts" came from Gate 7 (αSMA) vs Gate 8 (Vimentin+)?
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

def get_threshold(expr_df: pd.DataFrame, marker: str, percentile: float = 50) -> float:
    """Get threshold as percentile of non-zero values."""
    vals = expr_df[marker]
    nonzero = vals[vals > 0]
    if len(nonzero) > 0:
        return np.percentile(nonzero, percentile)
    return 0

def main():
    # Load protein data
    xenium_dir = Path('/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma')

    print("Loading Xenium protein data...")
    adata = sc.read_10x_h5(xenium_dir / 'cell_feature_matrix.h5', gex_only=False)
    protein_mask = adata.var['feature_types'] == 'Protein Expression'
    adata_protein = adata[:, protein_mask].copy()

    X = adata_protein.X.toarray() if hasattr(adata_protein.X, 'toarray') else adata_protein.X
    proteins = list(adata_protein.var_names)
    cell_ids = [str(x) for x in adata_protein.obs_names]

    expr_df = pd.DataFrame(X, index=cell_ids, columns=proteins)
    print(f"Loaded {len(expr_df)} cells x {len(proteins)} proteins")
    print(f"Proteins: {proteins}")

    # Calculate thresholds (same as create_protein_gt.py)
    CD3E_thresh = get_threshold(expr_df, 'CD3E', 50)
    CD68_thresh = get_threshold(expr_df, 'CD68', 50)
    CD31_thresh = get_threshold(expr_df, 'CD31', 50)
    PanCK_thresh = get_threshold(expr_df, 'PanCK', 25)
    ECad_thresh = get_threshold(expr_df, 'E-Cadherin', 90)
    alphaSMA_thresh = get_threshold(expr_df, 'alphaSMA', 75)
    Vim_thresh = get_threshold(expr_df, 'Vimentin', 50)

    print("\n" + "="*60)
    print("THRESHOLDS")
    print("="*60)
    print(f"PanCK (p25 of nonzero): {PanCK_thresh:.2f}")
    print(f"E-Cadherin (p90 of nonzero): {ECad_thresh:.2f}")
    print(f"alphaSMA (p75 of nonzero): {alphaSMA_thresh:.2f}")
    print(f"Vimentin (p50 of nonzero): {Vim_thresh:.2f}")

    # Define gates
    panck_pos = expr_df['PanCK'] > PanCK_thresh
    ecad_high = expr_df['E-Cadherin'] > ECad_thresh
    epithelial_markers = panck_pos | ecad_high

    vim_pos = expr_df['Vimentin'] > Vim_thresh
    asma_high = expr_df['alphaSMA'] > alphaSMA_thresh

    cd31_pos = expr_df['CD31'] > CD31_thresh
    cd68_pos = expr_df['CD68'] > CD68_thresh
    cd3e_pos = expr_df['CD3E'] > CD3E_thresh

    print("\n" + "="*60)
    print("A. VIMENTIN / EPITHELIAL OVERLAP")
    print("="*60)

    n_vim_pos = vim_pos.sum()
    n_epithelial = epithelial_markers.sum()
    n_overlap = (vim_pos & epithelial_markers).sum()

    print(f"Total cells: {len(expr_df)}")
    print(f"Vimentin+ cells: {n_vim_pos} ({n_vim_pos/len(expr_df)*100:.1f}%)")
    print(f"PanCK+ or E-Cad high cells: {n_epithelial} ({n_epithelial/len(expr_df)*100:.1f}%)")
    print(f"OVERLAP (VIM+ AND epithelial+): {n_overlap} ({n_overlap/len(expr_df)*100:.1f}%)")
    print(f"  -> {n_overlap/n_vim_pos*100:.1f}% of VIM+ cells are also epithelial+")
    print(f"  -> {n_overlap/n_epithelial*100:.1f}% of epithelial+ cells are also VIM+")

    # Check VIM expression levels in epithelial cells
    vim_in_epithelial = expr_df.loc[epithelial_markers, 'Vimentin']
    vim_in_non_epithelial = expr_df.loc[~epithelial_markers, 'Vimentin']

    print(f"\nVimentin expression:")
    print(f"  In epithelial+ cells: mean={vim_in_epithelial.mean():.2f}, median={vim_in_epithelial.median():.2f}")
    print(f"  In non-epithelial cells: mean={vim_in_non_epithelial.mean():.2f}, median={vim_in_non_epithelial.median():.2f}")

    print("\n" + "="*60)
    print("B. FIBROBLAST GATE BREAKDOWN")
    print("="*60)

    # Simulate the hierarchical gating to see what falls through to each fibroblast gate
    # First, mark cells that would be claimed by earlier gates

    # Load actual cell type assignments to get which cells are "Unknown" before Gate 7/8
    from create_protein_gt import classify_cells_by_protein

    # Actually, let's just re-run the logic manually
    # Track which cells are still "available" after gates 1-6

    CD20_thresh = get_threshold(expr_df, 'CD20', 25)
    CD4_thresh = get_threshold(expr_df, 'CD4', 50)
    CD8A_thresh = get_threshold(expr_df, 'CD8A', 50)

    cell_types = pd.Series('Unknown', index=expr_df.index)

    # Gate 1: B cells
    b_cells = expr_df['CD20'] > CD20_thresh
    cell_types[b_cells] = 'B cells'

    # Gate 2: CD4+ T cells
    t_cell_base = expr_df['CD3E'] > CD3E_thresh
    cd4_pos = expr_df['CD4'] > CD4_thresh
    cd8_neg = expr_df['CD8A'] < CD8A_thresh
    cd4_tcells = t_cell_base & cd4_pos & cd8_neg & (cell_types == 'Unknown')
    cell_types[cd4_tcells] = 'CD4+ T cells'

    # Gate 3: CD8+ T cells
    cd8_pos = expr_df['CD8A'] > CD8A_thresh
    cd8_tcells = t_cell_base & cd8_pos & (cell_types == 'Unknown')
    cell_types[cd8_tcells] = 'CD8+ T cells'

    # Gate 4: Macrophages
    cd68_neg = expr_df['CD68'] < CD68_thresh
    cd3e_neg = expr_df['CD3E'] < CD3E_thresh
    macrophages = cd68_pos & cd3e_neg & (cell_types == 'Unknown')
    cell_types[macrophages] = 'Macrophages'

    # Gate 5: Endothelial
    endothelial = cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[endothelial] = 'Endothelial'

    # Gate 6: Epithelial
    epithelial = epithelial_markers & (cell_types == 'Unknown')
    cell_types[epithelial] = 'Epithelial'

    # Now count what's available for Gate 7 and 8
    available_for_fib = cell_types == 'Unknown'
    n_available = available_for_fib.sum()

    # Gate 7: αSMA high fibroblasts
    gate7_eligible = asma_high & ~cd31_pos & cd68_neg & cd3e_neg & available_for_fib
    n_gate7 = gate7_eligible.sum()

    # After Gate 7, update availability
    remaining_after_gate7 = available_for_fib & ~gate7_eligible

    # Gate 8: Vimentin+ stromal
    gate8_eligible = vim_pos & remaining_after_gate7
    n_gate8 = gate8_eligible.sum()

    total_fibroblasts = n_gate7 + n_gate8

    print(f"Cells available after Gates 1-6: {n_available} ({n_available/len(expr_df)*100:.1f}%)")
    print(f"\nFibroblast breakdown:")
    print(f"  Gate 7 (αSMA high): {n_gate7} ({n_gate7/total_fibroblasts*100:.1f}% of fibroblasts)")
    print(f"  Gate 8 (VIM+ remaining): {n_gate8} ({n_gate8/total_fibroblasts*100:.1f}% of fibroblasts)")
    print(f"  TOTAL Fibroblasts: {total_fibroblasts}")

    # Check overlap: how many Gate 8 cells have epithelial markers?
    gate8_with_epithelial = (gate8_eligible & epithelial_markers).sum()
    print(f"\nGate 8 cells that are ALSO PanCK+ or E-Cad high: {gate8_with_epithelial}")
    print(f"  -> {gate8_with_epithelial/n_gate8*100:.1f}% of Gate 8 'fibroblasts' have epithelial markers!")

    # Remaining unknown
    still_unknown = remaining_after_gate7 & ~gate8_eligible
    print(f"\nStill Unknown after all gates: {still_unknown.sum()}")

    print("\n" + "="*60)
    print("C. VIMENTIN EXPRESSION BY FINAL CELL TYPE")
    print("="*60)

    # Update cell types with fibroblast gates
    cell_types[gate7_eligible] = 'Fibroblasts (αSMA)'
    cell_types[gate8_eligible] = 'Fibroblasts (VIM)'

    for ct in cell_types.unique():
        ct_mask = cell_types == ct
        vim_vals = expr_df.loc[ct_mask, 'Vimentin']
        pct_vim_pos = (vim_vals > Vim_thresh).mean() * 100
        print(f"{ct}: {ct_mask.sum()} cells, {pct_vim_pos:.1f}% VIM+, mean VIM={vim_vals.mean():.1f}")

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"""
KEY FINDINGS:
1. VIM/Epithelial overlap: {n_overlap} cells ({n_overlap/len(expr_df)*100:.1f}% of all cells)
   are positive for BOTH vimentin AND epithelial markers (PanCK/E-Cad)

2. Fibroblast composition:
   - {n_gate7/total_fibroblasts*100:.1f}% from Gate 7 (αSMA high) - true CAFs/myofibroblasts
   - {n_gate8/total_fibroblasts*100:.1f}% from Gate 8 (VIM+ catchall) - mixed population

3. Gate 8 contamination: {gate8_with_epithelial/n_gate8*100:.1f}% of VIM+ "fibroblasts"
   also express epithelial markers - likely EMT tumor cells!

RECOMMENDATION:
Consider removing Gate 8 (Vimentin+) from fibroblast definition in RCC.
αSMA alone is a more specific marker for cancer-associated fibroblasts.
""")

if __name__ == '__main__':
    main()
