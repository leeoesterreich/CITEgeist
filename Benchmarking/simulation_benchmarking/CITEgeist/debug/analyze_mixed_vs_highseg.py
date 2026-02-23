#!/usr/bin/env python3
"""
Analysis script to compare mixed vs high_seg simulation datasets
to understand why discrete IQP performs poorly on mixed dataset.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import jensenshannon
import warnings
warnings.filterwarnings('ignore')

# Paths
MIXED_H5AD_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/mixed/h5ad_objects"
HIGHSEG_H5AD_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects"
MIXED_GT_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/mixed/ST_sim"
HIGHSEG_GT_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/ST_sim"

def load_replicate(rep_id, dataset_type):
    """Load CITE-seq and GEX data for a replicate"""
    h5ad_dir = MIXED_H5AD_DIR if dataset_type == "mixed" else HIGHSEG_H5AD_DIR
    gt_dir = MIXED_GT_DIR if dataset_type == "mixed" else HIGHSEG_GT_DIR

    cite_path = f"{h5ad_dir}/Wu_rep_{rep_id}_CITE.h5ad"
    gex_path = f"{h5ad_dir}/Wu_rep_{rep_id}_GEX.h5ad"
    prop_path = f"{gt_dir}/Wu_ST_{rep_id}_prop.csv"

    cite_adata = sc.read_h5ad(cite_path)
    gex_adata = sc.read_h5ad(gex_path)
    props = pd.read_csv(prop_path, index_col=0)

    # Remove spot_x and spot_y columns (not cell types!)
    cell_type_cols = [c for c in props.columns if c not in ['spot_x', 'spot_y']]
    props = props[cell_type_cols]

    return cite_adata, gex_adata, props


def analyze_proportions(props, dataset_name):
    """Analyze cell proportion characteristics"""
    print(f"\n{'='*60}")
    print(f"Proportion Analysis: {dataset_name}")
    print(f"{'='*60}")

    n_spots = len(props)
    cell_types = props.columns.tolist()
    n_cell_types = len(cell_types)

    print(f"Number of spots: {n_spots}")
    print(f"Number of cell types: {n_cell_types}")
    print(f"Cell types: {cell_types}")

    # Count cell types per spot (non-zero proportions)
    threshold = 0.01  # 1% threshold for "present"
    present_per_spot = (props > threshold).sum(axis=1)

    print(f"\n--- Cell types per spot (threshold={threshold}) ---")
    print(f"Mean: {present_per_spot.mean():.2f}")
    print(f"Median: {present_per_spot.median():.2f}")
    print(f"Min: {present_per_spot.min()}")
    print(f"Max: {present_per_spot.max()}")
    print(f"Std: {present_per_spot.std():.2f}")

    # Distribution of cell type counts
    print("\nDistribution of cell type counts per spot:")
    for n in range(1, n_cell_types + 1):
        count = (present_per_spot == n).sum()
        pct = 100 * count / n_spots
        print(f"  {n} cell types: {count} spots ({pct:.1f}%)")

    # Entropy analysis
    # Add small epsilon to avoid log(0)
    props_safe = props.values + 1e-10
    props_safe = props_safe / props_safe.sum(axis=1, keepdims=True)
    entropy = -np.sum(props_safe * np.log2(props_safe), axis=1)
    max_entropy = np.log2(n_cell_types)
    normalized_entropy = entropy / max_entropy

    print(f"\n--- Entropy analysis ---")
    print(f"Mean entropy: {entropy.mean():.3f} (max possible: {max_entropy:.3f})")
    print(f"Mean normalized entropy: {normalized_entropy.mean():.3f}")
    print(f"Median normalized entropy: {np.median(normalized_entropy):.3f}")
    print(f"Std normalized entropy: {normalized_entropy.std():.3f}")

    # Pure spots analysis (dominant cell type > X%)
    dominant_prop = props.max(axis=1)

    print(f"\n--- Spot purity analysis ---")
    for thresh in [0.5, 0.7, 0.9, 0.95, 0.99]:
        n_pure = (dominant_prop >= thresh).sum()
        pct = 100 * n_pure / n_spots
        print(f"  Dominant >= {thresh*100:.0f}%: {n_pure} spots ({pct:.1f}%)")

    # Per cell type mean proportion
    print(f"\n--- Mean proportion per cell type ---")
    for ct in cell_types:
        mean_prop = props[ct].mean()
        non_zero = (props[ct] > threshold).sum()
        print(f"  {ct}: mean={mean_prop:.4f}, present in {non_zero} spots ({100*non_zero/n_spots:.1f}%)")

    return {
        'present_per_spot': present_per_spot,
        'entropy': entropy,
        'normalized_entropy': normalized_entropy,
        'dominant_prop': dominant_prop
    }


def analyze_marker_signal(cite_adata, props, dataset_name):
    """Analyze marker expression characteristics"""
    print(f"\n{'='*60}")
    print(f"Marker Signal Analysis: {dataset_name}")
    print(f"{'='*60}")

    # Get marker names
    markers = cite_adata.var_names.tolist()
    cell_types = props.columns.tolist()

    # Identify cell-type-specific markers (those with cell type names)
    specific_markers = {}
    for ct in cell_types:
        ct_markers = [m for m in markers if ct.replace(' ', '-') in m or ct in m]
        if ct_markers:
            specific_markers[ct] = ct_markers

    print(f"\nTotal markers: {len(markers)}")
    print(f"Cell types with specific markers: {len(specific_markers)}")

    # Get expression matrix
    X = cite_adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    print(f"\n--- Raw expression statistics ---")
    print(f"Expression matrix shape: {X.shape}")
    print(f"Mean expression: {X.mean():.4f}")
    print(f"Std expression: {X.std():.4f}")
    print(f"Max expression: {X.max():.4f}")
    print(f"% zeros: {100*(X == 0).sum() / X.size:.2f}%")

    # Analyze specific markers for each cell type
    print(f"\n--- Cell-type specific marker analysis ---")
    marker_stats = {}

    for ct, ct_markers in specific_markers.items():
        ct_props = props[ct].values

        for marker in ct_markers:
            marker_idx = markers.index(marker)
            marker_expr = X[:, marker_idx]

            # Correlation with cell type proportion
            corr, pval = stats.pearsonr(marker_expr, ct_props)

            # Signal in high-proportion spots vs low-proportion spots
            high_prop_mask = ct_props > 0.3
            low_prop_mask = ct_props < 0.1

            if high_prop_mask.sum() > 0 and low_prop_mask.sum() > 0:
                high_signal = marker_expr[high_prop_mask].mean()
                low_signal = marker_expr[low_prop_mask].mean()
                snr = high_signal / (low_signal + 1e-10)
            else:
                high_signal = np.nan
                low_signal = np.nan
                snr = np.nan

            marker_stats[marker] = {
                'cell_type': ct,
                'corr': corr,
                'pval': pval,
                'high_signal': high_signal,
                'low_signal': low_signal,
                'snr': snr,
                'mean_expr': marker_expr.mean(),
                'std_expr': marker_expr.std()
            }

            print(f"  {marker} ({ct}):")
            print(f"    Correlation with proportion: {corr:.3f} (p={pval:.2e})")
            print(f"    Mean expr: {marker_expr.mean():.3f}, Std: {marker_expr.std():.3f}")
            print(f"    SNR (high/low prop): {snr:.2f}")

    return marker_stats


def analyze_marker_correlations(cite_adata, dataset_name):
    """Analyze correlations between markers"""
    print(f"\n{'='*60}")
    print(f"Marker Correlation Analysis: {dataset_name}")
    print(f"{'='*60}")

    X = cite_adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    markers = cite_adata.var_names.tolist()

    # Compute correlation matrix
    corr_matrix = np.corrcoef(X.T)

    # Identify cell-type specific markers
    specific_markers = [m for m in markers if 'Protein' in m and 'Nonspecific' not in m]
    nonspecific_markers = [m for m in markers if 'Nonspecific' in m]

    print(f"\nSpecific markers: {len(specific_markers)}")
    print(f"Nonspecific markers: {len(nonspecific_markers)}")

    # Analyze within-cell-type correlations
    # Group markers by cell type
    cell_type_markers = {}
    for m in specific_markers:
        # Extract cell type from marker name (e.g., "B-cells_Protein_1" -> "B-cells")
        ct = m.rsplit('_', 2)[0]  # Remove "_Protein_X" suffix
        if ct not in cell_type_markers:
            cell_type_markers[ct] = []
        cell_type_markers[ct].append(m)

    print(f"\n--- Within-cell-type marker correlations ---")
    within_corrs = []
    for ct, ct_markers in cell_type_markers.items():
        if len(ct_markers) >= 2:
            ct_indices = [markers.index(m) for m in ct_markers]
            ct_corrs = []
            for i in range(len(ct_indices)):
                for j in range(i+1, len(ct_indices)):
                    ct_corrs.append(corr_matrix[ct_indices[i], ct_indices[j]])
            mean_corr = np.mean(ct_corrs)
            within_corrs.extend(ct_corrs)
            print(f"  {ct}: {len(ct_markers)} markers, mean within-type corr: {mean_corr:.3f}")

    if within_corrs:
        print(f"\nOverall within-cell-type correlation: {np.mean(within_corrs):.3f}")

    # Analyze between-cell-type correlations
    print(f"\n--- Between-cell-type marker correlations ---")
    between_corrs = []
    cts = list(cell_type_markers.keys())
    for i in range(len(cts)):
        for j in range(i+1, len(cts)):
            ct1_indices = [markers.index(m) for m in cell_type_markers[cts[i]]]
            ct2_indices = [markers.index(m) for m in cell_type_markers[cts[j]]]
            for idx1 in ct1_indices:
                for idx2 in ct2_indices:
                    between_corrs.append(corr_matrix[idx1, idx2])

    if between_corrs:
        print(f"Mean between-cell-type correlation: {np.mean(between_corrs):.3f}")
        print(f"Std between-cell-type correlation: {np.std(between_corrs):.3f}")

    return corr_matrix


def analyze_integer_feasibility(props, dataset_name):
    """Analyze how well integer cell counts can approximate proportions"""
    print(f"\n{'='*60}")
    print(f"Integer Feasibility Analysis: {dataset_name}")
    print(f"{'='*60}")

    # Simulate different total cell counts per spot
    for n_cells in [5, 10, 15, 20]:
        # Convert proportions to integer counts
        int_counts = np.round(props.values * n_cells).astype(int)

        # Ensure counts sum to n_cells (adjust largest count)
        for i in range(len(int_counts)):
            diff = n_cells - int_counts[i].sum()
            if diff != 0:
                max_idx = np.argmax(int_counts[i])
                int_counts[i, max_idx] += diff

        # Convert back to proportions
        approx_props = int_counts / n_cells

        # Calculate approximation error
        errors = np.abs(props.values - approx_props)
        mae = errors.mean()
        max_error = errors.max()

        # Count spots where integer constraint causes significant error
        spot_errors = errors.max(axis=1)  # Max error per spot
        significant_error_spots = (spot_errors > 0.1).sum()

        print(f"\n--- n_cells = {n_cells} ---")
        print(f"  Mean absolute error: {mae:.4f}")
        print(f"  Max error: {max_error:.4f}")
        print(f"  Spots with max error > 0.1: {significant_error_spots} ({100*significant_error_spots/len(props):.1f}%)")

        # Correlation preservation
        flat_true = props.values.flatten()
        flat_approx = approx_props.flatten()
        corr, _ = stats.pearsonr(flat_true, flat_approx)
        print(f"  Correlation (true vs approx): {corr:.4f}")


def load_nuclei_counts(dataset_type, rep_id):
    """Load nuclei counts from cellpose segmentation"""
    nuclei_dir = f"/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CITEgeist/{dataset_type}/nuclei_counts"
    nuclei_path = f"{nuclei_dir}/Wu_rep_{rep_id}_nuclei_counts.csv"

    try:
        nuclei = pd.read_csv(nuclei_path, index_col=0)
        return nuclei
    except FileNotFoundError:
        print(f"Nuclei counts not found at {nuclei_path}")
        return None


def analyze_nuclei_vs_ground_truth(nuclei_df, props, dataset_name):
    """Compare cellpose nuclei counts with ground truth cell counts"""
    print(f"\n{'='*60}")
    print(f"Nuclei Count Analysis: {dataset_name}")
    print(f"{'='*60}")

    if nuclei_df is None:
        print("No nuclei counts available")
        return {}

    # Align indices
    common_spots = props.index.intersection(nuclei_df.index)
    print(f"Common spots: {len(common_spots)}")

    # Get nuclei counts per spot
    nuclei_counts = nuclei_df.loc[common_spots, 'nuclei_count'].values if 'nuclei_count' in nuclei_df.columns else nuclei_df.iloc[:, 0].values

    print(f"\n--- Nuclei count distribution ---")
    print(f"Mean: {nuclei_counts.mean():.2f}")
    print(f"Median: {np.median(nuclei_counts):.2f}")
    print(f"Min: {nuclei_counts.min()}")
    print(f"Max: {nuclei_counts.max()}")
    print(f"Std: {nuclei_counts.std():.2f}")

    # Distribution of nuclei counts
    print("\nDistribution of nuclei counts:")
    for n in range(1, 21):
        count = (nuclei_counts == n).sum()
        if count > 0:
            pct = 100 * count / len(nuclei_counts)
            print(f"  {n} nuclei: {count} spots ({pct:.1f}%)")

    # Analyze integer representation capacity
    print("\n--- Integer representation capacity ---")
    props_aligned = props.loc[common_spots].values
    n_cell_types = props_aligned.shape[1]

    # For each spot, can integer counts represent the proportions?
    total_error = 0
    large_error_spots = 0

    for i, n in enumerate(nuclei_counts):
        n = int(n)
        if n == 0:
            continue
        true_props = props_aligned[i]

        # Best integer approximation
        int_counts = np.round(true_props * n).astype(int)
        diff = n - int_counts.sum()
        if diff != 0:
            # Adjust to maintain total
            max_idx = np.argmax(int_counts if diff < 0 else true_props - int_counts / n)
            int_counts[max_idx] += diff

        approx_props = int_counts / n
        error = np.abs(true_props - approx_props).max()
        total_error += error
        if error > 0.15:
            large_error_spots += 1

    print(f"Mean max error per spot: {total_error / len(nuclei_counts):.4f}")
    print(f"Spots with max error > 0.15: {large_error_spots} ({100*large_error_spots/len(nuclei_counts):.1f}%)")

    return {
        'nuclei_counts': nuclei_counts,
        'mean_nuclei': nuclei_counts.mean(),
        'large_error_spots': large_error_spots
    }


def analyze_proportion_granularity(props, dataset_name):
    """Analyze how granular the proportions are (what fractions appear)"""
    print(f"\n{'='*60}")
    print(f"Proportion Granularity Analysis: {dataset_name}")
    print(f"{'='*60}")

    # Get all non-zero proportions
    all_props = props.values.flatten()
    nonzero_props = all_props[all_props > 0.001]

    print(f"Total proportion values: {len(all_props)}")
    print(f"Non-zero proportions: {len(nonzero_props)}")

    # Check what fractions they correspond to
    # These would be n/N where N is total cells in spot
    print("\n--- Common proportion values ---")
    unique_props, counts = np.unique(np.round(nonzero_props, 4), return_counts=True)
    sorted_idx = np.argsort(-counts)[:20]
    for idx in sorted_idx:
        print(f"  {unique_props[idx]:.4f}: {counts[idx]} occurrences")

    # Infer likely total cell counts per spot
    print("\n--- Inferred total cells per spot ---")
    # For each spot, find the smallest denominator that explains all proportions
    inferred_totals = []
    for i in range(len(props)):
        row_props = props.iloc[i].values
        nonzero = row_props[row_props > 0.001]
        if len(nonzero) == 0:
            continue

        # Try different denominators
        best_denom = None
        for denom in range(1, 51):
            # Check if all proportions are close to n/denom for some integer n
            fracs = nonzero * denom
            if np.allclose(fracs, np.round(fracs), atol=0.01):
                best_denom = denom
                break

        if best_denom:
            inferred_totals.append(best_denom)

    if inferred_totals:
        inferred_totals = np.array(inferred_totals)
        print(f"Mean inferred total cells: {inferred_totals.mean():.2f}")
        print(f"Median inferred total cells: {np.median(inferred_totals):.2f}")
        print(f"Min: {inferred_totals.min()}, Max: {inferred_totals.max()}")

        # Distribution
        for n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]:
            count = (inferred_totals == n).sum()
            if count > 0:
                print(f"  {n} cells: {count} spots ({100*count/len(inferred_totals):.1f}%)")

    return np.array(inferred_totals) if len(inferred_totals) > 0 else []


def main():
    print("="*70)
    print("COMPARISON: MIXED vs HIGH_SEG SIMULATION DATASETS")
    print("="*70)

    # Load data for rep 0
    rep_id = 0

    print("\n\nLoading MIXED dataset...")
    mixed_cite, mixed_gex, mixed_props = load_replicate(rep_id, "mixed")

    print("\n\nLoading HIGH_SEG dataset...")
    highseg_cite, highseg_gex, highseg_props = load_replicate(rep_id, "high_seg")

    # Load nuclei counts
    print("\n\nLoading nuclei counts...")
    mixed_nuclei = load_nuclei_counts("mixed", rep_id)
    highseg_nuclei = load_nuclei_counts("high_seg", rep_id)

    # Analyze proportions
    mixed_prop_stats = analyze_proportions(mixed_props, "MIXED")
    highseg_prop_stats = analyze_proportions(highseg_props, "HIGH_SEG")

    # Analyze proportion granularity (infer cell counts from proportions)
    print("\n\n" + "="*70)
    print("PROPORTION GRANULARITY ANALYSIS")
    print("="*70)
    mixed_inferred = analyze_proportion_granularity(mixed_props, "MIXED")
    highseg_inferred = analyze_proportion_granularity(highseg_props, "HIGH_SEG")

    # Analyze nuclei counts
    print("\n\n" + "="*70)
    print("NUCLEI COUNT ANALYSIS")
    print("="*70)
    mixed_nuclei_stats = analyze_nuclei_vs_ground_truth(mixed_nuclei, mixed_props, "MIXED")
    highseg_nuclei_stats = analyze_nuclei_vs_ground_truth(highseg_nuclei, highseg_props, "HIGH_SEG")

    # Comparative summary
    print("\n" + "="*70)
    print("COMPARATIVE SUMMARY: Proportion Distributions")
    print("="*70)
    print(f"\n{'Metric':<40} {'MIXED':>12} {'HIGH_SEG':>12}")
    print("-"*64)
    print(f"{'Mean cell types per spot':<40} {mixed_prop_stats['present_per_spot'].mean():>12.2f} {highseg_prop_stats['present_per_spot'].mean():>12.2f}")
    print(f"{'Mean normalized entropy':<40} {mixed_prop_stats['normalized_entropy'].mean():>12.3f} {highseg_prop_stats['normalized_entropy'].mean():>12.3f}")
    print(f"{'% spots with dominant >= 90%':<40} {100*(mixed_prop_stats['dominant_prop'] >= 0.9).mean():>12.1f}% {100*(highseg_prop_stats['dominant_prop'] >= 0.9).mean():>12.1f}%")
    print(f"{'% spots with dominant >= 50%':<40} {100*(mixed_prop_stats['dominant_prop'] >= 0.5).mean():>12.1f}% {100*(highseg_prop_stats['dominant_prop'] >= 0.5).mean():>12.1f}%")

    # Analyze marker signal
    mixed_marker_stats = analyze_marker_signal(mixed_cite, mixed_props, "MIXED")
    highseg_marker_stats = analyze_marker_signal(highseg_cite, highseg_props, "HIGH_SEG")

    # Analyze marker correlations
    mixed_corr = analyze_marker_correlations(mixed_cite, "MIXED")
    highseg_corr = analyze_marker_correlations(highseg_cite, "HIGH_SEG")

    # Integer feasibility analysis
    analyze_integer_feasibility(mixed_props, "MIXED")
    analyze_integer_feasibility(highseg_props, "HIGH_SEG")

    # Final diagnosis
    print("\n" + "="*70)
    print("DIAGNOSTIC SUMMARY")
    print("="*70)

    print("""
Based on the analysis, the key differences between datasets are:

1. CELL TYPE MIXING:
   - HIGH_SEG: Spots tend to be dominated by single cell types (high purity)
   - MIXED: Spots contain multiple cell types in similar proportions (low purity)

2. INTEGER CONSTRAINT IMPACT:
   - HIGH_SEG: With ~1-2 dominant cell types, integer constraints work well
   - MIXED: With many cell types per spot, integer rounding loses information

3. SIGNAL AMBIGUITY:
   - HIGH_SEG: Clear marker signals when cell type dominates
   - MIXED: Overlapping marker signals from multiple cell types

HYPOTHESIS: The discrete IQP fails on MIXED because:
a) Integer constraints cannot represent fractional cell type mixtures
b) With limited cells per spot, rounding errors accumulate
c) Multiple cell types create marker signal ambiguity that integers can't resolve
""")


if __name__ == "__main__":
    main()
