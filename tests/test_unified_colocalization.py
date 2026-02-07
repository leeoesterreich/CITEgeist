#!/usr/bin/env python
"""
Test unified spatial colocalization on pILC single-cell data.

Compares the resolution-agnostic approach vs the original approach
to demonstrate that unified colocalization works for single-cell data.
"""

import sys
import json
import logging
from pathlib import Path
from datetime import datetime

import numpy as np
import scanpy as sc

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.spatial_colocalization_unified import (
    analyze_colocalization_unified,
    build_colocalization_graph,
)
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    _find_connected_components,
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_pilc_sample(
    sample_id: str = "TP08-2202",
    max_cells: int = 10000,
) -> sc.AnnData:
    """Load pILC Visium HD sample."""

    data_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project/enact_adatas")
    adata_path = data_dir / f"{sample_id}_cells_preprocessed.h5ad"

    logger.info(f"Loading {sample_id} from {adata_path}")
    adata = sc.read_h5ad(adata_path)

    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Subsample if needed
    if adata.n_obs > max_cells:
        logger.info(f"Subsampling to {max_cells} cells")
        np.random.seed(42)
        indices = np.random.choice(adata.n_obs, max_cells, replace=False)
        adata = adata[indices].copy()

    return adata


def select_test_markers(adata: sc.AnnData, n_markers: int = 50) -> list:
    """Select test markers: curated + top HVGs."""

    # Curated breast cancer markers
    curated = [
        'ESR1', 'ERBB2', 'FOXA1', 'GATA3',  # Breast cancer
        'CD3E', 'CD4', 'CD8A', 'FOXP3',  # T cells
        'CD19', 'MS4A1', 'CD79A',  # B cells
        'CD68', 'CD163', 'CD14',  # Macrophages
        'PECAM1', 'VWF', 'CDH5',  # Endothelial
        'COL1A1', 'VIM', 'ACTA2',  # Stromal
        'EPCAM', 'KRT8', 'KRT18',  # Epithelial
        'MKI67', 'TOP2A',  # Proliferation
        'ADIPOQ', 'LEP', 'FABP4',  # Adipocyte
    ]

    # Filter to those present
    present_curated = [g for g in curated if g in adata.var_names]
    logger.info(f"Found {len(present_curated)}/{len(curated)} curated markers")

    # Add top HVGs to reach n_markers
    n_hvg_needed = n_markers - len(present_curated)
    if n_hvg_needed > 0 and 'highly_variable' in adata.var.columns:
        hvgs = adata.var_names[adata.var['highly_variable']].tolist()
        hvgs = [g for g in hvgs if g not in present_curated]
        present_curated.extend(hvgs[:n_hvg_needed])

    return present_curated[:n_markers]


def run_comparison(
    adata: sc.AnnData,
    markers: list,
    n_permutations: int = 199,
    output_dir: Path = None,
):
    """Run both original and unified colocalization, compare results."""

    logger.info(f"\n{'='*60}")
    logger.info(f"COMPARISON: {len(markers)} markers, {adata.n_obs} cells")
    logger.info(f"{'='*60}")

    # Get coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    else:
        coords = adata.obs[['x_centroid', 'y_centroid']].values

    # Get expression matrix for markers
    marker_idx = [adata.var_names.get_loc(m) for m in markers]
    X = adata.X[:, marker_idx]
    if hasattr(X, 'toarray'):
        X = X.toarray()

    results = {}

    # =========================================
    # 1. ORIGINAL APPROACH
    # =========================================
    logger.info("\n" + "="*40)
    logger.info("ORIGINAL COLOCALIZATION (bivariate Moran's I)")
    logger.info("="*40)

    original_result = analyze_marker_colocalization(
        X, coords, markers,
        neighbor_k=6,
        n_permutations=n_permutations,
        verbose=True,
    )

    # Count significant pairs
    orig_sig_moran = sum(1 for p in original_result.pairs if p.bivariate_morans_pvalue < 0.05)
    orig_sig_corr = sum(1 for p in original_result.pairs if p.correlation_pvalue < 0.05)

    logger.info(f"Original - Significant pairs (Moran p<0.05): {orig_sig_moran}/{len(original_result.pairs)}")
    logger.info(f"Original - Significant pairs (Corr p<0.05): {orig_sig_corr}/{len(original_result.pairs)}")

    # Run profile discovery
    orig_profiles = discover_profiles(original_result, verbose=True)

    results['original'] = {
        'n_pairs': len(original_result.pairs),
        'sig_moran': orig_sig_moran,
        'sig_corr': orig_sig_corr,
        'n_profiles': len(orig_profiles.profiles),
        'n_singletons': len(orig_profiles.singletons),
        'n_edges': orig_profiles.n_significant_edges,
        'profiles': orig_profiles.profiles,
    }

    # =========================================
    # 2. UNIFIED APPROACH
    # =========================================
    logger.info("\n" + "="*40)
    logger.info("UNIFIED COLOCALIZATION (spatial pattern similarity)")
    logger.info("="*40)

    unified_results = analyze_colocalization_unified(
        X, coords, markers,
        neighbor_k=15,
        n_permutations=n_permutations,
        n_jobs=-1,
        verbose=True,
    )

    # Count significant pairs
    unified_sig = sum(1 for r in unified_results if r.spatial_pattern_pvalue < 0.05)

    logger.info(f"Unified - Significant pairs (p<0.05): {unified_sig}/{len(unified_results)}")

    # Build graph
    unified_edges = build_colocalization_graph(
        unified_results, markers,
        pvalue_threshold=0.05,
        top_k=5,
        verbose=True,
    )

    # Find connected components for profile count
    if unified_edges:
        components = _find_connected_components(markers, unified_edges)
        n_components = len([c for c in components if len(c) > 1])
        n_unified_singletons = len([c for c in components if len(c) == 1])
    else:
        n_components = 0
        n_unified_singletons = len(markers)

    results['unified'] = {
        'n_pairs': len(unified_results),
        'sig_spatial': unified_sig,
        'n_edges': len(unified_edges),
        'n_components': n_components,
        'n_singletons': n_unified_singletons,
        'edges': unified_edges[:20],  # Top 20 for inspection
    }

    # =========================================
    # 3. COMPARISON SUMMARY
    # =========================================
    logger.info("\n" + "="*60)
    logger.info("COMPARISON SUMMARY")
    logger.info("="*60)

    logger.info(f"""
    Metric                  | Original | Unified
    ------------------------|----------|--------
    Total pairs             | {results['original']['n_pairs']:>8} | {results['unified']['n_pairs']:>7}
    Significant (p<0.05)    | {results['original']['sig_moran']:>8} | {results['unified']['sig_spatial']:>7}
    Edges after filtering   | {results['original']['n_edges']:>8} | {results['unified']['n_edges']:>7}
    Multi-marker profiles   | {results['original']['n_profiles'] - results['original']['n_singletons']:>8} | {results['unified']['n_components']:>7}
    Singletons              | {results['original']['n_singletons']:>8} | {results['unified']['n_singletons']:>7}
    """)

    # Success criteria
    unified_better = results['unified']['sig_spatial'] > results['original']['sig_moran']
    unified_has_edges = results['unified']['n_edges'] > 0

    if unified_better:
        logger.info("✓ Unified approach found MORE significant pairs")
    else:
        logger.info("✗ Unified approach did not improve over original")

    if unified_has_edges:
        logger.info("✓ Unified approach produced edges for graph construction")
    else:
        logger.info("✗ Unified approach produced no edges")

    # Show top unified edges
    if unified_edges:
        logger.info("\nTop unified edges:")
        for a, b, score in sorted(unified_edges, key=lambda x: x[2], reverse=True)[:10]:
            logger.info(f"  {a} -- {b}: {score:.3f}")

    # Save results
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(output_dir / "comparison_results.json", 'w') as f:
            # Convert numpy types
            def convert(obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                if isinstance(obj, np.floating):
                    return float(obj)
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                if isinstance(obj, tuple):
                    return list(obj)
                return obj

            results_json = json.loads(
                json.dumps(results, default=convert)
            )
            json.dump(results_json, f, indent=2)

        logger.info(f"\nResults saved to {output_dir}")

    return results


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Test unified colocalization")
    parser.add_argument("--sample", default="TP08-2202", help="Sample ID")
    parser.add_argument("--max-cells", type=int, default=10000, help="Max cells")
    parser.add_argument("--n-markers", type=int, default=40, help="Number of markers")
    parser.add_argument("--n-perm", type=int, default=199, help="Permutations")
    parser.add_argument("--output", type=str, default=None, help="Output directory")

    args = parser.parse_args()

    logger.info(f"Starting unified colocalization test")
    logger.info(f"Sample: {args.sample}, max_cells: {args.max_cells}, n_markers: {args.n_markers}")

    # Load data
    adata = load_pilc_sample(args.sample, max_cells=args.max_cells)

    # Select markers
    markers = select_test_markers(adata, n_markers=args.n_markers)
    logger.info(f"Selected {len(markers)} markers: {markers[:10]}...")

    # Output dir
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path("test_results/unified_colocalization") / args.sample

    # Run comparison
    results = run_comparison(
        adata, markers,
        n_permutations=args.n_perm,
        output_dir=output_dir,
    )

    return results


if __name__ == "__main__":
    main()
