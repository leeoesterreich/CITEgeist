#!/usr/bin/env python
"""
Test unified colocalization on BOTH single-cell and spot-level data.

Compares scoring methods:
1. "weighted" (legacy): 0.5*spatial + 0.3*neigh + 0.2*expr
2. "geometric": sqrt(spatial * expr) - requires both to be good
3. "min": min(spatial, expr) - bottleneck by weaker
4. "harmonic": 2*s*e/(s+e) - strongly penalizes low values

Tests on:
- Single-cell: pILC Visium HD (TP08-2202)
- Spot-level: Xenium pseudo-Visium region
"""

import sys
import json
import logging
from pathlib import Path
from collections import defaultdict

import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.spatial_colocalization_unified import (
    analyze_colocalization_unified,
    build_colocalization_graph,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Curated markers for breast cancer analysis
CURATED_MARKERS = [
    # Luminal/Epithelial
    'ESR1', 'FOXA1', 'GATA3', 'EPCAM', 'KRT8', 'KRT18',
    # HER2
    'ERBB2',
    # T cells
    'CD3E', 'CD4', 'CD8A', 'FOXP3',
    # B cells
    'CD19', 'MS4A1', 'CD79A',
    # Macrophages
    'CD68', 'CD163', 'CD14',
    # Endothelial
    'PECAM1', 'VWF', 'CDH5',
    # Stromal
    'COL1A1', 'VIM', 'ACTA2',
    # Proliferation
    'MKI67', 'TOP2A',
]

# Expected biological groupings for validation (pILC single-cell)
EXPECTED_PROFILES = {
    'Luminal': {'ESR1', 'FOXA1', 'GATA3', 'EPCAM', 'KRT8', 'KRT18'},
    'T_cell': {'CD3E', 'CD4', 'CD8A'},
    'B_cell': {'CD19', 'MS4A1', 'CD79A'},
    'Macrophage': {'CD68', 'CD163', 'CD14'},
    'Endothelial': {'PECAM1', 'VWF', 'CDH5'},
    'Stromal': {'COL1A1', 'VIM', 'ACTA2'},
}

# Expected profiles for Xenium pseudo-Visium (different marker panel)
# Markers: PD-1, VISTA, PD-L1, LAG-3, CD16, GranzymeB, CD163, CD4, CD20, CD8A,
#          CD3E, CD138, HLA-DR, CD11c, CD68, CD45RA, PCNA, CD45RO, Ki-67,
#          Beta-catenin, CD31, PTEN, PanCK, Vimentin, alphaSMA, CD45, E-Cadherin
XENIUM_EXPECTED_PROFILES = {
    'T_cell': {'CD3E', 'CD4', 'CD8A', 'CD45'},
    'T_cell_activated': {'GranzymeB', 'CD45RO'},
    'T_cell_naive': {'CD45RA'},
    'B_cell': {'CD20', 'CD138'},  # CD20=B cells, CD138=plasma cells
    'Macrophage': {'CD68', 'CD163', 'CD11c'},
    'Immune_checkpoint': {'PD-1', 'PD-L1', 'LAG-3', 'VISTA'},
    'Epithelial': {'PanCK', 'E-Cadherin', 'Beta-catenin'},
    'Stromal': {'Vimentin', 'alphaSMA'},
    'Endothelial': {'CD31'},
    'Proliferation': {'Ki-67', 'PCNA'},
}


def load_single_cell_data(max_cells: int = 10000):
    """Load pILC single-cell data."""
    data_path = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project/enact_adatas/TP08-2202_cells_preprocessed.h5ad")

    logger.info(f"Loading single-cell data from {data_path}")
    adata = sc.read_h5ad(data_path)

    if adata.n_obs > max_cells:
        np.random.seed(42)
        indices = np.random.choice(adata.n_obs, max_cells, replace=False)
        adata = adata[indices].copy()

    logger.info(f"Single-cell: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata, "single_cell"


def load_spot_level_data(max_spots: int = 5000):
    """Load Xenium pseudo-Visium spot-level data with proper spatial coordinates."""
    # Use the Xenium pseudo-Visium CITE data which has actual spatial coords
    base_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt/h5ad_objects")

    if not base_dir.exists():
        logger.warning(f"Pseudo-Visium directory not found: {base_dir}")
        return None, None

    # Combine multiple regions for larger sample
    adata_list = []
    for region_id in range(5):  # 5 Xenium regions
        cite_path = base_dir / f"Xenium_region_{region_id}_CITE.h5ad"
        if cite_path.exists():
            adata_region = sc.read_h5ad(cite_path)
            adata_region.obs['region'] = f"region_{region_id}"
            adata_list.append(adata_region)
            logger.info(f"Loaded region {region_id}: {adata_region.n_obs} spots")

    if not adata_list:
        logger.warning("No Xenium CITE h5ad files found")
        return None, None

    # Concatenate regions
    import anndata as ad
    adata = ad.concat(adata_list, join='outer', label='region')

    # Ensure spatial coordinates exist
    if 'spatial' not in adata.obsm:
        logger.warning("No spatial coordinates in Xenium data")
        return None, None

    if adata.n_obs > max_spots:
        np.random.seed(42)
        indices = np.random.choice(adata.n_obs, max_spots, replace=False)
        adata = adata[indices].copy()

    logger.info(f"Spot-level: {adata.n_obs} spots, {adata.n_vars} markers")
    return adata, "spot_level"


def get_available_markers(adata, curated_list):
    """Get markers available in the dataset."""
    available = [m for m in curated_list if m in adata.var_names]
    logger.info(f"Found {len(available)}/{len(curated_list)} curated markers")
    return available


def extract_profiles_from_edges(edges, markers):
    """Extract connected components as profiles from edge list."""
    if not edges:
        return [[m] for m in markers]  # All singletons

    # Build adjacency
    adj = defaultdict(set)
    markers_with_edges = set()
    for a, b, _ in edges:
        adj[a].add(b)
        adj[b].add(a)
        markers_with_edges.add(a)
        markers_with_edges.add(b)

    # BFS for connected components
    visited = set()
    components = []

    for start in markers_with_edges:
        if start in visited:
            continue
        component = []
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            component.append(node)
            for neighbor in adj[node]:
                if neighbor not in visited:
                    queue.append(neighbor)
        components.append(sorted(component))

    # Add singletons
    for m in markers:
        if m not in visited:
            components.append([m])

    return sorted(components, key=len, reverse=True)


def evaluate_profile_quality(profiles, expected):
    """Evaluate how well discovered profiles match expected groupings."""

    scores = {}

    for expected_name, expected_markers in expected.items():
        best_jaccard = 0
        best_profile = None

        for profile in profiles:
            profile_set = set(profile)
            intersection = len(expected_markers & profile_set)
            union = len(expected_markers | profile_set)
            jaccard = intersection / union if union > 0 else 0

            if jaccard > best_jaccard:
                best_jaccard = jaccard
                best_profile = profile

        scores[expected_name] = {
            'jaccard': best_jaccard,
            'best_match': best_profile,
            'expected': list(expected_markers),
        }

    # Average Jaccard across expected profiles
    avg_jaccard = np.mean([s['jaccard'] for s in scores.values()])

    return scores, avg_jaccard


def run_comparison(adata, resolution_name, markers, output_dir, expected_profiles=None):
    """Run comparison of scoring methods on one dataset."""
    if expected_profiles is None:
        expected_profiles = EXPECTED_PROFILES

    logger.info(f"\n{'='*60}")
    logger.info(f"TESTING: {resolution_name} ({adata.n_obs} locations, {len(markers)} markers)")
    logger.info(f"{'='*60}")

    # Get coordinates and expression
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'X_spatial' in adata.obsm:
        coords = adata.obsm['X_spatial']
    elif 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
        coords = adata.obs[['x_centroid', 'y_centroid']].values
    elif 'array_row' in adata.obs.columns and 'array_col' in adata.obs.columns:
        coords = adata.obs[['array_row', 'array_col']].values
    else:
        # Try to find any coordinate columns
        coord_cols = [c for c in adata.obs.columns if any(x in c.lower() for x in ['x', 'y', 'row', 'col', 'coord'])]
        if len(coord_cols) >= 2:
            coords = adata.obs[coord_cols[:2]].values
        else:
            raise ValueError(f"Cannot find spatial coordinates in adata. Available obs columns: {list(adata.obs.columns)[:10]}")

    marker_idx = [adata.var_names.get_loc(m) for m in markers]
    X = adata.X[:, marker_idx]
    if hasattr(X, 'toarray'):
        X = X.toarray()

    results = {}
    # Include new methods: neighborhood handles sparse single-cell dropout
    scoring_methods = ["weighted", "ratio", "neighborhood", "spatial_neighborhood"]

    for method in scoring_methods:
        logger.info(f"\n--- Scoring method: {method} ---")

        # Run unified colocalization
        coloc_results = analyze_colocalization_unified(
            X, coords, markers,
            neighbor_k=15,
            n_permutations=199,
            n_jobs=-1,
            scoring_method=method,
            verbose=True,
        )

        # Build graph
        edges = build_colocalization_graph(
            coloc_results, markers,
            pvalue_threshold=0.05,
            top_k=5,
            verbose=True,
        )

        # Extract profiles
        profiles = extract_profiles_from_edges(edges, markers)

        # Evaluate against expected (passed as parameter)
        profile_scores, avg_jaccard = evaluate_profile_quality(profiles, expected_profiles)

        # Count statistics
        n_multi = len([p for p in profiles if len(p) > 1])
        n_singletons = len([p for p in profiles if len(p) == 1])

        # Store results
        results[method] = {
            'n_edges': len(edges),
            'n_profiles': len(profiles),
            'n_multi_marker': n_multi,
            'n_singletons': n_singletons,
            'avg_jaccard': avg_jaccard,
            'profile_scores': profile_scores,
            'profiles': profiles[:10],  # Top 10 for inspection
            'top_edges': edges[:15],
        }

        logger.info(f"  Edges: {len(edges)}, Profiles: {n_multi} multi + {n_singletons} single")
        logger.info(f"  Avg Jaccard vs expected: {avg_jaccard:.3f}")

        # Show profile matches
        for name, score in profile_scores.items():
            if score['jaccard'] > 0:
                logger.info(f"    {name}: Jaccard={score['jaccard']:.2f} -> {score['best_match']}")

    # Summary comparison
    logger.info(f"\n{'='*60}")
    logger.info(f"SUMMARY: {resolution_name}")
    logger.info(f"{'='*60}")

    header = f"{'Method':<12} | {'Edges':>6} | {'Multi':>5} | {'Single':>6} | {'AvgJaccard':>10}"
    logger.info(header)
    logger.info("-" * len(header))

    for method in scoring_methods:
        r = results[method]
        logger.info(
            f"{method:<12} | {r['n_edges']:>6} | {r['n_multi_marker']:>5} | "
            f"{r['n_singletons']:>6} | {r['avg_jaccard']:>10.3f}"
        )

    # Best method
    best_method = max(scoring_methods, key=lambda m: results[m]['avg_jaccard'])
    logger.info(f"\nBest method for {resolution_name}: {best_method} (Jaccard={results[best_method]['avg_jaccard']:.3f})")

    # Show best profiles
    logger.info(f"\nProfiles from best method ({best_method}):")
    for i, profile in enumerate(results[best_method]['profiles'][:5]):
        logger.info(f"  Profile {i+1}: {profile}")

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--max-cells", type=int, default=10000)
    parser.add_argument("--max-spots", type=int, default=5000)
    parser.add_argument("--output", type=str, default="test_results/unified_both_resolutions")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = {}

    # Test 1: Single-cell data
    logger.info("\n" + "="*70)
    logger.info("PART 1: SINGLE-CELL DATA (pILC Visium HD)")
    logger.info("="*70)

    adata_sc, res_name = load_single_cell_data(max_cells=args.max_cells)
    if adata_sc is not None:
        markers_sc = get_available_markers(adata_sc, CURATED_MARKERS)
        if len(markers_sc) >= 10:
            results_sc = run_comparison(
                adata_sc, "Single-Cell", markers_sc, output_dir,
                expected_profiles=EXPECTED_PROFILES
            )
            all_results['single_cell'] = results_sc

    # Test 2: Spot-level data
    logger.info("\n" + "="*70)
    logger.info("PART 2: SPOT-LEVEL DATA (Xenium pseudo-Visium)")
    logger.info("="*70)

    adata_spot, res_name = load_spot_level_data(max_spots=args.max_spots)
    if adata_spot is not None:
        # Use all available Xenium markers (different panel than pILC)
        markers_spot = list(adata_spot.var_names)
        logger.info(f"Using all {len(markers_spot)} Xenium markers")
        if len(markers_spot) >= 10:
            results_spot = run_comparison(
                adata_spot, "Spot-Level", markers_spot, output_dir,
                expected_profiles=XENIUM_EXPECTED_PROFILES
            )
            all_results['spot_level'] = results_spot

    # Cross-resolution comparison
    if 'single_cell' in all_results and 'spot_level' in all_results:
        logger.info("\n" + "="*70)
        logger.info("CROSS-RESOLUTION COMPARISON")
        logger.info("="*70)

        for method in ["weighted", "ratio", "neighborhood", "spatial_neighborhood"]:
            sc_jaccard = all_results['single_cell'][method]['avg_jaccard']
            spot_jaccard = all_results['spot_level'][method]['avg_jaccard']
            avg = (sc_jaccard + spot_jaccard) / 2
            logger.info(f"{method:<16}: SC={sc_jaccard:.3f}, Spot={spot_jaccard:.3f}, Avg={avg:.3f}")

        # Find best unified method
        best_unified = max(
            ["weighted", "ratio", "neighborhood", "spatial_neighborhood"],
            key=lambda m: (
                all_results['single_cell'][m]['avg_jaccard'] +
                all_results['spot_level'][m]['avg_jaccard']
            ) / 2
        )
        logger.info(f"\nBest unified method across resolutions: {best_unified}")

    # Save results
    def convert(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, tuple):
            return list(obj)
        if isinstance(obj, set):
            return list(obj)
        return obj

    with open(output_dir / "comparison_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=convert)

    logger.info(f"\nResults saved to {output_dir}")

    return all_results


if __name__ == "__main__":
    main()
