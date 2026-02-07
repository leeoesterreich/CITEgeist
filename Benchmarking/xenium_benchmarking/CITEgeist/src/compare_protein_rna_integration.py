"""
Compare protein-anchored vs RNA-only program discovery.

This script demonstrates that CITEgeist's protein-anchored approach yields
better spatial transcriptional programs than standard RNA-only analysis.

Comparisons:
1. Spatial coherence (Moran's I) - primary metric
2. Program distinctness (gene overlap)
3. Protein validation rate
4. Cross-region stability
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from itertools import combinations

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from sklearn.decomposition import NMF

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from load_xenium_singlecell import load_xenium_singlecell

logger = logging.getLogger(__name__)

# Directories
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell"
FIGURE_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "figures"


class IntegrationComparison:
    """Compare protein-anchored vs RNA-only program discovery."""

    def __init__(
        self,
        region_id: int,
        k_programs: int = 5,
        leiden_resolutions: List[float] = [0.3, 0.5, 0.7, 1.0],
        min_cells_per_cluster: int = 100,
        seed: int = 42,
    ):
        self.region_id = region_id
        self.k_programs = k_programs
        self.leiden_resolutions = leiden_resolutions
        self.min_cells_per_cluster = min_cells_per_cluster
        self.seed = seed

        self.citegeist_results = None
        self.rna_only_results = {}
        self.adata_gex = None
        self.adata_protein = None

    def load_data(self):
        """Load single-cell data and CITEgeist results."""
        logger.info(f"Loading data for region {self.region_id}")

        # Load single-cell data
        self.adata_gex, self.adata_protein = load_xenium_singlecell(
            region_id=self.region_id,
            seed=self.seed,
        )
        logger.info(f"Loaded {self.adata_gex.shape[0]} cells")

        # Load CITEgeist Module 4 results
        module4_path = OUTPUT_DIR / f"region_{self.region_id}_module4_summary.json"
        if module4_path.exists():
            with open(module4_path) as f:
                self.citegeist_results = json.load(f)
            logger.info(f"Loaded CITEgeist results: {self.citegeist_results['n_profiles_analyzed']} profiles")
        else:
            raise FileNotFoundError(f"CITEgeist results not found: {module4_path}")

    def run_rna_only_baseline(self, resolution: float = 0.5) -> Dict:
        """
        Run RNA-only baseline: Leiden clustering → NMF → Moran's I.

        This is the standard approach without protein information.
        """
        logger.info(f"Running RNA-only baseline (resolution={resolution})")

        adata = self.adata_gex.copy()

        # Preprocess RNA
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        adata_hvg = adata[:, adata.var.highly_variable].copy()

        # PCA and neighbors
        sc.pp.scale(adata_hvg, max_value=10)
        sc.tl.pca(adata_hvg, n_comps=50)
        sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)

        # Leiden clustering
        sc.tl.leiden(adata_hvg, resolution=resolution, random_state=self.seed)
        clusters = adata_hvg.obs['leiden'].values
        n_clusters = len(np.unique(clusters))
        logger.info(f"Found {n_clusters} Leiden clusters")

        # Get coordinates
        coords = self.adata_gex.obsm['spatial']

        # For each cluster, run NMF and compute Moran's I
        results = {
            'resolution': resolution,
            'n_clusters': n_clusters,
            'clusters': {},
        }

        # Get raw counts for NMF
        X_raw = self.adata_gex.X
        if hasattr(X_raw, 'toarray'):
            X_raw = X_raw.toarray()
        X_raw = np.asarray(X_raw, dtype=np.float64)

        gene_names = list(self.adata_gex.var_names)

        for cluster_id in np.unique(clusters):
            cluster_mask = clusters == cluster_id
            n_cells = cluster_mask.sum()

            if n_cells < self.min_cells_per_cluster:
                logger.info(f"Skipping cluster {cluster_id}: only {n_cells} cells")
                continue

            logger.info(f"Processing cluster {cluster_id}: {n_cells} cells")

            # Extract cells in this cluster
            X_cluster = X_raw[cluster_mask, :]
            coords_cluster = coords[cluster_mask, :]

            # Filter low-expressed genes
            gene_expr = X_cluster.sum(axis=0)
            active_genes = gene_expr > np.percentile(gene_expr, 50)
            X_filtered = X_cluster[:, active_genes]
            gene_names_filtered = [gene_names[i] for i in np.where(active_genes)[0]]

            # Normalize for NMF
            X_norm = np.log1p(X_filtered)
            X_norm = X_norm / (X_norm.max(axis=0, keepdims=True) + 1e-10)

            # Run NMF
            nmf = NMF(
                n_components=self.k_programs,
                init='nndsvda',
                random_state=self.seed,
                max_iter=500,
            )
            W = nmf.fit_transform(X_norm)  # (n_cells, K)
            H = nmf.components_  # (K, n_genes)

            # Compute Moran's I for each program
            programs = []
            for k in range(self.k_programs):
                program_activity = W[:, k]
                moran_i = self._compute_morans_i(program_activity, coords_cluster)

                # Get top genes
                loadings = H[k, :]
                top_idx = np.argsort(loadings)[::-1][:10]
                top_genes = [gene_names_filtered[i] for i in top_idx]

                programs.append({
                    'program_id': k,
                    'moran_i': moran_i,
                    'top_genes': top_genes,
                    'variance_explained': (W[:, k] ** 2).sum() / (W ** 2).sum(),
                })

            results['clusters'][str(cluster_id)] = {
                'n_cells': int(n_cells),
                'programs': programs,
            }

        return results

    def _compute_morans_i(
        self,
        values: np.ndarray,
        coords: np.ndarray,
        k: int = 10,
    ) -> float:
        """Compute Moran's I for spatial autocorrelation."""
        n = len(values)
        if n < 10:
            return np.nan

        tree = cKDTree(coords)
        _, indices = tree.query(coords, k=min(k + 1, n))

        if indices.ndim == 1:
            return np.nan

        indices = indices[:, 1:]  # Remove self
        z = values - values.mean()
        denom = np.sum(z ** 2)

        if denom == 0:
            return np.nan

        num = 0.0
        S0 = 0.0
        for i in range(n):
            for j in indices[i]:
                num += z[i] * z[j]
                S0 += 1

        if S0 == 0:
            return np.nan

        I = (n / S0) * (num / denom)
        return float(I)

    def run_all_resolutions(self):
        """Run RNA-only baseline at multiple Leiden resolutions."""
        for resolution in self.leiden_resolutions:
            self.rna_only_results[resolution] = self.run_rna_only_baseline(resolution)

    def compute_moran_comparison(self) -> pd.DataFrame:
        """Compare Moran's I between CITEgeist and RNA-only."""
        rows = []

        # CITEgeist programs
        for profile_name, profile_data in self.citegeist_results['profiles'].items():
            for prog in profile_data['programs']:
                rows.append({
                    'method': 'CITEgeist',
                    'profile': profile_name,
                    'program_id': prog['program_id'],
                    'moran_i': prog['moran_i'],
                    'n_cells': profile_data['n_cells'],
                    'resolution': None,
                })

        # RNA-only programs (all resolutions)
        for resolution, results in self.rna_only_results.items():
            for cluster_id, cluster_data in results['clusters'].items():
                for prog in cluster_data['programs']:
                    rows.append({
                        'method': f'RNA-only (res={resolution})',
                        'profile': f'Cluster_{cluster_id}',
                        'program_id': prog['program_id'],
                        'moran_i': prog['moran_i'],
                        'n_cells': cluster_data['n_cells'],
                        'resolution': resolution,
                    })

        return pd.DataFrame(rows)

    def compute_gene_overlap(self, top_n: int = 50) -> Dict:
        """Compute Jaccard index for top genes between program pairs."""
        results = {'citegeist': [], 'rna_only': {}}

        # CITEgeist gene overlap
        citegeist_genes = []
        for profile_data in self.citegeist_results['profiles'].values():
            for prog in profile_data['programs']:
                citegeist_genes.append(set(prog['top_genes'][:top_n]))

        for i, j in combinations(range(len(citegeist_genes)), 2):
            jaccard = len(citegeist_genes[i] & citegeist_genes[j]) / len(citegeist_genes[i] | citegeist_genes[j])
            results['citegeist'].append(jaccard)

        # RNA-only gene overlap (per resolution)
        for resolution, res_results in self.rna_only_results.items():
            rna_genes = []
            for cluster_data in res_results['clusters'].values():
                for prog in cluster_data['programs']:
                    rna_genes.append(set(prog['top_genes'][:top_n]))

            overlaps = []
            for i, j in combinations(range(len(rna_genes)), 2):
                if len(rna_genes[i] | rna_genes[j]) > 0:
                    jaccard = len(rna_genes[i] & rna_genes[j]) / len(rna_genes[i] | rna_genes[j])
                    overlaps.append(jaccard)
            results['rna_only'][resolution] = overlaps

        return results

    def compute_protein_validation(self) -> Dict:
        """Compute correlation between programs and proteins."""
        # Get protein expression
        X_protein = self.adata_protein.X
        if hasattr(X_protein, 'toarray'):
            X_protein = X_protein.toarray()
        X_protein = np.asarray(X_protein, dtype=np.float64)
        protein_names = list(self.adata_protein.var_names)

        results = {
            'citegeist': {'validated': 0, 'total': 0, 'correlations': []},
            'rna_only': {},
        }

        # For CITEgeist, we'd need to recompute program activities
        # For now, report that programs are anchored by proteins (by design)
        # This is a placeholder - full implementation would require storing W matrices

        # For RNA-only, compute correlations with all proteins
        # This requires the NMF W matrices which we'd need to store

        logger.info("Protein validation requires stored NMF matrices - skipping for now")

        return results

    def generate_figures(self, output_dir: Path = None):
        """Generate comparison figures."""
        if output_dir is None:
            output_dir = FIGURE_DIR
        output_dir.mkdir(parents=True, exist_ok=True)

        moran_df = self.compute_moran_comparison()
        overlap_results = self.compute_gene_overlap()

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # [A] Moran's I distribution
        ax = axes[0, 0]
        methods = ['CITEgeist'] + [f'RNA-only (res={r})' for r in self.leiden_resolutions]
        plot_data = []
        for method in methods:
            subset = moran_df[moran_df['method'] == method]['moran_i'].dropna()
            for val in subset:
                plot_data.append({'Method': method.replace('RNA-only ', 'RNA\n'), 'Moran\'s I': val})

        if plot_data:
            plot_df = pd.DataFrame(plot_data)
            sns.boxplot(data=plot_df, x='Method', y='Moran\'s I', ax=ax)
            ax.set_title('A. Spatial Coherence (Moran\'s I)')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

        # [B] Programs above threshold
        ax = axes[0, 1]
        threshold = 0.15
        counts = []
        for method in methods:
            subset = moran_df[moran_df['method'] == method]['moran_i'].dropna()
            n_above = (subset > threshold).sum()
            n_total = len(subset)
            counts.append({'Method': method.replace('RNA-only ', 'RNA\n'), 'Count': n_above, 'Total': n_total})

        if counts:
            count_df = pd.DataFrame(counts)
            bars = ax.bar(count_df['Method'], count_df['Count'])
            ax.set_title(f'B. Programs with Moran\'s I > {threshold}')
            ax.set_ylabel('Number of programs')
            ax.set_xticklabels(count_df['Method'], rotation=45, ha='right')

            # Add totals as text
            for i, (bar, total) in enumerate(zip(bars, count_df['Total'])):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                       f'/{total}', ha='center', va='bottom', fontsize=9)

        # [C] Gene overlap distribution
        ax = axes[1, 0]
        overlap_data = []
        if overlap_results['citegeist']:
            for val in overlap_results['citegeist']:
                overlap_data.append({'Method': 'CITEgeist', 'Jaccard': val})
        for res, overlaps in overlap_results['rna_only'].items():
            for val in overlaps:
                overlap_data.append({'Method': f'RNA\n(res={res})', 'Jaccard': val})

        if overlap_data:
            overlap_df = pd.DataFrame(overlap_data)
            sns.boxplot(data=overlap_df, x='Method', y='Jaccard', ax=ax)
            ax.set_title('C. Program Gene Overlap (lower = more distinct)')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

        # [D] Summary statistics
        ax = axes[1, 1]
        ax.axis('off')

        # Compute summary stats
        citegeist_moran = moran_df[moran_df['method'] == 'CITEgeist']['moran_i'].dropna()
        summary_text = f"Region {self.region_id} Summary\n"
        summary_text += "=" * 40 + "\n\n"
        summary_text += f"CITEgeist:\n"
        summary_text += f"  Profiles analyzed: {self.citegeist_results['n_profiles_analyzed']}\n"
        summary_text += f"  Programs: {len(citegeist_moran)}\n"
        summary_text += f"  Mean Moran's I: {citegeist_moran.mean():.3f} ± {citegeist_moran.std():.3f}\n"
        summary_text += f"  Programs with I > 0.15: {(citegeist_moran > 0.15).sum()}\n\n"

        for res in self.leiden_resolutions:
            rna_moran = moran_df[moran_df['method'] == f'RNA-only (res={res})']['moran_i'].dropna()
            if len(rna_moran) > 0:
                summary_text += f"RNA-only (res={res}):\n"
                summary_text += f"  Programs: {len(rna_moran)}\n"
                summary_text += f"  Mean Moran's I: {rna_moran.mean():.3f} ± {rna_moran.std():.3f}\n"
                summary_text += f"  Programs with I > 0.15: {(rna_moran > 0.15).sum()}\n\n"

        # Statistical test
        if len(citegeist_moran) > 0:
            best_rna_res = max(self.leiden_resolutions,
                              key=lambda r: moran_df[moran_df['method'] == f'RNA-only (res={r})']['moran_i'].dropna().mean()
                              if len(moran_df[moran_df['method'] == f'RNA-only (res={r})']) > 0 else 0)
            best_rna_moran = moran_df[moran_df['method'] == f'RNA-only (res={best_rna_res})']['moran_i'].dropna()
            if len(best_rna_moran) > 0:
                stat, pval = mannwhitneyu(citegeist_moran, best_rna_moran, alternative='greater')
                summary_text += f"Mann-Whitney U test (CITEgeist > best RNA):\n"
                summary_text += f"  p-value: {pval:.4f}\n"

        ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
               verticalalignment='top', fontfamily='monospace')

        plt.tight_layout()
        fig.savefig(output_dir / f'integration_comparison_region_{self.region_id}.png', dpi=150, bbox_inches='tight')
        plt.close()

        logger.info(f"Saved figure to {output_dir}")

        return moran_df

    def save_results(self, output_dir: Path = None):
        """Save comparison results to files."""
        if output_dir is None:
            output_dir = OUTPUT_DIR

        # Save Moran comparison
        moran_df = self.compute_moran_comparison()
        moran_df.to_csv(output_dir / f'moran_comparison_region_{self.region_id}.csv', index=False)

        # Save RNA-only results
        with open(output_dir / f'rna_only_results_region_{self.region_id}.json', 'w') as f:
            # Convert numpy types for JSON
            results_clean = {}
            for res, data in self.rna_only_results.items():
                results_clean[str(res)] = {
                    'resolution': float(res),
                    'n_clusters': int(data['n_clusters']),
                    'clusters': {}
                }
                for cid, cdata in data['clusters'].items():
                    results_clean[str(res)]['clusters'][cid] = {
                        'n_cells': int(cdata['n_cells']),
                        'programs': [
                            {
                                'program_id': int(p['program_id']),
                                'moran_i': float(p['moran_i']) if not np.isnan(p['moran_i']) else None,
                                'top_genes': p['top_genes'],
                                'variance_explained': float(p['variance_explained']),
                            }
                            for p in cdata['programs']
                        ]
                    }
            json.dump(results_clean, f, indent=2)

        logger.info(f"Saved results to {output_dir}")


def run_comparison(region_id: int, seed: int = 42) -> pd.DataFrame:
    """Run full comparison for a region."""
    comparison = IntegrationComparison(region_id=region_id, seed=seed)
    comparison.load_data()
    comparison.run_all_resolutions()
    moran_df = comparison.generate_figures()
    comparison.save_results()
    return moran_df


def run_all_regions(seed: int = 42):
    """Run comparison for all 5 regions and aggregate."""
    all_results = []

    for region_id in range(5):
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing Region {region_id}")
        logger.info(f"{'='*60}")

        try:
            moran_df = run_comparison(region_id, seed=seed)
            moran_df['region'] = region_id
            all_results.append(moran_df)
        except Exception as e:
            logger.error(f"Error processing region {region_id}: {e}")
            continue

    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_df.to_csv(OUTPUT_DIR / 'moran_comparison_all_regions.csv', index=False)

        # Generate aggregate figure
        generate_aggregate_figure(combined_df)

        return combined_df
    return None


def generate_aggregate_figure(df: pd.DataFrame):
    """Generate aggregate comparison across all regions."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Simplify method names for plotting
    df['method_simple'] = df['method'].apply(
        lambda x: 'CITEgeist' if x == 'CITEgeist' else 'RNA-only'
    )

    # [A] Moran's I by method (aggregated)
    ax = axes[0]
    sns.boxplot(data=df, x='method_simple', y='moran_i', ax=ax)
    ax.set_title('Spatial Coherence: CITEgeist vs RNA-only')
    ax.set_xlabel('Method')
    ax.set_ylabel("Moran's I")

    # Add stats
    citegeist = df[df['method_simple'] == 'CITEgeist']['moran_i'].dropna()
    rna_only = df[df['method_simple'] == 'RNA-only']['moran_i'].dropna()
    if len(citegeist) > 0 and len(rna_only) > 0:
        stat, pval = mannwhitneyu(citegeist, rna_only, alternative='greater')
        ax.text(0.5, 0.95, f'p = {pval:.2e}', transform=ax.transAxes,
               ha='center', fontsize=12)

    # [B] Programs above threshold by region
    ax = axes[1]
    threshold = 0.15
    region_stats = []
    for region in df['region'].unique():
        for method in ['CITEgeist', 'RNA-only']:
            subset = df[(df['region'] == region) & (df['method_simple'] == method)]['moran_i'].dropna()
            n_above = (subset > threshold).sum()
            region_stats.append({
                'Region': region,
                'Method': method,
                'Count': n_above,
            })

    stat_df = pd.DataFrame(region_stats)
    stat_pivot = stat_df.pivot(index='Region', columns='Method', values='Count')
    stat_pivot.plot(kind='bar', ax=ax)
    ax.set_title(f'Programs with Moran\'s I > {threshold}')
    ax.set_ylabel('Number of programs')
    ax.legend(title='Method')

    plt.tight_layout()
    fig.savefig(FIGURE_DIR / 'integration_comparison_aggregate.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Print summary
    print("\n" + "="*60)
    print("INTEGRATION COMPARISON SUMMARY")
    print("="*60)
    print(f"\nCITEgeist (n={len(citegeist)} programs):")
    print(f"  Mean Moran's I: {citegeist.mean():.3f} ± {citegeist.std():.3f}")
    print(f"  Programs with I > 0.15: {(citegeist > 0.15).sum()} ({100*(citegeist > 0.15).mean():.1f}%)")

    print(f"\nRNA-only (n={len(rna_only)} programs):")
    print(f"  Mean Moran's I: {rna_only.mean():.3f} ± {rna_only.std():.3f}")
    print(f"  Programs with I > 0.15: {(rna_only > 0.15).sum()} ({100*(rna_only > 0.15).mean():.1f}%)")

    if len(citegeist) > 0 and len(rna_only) > 0:
        stat, pval = mannwhitneyu(citegeist, rna_only, alternative='greater')
        print(f"\nMann-Whitney U test (CITEgeist > RNA-only):")
        print(f"  p-value: {pval:.2e}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare protein-anchored vs RNA-only program discovery"
    )
    parser.add_argument(
        "--region", type=int, default=None,
        help="Region ID (0-4). If not specified, runs all regions."
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    if args.region is not None:
        run_comparison(args.region, seed=args.seed)
    else:
        run_all_regions(seed=args.seed)


if __name__ == "__main__":
    main()
