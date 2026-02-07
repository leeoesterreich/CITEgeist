#!/usr/bin/env python
"""
Vignette 4: ESR1-D538G Mutation and Midkine Secretion Programs

This script runs the full MDK analysis pipeline:
1. Spatial program discovery in cancer cells
2. D538G region enrichment analysis
3. MDK contextual gene extraction
4. Bulk RNA-seq validation (GSE89888)
5. Integrated findings

Research Question: Why does MDK secretion increase in MCF7-D538G but not T47D-D538G?
"""

import os
import sys
import warnings
import logging
from pathlib import Path

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal

# Optional GEO download
try:
    import GEOparse
    HAS_GEOPARSE = True
except ImportError:
    HAS_GEOPARSE = False
    print("GEOparse not installed. GEO download will be skipped.")

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Add parent directory
sys.path.insert(0, os.path.abspath(".."))

from model import CitegeistModel
from model import (
    discover_programs_from_layers,
    analyze_program_regions,
    compare_programs_by_region,
    extract_program_context_genes,
)

# =============================================================================
# Configuration
# =============================================================================
# Gurobi license handled via module load gurobi/12.0.3
LICENSE_FILE = None  # Will use environment variable from module
DATA_FOLDER = "/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/"
SAMPLE_NAME = "HCC22-088-P4-S2_1i_rep"
SAMPLE_PATH = os.path.join(DATA_FOLDER, SAMPLE_NAME, "outs")
OUTPUT_DIR = Path("output_vignette4_mdk")
OUTPUT_DIR.mkdir(exist_ok=True)
FIGURES_DIR = OUTPUT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Cell type profiles
CELL_PROFILES = {
    "Cancer Cells": {"Major": ["EPCAM-1"], "Minor": ["SDC1-1", "KRT5-1"]},
    "Macrophages": {"Major": ["CD68-1"], "Minor": ["CD14-1"]},
    "CD4 T Cells": {"Major": ["CD3E-1", "CD4-1"]},
    "CD8 T Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "B Cells": {"Major": ["MS4A1-1", "CD19-1"]},
    "Endothelial Cells": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
}


def download_geo_rnaseq(output_path: Path):
    """Download and process GSE89888 RNA-seq data."""
    if not HAS_GEOPARSE:
        logger.warning("GEOparse not installed. Skipping GEO download.")
        return None, None

    logger.info("Downloading GSE89888 (bulk RNA-seq)...")

    geo_dir = output_path / "GSE89888"
    geo_dir.mkdir(exist_ok=True)

    # Download GEO series
    gse = GEOparse.get_GEO(geo="GSE89888", destdir=str(geo_dir))

    # Extract sample metadata
    samples = []
    for gsm_name, gsm in gse.gsms.items():
        meta = gsm.metadata

        # Parse characteristics
        chars = {}
        for char in meta.get('characteristics_ch1', []):
            if ':' in char:
                key, val = char.split(':', 1)
                chars[key.strip()] = val.strip()

        samples.append({
            'sample_id': gsm_name,
            'title': meta.get('title', [''])[0],
            'cell_line': chars.get('cell line', ''),
            'esr1_status': chars.get('esr1 status', chars.get('genotype', '')),
            'treatment': chars.get('treatment', chars.get('agent', '')),
        })

    meta_df = pd.DataFrame(samples)
    meta_df.to_csv(geo_dir / "sample_metadata.csv", index=False)
    logger.info(f"Saved metadata for {len(meta_df)} samples")

    # Try to get expression data from supplementary files
    # GSE89888 should have processed counts
    supp_files = gse.metadata.get('supplementary_file', [])
    logger.info(f"Supplementary files: {supp_files}")

    return geo_dir, meta_df


def load_rnaseq_from_geo(geo_dir: Path):
    """Load RNA-seq counts from GEO download."""

    # Look for counts file
    counts_files = list(geo_dir.glob("*counts*.txt")) + list(geo_dir.glob("*counts*.csv"))
    counts_files += list(geo_dir.glob("*matrix*.txt")) + list(geo_dir.glob("*FPKM*.txt"))

    if counts_files:
        logger.info(f"Found counts file: {counts_files[0]}")
        counts = pd.read_csv(counts_files[0], sep='\t', index_col=0)
        return counts

    # If no pre-processed counts, try to extract from GSM
    logger.warning("No counts file found. Checking for processed data in GSMs...")

    # For GSE89888, the data might be in a different format
    # Let's check what's available
    for f in geo_dir.iterdir():
        logger.info(f"  File: {f.name}")

    return None


def run_spatial_analysis():
    """Run Parts 1-4: Spatial program discovery and region analysis."""
    logger.info("="*60)
    logger.info("PART 1-2: LOADING DATA AND RUNNING CITEGEIST")
    logger.info("="*60)

    # Load spatial data
    adata = sq.read.visium(
        SAMPLE_PATH,
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False
    )
    logger.info(f"Loaded: {adata.shape[0]} spots x {adata.shape[1]} features")

    # Initialize and run CITEgeist
    model = CitegeistModel(
        sample_name=SAMPLE_NAME,
        adata=adata,
        output_folder=str(OUTPUT_DIR / "citegeist_output")
    )

    model.load_cell_profile_dict(CELL_PROFILES)
    model.split_adata()
    model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=25)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()
    model.preprocess_antibody()

    # Register Gurobi - use env variable if LICENSE_FILE is None
    if LICENSE_FILE:
        model.register_gurobi(LICENSE_FILE)
    else:
        # Find license from environment (set by module load)
        import os
        grb_license = os.environ.get("GRB_LICENSE_FILE")
        if grb_license:
            model.register_gurobi(grb_license)
        else:
            # Try default locations
            for path in ["/ihome/crc/install/gurobi/gurobi1203/gurobi.lic",
                        os.path.expanduser("~/gurobi.lic")]:
                if os.path.exists(path):
                    model.register_gurobi(path)
                    break
            else:
                logger.warning("No Gurobi license found, trying without explicit registration")

    # Run deconvolution (radius auto-detected from spatial coordinates)
    model.run_cell_proportion_model()
    model.run_cell_expression_pass1()

    # Get results
    model.append_proportions_to_adata(key="finetuned")
    model.append_gex_to_adata(pass_number=1)
    prop_gex_adata = model.get_adata()

    # Define D538G regions using basal cytokeratins
    logger.info("="*60)
    logger.info("DEFINING D538G MUTATION REGIONS")
    logger.info("="*60)

    paper_keratins = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT16", "KRT17"]
    available_keratins = [g for g in paper_keratins if g in prop_gex_adata.var_names]

    if available_keratins and "Cancer_Cells_genes_pass1" in prop_gex_adata.layers:
        keratin_idx = [prop_gex_adata.var_names.get_loc(g) for g in available_keratins]
        prop_gex_adata.obs["Basal_Cytokeratin"] = prop_gex_adata.layers["Cancer_Cells_genes_pass1"][:, keratin_idx].sum(axis=1)

        # Threshold at median of non-zero values
        nonzero = prop_gex_adata.obs["Basal_Cytokeratin"][prop_gex_adata.obs["Basal_Cytokeratin"] > 0]
        threshold = nonzero.median() if len(nonzero) > 0 else 0
        prop_gex_adata.obs["D538G_Mutation"] = prop_gex_adata.obs["Basal_Cytokeratin"] > threshold
    else:
        logger.warning("Keratin genes not found. Using spatial clustering for D538G annotation.")
        # Fallback: use spatial clustering on cancer cell proportions
        cancer_prop = prop_gex_adata.obs.get("Cancer Cells", np.zeros(prop_gex_adata.n_obs))
        prop_gex_adata.obs["D538G_Mutation"] = cancer_prop > cancer_prop.median()

    n_d538g_pos = prop_gex_adata.obs["D538G_Mutation"].sum()
    n_d538g_neg = (~prop_gex_adata.obs["D538G_Mutation"]).sum()
    logger.info(f"D538G+ spots: {n_d538g_pos}, D538G- spots: {n_d538g_neg}")

    # Discover programs
    logger.info("="*60)
    logger.info("PART 2: MODULE 4 PROGRAM DISCOVERY")
    logger.info("="*60)

    result = discover_programs_from_layers(
        adata=prop_gex_adata,
        layer_pattern="_genes_pass1",
        K_programs=5,
        top_n_genes=50,
    )

    logger.info(f"Discovered programs for {len(result.results_by_anchor)} cell types")

    # Get Cancer Cells result specifically
    cancer_result = result.results_by_anchor.get("Cancer Cells")
    if not cancer_result:
        # Try alternative names
        for key in result.results_by_anchor:
            if "Cancer" in key or "EPCAM" in key:
                cancer_result = result.results_by_anchor[key]
                break

    if cancer_result:
        for k, prog in enumerate(cancer_result.programs):
            logger.info(f"  Program {k}: Moran's I = {prog.spatial_moran_i:.3f}, "
                       f"top genes: {', '.join(prog.top_genes[:5])}")
    else:
        logger.warning("Cancer Cells not found in results")

    # Region enrichment analysis
    logger.info("="*60)
    logger.info("PART 3: REGION ENRICHMENT ANALYSIS")
    logger.info("="*60)

    if not cancer_result:
        logger.error("No Cancer Cells result to analyze")
        return None

    # Get the spot indices used in program discovery
    # The H matrix corresponds to spots_used in the cancer_result
    n_spots_in_H = cancer_result.H.shape[1]
    logger.info(f"H matrix has {n_spots_in_H} spots, adata has {prop_gex_adata.n_obs} spots")

    # Use the full adata for region analysis - H matrix should match
    # Note: discover_programs_from_layers uses all spots in the layer
    cancer_result = analyze_program_regions(
        result=cancer_result,
        adata=prop_gex_adata,
        region_column="D538G_Mutation",
        min_spots_per_region=20
    )

    # For downstream analysis, still filter to cancer cell spots
    cancer_mask = prop_gex_adata.obs["Cancer Cells"] > 0.1
    adata_cancer = prop_gex_adata[cancer_mask].copy()

    # Summary table
    region_summary = []
    for prog in cancer_result.programs:
        region_summary.append({
            "Program": prog.program_id,
            "D538G+ Activity": prog.region_enrichment.get("True", 0) if prog.region_enrichment else 0,
            "D538G- Activity": prog.region_enrichment.get("False", 0) if prog.region_enrichment else 0,
            "Specificity": prog.region_specificity,
            "P-value": prog.region_pvalue,
            "Enriched In": prog.enriched_region,
        })

    region_df = pd.DataFrame(region_summary)
    logger.info("\nRegion Enrichment Summary:")
    print(region_df.to_string())
    region_df.to_csv(OUTPUT_DIR / "region_enrichment_summary.csv", index=False)

    # Compare regions (use full adata to match H matrix)
    comparison_df = compare_programs_by_region(
        result=cancer_result,
        adata=prop_gex_adata,
        region_column="D538G_Mutation",
        region_a=True,
        region_b=False,
        top_n_genes=50
    )
    comparison_df.to_csv(OUTPUT_DIR / "program_region_comparison.csv", index=False)

    # Identify D538G-enriched programs
    d538g_enriched = comparison_df[(comparison_df["fold_change"] > 1.5) & (comparison_df["pvalue"] < 0.05)]
    logger.info(f"\nD538G-enriched programs: {len(d538g_enriched)}")

    # MDK context discovery
    logger.info("="*60)
    logger.info("PART 4: MDK SECRETION CONTEXT DISCOVERY")
    logger.info("="*60)

    # Check MDK loading in each program
    W = cancer_result.W
    mdk_loadings = []

    if "MDK" in cancer_result.gene_names:
        mdk_idx = cancer_result.gene_names.index("MDK")
        for k in range(len(cancer_result.programs)):
            mdk_loadings.append({
                "Program": k,
                "MDK_Loading": W[mdk_idx, k],
                "D538G_Enriched": k in d538g_enriched["program_id"].values if len(d538g_enriched) > 0 else False,
                "Top_Genes": ", ".join(cancer_result.programs[k].top_genes[:5])
            })

        mdk_df = pd.DataFrame(mdk_loadings).sort_values("MDK_Loading", ascending=False)
        logger.info("\nMDK Loading by Program:")
        print(mdk_df.to_string())
        mdk_df.to_csv(OUTPUT_DIR / "mdk_program_loadings.csv", index=False)

        # Identify secretion program
        if len(d538g_enriched) > 0:
            d538g_mdk = mdk_df[mdk_df["D538G_Enriched"]]
            if len(d538g_mdk) > 0:
                secretion_program_id = int(d538g_mdk.iloc[0]["Program"])
            else:
                secretion_program_id = int(mdk_df.iloc[0]["Program"])
        else:
            secretion_program_id = int(mdk_df.iloc[0]["Program"])

        logger.info(f"\nIdentified Secretion Program: {secretion_program_id}")
    else:
        logger.warning("MDK not found in gene list")
        secretion_program_id = 0
        mdk_df = pd.DataFrame()

    # Extract contextual genes
    context_genes = extract_program_context_genes(
        result=cancer_result,
        program_id=secretion_program_id,
        target_gene="MDK",
        top_n=100,
        exclude_target=True
    )

    context_df = pd.DataFrame(context_genes, columns=["Gene", "Loading"])
    logger.info(f"\nTop 20 Contextual Factors (co-loaded with MDK):")
    print(context_df.head(20).to_string())
    context_df.to_csv(OUTPUT_DIR / "mdk_contextual_genes.csv", index=False)

    # Check MDK receptors
    receptor_genes = ["SDC4", "NCL", "PTPRZ1", "LRP1", "GPC1"]
    receptor_in_context = [g for g, _ in context_genes if g in receptor_genes]
    logger.info(f"\nMDK receptors in context: {receptor_in_context}")

    return {
        'prop_gex_adata': prop_gex_adata,
        'result': cancer_result,
        'comparison_df': comparison_df,
        'd538g_enriched': d538g_enriched,
        'context_genes': context_genes,
        'context_df': context_df,
        'secretion_program_id': secretion_program_id,
        'mdk_df': mdk_df,
    }


def run_rnaseq_validation(context_genes, geo_dir: Path):
    """Run Parts 5-6: Bulk RNA-seq validation."""
    logger.info("="*60)
    logger.info("PART 5-6: BULK RNA-SEQ VALIDATION (GSE89888)")
    logger.info("="*60)

    # Load metadata
    meta_path = geo_dir / "sample_metadata.csv"
    if not meta_path.exists():
        logger.error("Metadata not found. Run download first.")
        return None

    meta_df = pd.read_csv(meta_path)
    logger.info(f"Loaded metadata for {len(meta_df)} samples")
    print(meta_df.groupby(['cell_line', 'esr1_status']).size())

    # Try to load counts
    counts = load_rnaseq_from_geo(geo_dir)

    if counts is None:
        # Try alternative: download from GEO supplementary
        logger.info("Attempting to download expression data from GEO...")

        try:
            import urllib.request
            # GSE89888 has FPKM data
            url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE89888&format=file&file=GSE89888%5Fgene%5Fexpression%5FFPKM%2Etxt%2Egz"
            dest = geo_dir / "GSE89888_gene_expression_FPKM.txt.gz"

            if not dest.exists():
                logger.info(f"Downloading to {dest}...")
                urllib.request.urlretrieve(url, dest)

            # Load
            import gzip
            with gzip.open(dest, 'rt') as f:
                counts = pd.read_csv(f, sep='\t', index_col=0)
            logger.info(f"Loaded expression data: {counts.shape}")

        except Exception as e:
            logger.error(f"Failed to download expression data: {e}")
            logger.info("Please manually download GSE89888 expression data")
            return None

    if counts is not None:
        # Standardize column names
        logger.info(f"Expression columns: {counts.columns.tolist()[:10]}...")

        # Parse sample info from column names
        # Expected format varies - need to inspect
        logger.info("Parsing sample metadata from expression data...")

        # Match columns to metadata
        # This may require manual adjustment based on actual column format

        return counts

    return None


def compute_interaction_effect(gene, counts_df, meta_df):
    """Compute interaction effect for a gene.

    Interaction = (D538G - WT in MCF7) - (D538G - WT in T47D)
    Positive = gene responds more to D538G in MCF7 than T47D
    """
    if gene not in counts_df.index:
        return None

    expr = counts_df.loc[gene]

    results = {}
    for cell_line in ["MCF7", "T47D"]:
        mask_wt = (meta_df["cell_line"] == cell_line) & (meta_df["esr1_status"].str.contains("WT", case=False, na=False))
        mask_mut = (meta_df["cell_line"] == cell_line) & (meta_df["esr1_status"].str.contains("D538G", case=False, na=False))

        wt_samples = meta_df[mask_wt]["sample_id"].tolist()
        mut_samples = meta_df[mask_mut]["sample_id"].tolist()

        wt_vals = expr[[s for s in wt_samples if s in expr.index]]
        mut_vals = expr[[s for s in mut_samples if s in expr.index]]

        if len(wt_vals) == 0 or len(mut_vals) == 0:
            return None

        results[f"{cell_line}_wt_mean"] = wt_vals.mean()
        results[f"{cell_line}_mut_mean"] = mut_vals.mean()
        results[f"{cell_line}_fc"] = (mut_vals.mean() + 1) / (wt_vals.mean() + 1)
        results[f"{cell_line}_diff"] = mut_vals.mean() - wt_vals.mean()

    # Interaction score
    if "MCF7_diff" in results and "T47D_diff" in results:
        results["interaction"] = results["MCF7_diff"] - results["T47D_diff"]

    return results


def generate_summary(spatial_results, rnaseq_results=None):
    """Generate integrated findings and summary."""
    logger.info("="*60)
    logger.info("PART 9: INTEGRATED FINDINGS")
    logger.info("="*60)

    context_df = spatial_results['context_df']
    d538g_enriched = spatial_results['d538g_enriched']

    print("\n" + "="*60)
    print("HIGH-CONFIDENCE PERMISSIVE FACTORS FOR MDK SECRETION")
    print("="*60)

    print(f"\n1. SPATIAL DISCOVERY:")
    print(f"   - Analyzed {len(spatial_results['result'].programs)} programs in cancer cells")
    print(f"   - Identified {len(d538g_enriched)} D538G-enriched programs")
    print(f"   - Extracted {len(spatial_results['context_genes'])} genes co-loaded with MDK")

    print(f"\n2. TOP MDK CONTEXTUAL GENES:")
    for i, (gene, loading) in enumerate(spatial_results['context_genes'][:15]):
        print(f"   {i+1:2d}. {gene:<12s} (loading: {loading:.4f})")

    # Biological interpretation
    print("\n" + "="*60)
    print("BIOLOGICAL INTERPRETATION")
    print("="*60)

    # Check for known pathways in context genes
    context_gene_list = [g for g, _ in spatial_results['context_genes'][:50]]

    # Known ESR1/ER signaling genes
    er_related = ["ESR1", "ESR2", "PGR", "GREB1", "TFF1", "CCND1", "MYC", "FOXA1", "XBP1"]
    er_in_context = [g for g in er_related if g in context_gene_list]

    # Known secretion/autocrine genes
    secretion_related = ["SDC4", "NCL", "PTPRZ1", "GPC1", "LRP1", "EGF", "EGFR", "ERBB2", "ERBB3"]
    secretion_in_context = [g for g in secretion_related if g in context_gene_list]

    # Transcription factors
    tfs = ["FOXA1", "GATA3", "RUNX1", "STAT1", "STAT3", "NFE2L2", "JUN", "FOS", "MYC"]
    tfs_in_context = [g for g in tfs if g in context_gene_list]

    print(f"\n   ER-related genes in MDK context: {er_in_context}")
    print(f"   Secretion pathway genes: {secretion_in_context}")
    print(f"   Transcription factors: {tfs_in_context}")

    print("\n" + "="*60)
    print("HYPOTHESIS: WHY MCF7 BUT NOT T47D?")
    print("="*60)

    print("""
   The contextual genes co-expressed with MDK in D538G+ regions may represent
   the "permissive context" that enables D538G-driven MDK secretion.

   POSSIBLE EXPLANATIONS:

   1. DIFFERENTIAL TF EXPRESSION
      - MCF7 and T47D have different baseline TF profiles
      - If contextual TFs are higher in MCF7, they may cooperate with D538G-ER
      - Check: FOXA1, GATA3 baseline expression differences

   2. EPIGENETIC DIFFERENCES
      - D538G-ER may only access MDK enhancers in permissive chromatin
      - MCF7 may have more accessible chromatin at MDK locus
      - Check: ChIP-seq for H3K27ac, ATAC-seq differences

   3. RECEPTOR AUTOCRINE LOOP
      - If MDK receptors (SDC4, NCL) are in context, autocrine signaling may
        amplify the initial D538G effect - but only if receptors are expressed
      - MCF7 may have higher baseline MDK receptor expression

   4. COFACTOR AVAILABILITY
      - D538G-ER recruits different coactivators than WT-ER
      - These cofactors may be differentially available in MCF7 vs T47D
""")

    # Save summary
    with open(OUTPUT_DIR / "biological_summary.txt", "w") as f:
        f.write("MDK Contextual Genes Analysis Summary\n")
        f.write("="*50 + "\n\n")
        f.write(f"Top contextual genes co-loaded with MDK:\n")
        for gene, loading in spatial_results['context_genes'][:20]:
            f.write(f"  {gene}: {loading:.4f}\n")
        f.write(f"\nER-related: {er_in_context}\n")
        f.write(f"Secretion pathway: {secretion_in_context}\n")
        f.write(f"Transcription factors: {tfs_in_context}\n")

    logger.info(f"\nResults saved to {OUTPUT_DIR}/")


def main():
    """Run full vignette analysis."""
    logger.info("Starting Vignette 4: ESR1-D538G MDK Analysis")

    # Part 1-4: Spatial analysis
    spatial_results = run_spatial_analysis()

    # Part 5-6: RNA-seq validation (download if needed)
    geo_dir = OUTPUT_DIR / "GSE89888"
    if not geo_dir.exists():
        try:
            geo_dir, meta_df = download_geo_rnaseq(OUTPUT_DIR)
        except Exception as e:
            logger.warning(f"Could not download GEO data: {e}")
            geo_dir = None

    rnaseq_results = None
    if geo_dir and geo_dir.exists():
        rnaseq_results = run_rnaseq_validation(
            spatial_results['context_genes'],
            geo_dir
        )

    # Part 9: Generate summary
    generate_summary(spatial_results, rnaseq_results)

    logger.info("Analysis complete!")
    return spatial_results


if __name__ == "__main__":
    main()
