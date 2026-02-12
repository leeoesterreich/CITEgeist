#!/usr/bin/env python
"""
Analyze Treatment Response: Biopsy (S1) vs Surgical (S2) Comparison.

This script uses Modules 4, 4b, and 5 to discover:
1. Gene expression programs in each sample
2. Spatial relationships between programs
3. Which programs/relationships are conserved vs treatment-specific

Biological questions:
- What programs are UP/DOWN after treatment?
- Do cell type spatial relationships change?
- What tumor signatures persist in residual disease?
- Which responses are patient-specific vs shared?

Usage:
    python examples/analyze_treatment_response.py --mode full
    python examples/analyze_treatment_response.py --mode integration-only
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model import (
    # Core model
    CitegeistModel,
    # Module 4
    discover_programs_from_layers,
    AnchoredProgramDiscoveryResult,
    # Module 4b
    analyze_program_relationships,
    BivariateProgramResult,
    # Module 5
    integrate_samples,
    save_integration_results,
    IntegrationResult,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

DATA_DIR = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
OUTPUT_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output")

# Patient sample mapping
PATIENT_SAMPLES = {
    "P1": {"biopsy": "HCC22-088-P1-S1", "surgical": "HCC22-088-P1-S2"},
    "P2": {"biopsy": "HCC22-088-P2-S1", "surgical": "HCC22-088-P2-S2"},
    "P3": {"biopsy": "HCC22-088-P3-S1_A", "surgical": "HCC22-088-P3-S2"},
    "P4": {"biopsy": "HCC22-088-P4-S1", "surgical": "HCC22-088-P4-S2_1i_rep"},
    "P5": {"biopsy": "HCC22-088-P5-S1", "surgical": "HCC22-088-P5-S2_F_rep"},
    "P6": {"biopsy": "HCC22-088-P6-S1", "surgical": "HCC22-088-P6-S2_D"},
}


# =============================================================================
# Helper Functions
# =============================================================================


def load_sample_data(sample_name: str) -> Tuple[sc.AnnData, sc.AnnData]:
    """Load GEX and antibody data for a sample."""
    sample_dir = DATA_DIR / sample_name / "outs"

    if not sample_dir.exists():
        raise FileNotFoundError(f"Sample directory not found: {sample_dir}")

    logger.info(f"Loading {sample_name}...")

    adata = sq.read.visium(
        str(sample_dir),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False,
    )

    # Split into GEX and antibody
    feature_types = adata.var['feature_types']
    gex_mask = feature_types == 'Gene Expression'
    ab_mask = feature_types == 'Antibody Capture'

    adata_gex = adata[:, gex_mask].copy()
    adata_cite = adata[:, ab_mask].copy()

    adata_gex.var_names_make_unique()

    logger.info(f"  GEX: {adata_gex.shape[1]} genes, Antibody: {adata_cite.shape[1]} proteins")

    return adata_gex, adata_cite


def load_deconvolved_layers(adata_gex: sc.AnnData, sample_name: str) -> sc.AnnData:
    """Load Module 3 deconvolved layers into adata."""
    layers_dir = OUTPUT_DIR / f"{sample_name}_pass1" / "layers" / "pass1"

    if not layers_dir.exists():
        logger.warning(f"No deconvolved layers found for {sample_name}")
        return adata_gex

    layer_files = list(layers_dir.glob("*_layer_pass1.csv"))

    for layer_file in layer_files:
        cell_type = layer_file.stem.replace("_layer_pass1", "")
        df = pd.read_csv(layer_file, index_col=0)

        # Align to adata
        common_spots = list(adata_gex.obs_names.intersection(df.index))
        if len(common_spots) == 0:
            continue

        layer_matrix = np.zeros((adata_gex.shape[0], adata_gex.shape[1]))

        adata_gene_to_idx = {g: i for i, g in enumerate(adata_gex.var_names)}
        common_genes = [g for g in df.columns if g in adata_gene_to_idx]

        if len(common_genes) == 0:
            continue

        spot_idx = np.array([adata_gex.obs_names.get_loc(s) for s in common_spots])
        gene_idx = np.array([adata_gene_to_idx[g] for g in common_genes])

        layer_data = df.loc[common_spots, common_genes].values
        for i, s_idx in enumerate(spot_idx):
            layer_matrix[s_idx, gene_idx] = layer_data[i, :]

        layer_name = f"{cell_type}_genes_pass1"
        adata_gex.layers[layer_name] = layer_matrix

    n_layers = len([l for l in adata_gex.layers.keys() if '_genes_pass1' in l])
    logger.info(f"  Loaded {n_layers} deconvolved layers")

    return adata_gex


def run_module4_for_sample(
    sample_name: str,
    force_rerun: bool = False,
) -> AnchoredProgramDiscoveryResult:
    """Run Module 4 for a single sample."""
    output_file = OUTPUT_DIR / f"{sample_name}_module4_v3_programs.csv"

    if output_file.exists() and not force_rerun:
        logger.info(f"Module 4 output exists for {sample_name}, loading...")
        # Load from CSV (simplified - full W matrix not available)
        df = pd.read_csv(output_file)
        # For full analysis, would need to re-run or store full results
        return None

    # Load data
    adata_gex, adata_cite = load_sample_data(sample_name)
    adata_gex = load_deconvolved_layers(adata_gex, sample_name)

    # Check if layers exist
    layer_names = [l for l in adata_gex.layers.keys() if '_genes_pass1' in l]
    if len(layer_names) == 0:
        logger.warning(f"No deconvolved layers for {sample_name}, skipping Module 4")
        return None

    # Load cell type proportions
    prop_file = OUTPUT_DIR / f"{sample_name}_cell_prop_finetuned_results.csv"
    cell_type_proportions = None
    if prop_file.exists():
        props_raw = pd.read_csv(prop_file, index_col=0)
        cell_type_proportions = pd.DataFrame(
            0.0, index=adata_gex.obs_names, columns=props_raw.columns
        )
        common = adata_gex.obs_names.intersection(props_raw.index)
        cell_type_proportions.loc[common, :] = props_raw.loc[common, :]

    # Run Module 4
    logger.info(f"Running Module 4 for {sample_name}...")
    result = discover_programs_from_layers(
        adata=adata_gex,
        layer_pattern="_genes_pass1",
        cell_type_proportions=cell_type_proportions,
        protein_adata=adata_cite,
        K_programs=3,
        lambda_spatial=0.1,
        lambda_sparsity=0.01,
        validate_with_proteins=True,
        detect_subpopulations=True,
        n_subpopulations=3,
        random_state=42,
    )

    # Save results
    result.to_dataframe().to_csv(output_file)
    logger.info(f"Saved Module 4 results to {output_file}")

    return result


def run_module4b_for_sample(
    sample_name: str,
    module4_result: AnchoredProgramDiscoveryResult,
    force_rerun: bool = False,
) -> BivariateProgramResult:
    """Run Module 4b for a single sample."""
    output_file = OUTPUT_DIR / f"{sample_name}_module4b_relationships.csv"

    if output_file.exists() and not force_rerun:
        logger.info(f"Module 4b output exists for {sample_name}")
        return None

    if module4_result is None:
        logger.warning(f"No Module 4 result for {sample_name}, skipping Module 4b")
        return None

    # Load spatial coordinates
    adata_gex, _ = load_sample_data(sample_name)

    # Run Module 4b
    logger.info(f"Running Module 4b for {sample_name}...")
    result = analyze_program_relationships(
        module4_result=module4_result,
        spatial_coords=adata_gex.obsm['spatial'],
        n_neighbors=6,
        n_permutations=999,
        fdr_threshold=0.1,
    )

    # Save results
    result.to_dataframe().to_csv(output_file)
    logger.info(f"Saved Module 4b results to {output_file}")

    return result


# =============================================================================
# Module 5 Analysis
# =============================================================================


def run_cross_sample_integration(
    sample_groups: Dict[str, List[str]],
    analysis_name: str,
) -> IntegrationResult:
    """
    Run Module 5 integration on specified sample groups.

    Args:
        sample_groups: Dict mapping group name to list of sample names
        analysis_name: Name for this analysis (used in output files)
    """
    from CITEgeist.model.cross_sample_integration import (
        load_module4_from_csv,
        load_module4b_from_csv,
    )

    logger.info(f"\n{'='*70}")
    logger.info(f"Module 5 Integration: {analysis_name}")
    logger.info(f"{'='*70}")

    # Collect all samples
    all_samples = []
    sample_labels = {}
    for group, samples in sample_groups.items():
        for s in samples:
            all_samples.append(s)
            sample_labels[s] = group

    logger.info(f"Integrating {len(all_samples)} samples across {len(sample_groups)} groups")

    # Load Module 4 results
    module4_results = {}
    for sample in all_samples:
        m4_file = OUTPUT_DIR / f"{sample}_module4_v3_programs.csv"
        if not m4_file.exists():
            m4_file = OUTPUT_DIR / f"{sample}_module4_programs.csv"

        if m4_file.exists():
            results_by_anchor = load_module4_from_csv(m4_file, sample)
            from CITEgeist.model import AnchoredProgramDiscoveryResult
            module4_results[sample] = AnchoredProgramDiscoveryResult(
                results_by_anchor=results_by_anchor,
                n_anchors=len(results_by_anchor),
                total_programs=sum(len(r.programs) for r in results_by_anchor.values()),
                profile_discovery_result=None,
                parameters={"source_file": str(m4_file), "group": sample_labels[sample]},
            )
            logger.info(f"  Loaded {sample}: {len(results_by_anchor)} anchors")
        else:
            logger.warning(f"  Missing Module 4 for {sample}")

    if len(module4_results) < 2:
        logger.error("Need at least 2 samples for integration")
        return None

    # Load Module 4b results
    module4b_results = {}
    for sample in all_samples:
        m4b_file = OUTPUT_DIR / f"{sample}_module4b_relationships.csv"
        if m4b_file.exists():
            module4b_results[sample] = load_module4b_from_csv(m4b_file)
            logger.info(f"  Loaded {sample} relationships: {module4b_results[sample].n_pairs_tested} pairs")

    # Run integration
    result = integrate_samples(
        module4_results=module4_results,
        module4b_results=module4b_results if module4b_results else None,
        n_components=30,
        n_clusters=min(50, len(module4_results) * 10),
        theta=2.0,
        similarity_threshold=0.6,
        max_iter=30,
    )

    # Add group labels to metadata
    result.parameters["sample_groups"] = sample_groups
    result.parameters["sample_labels"] = sample_labels

    # Save results
    output_prefix = f"module5_{analysis_name}"
    saved = save_integration_results(result, OUTPUT_DIR, output_prefix)

    return result


def analyze_treatment_changes(
    biopsy_result: IntegrationResult,
    surgical_result: IntegrationResult,
) -> pd.DataFrame:
    """
    Compare biopsy vs surgical integration results to find treatment-specific changes.

    Returns DataFrame of programs and their conservation in each condition.
    """
    # Get program summaries
    biopsy_progs = {p.program_id: p for p in biopsy_result.aligned_programs}
    surgical_progs = {p.program_id: p for p in surgical_result.aligned_programs}

    # Compare top genes to find similar programs
    records = []

    for prog in biopsy_result.aligned_programs:
        top_genes_set = set(prog.top_genes[:10])

        # Find matching program in surgical
        best_match = None
        best_overlap = 0
        for s_prog in surgical_result.aligned_programs:
            overlap = len(top_genes_set & set(s_prog.top_genes[:10]))
            if overlap > best_overlap:
                best_overlap = overlap
                best_match = s_prog

        records.append({
            "biopsy_program": prog.program_id,
            "cell_type": prog.cell_type,
            "biopsy_conservation": prog.conservation_score,
            "surgical_match": best_match.program_id if best_match else None,
            "surgical_conservation": best_match.conservation_score if best_match else 0,
            "gene_overlap": best_overlap,
            "top_genes": ", ".join(prog.top_genes[:5]),
            "change": "lost" if best_overlap < 3 else "conserved" if best_match else "unknown",
        })

    # Find surgical-specific programs
    matched_surgical = {r["surgical_match"] for r in records if r["surgical_match"]}
    for prog in surgical_result.aligned_programs:
        if prog.program_id not in matched_surgical:
            records.append({
                "biopsy_program": None,
                "cell_type": prog.cell_type,
                "biopsy_conservation": 0,
                "surgical_match": prog.program_id,
                "surgical_conservation": prog.conservation_score,
                "gene_overlap": 0,
                "top_genes": ", ".join(prog.top_genes[:5]),
                "change": "treatment_induced",
            })

    return pd.DataFrame(records)


# =============================================================================
# Main Analysis
# =============================================================================


def main():
    parser = argparse.ArgumentParser(description="Analyze treatment response with Modules 4-5")
    parser.add_argument(
        "--mode",
        choices=["full", "integration-only", "summary"],
        default="summary",
        help="Analysis mode"
    )
    parser.add_argument(
        "--force-rerun",
        action="store_true",
        help="Force re-run of Module 4/4b even if outputs exist"
    )
    args = parser.parse_args()

    logger.info("=" * 70)
    logger.info("Treatment Response Analysis: Biopsy (S1) vs Surgical (S2)")
    logger.info("=" * 70)

    # Show available data
    print("\nPatient Samples:")
    print("-" * 50)
    for patient, samples in PATIENT_SAMPLES.items():
        biopsy_exists = (DATA_DIR / samples["biopsy"] / "outs").exists()
        surgical_exists = (DATA_DIR / samples["surgical"] / "outs").exists()
        print(f"  {patient}: Biopsy={samples['biopsy']} ({'✓' if biopsy_exists else '✗'}), "
              f"Surgical={samples['surgical']} ({'✓' if surgical_exists else '✗'})")

    if args.mode == "summary":
        # Just show what's available
        print("\nExisting Module 4 outputs:")
        for f in OUTPUT_DIR.glob("*module4*programs.csv"):
            print(f"  {f.name}")

        print("\nExisting Module 4b outputs:")
        for f in OUTPUT_DIR.glob("*module4b*.csv"):
            print(f"  {f.name}")

        print("\nTo run full analysis: python examples/analyze_treatment_response.py --mode full")
        print("To run integration only: python examples/analyze_treatment_response.py --mode integration-only")
        return

    if args.mode == "full":
        # Run Module 4 and 4b for all samples
        module4_results = {}
        module4b_results = {}

        for patient, samples in PATIENT_SAMPLES.items():
            for sample_type, sample_name in samples.items():
                logger.info(f"\n{'='*50}")
                logger.info(f"Processing {patient} {sample_type}: {sample_name}")
                logger.info(f"{'='*50}")

                try:
                    m4 = run_module4_for_sample(sample_name, args.force_rerun)
                    if m4:
                        module4_results[sample_name] = m4
                        m4b = run_module4b_for_sample(sample_name, m4, args.force_rerun)
                        if m4b:
                            module4b_results[sample_name] = m4b
                except Exception as e:
                    logger.error(f"Error processing {sample_name}: {e}")

    # Run Module 5 integrations
    logger.info("\n" + "=" * 70)
    logger.info("Running Module 5 Cross-Sample Integration")
    logger.info("=" * 70)

    # Get available samples
    biopsy_samples = [s["biopsy"] for s in PATIENT_SAMPLES.values()
                      if (OUTPUT_DIR / f"{s['biopsy']}_module4_programs.csv").exists() or
                         (OUTPUT_DIR / f"{s['biopsy']}_module4_v3_programs.csv").exists()]

    surgical_samples = [s["surgical"] for s in PATIENT_SAMPLES.values()
                        if (OUTPUT_DIR / f"{s['surgical']}_module4_programs.csv").exists() or
                           (OUTPUT_DIR / f"{s['surgical']}_module4_v3_programs.csv").exists()]

    logger.info(f"Biopsy samples with Module 4: {biopsy_samples}")
    logger.info(f"Surgical samples with Module 4: {surgical_samples}")

    # Integration 1: All biopsies together
    if len(biopsy_samples) >= 2:
        logger.info("\n--- Integration 1: All Biopsies ---")
        biopsy_result = run_cross_sample_integration(
            {"biopsy": biopsy_samples},
            "all_biopsies"
        )
    else:
        logger.warning("Not enough biopsy samples for integration")
        biopsy_result = None

    # Integration 2: All surgicals together
    if len(surgical_samples) >= 2:
        logger.info("\n--- Integration 2: All Surgical Specimens ---")
        surgical_result = run_cross_sample_integration(
            {"surgical": surgical_samples},
            "all_surgical"
        )
    else:
        logger.warning("Not enough surgical samples for integration")
        surgical_result = None

    # Integration 3: Combined with group labels
    all_samples = biopsy_samples + surgical_samples
    if len(all_samples) >= 3:
        logger.info("\n--- Integration 3: Biopsy vs Surgical Combined ---")
        combined_result = run_cross_sample_integration(
            {"biopsy": biopsy_samples, "surgical": surgical_samples},
            "biopsy_vs_surgical"
        )

    # Analyze treatment changes
    if biopsy_result and surgical_result:
        logger.info("\n" + "=" * 70)
        logger.info("Treatment Response Summary")
        logger.info("=" * 70)

        changes_df = analyze_treatment_changes(biopsy_result, surgical_result)
        changes_df.to_csv(OUTPUT_DIR / "treatment_response_summary.csv", index=False)

        print("\nProgram Changes with Treatment:")
        print(changes_df.to_string())

        # Summary statistics
        n_conserved = len(changes_df[changes_df["change"] == "conserved"])
        n_lost = len(changes_df[changes_df["change"] == "lost"])
        n_induced = len(changes_df[changes_df["change"] == "treatment_induced"])

        print(f"\nSummary:")
        print(f"  Conserved programs: {n_conserved}")
        print(f"  Lost after treatment: {n_lost}")
        print(f"  Treatment-induced: {n_induced}")

    logger.info("\nAnalysis complete! Check output directory for results.")


if __name__ == "__main__":
    main()
