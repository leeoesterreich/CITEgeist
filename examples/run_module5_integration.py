#!/usr/bin/env python
"""
Module 5 Cross-Sample Integration Runner.

Integrates gene expression programs across all 12 patient samples using
Harmony-style batch correction, producing conserved programs and relationships.

Usage:
    python run_module5_integration.py --output-dir output/module5_integration

    # With custom parameters:
    python run_module5_integration.py --n-components 30 --theta 2.0
"""
import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import networkx as nx
import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.anchored_program_discovery import (
    AnchoredProgramDiscoveryResult,
    AnchoredProgramResult,
    BivariateProgramResult,
    ProgramPairRelationship,
    SpatialProgram,
)
from CITEgeist.model.cross_sample_integration import (
    IntegrationResult,
    integrate_samples,
    save_integration_results,
)

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# Data locations
MODULE4_OUTPUT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/module4_unified")
MODULE4B_OUTPUT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/module4b_unified")

# All 12 patient samples (deduplicated: P4-S2_1i_rep replaces P4-S2, P5-S2_F_rep replaces P5-S2)
ALL_SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

# Sample metadata for response analysis
SAMPLE_METADATA = {
    "HCC22-088-P1-S1": {"patient": "P1", "timepoint": "S1", "response": "Progressor"},
    "HCC22-088-P1-S2": {"patient": "P1", "timepoint": "S2", "response": "Progressor"},
    "HCC22-088-P2-S1": {"patient": "P2", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P2-S2": {"patient": "P2", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P3-S1_A": {"patient": "P3", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P3-S2": {"patient": "P3", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P4-S1": {"patient": "P4", "timepoint": "S1", "response": "Progressor"},
    "HCC22-088-P4-S2_1i_rep": {"patient": "P4", "timepoint": "S2", "response": "Progressor"},
    "HCC22-088-P5-S1": {"patient": "P5", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P5-S2_F_rep": {"patient": "P5", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P6-S1": {"patient": "P6", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P6-S2_D": {"patient": "P6", "timepoint": "S2", "response": "Responder"},
}


def load_module4_json(sample_name: str) -> Optional[AnchoredProgramDiscoveryResult]:
    """
    Load Module 4 results from JSON file.

    Args:
        sample_name: Sample identifier

    Returns:
        AnchoredProgramDiscoveryResult or None if not found
    """
    json_file = MODULE4_OUTPUT / f"{sample_name}_module4_discovery.json"
    h_file = MODULE4_OUTPUT / f"{sample_name}_anchored_H.npy"

    if not json_file.exists():
        logger.warning(f"Module 4 JSON not found: {json_file}")
        return None

    with open(json_file) as f:
        data = json.load(f)

    if "anchored" not in data:
        logger.warning(f"No anchored results in {json_file}")
        return None

    # Load H matrices if available
    H_matrices = {}
    if h_file.exists():
        H_matrices = np.load(h_file, allow_pickle=True).item()

    # Convert JSON to AnchoredProgramResult objects
    results_by_anchor = {}

    for anchor_name, anchor_data in data["anchored"]["anchors"].items():
        programs = []

        for prog_data in anchor_data["programs"]:
            top_genes = prog_data["top_genes"]
            # Create approximate gene loadings from rank
            gene_loadings = {g: 1.0 / (i + 1) for i, g in enumerate(top_genes)}

            programs.append(SpatialProgram(
                program_id=prog_data["program_id"],
                top_genes=top_genes,
                gene_loadings=gene_loadings,
                variance_explained=prog_data["variance_explained"],
                spatial_moran_i=prog_data["spatial_moran_i"],
                spatial_moran_pvalue=prog_data["spatial_moran_pvalue"],
                mean_activity=prog_data["mean_activity"],
                active_spots_fraction=0.0,  # Not stored in JSON
            ))

        # Create W matrix from gene loadings
        all_genes = list(set(g for p in programs for g in p.top_genes))
        n_progs = len(programs)
        W = np.zeros((len(all_genes), n_progs))

        for k, prog in enumerate(programs):
            for g, loading in prog.gene_loadings.items():
                if g in all_genes:
                    W[all_genes.index(g), k] = loading

        # Get H matrix from saved file or use placeholder
        H = H_matrices.get(anchor_name, np.zeros((n_progs, 1)))

        results_by_anchor[anchor_name] = AnchoredProgramResult(
            anchor_name=anchor_name,
            anchor_proteins=[],
            programs=programs,
            W=W,
            H=H,
            gene_names=all_genes,
            protein_correlations=pd.DataFrame(),
            reconstruction_error=0.0,
            n_spots_used=anchor_data["n_spots_used"],
            parameters={"sample_name": sample_name},
        )

    return AnchoredProgramDiscoveryResult(
        results_by_anchor=results_by_anchor,
        n_anchors=len(results_by_anchor),
        total_programs=sum(len(r.programs) for r in results_by_anchor.values()),
        profile_discovery_result=None,
        parameters={"source_file": str(json_file)},
    )


def load_all_module4_results() -> Dict[str, AnchoredProgramDiscoveryResult]:
    """Load Module 4 results for all samples."""
    results = {}

    for sample_name in ALL_SAMPLES:
        result = load_module4_json(sample_name)
        if result is not None:
            results[sample_name] = result
            logger.info(f"Loaded {sample_name}: {result.n_anchors} anchors, {result.total_programs} programs")

    logger.info(f"Loaded Module 4 results for {len(results)}/{len(ALL_SAMPLES)} samples")
    return results


def load_module4b_csv(sample_name: str) -> Optional[BivariateProgramResult]:
    """
    Load Module 4b results from CSV file.

    Args:
        sample_name: Sample identifier

    Returns:
        BivariateProgramResult or None if not found
    """
    csv_file = MODULE4B_OUTPUT / f"{sample_name}_module4b_relationships.csv"

    if not csv_file.exists():
        logger.warning(f"Module 4b CSV not found: {csv_file}")
        return None

    df = pd.read_csv(csv_file)

    if len(df) == 0:
        logger.info(f"No relationships in {csv_file}")
        return BivariateProgramResult(
            significant_pairs=[],
            all_pairs=[],
            n_programs_total=0,
            n_pairs_tested=0,
            n_significant=0,
            fdr_threshold=0.05,
        )

    all_pairs = []
    for _, row in df.iterrows():
        pair = ProgramPairRelationship(
            anchor1=row["anchor1"],
            program1_id=int(row["program1_id"]),
            anchor2=row["anchor2"],
            program2_id=int(row["program2_id"]),
            bivariate_morans_i=row["bivariate_morans_i"],
            bivariate_pvalue=row.get("bivariate_pvalue", 0.05),
            pearson_correlation=row.get("pearson_correlation", 0.0),
            pearson_pvalue=row.get("pearson_pvalue", 1.0),
            n_spots_used=int(row.get("n_spots_used", 0)),
            top_genes_overlap=row.get("top_genes_overlap", "").split(",") if pd.notna(row.get("top_genes_overlap", "")) else [],
            relationship_type=row.get("relationship_type", "unknown"),
        )
        all_pairs.append(pair)

    significant = [p for p in all_pairs if p.bivariate_pvalue < 0.05]

    return BivariateProgramResult(
        significant_pairs=significant,
        all_pairs=all_pairs,
        n_programs_total=len(set(
            [(p.anchor1, p.program1_id) for p in all_pairs] +
            [(p.anchor2, p.program2_id) for p in all_pairs]
        )),
        n_pairs_tested=len(all_pairs),
        n_significant=len(significant),
        fdr_threshold=0.05,
    )


def load_all_module4b_results() -> Dict[str, BivariateProgramResult]:
    """Load Module 4b results for all samples."""
    results = {}

    for sample_name in ALL_SAMPLES:
        result = load_module4b_csv(sample_name)
        if result is not None:
            results[sample_name] = result
            logger.info(f"Loaded {sample_name}: {result.n_significant} significant pairs")

    logger.info(f"Loaded Module 4b results for {len(results)}/{len(ALL_SAMPLES)} samples")
    return results


def analyze_response_patterns(
    result: IntegrationResult,
    metadata: Dict,
) -> Dict:
    """
    Analyze programs by treatment response (Responder vs Progressor).

    Returns summary of which programs are enriched in each group.
    """
    analysis = {
        "responder_enriched": [],
        "progressor_enriched": [],
        "timepoint_differential": [],
    }

    responder_samples = [s for s, m in metadata.items() if m["response"] == "Responder"]
    progressor_samples = [s for s, m in metadata.items() if m["response"] == "Progressor"]

    for prog in result.aligned_programs:
        samples_present = set(prog.samples_present)

        resp_count = len(samples_present & set(responder_samples))
        prog_count = len(samples_present & set(progressor_samples))

        resp_frac = resp_count / len(responder_samples) if responder_samples else 0
        prog_frac = prog_count / len(progressor_samples) if progressor_samples else 0

        # Enriched if >2x more frequent in one group
        if resp_frac > 0 and prog_frac > 0:
            if resp_frac / prog_frac > 2:
                analysis["responder_enriched"].append({
                    "program_id": prog.program_id,
                    "cell_type": prog.cell_type,
                    "responder_frac": resp_frac,
                    "progressor_frac": prog_frac,
                    "top_genes": prog.top_genes[:5],
                })
            elif prog_frac / resp_frac > 2:
                analysis["progressor_enriched"].append({
                    "program_id": prog.program_id,
                    "cell_type": prog.cell_type,
                    "responder_frac": resp_frac,
                    "progressor_frac": prog_frac,
                    "top_genes": prog.top_genes[:5],
                })

    return analysis


def run_module5(
    output_dir: Path,
    n_components: int = 30,
    n_clusters: int = 50,
    theta: float = 2.0,
    similarity_threshold: float = 0.7,
    seed: int = 42,
) -> IntegrationResult:
    """
    Run Module 5 cross-sample integration.

    Args:
        output_dir: Output directory
        n_components: PCA components for Harmony
        n_clusters: Soft clusters for Harmony
        theta: Diversity penalty
        similarity_threshold: Program matching threshold
        seed: Random seed

    Returns:
        IntegrationResult
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all Module 4 results
    logger.info("Loading Module 4 results...")
    module4_results = load_all_module4_results()

    if len(module4_results) < 2:
        raise ValueError(f"Need at least 2 samples, found {len(module4_results)}")

    # Load Module 4b results (bivariate relationships)
    logger.info("Loading Module 4b results...")
    module4b_results = load_all_module4b_results()
    logger.info(f"Loaded {len(module4b_results)} samples with bivariate data")

    # Run integration
    logger.info("Running cross-sample integration...")
    result = integrate_samples(
        module4_results=module4_results,
        module4b_results=module4b_results if len(module4b_results) > 0 else None,
        n_components=n_components,
        n_clusters=n_clusters,
        theta=theta,
        similarity_threshold=similarity_threshold,
        random_state=seed,
        build_network=True,
    )

    # Save results
    logger.info("Saving results...")
    saved_files = save_integration_results(
        result=result,
        output_dir=output_dir,
        prefix="module5_unified",
    )

    # Response analysis
    logger.info("Analyzing response patterns...")
    response_analysis = analyze_response_patterns(result, SAMPLE_METADATA)

    analysis_file = output_dir / "module5_response_analysis.json"
    with open(analysis_file, "w") as f:
        json.dump(response_analysis, f, indent=2)

    # Save summary
    summary = {
        "timestamp": datetime.now().isoformat(),
        "n_samples": len(module4_results),
        "samples": list(module4_results.keys()),
        "n_aligned_programs": len(result.aligned_programs),
        "n_conserved_relationships": len(result.conserved_relationships),
        "n_highly_conserved_programs": len([p for p in result.aligned_programs if p.conservation_score > 0.5]),
        "harmony_converged": result.harmony_converged,
        "n_iterations": result.n_iterations,
        "parameters": {
            "n_components": n_components,
            "n_clusters": n_clusters,
            "theta": theta,
            "similarity_threshold": similarity_threshold,
        },
        "response_analysis": {
            "n_responder_enriched": len(response_analysis["responder_enriched"]),
            "n_progressor_enriched": len(response_analysis["progressor_enriched"]),
        },
        "saved_files": {k: str(v) for k, v in saved_files.items()},
    }

    summary_file = output_dir / "module5_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Saved summary to {summary_file}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run Module 5 cross-sample integration"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/module5_integration",
        help="Output directory for results",
    )
    parser.add_argument("--n-components", type=int, default=30, help="PCA components")
    parser.add_argument("--n-clusters", type=int, default=50, help="Harmony clusters")
    parser.add_argument("--theta", type=float, default=2.0, help="Diversity penalty")
    parser.add_argument("--similarity-threshold", type=float, default=0.7, help="Program matching threshold")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    # Run Module 5
    result = run_module5(
        output_dir=Path(args.output_dir),
        n_components=args.n_components,
        n_clusters=args.n_clusters,
        theta=args.theta,
        similarity_threshold=args.similarity_threshold,
        seed=args.seed,
    )

    # Print summary
    print(f"\n{'='*60}")
    print("Module 5: Cross-Sample Integration Complete")
    print(f"{'='*60}")
    print(result.summary())
    print(f"\nOutput: {args.output_dir}")


if __name__ == "__main__":
    main()
