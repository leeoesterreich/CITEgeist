#!/usr/bin/env python
"""
Evaluate single-cell demonstration results.

Validates discovered profiles against expected RCC biology and runs
gene set enrichment analysis on discovered programs.

Usage:
    python evaluate_singlecell.py --mode full
    python evaluate_singlecell.py --mode quadrant --quadrant-id 0
    python evaluate_singlecell.py --mode all  # Aggregate all results
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell_demonstration"

# Expected profiles based on RCC literature
EXPECTED_PROFILES = {
    "T_helper": {"required": {"CD3E", "CD4"}, "optional": {"CD45RO"}},
    "T_cytotoxic": {"required": {"CD3E", "CD8A"}, "optional": {"GranzymeB"}},
    "Macrophage_M2": {"required": {"CD68", "CD163"}, "optional": {"HLA-DR"}},
    "B_cell": {"required": {"CD20"}, "optional": {"CD45RA"}},
    "Epithelial": {"required": {"PanCK"}, "optional": {"E-Cadherin", "Beta-catenin"}},
    "Stromal_CAF": {"required": {"Vimentin", "alphaSMA"}, "optional": set()},
    "Endothelial": {"required": {"CD31"}, "optional": {"Vimentin"}},
}


def match_profiles(
    discovered: Dict[str, Dict],
    expected: Dict[str, Dict],
) -> Dict:
    """
    Match discovered profiles to expected profiles.

    Returns metrics on profile recovery.
    """
    matches = []
    matched_expected = set()
    matched_discovered = set()

    for disc_name, disc_profile in discovered.items():
        disc_markers = set(disc_profile["markers"])

        best_match = None
        best_score = 0

        for exp_name, exp_profile in expected.items():
            required = exp_profile["required"]
            optional = exp_profile["optional"]

            # Check if required markers are present
            required_match = len(required & disc_markers) / len(required) if required else 0
            optional_match = len(optional & disc_markers) / len(optional) if optional else 0

            # Score: required markers count more
            score = 0.7 * required_match + 0.3 * optional_match

            if score > best_score and required_match >= 0.5:  # At least half of required
                best_score = score
                best_match = exp_name

        if best_match:
            matches.append({
                "discovered": disc_name,
                "expected": best_match,
                "score": best_score,
                "discovered_markers": list(disc_markers),
                "expected_required": list(expected[best_match]["required"]),
                "expected_optional": list(expected[best_match]["optional"]),
            })
            matched_expected.add(best_match)
            matched_discovered.add(disc_name)

    # Find unmatched
    unmatched_expected = set(expected.keys()) - matched_expected
    unmatched_discovered = set(discovered.keys()) - matched_discovered

    # Compute metrics
    precision = len(matched_discovered) / len(discovered) if discovered else 0
    recall = len(matched_expected) / len(expected) if expected else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {
        "matches": matches,
        "unmatched_expected": list(unmatched_expected),
        "unmatched_discovered": list(unmatched_discovered),
        "metrics": {
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "n_matched": len(matches),
            "n_expected": len(expected),
            "n_discovered": len(discovered),
        },
    }


def evaluate_programs(programs_dir: Path) -> Dict:
    """
    Evaluate discovered programs.

    For now, summarizes program statistics. GSEA can be added later.
    """
    results = {}

    for program_file in programs_dir.glob("*_programs.json"):
        cell_type = program_file.stem.replace("_programs", "")

        with open(program_file) as f:
            data = json.load(f)

        programs = data.get("programs", [])

        results[cell_type] = {
            "n_programs": len(programs),
            "n_cells": data.get("n_cells", 0),
            "programs": [
                {
                    "program_id": p["program_id"],
                    "top_5_genes": p["top_genes"][:5],
                    "morans_i": p["morans_i"],
                    "spatially_coherent": p["morans_i"] > 0.1,
                }
                for p in programs
            ],
            "mean_morans_i": np.mean([p["morans_i"] for p in programs]) if programs else 0,
            "n_spatially_coherent": sum(1 for p in programs if p["morans_i"] > 0.1),
        }

    return results


def generate_discovery_catalog(
    profile_eval: Dict,
    program_eval: Dict,
    output_dir: Path,
) -> str:
    """Generate markdown catalog of discoveries."""
    lines = [
        "# Single-Cell Discovery Catalog",
        "",
        "## Profile Discovery Summary",
        "",
        f"- **Profiles discovered:** {profile_eval['metrics']['n_discovered']}",
        f"- **Expected profiles recovered:** {profile_eval['metrics']['n_matched']}/{profile_eval['metrics']['n_expected']}",
        f"- **Precision:** {profile_eval['metrics']['precision']:.2f}",
        f"- **Recall:** {profile_eval['metrics']['recall']:.2f}",
        f"- **F1 Score:** {profile_eval['metrics']['f1']:.2f}",
        "",
        "### Matched Profiles",
        "",
    ]

    for match in profile_eval["matches"]:
        lines.append(f"- **{match['discovered']}** -> {match['expected']} (score: {match['score']:.2f})")
        lines.append(f"  - Markers: {', '.join(match['discovered_markers'])}")

    if profile_eval["unmatched_discovered"]:
        lines.extend([
            "",
            "### Novel Profiles (Not in Expected)",
            "",
        ])
        for name in profile_eval["unmatched_discovered"]:
            lines.append(f"- **{name}** - potential novel cell type or state")

    if profile_eval["unmatched_expected"]:
        lines.extend([
            "",
            "### Missing Expected Profiles",
            "",
        ])
        for name in profile_eval["unmatched_expected"]:
            lines.append(f"- {name}")

    lines.extend([
        "",
        "## Program Discovery Summary",
        "",
    ])

    for cell_type, data in program_eval.items():
        lines.extend([
            f"### {cell_type}",
            "",
            f"- Cells: {data['n_cells']:,}",
            f"- Programs: {data['n_programs']}",
            f"- Spatially coherent (Moran's I > 0.1): {data['n_spatially_coherent']}",
            "",
        ])

        for prog in data["programs"]:
            coherent = "Y" if prog["spatially_coherent"] else "x"
            lines.append(f"  - Program {prog['program_id']}: {', '.join(prog['top_5_genes'])} (I={prog['morans_i']:.3f} {coherent})")
        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Evaluate single-cell results")
    parser.add_argument(
        "--mode",
        choices=["full", "quadrant", "all"],
        required=True,
        help="Evaluate full, single quadrant, or aggregate all",
    )
    parser.add_argument(
        "--quadrant-id",
        type=int,
        choices=[0, 1, 2, 3],
        help="Quadrant ID (required if mode=quadrant)",
    )
    args = parser.parse_args()

    if args.mode == "quadrant" and args.quadrant_id is None:
        parser.error("--quadrant-id required when mode=quadrant")

    # Determine directories to evaluate
    if args.mode == "all":
        dirs_to_eval = [OUTPUT_BASE / "full"]
        dirs_to_eval.extend(OUTPUT_BASE / "quadrants" / f"Q{i}" for i in range(4))
        dirs_to_eval = [d for d in dirs_to_eval if d.exists()]
    elif args.mode == "full":
        dirs_to_eval = [OUTPUT_BASE / "full"]
    else:
        dirs_to_eval = [OUTPUT_BASE / "quadrants" / f"Q{args.quadrant_id}"]

    # Evaluate each directory
    all_results = {}
    for eval_dir in dirs_to_eval:
        if not eval_dir.exists():
            logger.warning(f"Directory not found: {eval_dir}")
            continue

        logger.info(f"Evaluating: {eval_dir}")

        # Load discovered profiles
        profiles_path = eval_dir / "module2c_profiles_selected.json"
        if not profiles_path.exists():
            logger.warning(f"No profiles found in {eval_dir}")
            continue

        with open(profiles_path) as f:
            discovered_profiles = json.load(f)

        # Evaluate profiles
        profile_eval = match_profiles(discovered_profiles, EXPECTED_PROFILES)

        # Evaluate programs
        programs_dir = eval_dir / "module4_programs"
        program_eval = evaluate_programs(programs_dir) if programs_dir.exists() else {}

        # Generate catalog
        catalog = generate_discovery_catalog(profile_eval, program_eval, eval_dir)

        # Save results
        eval_output_dir = eval_dir / "evaluation" if args.mode != "all" else OUTPUT_BASE / "evaluation"
        eval_output_dir.mkdir(parents=True, exist_ok=True)

        result_name = eval_dir.name if args.mode != "all" else "combined"

        with open(eval_output_dir / f"profile_validation_{result_name}.json", "w") as f:
            json.dump(profile_eval, f, indent=2)

        with open(eval_output_dir / f"program_evaluation_{result_name}.json", "w") as f:
            json.dump(program_eval, f, indent=2)

        with open(eval_output_dir / f"discovery_catalog_{result_name}.md", "w") as f:
            f.write(catalog)

        all_results[str(eval_dir)] = {
            "profile_eval": profile_eval,
            "program_eval": program_eval,
        }

        logger.info(f"Profile F1: {profile_eval['metrics']['f1']:.2f}")

    # Summary
    if args.mode == "all":
        summary_path = OUTPUT_BASE / "evaluation" / "validation_summary.json"
        with open(summary_path, "w") as f:
            json.dump(all_results, f, indent=2)
        logger.info(f"Summary saved to {summary_path}")

    logger.info("Evaluation complete!")


if __name__ == "__main__":
    main()
