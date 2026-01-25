#!/usr/bin/env python
"""
run_pipeline.py

Master orchestrator for MDK saturation analysis pipeline.
Runs all scripts in sequence, stopping on failure.

Usage:
    python run_pipeline.py              # Run full pipeline
    python run_pipeline.py --from 03    # Start from script 03
    python run_pipeline.py --only 02    # Run only script 02
"""

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent / "scripts"

SCRIPTS = [
    "01_summarize_spatial_finding.py",
    "02_analyze_chaperone_expression.py",
    "03_analyze_er_binding_changes.py",
    "04_quantify_saturation.py",
    "05_analyze_foxa1_perturbations.py",
    "06_cross_validate.py",
    "07_generate_report.py",
]


def run_script(script_name):
    """Run a single script and return success status."""
    script_path = SCRIPT_DIR / script_name
    print(f"\n{'='*60}")
    print(f"RUNNING: {script_name}")
    print('='*60)

    result = subprocess.run(
        ["python", str(script_path)],
        cwd=str(SCRIPT_DIR.parent)
    )

    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description="Run MDK saturation analysis pipeline")
    parser.add_argument("--from", dest="start_from", type=int,
                       help="Start from script number (e.g., --from 3)")
    parser.add_argument("--only", type=int,
                       help="Run only this script number (e.g., --only 2)")
    args = parser.parse_args()

    scripts_to_run = SCRIPTS

    if args.only:
        scripts_to_run = [s for s in SCRIPTS if s.startswith(f"{args.only:02d}_")]
        if not scripts_to_run:
            sys.exit(f"No script found with number {args.only}")
    elif args.start_from:
        scripts_to_run = [s for s in SCRIPTS if int(s[:2]) >= args.start_from]

    print("MDK SATURATION ANALYSIS PIPELINE")
    print("=" * 60)
    print(f"Scripts to run: {len(scripts_to_run)}")

    for script in scripts_to_run:
        success = run_script(script)
        if not success:
            print(f"\nPIPELINE FAILED at {script}")
            sys.exit(1)

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print("\nOutputs:")
    print("  - Tables: outputs/tables/")
    print("  - Figures: outputs/figures/")
    print("  - Report: outputs/report.md")

    return 0


if __name__ == "__main__":
    sys.exit(main())
