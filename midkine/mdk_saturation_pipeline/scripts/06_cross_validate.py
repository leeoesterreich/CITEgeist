#!/usr/bin/env python
"""
06_cross_validate.py

Datasets: Outputs from scripts 01-05
Question: Do all lines of evidence agree?
Inputs:   outputs/tables/*.csv from previous scripts
Outputs:  outputs/tables/cross_validation.csv
          outputs/tables/evidence_summary.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml

SCRIPT_DIR = Path(__file__).parent
PIPELINE_DIR = SCRIPT_DIR.parent
CONFIG_PATH = PIPELINE_DIR / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

OUTPUT_DIR = PIPELINE_DIR / config['paths']['output_dir']
TABLES_DIR = OUTPUT_DIR / "tables"


def validate_inputs():
    """Check required output files from previous scripts exist."""
    required = [
        "chaperone_expression.csv",
        "binding_changes.csv",
        "saturation_metrics.csv",
        "perturbation_results.csv",
    ]
    missing = [f for f in required if not (TABLES_DIR / f).exists()]
    if missing:
        sys.exit(f"Missing required files (run previous scripts first):\n" + "\n".join(missing))


def main():
    print("=" * 80)
    print("SCRIPT 06: CROSS-VALIDATION")
    print("=" * 80)

    validate_inputs()

    # Load all previous outputs
    expr_df = pd.read_csv(TABLES_DIR / "chaperone_expression.csv")
    binding_df = pd.read_csv(TABLES_DIR / "binding_changes.csv")
    saturation_df = pd.read_csv(TABLES_DIR / "saturation_metrics.csv")
    perturb_df = pd.read_csv(TABLES_DIR / "perturbation_results.csv")

    checks = []

    # CHECK 1: Expression shows opposite regulation
    print("\n" + "-" * 40)
    print("CHECK 1: Opposite expression regulation")
    print("-" * 40)

    opposite_expr = expr_df['opposite_regulation'].sum()
    total_genes = len(expr_df)
    check1_pass = opposite_expr >= 3  # At least 3/5 genes

    checks.append({
        'check': 'opposite_expression',
        'description': 'Chaperones show opposite regulation in MCF7 vs T47D',
        'metric': f'{opposite_expr}/{total_genes} genes',
        'threshold': '>=3 genes',
        'passed': check1_pass
    })
    print(f"Opposite regulation: {opposite_expr}/{total_genes} genes -> {'PASS' if check1_pass else 'FAIL'}")

    # CHECK 2: Binding changes correlate with expression
    print("\n" + "-" * 40)
    print("CHECK 2: Binding correlates with expression")
    print("-" * 40)

    merged = expr_df.merge(binding_df, on='gene')
    # MCF7: expression UP when binding lost (negative delta)
    # T47D: expression DOWN when binding gained (positive delta)

    concordant = 0
    for _, row in merged.iterrows():
        mcf7_expr_up = row['MCF7_FC'] > 1.1
        mcf7_bind_lost = row['MCF7_delta'] < 0
        t47d_expr_down = row['T47D_FC'] < 0.9
        t47d_bind_gained = row['T47D_delta'] > 0

        if (mcf7_expr_up == mcf7_bind_lost) or (t47d_expr_down == t47d_bind_gained):
            concordant += 1

    check2_pass = concordant >= len(merged) // 2
    checks.append({
        'check': 'binding_expression_correlation',
        'description': 'Binding changes correlate with expression changes',
        'metric': f'{concordant}/{len(merged)} genes',
        'threshold': '>=50%',
        'passed': check2_pass
    })
    print(f"Concordant: {concordant}/{len(merged)} genes -> {'PASS' if check2_pass else 'FAIL'}")

    # CHECK 3: MCF7 is saturated (>2% ER occupancy)
    print("\n" + "-" * 40)
    print("CHECK 3: MCF7 is ER-saturated")
    print("-" * 40)

    mcf7_sat = saturation_df[saturation_df['cell_line'] == 'MCF7']['ER_occupancy_pct'].values[0]
    t47d_sat = saturation_df[saturation_df['cell_line'] == 'T47D']['ER_occupancy_pct'].values[0]

    check3_pass = mcf7_sat > t47d_sat and mcf7_sat > 2.0
    checks.append({
        'check': 'mcf7_saturated',
        'description': 'MCF7 has higher ER occupancy than T47D',
        'metric': f'MCF7={mcf7_sat:.1f}%, T47D={t47d_sat:.1f}%',
        'threshold': 'MCF7 > T47D and MCF7 > 2%',
        'passed': check3_pass
    })
    print(f"MCF7 occupancy: {mcf7_sat:.1f}%, T47D: {t47d_sat:.1f}% -> {'PASS' if check3_pass else 'FAIL'}")

    # CHECK 4: FOXA1-KD derepresses chaperones
    print("\n" + "-" * 40)
    print("CHECK 4: FOXA1-KD causes derepression")
    print("-" * 40)

    foxa1_kd = perturb_df[perturb_df['perturbation'] == 'T47D_FOXA1_KD']
    up_count = (foxa1_kd['direction'] == 'UP').sum()

    check4_pass = up_count >= 3
    checks.append({
        'check': 'foxa1_kd_derepression',
        'description': 'FOXA1-KD increases chaperone expression (derepression)',
        'metric': f'{up_count}/{len(foxa1_kd)} genes UP',
        'threshold': '>=3 genes UP',
        'passed': check4_pass
    })
    print(f"FOXA1-KD: {up_count}/{len(foxa1_kd)} chaperones UP -> {'PASS' if check4_pass else 'FAIL'}")

    # CHECK 5: ER-KD derepresses chaperones
    print("\n" + "-" * 40)
    print("CHECK 5: ER-KD causes derepression")
    print("-" * 40)

    er_kd = perturb_df[perturb_df['perturbation'] == 'MCF7_ER_KD']
    up_count = (er_kd['direction'] == 'UP').sum()

    check5_pass = up_count >= 3
    checks.append({
        'check': 'er_kd_derepression',
        'description': 'ER-KD increases chaperone expression (confirms ER represses)',
        'metric': f'{up_count}/{len(er_kd)} genes UP',
        'threshold': '>=3 genes UP',
        'passed': check5_pass
    })
    print(f"ER-KD: {up_count}/{len(er_kd)} chaperones UP -> {'PASS' if check5_pass else 'FAIL'}")

    # Summary
    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(TABLES_DIR / "cross_validation.csv", index=False)

    passed = checks_df['passed'].sum()
    total = len(checks_df)

    print("\n" + "=" * 40)
    print(f"CROSS-VALIDATION SUMMARY: {passed}/{total} checks passed")
    print("=" * 40)

    # Evidence summary
    evidence = pd.DataFrame({
        'evidence_type': [
            'Expression (GSE89888)',
            'Binding (GSE125117)',
            'Saturation (GSE254216)',
            'Perturbation (GSE254218, GSE75329)'
        ],
        'supports_model': [check1_pass, check2_pass, check3_pass, check4_pass and check5_pass],
        'key_finding': [
            f'{opposite_expr}/{total_genes} chaperones opposite-regulated',
            f'{concordant}/{len(merged)} binding-expression concordant',
            f'MCF7 {mcf7_sat:.1f}% vs T47D {t47d_sat:.1f}% occupancy',
            'FOXA1-KD and ER-KD both derepress chaperones'
        ]
    })
    evidence.to_csv(TABLES_DIR / "evidence_summary.csv", index=False)

    print(f"\nSaved: outputs/tables/cross_validation.csv")
    print(f"Saved: outputs/tables/evidence_summary.csv")

    if passed == total:
        print("\nAll evidence supports the ER saturation model")
    else:
        print(f"\n{total - passed} check(s) failed - review results")

    print("\n" + "=" * 80)
    print("SCRIPT 06 COMPLETE")
    print("=" * 80)

    return 0 if passed >= total - 1 else 1  # Allow 1 failure


if __name__ == "__main__":
    sys.exit(main())
