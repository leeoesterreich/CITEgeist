#!/usr/bin/env python
"""Generate parameter grid CSVs for CITEgeist Xenium grid search."""

import itertools
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent

PROPS_GRID = {
    "alpha_elastic": [0.0, 0.3, 0.5, 0.7, 1.0],
    "lambda_reg": [0.1, 0.5, 1.0, 2.0, 5.0],
    "lambda_laplacian": [0.0, 0.01, 0.03, 0.1, 0.3],
    "detection_gate_ub": [0.005, 0.01, 0.02, 0.05],
    "nuclei_prior_lambda": [0.0, 0.5, 1.0, 2.0, 5.0],
}

SACE_GRID = {
    "max_iter": [1, 2, 5, 10, 20],
    "n_0": [1.0, 5.0, 10.0, 20.0, 50.0],
    "bandwidth": ["None", "50", "100", "200", "500"],
    "damping_eta": [0.5, 0.7, 1.0],
    "eps": [1e-8, 1e-6, 1e-4],
}


def generate_combos(grid: dict, output_path: Path) -> pd.DataFrame:
    keys = list(grid.keys())
    values = list(grid.values())
    rows = []
    for combo_id, combo in enumerate(itertools.product(*values)):
        row = {"combo_id": combo_id}
        for k, v in zip(keys, combo):
            row[k] = v
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"Wrote {len(df)} combos to {output_path}")
    return df


def main():
    props_df = generate_combos(PROPS_GRID, SCRIPT_DIR / "combos_props.csv")
    sace_df = generate_combos(SACE_GRID, SCRIPT_DIR / "combos_sace.csv")
    print(f"Props: {len(props_df)} combos ({' x '.join(str(len(v)) for v in PROPS_GRID.values())})")
    print(f"SACE:  {len(sace_df)} combos ({' x '.join(str(len(v)) for v in SACE_GRID.values())})")


if __name__ == "__main__":
    main()
