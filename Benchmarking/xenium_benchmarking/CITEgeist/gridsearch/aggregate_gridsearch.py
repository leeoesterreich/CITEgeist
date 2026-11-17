#!/usr/bin/env python -u
"""Aggregate grid search results: robustness ranking + sensitivity analysis.

Usage:
    python aggregate_gridsearch.py --stage props --results-dir results/props/ --grid-file combos_props.csv
    python aggregate_gridsearch.py --stage sace --results-dir results/sace/ --grid-file combos_sace.csv
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

N_REGIONS = 5
N_BOOTSTRAP = 1000
INSENSITIVE_THRESHOLD = 0.01
INTERACTION_THRESHOLD = 0.02

STAGE_CONFIG = {
    "props": {
        "primary_metric": "overall_pearson_r",
        "secondary_metric": "mean_per_type_r",
        "param_cols": ["alpha_elastic", "lambda_reg", "lambda_laplacian", "detection_gate_ub", "nuclei_prior_lambda"],
    },
    "sace": {
        "primary_metric": "spot_layer_overall_r",
        "secondary_metric": "percell_mean_gene_r",
        "param_cols": ["max_iter", "n_0", "bandwidth", "damping_eta", "eps"],
    },
}


def load_results(results_dir: Path, grid_df: pd.DataFrame):
    """Load all result JSONs into a DataFrame."""
    rows = []
    missing = []
    corrupt = []
    for combo_id in range(len(grid_df)):
        for region_id in range(N_REGIONS):
            path = results_dir / str(combo_id) / f"region_{region_id}.json"
            if not path.exists():
                missing.append((combo_id, region_id))
                continue
            try:
                with open(path) as f:
                    data = json.load(f)
                row = {"combo_id": combo_id, "region_id": region_id}
                row.update(data.get("params", {}))
                metrics = data.get("metrics", {})
                row.update(metrics)
                if "per_type_r" in metrics and isinstance(metrics["per_type_r"], dict):
                    row["mean_per_type_r"] = float(np.mean(list(metrics["per_type_r"].values())))
                row["runtime_seconds"] = data.get("runtime_seconds", 0)
                rows.append(row)
            except (json.JSONDecodeError, KeyError) as e:
                corrupt.append((combo_id, region_id, str(e)))

    if missing:
        logger.warning("%d missing result files", len(missing))
    if corrupt:
        logger.warning("%d corrupt result files", len(corrupt))

    return pd.DataFrame(rows), missing, corrupt


def compute_robustness(df: pd.DataFrame, metric: str, penalty: float = 1.0):
    """Compute robustness score per combo: mean - penalty * std."""
    agg = df.groupby("combo_id")[metric].agg(["mean", "std"]).reset_index()
    agg.columns = ["combo_id", f"{metric}_mean", f"{metric}_std"]
    agg[f"{metric}_std"] = agg[f"{metric}_std"].fillna(0)
    agg[f"{metric}_robustness"] = agg[f"{metric}_mean"] - penalty * agg[f"{metric}_std"]
    return agg


def sensitivity_analysis(df: pd.DataFrame, param_cols: list, metric: str, output_dir: Path):
    """Compute marginal effects and bootstrap CIs for each parameter."""
    results = []
    for param in param_cols:
        marginals = df.groupby(param)[metric].mean()
        marginal_range = float(marginals.max() - marginals.min())
        best_val = marginals.idxmax()

        boot_diffs = []
        for _ in range(N_BOOTSTRAP):
            boot_idx = np.random.choice(len(df), size=len(df), replace=True)
            boot_df = df.iloc[boot_idx]
            boot_marginals = boot_df.groupby(param)[metric].mean()
            if len(boot_marginals) >= 2:
                boot_diffs.append(float(boot_marginals.max() - boot_marginals.min()))
        ci_lower = float(np.percentile(boot_diffs, 2.5)) if boot_diffs else 0.0
        ci_upper = float(np.percentile(boot_diffs, 97.5)) if boot_diffs else 0.0

        insensitive = marginal_range < INSENSITIVE_THRESHOLD and ci_lower <= 0.0

        results.append(
            {
                "parameter": param,
                "marginal_range": marginal_range,
                "best_value": best_val,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
                "insensitive": insensitive,
                "marginal_values": {str(k): float(v) for k, v in marginals.items()},
            }
        )

    sens_df = pd.DataFrame(results)
    sens_df.to_csv(output_dir / "sensitivity.csv", index=False)
    return sens_df


def interaction_analysis(df: pd.DataFrame, param_cols: list, metric: str, output_dir: Path):
    """Compute pairwise interaction effects."""
    results = []
    for i, p1 in enumerate(param_cols):
        for p2 in param_cols[i + 1 :]:
            grand_mean = df[metric].mean()
            interaction_effects = []
            for v1 in df[p1].unique():
                for v2 in df[p2].unique():
                    mask = (df[p1] == v1) & (df[p2] == v2)
                    if mask.sum() == 0:
                        continue
                    cell_mean = df.loc[mask, metric].mean()
                    row_mean = df.loc[df[p1] == v1, metric].mean()
                    col_mean = df.loc[df[p2] == v2, metric].mean()
                    interaction = cell_mean - row_mean - col_mean + grand_mean
                    interaction_effects.append(abs(interaction))

            max_interaction = max(interaction_effects) if interaction_effects else 0.0
            if max_interaction > INTERACTION_THRESHOLD:
                results.append(
                    {
                        "param_1": p1,
                        "param_2": p2,
                        "max_interaction_effect": max_interaction,
                    }
                )

    inter_df = pd.DataFrame(results)
    inter_df.to_csv(output_dir / "interactions.csv", index=False)
    return inter_df


def main():
    parser = argparse.ArgumentParser(description="Aggregate grid search results")
    parser.add_argument("--stage", choices=["props", "sace"], required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--grid-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.results_dir.parent
    args.output_dir.mkdir(parents=True, exist_ok=True)

    config = STAGE_CONFIG[args.stage]
    grid_df = pd.read_csv(args.grid_file)

    df, missing, corrupt = load_results(args.results_dir, grid_df)
    logger.info("Loaded %d results (%d missing, %d corrupt)", len(df), len(missing), len(corrupt))

    if len(df) == 0:
        logger.error("No results found!")
        sys.exit(1)

    combo_counts = df.groupby("combo_id").size()
    complete_combos = combo_counts[combo_counts == N_REGIONS].index
    incomplete = combo_counts[combo_counts < N_REGIONS]
    if len(incomplete) > 0:
        logger.warning("%d incomplete combos excluded from ranking", len(incomplete))
    df_complete = df[df["combo_id"].isin(complete_combos)].copy()

    primary = config["primary_metric"]
    secondary = config["secondary_metric"]
    rob_primary = compute_robustness(df_complete, primary)
    rob_secondary = compute_robustness(df_complete, secondary)

    summary = rob_primary.merge(rob_secondary, on="combo_id")
    summary = summary.merge(grid_df, on="combo_id")
    summary = summary.sort_values(
        [f"{primary}_robustness", f"{secondary}_robustness"],
        ascending=False,
    ).reset_index(drop=True)
    summary.to_csv(args.output_dir / f"{args.stage}_summary.csv", index=False)

    best_id = int(summary.iloc[0]["combo_id"])
    best_row = summary.iloc[0].to_dict()
    best_row = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in best_row.items()}
    best_regions = df_complete[df_complete["combo_id"] == best_id].to_dict("records")
    best_config = {"best_combo": best_row, "per_region": best_regions}
    with open(args.output_dir / f"{args.stage}_best_config.json", "w") as f:
        json.dump(best_config, f, indent=2, default=str)

    sens_df = sensitivity_analysis(df_complete, config["param_cols"], primary, args.output_dir)
    inter_df = interaction_analysis(df_complete, config["param_cols"], primary, args.output_dir)

    insensitive = sens_df[sens_df["insensitive"]]["parameter"].tolist()
    with open(args.output_dir / f"{args.stage}_insensitive_params.txt", "w") as f:
        for p in insensitive:
            row = sens_df[sens_df["parameter"] == p].iloc[0]
            f.write(
                f"{p}: marginal_range={row['marginal_range']:.4f}, "
                f"CI=[{row['ci_lower']:.4f}, {row['ci_upper']:.4f}]\n"
            )

    print(f"\n{'='*70}")
    print(f"GRID SEARCH RESULTS: {args.stage.upper()}")
    print(f"{'='*70}")
    print(f"Total combos: {len(grid_df)}")
    print(f"Complete combos: {len(complete_combos)} ({len(incomplete)} incomplete)")
    print(f"\nTOP 10 CONFIGS (by {primary} robustness):")
    print(f"{'='*70}")
    top10 = summary.head(10)
    for _, row in top10.iterrows():
        param_str = ", ".join(f"{p}={row[p]}" for p in config["param_cols"])
        print(
            f"  combo {int(row['combo_id']):4d}: "
            f"mean={row[f'{primary}_mean']:.4f} "
            f"std={row[f'{primary}_std']:.4f} "
            f"rob={row[f'{primary}_robustness']:.4f} "
            f"| {param_str}"
        )

    print(f"\nSENSITIVITY ANALYSIS:")
    print(f"{'='*70}")
    for _, row in sens_df.iterrows():
        flag = " ** INSENSITIVE **" if row["insensitive"] else ""
        print(
            f"  {row['parameter']:25s}: range={row['marginal_range']:.4f} "
            f"CI=[{row['ci_lower']:.4f}, {row['ci_upper']:.4f}] "
            f"best={row['best_value']}{flag}"
        )

    if len(inter_df) > 0:
        print(f"\nSIGNIFICANT INTERACTIONS (>{INTERACTION_THRESHOLD:.2f}):")
        for _, row in inter_df.iterrows():
            print(f"  {row['param_1']} x {row['param_2']}: " f"max_effect={row['max_interaction_effect']:.4f}")

    if insensitive:
        print(f"\nINSENSITIVE PARAMETERS (candidates for removal):")
        for p in insensitive:
            print(f"  - {p}")
    else:
        print(f"\nNo parameters flagged as insensitive (all marginal ranges >= {INSENSITIVE_THRESHOLD})")


if __name__ == "__main__":
    main()
