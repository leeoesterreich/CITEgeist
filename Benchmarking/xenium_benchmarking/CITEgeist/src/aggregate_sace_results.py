#!/usr/bin/env python
"""Aggregate SACE benchmark results across all 5 Xenium pseudo-Visium regions.

Reads per-region JSON files produced by run_sace_benchmark.py from
results_sace/sace_region_{0-4}.json and prints:

  1. Summary table: Method | percell_r (mean±std) | spot_layer_r (mean±std)
  2. Per-type breakdown for sace_full (percell_per_type_r)
  3. SACE improvement over each baseline (delta ± std)
  4. Saves aggregated JSON to results_sace/sace_aggregate.json

Missing region files are reported and skipped; stats are computed over
available regions only.  Handles two JSON formats:

  - New format (run_sace_benchmark.py): top-level method keys
    ``equal_split``, ``type_reference``, ``sace_iter0``, ``sace_full``,
    each with ``percell_mean_gene_r``, ``percell_median_gene_r``,
    ``spot_layer_overall_r``, ``percell_per_type_r``, ``spot_layer_per_type_r``.

  - Old format (single-method): top-level keys ``mean_gene_r``,
    ``median_gene_r``, ``per_type_mean_r`` (treated as ``sace_full`` only).

Dependencies: json, pathlib, numpy only.
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import json
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RESULTS_DIR = Path(
    str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "results_sace")
)
N_REGIONS = 5
METHODS = ["equal_split", "type_reference", "sace_iter0", "sace_full"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fmt(values):
    """Return 'mean±std' string from a list of floats, or 'N/A' if empty."""
    if not values:
        return "N/A"
    arr = np.array(values, dtype=float)
    return f"{arr.mean():.4f}±{arr.std():.4f}"


def _normalise_record(data, region_id):
    """Return a dict mapping method_name -> metric_dict.

    Handles both old (single sace_full) and new (multi-method) JSON layouts.

    Args:
        data: Parsed JSON dict.
        region_id: Integer region identifier (for logging only).

    Returns:
        Dict[str, dict] with at least the method keys present in ``data``.
    """
    # Detect old format by absence of any of the method keys
    if not any(m in data for m in METHODS):
        # Old format: keys are percell metrics directly
        record = {
            "sace_full": {
                "percell_mean_gene_r": data.get("mean_gene_r", float("nan")),
                "percell_median_gene_r": data.get("median_gene_r", float("nan")),
                "spot_layer_overall_r": float("nan"),
                "percell_per_type_r": data.get("per_type_mean_r", {}),
                "spot_layer_per_type_r": {},
            }
        }
        return record

    # New format
    record = {}
    for method in METHODS:
        if method in data:
            md = data[method]
            record[method] = {
                "percell_mean_gene_r": md.get("percell_mean_gene_r", float("nan")),
                "percell_median_gene_r": md.get("percell_median_gene_r", float("nan")),
                "spot_layer_overall_r": md.get("spot_layer_overall_r", float("nan")),
                "percell_per_type_r": md.get("percell_per_type_r", {}),
                "spot_layer_per_type_r": md.get("spot_layer_per_type_r", {}),
            }
    return record


# ---------------------------------------------------------------------------
# Load all regions
# ---------------------------------------------------------------------------


def load_all_regions():
    """Load and normalise results for regions 0-4.

    Returns:
        Tuple of (region_records, missing_regions) where region_records is a
        list of (region_id, normalised_dict) and missing_regions is a list of
        ints for regions whose JSON files were not found.
    """
    region_records = []
    missing_regions = []

    for rid in range(N_REGIONS):
        path = RESULTS_DIR / f"sace_region_{rid}.json"
        if not path.exists():
            print(f"[WARNING] Missing region file: {path}")
            missing_regions.append(rid)
            continue
        with open(path) as fh:
            data = json.load(fh)
        record = _normalise_record(data, rid)
        region_records.append((rid, record))

    return region_records, missing_regions


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def aggregate_metrics(region_records):
    """Aggregate per-region metrics into per-method lists.

    Args:
        region_records: List of (region_id, normalised_dict) tuples.

    Returns:
        Dict[method -> Dict[metric -> list_of_floats_across_regions]].
    """
    agg = {
        m: {
            "percell_mean_gene_r": [],
            "percell_median_gene_r": [],
            "spot_layer_overall_r": [],
            "percell_per_type_r": {},  # type -> list
            "spot_layer_per_type_r": {},  # type -> list
        }
        for m in METHODS
    }

    for _rid, record in region_records:
        for method, metrics in record.items():
            if method not in agg:
                continue
            bucket = agg[method]
            for scalar_key in ("percell_mean_gene_r", "percell_median_gene_r", "spot_layer_overall_r"):
                val = metrics.get(scalar_key, float("nan"))
                if np.isfinite(val):
                    bucket[scalar_key].append(val)

            # Per-type percell
            for ct, val in metrics.get("percell_per_type_r", {}).items():
                if np.isfinite(val):
                    bucket["percell_per_type_r"].setdefault(ct, []).append(val)

            # Per-type spot layer
            for ct, val in metrics.get("spot_layer_per_type_r", {}).items():
                if np.isfinite(val):
                    bucket["spot_layer_per_type_r"].setdefault(ct, []).append(val)

    return agg


# ---------------------------------------------------------------------------
# Print reports
# ---------------------------------------------------------------------------


def print_summary_table(agg):
    """Print summary table: Method | percell_r (mean±std) | spot_layer_r (mean±std).

    Args:
        agg: Aggregated metrics dict from ``aggregate_metrics``.
    """
    header = f"{'Method':<18}  {'percell_r (mean±std)':>22}  {'spot_layer_r (mean±std)':>24}"
    sep = "-" * len(header)
    print("\n=== Summary Table ===")
    print(header)
    print(sep)
    for method in METHODS:
        bucket = agg[method]
        pc = _fmt(bucket["percell_mean_gene_r"])
        sl = _fmt(bucket["spot_layer_overall_r"])
        n = len(bucket["percell_mean_gene_r"])
        print(f"{method:<18}  {pc:>22}  {sl:>24}   (n={n} regions)")
    print()


def print_pertype_breakdown(agg):
    """Print per-type percell_r breakdown for sace_full.

    Args:
        agg: Aggregated metrics dict from ``aggregate_metrics``.
    """
    print("=== Per-Type Breakdown (sace_full, percell_per_type_r) ===")
    per_type = agg["sace_full"]["percell_per_type_r"]
    if not per_type:
        print("  No per-type data available.\n")
        return
    cell_types = sorted(per_type.keys())
    header = f"  {'Cell Type':<22}  {'mean±std':>18}  {'n regions':>10}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for ct in cell_types:
        vals = per_type[ct]
        print(f"  {ct:<22}  {_fmt(vals):>18}  {len(vals):>10}")
    print()


def print_delta_table(agg):
    """Print SACE improvement (delta) over each baseline.

    Reports delta = sace_full - baseline for percell_r and spot_layer_r,
    computed per region then summarised as mean±std.

    Args:
        agg: Aggregated metrics dict from ``aggregate_metrics``.
    """
    print("=== SACE Improvement Over Baselines (percell_mean_gene_r) ===")
    sace_vals = agg["sace_full"]["percell_mean_gene_r"]
    if not sace_vals:
        print("  No sace_full data available.\n")
        return

    # We need per-region deltas, so re-use the scalar lists indexed by
    # position (both lists are in the same region order if n equal)
    baselines = [m for m in METHODS if m != "sace_full"]
    header = f"  {'Baseline':<18}  {'delta percell_r':>18}  {'delta spot_layer_r':>20}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for baseline in baselines:
        b_pc = agg[baseline]["percell_mean_gene_r"]
        s_pc = agg["sace_full"]["percell_mean_gene_r"]
        b_sl = agg[baseline]["spot_layer_overall_r"]
        s_sl = agg["sace_full"]["spot_layer_overall_r"]

        n_pc = min(len(b_pc), len(s_pc))
        n_sl = min(len(b_sl), len(s_sl))

        if n_pc > 0:
            delta_pc_vals = np.array(s_pc[:n_pc]) - np.array(b_pc[:n_pc])
            delta_pc_str = _fmt(delta_pc_vals.tolist())
        else:
            delta_pc_str = "N/A"

        if n_sl > 0:
            delta_sl_vals = np.array(s_sl[:n_sl]) - np.array(b_sl[:n_sl])
            delta_sl_str = _fmt(delta_sl_vals.tolist())
        else:
            delta_sl_str = "N/A"

        print(f"  {baseline:<18}  {delta_pc_str:>18}  {delta_sl_str:>20}")
    print()


# ---------------------------------------------------------------------------
# Build aggregate JSON
# ---------------------------------------------------------------------------


def build_aggregate_json(region_records, missing_regions, agg):
    """Build the aggregated result dict for JSON serialisation.

    Args:
        region_records: List of (region_id, normalised_dict) tuples.
        missing_regions: List of missing region IDs.
        agg: Aggregated metrics dict from ``aggregate_metrics``.

    Returns:
        Dict suitable for json.dump.
    """

    def _summarise(values):
        if not values:
            return {"mean": None, "std": None, "n": 0, "values": []}
        arr = np.array(values, dtype=float)
        return {
            "mean": float(arr.mean()),
            "std": float(arr.std()),
            "n": len(values),
            "values": [float(v) for v in values],
        }

    methods_summary = {}
    for method in METHODS:
        bucket = agg[method]
        per_type_pc = {ct: _summarise(vals) for ct, vals in bucket["percell_per_type_r"].items()}
        per_type_sl = {ct: _summarise(vals) for ct, vals in bucket["spot_layer_per_type_r"].items()}
        methods_summary[method] = {
            "percell_mean_gene_r": _summarise(bucket["percell_mean_gene_r"]),
            "percell_median_gene_r": _summarise(bucket["percell_median_gene_r"]),
            "spot_layer_overall_r": _summarise(bucket["spot_layer_overall_r"]),
            "percell_per_type_r": per_type_pc,
            "spot_layer_per_type_r": per_type_sl,
        }

    return {
        "n_regions_total": N_REGIONS,
        "regions_loaded": [rid for rid, _ in region_records],
        "regions_missing": missing_regions,
        "methods": methods_summary,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    """Load region JSONs, print reports, and save aggregated JSON."""
    print(f"Loading SACE results from: {RESULTS_DIR}")
    region_records, missing_regions = load_all_regions()

    if not region_records:
        print("[ERROR] No region files found. Nothing to aggregate.", file=sys.stderr)
        sys.exit(1)

    print(
        f"Loaded {len(region_records)}/{N_REGIONS} regions "
        f"(missing: {missing_regions if missing_regions else 'none'})"
    )

    agg = aggregate_metrics(region_records)

    print_summary_table(agg)
    print_pertype_breakdown(agg)
    print_delta_table(agg)

    # Save aggregated JSON
    out_path = RESULTS_DIR / "sace_aggregate.json"
    aggregate = build_aggregate_json(region_records, missing_regions, agg)
    with open(out_path, "w") as fh:
        json.dump(aggregate, fh, indent=2)
    print(f"Aggregated results saved to: {out_path}")


if __name__ == "__main__":
    main()
