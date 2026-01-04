#!/bin/bash
#SBATCH --job-name=xenium_integrate      # Job name
#SBATCH --output=slurm_log/%x_%j.out     # Standard output
#SBATCH --error=slurm_log/%x_%j.err      # Standard error
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=4                # CPU cores per task
#SBATCH --mem=64G                        # Memory per node
#SBATCH --time=02:00:00                  # Time limit
#SBATCH --cluster=htc                    # Cluster
#SBATCH --partition=htc                  # Partition
# Dependency removed - run manually after array job completes

# ============================================================================
# Module 5: Cross-Sample Integration + Final Benchmarking Report
# ============================================================================
# Run this AFTER all regions have been processed
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
OUTPUT_DIR="${CITEGEIST_DIR}/output"
METRICS_DIR="${BENCH_DIR}/metrics"

mkdir -p "${METRICS_DIR}"
mkdir -p "${CITEGEIST_DIR}/slurm/slurm_log"

echo "============================================"
echo "Module 5: Cross-Sample Integration"
echo "============================================"
echo "Start time: $(date)"
echo ""

# Activate environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

# Run integration and aggregation
python << 'EOF'
import sys
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

OUTPUT_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output"
METRICS_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/metrics"

# ============================================================================
# Aggregate Results from All Regions
# ============================================================================
print("\n" + "="*60)
print("Aggregating Results from All Regions")
print("="*60)

all_results = []
for region_id in range(5):
    result_file = OUTPUT_DIR / f"Xenium_region_{region_id}" / "pipeline_results.json"
    if result_file.exists():
        with open(result_file) as f:
            results = json.load(f)
            all_results.append(results)
            print(f"Region {region_id}: Loaded")
    else:
        print(f"Region {region_id}: NOT FOUND")

if not all_results:
    print("ERROR: No results found!")
    sys.exit(1)

# ============================================================================
# Aggregate Proportion Metrics
# ============================================================================
print("\n" + "="*60)
print("Proportion Deconvolution Metrics")
print("="*60)

prop_metrics = [r["metrics"]["proportions"] for r in all_results if "proportions" in r.get("metrics", {})]

if prop_metrics:
    agg_props = {
        "JSD_median": {
            "mean": np.mean([m["JSD_median"] for m in prop_metrics]),
            "std": np.std([m["JSD_median"] for m in prop_metrics]),
            "values": [m["JSD_median"] for m in prop_metrics],
        },
        "RMSE_global": {
            "mean": np.mean([m["RMSE_global"] for m in prop_metrics]),
            "std": np.std([m["RMSE_global"] for m in prop_metrics]),
            "values": [m["RMSE_global"] for m in prop_metrics],
        },
        "MAE_global": {
            "mean": np.mean([m["MAE_global"] for m in prop_metrics]),
            "std": np.std([m["MAE_global"] for m in prop_metrics]),
            "values": [m["MAE_global"] for m in prop_metrics],
        },
        "Pearson_r": {
            "mean": np.mean([m["Pearson_r"] for m in prop_metrics]),
            "std": np.std([m["Pearson_r"] for m in prop_metrics]),
            "values": [m["Pearson_r"] for m in prop_metrics],
        },
    }

    print(f"\nCell Proportion Benchmarking (n={len(prop_metrics)} regions):")
    print(f"  JSD (median):  {agg_props['JSD_median']['mean']:.4f} ± {agg_props['JSD_median']['std']:.4f}")
    print(f"  RMSE (global): {agg_props['RMSE_global']['mean']:.4f} ± {agg_props['RMSE_global']['std']:.4f}")
    print(f"  MAE (global):  {agg_props['MAE_global']['mean']:.4f} ± {agg_props['MAE_global']['std']:.4f}")
    print(f"  Pearson r:     {agg_props['Pearson_r']['mean']:.4f} ± {agg_props['Pearson_r']['std']:.4f}")

# ============================================================================
# Aggregate GEX Metrics
# ============================================================================
print("\n" + "="*60)
print("Gene Expression Deconvolution Metrics")
print("="*60)

gex_metrics = [r["metrics"]["gex"] for r in all_results if "gex" in r.get("metrics", {}) and "error" not in r["metrics"]["gex"]]

if gex_metrics:
    agg_gex = {
        "RMSE_mean": {
            "mean": np.mean([m["RMSE_mean"] for m in gex_metrics]),
            "std": np.std([m["RMSE_mean"] for m in gex_metrics]),
            "values": [m["RMSE_mean"] for m in gex_metrics],
        },
        "NRMSE_mean": {
            "mean": np.mean([m["NRMSE_mean"] for m in gex_metrics if not np.isnan(m["NRMSE_mean"])]),
            "std": np.std([m["NRMSE_mean"] for m in gex_metrics if not np.isnan(m["NRMSE_mean"])]),
            "values": [m["NRMSE_mean"] for m in gex_metrics],
        },
        "MAE_mean": {
            "mean": np.mean([m["MAE_mean"] for m in gex_metrics]),
            "std": np.std([m["MAE_mean"] for m in gex_metrics]),
            "values": [m["MAE_mean"] for m in gex_metrics],
        },
        "Pearson_r_mean": {
            "mean": np.mean([m["Pearson_r_mean"] for m in gex_metrics if not np.isnan(m["Pearson_r_mean"])]),
            "std": np.std([m["Pearson_r_mean"] for m in gex_metrics if not np.isnan(m["Pearson_r_mean"])]),
            "values": [m["Pearson_r_mean"] for m in gex_metrics],
        },
    }

    print(f"\nGEX Deconvolution Benchmarking (n={len(gex_metrics)} regions):")
    print(f"  RMSE (mean):   {agg_gex['RMSE_mean']['mean']:.4f} ± {agg_gex['RMSE_mean']['std']:.4f}")
    print(f"  NRMSE (mean):  {agg_gex['NRMSE_mean']['mean']:.4f} ± {agg_gex['NRMSE_mean']['std']:.4f}")
    print(f"  MAE (mean):    {agg_gex['MAE_mean']['mean']:.4f} ± {agg_gex['MAE_mean']['std']:.4f}")
    print(f"  Pearson r:     {agg_gex['Pearson_r_mean']['mean']:.4f} ± {agg_gex['Pearson_r_mean']['std']:.4f}")
else:
    agg_gex = None
    print("\nNo GEX metrics available (ground truth may not be generated)")

# ============================================================================
# Timing Summary
# ============================================================================
print("\n" + "="*60)
print("Timing Summary")
print("="*60)

for r in all_results:
    region_id = r["region_id"]
    timings = r.get("timings", {})
    total = r.get("total_time", 0)
    print(f"Region {region_id}: {total:.1f}s ({total/60:.1f} min)")
    for step, t in timings.items():
        print(f"  - {step}: {t:.1f}s")

total_time = sum(r.get("total_time", 0) for r in all_results)
print(f"\nTotal pipeline time: {total_time:.1f}s ({total_time/60:.1f} min)")

# ============================================================================
# Module 5: Cross-Sample Integration
# ============================================================================
print("\n" + "="*60)
print("Module 5: Cross-Sample Integration")
print("="*60)

try:
    from model import (
        load_multi_sample_results,
        integrate_samples,
        integrate_programs_harmony,
        match_programs_across_samples,
        save_integration_results,
    )

    sample_dirs = [OUTPUT_DIR / f"Xenium_region_{i}" for i in range(5)]
    sample_dirs = [d for d in sample_dirs if d.exists()]

    print(f"Found {len(sample_dirs)} sample directories")

    if len(sample_dirs) >= 2:
        integration_dir = OUTPUT_DIR / "integration"
        integration_dir.mkdir(exist_ok=True)

        # Try to load and integrate
        try:
            integration_result = integrate_samples(
                sample_dirs=[str(d) for d in sample_dirs],
                output_dir=str(integration_dir),
            )
            print("Integration completed successfully")

            if hasattr(integration_result, 'conserved_programs'):
                print(f"Conserved programs: {len(integration_result.conserved_programs)}")
        except Exception as e:
            print(f"Integration step failed: {e}")
            print("(This may be expected if Module 4 results are not available)")
    else:
        print("Not enough samples for integration")

except ImportError as e:
    print(f"Module 5 imports not available: {e}")
except Exception as e:
    print(f"Module 5 failed: {e}")

# ============================================================================
# Save Final Report
# ============================================================================
print("\n" + "="*60)
print("Saving Final Report")
print("="*60)

final_report = {
    "n_regions": len(all_results),
    "total_spots": sum(r.get("n_spots", 0) for r in all_results),
    "proportion_metrics": agg_props if prop_metrics else None,
    "gex_metrics": agg_gex if gex_metrics else None,
    "timing": {
        "total_seconds": total_time,
        "total_minutes": total_time / 60,
        "per_region": {r["region_id"]: r.get("total_time", 0) for r in all_results},
    },
    "regions": all_results,
}

with open(METRICS_DIR / "final_benchmark_report.json", "w") as f:
    def convert(obj):
        if isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj
    json.dump(convert(final_report), f, indent=2)

# Create summary CSV
summary_df = pd.DataFrame([
    {
        "region_id": r["region_id"],
        "n_spots": r.get("n_spots", 0),
        "n_genes": r.get("n_genes_filtered", 0),
        # Proportion metrics
        "prop_JSD_median": r["metrics"]["proportions"]["JSD_median"] if "proportions" in r.get("metrics", {}) else None,
        "prop_RMSE_global": r["metrics"]["proportions"]["RMSE_global"] if "proportions" in r.get("metrics", {}) else None,
        "prop_MAE_global": r["metrics"]["proportions"]["MAE_global"] if "proportions" in r.get("metrics", {}) else None,
        "prop_Pearson_r": r["metrics"]["proportions"]["Pearson_r"] if "proportions" in r.get("metrics", {}) else None,
        # GEX metrics
        "gex_RMSE_mean": r["metrics"]["gex"]["RMSE_mean"] if "gex" in r.get("metrics", {}) and "error" not in r["metrics"]["gex"] else None,
        "gex_NRMSE_mean": r["metrics"]["gex"]["NRMSE_mean"] if "gex" in r.get("metrics", {}) and "error" not in r["metrics"]["gex"] else None,
        "gex_MAE_mean": r["metrics"]["gex"]["MAE_mean"] if "gex" in r.get("metrics", {}) and "error" not in r["metrics"]["gex"] else None,
        "gex_Pearson_r": r["metrics"]["gex"]["Pearson_r_mean"] if "gex" in r.get("metrics", {}) and "error" not in r["metrics"]["gex"] else None,
        "runtime_sec": r.get("total_time", 0),
    }
    for r in all_results
])
summary_df.to_csv(METRICS_DIR / "benchmark_summary.csv", index=False)

print(f"Report saved to: {METRICS_DIR / 'final_benchmark_report.json'}")
print(f"Summary saved to: {METRICS_DIR / 'benchmark_summary.csv'}")

print("\n" + "="*60)
print("XENIUM BENCHMARKING COMPLETE")
print("="*60)
EOF

echo ""
echo "============================================"
echo "Integration complete at $(date)"
echo "============================================"
