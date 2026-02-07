"""
Quick check: which of the 27 proteins best discriminate Fibroblasts from Epithelial?
"""
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from load_xenium_singlecell import load_xenium_singlecell

GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_type_assignments.csv"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Load region 2, 5k cells
adata_gex, adata_protein = load_xenium_singlecell(region_id=2, max_cells=5000, seed=42)

gt_df = pd.read_csv(GT_PATH, index_col=0)
common = sorted(set(adata_protein.obs_names) & set(gt_df.index))
gt = gt_df.loc[common, "cell_type"]

X = adata_protein.X
if hasattr(X, "toarray"):
    X = X.toarray()
X = np.asarray(X, dtype=np.float64)

protein_names = list(adata_protein.var_names)
cell_names = list(adata_protein.obs_names)
cell_to_idx = {c: i for i, c in enumerate(cell_names)}

logger.info(f"Full panel ({len(protein_names)} proteins): {protein_names}")

# Get indices for each GT type
type_indices = {}
for gt_type in sorted(gt.unique()):
    cells = [c for c in common if gt[c] == gt_type]
    indices = [cell_to_idx[c] for c in cells]
    type_indices[gt_type] = indices

# Print mean expression per type per protein
logger.info(f"\n{'Protein':<15}" + "".join(f"{t:>14}" for t in sorted(type_indices.keys())))
logger.info(f"{'N cells':<15}" + "".join(f"{len(v):>14}" for v in [type_indices[t] for t in sorted(type_indices.keys())]))
logger.info("-" * (15 + 14 * len(type_indices)))

for p_idx, p_name in enumerate(protein_names):
    row = f"{p_name:<15}"
    for gt_type in sorted(type_indices.keys()):
        idx = type_indices[gt_type]
        mean_val = X[idx, p_idx].mean() if idx else 0
        row += f"{mean_val:>14.4f}"
    logger.info(row)

# Specifically: Fibroblast vs Epithelial discrimination
logger.info("\n" + "=" * 80)
logger.info("FIBROBLAST vs EPITHELIAL discrimination")
logger.info("=" * 80)

fib_idx = type_indices.get("Fibroblasts", [])
epi_idx = type_indices.get("Epithelial", [])

if fib_idx and epi_idx:
    logger.info(f"N Fibroblasts: {len(fib_idx)}, N Epithelial: {len(epi_idx)}")
    logger.info(f"\n{'Protein':<15} {'Fib mean':>10} {'Epi mean':>10} {'Ratio F/E':>10} {'Ratio E/F':>10} {'|diff|':>10}")
    logger.info("-" * 65)

    results = []
    for p_idx, p_name in enumerate(protein_names):
        fib_mean = X[fib_idx, p_idx].mean()
        epi_mean = X[epi_idx, p_idx].mean()
        diff = abs(fib_mean - epi_mean)
        ratio_fe = fib_mean / (epi_mean + 1e-10)
        ratio_ef = epi_mean / (fib_mean + 1e-10)
        results.append((p_name, fib_mean, epi_mean, ratio_fe, ratio_ef, diff))

    # Sort by absolute difference
    results.sort(key=lambda x: -x[5])
    for p_name, fib_mean, epi_mean, ratio_fe, ratio_ef, diff in results:
        logger.info(f"{p_name:<15} {fib_mean:>10.4f} {epi_mean:>10.4f} {ratio_fe:>10.2f} {ratio_ef:>10.2f} {diff:>10.4f}")
