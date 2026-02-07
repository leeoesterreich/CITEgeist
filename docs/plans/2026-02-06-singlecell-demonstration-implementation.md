# Single-Cell Demonstration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Run CITEgeist Modules 1→2→4 on Xenium RCC single-cell data (465K cells) to demonstrate method applicability at single-cell resolution, validate against known RCC biology, and discover novel spatial programs.

**Architecture:** Four new scripts orchestrate the pipeline: (1) `run_singlecell_module12.py` for Modules 1-2, (2) `run_singlecell_module4.py` for Module 4, (3) `evaluate_singlecell.py` for validation/GSEA, (4) `generate_figures.py` for publication figures. SLURM scripts run full dataset and quadrant jobs in parallel.

**Tech Stack:** Python 3.10, scanpy, numpy, scipy, sklearn (NMF), gseapy (GSEA), matplotlib/seaborn (figures). Uses existing CITEgeist model modules.

---

## Task 1: Create Output Directory Structure

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/output_singlecell_demonstration/.gitkeep`

**Step 1: Create directory structure**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/singlecell-demonstration

mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/output_singlecell_demonstration/{full,quadrants/{Q0,Q1,Q2,Q3},evaluation,figures}
touch Benchmarking/xenium_benchmarking/CITEgeist/output_singlecell_demonstration/.gitkeep
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/output_singlecell_demonstration/
git commit -m "chore: create output directory structure for single-cell demonstration"
```

---

## Task 2: Create Module 1-2 Runner Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_module12.py`

**Step 1: Write the script**

```python
#!/usr/bin/env python
"""
Run Modules 1 and 2 on Xenium single-cell data.

Executes marker interest detection (Module 1) and profile discovery (Module 2)
at single-cell resolution. Supports full dataset or quadrant subsets.

Usage:
    python run_singlecell_module12.py --mode full
    python run_singlecell_module12.py --mode quadrant --quadrant-id 0
"""

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles_continuous,
    select_profiles,
)
from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Output directory
OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell_demonstration"


def get_quadrant_bounds(coords: np.ndarray, quadrant_id: int) -> Tuple[float, float, float, float]:
    """
    Get bounds for a spatial quadrant.

    Quadrants:
        Q0: bottom-left, Q1: bottom-right, Q2: top-left, Q3: top-right
    """
    x_mid = (coords[:, 0].min() + coords[:, 0].max()) / 2
    y_mid = (coords[:, 1].min() + coords[:, 1].max()) / 2
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

    bounds = {
        0: (x_min, x_mid, y_min, y_mid),  # bottom-left
        1: (x_mid, x_max, y_min, y_mid),  # bottom-right
        2: (x_min, x_mid, y_mid, y_max),  # top-left
        3: (x_mid, x_max, y_mid, y_max),  # top-right
    }
    return bounds[quadrant_id]


def load_data(
    mode: str,
    quadrant_id: Optional[int] = None,
) -> Tuple[sc.AnnData, sc.AnnData, str]:
    """
    Load Xenium data for full or quadrant mode.

    Returns:
        Tuple of (adata_gex, adata_protein, output_subdir)
    """
    if mode == "full":
        logger.info("Loading FULL dataset (all cells)")
        adata_gex, adata_protein = load_xenium_singlecell()
        output_subdir = "full"
    else:
        logger.info(f"Loading QUADRANT {quadrant_id}")
        # First load full to get bounds
        adata_full, _ = load_xenium_singlecell(max_cells=1000)  # Quick load for bounds
        full_coords = adata_full.obsm["spatial"]
        bounds = get_quadrant_bounds(full_coords, quadrant_id)

        # Now load with bounds
        adata_gex, adata_protein = load_xenium_singlecell(region_bounds=bounds)
        output_subdir = f"quadrants/Q{quadrant_id}"

    logger.info(f"Loaded {adata_gex.shape[0]:,} cells × {adata_gex.shape[1]} genes")
    logger.info(f"Proteins: {list(adata_protein.var_names)}")

    return adata_gex, adata_protein, output_subdir


def run_module1(
    X_protein: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    output_dir: Path,
) -> Dict:
    """Run Module 1: Marker Interest Detection."""
    logger.info("=" * 70)
    logger.info("MODULE 1: Marker Interest Detection")
    logger.info("=" * 70)

    # At single-cell resolution, use more neighbors for Moran's I
    result = identify_interesting_markers(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        kurtosis_threshold=2.0,
        morans_threshold=0.1,
        morans_k=15,  # More neighbors at single-cell resolution
        morans_n_perm=199,
        verbose=True,
    )

    # Build output
    output = {
        "n_markers_total": len(marker_names),
        "n_interesting": len(result.interesting_markers),
        "n_boring": len(result.boring_markers),
        "interesting_markers": result.interesting_markers,
        "boring_markers": result.boring_markers,
        "kurtosis_threshold": float(result.kurtosis_threshold),
        "morans_threshold": float(result.morans_threshold),
        "marker_details": [
            {
                "marker": m.name,
                "interest_score": float(m.interest_score),
                "kurtosis": float(m.kurtosis),
                "gmm_snr": float(m.gmm_snr),
                "morans_i": float(m.morans_i),
                "morans_i_pvalue": float(m.morans_i_pvalue),
                "passed_kurtosis": bool(m.passed_kurtosis),
                "passed_gmm": bool(m.passed_gmm),
                "passed_morans": bool(m.passed_morans),
            }
            for m in result.markers
        ],
    }

    # Save outputs
    with open(output_dir / "module1_marker_interest.json", "w") as f:
        json.dump(output, f, indent=2)

    df = result.to_dataframe()
    df.to_csv(output_dir / "module1_marker_interest.csv", index=False)

    logger.info(f"Interesting markers: {len(result.interesting_markers)}/{len(marker_names)}")
    logger.info(f"Markers: {result.interesting_markers}")

    return output, result


def run_module2(
    X_protein: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    interesting_markers: List[str],
    signal_masks: np.ndarray,
    signal_mask_names: List[str],
    output_dir: Path,
) -> Dict:
    """Run Module 2: Profile Discovery."""
    logger.info("=" * 70)
    logger.info("MODULE 2: Profile Discovery")
    logger.info("=" * 70)

    # Filter to interesting markers
    marker_indices = [marker_names.index(m) for m in interesting_markers]
    X_filtered = X_protein[:, marker_indices]

    # Module 2a: Colocalization analysis
    logger.info("Module 2a: Colocalization Analysis")
    coloc_result = analyze_marker_colocalization(
        X=X_filtered,
        coords=coords,
        marker_names=interesting_markers,
        signal_masks=signal_masks,
        signal_mask_marker_names=signal_mask_names,
        neighbor_k=15,  # More neighbors at single-cell resolution
        n_perm=199,
        n_jobs=8,
    )

    # Save colocalization results
    coloc_df = coloc_result.to_dataframe()
    coloc_df.to_csv(output_dir / "module2a_colocalization.csv", index=False)

    # Module 2b: Profile discovery
    logger.info("Module 2b: Profile Discovery")
    profile_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        X=X_filtered,
        coords=coords,
        marker_names=interesting_markers,
        signal_masks=signal_masks,
        signal_mask_marker_names=signal_mask_names,
        min_cluster_size=2,
        bivariate_pvalue_threshold=0.05,
    )

    # Save raw profiles
    raw_profiles = {
        name: {"markers": list(markers)}
        for name, markers in profile_result.profiles.items()
    }
    with open(output_dir / "module2b_profiles_raw.json", "w") as f:
        json.dump(raw_profiles, f, indent=2)

    # Module 2c: Profile selection
    logger.info("Module 2c: Profile Selection")
    selected_profiles = select_profiles(
        profile_result=profile_result,
        X=X_filtered,
        marker_names=interesting_markers,
        max_profiles=10,
        min_reconstruction_improvement=0.01,
    )

    # Save selected profiles
    selected_output = {
        name: {"markers": list(markers)}
        for name, markers in selected_profiles.items()
    }
    with open(output_dir / "module2c_profiles_selected.json", "w") as f:
        json.dump(selected_output, f, indent=2)

    logger.info(f"Discovered {len(selected_profiles)} profiles:")
    for name, markers in selected_profiles.items():
        logger.info(f"  {name}: {list(markers)}")

    return {
        "n_profiles_raw": len(profile_result.profiles),
        "n_profiles_selected": len(selected_profiles),
        "profiles": selected_output,
    }


def main():
    parser = argparse.ArgumentParser(description="Run Modules 1-2 on single-cell data")
    parser.add_argument(
        "--mode",
        choices=["full", "quadrant"],
        required=True,
        help="Run on full dataset or single quadrant",
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

    # Load data
    adata_gex, adata_protein, output_subdir = load_data(args.mode, args.quadrant_id)

    # Setup output directory
    output_dir = OUTPUT_BASE / output_subdir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract arrays
    X_protein = np.asarray(adata_protein.X.todense() if hasattr(adata_protein.X, "todense") else adata_protein.X)
    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    # Save data summary
    summary = {
        "mode": args.mode,
        "quadrant_id": args.quadrant_id,
        "n_cells": int(adata_gex.shape[0]),
        "n_genes": int(adata_gex.shape[1]),
        "n_proteins": int(adata_protein.shape[1]),
        "protein_names": marker_names,
        "spatial_extent": {
            "x_min": float(coords[:, 0].min()),
            "x_max": float(coords[:, 0].max()),
            "y_min": float(coords[:, 1].min()),
            "y_max": float(coords[:, 1].max()),
        },
    }
    with open(output_dir / "data_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Run Module 1
    module1_output, module1_result = run_module1(
        X_protein=X_protein,
        coords=coords,
        marker_names=marker_names,
        output_dir=output_dir,
    )

    # Run Module 2
    module2_output = run_module2(
        X_protein=X_protein,
        coords=coords,
        marker_names=marker_names,
        interesting_markers=module1_result.interesting_markers,
        signal_masks=module1_result.signal_masks,
        signal_mask_names=module1_result.signal_mask_marker_names,
        output_dir=output_dir,
    )

    # Save combined summary
    combined_summary = {
        "data": summary,
        "module1": module1_output,
        "module2": module2_output,
    }
    with open(output_dir / "module12_summary.json", "w") as f:
        json.dump(combined_summary, f, indent=2)

    logger.info("=" * 70)
    logger.info("MODULE 1-2 COMPLETE")
    logger.info(f"Output saved to: {output_dir}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_module12.py
git commit -m "feat: add Module 1-2 runner for single-cell demonstration"
```

---

## Task 3: Create Module 4 Runner Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_module4.py`

**Step 1: Write the script**

```python
#!/usr/bin/env python
"""
Run Module 4 on Xenium single-cell data using discovered profiles.

At single-cell resolution, assigns cells to profiles (no deconvolution needed),
then runs NMF to discover spatial transcriptomic programs within each cell type.

Usage:
    python run_singlecell_module4.py --mode full
    python run_singlecell_module4.py --mode quadrant --quadrant-id 0
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Output directory
OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell_demonstration"


def get_quadrant_bounds(coords: np.ndarray, quadrant_id: int) -> Tuple[float, float, float, float]:
    """Get bounds for a spatial quadrant."""
    x_mid = (coords[:, 0].min() + coords[:, 0].max()) / 2
    y_mid = (coords[:, 1].min() + coords[:, 1].max()) / 2
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

    bounds = {
        0: (x_min, x_mid, y_min, y_mid),
        1: (x_mid, x_max, y_min, y_mid),
        2: (x_min, x_mid, y_mid, y_max),
        3: (x_mid, x_max, y_mid, y_max),
    }
    return bounds[quadrant_id]


def assign_cells_to_profiles(
    X_protein: np.ndarray,
    marker_names: List[str],
    profiles: Dict[str, Dict],
    method: str = "soft",
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Assign cells to discovered profiles based on protein expression.

    Args:
        X_protein: Protein expression matrix (n_cells × n_proteins)
        marker_names: Protein names
        profiles: Dict of profile_name -> {"markers": [list of markers]}
        method: "soft" for probability assignment, "hard" for argmax

    Returns:
        assignments: (n_cells,) array of profile indices (hard) or
                     (n_cells, n_profiles) probability matrix (soft)
        profile_names: List of profile names in order
    """
    profile_names = list(profiles.keys())
    n_cells = X_protein.shape[0]
    n_profiles = len(profile_names)

    # Compute profile scores for each cell
    scores = np.zeros((n_cells, n_profiles))

    for i, (name, profile) in enumerate(profiles.items()):
        markers = profile["markers"]
        marker_indices = [marker_names.index(m) for m in markers if m in marker_names]

        if marker_indices:
            # Average normalized expression of profile markers
            profile_expr = X_protein[:, marker_indices]
            # Normalize per marker, then average
            profile_expr_norm = profile_expr / (profile_expr.max(axis=0, keepdims=True) + 1e-10)
            scores[:, i] = profile_expr_norm.mean(axis=1)

    if method == "soft":
        # Softmax to get probabilities
        scores_exp = np.exp(scores - scores.max(axis=1, keepdims=True))
        assignments = scores_exp / scores_exp.sum(axis=1, keepdims=True)
    else:
        # Hard assignment
        assignments = scores.argmax(axis=1)

    return assignments, scores, profile_names


def compute_morans_i(
    values: np.ndarray,
    coords: np.ndarray,
    k: int = 15,
) -> float:
    """Compute Moran's I for spatial autocorrelation."""
    n = len(values)
    if n < k + 1:
        return 0.0

    # Build KNN graph
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1)

    # Compute Moran's I
    z = values - values.mean()
    z_sq_sum = (z ** 2).sum()

    if z_sq_sum == 0:
        return 0.0

    # Sum of cross-products for neighbors
    cross_sum = 0.0
    for i in range(n):
        for j in indices[i, 1:]:  # Skip self
            cross_sum += z[i] * z[j]

    # Total weight (k neighbors per cell)
    W = n * k

    morans_i = (n / W) * (cross_sum / z_sq_sum)
    return float(morans_i)


def run_nmf_for_celltype(
    adata_gex: sc.AnnData,
    cell_mask: np.ndarray,
    n_programs: int = 5,
    min_cells: int = 500,
) -> Optional[Dict]:
    """
    Run NMF on cells of a specific type to discover programs.

    Returns:
        Dict with W (gene loadings), H (cell scores), top genes per program, Moran's I
    """
    n_cells = cell_mask.sum()
    if n_cells < min_cells:
        logger.warning(f"Only {n_cells} cells, skipping (min: {min_cells})")
        return None

    # Extract expression for this cell type
    X = adata_gex[cell_mask, :].X
    if hasattr(X, "todense"):
        X = np.asarray(X.todense())

    # Filter low-variance genes
    gene_var = X.var(axis=0)
    top_genes_mask = gene_var > np.percentile(gene_var, 50)
    X_filtered = X[:, top_genes_mask]
    gene_names_filtered = adata_gex.var_names[top_genes_mask].tolist()

    # Normalize
    X_norm = X_filtered / (X_filtered.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    X_norm = np.log1p(X_norm)

    # Make non-negative
    X_norm = X_norm - X_norm.min()

    # Run NMF
    logger.info(f"Running NMF with K={n_programs} on {n_cells} cells × {X_norm.shape[1]} genes")
    model = NMF(n_components=n_programs, init="nndsvd", random_state=42, max_iter=500)
    H = model.fit_transform(X_norm)  # (n_cells, n_programs)
    W = model.components_.T  # (n_genes, n_programs)

    # Get coordinates for this cell type
    coords = adata_gex[cell_mask, :].obsm["spatial"]

    # Compute top genes and Moran's I per program
    programs = []
    for prog_idx in range(n_programs):
        # Top genes
        gene_loadings = W[:, prog_idx]
        top_indices = np.argsort(gene_loadings)[::-1][:20]
        top_genes = [gene_names_filtered[i] for i in top_indices]
        top_loadings = gene_loadings[top_indices].tolist()

        # Moran's I for spatial coherence
        cell_scores = H[:, prog_idx]
        morans_i = compute_morans_i(cell_scores, coords, k=15)

        programs.append({
            "program_id": prog_idx,
            "top_genes": top_genes,
            "top_loadings": top_loadings,
            "morans_i": morans_i,
            "variance_explained": float((cell_scores ** 2).sum() / (H ** 2).sum()),
        })

    return {
        "n_cells": int(n_cells),
        "n_genes": int(X_norm.shape[1]),
        "n_programs": n_programs,
        "programs": programs,
        "W_shape": list(W.shape),
        "H_shape": list(H.shape),
        "reconstruction_error": float(model.reconstruction_err_),
    }


def main():
    parser = argparse.ArgumentParser(description="Run Module 4 on single-cell data")
    parser.add_argument(
        "--mode",
        choices=["full", "quadrant"],
        required=True,
        help="Run on full dataset or single quadrant",
    )
    parser.add_argument(
        "--quadrant-id",
        type=int,
        choices=[0, 1, 2, 3],
        help="Quadrant ID (required if mode=quadrant)",
    )
    parser.add_argument(
        "--n-programs",
        type=int,
        default=5,
        help="Number of NMF programs per cell type (default: 5)",
    )
    args = parser.parse_args()

    if args.mode == "quadrant" and args.quadrant_id is None:
        parser.error("--quadrant-id required when mode=quadrant")

    # Determine paths
    if args.mode == "full":
        output_subdir = "full"
    else:
        output_subdir = f"quadrants/Q{args.quadrant_id}"

    output_dir = OUTPUT_BASE / output_subdir

    # Load Module 2 profiles
    profiles_path = output_dir / "module2c_profiles_selected.json"
    if not profiles_path.exists():
        logger.error(f"Profiles not found at {profiles_path}. Run Module 1-2 first.")
        sys.exit(1)

    with open(profiles_path) as f:
        profiles = json.load(f)

    logger.info(f"Loaded {len(profiles)} profiles from Module 2")

    # Load data
    if args.mode == "full":
        adata_gex, adata_protein = load_xenium_singlecell()
    else:
        adata_full, _ = load_xenium_singlecell(max_cells=1000)
        bounds = get_quadrant_bounds(adata_full.obsm["spatial"], args.quadrant_id)
        adata_gex, adata_protein = load_xenium_singlecell(region_bounds=bounds)

    logger.info(f"Loaded {adata_gex.shape[0]:,} cells")

    # Extract arrays
    X_protein = np.asarray(adata_protein.X.todense() if hasattr(adata_protein.X, "todense") else adata_protein.X)
    marker_names = list(adata_protein.var_names)

    # Assign cells to profiles
    logger.info("Assigning cells to profiles...")
    assignments, scores, profile_names = assign_cells_to_profiles(
        X_protein=X_protein,
        marker_names=marker_names,
        profiles=profiles,
        method="hard",
    )

    # Save assignments
    assignments_df = pd.DataFrame({
        "cell_id": adata_gex.obs_names,
        "profile": [profile_names[i] for i in assignments],
        "profile_score": scores.max(axis=1),
    })
    assignments_df.to_csv(output_dir / "cell_assignments.csv", index=False)

    # Count cells per profile
    for i, name in enumerate(profile_names):
        n = (assignments == i).sum()
        logger.info(f"  {name}: {n:,} cells ({100*n/len(assignments):.1f}%)")

    # Create output directory for programs
    programs_dir = output_dir / "module4_programs"
    programs_dir.mkdir(exist_ok=True)

    # Run NMF for each cell type
    all_results = {}
    for i, profile_name in enumerate(profile_names):
        logger.info("=" * 70)
        logger.info(f"Processing: {profile_name}")
        logger.info("=" * 70)

        cell_mask = assignments == i
        result = run_nmf_for_celltype(
            adata_gex=adata_gex,
            cell_mask=cell_mask,
            n_programs=args.n_programs,
            min_cells=500,
        )

        if result is not None:
            all_results[profile_name] = result

            # Save per-cell-type results
            with open(programs_dir / f"{profile_name.replace(' ', '_')}_programs.json", "w") as f:
                json.dump(result, f, indent=2)

    # Save summary
    summary = {
        "mode": args.mode,
        "quadrant_id": args.quadrant_id,
        "n_cells_total": int(adata_gex.shape[0]),
        "n_profiles": len(profile_names),
        "profile_names": profile_names,
        "n_programs_per_type": args.n_programs,
        "results": {
            name: {
                "n_cells": r["n_cells"],
                "n_programs": r["n_programs"],
                "mean_morans_i": np.mean([p["morans_i"] for p in r["programs"]]),
            }
            for name, r in all_results.items()
        },
    }
    with open(output_dir / "module4_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("=" * 70)
    logger.info("MODULE 4 COMPLETE")
    logger.info(f"Output saved to: {output_dir}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_module4.py
git commit -m "feat: add Module 4 runner for single-cell demonstration"
```

---

## Task 4: Create Evaluation Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_singlecell.py`

**Step 1: Write the script**

```python
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
            all_markers = required | optional

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
        lines.append(f"- **{match['discovered']}** → {match['expected']} (score: {match['score']:.2f})")
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
            coherent = "✓" if prog["spatially_coherent"] else "✗"
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
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_singlecell.py
git commit -m "feat: add evaluation script for single-cell demonstration"
```

---

## Task 5: Create Figure Generation Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/generate_figures.py`

**Step 1: Write the script**

```python
#!/usr/bin/env python
"""
Generate publication figures for single-cell demonstration.

Creates:
1. Profile discovery heatmap
2. Spatial cell type map
3. Program discovery summary
4. Spatial program visualization
5. Validation summary

Usage:
    python generate_figures.py --mode full
    python generate_figures.py --mode quadrant --quadrant-id 0
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell_demonstration"

# Set style
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 10


def get_quadrant_bounds(coords: np.ndarray, quadrant_id: int):
    """Get bounds for a spatial quadrant."""
    x_mid = (coords[:, 0].min() + coords[:, 0].max()) / 2
    y_mid = (coords[:, 1].min() + coords[:, 1].max()) / 2
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

    bounds = {
        0: (x_min, x_mid, y_min, y_mid),
        1: (x_mid, x_max, y_min, y_mid),
        2: (x_min, x_mid, y_mid, y_max),
        3: (x_mid, x_max, y_mid, y_max),
    }
    return bounds[quadrant_id]


def plot_profile_heatmap(profiles: Dict, output_path: Path):
    """Plot heatmap of discovered profiles."""
    # Get all markers
    all_markers = set()
    for profile in profiles.values():
        all_markers.update(profile["markers"])
    all_markers = sorted(all_markers)

    # Build matrix
    profile_names = list(profiles.keys())
    matrix = np.zeros((len(profile_names), len(all_markers)))

    for i, (name, profile) in enumerate(profiles.items()):
        for marker in profile["markers"]:
            j = all_markers.index(marker)
            matrix[i, j] = 1

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(
        matrix,
        xticklabels=all_markers,
        yticklabels=profile_names,
        cmap="Blues",
        cbar_kws={"label": "Marker Present"},
        ax=ax,
    )
    ax.set_xlabel("Protein Markers")
    ax.set_ylabel("Discovered Profiles")
    ax.set_title("Discovered Cell Type Profiles")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_spatial_celltypes(
    assignments_df: pd.DataFrame,
    coords: np.ndarray,
    output_path: Path,
    max_points: int = 50000,
):
    """Plot spatial map of cell type assignments."""
    # Subsample if needed
    if len(assignments_df) > max_points:
        idx = np.random.choice(len(assignments_df), max_points, replace=False)
        assignments_df = assignments_df.iloc[idx]
        coords = coords[idx]

    # Get unique profiles
    profiles = assignments_df["profile"].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(profiles)))
    color_map = {p: c for p, c in zip(profiles, colors)}

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))

    for profile in profiles:
        mask = assignments_df["profile"] == profile
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=[color_map[profile]],
            label=profile,
            s=1,
            alpha=0.5,
        )

    ax.set_xlabel("X (microns)")
    ax.set_ylabel("Y (microns)")
    ax.set_title("Spatial Cell Type Distribution")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", markerscale=5)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_program_summary(programs_dir: Path, output_path: Path):
    """Plot summary of discovered programs per cell type."""
    # Load all program files
    data = []
    for program_file in sorted(programs_dir.glob("*_programs.json")):
        cell_type = program_file.stem.replace("_programs", "").replace("_", " ")

        with open(program_file) as f:
            prog_data = json.load(f)

        for prog in prog_data.get("programs", []):
            data.append({
                "Cell Type": cell_type,
                "Program": f"P{prog['program_id']}",
                "Moran's I": prog["morans_i"],
                "Top Genes": ", ".join(prog["top_genes"][:3]),
            })

    if not data:
        logger.warning("No program data found")
        return

    df = pd.DataFrame(data)

    # Plot Moran's I by cell type
    fig, ax = plt.subplots(figsize=(10, 6))

    cell_types = df["Cell Type"].unique()
    x = np.arange(len(cell_types))
    width = 0.15

    programs = sorted(df["Program"].unique())
    for i, prog in enumerate(programs):
        prog_data = df[df["Program"] == prog]
        values = [prog_data[prog_data["Cell Type"] == ct]["Moran's I"].values[0]
                  if ct in prog_data["Cell Type"].values else 0
                  for ct in cell_types]
        ax.bar(x + i * width, values, width, label=prog)

    ax.axhline(y=0.1, color="red", linestyle="--", alpha=0.5, label="Coherence threshold")
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Moran's I (Spatial Coherence)")
    ax.set_title("Spatial Coherence of Discovered Programs")
    ax.set_xticks(x + width * (len(programs) - 1) / 2)
    ax.set_xticklabels(cell_types, rotation=45, ha="right")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_validation_summary(eval_path: Path, output_path: Path):
    """Plot validation results."""
    with open(eval_path) as f:
        eval_data = json.load(f)

    metrics = eval_data["metrics"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Metrics bar chart
    ax1 = axes[0]
    metric_names = ["Precision", "Recall", "F1"]
    values = [metrics["precision"], metrics["recall"], metrics["f1"]]
    colors = ["#2ecc71", "#3498db", "#9b59b6"]

    bars = ax1.bar(metric_names, values, color=colors)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel("Score")
    ax1.set_title("Profile Recovery Metrics")

    for bar, val in zip(bars, values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f"{val:.2f}", ha="center", va="bottom")

    # Matched vs unmatched
    ax2 = axes[1]
    matched = metrics["n_matched"]
    unmatched_exp = len(eval_data["unmatched_expected"])
    unmatched_disc = len(eval_data["unmatched_discovered"])

    categories = ["Matched", "Missing\n(Expected)", "Novel\n(Discovered)"]
    values = [matched, unmatched_exp, unmatched_disc]
    colors = ["#2ecc71", "#e74c3c", "#f39c12"]

    ax2.bar(categories, values, color=colors)
    ax2.set_ylabel("Count")
    ax2.set_title("Profile Matching Summary")

    for i, val in enumerate(values):
        ax2.text(i, val + 0.1, str(val), ha="center", va="bottom")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate figures for single-cell results")
    parser.add_argument(
        "--mode",
        choices=["full", "quadrant"],
        required=True,
        help="Generate figures for full or quadrant results",
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

    # Determine paths
    if args.mode == "full":
        data_dir = OUTPUT_BASE / "full"
    else:
        data_dir = OUTPUT_BASE / "quadrants" / f"Q{args.quadrant_id}"

    figures_dir = OUTPUT_BASE / "figures"
    figures_dir.mkdir(exist_ok=True)

    prefix = "full" if args.mode == "full" else f"Q{args.quadrant_id}"

    # Load profiles
    profiles_path = data_dir / "module2c_profiles_selected.json"
    if profiles_path.exists():
        with open(profiles_path) as f:
            profiles = json.load(f)
        plot_profile_heatmap(profiles, figures_dir / f"fig_profiles_{prefix}.png")

    # Load cell assignments and coordinates
    assignments_path = data_dir / "cell_assignments.csv"
    if assignments_path.exists():
        assignments_df = pd.read_csv(assignments_path)

        # Load coordinates
        if args.mode == "full":
            _, adata_protein = load_xenium_singlecell(max_cells=50000)
        else:
            adata_full, _ = load_xenium_singlecell(max_cells=1000)
            bounds = get_quadrant_bounds(adata_full.obsm["spatial"], args.quadrant_id)
            _, adata_protein = load_xenium_singlecell(region_bounds=bounds, max_cells=50000)

        coords = adata_protein.obsm["spatial"]

        # Align if subsampled
        if len(assignments_df) != len(coords):
            # Use subset
            n = min(len(assignments_df), len(coords))
            assignments_df = assignments_df.iloc[:n]
            coords = coords[:n]

        plot_spatial_celltypes(assignments_df, coords, figures_dir / f"fig_spatial_celltypes_{prefix}.png")

    # Plot program summary
    programs_dir = data_dir / "module4_programs"
    if programs_dir.exists():
        plot_program_summary(programs_dir, figures_dir / f"fig_programs_{prefix}.png")

    # Plot validation summary
    eval_dir = data_dir / "evaluation" if (data_dir / "evaluation").exists() else OUTPUT_BASE / "evaluation"
    eval_path = eval_dir / f"profile_validation_{data_dir.name}.json"
    if eval_path.exists():
        plot_validation_summary(eval_path, figures_dir / f"fig_validation_{prefix}.png")

    logger.info(f"Figures saved to: {figures_dir}")


if __name__ == "__main__":
    main()
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/generate_figures.py
git commit -m "feat: add figure generation script for single-cell demonstration"
```

---

## Task 6: Create SLURM Scripts

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_singlecell_full.sh`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_singlecell_quadrants.sh`

**Step 1: Write full dataset script**

```bash
#!/bin/bash
#SBATCH --job-name=singlecell_full
#SBATCH --output=slurm_log/singlecell_full_%j.out
#SBATCH --error=slurm_log/singlecell_full_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Single-cell demonstration: Full dataset (465K cells)
# Runs Modules 1-2 then Module 4 sequentially

set -e

echo "=========================================="
echo "Single-Cell Demonstration: Full Dataset"
echo "Started: $(date)"
echo "=========================================="

# Load environment
module load python/ondemand-jupyter-python3.11

# Change to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/src

# Run Module 1-2
echo ""
echo "Running Modules 1-2..."
python run_singlecell_module12.py --mode full

# Run Module 4
echo ""
echo "Running Module 4..."
python run_singlecell_module4.py --mode full

# Run evaluation
echo ""
echo "Running evaluation..."
python evaluate_singlecell.py --mode full

# Generate figures
echo ""
echo "Generating figures..."
python generate_figures.py --mode full

echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "=========================================="
```

**Step 2: Write quadrant array script**

```bash
#!/bin/bash
#SBATCH --job-name=singlecell_quad
#SBATCH --output=slurm_log/singlecell_quad_%A_%a.out
#SBATCH --error=slurm_log/singlecell_quad_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --array=0-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Single-cell demonstration: Quadrants (array job)
# Each array task processes one quadrant (~116K cells)

set -e

QUADRANT_ID=$SLURM_ARRAY_TASK_ID

echo "=========================================="
echo "Single-Cell Demonstration: Quadrant ${QUADRANT_ID}"
echo "Started: $(date)"
echo "=========================================="

# Load environment
module load python/ondemand-jupyter-python3.11

# Change to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/src

# Run Module 1-2
echo ""
echo "Running Modules 1-2..."
python run_singlecell_module12.py --mode quadrant --quadrant-id ${QUADRANT_ID}

# Run Module 4
echo ""
echo "Running Module 4..."
python run_singlecell_module4.py --mode quadrant --quadrant-id ${QUADRANT_ID}

# Run evaluation
echo ""
echo "Running evaluation..."
python evaluate_singlecell.py --mode quadrant --quadrant-id ${QUADRANT_ID}

# Generate figures
echo ""
echo "Generating figures..."
python generate_figures.py --mode quadrant --quadrant-id ${QUADRANT_ID}

echo ""
echo "=========================================="
echo "Quadrant ${QUADRANT_ID} Completed: $(date)"
echo "=========================================="
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_singlecell_full.sh
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_singlecell_quadrants.sh
git commit -m "feat: add SLURM scripts for single-cell demonstration"
```

---

## Task 7: Test Quadrant Data Loading

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/load_xenium_singlecell.py`

**Step 1: Add quadrant support to loader**

Add after the existing `get_region_bounds` function (around line 170):

```python
def get_quadrant_bounds(
    data_dir: str = XENIUM_DATA_DIR,
) -> List[Tuple[float, float, float, float]]:
    """
    Get quadrant boundaries for the Xenium dataset.

    Quadrants are defined by X/Y midpoint:
        Q0: bottom-left, Q1: bottom-right, Q2: top-left, Q3: top-right

    Returns:
        List of (x_min, x_max, y_min, y_max) tuples for each quadrant.
    """
    try:
        from load_xenium import load_xenium_data
    except ImportError:
        from Benchmarking.xenium_pseudovisium.src.load_xenium import load_xenium_data

    adata = load_xenium_data(data_dir, min_counts=0)
    coords = adata.obsm["spatial"]

    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    x_mid = (x_min + x_max) / 2
    y_mid = (y_min + y_max) / 2

    return [
        (x_min, x_mid, y_min, y_mid),  # Q0: bottom-left
        (x_mid, x_max, y_min, y_mid),  # Q1: bottom-right
        (x_min, x_mid, y_mid, y_max),  # Q2: top-left
        (x_mid, x_max, y_mid, y_max),  # Q3: top-right
    ]
```

**Step 2: Test the loader**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/singlecell-demonstration
module load python/ondemand-jupyter-python3.11
python -c "
from Benchmarking.xenium_benchmarking.CITEgeist.src.load_xenium_singlecell import get_quadrant_bounds
bounds = get_quadrant_bounds()
for i, b in enumerate(bounds):
    print(f'Q{i}: x=[{b[0]:.0f}, {b[1]:.0f}], y=[{b[2]:.0f}, {b[3]:.0f}]')
"
```

Expected output:
```
Q0: x=[2, 5737], y=[3, 2904]
Q1: x=[5737, 11473], y=[3, 2904]
Q2: x=[2, 5737], y=[2904, 5805]
Q3: x=[5737, 11473], y=[2904, 5805]
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/load_xenium_singlecell.py
git commit -m "feat: add quadrant bounds function to single-cell loader"
```

---

## Task 8: Create slurm_log Directory

**Step 1: Create directory**

```bash
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log
```

**Step 2: Verify scripts are executable**

```bash
chmod +x Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_singlecell_*.sh
```

---

## Task 9: Run Quadrant Jobs First (Faster Validation)

**Step 1: Submit quadrant jobs**

```bash
cd Benchmarking/xenium_benchmarking/CITEgeist/slurm
sbatch run_singlecell_quadrants.sh
```

Expected: Job array with 4 tasks (0-3) submitted.

**Step 2: Monitor progress**

```bash
squeue -u $USER
```

---

## Task 10: Run Full Dataset Job

**Step 1: Submit full job (after quadrants validate)**

```bash
cd Benchmarking/xenium_benchmarking/CITEgeist/slurm
sbatch run_singlecell_full.sh
```

---

## Task 11: Aggregate Results

**After all jobs complete:**

**Step 1: Run aggregated evaluation**

```bash
cd Benchmarking/xenium_benchmarking/CITEgeist/src
python evaluate_singlecell.py --mode all
```

**Step 2: Review outputs**

```bash
ls -la ../output_singlecell_demonstration/
ls -la ../output_singlecell_demonstration/figures/
cat ../output_singlecell_demonstration/evaluation/validation_summary.json
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Create output directories | output_singlecell_demonstration/ |
| 2 | Module 1-2 runner | run_singlecell_module12.py |
| 3 | Module 4 runner | run_singlecell_module4.py |
| 4 | Evaluation script | evaluate_singlecell.py |
| 5 | Figure generation | generate_figures.py |
| 6 | SLURM scripts | run_singlecell_full.sh, run_singlecell_quadrants.sh |
| 7 | Add quadrant support | load_xenium_singlecell.py |
| 8 | Create slurm_log | slurm/slurm_log/ |
| 9 | Run quadrant jobs | sbatch array job |
| 10 | Run full dataset | sbatch full job |
| 11 | Aggregate results | evaluate_singlecell.py --mode all |
