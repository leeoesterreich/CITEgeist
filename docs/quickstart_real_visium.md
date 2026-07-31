# Quickstart: Real Visium Data → Single-Cell Output

This guide walks through the complete CITEgeist pipeline on real Visium CITE-seq data, from SpaceRanger output to per-cell gene expression.

## What You Need

| Input | Description |
|-------|-------------|
| SpaceRanger output | Directory with `filtered_feature_bc_matrix.h5` and `spatial/` folder |
| H&E image | `spatial/tissue_fullres_image.tif` (or CytAssist/hires fallback) |
| GPU node | Required for Module 3 (cuOPT QP solver) and StarDist segmentation |

No external scRNA-seq reference is needed — CITEgeist builds cell-type profiles directly from the same-slide CITE-seq antibody panel.

## Pipeline Overview

```
SpaceRanger output
       │
       ▼
┌─────────────────────┐
│ Module 1-2:         │  Identify spatially interesting markers,
│ Profile Discovery   │  build cell-type protein profiles
└─────────┬───────────┘
          │ cell_profile_dict (JSON)
          ▼
┌─────────────────────┐
│ Module 3:           │  GPU-accelerated QP deconvolution →
│ Proportions         │  spot-level cell-type proportions
└─────────┬───────────┘
          │ proportions CSV
          ▼
┌─────────────────────┐
│ Module 3-post:      │  StarDist nuclei segmentation →
│ Cell Assignment     │  Hungarian assignment → SACE per-cell GEX
└─────────┬───────────┘
          │
          ▼
   Per-cell AnnData (.h5ad)
   • .X = gene expression per nucleus (count-scale allocations)
   • .obs["cell_type"] = assigned type
   • .obsm["spatial"] = (x, y) coordinates
```

---

## Step 0: Environment Setup

```bash
conda activate CITEgeist_env   # Case-sensitive! Not citegeist_env
```

All heavy computation should be submitted via SLURM (`sbatch`), not run on login nodes.

---

## Step 1: Profile Discovery (Modules 1–2)

This step identifies which antibody markers are spatially informative and groups them into cell-type profiles. Run once per cohort or tissue type — the output profiles are reusable across samples.

The discovery API is functional (not method-based on the model class):

```python
"""run_discovery.py — Module 1-2 profile discovery for one sample."""
import sys
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import scipy.sparse
import squidpy as sq

from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    select_profiles,
)
from CITEgeist.model.discovery.spatial_colocalization import discover_profiles_continuous

# --- Configuration ---
SAMPLE_NAME = "sample_P1_S1"
SPACERANGER_DIR = Path("/path/to/spaceranger/output/sample_P1_S1")
OUTPUT_DIR = Path("output/module12_discovery")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load from SpaceRanger (gex_only=False to include antibodies)
adata = sq.read.visium(
    str(SPACERANGER_DIR),
    counts_file="filtered_feature_bc_matrix.h5",
    load_images=True,
    gex_only=False,
)

# Extract antibody matrix (antibody features are typically the last N features)
# Real data has antibody names ending in "-1" (e.g., "CD68-1")
antibody_mask = adata.var_names.str.endswith("-1")
ab_adata = adata[:, antibody_mask]
X = ab_adata.X.toarray() if scipy.sparse.issparse(ab_adata.X) else np.array(ab_adata.X)
coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
marker_names = list(ab_adata.var_names)

# Module 1: Score markers for spatial interestingness
# (kurtosis, GMM SNR, Moran's I)
m1_result = identify_interesting_markers(
    X=X, coords=coords, marker_names=marker_names,
    morans_k=8, morans_n_perm=99,
)
print(f"Interesting markers: {m1_result.interesting_markers}")

# Module 2a: Pairwise spatial colocalization
coloc_result = analyze_marker_colocalization(
    X=X, coords=coords, marker_names=marker_names,
    markers_to_analyze=m1_result.interesting_markers,
    neighbor_k=6, n_permutations=999,
)

# Module 2b: Discover profiles from colocalization graph
discovery_result = discover_profiles_continuous(
    colocalization_result=coloc_result,
    top_k=3,
)

# Module 2c: Select profiles by spatial + protein variance explained
selection_result = select_profiles(
    X=X, coords=coords, marker_names=marker_names,
    profiles=discovery_result.profiles,
    _interesting_markers=m1_result.interesting_markers,
    _colocalization_result=coloc_result,
    min_spatial_explained=0.90,
    min_protein_explained=0.90,
)

# Save results
results = {
    "sample_name": SAMPLE_NAME,
    "interesting_markers": m1_result.interesting_markers,
    "discovered_profiles": [
        {"profile_id": i, "markers": list(p)}
        for i, p in enumerate(selection_result.selected_profiles)
    ],
}
with open(OUTPUT_DIR / f"{SAMPLE_NAME}_module12_discovery.json", "w") as f:
    json.dump(results, f, indent=2)
```

**Output:** A JSON with discovered marker groupings (raw profiles, not yet labeled):

```json
{
  "discovered_profiles": [
    {"profile_id": 0, "markers": ["CD68-1", "CD163-1"]},
    {"profile_id": 1, "markers": ["CD3E-1", "CD8A-1"]},
    {"profile_id": 2, "markers": ["PECAM1-1"]}
  ]
}
```

**Important:** The discovered profiles are unlabeled marker sets. You must assign biological cell-type names to create the `cell_profile_dict` format needed for Module 3:

```python
# Convert discovered profiles to the format expected by load_cell_profile_dict()
# This requires biological interpretation of each marker group
cell_profile_dict = {
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "Endothelial": {"Major": ["PECAM1-1"]},
    # ... assign a cell type name to each discovered marker group
}
```

> **Tip:** If you already have a curated profile (e.g., from prior CITE-seq knowledge), you can skip Module 1-2 entirely and provide the dict directly in Step 2. Marker names must include the `-1` antibody suffix for real data.

---

## Step 2: Spot-Level Proportions (Module 3)

This step estimates cell-type proportions per spot using GPU-accelerated quadratic programming.

```python
"""run_proportions.py — Module 3 QP on a single sample. Requires GPU node."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import squidpy as sq
from CITEgeist.model.citegeist_model import CitegeistModel

# --- Configuration ---
SAMPLE_NAME = "sample_P1_S1"
SPACERANGER_DIR = Path("/path/to/spaceranger/output/sample_P1_S1")
OUTPUT_DIR = Path("output/module3")
MIN_COUNTS = 100  # 100 for biopsy, 25 for surgical (low-depth) samples

# Cell-type profile dict (from Step 1 or curated)
CELL_PROFILE_DICT = {
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1", "CD3E-1"]},
    "Cancer_Luminal": {"Major": ["EPCAM-1"]},
    "Cancer_Basal": {"Major": ["KRT5-1", "SDC1-1", "EPCAM-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1", "HLA-DRA-1"]},
}

# Load Visium data
adata = sq.read.visium(
    str(SPACERANGER_DIR),
    counts_file="filtered_feature_bc_matrix.h5",
    load_images=True,
    gex_only=False,
)

# Initialize
model = CitegeistModel(
    sample_name=SAMPLE_NAME,
    adata=adata,
    output_folder=str(OUTPUT_DIR),
    simulation=False,
)

# Preprocessing
model.split_adata()
model.filter_gex(
    nonzero_percentage=0.01,
    mean_expression_threshold=1.1,
    min_counts=MIN_COUNTS,
)
model.copy_gex_to_protein_adata()
model.preprocess_gex(target_sum=10000)
model.preprocess_antibody()

# Load cell-type profiles
model.load_cell_profile_dict(CELL_PROFILE_DICT)

# Run cell-type proportions (cuOPT QP, requires GPU)
model.run_cell_proportion_model(
    method="qp",
    use_detection_gating=True,
    use_gex_detection=True,      # GEX-informed detection refinement
    gex_detection_k=10,
    validation_warn_only=True,
)

# Append proportions to adata and save
model.append_proportions_to_adata(key="finetuned")
final_adata = model.get_adata()
final_adata.write_h5ad(str(OUTPUT_DIR / f"{SAMPLE_NAME}_module3_results.h5ad"))
```

**Key outputs:**
- **Proportions CSV** — `{SAMPLE_NAME}_cell_prop_finetuned_results.csv` (spots x cell types), auto-saved
- **AnnData** with proportion columns appended to `adata.obs`

At this point you have spot-level deconvolution. Continue to Step 3 for single-cell resolution.

> **Note:** Gene expression deconvolution (SACE) happens *after* cell assignment in Step 3, not at the spot level. SACE allocates each spot's GEX counts to the individual nuclei assigned within that spot.

---

## Step 3: Single-Cell Assignment + GEX (Module 3-post)

This step segments individual nuclei from the H&E image, assigns each nucleus a cell type from the spot-level proportions (Hungarian matching), and allocates gene expression per cell via SACE.

The production script is `examples/scripts/run_single_cell_assignment.py`, which handles this in stages. Below is the conceptual flow:

```python
"""run_assignment.py — Nuclei segmentation + cell assignment + GEX. Requires GPU."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import squidpy as sq
from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.morphology.segmentation import (
    compute_spot_nuclei_counts_patchwise,
    estimate_pixel_size_um,
)

# --- Configuration ---
SAMPLE_NAME = "sample_P1_S1"
SPACERANGER_DIR = Path("/path/to/spaceranger/output/sample_P1_S1")
MODULE3_OUTPUT = Path("output/module3")
OUTPUT_DIR = Path("output/single_cell")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Same profile dict used in Step 2
CELL_PROFILE_DICT = {
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    # ... (same as Step 2)
}

# =========================================================================
# Stage 1: Nuclei segmentation (StarDist patchwise on H&E)
# =========================================================================
adata = sq.read.visium(
    str(SPACERANGER_DIR), counts_file="filtered_feature_bc_matrix.h5",
    load_images=True, gex_only=False,
)
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
fullres_img = np.array(
    Image.open(SPACERANGER_DIR / "spatial" / "tissue_fullres_image.tif").convert("RGB")
)
px_size = estimate_pixel_size_um(adata)

result = compute_spot_nuclei_counts_patchwise(
    fullres_image=fullres_img,
    spot_centers_fullres=np.asarray(adata.obsm["spatial"]),
    spot_names=adata.obs_names,
    pixel_size_um=px_size,
    modality="he",
)
# Save nucleus-spot mapping (columns: nucleus_id, spot_barcode, centroid_x, centroid_y)
nuc_map = result.nucleus_spot_map
nuc_map.to_csv(OUTPUT_DIR / "nucleus_spot_mapping.csv", index=False)
del fullres_img

# =========================================================================
# Stage 2: Hungarian cell assignment + SACE GEX
# =========================================================================
# Load pre-computed Module 3 proportions
prop_df = pd.read_csv(
    MODULE3_OUTPUT / f"{SAMPLE_NAME}_cell_prop_finetuned_results.csv", index_col=0
)

# Re-initialize model for GEX access
model = CitegeistModel(
    sample_name=SAMPLE_NAME, adata=adata,
    output_folder=str(OUTPUT_DIR), simulation=False,
)
model.split_adata()
model.filter_gex(min_counts=100)
model.preprocess_gex()
model.preprocess_antibody()
model.load_cell_profile_dict(CELL_PROFILE_DICT)

# Inject pre-computed proportions
model.results["cell_prop"] = prop_df

# Build cell_to_spot mapping (positional indices into prop_df)
spot_order = {bc: i for i, bc in enumerate(prop_df.index)}
nuc_map = nuc_map[nuc_map["spot_barcode"].isin(spot_order)].copy()
cell_ids = nuc_map["nucleus_id"].astype(str).values
cell_to_spot = nuc_map["spot_barcode"].map(spot_order).values.astype(int)
nuclei_counts = nuc_map.groupby("spot_barcode")["nucleus_id"].count().reindex(
    prop_df.index, fill_value=0
)

# Hungarian assignment: discretize spot proportions into per-nucleus types.
# (Bayesian mode is available via assignment_method="bayesian" only when a
# per-nucleus (C, n_types) morphology-score array is passed as
# morph_scores_precomputed; the post-E3 build ships no score producer.)
assignments_df = model.assign_cells(
    nuclei_counts=nuclei_counts,
    cell_to_spot=cell_to_spot,
    cell_ids=cell_ids,
    assignment_method="hungarian",
    detection_mask=np.ones((len(prop_df), len(prop_df.columns)), dtype=bool),
)

# SACE per-cell GEX allocation. output_mode="single_cell" auto-orchestrates
# StarDist + assignment + SACE projection behind an H&E resolution gate
# (requires <=1.0 um/px), so the assign_cells call above is illustrative.
spot_type_gex, cell_adata, _ = model.run_sace_allocation(
    output_mode="single_cell",
)

# Save final single-cell AnnData
cell_adata.write_h5ad(OUTPUT_DIR / f"{SAMPLE_NAME}_single_cell.h5ad")
print(f"Done! {cell_adata.n_obs} cells × {cell_adata.n_vars} genes")
print(f"Cell types: {cell_adata.obs['cell_type'].value_counts().to_dict()}")
```

> **In practice**, use the production script directly:
> ```bash
> python examples/scripts/run_single_cell_assignment.py --sample sample_P1_S1 --stages 1,5,6
> ```
> It handles all the data loading, alignment, and edge cases (NaN coords, missing embeddings, image fallbacks, etc.).

**Final output:** `{SAMPLE_NAME}_single_cell.h5ad`

| Field | Contents |
|-------|----------|
| `.X` | Per-cell gene expression (count-scale SACE allocations) |
| `.obs["cell_type"]` | Assigned cell type |
| `.obs["spot_barcode"]` | Parent Visium spot |
| `.obs["cell_id"]` | Unique cell/nucleus identifier |
| `.obsm["spatial"]` | (x, y) pixel coordinates |

This AnnData is compatible with standard single-cell tools (scanpy, scVI, CellxGene, etc.).

---

## SLURM Submission

All steps require significant compute. Example sbatch wrappers:

```bash
#!/bin/bash
#SBATCH --job-name=citegeist_m3
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your_email@institution.edu
#SBATCH --output=slurm_%j.log

module load cuda/12.1
conda activate CITEgeist_env

cd /path/to/CITEgeist/CITEgeist

# Step 2: Proportions
python ../examples/scripts/run_module3_unified.py --sample sample_P1_S1 \
    --output-dir ../output/module3_unified

# Step 3: Single-cell assignment + GEX
python ../examples/scripts/run_single_cell_assignment.py --sample sample_P1_S1 \
    --stages 1,5,6
```

For GPU jobs on CRC, use `gpu-race-submit.sh` instead of direct `sbatch` to race across available GPU partitions.

---

## Quick Reference: Key Parameters

| Parameter | Default | When to Change |
|-----------|---------|----------------|
| `min_counts` | 100 | Lower to 25 for surgical/low-depth samples |
| `method` | `"qp"` | Only option for production (cuOPT) |
| `use_gex_detection` | `True` | Disable if antibody panel is very sparse |
| `assignment_method` | `"hungarian"` | Use `"bayesian"` only with precomputed per-nucleus morphology scores |
| `max_iter` (SACE) | `1` | Do not increase — EM overfits beyond iteration 1 |

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| cuOPT import fails | You're on a CPU node — must run on GPU with `--gres=gpu:1` |
| `obsm["spatial"]` missing | Use `load_images=True` in `sq.read.visium()`, or inject manually from `tissue_positions.csv` |
| Antibody names don't match profile | Real data has `-1` suffix (e.g., `CD68-1`); simulated data does not |
| NaN coordinates | Surgical samples may have NaN spatial coords — filter explicitly after `split_adata()` |
| StarDist OOM | Use `resolution_mode="hires"` or reduce patch size |
| HDF5 locking errors | Set `HDF5_USE_FILE_LOCKING=FALSE` in sbatch environment |
| No fullres TIF | Pipeline falls back to CytAssist or hires image (see `_resolve_he_image` in production script) |

---

## What's Next?

With your per-cell AnnData, you can:

- **Spatial programs** (Module 4): Discover gene programs anchored to specific cell types
- **Cross-sample integration** (Module 5): Align programs across patients
- **Downstream analysis**: Use standard scanpy workflows (clustering, DE, trajectory) on the single-cell output
- **Visualization**: Plot cell-type maps, pie-chart spots, or spatial gene expression

See `examples/scripts/` for downstream analysis examples.
