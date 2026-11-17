"""Expand spot-level deconvolution into pseudo-single-cell AnnData.

For each spot, allocates CELLS_PER_SPOT (5) pseudo-cells based on
proportions from C2L or Tangram. Each pseudo-cell gets the GEX from
its cell type's layer, jittered spatial coordinates within the spot
radius, and a cell_type / spot_barcode annotation.
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import CELLS_PER_SPOT, OUTPUT_ROOT, SAMPLES, SPOT_RADIUS_UM


def largest_remainder_allocate(proportions: np.ndarray, total: int) -> np.ndarray:
    """Allocate integer counts summing to *total* using largest-remainder method.

    Parameters
    ----------
    proportions : 1-D array of non-negative floats (should sum to ~1)
    total : desired integer total

    Returns
    -------
    1-D int array with same length, summing to *total*.
    """
    scaled = proportions * total
    floored = np.floor(scaled).astype(int)
    remainders = scaled - floored
    deficit = total - floored.sum()
    # Distribute remaining slots to the types with the largest remainders
    if deficit > 0:
        indices = np.argsort(-remainders)[:deficit]
        floored[indices] += 1
    return floored


def expand_sample(method: str, sample: str) -> None:
    """Build pseudo-single-cell AnnData for one method+sample pair."""
    sample_dir = OUTPUT_ROOT / method / sample

    # Handle naming variations across deconvolution scripts
    prop_candidates = [
        sample_dir / f"{sample}_proportions.csv",
        sample_dir / "proportions.csv",
    ]
    gex_candidates = [
        sample_dir / f"{sample}_gex_layers.h5ad",
        sample_dir / "gex_layers.h5ad",
        sample_dir / "c2l_gex_layers.h5ad",
    ]
    prop_path = next((p for p in prop_candidates if p.exists()), None)
    gex_path = next((p for p in gex_candidates if p.exists()), None)

    if prop_path is None:
        raise FileNotFoundError(f"Proportions not found in {sample_dir} (tried: {[p.name for p in prop_candidates]})")
    if gex_path is None:
        raise FileNotFoundError(f"GEX layers not found in {sample_dir} (tried: {[p.name for p in gex_candidates]})")

    # ------------------------------------------------------------------
    # Load inputs
    # ------------------------------------------------------------------
    prop_df = pd.read_csv(prop_path, index_col=0)  # spots x cell types
    adata_layers = ad.read_h5ad(gex_path)  # spots x genes, layers per type

    # Strip C2L prefixes from column/layer names to get clean cell type names
    C2L_PROP_PREFIX = "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"
    C2L_LAYER_PREFIX = "means_per_cluster_mu_fg_"

    def strip_prefix(name):
        for pfx in [C2L_PROP_PREFIX, C2L_LAYER_PREFIX]:
            if name.startswith(pfx):
                return name[len(pfx) :]
        return name

    # Build mappings: clean_name → original_name
    layer_map = {strip_prefix(k): k for k in adata_layers.layers.keys()}
    prop_map = {strip_prefix(c): c for c in prop_df.columns}

    common_types = sorted(set(layer_map.keys()) & set(prop_map.keys()))
    if not common_types:
        raise ValueError(
            f"No overlapping cell types between proportions ({set(prop_map.keys())}) "
            f"and layers ({set(layer_map.keys())})"
        )
    print(f"  Using {len(common_types)} common cell types: {common_types}")

    # Subset and renormalize proportions to common types
    prop_sub = prop_df[[prop_map[ct] for ct in common_types]].copy()
    prop_sub.columns = common_types
    row_sums = prop_sub.sum(axis=1)
    row_sums = row_sums.replace(0, 1)  # avoid division by zero
    prop_sub = prop_sub.div(row_sums, axis=0)

    # Intersect spots
    common_spots = prop_sub.index.intersection(adata_layers.obs_names)
    if len(common_spots) == 0:
        raise ValueError("No overlapping spot barcodes between proportions and GEX layers")
    prop_sub = prop_sub.loc[common_spots]
    adata_layers = adata_layers[common_spots]
    print(f"  {len(common_spots)} spots with both proportions and GEX")

    # Spatial coordinates
    if "spatial" not in adata_layers.obsm:
        raise ValueError("GEX layers AnnData missing .obsm['spatial']")
    spatial = adata_layers.obsm["spatial"]  # (n_spots, 2)

    # ------------------------------------------------------------------
    # Allocate pseudo-cells
    # ------------------------------------------------------------------
    import hashlib

    seed_bytes = hashlib.sha256(f"{method}_{sample}".encode()).digest()[:4]
    seed = int.from_bytes(seed_bytes, "big")
    rng = np.random.default_rng(seed=seed)
    n_genes = adata_layers.n_vars
    gene_names = adata_layers.var_names

    # Pre-extract layer matrices (dense for simpler slicing)
    layer_matrices = {}
    for ct in common_types:
        mat = adata_layers.layers[layer_map[ct]]
        if sp.issparse(mat):
            mat = mat.toarray()
        layer_matrices[ct] = mat  # (n_spots, n_genes)

    # Build lists for output
    all_X = []
    all_cell_type = []
    all_barcode = []
    all_spatial = []

    prop_values = prop_sub.values  # (n_spots, n_types)

    for i in range(len(common_spots)):
        props_i = prop_values[i]
        counts = largest_remainder_allocate(props_i, CELLS_PER_SPOT)
        spot_xy = spatial[i]  # (2,)

        for j, ct in enumerate(common_types):
            n_cells = counts[j]
            if n_cells < 1:
                continue
            gex_vec = layer_matrices[ct][i]  # (n_genes,)
            for _ in range(n_cells):
                # Uniform jitter within spot disc
                angle = rng.uniform(0, 2 * np.pi)
                radius = SPOT_RADIUS_UM * np.sqrt(rng.uniform(0, 1))
                jitter = np.array([radius * np.cos(angle), radius * np.sin(angle)])
                all_X.append(gex_vec)
                all_cell_type.append(ct)
                all_barcode.append(common_spots[i])
                all_spatial.append(spot_xy + jitter)

    # ------------------------------------------------------------------
    # Assemble AnnData
    # ------------------------------------------------------------------
    X = np.vstack(all_X)
    obs = pd.DataFrame(
        {
            "cell_type": all_cell_type,
            "spot_barcode": all_barcode,
        }
    )
    obs.index = [f"cell_{k}" for k in range(len(obs))]
    spatial_out = np.vstack(all_spatial)

    adata_out = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=gene_names),
    )
    adata_out.obsm["spatial"] = spatial_out

    out_path = sample_dir / f"{sample}_single_cell.h5ad"
    adata_out.write_h5ad(out_path)
    print(f"  Saved {adata_out.n_obs} pseudo-cells to {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Expand spot deconvolution to pseudo-single-cell AnnData")
    parser.add_argument(
        "--method",
        required=True,
        choices=["cell2location", "tangram"],
        help="Deconvolution method",
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample ID (must be in constants.SAMPLES)",
    )
    args = parser.parse_args()

    if args.sample not in SAMPLES:
        raise ValueError(f"Unknown sample {args.sample!r}; expected one of {SAMPLES}")

    print(f"Expanding pseudo-cells: method={args.method}, sample={args.sample}")
    expand_sample(args.method, args.sample)
    print("Done.")


if __name__ == "__main__":
    main()
