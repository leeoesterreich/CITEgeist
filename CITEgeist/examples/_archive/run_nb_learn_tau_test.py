#!/usr/bin/env python3
"""
Quick test: learned tau on P4-S2_1i_rep.
Reports learned tau, cancer std, r(EPCAM).
Single start for speed.
"""
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

from CITEgeist.model.nb_deconvolution import run_nb_optimization
from CITEgeist.model.nb_initialization import (
    build_lam_prior_sigma_matrix,
    initialize_background,
    initialize_dispersion,
    initialize_lambda,
)

DATA_ROOT = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")

PATIENT_CELL_TYPES = [
    "Cancer", "Macrophages", "CD8_T_Cells", "CD4_T_Cells",
    "B_Cells", "Endothelial", "Fibroblasts", "Monocytes", "Dendritic_Cells",
]

PATIENT_MARKER_TYPE_TABLE = {
    "EPCAM":   [("Cancer", "strong")],
    "CD68":    [("Macrophages", "strong")],
    "CD163":   [("Macrophages", "strong")],
    "CD3E":    [("CD8_T_Cells", "strong"), ("CD4_T_Cells", "strong")],
    "CD8A":    [("CD8_T_Cells", "strong")],
    "CD4":     [("CD4_T_Cells", "strong")],
    "CD19":    [("B_Cells", "strong")],
    "PECAM1":  [("Endothelial", "strong")],
    "ACTA2":   [("Fibroblasts", "strong")],
    "CD14":    [("Monocytes", "strong")],
    "ITGAX":   [("Dendritic_Cells", "strong"), ("Macrophages", "soft")],
    "HLA-DRA": [("Dendritic_Cells", "strong"), ("B_Cells", "soft"), ("Macrophages", "soft")],
}


def main():
    import squidpy as sq
    import scipy.sparse

    sample = "HCC22-088-P4-S2_1i_rep"
    visium_path = DATA_ROOT / sample / "outs"

    print(f"Loading {sample}...")
    adata = sq.read.visium(str(visium_path), counts_file='filtered_feature_bc_matrix.h5',
                           load_images=True, gex_only=False)
    adata.var_names_make_unique()

    is_ab = adata.var['feature_types'] == 'Antibody Capture'
    ab_data = adata[:, is_ab].copy()
    ab_X = ab_data.X.toarray() if scipy.sparse.issparse(ab_data.X) else np.array(ab_data.X)
    ab_df = pd.DataFrame(ab_X, index=ab_data.obs_names, columns=ab_data.var_names)
    epcam_raw = ab_df['EPCAM-1'] if 'EPCAM-1' in ab_df.columns else None

    # Build marker config
    type_names = list(PATIENT_CELL_TYPES)
    T = len(type_names)

    def canonical(name):
        base = name.replace("-1", "").replace("-", "")
        remap = {"CD31": "PECAM1", "aSMA": "ACTA2", "CD11c": "ITGAX", "CD20": "CD19", "panCK": None}
        if base in remap:
            return remap[base]
        if base in PATIENT_MARKER_TYPE_TABLE:
            return base
        return name

    markers = []
    raw_indices = []
    seen = set()
    for idx, raw_name in enumerate(ab_data.var_names):
        can = canonical(raw_name)
        if can is None:
            continue
        if can in PATIENT_MARKER_TYPE_TABLE and can not in seen:
            markers.append(can)
            raw_indices.append(idx)
            seen.add(can)

    M = len(markers)
    active_mask = np.zeros((T, M), dtype=bool)
    for m_idx, marker in enumerate(markers):
        for type_name, _strength in PATIENT_MARKER_TYPE_TABLE[marker]:
            if type_name in type_names:
                t_idx = type_names.index(type_name)
                active_mask[t_idx, m_idx] = True

    S = ab_X[:, raw_indices]
    I = S.shape[0]
    N = np.ones(I)

    # QP warm start
    qp_file = REPO_ROOT / "output/module3_nb" / sample / f"{sample}_cell_prop_global_results.csv"
    qp_props = pd.read_csv(qp_file, index_col=0)
    p_init = np.zeros((I, T))
    for t_idx, ct in enumerate(type_names):
        if ct in qp_props.columns:
            p_init[:, t_idx] = qp_props.loc[adata.obs_names[:I], ct].values
    row_sums = p_init.sum(axis=1, keepdims=True)
    p_init = np.where(row_sums > 0, p_init / row_sums, 1.0 / T)

    median_N = float(np.median(N))
    lam_init = initialize_lambda(S, active_mask, median_N)
    b_init = initialize_background(S)
    r_init = initialize_dispersion(S)

    lam_sigma_matrix = build_lam_prior_sigma_matrix(
        active_mask, markers=markers, type_names=type_names,
        marker_type_table=PATIENT_MARKER_TYPE_TABLE,
        sigma_exclusive=5.0, sigma_shared=2.0, sigma_inactive=1.0,
    )

    # Test configs
    configs = [
        {"name": "baseline (tau=1.0 fixed)", "tau": 1.0, "learn_tau": False},
        {"name": "learn_tau (init=1.0)", "tau": 1.0, "learn_tau": True, "tau_prior_sigma": 1.0},
        {"name": "learn_tau (init=0.5)", "tau": 0.5, "learn_tau": True, "tau_prior_sigma": 1.0},
        {"name": "learn_tau weak prior (init=1.0, sigma=2.0)", "tau": 1.0, "learn_tau": True, "tau_prior_sigma": 2.0},
    ]

    for cfg in configs:
        print(f"\n{'='*60}")
        print(f"Config: {cfg['name']}")
        print(f"{'='*60}")

        t0 = time.time()
        result = run_nb_optimization(
            S=S, N=N, active_mask=active_mask, n_iter=100,
            lam_init=lam_init, r_init=r_init, b_init=b_init, p_init=p_init,
            lam_prior_sigma=2.0,
            lam_prior_sigma_matrix=lam_sigma_matrix,
            lr_global=0.02, lr_local=0.05, patience=5, n_starts=1,
            tau=cfg['tau'],
            learn_tau=cfg.get('learn_tau', False),
            tau_prior_sigma=cfg.get('tau_prior_sigma', 1.0),
            device="cpu",
        )
        elapsed = time.time() - t0

        props = pd.DataFrame(result["proportions"], index=adata.obs_names[:I], columns=type_names)
        cancer = props["Cancer"]

        learned_tau = result.get("tau", cfg["tau"])
        r_epcam = stats.pearsonr(cancer, epcam_raw.loc[cancer.index])[0] if epcam_raw is not None else np.nan

        print(f"  Learned tau:   {learned_tau:.4f}")
        print(f"  Cancer std:    {cancer.std():.4f}")
        print(f"  Cancer >0.5:   {(cancer > 0.5).sum()}")
        print(f"  Cancer >0.3:   {(cancer > 0.3).sum()}")
        print(f"  r(EPCAM raw):  {r_epcam:.3f}")
        print(f"  Iters:         {len(result['held_out_nll'])}")
        print(f"  Time:          {elapsed:.0f}s")


if __name__ == "__main__":
    main()
