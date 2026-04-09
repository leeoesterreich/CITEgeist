#!/usr/bin/env python3
"""
Full fix test: QP without Laplacian → NB with learn_tau + adaptive sigma.
Runs on P4-S2_1i_rep, reports cancer std, r(EPCAM), learned tau.

Tests two conditions:
  A) Old: QP(Laplacian=0.1) → NB(tau=1.0 fixed) [baseline]
  B) New: QP(Laplacian=0.0) → NB(learn_tau=True, adaptive sigma)
"""
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import scipy.sparse

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

import squidpy as sq
from model import CitegeistModel
from model.unified_config import CELL_PROFILES_NESTED
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


def run_qp_warm_start(sample, lambda_laplacian):
    """Run QP on a sample with specified Laplacian strength. Returns proportions."""
    visium_path = DATA_ROOT / sample / "outs"
    adata = sq.read.visium(
        str(visium_path),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False,
    )
    adata.var_names_make_unique()

    model = CitegeistModel(
        sample_name=sample, adata=adata,
        output_folder=str(REPO_ROOT / "output" / "nb_fix_test"),
    )
    model.split_adata()
    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(CELL_PROFILES_NESTED)

    model.run_cell_proportion_model(
        validation_warn_only=True,
        skip_finetuning=True,
        lambda_laplacian=lambda_laplacian,
    )

    return model.results["cell_prop"], model


def build_nb_inputs(model, sample):
    """Extract S, N, active_mask, p_init from model for NB."""
    ab_adata = model.antibody_capture_adata
    ab_X = ab_adata.X.toarray() if scipy.sparse.issparse(ab_adata.X) else np.array(ab_adata.X)

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
    for idx, raw_name in enumerate(ab_adata.var_names):
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

    return S, N, active_mask, markers, type_names


def run_nb_phase(S, N, active_mask, markers, type_names, qp_props, spot_names,
                 learn_tau, tau_init, lam_sigma_matrix):
    """Run NB optimization with given settings."""
    T = len(type_names)
    I = S.shape[0]

    p_init = np.zeros((I, T))
    for t_idx, ct in enumerate(type_names):
        if ct in qp_props.columns:
            p_init[:, t_idx] = qp_props.loc[spot_names[:I], ct].values
    row_sums = p_init.sum(axis=1, keepdims=True)
    p_init = np.where(row_sums > 0, p_init / row_sums, 1.0 / T)

    median_N = float(np.median(N))
    lam_init = initialize_lambda(S, active_mask, median_N)
    b_init = initialize_background(S)
    r_init = initialize_dispersion(S)

    result = run_nb_optimization(
        S=S, N=N, active_mask=active_mask, n_iter=100,
        lam_init=lam_init, r_init=r_init, b_init=b_init, p_init=p_init,
        lam_prior_sigma=2.0,
        lam_prior_sigma_matrix=lam_sigma_matrix,
        lr_global=0.02, lr_local=0.05, patience=5, n_starts=1,
        tau=tau_init,
        learn_tau=learn_tau,
        tau_prior_sigma=1.0,
        device="cpu",
    )
    return result


def evaluate(result, type_names, spot_names, epcam_raw):
    """Compute metrics on NB result."""
    props = pd.DataFrame(result["proportions"], index=spot_names, columns=type_names)
    cancer = props["Cancer"]

    common = cancer.index.intersection(epcam_raw.index)
    r_epcam, _ = stats.pearsonr(cancer.loc[common], epcam_raw.loc[common])

    return {
        'cancer_std': cancer.std(),
        'cancer_gt05': (cancer > 0.5).sum(),
        'cancer_gt03': (cancer > 0.3).sum(),
        'r_epcam_raw': r_epcam,
        'tau': result.get('tau', None),
        'n_iters': len(result['held_out_nll']),
    }


def main():
    sample = "HCC22-088-P4-S2_1i_rep"
    print(f"{'='*70}")
    print(f"FULL FIX TEST: {sample}")
    print(f"{'='*70}")

    # Get raw EPCAM for validation
    visium_path = DATA_ROOT / sample / "outs"
    adata_raw = sq.read.visium(str(visium_path), counts_file='filtered_feature_bc_matrix.h5',
                                load_images=True, gex_only=False)
    adata_raw.var_names_make_unique()
    is_ab = adata_raw.var['feature_types'] == 'Antibody Capture'
    ab_raw = adata_raw[:, is_ab].copy()
    ab_X = ab_raw.X.toarray() if scipy.sparse.issparse(ab_raw.X) else np.array(ab_raw.X)
    ab_df = pd.DataFrame(ab_X, index=ab_raw.obs_names, columns=ab_raw.var_names)
    epcam_raw = ab_df['EPCAM-1']

    # Build adaptive sigma matrix (same for both conditions)
    type_names = list(PATIENT_CELL_TYPES)
    # We'll build it after we get active_mask from the first run

    # === Condition A: Old baseline — reuse existing QP output (with Laplacian) ===
    print(f"\n{'='*60}")
    print("A) BASELINE: QP(Laplacian=0.1) → NB(tau=1.0 fixed)")
    print(f"{'='*60}")
    t0 = time.time()

    # Reuse QP from the previous NB run (already has Laplacian=0.1)
    qp_file_A = REPO_ROOT / "output/module3_nb" / sample / f"{sample}_cell_prop_global_results.csv"
    qp_props_A = pd.read_csv(qp_file_A, index_col=0)
    print(f"  Reusing existing QP from {qp_file_A.name}")

    # Still need model for ab data extraction — load once, reuse
    visium_path_model = DATA_ROOT / sample / "outs"
    adata_model = sq.read.visium(str(visium_path_model), counts_file='filtered_feature_bc_matrix.h5',
                                  load_images=True, gex_only=False)
    adata_model.var_names_make_unique()
    model_A = CitegeistModel(sample_name=sample, adata=adata_model,
                              output_folder=str(REPO_ROOT / "output" / "nb_fix_test"))
    model_A.split_adata()

    S, N, active_mask, markers, type_names = build_nb_inputs(model_A, sample)
    spot_names = model_A.antibody_capture_adata.obs_names

    # No adaptive sigma for baseline
    result_A = run_nb_phase(S, N, active_mask, markers, type_names, qp_props_A,
                            spot_names, learn_tau=False, tau_init=1.0,
                            lam_sigma_matrix=None)
    metrics_A = evaluate(result_A, type_names, spot_names, epcam_raw)
    elapsed_A = time.time() - t0
    print(f"  Cancer std:    {metrics_A['cancer_std']:.4f}")
    print(f"  Cancer >0.5:   {metrics_A['cancer_gt05']}")
    print(f"  Cancer >0.3:   {metrics_A['cancer_gt03']}")
    print(f"  r(EPCAM raw):  {metrics_A['r_epcam_raw']:.3f}")
    print(f"  tau:           {metrics_A['tau']}")
    print(f"  Time:          {elapsed_A:.0f}s")

    # === Condition B: Full fix (Laplacian=0.0, learn_tau, adaptive sigma) ===
    print(f"\n{'='*60}")
    print("B) FULL FIX: QP(Laplacian=0.0) → NB(learn_tau=True, adaptive sigma)")
    print(f"{'='*60}")
    t0 = time.time()
    qp_props_B, model_B = run_qp_warm_start(sample, lambda_laplacian=0.0)
    S_B, N_B, active_mask_B, markers_B, type_names_B = build_nb_inputs(model_B, sample)
    spot_names_B = model_B.antibody_capture_adata.obs_names

    lam_sigma_matrix = build_lam_prior_sigma_matrix(
        active_mask_B, markers=markers_B, type_names=type_names_B,
        marker_type_table=PATIENT_MARKER_TYPE_TABLE,
        sigma_exclusive=5.0, sigma_shared=2.0, sigma_inactive=1.0,
    )

    result_B = run_nb_phase(S_B, N_B, active_mask_B, markers_B, type_names_B, qp_props_B,
                            spot_names_B, learn_tau=True, tau_init=1.0,
                            lam_sigma_matrix=lam_sigma_matrix)
    metrics_B = evaluate(result_B, type_names_B, spot_names_B, epcam_raw)
    elapsed_B = time.time() - t0
    print(f"  Cancer std:    {metrics_B['cancer_std']:.4f}")
    print(f"  Cancer >0.5:   {metrics_B['cancer_gt05']}")
    print(f"  Cancer >0.3:   {metrics_B['cancer_gt03']}")
    print(f"  r(EPCAM raw):  {metrics_B['r_epcam_raw']:.3f}")
    print(f"  Learned tau:   {metrics_B['tau']:.4f}")
    print(f"  Time:          {elapsed_B:.0f}s")

    # === Condition C: No Laplacian only (tau=1.0 fixed, no adaptive sigma) ===
    print(f"\n{'='*60}")
    print("C) ABLATION: QP(Laplacian=0.0) → NB(tau=1.0 fixed, no adaptive sigma)")
    print(f"{'='*60}")
    t0 = time.time()
    # Reuse qp_props_B from above
    result_C = run_nb_phase(S_B, N_B, active_mask_B, markers_B, type_names_B, qp_props_B,
                            spot_names_B, learn_tau=False, tau_init=1.0,
                            lam_sigma_matrix=None)
    metrics_C = evaluate(result_C, type_names_B, spot_names_B, epcam_raw)
    elapsed_C = time.time() - t0
    print(f"  Cancer std:    {metrics_C['cancer_std']:.4f}")
    print(f"  Cancer >0.5:   {metrics_C['cancer_gt05']}")
    print(f"  Cancer >0.3:   {metrics_C['cancer_gt03']}")
    print(f"  r(EPCAM raw):  {metrics_C['r_epcam_raw']:.3f}")
    print(f"  tau:           {metrics_C['tau']}")
    print(f"  Time:          {elapsed_C:.0f}s")

    # === Summary ===
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"  {'Condition':<55} {'std':>7} {'r(EPCAM)':>9} {'tau':>6}")
    print(f"  {'-'*55} {'-'*7} {'-'*9} {'-'*6}")
    print(f"  {'A) QP(Lap=0.1) + NB(tau=1.0)':<55} {metrics_A['cancer_std']:>7.4f} {metrics_A['r_epcam_raw']:>9.3f} {'1.000':>6}")
    print(f"  {'B) QP(Lap=0.0) + NB(learn_tau, adaptive sigma)':<55} {metrics_B['cancer_std']:>7.4f} {metrics_B['r_epcam_raw']:>9.3f} {metrics_B['tau']:>6.3f}")
    print(f"  {'C) QP(Lap=0.0) + NB(tau=1.0, no sigma)':<55} {metrics_C['cancer_std']:>7.4f} {metrics_C['r_epcam_raw']:>9.3f} {'1.000':>6}")


if __name__ == "__main__":
    main()
