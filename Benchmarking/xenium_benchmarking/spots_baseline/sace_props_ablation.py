"""SACE proportion-source ablation on Xenium RCC (7-type): does SACE's per-cell-type GEX depend on the QP proportions?

One code path, three arms; ONLY the proportion input to SACE differs:
  qp    : canonical cuOPT QP proportions (identical call to CITEgeist/src/benchmark_sace_xenium.py; needs a GPU)
  qp_saved : the saved QP proportions of the protein_sensitivity baseline (CITEgeist_7type_results.json _inputs.run_dirs;
            a DIFFERENT QP call: lambda_laplacian=0 + nuclei_counts) -- CPU only, no re-solve
  spots : SPOTS reading (a) marker-share proportions (spots_baseline.spots_composition, ACHIEVABLE_7 Major markers)
  gt    : SingleR 7-type ground-truth proportions (oracle)
  marker_share : v41 CITEgeist marker-share props (output_marker_share_7type, env MAIN); output is FLAT (OUT_DIR/Xenium_region_N)
Everything else (filter_gex(min_counts=25), preprocess, Major-only profile dict, init_method="marker", 5 cells/spot
largest-remainder rounding, run_sace, per-type layer export) is the canonical benchmark_sace_xenium.run_region flow.

Detection mask: the canonical flow does NOT pass the QP detection mask to SACE (detection gating acts inside the QP only);
SACE receives just the (spots x 7) proportion matrix, so no extra handling is needed for SPOTS/GT. A zero proportion
already means "no cells of this type" via largest-remainder rounding. Spots with an all-zero proportion row (possible for
GT = spots with no known cell) are left as-is (identical to what the canonical flow does with a zero QP row); their count
is recorded in meta.json.

Usage: sace_props_ablation.py ARM REGION OUT_DIR [PROPS_DIR [PROFILE]]
       -> OUT_DIR/ARM/Xenium_region_REGION/{layers/*.csv, props.csv, meta.json}
  marker_share arm only: PROPS_DIR (default <MAIN>/.../CITEgeist/output_marker_share_7type) holds Xenium_region_N/predictions/baseline.csv.gz;
  PROFILE in {achievable (default), ms_v41} selects the dict used BOTH for the types/Major-only model profile and the SACE marker init.

Cell mode (run_region(cell_mode=...), 6th CLI arg):
  five     (default, reproduces every old arm): each spot is given a fixed 5 pseudo-cells by largest-remainder rounding of the
           proportions, so a type that rounds to 0 cells in a spot has no row in its layer and its mass is moved onto the
           other types of that spot (see rna_refinement/cell_count_DECISION_RULE.md).
  no_round : `run_sace_layers` directly, no cell rounding; the layer of type t is E_x[s,t,:] for every spot with nonzero
           mass, so no type loses its mass. Same inputs, same marker init, same seeds as `five`.
`cells_per_spot=5` is a fixed placeholder for the spot-level layers, NOT a nuclei count. Single-cell output
(`output_mode="single_cell"`) is a separate path that still rounds to the real nuclei counts.
"""
import json
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
GT_PROP = XB / "ground_truth_singler" / "ground_truth_7type" / "ground_truth"
ARMS = ("qp", "qp_saved", "spots", "gt", "marker_share")
CELL_MODES = ("five", "no_round")
FLAT_ARMS = ("marker_share",)  # these write OUT_DIR/Xenium_region_N directly (no <arm>/ level)


def align_props(df, spot_names, types):
    """Reindex to (spot_names x types), missing -> 0, then row-normalize (zero rows stay zero). Same tail as the canonical flow."""
    out = df.reindex(index=list(spot_names), columns=list(types)).fillna(0.0).astype(float)
    return out.div(out.sum(axis=1).replace(0, 1), axis=0)


def spots_arm_props(adata_cite, profile, spot_names):
    """SPOTS reading (a) on the RAW CITE counts (the .X of the h5ad, as in spots_baseline.py), aligned to the SACE spots."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from spots_baseline import spots_composition

    X = adata_cite.X.toarray() if hasattr(adata_cite.X, "toarray") else np.asarray(adata_cite.X)
    df, _ = spots_composition(X, adata_cite.var_names, profile, "a")
    df.index = adata_cite.obs_names
    return align_props(df, spot_names, list(profile.keys()))


def qp_saved_props(region, spot_names, types):
    rd = json.loads((XB / "evaluation/results/method_comparison_singler_7type/CITEgeist_7type_results.json").read_text())
    path = REPO / rd["_inputs"]["run_dirs"][str(region)] / "predictions" / "baseline.csv.gz"
    df = pd.read_csv(path, index_col=0).drop(columns=["recon_error"], errors="ignore")
    return align_props(df[list(types)], spot_names, types)


def get_profile(name):
    """SACE profile dict by name (7-type): achievable (default, v41) or ms_v41 (v41p: B Major = CD20 + CD138)."""
    sys.path.insert(0, str(XB))
    import benchmark_constants as C

    return {"achievable": C.ACHIEVABLE_7_CELL_PROFILE_DICT, "ms_v41": C.MS_V41_7_CELL_PROFILE_DICT}[name]


def marker_share_props(region, spot_names, types, props_dir=None):
    """Marker-share proportions written by write_marker_share_preds.py (env MAIN = repo root holding the gitignored predictions)."""
    main = Path(os.environ.get("MAIN", REPO))
    root = Path(props_dir) if props_dir else main / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_marker_share_7type"
    path = root / f"Xenium_region_{region}" / "predictions" / "baseline.csv.gz"
    df = pd.read_csv(path, index_col=0).drop(columns=["recon_error"], errors="ignore")
    return align_props(df[list(types)], spot_names, types)


def arm_region_dir(out_dir, arm, sample):
    return Path(out_dir) / sample if arm in FLAT_ARMS else Path(out_dir) / arm / sample


def gt_arm_props(region, spot_names, types):
    gt = pd.read_csv(GT_PROP / f"Xenium_region_{region}_prop.csv", index_col=0)
    return align_props(gt[list(types)], spot_names, types)


def no_round_layers(E_x, types, spot_names, gene_names):
    """{type: DataFrame(spots x genes) or None}: layer of type t = E_x[s,t,:] for spots with nonzero mass (no cell rounding)."""
    out = {}
    for i, t in enumerate(types):
        m = E_x[:, i, :].sum(axis=1) > 0
        out[t] = pd.DataFrame(E_x[m, i, :], index=[spot_names[j] for j in np.where(m)[0]], columns=gene_names) if m.any() else None
    return out


def run_region(arm, region, out_dir, props_dir=None, profile="achievable", cell_mode="five"):
    assert arm in ARMS, arm
    assert cell_mode in CELL_MODES, cell_mode
    assert profile == "achievable" or arm == "marker_share", "non-default profile only for the marker_share arm"
    import scanpy as sc

    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(XB))
    sys.path.insert(0, str(XB / "CITEgeist" / "src"))
    from benchmark_sace_xenium import TYPE_TO_FILENAME, build_cell_spot_map

    from CITEgeist.model.citegeist_model import CitegeistModel
    from CITEgeist.model.gex.sace_gex import run_sace, run_sace_layers

    PROFILE = get_profile(profile)
    sample = f"Xenium_region_{region}"
    types = list(PROFILE.keys())
    adata_gex = sc.read_h5ad(DATA / f"{sample}_GEX.h5ad")
    adata_cite = sc.read_h5ad(DATA / f"{sample}_CITE.h5ad")
    with tempfile.TemporaryDirectory() as tmpdir:
        model = CitegeistModel(
            sample_name=sample,
            output_folder=tmpdir,
            simulation=True,
            gene_expression_adata=adata_gex.copy(),
            antibody_capture_adata=adata_cite.copy(),
        )
        model.filter_gex(min_counts=25)
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()
        major_only = {ct: {"Major": info["Major"], "Minor": []} for ct, info in PROFILE.items()}
        model.load_cell_profile_dict(major_only)
        spot_names = list(model.gene_expression_adata.obs_names)

        # ---- the ONLY arm-specific step: where the proportions come from ----
        if arm == "qp":
            model.run_cell_proportion_model(method="qp", lambda_reg=1.0, use_detection_gating=True, detection_gate_ub=0.01)
            props = model.results["cell_prop"].copy().drop(columns=["recon_error"], errors="ignore")
            for ct in types:
                if ct not in props.columns:
                    props[ct] = 0.0
            props = align_props(props[types], spot_names, types)
        elif arm == "qp_saved":
            props = qp_saved_props(region, spot_names, types)
        elif arm == "marker_share":
            props = marker_share_props(region, spot_names, types, props_dir)
        elif arm == "spots":
            props = spots_arm_props(adata_cite, PROFILE, spot_names)
        else:
            props = gt_arm_props(region, spot_names, types)
        # ---------------------------------------------------------------------

        raw = model.gene_expression_adata
        X_raw = raw.layers["raw_counts"] if "raw_counts" in raw.layers else raw.X
        X_raw = np.asarray(X_raw.toarray() if hasattr(X_raw, "toarray") else X_raw, dtype=np.float64)
        gene_names = list(raw.var_names)
        if "spatial" in raw.obsm:
            coords = np.asarray(raw.obsm["spatial"], dtype=np.float64)
        else:
            coords = np.zeros((len(spot_names), 2))
        cite_raw = np.asarray(model.antibody_capture_adata.layers.get("raw_counts", model.antibody_capture_adata.X))
        sace_kwargs = dict(
            spot_counts=X_raw,
            proportions=props,
            spot_coords=coords,
            gene_names=gene_names,
            init_method="marker",
            antibody_data=cite_raw,
            antibody_names=list(model.antibody_capture_adata.var_names),
            cell_profile_dict=PROFILE,
        )
        reg_dir = arm_region_dir(out_dir, arm, sample)
        (reg_dir / "layers").mkdir(parents=True, exist_ok=True)
        if cell_mode == "no_round":
            _, E_x, _, diag = run_sace_layers(**sace_kwargs)
            for t, lay in no_round_layers(E_x, types, spot_names, gene_names).items():
                if lay is None:
                    lay = pd.DataFrame(index=pd.Index([], name="spot_id"), columns=gene_names)
                else:
                    lay.index.name = "spot_id"
                lay.to_csv(reg_dir / "layers" / f"{TYPE_TO_FILENAME[t]}_layer.csv")
            n_cells = None
        else:
            cell_assignments, cell_spot_map = build_cell_spot_map(props, spot_coords=coords, cells_per_spot=5)
            _, sace_adata, diag = run_sace(cell_assignments=cell_assignments, cell_spot_map=cell_spot_map, **sace_kwargs)
            sx = sace_adata.X
            sace_X = sx.toarray() if hasattr(sx, "toarray") else np.asarray(sx, dtype=np.float64)
            cid_to_spot = dict(zip(cell_spot_map["cell_id"].values, cell_spot_map["spot_barcode"].values))
            obs_ids = np.array(sace_adata.obs_names)
            for t in types:
                mask = sace_adata.obs["cell_type"].values == t
                if mask.sum() == 0:
                    # No cells of this type assigned (count-constrained SACE rounds a tiny
                    # proportion to 0). Write an EMPTY layer so the evaluator records
                    # insufficient_overlap (allowed) instead of layer_missing (fatal): the GT
                    # still defines the type, so this is a genuine no-overlap block, not a
                    # missing prediction file.
                    empty = pd.DataFrame(index=pd.Index([], name="spot_id"), columns=gene_names)
                    empty.to_csv(reg_dir / "layers" / f"{TYPE_TO_FILENAME[t]}_layer.csv")
                    continue
                spot_rows = defaultdict(list)
                for i, cid in enumerate(obs_ids[mask]):
                    spot_rows[cid_to_spot[cid]].append(i)
                te = sace_X[mask, :]
                spots = sorted(spot_rows)
                mat = np.stack([te[spot_rows[s]].sum(axis=0) for s in spots])
                lay = pd.DataFrame(mat, index=spots, columns=gene_names)
                lay.index.name = "spot_id"
                lay.to_csv(reg_dir / "layers" / f"{TYPE_TO_FILENAME[t]}_layer.csv")
            n_cells = len(cell_assignments)
        props.to_csv(reg_dir / "props.csv")
        meta = {
            "arm": arm,
            "region": region,
            "profile": profile,
            "cell_mode": cell_mode,
            "props_dir": str(props_dir) if props_dir else None,
            "n_spots": len(spot_names),
            "n_zero_prop_rows": int((props.sum(axis=1) == 0).sum()),
            "n_cells": n_cells,
            "log_likelihood": float(diag.get("log_likelihood", float("nan"))),
            "mean_prop_by_type": props.mean().round(4).to_dict(),
        }
        (reg_dir / "meta.json").write_text(json.dumps(meta, indent=2))
        print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    run_region(
        sys.argv[1],
        int(sys.argv[2]),
        sys.argv[3],
        *(sys.argv[4:5] or [None]),
        *(sys.argv[5:6] or ["achievable"]),
        *(sys.argv[6:7] or ["five"]),
    )
