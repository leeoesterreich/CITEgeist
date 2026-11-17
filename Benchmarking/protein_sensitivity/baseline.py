"""Production call sequences for the v40 stress test (ported from recovered scripts).

Heavy imports (CitegeistModel, SACE, evaluators) are done lazily inside the
functions that need them so that this module imports cheaply on CPU.
"""
from __future__ import annotations

import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))
from benchmark_constants import (  # noqa: E402
    ACHIEVABLE_6_CELL_PROFILE_DICT,
    ACHIEVABLE_7_CELL_PROFILE_DICT,
)

XENIUM_DATA = REPO / "Benchmarking/xenium_pseudovisium/data_protein_gt"
XB = REPO / "Benchmarking/xenium_benchmarking"
GT_DIR = {
    6: XB / "ground_truth_singler/ground_truth_6type/ground_truth",
    7: XB / "ground_truth_singler/ground_truth_7type/ground_truth",
}
GEX_GT_DIR = XB / "ground_truth_singler/ground_truth_7type/ground_truth_gex"
PROFILE = {6: ACHIEVABLE_6_CELL_PROFILE_DICT, 7: ACHIEVABLE_7_CELL_PROFILE_DICT}
DEFAULT_SIM_ROOT = REPO / "Benchmarking/simulation_benchmarking/replicates"


def resolve_sim_root(env) -> Path:
    """Sim replicate root: env SIM_ROOT if set, else the canonical `replicates` dir (unchanged).

    An explicit SIM_ROOT must exist and hold both condition dirs (high_seg, mixed), else ValueError. The default is
    returned unchecked, exactly as before (it is untracked, so a clean worktree legitimately lacks it at import time).
    """
    raw = env.get("SIM_ROOT")
    if not raw:
        return DEFAULT_SIM_ROOT
    root = Path(raw)
    if not root.is_dir():
        raise ValueError(f"SIM_ROOT {root} is not an existing directory")
    missing = [c for c in ("high_seg", "mixed") if not (root / c).is_dir()]
    if missing:
        raise ValueError(f"SIM_ROOT {root} lacks condition dir(s): {missing}")
    return root


def sim_unit_paths(root: Path, condition: str, rep: int) -> dict:
    """Per-unit sim inputs under `root` (used by eval_unit.py and tests)."""
    root = Path(root)
    return {
        "gt_prop": root / condition / "ST_sim" / f"Wu_ST_{rep}_prop.csv",
        "gt_gex_layers": root / condition / "ST_GEX_sim" / f"sample_{rep}" / "layers",
        "gex_h5ad": root / condition / "h5ad_objects" / f"Wu_rep_{rep}_GEX.h5ad",
        "cite_h5ad": root / condition / "h5ad_objects" / f"Wu_rep_{rep}_CITE.h5ad",
    }


SIM_ROOT = resolve_sim_root(os.environ)

_XENIUM_SRC = XB / "CITEgeist" / "src"
_XENIUM_EVAL_SRC = XB / "evaluation" / "src"
_SIM_SRC = REPO / "Benchmarking" / "simulation_benchmarking" / "src"
_SIM_SACE_SRC = REPO / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "src"


def _add_path(p: Path) -> None:
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def load_xenium(region: int):
    d = XENIUM_DATA / "h5ad_objects"
    return (
        sc.read_h5ad(d / f"Xenium_region_{region}_GEX.h5ad"),
        sc.read_h5ad(d / f"Xenium_region_{region}_CITE.h5ad"),
    )


def load_sim(condition: str, rep: int):
    d = SIM_ROOT / condition / "h5ad_objects"
    return sc.read_h5ad(d / f"Wu_rep_{rep}_GEX.h5ad"), sc.read_h5ad(d / f"Wu_rep_{rep}_CITE.h5ad")


def sim_profile() -> dict:
    """Return SIMULATION_CELL_PROFILE_DICT (9 types, 2 Major markers each).

    Imported from benchmark_sace_simulation.py (module-level code there is
    only constants and logging.basicConfig).
    """
    _add_path(_SIM_SACE_SRC)
    from benchmark_sace_simulation import SIMULATION_CELL_PROFILE_DICT

    return SIMULATION_CELL_PROFILE_DICT


def align_prediction(pred: pd.DataFrame, types, renormalize: bool = True) -> pd.DataFrame:
    """Reindex to `types` (missing -> 0, extras dropped). With renormalize=True
    (Xenium path) also clip at 0 and row-normalize nonzero rows; with False
    (sim path, which the published pipeline does not renormalize) only
    reindex + zero-fill, leaving row sums untouched.
    """
    out = pred.reindex(columns=list(types)).fillna(0.0)
    if not renormalize:
        return out
    out = out.clip(lower=0)
    s = out.sum(1)
    out.loc[s > 0] = out.loc[s > 0].div(s[s > 0], axis=0)
    return out


def load_nuclei_counts(region: int) -> pd.Series:
    """Per-spot nuclei counts from the Xenium cell-to-spot mapping.

    Ported verbatim from recovered/run_qp_7type_benchmark.py:59-67
    (get_nuclei_counts); identical in run_qp_only_benchmark.py.
    """
    mapping = pd.read_csv(XENIUM_DATA / "cell_to_spot_mapping.csv", index_col=0)
    cell_types = pd.read_csv(XENIUM_DATA / "cell_type_assignments.csv", index_col=0)

    reg = mapping[(mapping["region_id"] == region) & (mapping["spot_idx"] != -1)]
    known = reg.merge(cell_types, left_index=True, right_index=True)
    known = known[known["cell_type"] != "Unknown"]
    return known.groupby("spot_id").size()


def _resolve_phase1(phase1=None) -> str:
    """Proportion method for Phase 1: explicit arg, else env PHASE1, else "qp" (the v40 pinned call).

    "marker_share_bright" = marker share with the brightness method named by env BRIGHTNESS_METHOD.
    """
    phase1 = phase1 if phase1 is not None else os.environ.get("PHASE1", "qp")
    if phase1 not in ("qp", "marker_share", "marker_share_bright"):
        raise ValueError(f"PHASE1 must be 'qp', 'marker_share' or 'marker_share_bright', got {phase1!r}")
    return phase1


def _bright_kwargs(panel_set: str) -> dict:
    """run_cell_proportion_model kwargs for PHASE1=marker_share_bright (method from env BRIGHTNESS_METHOD).

    The RNA estimation panel is the dataset's fixed panel, re-keyed to profile type names.
    """
    method = os.environ.get("BRIGHTNESS_METHOD")
    if not method:
        raise ValueError("PHASE1=marker_share_bright needs env BRIGHTNESS_METHOD (e.g. plain, gmm, rna_slope_adj)")
    from CITEgeist.model.deconvolution.brightness_estimators import (
        ESTIMATION_PANELS,
        PANEL_TO_PROFILE_TYPES,
    )

    panel = {
        t: list(ESTIMATION_PANELS[panel_set][key])
        for key, types in PANEL_TO_PROFILE_TYPES[panel_set].items()
        for t in types
    }
    return {"method": "marker_share", "brightness_method": method, "rna_panel": panel}


def _xenium_panel_set(profile) -> str:
    return "xenium7" if "CD4+ T cells" in profile else "xenium6"


_QP_ONLY_LOGGED: set = set()


def _log_qp_only_settings_once(setting) -> None:
    """marker_share path: log (once per param) the `setting` overrides that only tune the QP solver.

    `markers` is a data/profile change and IS honoured (resolve() promotes Minor to Major before
    load_cell_profile_dict). Everything else (solver kwargs, nuclei_source, detection-GMM variants, joint
    block) has no marker-share counterpart, so it is ignored, loudly.
    """
    if not setting or setting.get("axis") not in ("setting", "joint"):
        return
    param = setting.get("param")
    if param == "markers" or param in _QP_ONLY_LOGGED:
        return
    _QP_ONLY_LOGGED.add(param)
    import logging

    logging.warning(
        "PHASE1=marker_share: setting %r is a QP-solver/nuclei knob with no marker-share counterpart; "
        "ignored (run_xenium_qp results are identical to baseline for this condition).",
        param,
    )


def run_xenium_qp(adata_gex, adata_cite, profile, region, setting=None, phase1=None, variant=None) -> pd.DataFrame:
    """Xenium proportion QP (spots x profile types, row-normalized).

    Ported from recovered/run_qp_7type_benchmark.py:111-152 (and the identical
    6-type recovered/run_qp_only_benchmark.py): CitegeistModel(simulation=True)
    -> filter_gex(min_counts=25) -> preprocess_gex(target_sum=10000) ->
    preprocess_antibody() -> load_cell_profile_dict(profile) ->
    run_cell_proportion_model(method="qp", lambda_laplacian=0.0, nuclei_counts=...). Output
    is rows renormalized with sum.replace(0, 1), exactly as the recovered
    scripts do (no reindex-with-zero-fill: the recovered evaluation
    intersects spots with GT, so we keep the model's spots).

    `setting` (settings-sensitivity axis, None/baseline = the published call above) is
    {"axis": "setting", "param", "level"}; settings_conditions.resolve turns it into a
    profile and run_cell_proportion_model kwargs (lambda_laplacian=0.0 + nuclei_counts always
    present, as published).
    """
    from CITEgeist.model.citegeist_model import CitegeistModel
    from protein_sensitivity import settings_conditions as SC

    import contextlib

    profile, qp_kwargs = SC.resolve(setting, region, profile, load_nuclei_counts)
    # setting=None (stress-test path) is untouched. With a setting: benchmark-only protein-GMM variant (round 2)
    # applied for this condition only, plus a passive threshold recorder for diagnostics.
    recorder = SC.ThresholdRecorder() if setting is not None else None
    det_cm = (
        SC.detection_context(SC.detection_variant(setting), recorder) if setting is not None else contextlib.nullcontext()
    )
    diagnostics = None
    sample_name = f"Xenium_region_{region}"
    with tempfile.TemporaryDirectory() as tmpdir, det_cm:
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=tmpdir,
            simulation=True,
            gene_expression_adata=adata_gex.copy(),
            antibody_capture_adata=adata_cite.copy(),
        )
        model.filter_gex(min_counts=25, keep_unfiltered=_resolve_phase1(phase1) == "marker_share_bright")
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()
        model.load_cell_profile_dict(profile)
        if _resolve_phase1(phase1) == "marker_share_bright":
            _log_qp_only_settings_once(setting)
            model.run_cell_proportion_model(**_bright_kwargs(_xenium_panel_set(profile)))
        elif _resolve_phase1(phase1) == "marker_share":
            _log_qp_only_settings_once(setting)
            if variant is not None and variant.brightness:
                model.run_cell_proportion_model(method="marker_share", brightness_correction=True)
            else:
                model.run_cell_proportion_model(method="marker_share")
        else:
            model.run_cell_proportion_model(method="qp", **qp_kwargs)
        props = model.results["cell_prop"].copy()
        if setting is not None:
            diagnostics = SC.collect_diagnostics(model.results, qp_kwargs, recorder)
        if variant is not None and variant.brightness:
            diagnostics = {
                **(diagnostics or {}),
                "brightness": model.results.get("marker_brightness_diagnostics"),
                "refit_on": "perturbed_raw_counts",
            }
    props = props.drop(columns=["recon_error"], errors="ignore")
    types = list(profile.keys())
    for ct in types:
        if ct not in props.columns:
            props[ct] = 0.0
    props = props[types]
    out = props.div(props.sum(axis=1).replace(0, 1), axis=0)
    if diagnostics is not None:
        out.attrs["diagnostics"] = diagnostics
    return out


def run_xenium_sace(adata_gex, adata_cite, profile, region, phase1=None, variant=None) -> dict:
    """Xenium SACE layers: {TYPE_TO_FILENAME stem: spots x genes DataFrame}.

    Ported from benchmark_sace_xenium.py:145-281. Differences (all
    parameterisation, not method): the profile is the argument instead of the
    fixed 7-type dict; the all-types list is the 7-type dict's keys (GEX GT is
    7-type); output is returned instead of written. Only meaningful for
    7-type profiles.
    """
    _add_path(_XENIUM_SRC)
    from benchmark_sace_xenium import TYPE_TO_FILENAME, build_cell_spot_map

    from CITEgeist.model.citegeist_model import CitegeistModel
    from CITEgeist.model.gex.sace_gex import run_sace

    sample = f"Xenium_region_{region}"
    all_types = list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys())
    with tempfile.TemporaryDirectory() as tmpdir:
        model = CitegeistModel(
            sample_name=sample,
            output_folder=tmpdir,
            simulation=True,
            gene_expression_adata=adata_gex.copy(),
            antibody_capture_adata=adata_cite.copy(),
        )
        model.filter_gex(min_counts=25, keep_unfiltered=_resolve_phase1(phase1) == "marker_share_bright")
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()

        major_only_dict = {ct: {"Major": info.get("Major", []), "Minor": []} for ct, info in profile.items()}
        model.load_cell_profile_dict(major_only_dict)
        if _resolve_phase1(phase1) == "marker_share_bright":
            model.run_cell_proportion_model(**_bright_kwargs(_xenium_panel_set(profile)))
        elif _resolve_phase1(phase1) == "marker_share":
            if variant is not None and variant.brightness:
                model.run_cell_proportion_model(method="marker_share", brightness_correction=True)
            else:
                model.run_cell_proportion_model(method="marker_share")
        else:
            model.run_cell_proportion_model(
                method="qp",
                lambda_reg=1.0,
                use_detection_gating=True,
                detection_gate_ub=0.01,
            )

        props = model.results["cell_prop"].copy()
        if "recon_error" in props.columns:
            props = props.drop(columns=["recon_error"])
        for ct in all_types:
            if ct not in props.columns:
                props[ct] = 0.0
        props = props[all_types]
        spot_names = list(model.gene_expression_adata.obs_names)
        props = props.reindex(spot_names, fill_value=0.0)
        props = props.div(props.sum(axis=1).replace(0, 1), axis=0)

        raw_adata = model.gene_expression_adata
        if "raw_counts" in raw_adata.layers:
            X_raw = raw_adata.layers["raw_counts"]
        else:
            X_raw = raw_adata.X
        if hasattr(X_raw, "toarray"):
            X_raw = X_raw.toarray()
        X_raw = np.asarray(X_raw, dtype=np.float64)
        gene_names = list(raw_adata.var_names)
        if "spatial" in raw_adata.obsm:
            coords = np.asarray(raw_adata.obsm["spatial"], dtype=np.float64)
        else:
            coords = np.zeros((len(spot_names), 2), dtype=np.float64)

        cell_assignments, cell_spot_map = build_cell_spot_map(props, spot_coords=coords, cells_per_spot=5)
        cite_raw = np.asarray(model.antibody_capture_adata.layers.get("raw_counts", model.antibody_capture_adata.X))
        cite_names = list(model.antibody_capture_adata.var_names)

        _, sace_adata, _diagnostics = run_sace(
            spot_counts=X_raw,
            proportions=props,
            cell_assignments=cell_assignments,
            cell_spot_map=cell_spot_map,
            spot_coords=coords,
            gene_names=gene_names,
            init_method="marker",
            antibody_data=cite_raw,
            antibody_names=cite_names,
            cell_profile_dict=profile,
        )

        sace_X = sace_adata.X.toarray() if hasattr(sace_adata.X, "toarray") else np.asarray(sace_adata.X)
        sace_X = np.asarray(sace_X, dtype=np.float64)
        cid_to_spot = dict(zip(cell_spot_map["cell_id"].values, cell_spot_map["spot_barcode"].values))
        layers = {}
        for type_name in all_types:
            if "cell_type" in sace_adata.obs.columns:
                mask = sace_adata.obs["cell_type"].values == type_name
            else:
                mask = np.array([cell_assignments.get(cid, "") == type_name for cid in sace_adata.obs_names])
            if mask.sum() == 0:
                continue
            type_expr = sace_X[mask, :]
            cell_ids_in_type = list(np.array(sace_adata.obs_names)[mask])
            spot_to_rows = defaultdict(list)
            for row_i, cell_id in enumerate(cell_ids_in_type):
                spot_name = cid_to_spot.get(cell_id)
                if spot_name is not None:
                    spot_to_rows[spot_name].append(row_i)
            spot_list = sorted(spot_to_rows.keys())
            if not spot_list:
                continue
            spot_matrix = np.zeros((len(spot_list), len(gene_names)), dtype=np.float64)
            for si, spot in enumerate(spot_list):
                spot_matrix[si] = type_expr[spot_to_rows[spot], :].sum(axis=0)
            layers[TYPE_TO_FILENAME[type_name]] = pd.DataFrame(spot_matrix, index=spot_list, columns=gene_names)
    return layers


def run_sim(adata_gex, adata_cite, profile, sace: bool, phase1=None):
    """Simulation QP (+ optional SACE layers). Returns (props_df, layers|None).

    Ported from benchmark_sace_simulation.py:356-433 (non-cached branch): no
    filter_gex; model.cell_profile_dict = profile; run_cell_proportion_model(method="qp")
    with defaults; props reindexed to GEX obs_names with NaN rows -> uniform;
    results["cell_prop"] = props; then run_sace_allocation(output_mode="layers",
    init_method="marker", integer_allocation=False) and spot_type_gex_to_layers.
    The synthetic cell-spot map built in the source is unused by
    run_sace_allocation and is omitted.
    """
    _add_path(_SIM_SACE_SRC)
    from benchmark_sace_simulation import spot_type_gex_to_layers

    from CITEgeist.model.citegeist_model import CitegeistModel

    with tempfile.TemporaryDirectory() as tmpdir:
        model = CitegeistModel(
            sample_name="Wu_rep",
            output_folder=tmpdir,
            simulation=True,
            gene_expression_adata=adata_gex.copy(),  # CitegeistModel keeps references and preprocessing is in place
            antibody_capture_adata=adata_cite.copy(),
        )
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()
        model.cell_profile_dict = profile
        raw_phase1 = phase1 if phase1 is not None else os.environ.get("PHASE1", "qp")
        if raw_phase1 == "marker_share_v41p":  # v41p = marker share + per-sample GMM brightness correction (sim only)
            global_props_df, _finetuned = model.run_cell_proportion_model(
                method="marker_share", brightness_correction=True
            )
        elif raw_phase1 == "marker_share_bright":
            global_props_df, _finetuned = model.run_cell_proportion_model(**_bright_kwargs("sim"))
        else:
            global_props_df, _finetuned = model.run_cell_proportion_model(method=_resolve_phase1(phase1))
        brightness_diag = model.results.get("marker_brightness_diagnostics")
        props_df = global_props_df.drop(columns=["recon_error"], errors="ignore")
        props_df = props_df.reindex(model.gene_expression_adata.obs_names)
        props_df = props_df.fillna(1.0 / props_df.shape[1])
        model.results["cell_prop"] = props_df
        if brightness_diag is not None:
            props_df.attrs["diagnostics"] = {"brightness": brightness_diag, "refit_on": "sample_raw_counts"}
        if not sace:
            return props_df, None
        spot_type_gex, _diag = model.run_sace_allocation(
            output_mode="layers", init_method="marker", integer_allocation=False
        )
        layers = spot_type_gex_to_layers(
            spot_type_gex,
            list(model.gene_expression_adata.obs_names),
            list(model.gene_expression_adata.var_names),
            list(props_df.columns),
        )
    return props_df, layers


# ---------------------------------------------------------------------------
# Metric wrappers (return dicts shaped like the fakes in the runner tests)
# ---------------------------------------------------------------------------


def xenium_prop_metrics(pred: pd.DataFrame, region: int, n_types: int, gt_root=None) -> dict:
    """Proportion metrics vs SingleR GT via evaluate_benchmark.calculate_metrics.

    Spots are intersected with GT (as in recovered evaluate_against_gt).
    Maps overall_pearson_r/overall_rmse/overall_mae/mean_jsd to
    pearson_r/rmse/mae/jsd; per_type = per-type Pearson r.
    """
    _add_path(_XENIUM_EVAL_SRC)
    from evaluate_benchmark import calculate_metrics

    gt = pd.read_csv(
        (GT_DIR[n_types] if gt_root is None else Path(gt_root)) / f"Xenium_region_{region}_prop.csv", index_col=0
    )
    types = list(PROFILE[n_types].keys())
    aligned = align_prediction(pred, types)
    common = aligned.index.intersection(gt.index)
    m = calculate_metrics(gt.loc[common, types], aligned.loc[common], types)
    return {
        "pearson_r": float(m["overall_pearson_r"]),
        "rmse": float(m["overall_rmse"]),
        "mae": float(m["overall_mae"]),
        "jsd": float(m["mean_jsd"]),
        "per_type": {ct: float(m[f"{ct}_pearson_r"]) for ct in types if f"{ct}_pearson_r" in m},
    }


def xenium_gex_metrics(layers: dict, region: int, gex_gt_dir=None) -> dict:
    """GEX metrics vs GT layers; mean over cell types of flattened_r / nrmse.

    Per-type prep ported from compare_xenium_gex_methods.py:181-250
    (uppercase genes, shared spots/genes, GT-expressed genes), then
    compute_metrics(pred, gt).
    """
    _add_path(_XENIUM_EVAL_SRC)
    from compare_xenium_gex_methods import compute_metrics, load_csv_spots_genes

    flat, nrmse = [], []
    gex_gt_root = GEX_GT_DIR if gex_gt_dir is None else Path(gex_gt_dir)
    for stem, pred_df in layers.items():
        gt_path = gex_gt_root / f"Xenium_region_{region}" / f"{stem}_GT.csv"
        if not gt_path.exists():
            continue
        gt_df = load_csv_spots_genes(gt_path)
        pred_df = pred_df.copy()
        gt_df.columns = gt_df.columns.str.upper()
        pred_df.columns = pred_df.columns.str.upper()
        gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
        pred_df = pred_df.loc[:, ~pred_df.columns.duplicated()]
        spots = sorted(set(gt_df.index) & set(pred_df.index))
        genes = sorted(set(gt_df.columns) & set(pred_df.columns))
        if len(spots) < 5 or len(genes) < 5:
            continue
        g = gt_df.loc[spots, genes].values.astype(np.float64)
        p = pred_df.loc[spots, genes].values.astype(np.float64)
        keep = g.sum(axis=0) > 0
        if keep.sum() < 5:
            continue
        m = compute_metrics(p[:, keep], g[:, keep])
        flat.append(m["flattened_r"])
        nrmse.append(m["nrmse"])
    if not flat:
        return {"gex_flat_r": float("nan"), "gex_nrmse": float("nan")}
    return {"gex_flat_r": float(np.nanmean(flat)), "gex_nrmse": float(np.nanmean(nrmse))}


def sim_prop_metrics(pred: pd.DataFrame, condition: str, rep: int) -> dict:
    """Simulation proportion metrics via benchmarking_spot_deconv.benchmark_performance.

    Mapping: sim `corr` -> pearson_r, `Sum_RMSE` -> rmse, `Sum_MAE` -> mae,
    median `JSD` -> jsd. per_type = per-type Pearson r. GT is
    ST_sim/Wu_ST_{rep}_prop.csv (spot_x/spot_y dropped; sorted rows/cols as
    in the module's main()).
    """
    _add_path(_SIM_SRC)
    from benchmarking_spot_deconv import benchmark_performance

    gt = pd.read_csv(SIM_ROOT / condition / "ST_sim" / f"Wu_ST_{rep}_prop.csv", index_col=0)
    gt = gt.drop(columns=["spot_x", "spot_y"], errors="ignore").sort_index().sort_index(axis=1)
    cols = list(gt.columns)
    aligned = align_prediction(pred, cols, renormalize=False).reindex(gt.index).fillna(0.0)
    res = benchmark_performance(gt.values, aligned.values, cols)
    per_type = {}
    for c in cols:
        a, b = gt[c].values, aligned[c].values
        per_type[c] = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float("nan")
    return {
        "pearson_r": float(res["corr"]),
        "rmse": float(res["Sum_RMSE"]),
        "mae": float(res["Sum_MAE"]),
        "jsd": float(res["JSD"]),
        "per_type": per_type,
    }


def sim_gex_metrics(layers: dict, condition: str, rep: int) -> dict:
    """Simulation GEX metrics via benchmarking_gex.calculate_rmse (normalize="range").

    Layers are written as {type}_layer.csv into a temp dir (the layout
    calculate_rmse expects). gex_flat_r = mean_pearson_r, gex_nrmse = average_nrmse.
    """
    _add_path(_SIM_SRC)
    from benchmarking_gex import calculate_rmse

    gt_dir = SIM_ROOT / condition / "ST_GEX_sim" / f"sample_{rep}" / "layers"
    with tempfile.TemporaryDirectory() as tmpdir:
        for ct, df in layers.items():
            df.to_csv(Path(tmpdir) / f"{ct}_layer.csv")
        m = calculate_rmse(ground_truth_dir=str(gt_dir), predictions_dir=tmpdir, normalize="range")
    return {"gex_flat_r": float(m["mean_pearson_r"]), "gex_nrmse": float(m["average_nrmse"])}
