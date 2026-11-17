"""Task A8: QP proportions + SACE `layers` runs (real, init=prop, three nulls) for one specimen.

Clinical config follows examples/scripts/run_cuopt_qp_patient.py (10-type major-marker dictionary,
lambda_reg=1.0, lambda_laplacian=0.03, detection gating 0.01), with the real GSE289326 data root.
Writes per specimen: props.csv, layers_<variant>.parquet (types x genes layer means), run_meta.json, .done
"""
import argparse
import hashlib
import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE.parent))
from ora_check import ora_lib as L  # noqa: E402

PROTOCOL = json.loads((HERE / "ora_protocol.json").read_text())

# Same dictionary as examples/scripts/run_cuopt_qp_patient.py (major markers only).
PROFILE_DICT = {
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "CD8_T_Cells": {"Major": ["CD8A-1", "CD3E-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1", "CD3E-1"]},
    "Cancer_Luminal": {"Major": ["EPCAM-1"]},
    "Cancer_Basal": {"Major": ["KRT5-1", "SDC1-1", "EPCAM-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1", "HLA-DRA-1"]},
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("a8")


def load_model(sample: str, out_dir: Path, run_qp: bool = True, profile: dict = None):
    import squidpy as sq

    from CITEgeist.model import CitegeistModel

    patient_dir = Path(os.environ.get("CITEGEIST_DATA_ROOT", "data")) / PROTOCOL["data_root"] / sample / "outs"
    adata = sq.read.visium(
        str(patient_dir), counts_file="filtered_feature_bc_matrix.h5", load_images=True, gex_only=False
    )
    model = CitegeistModel(sample_name=sample, adata=adata, output_folder=str(out_dir))
    model.split_adata()
    for attr in ("gene_expression_adata", "antibody_capture_adata"):
        ad = getattr(model, attr)
        if ad is not None and "spatial" in ad.obsm:
            valid = np.all(np.isfinite(ad.obsm["spatial"]), axis=1)
            if (~valid).sum() > 0:
                log.warning("Filtering %d spots with NaN coords from %s", (~valid).sum(), attr)
                setattr(model, attr, ad[valid].copy())
    min_counts = PROTOCOL["clinical_config"]["min_counts"]["S1" if "S1" in sample else "S2"]
    model.filter_gex(min_counts=min_counts)
    model.copy_gex_to_protein_adata()  # restrict antibody spots to the filtered GEX spots (as in the original pipeline)
    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(PROFILE_DICT if profile is None else profile)
    if not run_qp:
        return model
    qp = PROTOCOL["clinical_config"]["qp"]
    model.run_cell_proportion_model(
        method=qp["method"],
        lambda_reg=qp["lambda_reg"],
        lambda_laplacian=qp["lambda_laplacian"],
        use_detection_gating=qp["use_detection_gating"],
        detection_gate_ub=qp["detection_gate_ub"],
    )
    return model


def variant_inputs(name: str, real_props: pd.DataFrame, seed: int, profile: dict = None):
    """Return (props, dict, init_method) for a variant."""
    profile = PROFILE_DICT if profile is None else profile
    init = "prop" if name.endswith("_prop") else "marker"
    if name.startswith("real"):
        return real_props, profile, init
    if name.startswith("spot_permuted"):
        return L.spot_permuted(real_props, seed)[0], profile, init
    if name.startswith("uniform"):
        return L.uniform_props(real_props), profile, init
    if name.startswith("scrambled"):
        return real_props, L.scrambled_dict(profile, seed), init
    if name.startswith("blend"):  # blend<W>_marker: (1-W/100) * real + W/100 * uniform, rows renormalized (ORA sense check)
        w = int(name[len("blend") :].split("_")[0]) / 100.0
        q = (1.0 - w) * real_props + w * L.uniform_props(real_props)
        tot = q.sum(axis=1).replace(0, 1.0)
        return q.div(tot, axis=0), profile, init
    raise ValueError(name)


def run_variants(model, real_props, profile, seed, out_dir, gene_names, variants):
    """SACE `layers` for each variant; writes layers_<variant>.parquet; returns the per-variant meta."""
    out = {}
    for name in variants:
        props, cpd, init = variant_inputs(name, real_props, seed, profile)
        model.results["cell_prop"] = props
        model.results.pop("gene_expression_pass1", None)
        model.cell_profile_dict = cpd
        spot_type_gex, diag = model.run_sace_allocation(output_mode="layers", init_method=init)
        acc = np.zeros((props.shape[1], len(gene_names)))
        for s in spot_type_gex:
            acc += spot_type_gex[s]
        means = acc / len(spot_type_gex)
        pd.DataFrame(means, index=props.columns, columns=gene_names).to_parquet(out_dir / f"layers_{name}.parquet")
        out[name] = {
            "init_method": init,
            "log_likelihood": diag.get("log_likelihood"),
            "total_mean_count": float(means.sum()),
        }
        log.info("variant %s done (init=%s)", name, init)
        del spot_type_gex, acc, means
    model.results["cell_prop"] = real_props
    return out


def load_v41_profile() -> dict:
    """BREAST_PROFILE of CITEgeist/examples/breast_dictionary.py (single source), loaded by path."""
    import importlib.util

    f = REPO / "CITEgeist" / "examples" / "breast_dictionary.py"
    spec = importlib.util.spec_from_file_location("_breast_dictionary_single_source", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.BREAST_PROFILE


def run_arms(a, sample, seed, t0):
    """Supplied-proportion mode (no QP): every ``--arm NAME=CSV`` is run on the strict barcode intersection of the arms.

    Writes <out-root>/<arm>/sace_layers/<sample>/{props.csv, layers_<variant>.parquet, run_meta.json, .done}.
    """
    from ora_check import ora_v2_lib as V

    profile = load_v41_profile() if a.profile == "v41" else PROFILE_DICT
    arm_props = {}
    for spec in a.arm:
        name, path = spec.split("=", 1)
        p = pd.read_csv(path.replace("{sample}", sample), index_col=0)
        arm_props[name] = p.drop(columns="recon_error") if "recon_error" in p.columns else p
    types = list(profile)
    for name, p in arm_props.items():
        assert list(p.columns) == types, f"arm {name}: columns {list(p.columns)} != dictionary types {types}"
    Path(a.out_root).mkdir(parents=True, exist_ok=True)
    model = load_model(sample, Path(a.out_root), run_qp=False, profile=profile)
    barcodes = list(model.gene_expression_adata.obs_names)
    common = V.align_arms(arm_props, barcodes, reference_arm=a.subset_arm)
    n_model = len(barcodes)
    model.gene_expression_adata = model.gene_expression_adata[common].copy()
    model.antibody_capture_adata = model.antibody_capture_adata[common].copy()
    assert list(model.gene_expression_adata.obs_names) == common
    gene_names = list(model.gene_expression_adata.var_names)
    variants = a.variants.split(",") if a.variants else list(PROTOCOL["variants"])
    for name, p in arm_props.items():
        if a.write_arms and name not in a.write_arms.split(","):
            continue
        real = p.loc[common]
        assert real.index.equals(model.gene_expression_adata.obs_names)
        out_dir = Path(a.out_root) / name / "sace_layers" / sample
        out_dir.mkdir(parents=True, exist_ok=True)
        real.to_csv(out_dir / "props.csv")
        meta = {
            "sample": sample,
            "arm": name,
            "seed": seed,
            "n_spots": int(len(common)),
            "n_arm_table_spots": int(len(p)),
            "n_model_spots": n_model,
            "n_arm_spots_dropped": int(len(p) - len(common)),
            "types": types,
            "profile": a.profile,
        }
        meta["variants"] = run_variants(model, real, profile, seed, out_dir, gene_names, variants)
        meta["runtime_s"] = time.time() - t0
        meta["protocol_sha256"] = hashlib.sha256((HERE / "ora_protocol.json").read_bytes()).hexdigest()
        (out_dir / "run_meta.json").write_text(json.dumps(meta, indent=2))
        (out_dir / ".done").touch()
        log.info("DONE arm=%s %s n_spots=%d", name, sample, len(common))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-index", type=int, required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--race-ns", default="")
    ap.add_argument("--arm", action="append", default=[], help="NAME=CSV (use {sample}); supplied proportions, skips QP")
    ap.add_argument("--subset-arm", default=None, help="arm whose barcodes must all be present in every other arm")
    ap.add_argument("--profile", choices=["legacy", "v41"], default="legacy")
    ap.add_argument("--write-arms", default="", help="comma list of arms to write (default all); alignment still uses every --arm")
    ap.add_argument("--variants", default="", help="comma list (default: all protocol variants)")
    a = ap.parse_args()

    sample = PROTOCOL["specimens"][a.sample_index]
    seed = PROTOCOL["seeds"]["base"] + a.sample_index
    if a.arm:
        log.info("sample=%s seed=%d arms=%s out=%s", sample, seed, a.arm, a.out_root)
        run_arms(a, sample, seed, time.time())
        return
    out_dir = Path(a.out_root) / (f"{sample}_race_{a.race_ns}" if a.race_ns else sample)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    log.info("sample=%s seed=%d out=%s", sample, seed, out_dir)

    model = load_model(sample, out_dir)
    real_props = model.results["cell_prop"].copy()
    if "recon_error" in real_props.columns:
        real_props = real_props.drop(columns="recon_error")
    real_props.to_csv(out_dir / "props.csv")
    gene_names = list(model.gene_expression_adata.var_names)
    log.info("QP done: %d spots x %d types, %d genes", *real_props.shape, len(gene_names))

    meta = {
        "sample": sample,
        "seed": seed,
        "n_spots": int(real_props.shape[0]),
        "types": list(real_props.columns),
    }
    meta["variants"] = run_variants(model, real_props, None, seed, out_dir, gene_names, list(PROTOCOL["variants"]))

    meta["runtime_s"] = time.time() - t0
    meta["protocol_sha256"] = hashlib.sha256((HERE / "ora_protocol.json").read_bytes()).hexdigest()
    (out_dir / "run_meta.json").write_text(json.dumps(meta, indent=2))
    (out_dir / ".done").touch()
    log.info("DONE %s in %.0fs", sample, meta["runtime_s"])


if __name__ == "__main__":
    main()
