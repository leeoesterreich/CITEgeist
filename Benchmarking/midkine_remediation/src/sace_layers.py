"""A12.6 producer: SACE layers (spot x type allocated counts) for one module3_7type specimen.

    model.run_sace_allocation(output_mode="layers", init_method="marker", integer_allocation=False)

Same model preparation as ``examples/scripts/run_module3_unified.py`` (split, spot filter ``min_counts`` 100 for -S1 / 25 for -S2,
copy_gex_to_protein_adata, preprocess_gex, preprocess_antibody, the 7-type marker dictionary, stored finetuned proportions) EXCEPT the gene
filter: here ``filter_gex(nonzero_percentage=0.0, mean_expression_threshold=-1.0)`` keeps EVERY gene. Reason: the module3 gene filter (>= 1% nonzero
spots AND mean over nonzero spots > 1) removes KRT6A in all 12 specimens (and most of the other keratins in several), and the protocol's
``k_share`` needs all six keratins and MDK in the layers (``compartments.MissingGeneError``). Gene-set choice is outcome-blind.

Per specimen output (written to a temp dir, atomically renamed; non-empty target refused; namespace <protocol_sha12>/<run_id>):
    layer_<Type>.npz   scipy CSR float64, spots x genes (expected allocated counts)
    genes.txt, barcodes.txt, coords.csv (barcode, x, y), spot_umi.csv (barcode, raw GEX UMI over all genes), nuclei.csv (B2 stamped StarDist counts, read-only)
    manifest.json      conservation gate, input md5s, protocol sha256, git head, seconds
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(HERE.parent))
import chain_guard as PF  # noqa: E402
import provenance  # noqa: E402
import protocol_io as PIO  # noqa: E402
import nuclei_artifact as NA  # noqa: E402
from breast_dict_loader import profile_for_protocol  # noqa: E402

CONSERVATION_MAX_REL = 1e-4


def md5_file(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def conservation(E_sparse_by_type: dict, Y: sp.spmatrix) -> dict:
    """Gate: sum over types of allocated counts equals the RAW counts entrywise; relative error |sum-Y| / max(Y, 1).

    ``Y`` is read independently from the Visium h5 (not from the model). Fails closed on a shape mismatch.
    """
    total = None
    for E in E_sparse_by_type.values():
        if E.shape != Y.shape:
            return {"passed": False, "reason": f"shape {E.shape} != raw {Y.shape}", "max_rel_error": float("inf")}
        total = E.copy() if total is None else total + E
    diff = (total - Y).tocoo()
    yv = np.asarray(Y[diff.row, diff.col]).ravel() if diff.nnz else np.array([])
    rel = np.abs(diff.data) / np.maximum(yv, 1.0) if diff.nnz else np.array([0.0])
    mx = float(rel.max()) if rel.size else 0.0
    return {"passed": bool(mx <= CONSERVATION_MAX_REL), "max_rel_error": mx, "limit": CONSERVATION_MAX_REL, "reason": ""}


def read_raw_counts(h5_path: Path, barcodes, genes) -> sp.csr_matrix:
    """Independent raw-count read: spots x genes (Gene Expression features only), rows/cols in the given order."""
    import h5py

    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        bc = np.array([b.decode() for b in m["barcodes"][:]])
        ft = m["features"]["feature_type"][:]
        names = np.array([b.decode() for b in m["features"]["name"][:]])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])).tocsr()   # genes x spots
    ge = np.where(ft == b"Gene Expression")[0]
    X, names = X[ge], names[ge]
    # duplicated symbols are made unique by the model as name-1, name-2: reproduce for the lookup
    seen, uniq = {}, []
    for n in names:
        k = seen.get(n, 0)
        uniq.append(n if k == 0 else f"{n}-{k}")
        seen[n] = k + 1
    gidx = pd.Index(uniq).get_indexer(list(genes))
    bidx = pd.Index(bc).get_indexer(list(barcodes))
    if (gidx < 0).any() or (bidx < 0).any():
        raise KeyError("barcodes/genes of the model are absent from the raw h5")
    return X[gidx][:, bidx].T.tocsr().astype(np.float64)


def produce(specimen, data_root, module3_dir, nuclei_dir, out_root, run_id, protocol_path):  # pragma: no cover - real data, sbatch only
    NA.refuse_standin(nuclei_dir)  # first: no v2 stage reads output/morphology_assignment_v3
    import squidpy as sq

    import rebuild_single_cell as RB
    from CITEgeist.model.citegeist_model import CitegeistModel

    proto = PIO.load_protocol(protocol_path)
    if tuple(Path(out_root).parts[-2:]) != (proto["_protocol_sha12"], run_id):
        raise PIO.OutputTargetError(f"out-root must end in <{proto['_protocol_sha12']}>/<{run_id}>: {out_root}")
    target = Path(out_root) / specimen
    PIO.check_output_target(target, proto["_protocol_sha12"], run_id)
    t0 = time.time()
    ent = RB.SPECIMENS[specimen]
    outs = Path(data_root) / ent["visium_dir"]
    qp_path = Path(module3_dir) / ent["qp_csv"]
    prop_df = pd.read_csv(qp_path, index_col=0)
    assert list(prop_df.columns) == proto["inputs"]["types"], "type order differs from the protocol"
    b2 = NA.load_nuclei(nuclei_dir, specimen, prop_df.index, protocol_sha12=proto["_protocol_sha12"])  # stamp + exact barcode equality
    h5 = outs / "filtered_feature_bc_matrix.h5"

    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    tmp = out_root / f".tmp_{specimen}_{os.getpid()}"
    tmp.mkdir()

    adata = sq.read.visium(str(outs), counts_file="filtered_feature_bc_matrix.h5", load_images=True, gex_only=False)
    if "spatial" in adata.obsm:
        ok = np.all(np.isfinite(adata.obsm["spatial"]), axis=1)
        if (~ok).any():
            adata = adata[ok].copy()
    model = CitegeistModel(sample_name=specimen, adata=adata, output_folder=str(tmp))
    model.split_adata()
    model.filter_gex(nonzero_percentage=0.0, mean_expression_threshold=-1.0, min_counts=100 if "-S1" in specimen else 25)   # spot filter only, ALL genes kept
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(profile_for_protocol(proto))
    common = prop_df.index.intersection(model.gene_expression_adata.obs_names)
    model.gene_expression_adata = model.gene_expression_adata[common].copy()
    model.antibody_capture_adata = model.antibody_capture_adata[common].copy()
    prop_df = prop_df.loc[common]
    # n_spots vs expected_n_spots_qp_file is reported in the manifest (the spot filter is the module3 rule)
    model.results["cell_prop"] = prop_df
    assert "raw_counts" in model.gene_expression_adata.layers, "SACE would fall back to normalised .X"

    spot_type_gex, diag = model.run_sace_allocation(output_mode="layers", init_method="marker", integer_allocation=False)
    genes = list(model.gene_expression_adata.var_names)
    barcodes = list(model.gene_expression_adata.obs_names)
    types = list(prop_df.columns)
    N, T, G = len(barcodes), len(types), len(genes)
    E_by_type = {}
    for ti, t in enumerate(types):
        dense = np.stack([spot_type_gex[s][ti] for s in range(N)])      # N x G
        E_by_type[t] = sp.csr_matrix(dense)
        del dense
    Y = read_raw_counts(h5, barcodes, genes)
    cons = conservation(E_by_type, Y)
    for t, E in E_by_type.items():
        sp.save_npz(tmp / f"layer_{t}.npz", E, compressed=True)
    Path(tmp / "genes.txt").write_text("\n".join(genes) + "\n")
    Path(tmp / "barcodes.txt").write_text("\n".join(barcodes) + "\n")
    xy = np.asarray(model.gene_expression_adata.obsm["spatial"], float)
    pd.DataFrame({"barcode": barcodes, "x": xy[:, 0], "y": xy[:, 1]}).to_csv(tmp / "coords.csv", index=False)
    pd.DataFrame({"barcode": barcodes, "spot_umi": np.asarray(Y.sum(axis=1)).ravel()}).to_csv(tmp / "spot_umi.csv", index=False)
    nuc_path = NA.nuclei_csv(nuclei_dir, specimen)
    b2.reindex(barcodes).rename("nuclei_count_raw").to_csv(tmp / "nuclei.csv", header=True)
    manifest = {
        "specimen": specimen, "run_id": run_id, "protocol_sha256": proto["_protocol_sha256"], "protocol_sha12": proto["_protocol_sha12"],
        "call": "run_sace_allocation(output_mode='layers', init_method='marker', integer_allocation=False); gene filter OFF (all genes), spot filter as module3",
        "n_spots": N, "n_types": T, "n_genes": G, "expected_n_spots_qp_file": proto["specimens"]["expected_n_spots_qp_file"][specimen],
        "types": types, "conservation": cons, "sace_diagnostics": {k: (float(v) if isinstance(v, (int, float, np.floating)) else str(v)) for k, v in diag.items()},
        "inputs_md5": {"qp_csv": md5_file(qp_path), "visium_h5": md5_file(h5), "tissue_positions": md5_file(outs / "spatial" / "tissue_positions.csv"), "b2_nuclei_csv": md5_file(nuc_path)},
        "git_head": provenance.code_commit(REPO_ROOT),
        "seconds": round(time.time() - t0, 1),
    }
    (tmp / "manifest.json").write_text(json.dumps(manifest, indent=2))
    if not cons["passed"]:
        failed = out_root / f".FAILED_{specimen}_{os.getpid()}"
        os.rename(tmp, failed)
        raise RuntimeError(f"conservation gate failed ({cons}); kept at {failed}")
    if target.exists():
        target.rmdir()
    os.rename(tmp, target)
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--specimen", required=True)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--module3-dir", type=Path, required=True, help="the B1 output dir (<specimen>_cell_prop_finetuned_results.csv)")
    ap.add_argument("--nuclei-dir", type=Path, required=True, help="the stamped B2 artifact dir (<root>/v2b2)")
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry("B4", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.module3_dir, nuclei_dir=a.nuclei_dir)
    produce(a.specimen, a.data_root, a.module3_dir, a.nuclei_dir, a.out_root, a.run_id, a.protocol)
    return 0


if __name__ == "__main__":
    sys.exit(main())
