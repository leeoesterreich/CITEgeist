"""Task 37 Step 1: write the fixed gene list for the fixed-gene Cell2Location reruns.

Setting ``full_training_5783`` (user decision, ``composition_arms/fixed_gene_set.json``): the list is the gene
index of the full-size v40 model's ``cell_type_signatures.csv`` (one gene per line, file order). The sidecar
records sha256, source path + sha256 and counts. The gate checks that, after the same MT filter and
intersection that ``Cell2Location/src/run_benchmark.py preprocess_spatial_data`` applies, the list leaves
exactly the genes of the old full-size inference (``Xenium_region_*_annotated.h5ad`` var_names) in all regions.
Refuses to overwrite. Login-node safe (reads only var names).
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
sys.path.insert(0, str(XB))

SETTING = "full_training_5783"
ALLOWED_SETTINGS = ("full_training_5783", "panel_180")
SRC_SIG = XB / "Cell2Location" / "reference_model" / "cell_type_signatures.csv"
OLD_RUN = XB / "Cell2Location" / "output_singler_6type"
GEX_DIR = XB.parent / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
OUT = XB / "composition_arms" / "fixed_genes_full_size.txt"
SETTING_FILE = XB / "composition_arms" / "fixed_gene_set.json"
N_REGIONS = 5


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def signature_genes(sig_csv) -> List[str]:
    """Gene index of the signature table (file order), rejecting empty/duplicated/blank names."""
    idx = pd.read_csv(sig_csv, index_col=0).index.astype(str).tolist()
    if not idx:
        raise ValueError(f"no genes in {sig_csv}")
    if any(not g.strip() for g in idx):
        raise ValueError("blank gene name in signatures")
    if len(set(idx)) != len(idx):
        raise ValueError("duplicate gene names in signatures")
    return idx


def panel_after_mt_filter(var_names: Sequence[str], symbols: Sequence[str] = None) -> List[str]:
    """Genes kept by run_benchmark.preprocess_spatial_data before the signature intersection."""
    symbols = list(var_names) if symbols is None else list(symbols)
    return [g for g, s in zip(var_names, symbols) if not str(s).startswith("MT-")]


def read_var(h5ad) -> Dict[str, list]:
    import anndata as ad

    a = ad.read_h5ad(h5ad, backed="r")
    try:
        names = [str(x) for x in a.var_names]
        sym = [str(x) for x in a.var["SYMBOL"]] if "SYMBOL" in a.var.columns else names
    finally:
        a.file.close()
    return {"names": names, "symbols": sym}


def read_annotated_vars(h5ad) -> List[str]:
    return read_var(h5ad)["names"]


def intersection_check(
    genes: Sequence[str], gex_h5ads: Dict[int, Path], old_annotated: Dict[int, Path]
) -> Dict[str, int]:
    """Per region: sorted(list ∩ panel-after-MT) must equal sorted(old annotated var_names). Raises otherwise."""
    if not gex_h5ads or set(gex_h5ads) != set(old_annotated):
        raise ValueError("region sets for GEX inputs and old annotated outputs differ or are empty")
    gset = set(genes)
    sizes = {}
    for r in sorted(gex_h5ads):
        v = read_var(gex_h5ads[r])
        kept = panel_after_mt_filter(v["names"], v["symbols"])
        got = sorted(gset.intersection(kept))
        want = sorted(read_annotated_vars(old_annotated[r]))
        if got != want:
            extra, miss = sorted(set(got) - set(want))[:5], sorted(set(want) - set(got))[:5]
            raise ValueError(
                f"region {r}: intersection {len(got)} != old inference {len(want)}; extra {extra} missing {miss}"
            )
        sizes[str(r)] = len(got)
    return sizes


def check_setting(setting_file=SETTING_FILE) -> str:
    d = json.loads(Path(setting_file).read_text())
    if d.get("value") != SETTING:
        raise RuntimeError(f"{setting_file} value {d.get('value')!r} != {SETTING!r}; this script builds only {SETTING}")
    return d["value"]


def build(
    sig_csv=SRC_SIG, out=OUT, gex_dir=GEX_DIR, old_run=OLD_RUN, setting_file=SETTING_FILE, prefix="Xenium"
) -> dict:
    out = Path(out)
    side = out.with_name(out.name + ".provenance.json")
    for p in (out, side):
        if p.exists():
            raise FileExistsError(f"refusing to overwrite {p}")
    check_setting(setting_file)
    genes = signature_genes(sig_csv)
    gex = {r: Path(gex_dir) / f"{prefix}_region_{r}_GEX.h5ad" for r in range(N_REGIONS)}
    old = {r: Path(old_run) / f"{prefix}_region_{r}" / f"{prefix}_region_{r}_annotated.h5ad" for r in range(N_REGIONS)}
    sizes = intersection_check(genes, gex, old)
    out.write_text("\n".join(genes) + "\n")
    lines = [ln for ln in out.read_text().splitlines() if ln]
    if lines != genes:
        raise RuntimeError("written list does not round-trip")
    prov = {
        "setting": SETTING,
        "n_genes": len(genes),
        "sha256": sha256_file(out),
        "source_path": str(Path(sig_csv).resolve()),
        "source_sha256": sha256_file(sig_csv),
        "source_n_rows": len(genes),
        "intersection_with_panel_after_mt_filter": sizes,
        "setting_file": str(Path(setting_file).resolve()),
        "setting_file_sha256": sha256_file(setting_file),
        "producer": "Benchmarking/xenium_benchmarking/gt_rebuild/build_fixed_gene_list.py",
    }
    side.write_text(json.dumps(prov, indent=2) + "\n")
    return prov


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
