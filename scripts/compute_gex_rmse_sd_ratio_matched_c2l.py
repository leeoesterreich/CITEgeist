#!/usr/bin/env python
"""RMSE/sigma on Cell2Location-matched genes (and spots) -- rnaprot / rnaprot_v41 labelsets only.

Same definition as scripts/compute_gex_rmse_sd_ratio.py (it is imported, not copied): per (region, type) block,
per gene g: rmse_g / sd_g on log1p expression (sd_g = sample std (ddof=1) of GT log1p over the block's spots, genes with sd 0
dropped, genes with no GT counts dropped), block_ratio = mean over genes, method scalar = mean over scored blocks.

Variant 1 "matched_genes": genes restricted to C2L's gene set (union over all C2L layer CSVs, upper-cased --
matched_gene_r.c2l_union_genes, identical to the matched-gene r). Roster: CITEgeist, CITEgeist_v41p, C2L, Tangram
(rnaprot_v41 adds CITEgeist_p138, the paper's v41 CITEgeist, and writes to evaluation/rnaprot_v41/).
Variant 2 "matched_genes_spots": additionally, for each pair (CITEgeist vs C2L, CITEgeist_v41p vs C2L), per block the
spots AND genes are the intersection of GT, the method's layer and C2L's layer, so both sides of a pair score the
identical genes x spots. Rows are named "<method>|vs_C2L" and "C2L|vs_<method>".

WRITES (new files only, stamped via scripts/provenance.py; never overwrites):
  evaluation/rnaprot/rmse_sd_ratio_matched_c2l.json, rmse_sd_ratio_matched_c2l_table.csv
"""

from __future__ import annotations

import logging
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import compute_gex_rmse_sd_ratio as base  # noqa: E402  (same module object the tests patch)

PRODUCER = "scripts/compute_gex_rmse_sd_ratio_matched_c2l.py"
OUT_JSON = "rmse_sd_ratio_matched_c2l.json"
OUT_CSV = "rmse_sd_ratio_matched_c2l_table.csv"
LABELSETS = ("rnaprot", "rnaprot_v41")
PAIRED = ("CITEgeist", "CITEgeist_v41p", "CITEgeist_p138", "CITEgeist_v2")
SPOT_RULE = (
    "per (region, type) block the scored spots are the intersection of the GT block, the method's layer and C2L's "
    "layer (genes likewise: GT x method x C2L-layer x C2L gene set); both members of a pair use identical genes x spots"
)
ALL_SPOTS_RULE = (
    "fair standard (RMSE_SIGMA_ALL_SPOTS): every block is scored on ALL GT spots for every method, a spot absent from a layer is an "
    "all-zero row; genes are GT x method layer x C2L gene set (and x C2L layer for the paired variant)"
)
logger = logging.getLogger(__name__)


def summarize(records: list[dict], aggregate: dict) -> dict:
    """Per-method scalar plus the blocks / genes / spots actually used."""
    ng = [r["n_genes_kept"] for r in records]
    return {
        "method_scalar": aggregate["method_scalar"],
        "n_blocks": aggregate["n_blocks"],
        "n_blocks_excluded": aggregate["n_blocks_excluded"],
        "n_blocks_gt_all_zero": aggregate["n_blocks_gt_all_zero"],
        "excluded": aggregate["excluded"],
        "n_spots_total": int(sum(r["n_spots"] for r in records)),
        "n_genes_kept_total": int(sum(ng)),
        "n_genes_kept_mean_per_block": float(np.mean(ng)) if ng else float("nan"),
        "n_genes_kept_min_per_block": int(min(ng)) if ng else 0,
        "n_genes_kept_max_per_block": int(max(ng)) if ng else 0,
    }


def per_type_means(records: list[dict]) -> dict:
    """cell_type -> mean over regions of block_ratio (scored blocks, NaN ignored)."""
    out = {}
    for t in base.CELL_TYPES:
        vals = [r["block_ratio"] for r in records if r["cell_type"] == t and not math.isnan(r["block_ratio"])]
        out[t] = {"mean_over_regions": float(np.mean(vals)) if vals else None, "n_regions": len(vals)}
    return out


def run_variant(roster: dict, gene_set: set, c2l_dir: Path, spot_match: bool) -> dict:
    """roster: {name: layer dir}. spot_match=False -> gene restriction only; True -> pairwise vs C2L (PAIRED only)."""
    jobs = {}
    if not spot_match:
        for name, d in roster.items():
            jobs[name] = (name, d, None)
    else:
        for name in PAIRED:
            if name in roster:
                jobs[f"{name}|vs_C2L"] = (name, roster[name], c2l_dir)
                jobs[f"C2L|vs_{name}"] = ("C2L", c2l_dir, roster[name])
    methods = {}
    for key, (name, d, partner) in jobs.items():
        if not Path(d).exists():
            raise base.MissingLayerError(f"expected method dir missing: method={name} dir={d}")
        records, agg = base.evaluate_method(name, Path(d), gene_set=gene_set, partner_dir=partner, all_spots=base.ALL_SPOTS)
        methods[key] = {
            "summary": summarize(records, agg),
            "per_type": per_type_means(records),
            "records": records,
            "method_dir": _rel(d),
        }
        logger.info("%s scalar=%.3f blocks=%d", key, agg["method_scalar"], agg["n_blocks"])
    return methods


def table(variants: dict) -> pd.DataFrame:
    rows = []
    for vname, methods in variants.items():
        df = base.blocks_table({k: {"records": v["records"], "aggregate": v["summary"]} for k, v in methods.items()})
        # blocks_table keys excluded rows off aggregate["excluded"]; add gene/spot counts for scored rows
        n_genes = {(k, r["region"], r["cell_type"]): r["n_genes_kept"] for k, v in methods.items() for r in v["records"]}
        df.insert(0, "variant", vname)
        df["n_genes_kept"] = [n_genes.get((m, r, c)) for m, r, c in zip(df["method"], df["region"], df["cell_type"])]
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def main() -> int:
    if base._LS not in LABELSETS:
        raise SystemExit(f"XEN_LABELSET must be one of {LABELSETS} (got {base._LS!r})")
    base.manifest.require_rnagt_signoff(labelset=base._LS)  # before reading or writing anything
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    sys.path.insert(0, str(base.BENCH_ROOT / "spots_baseline"))
    from matched_gene_r import c2l_union_genes  # noqa: E402  (the one C2L gene-set definition)
    from provenance import stamp_payload, write_csv_sidecar

    out_dir = base.output_dir(base._LS)
    json_path, csv_path = out_dir / OUT_JSON, out_dir / OUT_CSV
    for p in (json_path, csv_path):
        if p.exists():
            raise FileExistsError(f"refusing to overwrite {p}")
    roster = {k: v for k, v in base.ALL_METHODS.items() if k in PAIRED + ("C2L", "Tangram")}
    c2l_dir = Path(roster["C2L"])
    genes = {g.upper() for g in c2l_union_genes(c2l_dir)}
    assert genes, f"no C2L layer CSVs under {c2l_dir}"
    variants = {
        "matched_genes": run_variant(roster, genes, c2l_dir, spot_match=False),
        "matched_genes_spots": run_variant(roster, genes, c2l_dir, spot_match=True),
    }
    inputs = {}
    for ms in variants.values():
        for v in ms.values():
            for r in v["records"]:
                inputs[r["gt_path"]] = r["gt_sha256"]
                inputs[r["pred_path"]] = r["pred_sha256"]
    payload = {
        "_generated_at": datetime.now(timezone.utc).isoformat(),
        "labelset": base._LS,
        "definition": base_definition(),
        "gene_set": "union of upper-cased gene names over all Cell2Location layer CSVs (matched_gene_r.c2l_union_genes)",
        "n_c2l_genes": len(genes),
        "spot_match_rule": SPOT_RULE if not base.ALL_SPOTS else ALL_SPOTS_RULE,
        "spot_match_scope": f"only {', '.join(m for m in PAIRED if m in roster)} are paired with C2L; Tangram appears in matched_genes only",
        "regions": base.REGIONS,
        "cell_types": base.CELL_TYPES,
        "gt_gex_dir": str(base.GT_ROOT.relative_to(base.REPO_ROOT)),
        "variants": {
            vn: {
                "method_scalars": {k: _r3(v["summary"]["method_scalar"]) for k, v in ms.items()},
                "methods": {k: {kk: vv for kk, vv in v.items() if kk != "records"} for k, v in ms.items()},
            }
            for vn, ms in variants.items()
        },
        "input_sha256": inputs,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    base._atomic_write_json(json_path, stamp_payload(payload, PRODUCER))
    table(variants).to_csv(csv_path, index=False)
    write_csv_sidecar(csv_path, PRODUCER)
    logger.info("wrote %s and %s", json_path, csv_path)
    return 0


def _rel(d) -> str:
    try:
        return str(Path(d).relative_to(base.REPO_ROOT))
    except ValueError:
        return str(d)


def _r3(x):
    return round(x, 3) if not math.isnan(x) else None


def base_definition() -> str:
    return (
        "per-gene RMSE / per-gene GT std (log1p), mean across genes with sd>0 per (region, type) block, "
        "then mean across scored blocks; unchanged from compute_gex_rmse_sd_ratio.py, with genes (and for "
        "matched_genes_spots also spots) restricted as described"
    )


if __name__ == "__main__":
    sys.exit(main())
