#!/usr/bin/env python3
"""Distance of the nearest ER ChIP-seq peak to the MDK TSS, and the secretory-gene peak counts, from the GSE125117 peak BEDs.

Companion of build_er_chipseq_table.py (same BEDs, same TSS table, same window). Feeds the Discussion statements
"ER peaks within 10 kb of the TSS at two of eight secretory-pathway genes", "no peak within 10 kb of the MDK TSS" and
"in all four E2-treated conditions ER peaks lay about 18 kb upstream".

Checks (the script fails, never writes, if one breaks):
  * the md5 of every BED equals the one recorded in er_chipseq_inputs.json;
  * the TSS table equals the MANE Select TSS re-derived from the GENCODE v49 GTF (--gtf), for all genes of the table;
  * for every gene and condition, "a peak overlaps TSS +/- 10 kb" agrees with er_chipseq_peak_scores.csv (> 0).

Distances: TSS is 0-based (table tss - 1). For a peak [st, en), edge distance = bp from the TSS to the nearest edge of the peak
(0 if the peak covers the TSS); midpoint distance = (st + en) / 2 - TSS. Signed so that NEGATIVE = upstream (against the
direction of transcription of the gene).

Output (stamped with scripts/provenance.py): er_chipseq/er_chipseq_mdk_distance.json
Usage:  python -u build_er_chipseq_mdk_distance.py --gtf GENCODE_v49.gtf.gz [--bed-dir DIR] [--out PATH]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[2]  # repo root
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(PROJECT / "scripts"))
import build_er_chipseq_table as T  # noqa: E402
import provenance  # noqa: E402

OUT = T.OUT_DIR / "er_chipseq_mdk_distance.json"
PRODUCER = "Benchmarking/midkine_remediation/wetlab/build_er_chipseq_mdk_distance.py"
TARGET = "MDK"
SECRETORY = [g for g in T.FIG6A_GENES if g != TARGET]


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def overlaps_window(peaks: dict, chrom: str, tss0: int) -> list:
    lo, hi = tss0 - T.WINDOW, tss0 + T.WINDOW + 1  # same window as build_er_chipseq_table.max_score
    return [(st, en, sc) for st, en, sc in peaks.get(chrom, []) if st < hi and en > lo]


def nearest(peaks: dict, chrom: str, tss0: int, strand: str) -> dict | None:
    best = None
    for st, en, sc in peaks.get(chrom, []):
        if st <= tss0 < en:
            edge = 0
        elif en <= tss0:
            edge = tss0 - en  # peak lies at lower coordinates
        else:
            edge = st - tss0
        if best is None or edge < best[0]:
            best = (edge, st, en, sc)
    if best is None:
        return None
    edge, st, en, sc = best
    mid = (st + en) / 2 - tss0
    sign = 1 if strand == "+" else -1  # upstream = lower coordinate on the + strand, higher on the - strand
    side = 0 if edge == 0 else (1 if (st > tss0) else -1)
    return {
        "chrom": chrom,
        "start": st,
        "end": en,
        "score": sc,
        "edge_distance_bp": sign * side * edge,
        "midpoint_distance_bp": sign * mid,
        "side": "covers_tss" if edge == 0 else ("downstream" if sign * side > 0 else "upstream"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed-dir", type=Path, default=T.DEFAULT_BED_DIR)
    ap.add_argument("--gtf", type=Path, required=True, help="GENCODE v49 GTF(.gz) to verify the TSS table against")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    recorded = json.loads(T.INPUTS_JSON.read_text())["peak_files_md5"]
    peaks, md5s = {}, {}
    for cond, fname in T.BED_FILES.items():
        path = args.bed_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"Missing peak file: {path}")
        md5s[fname] = md5(path)
        if md5s[fname] != recorded.get(fname):
            raise ValueError(f"{fname}: md5 {md5s[fname]} != {recorded.get(fname)} recorded in er_chipseq_inputs.json")
        peaks[cond] = T.load_bed(path)

    tss = pd.read_csv(T.TSS_TSV, sep="\t").set_index("gene")
    derived = T.derive_tss_table(args.gtf, list(tss.index)).set_index("gene")
    for g in tss.index:
        a, b = tss.loc[g], derived.loc[g]
        if (a["chrom"], int(a["tss"]), a["strand"], a["transcript"]) != (
            b["chrom"],
            int(b["tss"]),
            b["strand"],
            b["transcript"],
        ):
            raise ValueError(f"{g}: TSS table {tuple(a)} != GTF-derived {tuple(b)}")

    scores = pd.read_csv(T.TABLE_CSV).set_index("gene")
    per_gene, with_peak = {}, []
    for g in SECRETORY + [TARGET]:
        chrom, tss0 = tss.loc[g, "chrom"], int(tss.loc[g, "tss"]) - 1
        hits = {c: overlaps_window(peaks[c], chrom, tss0) for c in T.BED_FILES}
        for c in T.BED_FILES:
            if bool(hits[c]) != bool(scores.loc[g, c] > 0):
                raise ValueError(
                    f"{g} {c}: window overlap {bool(hits[c])} disagrees with er_chipseq_peak_scores.csv ({scores.loc[g, c]})"
                )
        conds = [c for c in T.BED_FILES if hits[c]]
        per_gene[g] = {"conditions_with_peak_within_10kb": conds, "n_conditions": len(conds)}
        if g != TARGET and conds:
            with_peak.append(g)

    mdk = tss.loc[TARGET]
    mdk_tss0 = int(mdk["tss"]) - 1
    per_cond = {}
    for c in T.BED_FILES:
        near = nearest(peaks[c], mdk["chrom"], mdk_tss0, mdk["strand"])
        per_cond[c] = {"n_peaks_in_file": sum(len(v) for v in peaks[c].values()), "nearest_peak": near}
    e2 = [c for c in T.BED_FILES if c.endswith("_E2")]
    e2_edges = [per_cond[c]["nearest_peak"]["edge_distance_bp"] for c in e2]
    payload = {
        "genome_build": T.GENOME_BUILD,
        "window_bp_each_side_of_tss": T.WINDOW,
        "peak_files_md5": md5s,
        "tss_verified_against": {
            "file": args.gtf.name,
            "md5": md5(args.gtf),
            "genes": list(tss.index),
            "rule": "MANE Select transcript start, strand-aware",
        },
        "mdk_tss": {
            "chrom": mdk["chrom"],
            "tss": int(mdk["tss"]),
            "strand": mdk["strand"],
            "transcript": mdk["transcript"],
        },
        "secretory_genes": SECRETORY,
        "n_secretory_genes": len(SECRETORY),
        "secretory_genes_with_peak_within_10kb": with_peak,
        "n_secretory_genes_with_peak_within_10kb": len(with_peak),
        "n_secretory_genes_without_peak_within_10kb": len(SECRETORY) - len(with_peak),
        "per_gene": per_gene,
        "mdk_n_conditions_with_peak_within_10kb": per_gene[TARGET]["n_conditions"],
        "mdk_nearest_peak_by_condition": per_cond,
        "mdk_e2_conditions": e2,
        "mdk_e2_nearest_edge_distance_bp": dict(zip(e2, e2_edges)),
        "mdk_e2_all_upstream": all(per_cond[c]["nearest_peak"]["side"] == "upstream" for c in e2),
        "mdk_e2_nearest_edge_min_abs_bp": min(abs(x) for x in e2_edges),
        "mdk_e2_nearest_edge_max_abs_bp": max(abs(x) for x in e2_edges),
    }
    provenance.write_stamped_json(args.out, payload, PRODUCER)
    print(f"Wrote {args.out}")
    print(
        json.dumps(
            {k: payload[k] for k in payload if k.startswith(("n_secretory", "secretory_genes_with", "mdk_"))}, indent=1
        )[:3000]
    )


if __name__ == "__main__":
    main()
