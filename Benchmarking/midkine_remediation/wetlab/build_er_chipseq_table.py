#!/usr/bin/env python3
"""Build the ER ChIP-seq peak-score table behind Figure 6A (draw_panel_H).

Method (identical to Supplementary Figure S9B, generate_supp_mdk_intensity_chipseq.py):
for each gene and condition, the value is the largest MACS peak score (BED
column 5) among peaks overlapping the window TSS +/- 10 kb, capped at 3100.
0 means no peak overlapped the window. No other score threshold is applied.

Genome build: hg38 (GRCh38). The GSE125117 peak BEDs carry GRCh38-only contig
names (e.g. chr17_GL000205v2_random, chr22_KI270733v1_random), and ER peaks sit
on the canonical E2 targets only when TSS coordinates are hg38 (TFF1 scores 3100
in MCF7 WT +E2, 0 in vehicle). TSS = start of the GENCODE v49 MANE Select
transcript (strand-aware), stored in er_chipseq_tss_hg38.tsv. Pass --gtf to
re-derive that table from a GENCODE v49 GTF.

Inputs : 8 GSE125117 peak BEDs (GSM3563750/51/56/57/59/60/65/66), looked up in
         --bed-dir (default er_chipseq/beds next to this script; download the BEDs from GEO GSE125117). MD5s are recorded in
         er_chipseq_inputs.json.
Outputs: er_chipseq_peak_scores.csv   (Figure 6A genes x 8 conditions)
         er_chipseq_inputs.json       (BED md5s, window, cap, genome build)

Usage:
    python build_er_chipseq_table.py [--bed-dir DIR] [--gtf GTF.gz]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parents[2]  # repo root
OUT_DIR = SCRIPT_DIR / "er_chipseq"
TSS_TSV = OUT_DIR / "er_chipseq_tss_hg38.tsv"
TABLE_CSV = OUT_DIR / "er_chipseq_peak_scores.csv"
INPUTS_JSON = OUT_DIR / "er_chipseq_inputs.json"
DEFAULT_BED_DIR = OUT_DIR / "beds"

WINDOW = 10000  # TSS +/- 10 kb
SCORE_CAP = 3100.0  # GEO MACS score ceiling
GENOME_BUILD = "hg38"

# Figure 6A genes (row order is re-sorted by draw_panel_H) and S9B-only genes.
FIG6A_GENES = ["HSP90B1", "HSPA5", "XBP1", "ATF6", "ATF4", "CALR", "CANX", "PDIA4", "MDK"]
EXTRA_GENES = ["TFF1", "GREB1", "PGR", "SEC61B", "ERN1"]  # used by Supplementary S9B

# Column order matches the table draw_panel_H has always read.
BED_FILES = {
    "MCF7_WT_Veh": "GSM3563750_MCF7_WT_Vehicle_peaks.bed.gz",
    "MCF7_WT_E2": "GSM3563751_MCF7_WT_E2_peaks.bed.gz",
    "MCF7_D538G_Veh": "GSM3563756_MCF7_D538G_Vehicle_peaks.bed.gz",
    "MCF7_D538G_E2": "GSM3563757_MCF7_D538G_E2_peaks.bed.gz",
    "T47D_WT_Veh": "GSM3563759_T47D_WT_Vehicle_peaks.bed.gz",
    "T47D_WT_E2": "GSM3563760_T47D_WT_E2_peaks.bed.gz",
    "T47D_D538G_Veh": "GSM3563765_T47D_D538G_Vehicle_peaks.bed.gz",
    "T47D_D538G_E2": "GSM3563766_T47D_D538G_E2_peaks.bed.gz",
}


def derive_tss_table(gtf_gz: Path, genes: list[str]) -> pd.DataFrame:
    """TSS of the MANE Select transcript of each gene (GENCODE GTF, 1-based)."""
    rows: dict[str, dict] = {}
    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] != "transcript" or 'tag "MANE_Select"' not in p[8]:
                continue
            gene = re.search(r'gene_name "([^"]+)"', p[8]).group(1)
            if gene not in genes:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
            tss = int(p[3]) if p[6] == "+" else int(p[4])
            rows[gene] = {"gene": gene, "chrom": p[0], "tss": tss, "strand": p[6], "transcript": tid}
    missing = [g for g in genes if g not in rows]
    if missing:
        raise ValueError(f"no MANE Select transcript found for: {missing}")
    return pd.DataFrame([rows[g] for g in genes])


def load_bed(path: Path) -> dict[str, list[tuple[int, int, float]]]:
    by_chrom: dict[str, list[tuple[int, int, float]]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            by_chrom.setdefault(p[0], []).append((int(p[1]), int(p[2]), float(p[4])))
    return by_chrom


def max_score(peaks: dict, chrom: str, start: int, end: int) -> float:
    """Largest peak score among BED peaks (0-based, half-open) overlapping [start, end)."""
    best = 0.0
    for st, en, sc in peaks.get(chrom, []):
        if st < end and en > start and sc > best:
            best = sc
    return best


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed-dir", type=Path, default=DEFAULT_BED_DIR)
    ap.add_argument("--gtf", type=Path, default=None, help="GENCODE v49 GTF(.gz): re-derive the TSS table")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.gtf is not None:
        tss = derive_tss_table(args.gtf, FIG6A_GENES + EXTRA_GENES)
        tss.to_csv(TSS_TSV, sep="\t", index=False)
        print(f"Wrote {TSS_TSV}")
    tss = pd.read_csv(TSS_TSV, sep="\t").set_index("gene")

    peaks, md5 = {}, {}
    for cond, fname in BED_FILES.items():
        path = args.bed_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"Missing peak file: {path}")
        peaks[cond] = load_bed(path)
        md5[fname] = hashlib.md5(path.read_bytes()).hexdigest()

    rows = []
    for gene in FIG6A_GENES:
        chrom, tss0 = tss.loc[gene, "chrom"], int(tss.loc[gene, "tss"]) - 1  # 1-based -> 0-based
        row = {"gene": gene}
        for cond in BED_FILES:
            row[cond] = min(max_score(peaks[cond], chrom, tss0 - WINDOW, tss0 + WINDOW + 1), SCORE_CAP)
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(TABLE_CSV, index=False, float_format="%g")
    INPUTS_JSON.write_text(
        json.dumps(
            {
                "genome_build": GENOME_BUILD,
                "tss_source": "GENCODE v49 MANE Select transcript start (strand-aware)",
                "window_bp_each_side_of_tss": WINDOW,
                "score_cap": SCORE_CAP,
                "peak_files_md5": md5,
            },
            indent=2,
        )
        + "\n"
    )
    print(table.to_string(index=False))
    print(f"Wrote {TABLE_CSV}")


if __name__ == "__main__":
    main()
