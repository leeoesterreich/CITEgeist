# Wet-lab analyses for the midkine figures

Scripts and small input tables behind the wet-lab numbers in the manuscript. Run every command from this folder.
Outputs are written next to the scripts and are not shipped.

## IF quantification (Figure 6C, immunofluorescence)

- `analyze_hw120_if.py` quantifies midkine immunofluorescence in MCF7 and T47D, ESR1 wild type against D538G. Steps: max projection of each z-stack, rolling-ball background subtraction, Cellpose nuclear segmentation (Cellpose 4.0.9), nuclear, perinuclear and cytoplasmic compartments, per-cell intensity, two-way ANOVA with Tukey HSD, figures.
- Input: `--data-root`, the folder of raw z-stack TIFFs with one sub-folder per condition. The raw images are not part of this repository.
- Output: `--output-dir` (default `if_analysis/`), with the per-cell table, the statistics and the figures.
- Run: `IF_DATA_ROOT=/path/to/raw_zstacks gpu-race-submit.sh sbatch_hw120_if.sh`, or call the script directly with `--data-root` (`--dry-run` only lists the images). Environment: Python with cellpose 4.0.9, scikit-image, tifffile, seaborn, scipy and pandas.

## ELISA (Figure 6D)

- `elisa_table.py` fits a four-parameter logistic standard curve per plate, subtracts the blanks, divides by protein per microlitre and compares arms with Welch tests.
- Input: `elisa_plate2_wells.csv`, the raw 450 nm readings of plate 2.
- Run: `python elisa_table.py`.

## ddPCR

- `ddpcr_table.py` turns the droplet counts into the copy numbers and thresholds used in the text (P4-S2 DNA and cDNA, control wells).
- Input: `ddpcr_wells.csv`, a de-identified table copied from the instrument export.
- Run: `python ddpcr_table.py`.

## ER ChIP-seq (Figure 6A)

- `build_er_chipseq_table.py` builds the peak-score table (largest peak score within 10 kb of each transcription start site, capped at 3100). `build_er_chipseq_mdk_distance.py` then records the distance from the nearest ER peak to the MDK start site.
- Inputs: the eight GSE125117 peak BED files, downloaded from GEO into `er_chipseq/beds/` (the file names are listed in `build_er_chipseq_table.py`); `er_chipseq/er_chipseq_tss_hg38.tsv`, the shipped transcription start sites (GENCODE v49 MANE Select); and the GENCODE v49 annotation GTF for the re-derivation check.
- Run: `GTF=/path/to/gencode.v49.annotation.gtf.gz sbatch sbatch_er_chipseq_mdk_distance.sh`, or run the two scripts in that order with `--gtf`. The distance script stops without writing if a BED checksum, the start-site table or the window overlap disagrees with the table script.
