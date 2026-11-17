#!/bin/bash
# ER ChIP-seq peak table and the MDK peak-distance record (Figure 6A and the Discussion statements about ER peaks near MDK).
# Usage: GTF=/path/to/gencode.v49.annotation.gtf.gz sbatch sbatch_er_chipseq_mdk_distance.sh
# First download the eight GSE125117 peak BEDs into er_chipseq/beds (names in build_er_chipseq_table.py).
#SBATCH --job-name=er_chipseq_mdk
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=er_chipseq_mdk_%j.out
#SBATCH --error=er_chipseq_mdk_%j.err

set -euo pipefail

: "${GTF:?set GTF to the GENCODE v49 annotation GTF (gz)}"
HERE="${HERE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"  # submit from this folder, or set HERE_DIR
python -u "$HERE/build_er_chipseq_table.py" --gtf "$GTF"
python -u "$HERE/build_er_chipseq_mdk_distance.py" --gtf "$GTF"
