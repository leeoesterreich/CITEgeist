#!/usr/bin/env python
"""
13_discriminate_hsp90b1_mechanism.py

HSP90B1 MECHANISM DISCRIMINATION: FOLDING vs LRP1 TRAPPING

Purpose:
Given that HSP90B1 is differentially regulated between MCF7 and T47D
(established by scripts 02-10), determine HOW HSP90B1 influences MDK
secretion through hypothesis-driven gene signature analysis.

Two mechanistic models:
  Model A (Direct folding/export): HSP90B1 as ER chaperone improves MDK
    folding capacity, reducing ER retention and ERAD.
  Model B (LRP1 premature binding trap): HSP90B1 matures LRP1, which can
    trap MDK intracellularly during biosynthesis.

Method: 8 diagnostic tests using existing pipeline datasets, producing
an evidence scorecard with effect sizes and p-values.

Datasets:
  GSE89888  - RNA-seq, MCF7/T47D WT vs D538G (primary)
  GSE75329  - RNA-seq, MCF7 perturbations (validation)
  GSE125117 - ER ChIP-seq (binding at LRP1)
  GSE254216 - ATAC-seq (accessibility at LRP1)

Outputs:
  outputs/tables/mechanism_discrimination_scorecard.csv
  outputs/tables/mechanism_gene_signatures.csv
  outputs/figures/fig13_mechanism_discrimination.png
"""

import sys
import gzip
from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd
import numpy as np
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).parent
PIPELINE_DIR = SCRIPT_DIR.parent
CONFIG_PATH = PIPELINE_DIR / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

DATA_DIR = PIPELINE_DIR / config['paths']['data_dir']
OUTPUT_DIR = PIPELINE_DIR / config['paths']['output_dir']
GENE_COORDS_HG38 = {k: tuple(v) for k, v in config['gene_coordinates_hg38'].items()}
SIGNATURES = config['mechanism_discrimination']

# Sample groups for GSE89888
GROUPS = {
    'MCF7_WT': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
    'MCF7_D538G': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
    'T47D_WT': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
    'T47D_D538G': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
}

# ELISA MDK secretion fold changes (D538G/WT)
ELISA_MDK = {'MCF7': 1.83, 'T47D': 0.38}


@dataclass
class TestResult:
    test_id: int
    test_name: str
    model_tested: str          # "A", "B", or "both"
    verdict: str               # "Supports A", "Supports B", "Ambiguous"
    effect_size: float         # Cohen's d or Pearson r
    p_value: float
    details: str               # One-sentence interpretation
    gene_level_data: pd.DataFrame = field(default_factory=pd.DataFrame)


def load_tpm():
    """Load GSE89888 TPM data with gene symbol mapping."""
    filepath = DATA_DIR / config['datasets']['GSE89888']
    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(df.index.astype(str).tolist(), scopes='entrezgene',
                              fields='symbol', species='human', returnall=True)
        id_to_sym = {str(h['query']): h['symbol'] for h in results['out'] if 'symbol' in h}
        df.index = df.index.astype(str).map(lambda x: id_to_sym.get(x, x))
    except Exception as e:
        print(f"Warning: Gene mapping failed: {e}")
    return df


def load_peaks(filepath):
    """Load BED peak file with signal score."""
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    peak = {
                        'chr': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'name': parts[3] if len(parts) > 3 else '',
                        'score': float(parts[4]) if len(parts) > 4 else 0.0,
                    }
                    peaks.append(peak)
                except ValueError:
                    continue
    return pd.DataFrame(peaks)


def count_peaks_in_region(peaks_df, chrom, start, end):
    """Count peaks overlapping a genomic region."""
    if peaks_df is None or len(peaks_df) == 0:
        return 0
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    overlapping = chr_peaks[(chr_peaks['start'] < end) & (chr_peaks['end'] > start)]
    return len(overlapping)


def get_peak_signal_in_region(peaks_df, chrom, start, end):
    """Get max peak signal in a genomic region."""
    if peaks_df is None or len(peaks_df) == 0:
        return 0.0
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    overlapping = chr_peaks[(chr_peaks['start'] < end) & (chr_peaks['end'] > start)]
    if len(overlapping) == 0:
        return 0.0
    return overlapping['score'].max()


def calculate_cohens_d(group1, group2):
    """Calculate Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0
    return (np.mean(group1) - np.mean(group2)) / pooled_std


def get_gene_fc_and_stats(tpm, gene):
    """Get fold changes and stats for a gene across cell lines."""
    if gene not in tpm.index:
        return None

    mcf7_wt = tpm.loc[gene, GROUPS['MCF7_WT']].values.astype(float)
    mcf7_d538g = tpm.loc[gene, GROUPS['MCF7_D538G']].values.astype(float)
    t47d_wt = tpm.loc[gene, GROUPS['T47D_WT']].values.astype(float)
    t47d_d538g = tpm.loc[gene, GROUPS['T47D_D538G']].values.astype(float)

    mcf7_wt_mean = np.mean(mcf7_wt)
    t47d_wt_mean = np.mean(t47d_wt)

    if mcf7_wt_mean <= 0 or t47d_wt_mean <= 0:
        return None

    mcf7_fc = np.mean(mcf7_d538g) / mcf7_wt_mean
    t47d_fc = np.mean(t47d_d538g) / t47d_wt_mean

    _, mcf7_pval = stats.ttest_ind(mcf7_d538g, mcf7_wt, equal_var=False)
    _, t47d_pval = stats.ttest_ind(t47d_d538g, t47d_wt, equal_var=False)

    return {
        'gene': gene,
        'MCF7_WT_TPM': mcf7_wt_mean,
        'MCF7_D538G_TPM': np.mean(mcf7_d538g),
        'MCF7_FC': mcf7_fc,
        'MCF7_log2FC': np.log2(mcf7_fc),
        'MCF7_pval': float(mcf7_pval),
        'MCF7_cohens_d': calculate_cohens_d(mcf7_d538g, mcf7_wt),
        'T47D_WT_TPM': t47d_wt_mean,
        'T47D_D538G_TPM': np.mean(t47d_d538g),
        'T47D_FC': t47d_fc,
        'T47D_log2FC': np.log2(t47d_fc),
        'T47D_pval': float(t47d_pval),
        'T47D_cohens_d': calculate_cohens_d(t47d_d538g, t47d_wt),
    }


# ============================================================================
# DIAGNOSTIC TESTS
# ============================================================================

def run_test_1_erad(gene_data, hsp90b1_stats):
    """Test 1: When HSP90B1 UP (MCF7-D538G), do ERAD genes go DOWN?"""
    print("\n" + "-" * 60)
    print("TEST 1: ERAD Pathway Suppression")
    print("-" * 60)

    erad_genes = SIGNATURES['erad_genes']
    mcf7_fcs = []
    gene_rows = []

    for gene in erad_genes:
        if gene in gene_data:
            d = gene_data[gene]
            mcf7_fcs.append(d['MCF7_log2FC'])
            gene_rows.append(d)
            direction_mcf7 = "UP" if d['MCF7_FC'] > 1.1 else "DOWN" if d['MCF7_FC'] < 0.9 else "NC"
            direction_t47d = "UP" if d['T47D_FC'] > 1.1 else "DOWN" if d['T47D_FC'] < 0.9 else "NC"
            print(f"  {gene}: MCF7 FC={d['MCF7_FC']:.2f} ({direction_mcf7}), T47D FC={d['T47D_FC']:.2f} ({direction_t47d})")

    if len(mcf7_fcs) < 3:
        return TestResult(1, "ERAD Pathway Suppression", "A", "Ambiguous", 0.0, 1.0,
                         f"Only {len(mcf7_fcs)} ERAD genes found", pd.DataFrame(gene_rows))

    mcf7_fcs = np.array(mcf7_fcs)
    t_stat, p_val = stats.ttest_1samp(mcf7_fcs, 0)
    mean_fc = np.mean(mcf7_fcs)
    effect = mean_fc / np.std(mcf7_fcs) if np.std(mcf7_fcs) > 0 else 0

    if p_val < 0.05 and mean_fc < 0:
        verdict = "Supports A"
        details = f"ERAD genes decrease in MCF7-D538G (mean log2FC={mean_fc:.3f}, p={p_val:.4f}), consistent with reduced misfolding"
    elif p_val < 0.05 and mean_fc > 0:
        verdict = "Supports B"
        details = f"ERAD genes increase in MCF7-D538G (mean log2FC={mean_fc:.3f}, p={p_val:.4f}), inconsistent with improved folding"
    else:
        verdict = "Ambiguous"
        details = f"ERAD genes show no significant change (mean log2FC={mean_fc:.3f}, p={p_val:.4f})"

    print(f"\n  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(1, "ERAD Pathway Suppression", "A", verdict, effect, p_val, details, pd.DataFrame(gene_rows))


def run_test_2_upr(gene_data, hsp90b1_stats):
    """Test 2: When HSP90B1 UP (MCF7-D538G), do UPR markers go DOWN?"""
    print("\n" + "-" * 60)
    print("TEST 2: UPR Stress Relief")
    print("-" * 60)

    upr_genes = SIGNATURES['upr_genes']
    mcf7_fcs = []
    gene_rows = []

    for gene in upr_genes:
        if gene in gene_data:
            d = gene_data[gene]
            mcf7_fcs.append(d['MCF7_log2FC'])
            gene_rows.append(d)
            direction = "UP" if d['MCF7_FC'] > 1.1 else "DOWN" if d['MCF7_FC'] < 0.9 else "NC"
            print(f"  {gene}: MCF7 FC={d['MCF7_FC']:.2f} ({direction})")

    if len(mcf7_fcs) < 3:
        return TestResult(2, "UPR Stress Relief", "A", "Ambiguous", 0.0, 1.0,
                         f"Only {len(mcf7_fcs)} UPR genes found", pd.DataFrame(gene_rows))

    mcf7_fcs = np.array(mcf7_fcs)
    t_stat, p_val = stats.ttest_1samp(mcf7_fcs, 0)
    mean_fc = np.mean(mcf7_fcs)
    effect = mean_fc / np.std(mcf7_fcs) if np.std(mcf7_fcs) > 0 else 0
    n_down = int(np.sum(mcf7_fcs < -np.log2(1.1)))
    n_up = int(np.sum(mcf7_fcs > np.log2(1.1)))
    n_nc = len(mcf7_fcs) - n_down - n_up

    if p_val < 0.05 and mean_fc < 0:
        verdict = "Supports A"
        details = f"UPR markers decrease in MCF7-D538G ({n_down} DOWN, {n_up} UP, {n_nc} NC; p={p_val:.4f})"
    elif p_val < 0.05 and mean_fc > 0:
        verdict = "Supports B"
        details = f"UPR markers increase in MCF7-D538G ({n_down} DOWN, {n_up} UP, {n_nc} NC; p={p_val:.4f})"
    else:
        verdict = "Ambiguous"
        details = f"UPR markers show no clear trend ({n_down} DOWN, {n_up} UP, {n_nc} NC; p={p_val:.4f})"

    print(f"\n  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(2, "UPR Stress Relief", "A", verdict, effect, p_val, details, pd.DataFrame(gene_rows))


def run_test_3_secretory(gene_data, hsp90b1_stats):
    """Test 3: Do other secretory chaperones co-regulate with HSP90B1?"""
    print("\n" + "-" * 60)
    print("TEST 3: Secretory Capacity Co-regulation")
    print("-" * 60)

    sec_genes = SIGNATURES['secretory_chaperones']
    hsp_mcf7_up = hsp90b1_stats['MCF7_FC'] > 1.0
    hsp_t47d_up = hsp90b1_stats['T47D_FC'] > 1.0

    mcf7_concordant = 0
    t47d_concordant = 0
    total = 0
    gene_rows = []

    for gene in sec_genes:
        if gene in gene_data:
            d = gene_data[gene]
            total += 1
            gene_rows.append(d)
            gene_mcf7_up = d['MCF7_FC'] > 1.0
            gene_t47d_up = d['T47D_FC'] > 1.0
            if gene_mcf7_up == hsp_mcf7_up:
                mcf7_concordant += 1
            if gene_t47d_up == hsp_t47d_up:
                t47d_concordant += 1
            print(f"  {gene}: MCF7 FC={d['MCF7_FC']:.2f} ({'same' if gene_mcf7_up == hsp_mcf7_up else 'opp'}), "
                  f"T47D FC={d['T47D_FC']:.2f} ({'same' if gene_t47d_up == hsp_t47d_up else 'opp'})")

    if total == 0:
        return TestResult(3, "Secretory Capacity Co-regulation", "A", "Ambiguous", 0.0, 1.0,
                         "No secretory chaperone genes found", pd.DataFrame(gene_rows))

    mcf7_frac = mcf7_concordant / total
    t47d_frac = t47d_concordant / total
    avg_frac = (mcf7_frac + t47d_frac) / 2

    concordant_total = mcf7_concordant + t47d_concordant
    p_val = float(stats.binomtest(concordant_total, 2 * total, 0.5, alternative='greater').pvalue)

    print(f"\n  MCF7 concordance: {mcf7_concordant}/{total} ({mcf7_frac:.1%})")
    print(f"  T47D concordance: {t47d_concordant}/{total} ({t47d_frac:.1%})")

    if avg_frac > 0.6 and p_val < 0.1:
        verdict = "Supports A"
        details = f"Secretory chaperones co-regulate with HSP90B1 ({avg_frac:.0%} concordant, p={p_val:.4f})"
    elif avg_frac < 0.4:
        verdict = "Supports B"
        details = f"Secretory chaperones do NOT co-regulate with HSP90B1 ({avg_frac:.0%} concordant)"
    else:
        verdict = "Ambiguous"
        details = f"Mixed co-regulation pattern ({avg_frac:.0%} concordant, p={p_val:.4f})"

    print(f"\n  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(3, "Secretory Capacity Co-regulation", "A", verdict, avg_frac, p_val, details, pd.DataFrame(gene_rows))


def run_test_4_clients(gene_data, hsp90b1_stats):
    """Test 4: Do known GRP94 clients show expression concordant with HSP90B1?"""
    print("\n" + "-" * 60)
    print("TEST 4: GRP94 Client Concordance")
    print("-" * 60)

    client_genes = SIGNATURES['grp94_clients']
    hsp_fcs = np.array([hsp90b1_stats['MCF7_log2FC'], hsp90b1_stats['T47D_log2FC']])

    client_fcs_mcf7 = []
    client_fcs_t47d = []
    gene_rows = []

    for gene in client_genes:
        if gene in gene_data:
            d = gene_data[gene]
            client_fcs_mcf7.append(d['MCF7_log2FC'])
            client_fcs_t47d.append(d['T47D_log2FC'])
            gene_rows.append(d)
            print(f"  {gene}: MCF7 FC={d['MCF7_FC']:.2f}, T47D FC={d['T47D_FC']:.2f}")

    if len(client_fcs_mcf7) < 2:
        return TestResult(4, "GRP94 Client Concordance", "A", "Ambiguous", 0.0, 1.0,
                         f"Only {len(client_fcs_mcf7)} GRP94 clients found", pd.DataFrame(gene_rows))

    concordant = 0
    for mcf7_fc, t47d_fc in zip(client_fcs_mcf7, client_fcs_t47d):
        client_pattern = np.array([mcf7_fc, t47d_fc])
        if np.corrcoef(hsp_fcs, client_pattern)[0, 1] > 0:
            concordant += 1

    total = len(client_fcs_mcf7)
    frac = concordant / total
    p_val = float(stats.binomtest(concordant, total, 0.5, alternative='greater').pvalue)

    print(f"\n  Concordant with HSP90B1 pattern: {concordant}/{total} ({frac:.0%})")

    if frac > 0.6 and p_val < 0.1:
        verdict = "Supports A"
        details = f"GRP94 clients follow HSP90B1 pattern ({concordant}/{total} concordant, p={p_val:.4f})"
    elif frac < 0.3:
        verdict = "Supports B"
        details = f"GRP94 clients do NOT follow HSP90B1 ({concordant}/{total} concordant)"
    else:
        verdict = "Ambiguous"
        details = f"Mixed client concordance ({concordant}/{total}, p={p_val:.4f})"

    print(f"\n  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(4, "GRP94 Client Concordance", "A", verdict, frac, p_val, details, pd.DataFrame(gene_rows))


def run_test_5_lrp1_expression(gene_data):
    """Test 5: Is LRP1 expressed and regulated in these cell lines?"""
    print("\n" + "-" * 60)
    print("TEST 5: LRP1 Expression and Regulation")
    print("-" * 60)

    if 'LRP1' not in gene_data:
        print("  LRP1 NOT FOUND in expression data")
        return TestResult(5, "LRP1 Expression", "B", "Supports A", 0.0, 1.0,
                         "LRP1 not detected in dataset, Model B implausible", pd.DataFrame())

    d = gene_data['LRP1']
    mcf7_wt_tpm = d['MCF7_WT_TPM']
    t47d_wt_tpm = d['T47D_WT_TPM']
    mcf7_fc = d['MCF7_FC']
    t47d_fc = d['T47D_FC']

    print(f"  LRP1 MCF7: WT={mcf7_wt_tpm:.1f} TPM, D538G FC={mcf7_fc:.2f}")
    print(f"  LRP1 T47D: WT={t47d_wt_tpm:.1f} TPM, D538G FC={t47d_fc:.2f}")

    expressed = mcf7_wt_tpm > 1 or t47d_wt_tpm > 1
    regulated = abs(d['MCF7_log2FC']) > np.log2(1.2) or abs(d['T47D_log2FC']) > np.log2(1.2)

    if not expressed:
        verdict = "Supports A"
        details = f"LRP1 below threshold (MCF7={mcf7_wt_tpm:.1f}, T47D={t47d_wt_tpm:.1f} TPM), Model B implausible"
        effect = 0.0
    elif expressed and regulated:
        verdict = "Supports B"
        details = f"LRP1 expressed (MCF7={mcf7_wt_tpm:.1f}, T47D={t47d_wt_tpm:.1f} TPM) and regulated by D538G"
        effect = max(abs(d['MCF7_cohens_d']), abs(d['T47D_cohens_d']))
    else:
        verdict = "Ambiguous"
        details = f"LRP1 expressed but not regulated (MCF7 FC={mcf7_fc:.2f}, T47D FC={t47d_fc:.2f})"
        effect = 0.0

    p_val = min(d['MCF7_pval'], d['T47D_pval'])
    print(f"\n  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(5, "LRP1 Expression", "B", verdict, effect, p_val, details, pd.DataFrame([d]))


def run_test_6_rap_gating(gene_data):
    """Test 6: Is LRPAP1 (RAP) highly expressed, suppressing premature binding?"""
    print("\n" + "-" * 60)
    print("TEST 6: RAP/LRPAP1 Gating")
    print("-" * 60)

    if 'LRPAP1' not in gene_data:
        print("  LRPAP1 NOT FOUND in expression data")
        return TestResult(6, "RAP/LRPAP1 Gating", "B", "Ambiguous", 0.0, 1.0,
                         "LRPAP1 not found in dataset", pd.DataFrame())

    d = gene_data['LRPAP1']
    mcf7_wt_tpm = d['MCF7_WT_TPM']
    t47d_wt_tpm = d['T47D_WT_TPM']
    mcf7_fc = d['MCF7_FC']
    t47d_fc = d['T47D_FC']

    print(f"  LRPAP1 MCF7: WT={mcf7_wt_tpm:.1f} TPM, D538G FC={mcf7_fc:.2f}")
    print(f"  LRPAP1 T47D: WT={t47d_wt_tpm:.1f} TPM, D538G FC={t47d_fc:.2f}")

    if 'LRP1' in gene_data:
        lrp1_tpm = max(gene_data['LRP1']['MCF7_WT_TPM'], gene_data['LRP1']['T47D_WT_TPM'])
        rap_to_lrp1 = max(mcf7_wt_tpm, t47d_wt_tpm) / lrp1_tpm if lrp1_tpm > 0 else float('inf')
        print(f"  LRPAP1/LRP1 ratio: {rap_to_lrp1:.1f}")

    highly_expressed = mcf7_wt_tpm > 10 and t47d_wt_tpm > 10
    stable = abs(d['MCF7_log2FC']) < np.log2(1.2) and abs(d['T47D_log2FC']) < np.log2(1.2)
    decreasing = d['MCF7_FC'] < 0.8 or d['T47D_FC'] < 0.8

    if highly_expressed and stable:
        verdict = "Supports A"
        details = f"LRPAP1 highly expressed ({mcf7_wt_tpm:.0f}/{t47d_wt_tpm:.0f} TPM) and stable, premature binding suppressed"
        effect = 0.0
    elif decreasing:
        verdict = "Supports B"
        details = f"LRPAP1 decreases (MCF7 FC={mcf7_fc:.2f}, T47D FC={t47d_fc:.2f}), trap may open"
        effect = max(abs(d['MCF7_cohens_d']), abs(d['T47D_cohens_d']))
    else:
        verdict = "Ambiguous"
        details = f"LRPAP1 inconclusive (MCF7={mcf7_wt_tpm:.0f} TPM FC={mcf7_fc:.2f}, T47D={t47d_wt_tpm:.0f} TPM FC={t47d_fc:.2f})"
        effect = 0.0

    p_val = min(d['MCF7_pval'], d['T47D_pval'])
    print(f"\n  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(6, "RAP/LRPAP1 Gating", "B", verdict, effect, p_val, details, pd.DataFrame([d]))


def run_test_7_coregulation(gene_data, hsp90b1_stats):
    """Test 7: Are LRP1 and HSP90B1 co-regulated? Is ER binding at LRP1?"""
    print("\n" + "-" * 60)
    print("TEST 7: LRP1-HSP90B1 Co-regulation + Epigenomics")
    print("-" * 60)

    gene_rows = []
    evidence_for_b = 0
    evidence_against_b = 0
    details_parts = []

    # Part A: Expression co-regulation
    if 'LRP1' in gene_data:
        lrp1 = gene_data['LRP1']
        gene_rows.append(lrp1)
        hsp_pattern = np.array([hsp90b1_stats['MCF7_log2FC'], hsp90b1_stats['T47D_log2FC']])
        lrp1_pattern = np.array([lrp1['MCF7_log2FC'], lrp1['T47D_log2FC']])

        if np.std(hsp_pattern) > 0 and np.std(lrp1_pattern) > 0:
            corr = np.corrcoef(hsp_pattern, lrp1_pattern)[0, 1]
        else:
            corr = 0

        print(f"  HSP90B1 pattern: MCF7={hsp90b1_stats['MCF7_log2FC']:.3f}, T47D={hsp90b1_stats['T47D_log2FC']:.3f}")
        print(f"  LRP1 pattern:    MCF7={lrp1['MCF7_log2FC']:.3f}, T47D={lrp1['T47D_log2FC']:.3f}")
        print(f"  Correlation: {corr:.2f}")

        if corr > 0.5:
            evidence_for_b += 1
            details_parts.append(f"LRP1-HSP90B1 positively correlated (r={corr:.2f})")
        else:
            evidence_against_b += 1
            details_parts.append(f"LRP1-HSP90B1 not positively correlated (r={corr:.2f})")
    else:
        evidence_against_b += 1
        details_parts.append("LRP1 not found in dataset")

    # Part B: ER ChIP-seq binding at LRP1 locus
    print("\n  Checking ER binding at LRP1 locus (GSE125117)...")

    if 'LRP1' in GENE_COORDS_HG38:
        chrom, start, end = GENE_COORDS_HG38['LRP1']
        peak_files = {
            'MCF7_WT': DATA_DIR / config['datasets']['GSE125117_MCF7_WT'],
            'MCF7_D538G': DATA_DIR / config['datasets']['GSE125117_MCF7_D538G'],
            'T47D_WT': DATA_DIR / config['datasets']['GSE125117_T47D_WT'],
            'T47D_D538G': DATA_DIR / config['datasets']['GSE125117_T47D_D538G'],
        }

        binding_found = False
        for name, filepath in peak_files.items():
            if filepath.exists():
                peaks = load_peaks(filepath)
                count = count_peaks_in_region(peaks, chrom, start, end)
                signal = get_peak_signal_in_region(peaks, chrom, start, end)
                print(f"    {name}: {count} peaks (signal={signal:.1f})")
                if count > 0:
                    binding_found = True

        if binding_found:
            evidence_for_b += 1
            details_parts.append("ER binds at LRP1 locus")
        else:
            evidence_against_b += 1
            details_parts.append("No ER binding at LRP1 locus")

    # Part C: ATAC-seq accessibility at LRP1
    print("\n  Checking ATAC accessibility at LRP1 (GSE254216)...")

    if 'LRP1' in GENE_COORDS_HG38:
        chrom, start, end = GENE_COORDS_HG38['LRP1']
        atac_files = {
            'MCF7': DATA_DIR / config['datasets']['GSE254216_MCF7'],
            'T47D': DATA_DIR / config['datasets']['GSE254216_T47D'],
        }
        for name, filepath in atac_files.items():
            if filepath.exists():
                atac_peaks = load_peaks(filepath)
                count = count_peaks_in_region(atac_peaks, chrom, start, end)
                print(f"    {name} ATAC at LRP1: {count} peaks")

    total_evidence = evidence_for_b + evidence_against_b
    b_fraction = evidence_for_b / total_evidence if total_evidence > 0 else 0

    if evidence_for_b >= 2:
        verdict = "Supports B"
    elif evidence_against_b >= 2:
        verdict = "Supports A"
    else:
        verdict = "Ambiguous"

    details = "; ".join(details_parts)
    print(f"\n  Evidence for B: {evidence_for_b}, against: {evidence_against_b}")
    print(f"  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(7, "LRP1-HSP90B1 Co-regulation", "B", verdict, b_fraction, 1.0, details, pd.DataFrame(gene_rows))


def run_test_8_consistency(previous_results):
    """Test 8: Does the leading model explain both cell lines consistently?"""
    print("\n" + "-" * 60)
    print("TEST 8: Cross-Cell-Line Consistency")
    print("-" * 60)

    a_count = sum(1 for r in previous_results if r.verdict == "Supports A")
    b_count = sum(1 for r in previous_results if r.verdict == "Supports B")
    print(f"  Tests 1-7 tally: Model A={a_count}, Model B={b_count}")

    if a_count > b_count:
        leading = "A"
        consistent = True
        explanation = ("Model A: HSP90B1 UP (MCF7) -> improved folding -> MDK secretion UP; "
                      "HSP90B1 DOWN (T47D) -> reduced folding -> MDK secretion DOWN. "
                      "Both lines explained by same mechanism.")
    elif b_count > a_count:
        leading = "B"
        consistent = True
        explanation = ("Model B (rescue): HSP90B1 UP (MCF7) -> better LRP1 chaperoning -> "
                      "less premature trapping -> MDK secretion UP; "
                      "HSP90B1 DOWN (T47D) -> worse chaperoning -> more trapping -> MDK secretion DOWN.")
    else:
        leading = "tie"
        consistent = False
        explanation = "Equal evidence for both models; cannot determine consistency."

    print(f"  Leading model: {leading}")
    print(f"  Consistent across cell lines: {consistent}")

    if consistent and leading != "tie":
        verdict = f"Supports {leading}"
        details = explanation
    else:
        verdict = "Ambiguous"
        details = explanation

    print(f"\n  VERDICT: {verdict}")
    print(f"  {details}")
    return TestResult(8, "Cross-Cell-Line Consistency", "both", verdict,
                     max(a_count, b_count) / 7 if a_count + b_count > 0 else 0,
                     1.0, details)


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def compile_scorecard(results):
    """Compile all test results into scorecard DataFrame and save CSV."""
    rows = []
    for r in results:
        rows.append({
            'test_id': r.test_id,
            'test_name': r.test_name,
            'model_tested': r.model_tested,
            'verdict': r.verdict,
            'effect_size': r.effect_size,
            'p_value': r.p_value,
            'details': r.details,
        })

    scorecard = pd.DataFrame(rows)

    # Apply Benjamini-Hochberg FDR correction across all tests
    _, fdr_corrected, _, _ = multipletests(scorecard['p_value'], method='fdr_bh')
    scorecard['p_value_BH'] = fdr_corrected

    scorecard.to_csv(OUTPUT_DIR / "tables" / "mechanism_discrimination_scorecard.csv", index=False)
    print(f"\nSaved scorecard to: outputs/tables/mechanism_discrimination_scorecard.csv")
    return scorecard


def plot_figure(tpm, gene_data, hsp90b1_stats, scorecard, results):
    """Generate 3x3 multi-panel figure."""
    print("\nGenerating fig13: Mechanism Discrimination...")

    fig = plt.figure(figsize=(20, 16))

    # ---- Panel A: ERAD gene heatmap ----
    ax_a = fig.add_subplot(3, 3, 1)
    erad_genes = [g for g in SIGNATURES['erad_genes'] if g in gene_data]
    if erad_genes:
        heatmap_data = pd.DataFrame({
            'MCF7 log2FC': [gene_data[g]['MCF7_log2FC'] for g in erad_genes],
            'T47D log2FC': [gene_data[g]['T47D_log2FC'] for g in erad_genes],
        }, index=erad_genes)
        annot_data = pd.DataFrame({
            'MCF7 log2FC': [f"{gene_data[g]['MCF7_FC']:.2f}" for g in erad_genes],
            'T47D log2FC': [f"{gene_data[g]['T47D_FC']:.2f}" for g in erad_genes],
        }, index=erad_genes)
        sns.heatmap(heatmap_data, annot=annot_data, fmt='', cmap='RdBu_r', center=0,
                    ax=ax_a, cbar_kws={'label': 'log2FC', 'shrink': 0.8})
    ax_a.set_title('A. ERAD Genes (Test 1)\nD538G/WT Fold Change', fontsize=10, fontweight='bold')

    # ---- Panel B: UPR marker bars ----
    ax_b = fig.add_subplot(3, 3, 2)
    upr_genes = [g for g in SIGNATURES['upr_genes'] if g in gene_data]
    if upr_genes:
        x = np.arange(len(upr_genes))
        width = 0.35
        mcf7_fcs = [gene_data[g]['MCF7_log2FC'] for g in upr_genes]
        t47d_fcs = [gene_data[g]['T47D_log2FC'] for g in upr_genes]
        ax_b.bar(x - width/2, mcf7_fcs, width, label='MCF7', color='steelblue', edgecolor='black')
        ax_b.bar(x + width/2, t47d_fcs, width, label='T47D', color='coral', edgecolor='black')
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(upr_genes, rotation=45, ha='right')
        ax_b.axhline(0, color='black', linewidth=0.5)
        ax_b.set_ylabel('log2(D538G/WT)')
        ax_b.legend(fontsize=8)
    ax_b.set_title('B. UPR Markers (Test 2)', fontsize=10, fontweight='bold')

    # ---- Panel C: Secretory capacity dot plot ----
    ax_c = fig.add_subplot(3, 3, 3)
    sec_genes = [g for g in SIGNATURES['secretory_chaperones'] if g in gene_data]
    all_sec = sec_genes + (['HSP90B1'] if 'HSP90B1' in gene_data else [])
    if all_sec:
        mcf7_vals = [gene_data[g]['MCF7_log2FC'] for g in all_sec]
        t47d_vals = [gene_data[g]['T47D_log2FC'] for g in all_sec]
        colors = ['gold' if g == 'HSP90B1' else 'steelblue' for g in all_sec]
        ax_c.scatter(mcf7_vals, t47d_vals, c=colors, s=80, edgecolors='black', zorder=3)
        for g, mx, ty in zip(all_sec, mcf7_vals, t47d_vals):
            ax_c.annotate(g, (mx, ty), fontsize=7, textcoords='offset points', xytext=(4, 4))
        ax_c.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax_c.axvline(0, color='gray', linestyle='--', alpha=0.5)
        ax_c.set_xlabel('MCF7 log2FC')
        ax_c.set_ylabel('T47D log2FC')
    ax_c.set_title('C. Secretory Chaperones (Test 3)\n(Gold=HSP90B1)', fontsize=10, fontweight='bold')

    # ---- Panel D: GRP94 clients ----
    ax_d = fig.add_subplot(3, 3, 4)
    client_genes = [g for g in SIGNATURES['grp94_clients'] if g in gene_data]
    if client_genes:
        y = np.arange(len(client_genes))
        mcf7_fcs = [gene_data[g]['MCF7_log2FC'] for g in client_genes]
        t47d_fcs = [gene_data[g]['T47D_log2FC'] for g in client_genes]
        ax_d.barh(y - 0.15, mcf7_fcs, 0.3, label='MCF7', color='steelblue', edgecolor='black')
        ax_d.barh(y + 0.15, t47d_fcs, 0.3, label='T47D', color='coral', edgecolor='black')
        ax_d.set_yticks(y)
        ax_d.set_yticklabels(client_genes)
        ax_d.axvline(0, color='black', linewidth=0.5)
        ax_d.set_xlabel('log2(D538G/WT)')
        ax_d.legend(fontsize=8)
    else:
        ax_d.text(0.5, 0.5, 'No GRP94 clients\nexpressed in dataset',
                 ha='center', va='center', transform=ax_d.transAxes, fontsize=11)
    ax_d.set_title('D. GRP94 Clients (Test 4)', fontsize=10, fontweight='bold')

    # ---- Panel E: LRP1/LRPAP1 expression bars ----
    ax_e = fig.add_subplot(3, 3, 5)
    lrp_genes = [g for g in SIGNATURES['lrp1_axis'] if g in gene_data]
    if lrp_genes:
        x = np.arange(len(lrp_genes))
        width = 0.2
        mcf7_wt = [gene_data[g]['MCF7_WT_TPM'] for g in lrp_genes]
        mcf7_d538g = [gene_data[g]['MCF7_D538G_TPM'] for g in lrp_genes]
        t47d_wt = [gene_data[g]['T47D_WT_TPM'] for g in lrp_genes]
        t47d_d538g = [gene_data[g]['T47D_D538G_TPM'] for g in lrp_genes]
        ax_e.bar(x - 1.5*width, mcf7_wt, width, label='MCF7 WT', color='steelblue', alpha=0.6, edgecolor='black')
        ax_e.bar(x - 0.5*width, mcf7_d538g, width, label='MCF7 D538G', color='steelblue', edgecolor='black')
        ax_e.bar(x + 0.5*width, t47d_wt, width, label='T47D WT', color='coral', alpha=0.6, edgecolor='black')
        ax_e.bar(x + 1.5*width, t47d_d538g, width, label='T47D D538G', color='coral', edgecolor='black')
        ax_e.set_xticks(x)
        ax_e.set_xticklabels(lrp_genes)
        ax_e.set_ylabel('TPM')
        ax_e.legend(fontsize=7, ncol=2)
    else:
        ax_e.text(0.5, 0.5, 'LRP1/LRPAP1 not found\nin dataset',
                 ha='center', va='center', transform=ax_e.transAxes, fontsize=11)
    ax_e.set_title('E. LRP1/LRPAP1 Expression (Tests 5-6)', fontsize=10, fontweight='bold')

    # ---- Panel F: LRP1 locus epigenomic strip ----
    ax_f = fig.add_subplot(3, 3, 6)
    if 'LRP1' in GENE_COORDS_HG38:
        chrom, start, end = GENE_COORDS_HG38['LRP1']
        conditions = ['MCF7_WT', 'MCF7_D538G', 'T47D_WT', 'T47D_D538G']
        chipseq_keys = ['GSE125117_MCF7_WT', 'GSE125117_MCF7_D538G',
                       'GSE125117_T47D_WT', 'GSE125117_T47D_D538G']
        peak_counts = []
        signals = []
        for key in chipseq_keys:
            filepath = DATA_DIR / config['datasets'][key]
            if filepath.exists():
                peaks = load_peaks(filepath)
                peak_counts.append(count_peaks_in_region(peaks, chrom, start, end))
                signals.append(get_peak_signal_in_region(peaks, chrom, start, end))
            else:
                peak_counts.append(0)
                signals.append(0)

        colors_bar = ['steelblue', 'steelblue', 'coral', 'coral']
        alphas = [0.5, 1.0, 0.5, 1.0]
        for i in range(len(conditions)):
            ax_f.bar(i, signals[i], color=colors_bar[i], alpha=alphas[i], edgecolor='black')
            ax_f.text(i, signals[i] + 1, f'{peak_counts[i]}pk', ha='center', fontsize=8)
        ax_f.set_xticks(range(len(conditions)))
        ax_f.set_xticklabels(['MCF7\nWT', 'MCF7\nD538G', 'T47D\nWT', 'T47D\nD538G'], fontsize=8)
        ax_f.set_ylabel('ER ChIP-seq Signal')
    ax_f.set_title('F. ER Binding at LRP1 Locus (Test 7)\n(GSE125117)', fontsize=10, fontweight='bold')

    # ---- Panel G: Correlation matrix ----
    ax_g = fig.add_subplot(3, 3, 7)
    corr_genes = ['HSP90B1', 'MDK', 'LRP1', 'LRPAP1']
    corr_genes_present = [g for g in corr_genes if g in gene_data]
    if len(corr_genes_present) >= 2:
        fc_data = {}
        for g in corr_genes_present:
            fc_data[g] = [gene_data[g]['MCF7_log2FC'], gene_data[g]['T47D_log2FC']]
        fc_df = pd.DataFrame(fc_data)
        corr_matrix = fc_df.corr()
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                    ax=ax_g, vmin=-1, vmax=1, square=True)
    else:
        ax_g.text(0.5, 0.5, 'Insufficient genes\nfor correlation',
                 ha='center', va='center', transform=ax_g.transAxes, fontsize=11)
    ax_g.set_title('G. Log2FC Correlation Matrix (Test 7)', fontsize=10, fontweight='bold')

    # ---- Panel H: Consistency diagram ----
    ax_h = fig.add_subplot(3, 3, 8)
    ax_h.axis('off')
    test8 = results[7] if len(results) >= 8 else None
    consistency_text = ("CELL LINE CONSISTENCY CHECK (Test 8)\n"
                       "=" * 45 + "\n\n"
                       "Known phenotypes:\n"
                       "  MCF7-D538G:  HSP90B1 UP   -> MDK secretion UP (+83%)\n"
                       "  T47D-D538G:  HSP90B1 DOWN -> MDK secretion DOWN (-62%)\n\n"
                       "Model A (Direct Folding):\n"
                       "  MCF7:  HSP90B1 UP -> folding capacity UP -> MDK export UP\n"
                       "  T47D:  HSP90B1 DN -> folding capacity DN -> MDK export DN\n"
                       "  -> Consistent across both lines\n\n"
                       "Model B (LRP1 Trapping - Rescue):\n"
                       "  MCF7:  HSP90B1 UP -> LRP1 chaperoning UP -> less trapping\n"
                       "  T47D:  HSP90B1 DN -> LRP1 chaperoning DN -> more trapping\n"
                       "  -> Also consistent (rescue sub-model)\n\n"
                       f"Leading model: {'A' if test8 and 'A' in test8.verdict else 'B' if test8 and 'B' in test8.verdict else 'Tie'}")
    ax_h.text(0.02, 0.98, consistency_text, transform=ax_h.transAxes, fontsize=8.5,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # ---- Panel I: Scorecard summary ----
    ax_i = fig.add_subplot(3, 3, 9)
    ax_i.axis('off')
    table_data = []
    row_colors = []
    for _, row in scorecard.iterrows():
        p_str = f"{row['p_value']:.2e}" if row['p_value'] < 0.05 else f"{row['p_value']:.2f}"
        table_data.append([
            f"T{int(row['test_id'])}",
            row['test_name'][:28],
            row['verdict'],
            f"{row['effect_size']:.2f}",
            p_str,
        ])
        if 'A' in row['verdict']:
            row_colors.append('#cce5ff')
        elif 'B' in row['verdict']:
            row_colors.append('#ffe5cc')
        else:
            row_colors.append('#e0e0e0')

    table = ax_i.table(
        cellText=table_data,
        colLabels=['#', 'Test', 'Verdict', 'Effect', 'p-value'],
        cellLoc='center', loc='center',
        colWidths=[0.06, 0.38, 0.22, 0.14, 0.20],
    )
    for i, color in enumerate(row_colors):
        for j in range(5):
            table[i + 1, j].set_facecolor(color)
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.4)

    a_count = sum(1 for r in results if r.verdict == "Supports A")
    b_count = sum(1 for r in results if r.verdict == "Supports B")
    ax_i.set_title(f'I. Evidence Scorecard\nModel A: {a_count}/8 | Model B: {b_count}/8',
                   fontsize=10, fontweight='bold')

    plt.tight_layout()
    fig_path = OUTPUT_DIR / "figures" / "fig13_mechanism_discrimination.png"
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / "figures" / "fig13_mechanism_discrimination.pdf", bbox_inches='tight')
    plt.close()
    print(f"Saved figure: {fig_path}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 80)
    print("SCRIPT 13: HSP90B1 MECHANISM DISCRIMINATION")
    print("  Model A: Direct folding/export (ER chaperone capacity)")
    print("  Model B: LRP1 premature binding trap")
    print("=" * 80)

    # Load data
    print("\nLoading GSE89888 RNA-seq data...")
    tpm = load_tpm()
    all_samples = [s for g in GROUPS.values() for s in g if s in tpm.columns]
    tpm = tpm[tpm[all_samples].max(axis=1) > 0.1]
    print(f"Loaded {len(tpm)} genes")

    # Get HSP90B1 baseline
    hsp90b1_stats = get_gene_fc_and_stats(tpm, 'HSP90B1')
    if hsp90b1_stats is None:
        sys.exit("ERROR: HSP90B1 not found in expression data")

    print(f"\nHSP90B1 reference:")
    print(f"  MCF7: FC={hsp90b1_stats['MCF7_FC']:.2f} (D538G {'UP' if hsp90b1_stats['MCF7_FC'] > 1 else 'DOWN'})")
    print(f"  T47D: FC={hsp90b1_stats['T47D_FC']:.2f} (D538G {'UP' if hsp90b1_stats['T47D_FC'] > 1 else 'DOWN'})")

    # Compute stats for all signature genes
    print("\nComputing expression for all signature genes...")
    all_sig_genes = set()
    for gene_list in SIGNATURES.values():
        all_sig_genes.update(gene_list)
    all_sig_genes.add('HSP90B1')
    all_sig_genes.add('MDK')

    gene_data = {}
    for gene in sorted(all_sig_genes):
        result = get_gene_fc_and_stats(tpm, gene)
        if result is not None:
            gene_data[gene] = result

    print(f"Found {len(gene_data)}/{len(all_sig_genes)} signature genes in dataset")
    missing = all_sig_genes - set(gene_data.keys())
    if missing:
        print(f"Missing: {', '.join(sorted(missing))}")

    # Run tests
    results = []
    results.append(run_test_1_erad(gene_data, hsp90b1_stats))
    results.append(run_test_2_upr(gene_data, hsp90b1_stats))
    results.append(run_test_3_secretory(gene_data, hsp90b1_stats))
    results.append(run_test_4_clients(gene_data, hsp90b1_stats))
    results.append(run_test_5_lrp1_expression(gene_data))
    results.append(run_test_6_rap_gating(gene_data))
    results.append(run_test_7_coregulation(gene_data, hsp90b1_stats))
    results.append(run_test_8_consistency(results))

    # Compile scorecard
    scorecard = compile_scorecard(results)

    # Save gene signatures table
    sig_df = pd.DataFrame(list(gene_data.values()))
    sig_df.to_csv(OUTPUT_DIR / "tables" / "mechanism_gene_signatures.csv", index=False)

    # Generate figure
    plot_figure(tpm, gene_data, hsp90b1_stats, scorecard, results)

    # Summary
    print("\n" + "=" * 80)
    print("SCORECARD SUMMARY")
    print("=" * 80)
    a_count = sum(1 for r in results if r.verdict == "Supports A")
    b_count = sum(1 for r in results if r.verdict == "Supports B")
    amb_count = sum(1 for r in results if r.verdict == "Ambiguous")
    print(f"  Model A (Direct folding): {a_count}/8 tests")
    print(f"  Model B (LRP1 trapping):  {b_count}/8 tests")
    print(f"  Ambiguous:                {amb_count}/8 tests")

    print("\n" + "=" * 80)
    print("SCRIPT 13 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
