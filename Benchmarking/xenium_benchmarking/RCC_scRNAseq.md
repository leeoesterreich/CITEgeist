# RCC Single-Cell RNA-seq Reference Datasets for Spatial Deconvolution

## Purpose
Reference datasets for benchmarking spatial deconvolution methods (Cell2Location, RCTD, Tangram, Seurat) on the Xenium-derived pseudo-Visium test dataset.

---

## Top Recommended Datasets

| Dataset | Cells | Samples | Best For | Source |
|---------|-------|---------|----------|--------|
| **GSE156632** | 54,776 | 7 KIRC + 5 normal | Comprehensive TME | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156632) |
| **GSE159115** | ~75,000 | Tumor + normal + PBMC | Immune profiling | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159115) |
| **Sanger Kidney Atlas** | Large | Multi-patient | Pre-processed h5ad | [Mendeley](https://doi.org/10.17632/g67bkbnhhg.1) |
| **GSE121636** | Medium | 3 ccRCC patients | TIL analysis | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121636) |

---

## All Available GEO Datasets

### Primary Recommendations
- **GSE156632**: 54,776 cells from 7 KIRC + 5 adjacent normal samples (Kim et al.)
- **GSE159115**: Tumor, normal adjacent tissue, and peripheral blood
- **GSE210042**: ccRCC tissues and benign kidney
- **GSE178481**: ccRCC single-cell data

### Additional Datasets
- **GSE237425**: Tumor and matched normal endothelial cells from ccRCC - reference-quality atlas
- **GSE237427**: RNAseq data for kidney cancer endothelial cells
- **GSE121638**: Borcherding et al. kidney cancer dataset
- **GSE154763**: Cheng et al. kidney cancer dataset
- **GSE152938**: 4 tumor samples (2 ccRCC, 1 ChRCC, 1 pRCC)
- **GSE202374**: Single-cell RNA sequencing data (Chow J)
- **GSE224630**: Processed scRNA-seq data (Saout JR)
- **GSE207493**: Yu Z's single-cell RNA sequencing data
- **GSE139555**: Peripheral blood, adjacent tissues, and tumor tissues from 3 RCC patients
- **GSE171306**: Bilateral RCC - 3,575 cells (left) + 3,568 cells (right)
- **GSE73122**: Intratumoral heterogeneity of primary RCC and lung metastasis

### European Genome-Phenome Archive
- **EGAD00001008029**: Whole-exome sequencing
- **EGAD00001008030**: Single cell RNA sequencing
- **EGAD00001008781**: Spatial transcriptomic data

---

## Cell Types in ccRCC Tumor Microenvironment

Typical cell types identified in ccRCC scRNA-seq studies:

1. **Tumor cells** (proximal tubule origin)
2. **CD8+ T cells**
3. **CD4+ T cells**
4. **Macrophages / Myeloid cells**
5. **B cells**
6. **Plasma cells**
7. **NK cells**
8. **Endothelial cells**
9. **Fibroblasts / CAFs**
10. **Dendritic cells**

---

## Key Publications

### 2024-2025 Studies
1. **Metabolic heterogeneity in ccRCC** (J Transl Med, Feb 2024)
   - Gathered 9 major scRNA-seq databases, 195 samples
   - Combined with spatial transcriptomics
   - [Link](https://link.springer.com/article/10.1186/s12967-024-04848-x)

2. **Single-cell multiomics in ccRCC** (Cell Discovery, 2022)
   - scRNA-seq + scATAC-seq
   - [Link](https://www.nature.com/articles/s41421-022-00415-0)

3. **Mapping the TME in ccRCC** (Front Genet, 2023)
   - Identified 11 cell types
   - Focus on tumor cells, myeloid, T cells, fibroblasts, ECs
   - [Link](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2023.1207233/full)

### Foundational Studies
- **Young et al. (Science, 2018)**: 72,501 single-cell transcriptomes of human renal tumors
- **Zhang et al. (PNAS, 2021)**: Cell atlases predicting cells of origin for 10+ RCC subtypes

---

## Deconvolution Method Benchmarks

From [Frontiers benchmarking study](https://www.frontiersin.org/articles/10.3389/fbinf.2024.1352594/full):

| Method | Accuracy | Speed | Resources | Notes |
|--------|----------|-------|-----------|-------|
| **Cell2Location** | Highest | Slow | Heavy | Best overall accuracy |
| **RCTD** | High | Fast | Light | Best speed/accuracy balance |
| **Tangram** | High | Medium | Medium | Nearly as accurate as Cell2Location |
| **Seurat** | Medium | Fast | Light | Simple workflow |

Key findings:
- Cell2Location and RCTD are top-performing methods
- Choice of reference dataset has largest impact on predictions
- RCTD shows best performance in CVD/CKD samples
- Both RCTD and Cell2Location reduce platform effects

---

## Data Resources

### Web Portals
- **Sanger Kidney Cancer Portal**: https://www.sanger.ac.uk/project/microenvironment-of-kidney-cancer
- **Interactive RCC Endothelial Atlas**: https://yxu2.shinyapps.io/shinyapp/
- **Kidney Interactive Transcriptomics**: https://humphreyslab.com/SingleCell/

### Code Repositories
- RCC Endothelial Analysis: https://github.com/YuexinXu/RCC_Endo

### Pre-processed Data
- H5AD objects available at Mendeley: https://doi.org/10.17632/g67bkbnhhg.1

---

## Recommended Workflow

1. **Download GSE156632** as primary reference (54,776 cells, good cell type diversity)
2. **Process with Scanpy** to create reference signatures
3. **Run deconvolution methods**:
   - Cell2Location
   - RCTD
   - Tangram
   - Seurat (TransferData)
4. **Benchmark against Xenium ground truth** using same metrics as CITEgeist

---

*Last updated: 2025-12-27*
