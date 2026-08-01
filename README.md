# CITEgeist

[![GitHub stars](https://img.shields.io/github/stars/leeoesterreich/CITEgeist)](https://github.com/leeoesterreich/CITEgeist/stargazers)
[![License](https://img.shields.io/github/license/leeoesterreich/CITEgeist)](https://github.com/leeoesterreich/CITEgeist/blob/main/LICENSE)
[![Issues](https://img.shields.io/github/issues/leeoesterreich/CITEgeist)](https://github.com/leeoesterreich/CITEgeist/issues)
[![Python Version](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![CI](https://github.com/leeoesterreich/CITEgeist/actions/workflows/ci.yml/badge.svg)](https://github.com/leeoesterreich/CITEgeist/actions/workflows/ci.yml)
[![Code Quality](https://github.com/leeoesterreich/CITEgeist/actions/workflows/quality.yml/badge.svg)](https://github.com/leeoesterreich/CITEgeist/actions/workflows/quality.yml)

> Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends

CITEgeist is a comprehensive computational framework for analyzing spatial multi-omic data, with a focus on integrating CITE-seq and spatial transcriptomics. Our toolkit enables researchers to uncover spatial patterns in cellular organization and gene expression, providing deeper insights into tissue architecture and function.

## Key Features

- Reference-free cell-type deconvolution using same-slide CITE-seq antibody capture — no scRNA-seq reference required
- GPU-accelerated cell-type proportion estimation via quadratic programming (cuOPT)
- SACE per-cell gene expression deconvolution via single-pass Poisson-multinomial allocation
- Per-cell type assignment from spot proportions via StarDist nuclei segmentation. The default
  routes nucleus morphology scores (`morphology/nucleus_scores.py`, shipped) through a Bayesian
  posterior; count-constrained Hungarian matching remains available via `use_morphology=False`
- Spatial gene program discovery and cross-sample integration (Modules 4–5)
- Validated on dense tumor microenvironments including breast cancer and RCC clinical samples

## Documentation

- [Quick Start Guide](https://github.com/leeoesterreich/CITEgeist/blob/main/CITEgeist/README.md)
- [Running on Real Visium Data (end-to-end tutorial)](https://github.com/leeoesterreich/CITEgeist/blob/main/docs/quickstart_real_visium.md)
- [Benchmarking Methodology](https://github.com/leeoesterreich/CITEgeist/blob/main/docs/Benchmarking/README.md)
- [Example Scripts](https://github.com/leeoesterreich/CITEgeist/tree/main/examples/scripts)

## Data Availability

| Dataset | Used for | Source |
|---|---|---|
| Visium + CITE-seq clinical-trial samples (this study) | Midkine / ESR1-D538G analyses; competitor-GEX validation | NCBI GEO [GSE289326](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289326) |
| Wu et al. 2021 breast cancer scRNA-seq atlas | Reference for reference-based competitor methods; source cells for the simulation | NCBI GEO [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) |
| Xenium Human Kidney Dataset (renal cell carcinoma) | Real-tissue benchmark | [10x Genomics](https://www.10xgenomics.com/datasets), publicly available |

## Reproducibility

Example end-to-end pipelines for each module are in [`examples/scripts`](https://github.com/leeoesterreich/CITEgeist/tree/main/examples/scripts).
See the [end-to-end Visium tutorial](https://github.com/leeoesterreich/CITEgeist/blob/main/docs/quickstart_real_visium.md) to run CITEgeist on real data.

The curated analysis code behind the paper's reported numbers — with a map from each script to
the figure it produces — is in
[`examples/paper_analysis`](https://github.com/leeoesterreich/CITEgeist/tree/main/examples/paper_analysis).

## System Requirements

### Minimum Requirements
- Python 3.10
- 16GB RAM
- Multi-core processor
- NVIDIA GPU with 8GB+ VRAM (required for cuOPT QP solver, Module 3)
- Linux, macOS, or Windows 10 with WSL2

## Quick Installation

CITEgeist is published on PyPI:

```bash
pip install citegeist
```

> **Note — package name vs import name:** The installed distribution is named
> `citegeist` (lowercase), but the importable package is `CITEgeist` (capital C). After installation:
>
> ```python
> import CITEgeist                          # bare import (lightweight)
> from CITEgeist import CitegeistModel      # lazy re-export — heavy stack loaded on demand
> ```

> **GPU requirement:** CITEgeist's QP deconvolution (Module 3) requires an **NVIDIA GPU with
> CUDA** and the [cuOPT](https://developer.nvidia.com/cuopt-logistics-optimization) library.
> cuOPT is **not available on PyPI** — install it via NVIDIA RAPIDS or NGC
> (`pip install cuopt-cu12` from the NVIDIA index, or use a pre-built RAPIDS container).
> CPU-only environments **cannot run** the cell-type proportion solver.

To install from source instead:

```bash
git clone https://github.com/leeoesterreich/CITEgeist.git
cd CITEgeist
pip install -e .
```

For a development installation with test and lint dependencies:

```bash
pip install -e .[dev]
```

## Development Setup

### For Contributors

If you want to contribute to CITEgeist, follow these steps:

```bash
# Clone the repository
git clone https://github.com/leeoesterreich/CITEgeist.git
cd CITEgeist

# Create conda environment
conda env create -f CITEgeist_env.yml
conda activate CITEgeist_env

# Install in development mode with dev dependencies
pip install -e ".[dev]"
```

### Code Quality Tools

We use several tools to maintain code quality:

- **black**: Code formatting (120 char line length)
- **isort**: Import sorting
- **flake8**: Linting
- **mypy**: Type checking
- **pylint**: Static code analysis
- **pytest**: Testing framework

Run quality checks before committing:
```bash
# Format code
black CITEgeist/
isort CITEgeist/

# Run linting
flake8 CITEgeist/
pylint CITEgeist/

# Run tests
pytest
```

**CI enforcement (blocks sync to main):** black, isort, flake8, pytest/coverage

**Advisory (badge only, no merge gate):** pylint, mypy

See [CONTRIBUTING.md](https://github.com/leeoesterreich/CITEgeist/blob/main/CONTRIBUTING.md) for detailed contribution guidelines.

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](https://github.com/leeoesterreich/CITEgeist/blob/main/LICENSE) file for details. BSD 3-Clause is a permissive license: you may use, modify, and redistribute the code, including in proprietary work, provided you retain the copyright notice and do not use the authors' names to endorse derived products.

## Contact

- **Lab Website**: [Lee/Oesterreich Laboratory](https://leeoesterreich.org/)
- **Issues**: [GitHub Issues](https://github.com/leeoesterreich/CITEgeist/issues)
- **Email**: [Contact Us](mailto:leeav@upmc.edu)

## Citation

If you use CITEgeist in your research, please cite our paper:
(to be updated)
```bibtex
@article{ChangSchlegelCITEgeistCellularIndexing2025,
  title = {{{CITEgeist}}: {{Cellular Indexing}} of {{Transcriptomes}} and {{Epitopes}} for {{Guided Exploration}} of {{Intrinsic Spatial Trends}}},
  shorttitle = {{{CITEgeist}}},
  author = {Chang, Alexander Chih-Chieh and Schlegel, Brent T. and Carleton, Neil and McAulife, Priscilla F. and Oesterreich, Steffi and Schwartz, Russell and Lee, Adrian V.},
  date = {2025-02-17},
  eprinttype = {bioRxiv},
  eprintclass = {New Results},
  pages = {2025.02.15.638331},
  doi = {10.1101/2025.02.15.638331},
  url = {https://www.biorxiv.org/content/10.1101/2025.02.15.638331v1},
  urldate = {2025-02-17},
  abstract = {Spatial transcriptomics provides insights into tissue architecture by linking gene expression with spatial localization. Current deconvolution methods rely heavily on single-cell RNA sequencing (scRNA-seq) references, which are costly and often unavailable, mainly if the tissue under evaluation is limited, such as in a core biopsy specimen. We present a novel tool, CITEgeist, that deconvolutes spatial transcriptomics data using antibody capture from the same slide as the reference, directly leveraging cell surface protein measurements from the same tissue section. This approach circumvents the limitations of scRNA-seq as a reference, offering a cost-effective and biologically grounded alternative. Our method employs mathematical optimization to estimate cell type proportions and gene expression profiles, incorporating sparsity constraints for robustness and interpretability. Benchmarks against state-of-the-art deconvolution methods show improved accuracy in cell type resolution, particularly in dense tumor microenvironments, while maintaining computational efficiency. This antibody-based tool advances spatial transcriptomics by providing a scalable, accurate, and reference-independent solution for deconvolution in complex tissues. We validate this tool by using a combined approach of simulated data and clinical samples by applying CITEgeist to translational pre-treatment and post-treatment ER+ breast tumors from an ongoing clinical trial, emphasizing the applicability and robustness of CITEgeist.},
}


```

## Contributing

We welcome contributions! Please see our [CONTRIBUTING.md](https://github.com/leeoesterreich/CITEgeist/blob/main/CONTRIBUTING.md) for detailed guidelines on:
- Setting up your development environment
- Code quality standards
- Testing requirements
- Pull request process

## Support

For support, please:
- Open an issue on our [GitHub issue tracker](https://github.com/leeoesterreich/CITEgeist/issues)
- Check existing issues and discussions
- Read our documentation and examples

---
Copyright (c) 2025 Lee/Oesterreich Lab
