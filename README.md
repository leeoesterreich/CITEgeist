# CITEgeist

[![GitHub stars](https://img.shields.io/github/stars/leeoesterreich/CITEgeist)](https://github.com/leeoesterreich/CITEgeist/stargazers)
[![License](https://img.shields.io/github/license/leeoesterreich/CITEgeist)](https://github.com/leeoesterreich/CITEgeist/blob/main/LICENSE)
[![Issues](https://img.shields.io/github/issues/leeoesterreich/CITEgeist)](https://github.com/leeoesterreich/CITEgeist/issues)
[![Python Version](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)

> Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends

CITEgeist is a comprehensive computational framework for analyzing spatial multi-omic data, with a focus on integrating CITE-seq and spatial transcriptomics. Our toolkit enables researchers to uncover spatial patterns in cellular organization and gene expression, providing deeper insights into tissue architecture and function.

## 🚀 Key Features

- Integration of CITE-seq and spatial transcriptomics data
- Robust cell-type deconvolution
- Spatial pattern analysis
- Flexible regularization options
- Scalable to large datasets
- Comprehensive visualization tools

## 📚 Documentation

- [Quick Start Guide](CITEgeist/README.md)
- [Benchmarking Results](Benchmarking/README.md)
- [Example Notebooks](CITEgeist/examples)

## 💻 System Requirements

### Minimum Requirements
- Python 3.10
- 16GB RAM
- Multi-core processor
- Linux, macOS, or Windows 10 with WSL2

## 🔧 Quick Installation

You can install CITEgeist using pip:

```bash
pip install citegeist
```

For development installation:

```bash
git clone https://github.com/alee-x/CITEgeist.git
cd CITEgeist
pip install -e .[dev]
```

## 📜 License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details. This license ensures that all modifications and derivative works must also be open source.

## 📫 Contact

- **Lab Website**: [Lee/Oesterreich Laboratory](https://leeoesterreich.org/)
- **Issues**: [GitHub Issues](https://github.com/leeoesterreich/CITEgeist/issues)
- **Email**: [Contact Us](mailto:alc376@pitt.edu)

## 📝 Citation

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

## 📝 Citation

If you use CITEgeist in your research, please cite:

[Citation information to be added]

## 📝 Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

## 📝 Support

For support, please open an issue on our [GitHub issue tracker](https://github.com/alee-x/CITEgeist/issues).

---
Copyright (c) 2025 Lee/Oesterreich Lab

