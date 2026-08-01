# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the project root to sys.path
project_root = os.path.abspath("../../")
sys.path.insert(0, project_root)
sys.path.insert(0, os.path.join(project_root, "CITEgeist"))
sys.path.insert(0, os.path.join(project_root, "CITEgeist", "model"))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "CITEgeist"
copyright = "2025, Lee/Oesterreich Lab"
author = "Lee/Oesterreich Lab"
release = "2.0.1"
version = "2.0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
exclude_patterns = []

# The docs build does not install CITEgeist or its runtime stack: cuOPT is
# GPU-only and not on PyPI, and torch/opencv/stardist would exceed Read the
# Docs' build resources. autodoc imports each module for its docstrings, so
# every third-party package NOT in docs/requirements.txt must be mocked here.
# numpy, pandas, scipy and anndata are installed for real so that signatures
# render true annotations instead of mock reprs.
autodoc_mock_imports = [
    "csbdeep",
    "cuopt",
    "cv2",
    "esda",
    "gurobipy",
    "joblib",
    "libpysal",
    "matplotlib",
    "networkx",
    "PIL",
    "psutil",
    # NOT pyarrow: pandas parses pyarrow.__version__ with a regex at its own
    # import time, and a mock's __version__ is a Mock, not a string — mocking it
    # raises TypeError and breaks pandas entirely, which silently empties every
    # autodoc page while sphinx-build still exits 0. The package already guards
    # its pyarrow import (citegeist_model.py:19), so absence is handled.
    "scanpy",
    "skimage",
    "sklearn",
    "squidpy",
    "stardist",
    "statsmodels",
    "timm",
    "torch",
    "torch_geometric",
    "torch_scatter",
    "torchvision",
    "tqdm",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "navigation_depth": 4,
    "titles_only": False,
}

# GitHub context for edit page button
html_context = {
    "github_user": "leeoesterreich",
    "github_repo": "CITEgeist",
    "github_version": "main",
    "doc_path": "docs/source",
}

# -- Options for autosummary -------------------------------------------------
autosummary_generate = True
autosummary_imported_members = True

# -- Options for Napoleon ----------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True
