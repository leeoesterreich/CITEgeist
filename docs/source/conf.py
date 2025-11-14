# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the project root to sys.path
project_root = os.path.abspath('../../')
sys.path.insert(0, project_root)
sys.path.insert(0, os.path.join(project_root, 'CITEgeist'))
sys.path.insert(0, os.path.join(project_root, 'CITEgeist', 'model'))
print("!!! DEBUG: The project root path added is:", project_root)
print("!!! DEBUG: The full sys.path is:", sys.path)


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CITEgeist'
copyright = '2025, Lee/Oesterreich Lab'
author = 'Lee/Oesterreich Lab'
release = ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx_autodoc_typehints',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

html_theme_options = {
    "github_url": "ll ghttps://github.com/leeoesterreich/CITEgeist",
    "show_nav_level": 2,
    "navigation_depth": 4,
    "show_toc_level": 2,
    "navbar_align": "left",
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["navbar-icon-links"],
    "footer_start": ["copyright", "sphinx-version"],
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

