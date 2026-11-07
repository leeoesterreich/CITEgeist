# Changelog

## 0.2.0

One install that reproduces the ERT paper: v1 proportions plus SACE.

- `citegeist.core`: the v1 proportion code (Gurobi QP with radius-400 refinement). The six
  `.py` files are byte-identical to PyPI `citegeist==0.1.1`, so `from citegeist.core.citegeist_model
  import CitegeistModel` keeps working.
- `citegeist.sace`: SACE per-cell expression allocation, taken from the 2.0.0 release (tag `v2.0.0`).
  Only the modules SACE needs are included. The only source change is the import path
  (`CITEgeist.model.*` became `citegeist.sace.*`), plus removal of lazy-export entries whose
  module is not shipped. No second top-level `CITEgeist` package is installed.
- Packaging: requirements are the union of 0.1.1 and 2.0.0. `psutil` and `statsmodels`, which
  `citegeist.core` always imports, are now declared. Python >= 3.10.
- Pin the version: `pip install citegeist==0.2.0`. A plain `pip install citegeist` installs the
  newer 2.x line, which has a different layout (`CITEgeist.model`).

## 0.1.1

v1 proportion code on PyPI (Gurobi QP).

## 2.0.0

SACE release (package name `citegeist`, import name `CITEgeist`).
