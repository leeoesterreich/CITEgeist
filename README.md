# citegeist 0.2.0 (v1 proportions + SACE, for the ERT paper)

> Install this release with a pin: `pip install citegeist==0.2.0`.
> `pip install citegeist` without a pin installs the newer 2.x line, which has a different layout
> (the `CITEgeist.model` package) and does not contain the v1 proportion code.

CITEgeist is the Lee/Oesterreich lab framework for spatial multi-omics (CITE-seq antibody capture
with Visium). Release 0.2.0 puts two earlier releases in a single install so that the ERT paper
needs one dependency instead of two:

| Import | What it is | Same as |
|---|---|---|
| `citegeist.core` | v1 cell-type proportions: Gurobi quadratic program plus radius-400 refinement | PyPI `citegeist==0.1.1` (the six files are byte-identical) |
| `citegeist.sace` | SACE per-cell gene-expression allocation from spot proportions | `citegeist==2.0.0` (`CITEgeist.model`, only import paths changed) |

## Requirements

- Python >= 3.10
- `gurobipy` and a Gurobi licence for the v1 proportions (`citegeist.core`). SACE does not use Gurobi.
- The paper results were run with numpy 1.26.4, pandas 2.2.3, scanpy 1.10.4, squidpy 1.6.2 and
  gurobipy 11.0.2.

## Usage

Proportions (v1):

```python
from citegeist.core.citegeist_model import CitegeistModel
```

SACE (previously `from CITEgeist.model.citegeist_model import CitegeistModel` with citegeist 2.0.0):

```python
from citegeist.sace.citegeist_model import CitegeistModel
# ... same calls as in 2.0.0, for example
# spot_type_gex, diagnostics = model.run_sace_allocation(output_mode="layers", init_method="marker")
```

The SACE function code is `citegeist.sace.gex.sace_gex` (`run_sace_layers`, `run_sace`,
`project_sace_to_cells`). The lazy exports of `citegeist.sace` also work, for example
`from citegeist.sace import CitegeistModel, run_sace_layers`.

`citegeist.core.CitegeistModel` (v1) and `citegeist.sace.CitegeistModel` (2.0.0) are different classes.
Use the first for proportions and the second for SACE.

## Newer code

Cell-type proportions by SPOTS marker share with RNA refinement, the cuOPT solver, the QC tools
and the full 2.x documentation are in the 2.x releases of this package, not in 0.2.0.

## Licence

BSD 3-Clause, see `LICENSE`.
