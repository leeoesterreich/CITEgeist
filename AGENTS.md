# Repository Guidelines

## Project Structure & Module Organization
- Core package lives in `CITEgeist/model/` (`citegeist_model.py`, `gurobi_impl.py`, `utils.py`, `checkpoints.py`); import via `from CITEgeist.model import CitegeistModel`.
- Tests sit in `tests/` (unit + integration); fixture config in `tests/conftest.py`; coverage outputs to `htmlcov/` and `coverage.xml`.
- Example scripts and notebooks are under `examples/`; benchmarking code and reports live in `Benchmarking/` and `benchmarking_logs/`.
- Documentation lives in `docs/`; legacy or archived models reside in `model_archive/`. Use `output/` and `intermediate_files/` for generated artifacts and avoid committing them.

## Build, Test, and Development Commands
```bash
# Set up dev env
pip install -e .[dev]
pre-commit install

# Format and lint
black CITEgeist/ --line-length 120
isort CITEgeist/ --profile black
flake8 CITEgeist/ --max-line-length=120 --extend-ignore=E203,E501,W503
mypy CITEgeist/
pylint CITEgeist/

# Run tests (includes coverage by default via pyproject addopts)
pytest
pytest -m "not slow"           # skip slow cases
pytest --cov=CITEgeist/model --cov-report=html
```
Use the provided `CITEgeist_env.yml` if you prefer conda.

## Coding Style & Naming Conventions
- Python 3.10+, 4-space indentation, 120-char lines (Black enforced), trailing commas allowed for multi-line lists.
- Naming: `snake_case` for functions/vars, `PascalCase` for classes, `UPPER_SNAKE_CASE` for constants, `_private` for internals.
- Write docstrings (Google/NumPy style) explaining intent; cite algorithms when relevant. Keep imports sorted (isort) and type hints where clarity improves.

## Testing Guidelines
- Tests should live alongside peers in `tests/` and follow `test_*.py` / `test_*` naming.
- Pytest markers in use: `unit`, `integration`, `slow`, `requires_gurobi`, `requires_data`; prefer `-m "not slow"` for quick checks.
- CI enforces coverage ≥60%; target higher for critical modules (`checkpoints`, `utils`).
- Use fixtures from `conftest.py` for AnnData mocks and temp directories; avoid relying on real patient data.

## Commit & Pull Request Guidelines
- Commit messages: short imperative summary on the first line (≈50 chars), optional bullet body; example `Add feature: cell type deconvolution`.
- Before pushing, run `pre-commit run --all-files` and `pytest` to match CI.
- PRs should describe intent, list key changes, link related issues, and include screenshots or brief logs for user-facing or result changes. Keep branches up to date with `main` and respond to review feedback promptly.

## Security & Configuration Tips
- Do not commit licenses, credentials, or large datasets; Gurobi license files should stay local and tests requiring them are marked `requires_gurobi`.
- Generated outputs in `output/`, `intermediate_files/`, and benchmarking artifacts should remain untracked unless explicitly needed for documentation.
