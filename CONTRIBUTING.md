# Contributing to CITEgeist

Thank you for your interest in contributing to CITEgeist! This document provides guidelines and instructions for contributing to the project.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Code Quality Standards](#code-quality-standards)
- [Testing](#testing)
- [Pull Request Process](#pull-request-process)
- [Coding Conventions](#coding-conventions)

## Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please be respectful and constructive in all interactions.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/CITEgeist.git
   cd CITEgeist
   ```
3. **Add upstream remote**:
   ```bash
   git remote add upstream https://github.com/leeoesterreich/CITEgeist.git
   ```

## Development Setup

### Option 1: Conda Environment (Recommended)

```bash
# Create and activate conda environment
conda env create -f CITEgeist_env.yml
conda activate CITEgeist_env

# Install development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Option 2: pip + venv

```bash
# Create virtual environment
python3.10 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install package with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Verify Installation

```bash
# Run a quick test
python -c "from CITEgeist import CitegeistModel; print('Installation successful!')"
```

## Code Quality Standards

We use automated tools to maintain code quality. All contributions must pass these checks.

### Pre-commit Hooks

Pre-commit hooks automatically run when you commit. They check:
- Code formatting (black, isort)
- Linting (flake8)
- Type hints (mypy)
- Common issues (trailing whitespace, large files, etc.)

**Run manually on all files:**
```bash
pre-commit run --all-files
```

### Code Formatting

We use **black** for consistent formatting:
```bash
# Format your code
black CITEgeist/ --line-length=120

# Check without modifying
black --check CITEgeist/
```

### Import Sorting

We use **isort** to organize imports:
```bash
# Sort imports
isort CITEgeist/ --profile black

# Check without modifying
isort --check-only CITEgeist/
```

### Linting

We use **flake8** for linting:
```bash
flake8 CITEgeist/ --max-line-length=120 --extend-ignore=E203,E501,W503
```

### Type Checking

We use **mypy** for type checking:
```bash
mypy CITEgeist/ --ignore-missing-imports
```

### Running All Quality Checks

```bash
# Run all checks at once
black --check CITEgeist/
isort --check-only CITEgeist/
flake8 CITEgeist/
mypy CITEgeist/
pylint CITEgeist/
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=CITEgeist --cov-report=html

# Run specific test file
pytest CITEgeist/tests/test_citegeist_simulated.py

# Run tests matching a pattern
pytest -k "test_model"
```

### Writing Tests

- Place tests in `CITEgeist/tests/`
- Name test files `test_*.py`
- Name test functions `test_*`
- Use descriptive names that explain what is being tested

Example:
```python
def test_model_initialization():
    """Test that CitegeistModel initializes correctly."""
    model = CitegeistModel(sample_name="test")
    assert model.sample_name == "test"
```

## Pull Request Process

### 1. Create a Feature Branch

```bash
# Update your fork
git fetch upstream
git checkout main
git merge upstream/main

# Create feature branch
git checkout -b feature/your-feature-name
```

### 2. Make Your Changes

- Write clean, well-documented code
- Add tests for new functionality
- Update documentation as needed
- Follow our coding conventions

### 3. Test Your Changes

```bash
# Run tests
pytest

# Run quality checks
pre-commit run --all-files

# Or run specific checks
black --check CITEgeist/
isort --check-only CITEgeist/
flake8 CITEgeist/
```

### 4. Commit Your Changes

```bash
# Stage changes
git add .

# Commit (pre-commit hooks will run automatically)
git commit -m "Add feature: brief description"

# If pre-commit modifies files, review and commit again
git add .
git commit -m "Add feature: brief description"
```

### 5. Push and Create PR

```bash
# Push to your fork
git push origin feature/your-feature-name
```

Then create a Pull Request on GitHub:
- Provide a clear description of changes
- Reference any related issues
- Ensure all CI checks pass

### 6. Code Review

- Address reviewer feedback promptly
- Push additional commits to the same branch
- Keep the conversation constructive

## Coding Conventions

### Python Style

- **Line length**: 120 characters
- **Indentation**: 4 spaces (no tabs)
- **Naming conventions**:
  - Functions/variables: `snake_case`
  - Classes: `PascalCase`
  - Constants: `UPPER_SNAKE_CASE`
  - Private members: `_leading_underscore`

### Documentation

- **Docstrings**: Use Google or NumPy style
- **Comments**: Explain *why*, not *what*
- **Type hints**: Use where it improves clarity

Example:
```python
def optimize_cell_proportions(
    data: np.ndarray,
    reference: pd.DataFrame,
    regularization: float = 0.1
) -> Dict[str, np.ndarray]:
    """
    Optimize cell type proportions using Gurobi.

    Args:
        data: Spatial expression data (n_spots x n_genes).
        reference: Reference cell type profiles (n_genes x n_celltypes).
        regularization: L1 regularization strength.

    Returns:
        Dictionary containing optimized proportions and metadata.
    """
    # Implementation
    pass
```

### Scientific Code Specifics

- **Single-letter variables** are acceptable for mathematical notation (x, y, i, j, k)
- **Matrix/array dimensions** should be commented clearly
- **Algorithm references** should be cited in docstrings
- **Numerical stability** considerations should be documented

### Git Commit Messages

Good commit messages help track project history:

```
Add feature: cell type deconvolution with spatial regularization

- Implement Gurobi optimization for spatial constraints
- Add tests for boundary conditions
- Update documentation with usage examples

Fixes #123
```

Format:
- **First line**: Brief summary (50 chars or less)
- **Body**: Detailed explanation if needed
- **Footer**: Reference issues/PRs

## Project Structure

```
CITEgeist/
├── CITEgeist/          # Main package
│   ├── model/          # Core model implementation
│   ├── tests/          # Test files
│   └── examples/       # Example notebooks
├── Benchmarking/       # Benchmarking code
├── .github/            # GitHub workflows
│   └── workflows/      # CI/CD configurations
├── docs/               # Documentation
└── tests/              # Additional tests
```

## Questions?

- **Issues**: [GitHub Issues](https://github.com/leeoesterreich/CITEgeist/issues)
- **Discussions**: [GitHub Discussions](https://github.com/leeoesterreich/CITEgeist/discussions)
- **Email**: alc376@pitt.edu

## License

By contributing, you agree that your contributions will be licensed under the BSD-3-Clause License.

---

Thank you for contributing to CITEgeist! 🎉
