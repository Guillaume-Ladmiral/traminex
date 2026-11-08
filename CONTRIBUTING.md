# Contributing to Traminex

Thank you for your interest in contributing to traminex! This document provides guidelines and instructions for contributing to the project.

## 📋 Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Testing](#testing)
- [Code Style](#code-style)
- [Submitting Changes](#submitting-changes)
- [Documentation](#documentation)
- [Areas for Contribution](#areas-for-contribution)

---

## Code of Conduct

This project follows standard open-source community guidelines:
- Be respectful and inclusive
- Focus on constructive feedback
- Help create a welcoming environment for all contributors

---

## Getting Started

### Prerequisites

- Ubuntu 24.04
- Python 3.11 or 3.12 (the family-standard development environment uses 3.12)
- The `~/miniforge3/envs/datasci` environment, or a private virtual
  environment layered on it that leaves it unchanged (docs/DEV_BUILD.md,
  "Environment"), for the local gates
- Git
- Familiarity with sequence analysis (helpful but not required)

### Quick Start

1. **Fork the repository** on GitHub
2. **Clone your fork**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/traminex.git
   cd traminex
   ```
3. **Set up development environment** (see below)
4. **Create a branch** for your changes
5. **Make your changes** following our guidelines
6. **Submit a pull request**

---

## Development Setup

### Basic Installation

```bash
# Install in editable mode with all extras in the standard environment
~/miniforge3/envs/datasci/bin/python -m pip install -e '.[plot,perf,extras,tree,dev]'
```

### Optional: C++ Extensions

For performance-critical development, you can build optional C++ extensions:

```bash
# See docs/DEV_BUILD.md for detailed instructions
PY=~/miniforge3/envs/datasci/bin/python
PYBIND11_DIR=$($PY -m pybind11 --cmakedir)
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DTRAMINEX_BUILD_PYBIND11=ON \
  -DPython_EXECUTABLE="$PY" \
  -Dpybind11_DIR="$PYBIND11_DIR"
cmake --build build --parallel
```

### Verify Installation

```bash
# Run tests
make test-python
make test-pybind11

# Check imports
~/miniforge3/envs/datasci/bin/python -c "import traminex; print('OK')"
```

Release-hardening work must also verify artifacts outside the checkout. Build
one default sdist, one default pure wheel, and one explicit compiled wheel, then
run the maintained six-profile verifier described in
[`docs/DEV_BUILD.md`](docs/DEV_BUILD.md). It checks Python 3.11 floor and 3.12
current installs, every optional dependency boundary, installed import origin,
PEP 561 typing, pybind11/OpenMP identity, and fallback behavior.

---

## Making Changes

### Branch Naming

Use descriptive branch names:
- `feature/add-seqfoo-function` - New features
- `fix/issue-123-distance-bug` - Bug fixes
- `docs/update-distance-guide` - Documentation
- `refactor/split-stats-module` - Code refactoring

### Commit Messages

Write clear, descriptive commit messages:

```
Add seqfoo function for feature extraction

- Implement core algorithm from TraMineR seqfoo
- Add comprehensive tests with edge cases
- Update API_PARITY.md to reflect new function
- Add usage examples to docs/USAGE_new_apis.md

Closes #123
```

**Format**:
- First line: Brief summary (50 chars or less)
- Blank line
- Detailed description (wrapped at 72 chars)
- Reference related issues

---

## Testing

### Running Tests

```bash
# All tests
make test-python

# Build pybind11 with OpenMP and run its complete suite
make test-pybind11

# Specific test file
pytest tests/test_stats.py -v

# Specific test function
pytest tests/test_stats.py::test_seqtab_counts_and_props -v

# With coverage
TRAMINEX_BINDING=python ~/miniforge3/envs/datasci/bin/python \
  -m pytest --cov=traminex --cov-branch --cov-report=html
```

### Writing Tests

- Place tests in `tests/` directory
- Name test files `test_*.py`
- Name test functions `test_*`
- Use descriptive test names that explain what is being tested
- A warning fails the suite (`UserWarning`, which includes traminex's
  `CostDiagnosticMessage`, and `DeprecationWarning`; see `pyproject.toml`). A test
  that triggers one of R's diagnostics on purpose states it, with
  `pytest.warns(..., match=...)` or a
  `@pytest.mark.filterwarnings("ignore:<message regex>:<category>")` mark. A
  `DeprecationWarning` raised inside matplotlib's own modules does not fail it: at the
  declared floor, matplotlib 3.6.0 calls pyparsing names that pyparsing 3.3 deprecates
  (TX-REP-20261025-001)
- After adding, removing or renaming a file, stage the change (`git add` or `git rm`)
  and run `make sdist-manifest`: it rewrites the sdist list in `pyproject.toml` from
  the tracked files, and `tests/test_packaging_profiles.py` fails while the list is out
  of date

**Example**:

```python
def test_seqfoo_with_missing_states():
    """Test that seqfoo correctly handles missing states."""
    data = np.array([['A', '*', 'B'], ['A', 'B', '*']])
    sts = seqdef(data, nr='*')
    result = seqfoo(sts, include_missing=False)
    assert result.shape == (2, 2)  # Missing states excluded
```

### Test Coverage Guidelines

- New behavior must cover its material success, validation, and fallback branches
- Include edge cases: empty input, single element, large datasets
- Test parameter validation and error messages
- Test TraMineR parity where applicable

---

## Code Style

### Python Style

We follow **PEP 8** with some project-specific conventions:

```bash
# Format code with black
black traminex/ tests/

# Lint with ruff
ruff check traminex/ tests/ --fix

# Run pre-commit hooks
pre-commit run --all-files

# Check every configured package module
make typecheck

# Check the 226-name root stub, bidirectional module/signature symmetry,
# and that the stub tree type-checks under mypy
make stub-check
```

### Style Guidelines

1. **Line length**: 100 characters (set for black and ruff in `pyproject.toml`;
   black's own default is 88)
2. **Imports**:
   - Standard library first
   - Third-party packages second
   - Local imports last
   - Use `from __future__ import annotations` for type hints
3. **Type hints**: Use type hints for function signatures
4. **Docstrings**: NumPy-style docstrings for all public functions
5. **Function size**: a function stays within 200 lines and an approximate
   complexity of 50. The 22 longer ones are recorded with their sizes in
   `tests/test_function_size_limit.py`, which fails when one grows or shrinks
   without its recorded size changing; when you shrink one, record its new size
   there. The test cannot stop one commit that raises a recorded size together
   with the code: that shows in the diff as an edit to the table, and review must
   catch it.

**Example**:

```python
def seqfoo(
    x: STSList,
    *,
    param1: float = 1.0,
    include_missing: bool = False,
) -> np.ndarray:
    """Brief one-line description.

    More detailed description if needed. Explain what the function does,
    when to use it, and any important notes.

    Parameters
    ----------
    x : STSList
        Input sequence object
    param1 : float, default 1.0
        Description of parameter
    include_missing : bool, default False
        Whether to include missing/void states

    Returns
    -------
    np.ndarray
        Description of return value, including shape

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqfoo
    >>> data = np.array([['A', 'B', 'C']])
    >>> sts = seqdef(data)
    >>> result = seqfoo(sts)
    >>> print(result.shape)
    (1, 3)

    Notes
    -----
    Additional notes about implementation, algorithms, or references.
    """
    # Implementation here
    pass
```

### Naming Conventions

- **Functions**: `snake_case` (TraMineR-compatible names like `seqdef`, `seqdist`)
- **Classes**: `PascalCase` (`STSList`, `SeqAlign`)
- **Constants**: `UPPER_SNAKE_CASE`
- **Private functions**: `_leading_underscore`

---

## Submitting Changes

### Pull Request Process

1. **Update your branch**:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Ensure all tests pass**:
   ```bash
   make release-check
   make test-python
   make test-pybind11
   ```

3. **Update documentation**:
   - Add docstrings to new functions
   - Update relevant docs in `docs/`
   - Add examples if applicable
   - Update `CHANGELOG.md`

4. **Create pull request**:
   - Provide clear description of changes
   - Reference related issues
   - Include test results
   - Add screenshots for UI/plotting changes

5. **Code review**:
   - Address reviewer feedback
   - Update PR based on comments
   - Re-request review when ready

### Pull Request Template

```markdown
## Description
Brief description of changes

## Related Issues
Closes #123

## Changes Made
- Added seqfoo function
- Updated API_PARITY.md
- Added comprehensive tests

## Testing
- [ ] All existing tests pass
- [ ] Added new tests for new functionality
- [ ] Tested manually with sample data

## Documentation
- [ ] Updated docstrings
- [ ] Updated relevant guides in docs/
- [ ] Added usage examples
- [ ] Updated CHANGELOG.md

## Checklist
- [ ] Code follows style guidelines
- [ ] Tests added and passing
- [ ] Documentation updated
- [ ] No breaking changes (or documented if necessary)
```

---

## Documentation

### Documentation Structure

- **`README.md`**: Quick start, installation, basic examples
- **`docs/`**: Detailed guides and API documentation
- **Docstrings**: In-code documentation for all public functions
- **Examples**: Jupyter notebooks and Python scripts in `examples/`

### Writing Documentation

1. **Be clear and concise**: Avoid jargon where possible
2. **Include examples**: Show concrete usage
3. **Link to related docs**: Cross-reference relevant guides
4. **Keep it up-to-date**: Update docs when code changes

### Documentation Standards

- Use **Markdown** for all documentation files
- Follow **NumPy docstring style** for Python docstrings
- Include **runnable examples** in docstrings
- Add **cross-references** to TraMineR equivalents where applicable

---

## Areas for Contribution

### Priority Areas

Based on our [ROADMAP.md](ROADMAP.md), [docs/API_PARITY.md](docs/API_PARITY.md), and `reports/audit/traminex_api_audit.md`:

#### High-value work

`seqauto`, `seqimplic`, and the weighted single-group `seqsurv` path are already
implemented and R-golden tested. Current priorities are:

1. Complete the remaining explicitly documented `seqformat` conversion modes.
2. Extend `seqsurv` only where a new path can be pinned to the R oracle.
3. Keep all six isolated Ubuntu install profiles and the pybind11/OpenMP
   artifact profile reproducible.
4. Expand property tests, honest branch coverage, and benchmark agreement
   gates without inventing analytical expectations that lack R-oracle evidence.
5. Keep generated parity and API documentation synchronized with the canonical
   manifest and current R oracle.

#### Documentation Improvements

- Add more .qmd tutorial vignettes
- Translate TraMineR examples to traminex
- Create video walkthroughs
- Improve API reference completeness

#### Performance Optimization

- Profile bottlenecks outside the already-accelerated OM/LCS/HAM paths
- Add accelerated paths only with byte-identity or R-golden evidence
- Implement caching for expensive operations
- Sparse matrix support for large datasets

#### Testing & Quality

- Add more integration tests
- Increase integration and workflow coverage
- Add property-based tests with Hypothesis
- Performance regression tests

#### Examples & Use Cases

- Life course analysis workflows
- Event sequence mining examples
- Clustering best practices
- Real-world datasets and tutorials

### Good First Issues

Great for new contributors:

1. **Add missing examples** to docstrings
2. **Fix typos** in documentation
3. **Improve error messages** for better UX
4. **Add parameter validation** to existing functions
5. **Write tutorial notebooks** for specific use cases

---

## Project Structure

```
traminex/
├── traminex/              # Main package
│   ├── __init__.py        # Public API exports
│   ├── format.py          # Data format conversions
│   ├── sequences.py       # Core STSList class
│   ├── stats.py           # Statistics facade
│   ├── stats_basic.py     # Basic counts/distributions
│   ├── stats_complexity.py # Turbulence/complexity metrics
│   ├── stats_utils.py     # Helper/indicator utilities
│   ├── stats_p1.py        # Additional parity statistics
│   ├── stats_p2.py        # Additional parity statistics
│   ├── distances/         # Distance algorithms
│   │   ├── __init__.py
│   │   ├── costs.py       # Cost matrices
│   │   ├── omloc.py       # Position-dependent OM
│   │   └── python_kernels.py  # Distance implementations
│   ├── visualization/     # Plotting subsystem
│   └── ...
├── tests/                 # Test suite
├── docs/                  # Documentation
├── examples/              # Example scripts and notebooks
├── benchmarks/            # Performance benchmarks
└── scripts/               # Utility scripts
```

---

## Getting Help

- **Documentation**: Check [docs/README.md](docs/README.md)
- **Development**: See [docs/DOCS_POLICY.md](docs/DOCS_POLICY.md)
- **Questions**: Open a GitHub Discussion
- **Bugs**: Open a GitHub Issue with reproduction steps

---

## License

By contributing to traminex, you agree that your contributions will be licensed under the **GPL-3.0-or-later** license.

---

## Recognition

Contributors will be:
- Acknowledged in release notes
- Credited in relevant documentation

Thank you for contributing to traminex.

---

**Last Updated**: 2026-07-15
**Version**: 0.2.2
