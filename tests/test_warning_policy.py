"""A diagnostic a test does not state fails it (TX-AUD-20261019-068).

The suite emitted 269 warnings from 72 sites, mostly R's cost diagnostics from tests that
feed unusual costs on purpose, so a new, unexpected diagnostic was one more line in the
summary. ``pyproject.toml`` now turns ``UserWarning`` (traminex's ``CostDiagnosticMessage``
is one) and ``DeprecationWarning`` into errors; a test that meets a diagnostic by design
states it with ``pytest.warns`` or a ``filterwarnings`` mark naming the message.
"""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

PROBE = """
import numpy as np
import pytest

import traminex as tx


def _distance():
    x = tx.seqdef(np.array([list("AAB"), list("ABB")], dtype=object))
    return tx.seqdist(x, method="OM", subm=np.array([[0.0, 5.0], [5.0, 0.0]]), indel=1.0)

"""


def test_the_suite_turns_warnings_into_errors() -> None:
    options = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    filters = options["tool"]["pytest"]["ini_options"]["filterwarnings"]

    assert "error::UserWarning" in filters
    assert "error::DeprecationWarning" in filters
    # A deprecation raised inside matplotlib's own modules is matplotlib's (TX-REP-20261025-001).
    assert r"ignore::DeprecationWarning:matplotlib(\..*)?$" in filters


def _run(tmp_path: Path, test: str) -> subprocess.CompletedProcess:
    (tmp_path / "test_probe.py").write_text(PROBE + test, encoding="utf-8")
    return subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         "-c", str(ROOT / "pyproject.toml"), "--rootdir", str(tmp_path),
         str(tmp_path / "test_probe.py")],  # fmt: skip
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_an_unstated_cost_diagnostic_fails_the_test(tmp_path) -> None:
    result = _run(tmp_path, "\ndef test_unstated():\n    _distance()\n")

    assert result.returncode == 1, result.stdout + result.stderr
    assert "CostDiagnosticMessage" in result.stdout
    assert "greater than twice the indel cost" in result.stdout


def test_a_stated_cost_diagnostic_passes(tmp_path) -> None:
    result = _run(
        tmp_path,
        "\n@pytest.mark.filterwarnings(\n"
        '    r"ignore:\\s*\\[!!\\] At least one substitution cost greater than twice the indel cost"\n'
        '    ":traminex.distances.python_kernels.CostDiagnosticMessage"\n'
        ")\ndef test_marked():\n    _distance()\n\n\n"
        "def test_asserted():\n"
        "    from traminex.distances.python_kernels import CostDiagnosticMessage\n\n"
        '    with pytest.warns(CostDiagnosticMessage, match="greater than twice the indel cost"):\n'
        "        _distance()\n",
    )

    assert result.returncode == 0, result.stdout + result.stderr


# warnings.warn_explicit attributes the warning to the named module, as a library's own warn() call
# would; the configured filters then decide.
_DEPRECATION = (
    "\nimport warnings\n\n\ndef test_warned():\n"
    "    warnings.warn_explicit(\"'parseString' deprecated - use 'parse_string'\", DeprecationWarning, "
    "'probe.py', 1, module={module!r}, registry={{}})\n"
)


def test_a_deprecation_raised_inside_matplotlib_does_not_fail_a_test(tmp_path) -> None:
    # At the declared floors, pyparsing 3.3 deprecations raised inside matplotlib 3.6.0's own modules
    # stopped collection of 72 test modules (TX-REP-20261025-001).
    result = _run(tmp_path, _DEPRECATION.format(module="matplotlib._fontconfig_pattern"))

    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("module", ["traminex.distances", "matplotlib_inline.backend_inline"])
def test_a_deprecation_raised_elsewhere_still_fails_a_test(tmp_path, module: str) -> None:
    result = _run(tmp_path, _DEPRECATION.format(module=module))

    assert result.returncode == 1, result.stdout + result.stderr
    assert "DeprecationWarning" in result.stdout


def test_contributing_states_the_matplotlib_exception_with_its_finding() -> None:
    text = " ".join((ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8").split())

    assert (
        "matplotlib 3.6.0 calls pyparsing names that pyparsing 3.3 deprecates (TX-REP-20261025-001)"
        in text
    )
