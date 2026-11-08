"""TRAMINEX_BINDING and TRAMINEX_NUMBA are read the same way (TX-AUD-20261019-031, TX-AUD-20261011-049).

The comment on TRAMINEX_BINDING says it is validated "exactly like TRAMINEX_NUMBA", but the binding refused
" python " and "" while the Numba switch strips spaces and reads "" as "auto"; and an invalid TRAMINEX_NUMBA was
accepted silently when numba does not import, because the check sat after the import short-circuit. Each case
runs in a fresh interpreter, since both variables are read at import.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pytest

PROBE = (
    "import traminex.distances as d, traminex.distances.python_kernels as pk; "
    "print(d._impl_name(), pk._HAVE_NUMBA)"
)


def _run(code: str, **env: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-c", code], env={**os.environ, **env}, capture_output=True, text=True
    )


@pytest.mark.parametrize(
    ("value", "impl"), [(" python ", "python"), ("", None), ("PYTHON", "python")]
)
def test_binding_is_stripped_and_empty_means_auto(value: str, impl: str | None) -> None:
    result = _run(PROBE, TRAMINEX_BINDING=value)
    assert result.returncode == 0, result.stderr
    if impl is not None:
        assert result.stdout.split()[0] == impl
    else:
        assert result.stdout.split()[0] == _run(PROBE, TRAMINEX_BINDING="auto").stdout.split()[0]


def test_an_invalid_binding_is_still_refused() -> None:
    result = _run(PROBE, TRAMINEX_BINDING="bogus")
    assert result.returncode != 0 and "TRAMINEX_BINDING must be one of" in result.stderr


@pytest.mark.parametrize("numba_present", [True, False], ids=["numba-present", "numba-blocked"])
def test_an_invalid_numba_mode_is_refused_whether_or_not_numba_imports(numba_present: bool) -> None:
    prefix = "" if numba_present else "import sys; sys.modules['numba'] = None; "
    result = _run(prefix + PROBE, TRAMINEX_NUMBA="bogus")
    assert result.returncode != 0
    assert "TRAMINEX_NUMBA must be one of" in result.stderr
