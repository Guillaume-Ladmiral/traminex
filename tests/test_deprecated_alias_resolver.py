from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from traminex._deprecated import MISSING, resolve_deprecated_alias


def _resolve(canonical, alias, *, default=MISSING):
    return resolve_deprecated_alias(
        canonical,
        alias,
        canonical_name="canonical",
        alias_name="legacy",
        function_name="example",
        default=default,
    )


def test_canonical_only_and_required_omission_contracts() -> None:
    value = object()
    assert _resolve(value, MISSING) is value

    with pytest.raises(TypeError, match="missing required argument"):
        _resolve(MISSING, MISSING)


def test_alias_only_warns_at_the_caller_and_preserves_identity() -> None:
    value = np.array([1.0, np.nan])

    with pytest.warns(DeprecationWarning, match="legacy.*canonical") as caught:
        result = _resolve(MISSING, value)

    assert result is value
    # Resolved: through a symbolic link the code's file name and __file__ can name
    # the same file by two paths (TX-AUD-20261010-090 (a)).
    assert Path(caught[0].filename).resolve() == Path(__file__).resolve()


@pytest.mark.parametrize(
    ("canonical", "default", "alias"),
    [
        (float("nan"), float("nan"), 3.0),
        (np.array([1.0, np.nan]), np.array([1.0, np.nan]), np.array([2.0])),
        (
            {"values": np.array([1.0, np.nan]), "label": "x"},
            {"values": np.array([1.0, np.nan]), "label": "x"},
            {"values": np.array([2.0]), "label": "y"},
        ),
    ],
)
def test_explicit_canonical_equal_to_default_uses_alias(canonical, default, alias) -> None:
    canonical_before = deepcopy(canonical)
    default_before = deepcopy(default)
    alias_before = deepcopy(alias)

    with pytest.warns(DeprecationWarning):
        result = _resolve(canonical, alias, default=default)

    assert result is alias
    _assert_unchanged(canonical, canonical_before)
    _assert_unchanged(default, default_before)
    _assert_unchanged(alias, alias_before)


@pytest.mark.parametrize(
    ("canonical", "default", "alias"),
    [
        (2, 1, 3),
        (np.array([1.0, 2.0]), np.array([1.0, 3.0]), np.array([4.0])),
        (
            {"values": np.array([1.0]), "label": "x"},
            {"values": np.array([1.0]), "label": "y"},
            {"values": np.array([2.0]), "label": "z"},
        ),
    ],
)
def test_nondefault_canonical_and_alias_conflict(canonical, default, alias) -> None:
    with pytest.raises(ValueError, match="canonical.*legacy"):
        _resolve(canonical, alias, default=default)


def _assert_unchanged(actual, expected) -> None:
    if isinstance(actual, np.ndarray):
        np.testing.assert_array_equal(actual, expected)
        return
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for key in actual:
            _assert_unchanged(actual[key], expected[key])
        return
    if isinstance(actual, float) and np.isnan(actual):
        assert isinstance(expected, float) and np.isnan(expected)
        return
    assert actual == expected


def _small_sequences():
    import traminex as tx

    x = tx.seqdef(np.array([["A", "B", "B"], ["A", "A", "B"], ["B", "B", "A"]], dtype=object))
    return x, tx.seqdist(x, method="HAM")


def test_wrapped_public_function_attributes_alias_warning_to_its_caller() -> None:
    """A decorator between the caller and the resolver must not hide the warning.

    ``seqrep`` is wrapped by ``sequences._require_stslist``, so the resolver's
    fixed ``stacklevel=3`` attributed ``'trep' is deprecated for seqrep()`` to
    ``traminex/sequences.py``, where Python's default filters hide a
    ``DeprecationWarning`` (audit S7-07, 2026-09-22). R always shows it. The
    warning now points at the first frame outside the package.
    """
    import traminex as tx

    x, diss = _small_sequences()
    with pytest.warns(DeprecationWarning, match="'trep' is deprecated for seqrep") as caught:
        tx.seqrep(x, diss=diss, trep=0.5)

    # Resolved: through a symbolic link the code's file name and __file__ can name
    # the same file by two paths (TX-AUD-20261010-090 (a)).
    assert Path(caught[0].filename).resolve() == Path(__file__).resolve()


def test_alias_warning_of_a_wrapped_function_is_shown_by_default_filters() -> None:
    import os
    import subprocess
    import sys

    code = (
        "import numpy as np\n"
        "import traminex as tx\n"
        "x = tx.seqdef(np.array([['A','B','B'],['A','A','B'],['B','B','A']], dtype=object))\n"
        "tx.seqrep(x, diss=tx.seqdist(x, method='HAM'), trep=0.5)\n"
    )
    env = {key: value for key, value in os.environ.items() if key != "PYTHONWARNINGS"}
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).resolve().parents[1],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "'trep' is deprecated for seqrep()" in result.stderr, result.stderr
