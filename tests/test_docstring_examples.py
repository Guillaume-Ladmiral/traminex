"""Every docstring example in traminex runs and shows what it prints (TX-AUD-20261010-097).

The 2026-10-10 audit ran doctest over the package: 17 of 108 examples failed (names never
imported, undefined variables, NumPy 2 reprs) and nothing ran them. Each module's examples
run here as doctest runs them, in a temporary working directory, so an example that saves
a figure writes nothing into the checkout.
"""

from __future__ import annotations

import doctest
import importlib
import io
import os
import pkgutil

import pytest

import traminex

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

PACKAGE_DIR = os.path.dirname(traminex.__file__)
MODULES = ["traminex"] + sorted(
    module.name for module in pkgutil.walk_packages(traminex.__path__, "traminex.")
)


@pytest.mark.parametrize("name", MODULES)
def test_every_docstring_example_runs(name: str, tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    module = importlib.import_module(name)
    finder = doctest.DocTestFinder(exclude_empty=True)
    tests = [
        test
        for test in finder.find(module, name)
        if test.examples and (test.filename or "").startswith(PACKAGE_DIR)
    ]
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
    report = io.StringIO()
    for test in tests:
        runner.run(test, out=report.write)
    result = runner.summarize(verbose=False)
    assert result.failed == 0, report.getvalue()
