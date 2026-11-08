"""A missing optional dependency must name the extra that installs it (F-014).

``README.md`` advertises the extras by name, but a user who installed only the
core got ``ImportError: No module named 'matplotlib'`` at attribute resolution,
which never connects back to them.  The scikit-learn path was already exemplary
(``RuntimeError: scikit-learn is required for disstree_fit``); the lazy facade
resolver now matches that wording for every optional distribution.

The absence is produced in a subprocess with a real ``sys.meta_path`` blocker, so
the test exercises the actual import failure rather than a patched stand-in.
"""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

_BLOCKER = """
import sys

class _Block:
    def __init__(self, blocked):
        self.blocked = blocked

    def find_spec(self, name, path=None, target=None):
        root = name.split(".")[0]
        if root in self.blocked:
            raise ImportError("No module named %r" % name, name=root)
        return None

sys.meta_path.insert(0, _Block({BLOCKED!r}))
sys.path.insert(1, {ROOT!r})

import traminex

print("EXPORTS", len(traminex.__all__))
try:
    getattr(traminex, {EXPORT!r})
except Exception as exc:
    print("RAISED", type(exc).__name__, exc)
else:
    print("RESOLVED")
"""


def _resolve_with_blocked(blocked: set[str], export: str) -> str:
    script = _BLOCKER.format(BLOCKED=blocked, ROOT=str(ROOT), EXPORT=export)
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    )
    return result.stdout


@pytest.mark.parametrize("export", ["seqplot", "sequence_index_plot", "state_legend", "chronogram"])
def test_absent_matplotlib_names_the_plot_extra(export):
    output = _resolve_with_blocked({"matplotlib", "seaborn"}, export)

    assert "EXPORTS 226" in output, output
    assert "RAISED ImportError" in output, output
    assert "matplotlib is required for" in output, output
    assert "'plot' extra" in output, output
    assert "traminex[plot]" in output, output


_INTROSPECTION = """
import sys

class _Block:
    def __init__(self, blocked):
        self.blocked = blocked

    def find_spec(self, name, path=None, target=None):
        root = name.split(".")[0]
        if root in self.blocked:
            raise ImportError("No module named %r" % name, name=root)
        return None

sys.meta_path.insert(0, _Block({BLOCKED!r}))
sys.path.insert(1, {ROOT!r})

import inspect
import pydoc

import traminex

listed = set(dir(traminex))
unresolvable = set()
for name in traminex.__all__:
    try:
        getattr(traminex, name)
    except ImportError:
        unresolvable.add(name)
print("ALL", len(traminex.__all__))
print("UNRESOLVABLE", len(unresolvable))
print("UNRESOLVABLE_LISTED", sorted(unresolvable & listed))
print("RESOLVABLE_UNLISTED", sorted(set(traminex.__all__) - unresolvable - listed))
for label, call in (
    ("GETMEMBERS", lambda: inspect.getmembers(traminex)),
    ("RENDER_DOC", lambda: pydoc.render_doc(traminex)),
):
    try:
        call()
    except Exception as exc:
        print(label, "RAISED", type(exc).__name__, exc)
    else:
        print(label, "OK")
"""


@pytest.mark.parametrize(
    "blocked",
    [{"matplotlib", "seaborn"}, {"seaborn"}, {"numba"}, {"sklearn"}],
    ids=["matplotlib+seaborn", "seaborn", "numba", "sklearn"],
)
def test_introspection_skips_only_the_names_a_missing_extra_disables(blocked):
    """``help(traminex)`` and ``inspect.getmembers`` work on a core install (S7-05).

    Both iterate ``dir(traminex)`` and call ``getattr`` on every name, skipping
    only ``AttributeError``. Without matplotlib the 37 plotting names raised the
    informative ``ImportError`` above, so both failed on a core install (audit
    S7-05, 2026-09-22, reproduced in a real environment without matplotlib).
    ``dir()`` now omits exactly the names the missing distribution disables;
    ``__all__`` keeps all 226 and access still raises the ``ImportError``.
    """
    script = _INTROSPECTION.format(BLOCKED=blocked, ROOT=str(ROOT))
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, cwd=ROOT, check=True
    )
    output = result.stdout

    assert "ALL 226" in output, output
    assert "UNRESOLVABLE_LISTED []" in output, output
    assert "RESOLVABLE_UNLISTED []" in output, output
    assert "GETMEMBERS OK" in output, output
    assert "RENDER_DOC OK" in output, output
    if "matplotlib" in blocked:
        assert "UNRESOLVABLE 37" in output, output
    else:
        assert "UNRESOLVABLE 0" in output, output


def test_the_extra_map_matches_the_declared_optional_dependencies():
    """Every extra the resolver names must exist in pyproject."""
    from traminex import _OPTIONAL_DEPENDENCY_EXTRAS

    declared = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"][
        "optional-dependencies"
    ]
    unknown = sorted(set(_OPTIONAL_DEPENDENCY_EXTRAS.values()) - set(declared))
    assert not unknown, f"extras named by the resolver but not declared: {unknown}"


def test_an_unrelated_import_error_is_not_relabelled():
    """Only the declared optional distributions get the extras message."""
    output = _resolve_with_blocked({"this_module_does_not_exist"}, "seqplot")
    assert "RESOLVED" in output or "matplotlib" not in output, output


def test_the_sklearn_refusal_names_the_function_that_requires_it(monkeypatch):
    """``seqtree_fit`` must not blame the R-faithful ``seqtree`` (2026-09-20 T-04).

    ``seqtree_explore.seqtree_fit`` raised "scikit-learn is required for
    seqtree", naming a *different* public function -- the R-faithful weighted
    tree, which needs no scikit-learn at all (pinned by the companion test
    below).  A user without the ``extras`` profile was told the R-faithful route
    was unavailable when it is not.  The sibling
    ``disstree_explore.disstree_fit`` already names itself.
    """
    import numpy as np

    import traminex as tx
    from traminex import seqtree_explore

    # A real sequence object: the type of x is checked first (TX-AUD-20261019-036).
    sequences = tx.seqdef(np.array([["A", "B"], ["B", "A"]], dtype=object))
    # scikit-learn is imported where it is used (TX-AUD-20261011-034): block the import.
    monkeypatch.setitem(sys.modules, "sklearn", None)
    with pytest.raises(RuntimeError) as excinfo:
        seqtree_explore.seqtree_fit(sequences, [0, 1])
    assert str(excinfo.value) == "scikit-learn is required for seqtree_fit"


def test_the_r_faithful_seqtree_does_not_require_scikit_learn(monkeypatch):
    """``seqtree`` is the weighted TraMineR tree; it must run on the core deps."""
    import numpy as np

    import traminex as tx

    sequences = tx.seqdef(
        np.array([["A", "A", "B"], ["A", "B", "B"], ["B", "B", "A"], ["A", "A", "A"]], dtype=object)
    )
    predictors = pd.DataFrame({"g": ["x", "x", "y", "y"]})
    diss = tx.seqdist(sequences, method="LCS")

    monkeypatch.setitem(sys.modules, "sklearn", None)
    tree = tx.seqtree(sequences, predictors, diss=diss, R=1, pval=1.0, min_size=1)
    assert tree.root.info["n"] == 4


def test_hand_written_refusals_name_the_installable_distribution():
    """Refusals must say ``scikit-learn``, not the import name ``sklearn``.

    ``traminex/__init__.py`` records the intended wording: "The scikit-learn path
    was already exemplary -- ``RuntimeError: scikit-learn is required for
    disstree_fit`` -- and this matches that wording."  Three hand-written
    refusals used the *import* root instead, which is not what a user types to
    install it (2026-09-20 T-08).
    """
    import re

    offenders = []
    for path in sorted((ROOT / "traminex").rglob("*.py")):
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if re.search(r'"sklearn is required for', line):
                offenders.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()}")
    assert not offenders, "use the distribution name 'scikit-learn':\n" + "\n".join(offenders)


_SEABORN_CALLS = """
import sys

class _Block:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] == "seaborn":
            raise ImportError("No module named %r" % name, name="seaborn")
        return None

sys.meta_path.insert(0, _Block())
sys.path.insert(1, {ROOT!r})

import matplotlib

matplotlib.use("Agg")
import numpy as np

import traminex as tx
from traminex.plotting import transition_heatmap

x = tx.seqdef(np.array([["A", "B", "B"], ["A", "A", "B"], ["B", "B", "A"]], dtype=object))
calls = {{
    "seqfrequency_plot": lambda: tx.seqfrequency_plot(x),
    "duration_boxplot": lambda: tx.duration_boxplot(x, by_state=True),
    "chronogram": lambda: tx.chronogram(x, seaborn=True),
    "transition_heatmap": lambda: transition_heatmap(x, seaborn=True),
    "seqfrequency_plot_without_seaborn": lambda: tx.seqfrequency_plot(x, seaborn=False),
}}
for label, call in calls.items():
    try:
        call()
    except Exception as exc:
        print(label, "RAISED", type(exc).__name__, exc)
    else:
        print(label, "OK")
"""


def test_absent_seaborn_names_the_plot_extra_and_the_opt_out():
    """A seaborn-backed plot without seaborn must say how to fix it (S7-16).

    ``seqfrequency_plot`` defaults to ``seaborn=True``, and it and three other
    plots import seaborn inside the function, where the facade's extras message
    can never see the failure: a user with matplotlib but not seaborn got a bare
    ``No module named 'seaborn'`` (audit S7-16, 2026-09-22).
    """
    result = subprocess.run(
        [sys.executable, "-c", _SEABORN_CALLS.format(ROOT=str(ROOT))],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    )
    lines = dict(line.split(" ", 1) for line in result.stdout.splitlines())

    for label in ("seqfrequency_plot", "duration_boxplot", "chronogram", "transition_heatmap"):
        assert lines[label].startswith(f"RAISED ImportError seaborn is required for {label}"), lines
        assert "install the 'plot' extra (pip install 'traminex[plot]')" in lines[label], lines
        assert "seaborn=False" in lines[label], lines
    assert lines["seqfrequency_plot_without_seaborn"] == "OK", lines


def test_every_mapped_extra_is_reachable_from_a_real_failure_path():
    """Each ``_OPTIONAL_DEPENDENCY_EXTRAS`` key must be able to produce its message.

    The resolver sees only an ``ImportError`` raised while a lazy module is
    imported; scikit-learn and numba are imported inside ``try`` blocks (their
    features fail locally or fall back), so their entries could never fire
    (audit S7-16). seaborn's entry is read by the plots' own import helper.
    """
    from traminex import _OPTIONAL_DEPENDENCY_EXTRAS

    assert set(_OPTIONAL_DEPENDENCY_EXTRAS) == {"matplotlib", "seaborn"}


_IMPORT_PROBE = """
import os
import sys

before = {{key for key in os.environ if key.startswith("KMP_")}}
import {MODULE}  # noqa: E402

added = sorted({{key for key in os.environ if key.startswith("KMP_")}} - before)
print("SKLEARN", any(name == "sklearn" or name.startswith("sklearn.") for name in sys.modules))
print("KMP", added)
"""


@pytest.mark.parametrize("module", ["traminex.disstree", "traminex.seqtree", "traminex.cluster"])
def test_importing_a_module_loads_no_scikit_learn(module):
    """scikit-learn is imported where it is used (TX-AUD-20261011-034).

    Importing traminex.disstree, traminex.seqtree or traminex.cluster loaded
    scikit-learn at module level for the exploratory helpers, and scikit-learn
    sets KMP_DUPLICATE_LIB_OK and KMP_INIT_AT_FORK in the process environment.
    """
    result = subprocess.run(
        [sys.executable, "-c", _IMPORT_PROBE.format(MODULE=module)],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    )
    assert "SKLEARN False" in result.stdout, result.stdout
    assert "KMP []" in result.stdout, result.stdout


def test_the_exploratory_tree_plots_name_the_plot_extra(monkeypatch):
    """disstree_plot and seqtree_plot name the 'plot' extra without matplotlib (TX-AUD-20261011-033).

    They imported matplotlib.pyplot themselves and let a bare ModuleNotFoundError
    through, where the package facade names the extra to install.
    """
    pytest.importorskip("sklearn")
    import numpy as np

    import traminex as tx
    from traminex.disstree_explore import disstree_fit, disstree_plot
    from traminex.seqtree_explore import seqtree_fit, seqtree_plot

    diss_model = disstree_fit(np.array([[0.0, 1, 4], [1, 0, 4], [4, 4, 0]]), [0, 0, 1])
    sequences = tx.seqdef(
        np.array([list("AAB"), list("ABB"), list("BBA"), list("BAA")], dtype=object)
    )
    seq_model = seqtree_fit(sequences, np.array([0, 0, 1, 1]))
    monkeypatch.setitem(sys.modules, "matplotlib", None)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", None)
    for plot, model, name in (
        (disstree_plot, diss_model, "disstree_plot"),
        (seqtree_plot, seq_model, "seqtree_plot"),
    ):
        with pytest.raises(ImportError) as excinfo:
            plot(model)
        assert str(excinfo.value) == (
            f"matplotlib is required for {name}; install the 'plot' extra (pip install 'traminex[plot]')"
        )
