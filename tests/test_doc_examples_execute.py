"""Documented Python examples must run (audit S7-04, 2026-09-22).

What was wrong: ``docs/NEW_DISTANCE_METHODS.md`` showed
``seqdist(..., method="OMstran", otto=0.5, transindel="constant")`` and
``docs/USAGE_new_apis.md`` showed ``seqrep(X, criterion="density",
coverage=0.5)`` and ``tx.seqdist(X, method="OM")`` without a substitution
cost. Since ``seqdist`` requires ``subm`` for the OM family exactly where
TraMineR requires ``sm`` (``docs/DISTANCES.md``; ``seqdist.R`` stops with
"'sm' must be specified"), all three raised ``ValueError``. Nothing executed
the examples: the docs gate checks call signatures only, and Quarto does not
run plain Markdown ``python`` fences.

Source of truth: the documented ``seqdist`` contract above and R's
``seqrep``, which computes ``seqdist(seqdata, ...)`` from the forwarded
arguments when ``diss`` is not given (``seqrep.R``).

Each listed document's fences run cumulatively in one namespace, because
later fences reuse names defined by earlier ones. ``docs/TraMineR_examples.qmd``
is listed too: its examples are plain ``python`` fences, which Quarto displays
but does not run, so nothing else executed them (audit S7-21). ``docs/DISTANCES.md`` and
``CONTRIBUTING.md`` are not listed: they contain deliberate fragments (a
function-template example and a sketch over undefined ``corpus``/``a``/``b``).
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

_FENCE = re.compile(r"^([ \t]*)```(?:python|\{python[^}\n]*\})[^\n]*\n(.*?)^[ \t]*```", re.S | re.M)

RUNNABLE_DOCS = (
    "README.md",
    "docs/API_CONSISTENCY.md",
    "docs/DISTANCES_NORMALIZATION.md",
    "docs/NEW_DISTANCE_METHODS.md",
    "docs/OMLOC.md",
    "docs/USAGE_new_apis.md",
    "docs/TraMineR_examples.qmd",
)


def _fences(text: str) -> list[tuple[int, str]]:
    fences = []
    for match in _FENCE.finditer(text):
        indent, code = match.group(1), match.group(2)
        if indent:
            code = "\n".join(
                line[len(indent) :] if line.startswith(indent) else line
                for line in code.splitlines()
            )
        fences.append((text.count("\n", 0, match.start(2)) + 1, code))
    return fences


def _compile_fence(code: str, filename: str):
    """Compile one documented fence as a reader's own script would be compiled.

    ``dont_inherit=True``: without it the fence inherits this module's
    ``from __future__ import annotations`` and no annotation is evaluated
    (TX-AUD-20261010-063).
    """
    return compile(code, filename, "exec", dont_inherit=True)


@pytest.mark.parametrize("relative", RUNNABLE_DOCS)
def test_every_python_fence_in_the_document_runs(relative: str, tmp_path, monkeypatch) -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    monkeypatch.chdir(tmp_path)  # any file a snippet writes lands in tmp_path
    fences = _fences((ROOT / relative).read_text(encoding="utf-8"))
    assert fences, f"{relative} has no python fence"
    namespace: dict[str, object] = {"__name__": "__doc_example__"}
    for line, code in fences:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                exec(_compile_fence(code, f"{relative}:{line}"), namespace)
            except Exception as exc:  # noqa: BLE001 - report which fence failed
                pytest.fail(f"{relative} fence at line {line} raised {type(exc).__name__}: {exc}")


def test_a_fence_annotation_with_an_undefined_name_fails() -> None:
    """A fence is compiled without this module's ``from __future__ import annotations``.

    ``compile`` inherits the calling module's future flags unless told not to, so the
    gate evaluated no annotation: ``x: undefined_name = 1`` passed here and raises
    NameError in a reader's script (TX-AUD-20261010-063).
    """
    compiled = _compile_fence("x: undefined_name_in_annotation = 1\n", "fence.md:1")

    with pytest.raises(NameError, match="undefined_name_in_annotation"):
        exec(compiled, {"__name__": "__doc_example__"})
