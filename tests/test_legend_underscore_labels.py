"""Automatic legends keep labels that start with "_" (TX-REP-20261017-003).

Matplotlib's automatic legend leaves out an artist whose label starts with "_",
so a state or group named "_x" vanished from the legends of the distribution
plot, seqdiffplot and seqsurvplot (R's legends show every state), and an event
named "_x" from seqedplot's legend (batch-4 review M-4). The legends now
receive their handles and labels explicitly, as seqHtplot's do since batch 2.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

SEQS = tx.seqdef(
    np.array([list("_A_A"), list("A_AA"), list("AA__"), list("_AAA")], dtype=object),
    alphabet=["_", "A"],
)


def _legend_texts(ax) -> list[str]:
    legend = ax.get_legend()
    return [] if legend is None else [text.get_text() for text in legend.get_texts()]


def test_distribution_plot_legend_keeps_an_underscore_state() -> None:
    fig, ax = tx.state_distribution_plot(SEQS)
    assert "_" in _legend_texts(ax)
    plt.close(fig)


def test_seqsurvplot_legend_keeps_an_underscore_state() -> None:
    fig, ax = tx.seqsurvplot(SEQS)
    assert "_" in _legend_texts(ax)
    plt.close(fig)


def test_seqdiffplot_legend_keeps_an_underscore_group() -> None:
    fig, ax = tx.seqdiffplot(SEQS, ["_g", "_g", "h", "h"], stat="discrepancy", weighted=False)
    assert "_g" in _legend_texts(ax)
    plt.close(fig)


@pytest.mark.parametrize("kind", ["survival", "hazard"])
def test_seqedplot_legend_keeps_an_underscore_event(kind: str) -> None:
    events = tx.seqecreate(
        id=[1, 1, 2, 2, 3], timestamp=[1, 2, 1, 3, 2], event=["_x", "b", "b", "_x", "b"]
    )
    events.lengths = [4, 4, 4]
    fig, _ = tx.seqedplot(events, type=kind)
    texts = [text.get_text() for legend in fig.legends for text in legend.get_texts()]
    assert texts == ["_x", "b"]
    plt.close(fig)


def test_every_legend_goes_through_the_label_keeping_helper() -> None:
    """Explicit handles alone do not keep "_x" at the declared floor (TX-AUD-20261019-060).

    matplotlib 3.6.0, the declared floor, drops a label starting with "_" even when the
    handles and labels are given, so the tests above failed in the floors profile while
    passing on the newer matplotlib of the ordinary suite. Every legend is therefore made by
    ``visualization.utils._legend_with_labels``, which this census checks in the source:
    no ``.legend(...)`` call outside that helper.
    """
    import ast
    from pathlib import Path

    package = Path(tx.__file__).resolve().parent
    outside: list[str] = []
    for path in sorted(package.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        helpers = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name == "_legend_with_labels"
        ]
        inside = {id(call) for helper in helpers for call in ast.walk(helper)}
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "legend"
                and id(node) not in inside
            ):
                outside.append(f"{path.relative_to(package)}:{node.lineno}")
    assert outside == []
