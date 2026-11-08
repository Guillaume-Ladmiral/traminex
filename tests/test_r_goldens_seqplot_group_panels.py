"""seqplot draws one panel per group level for every type (TX-AUD-20261010-044).

R's seqplot takes ``group`` for every type and draws one panel per level,
titled by the levels (``main="auto"``), by "<main> - <level>" (one ``main``) or
by ``main`` itself (one title per level), seqplot.R:163-205. traminex's
dispatcher handed ``group`` to renderers that do not take it (TypeError). The
levels are factor()'s: R's internal group() pastes a space before each level and
orders numeric groups as text, a quirk traminex corrects as for seqHtplot
(TX-REP-20261017-002). Fixture: ``scripts/r_goldens/make_goldens_seqplot_group_titles.R``.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqplot_group_titles.csv",
    dtype=str,
    keep_default_na=False,
)
STRINGS = ["AABBC", "ABBBC", "BBAAC", "AABBB", "BBBAA"]
SEQS = tx.seqdef(np.array([list(s) for s in STRINGS], dtype=object), alphabet=["A", "B", "C"])
TEXT = ["b", "a", "b", "a", "b"]
CASES = {
    "d_auto": ("d", TEXT, "auto"),
    "i_auto": ("i", TEXT, "auto"),
    "f_auto": ("f", TEXT, "auto"),
    "mt_auto": ("mt", TEXT, "auto"),
    "ms_auto": ("ms", TEXT, "auto"),
    "d_one_main": ("d", TEXT, "X"),
    "d_main_per_level": ("d", TEXT, ["P", "Q"]),
    "d_numeric": ("d", [10, 2, 10, 2, 10], "auto"),
    "d_missing": ("d", ["b", None, "b", "a", "b"], "auto"),
    "d_one_sequence": ("d", ["b", "a", "b", "b", "b"], "auto"),
}


def _r(case: str, column: str) -> list[str]:
    rows = GOLDEN[GOLDEN["case"] == case]
    return rows.sort_values("position", key=lambda c: c.astype(int))[column].tolist()


def _signature(ax) -> tuple:
    patches = sorted(
        (
            round(p.get_x(), 9),
            round(p.get_y(), 9),
            round(p.get_width(), 9),
            round(p.get_height(), 9),
        )
        for p in ax.patches
    )
    images = [np.asarray(image.get_array()).tolist() for image in ax.get_images()]
    return patches, images


@pytest.mark.parametrize("case", sorted(CASES))
def test_one_panel_per_level_titled_as_r(case: str) -> None:
    kind, group, main = CASES[case]
    fig, axes = tx.seqplot(SEQS, type=kind, group=group, main=main)
    assert [ax.get_title() for ax in axes] == _r(case, "factor_title")
    plt.close(fig)


@pytest.mark.parametrize("case", sorted(CASES))
def test_each_panel_is_the_plot_of_its_subset(case: str) -> None:
    kind, group, main = CASES[case]
    fig, axes = tx.seqplot(SEQS, type=kind, group=group, main=main)
    for ax, level in zip(axes, _r(case, "factor_level")):
        rows = [i for i, value in enumerate(group) if value is not None and str(value) == level]
        subset = tx.seqdef(
            np.array([list(STRINGS[i]) for i in rows], dtype=object), alphabet=["A", "B", "C"]
        )
        single_fig, single_ax = tx.seqplot(subset, type=kind)
        assert _signature(ax) == _signature(single_ax)
        plt.close(single_fig)
    plt.close(fig)


def test_r_titles_carry_its_group_quirk() -> None:
    assert _r("d_auto", "r_title") == [" a", " b"]
    assert _r("d_numeric", "r_title") == [" 10", " 2"]
    assert _r("d_numeric", "factor_title") == ["2", "10"]


@pytest.mark.parametrize(
    "function, kind",
    [(tx.seqfplot, "f"), (tx.seqmtplot, "mt"), (tx.seqmsplot, "ms"), (tx.seqdHplot, "dH")],
)
def test_the_per_type_functions_take_group(function, kind) -> None:
    fig, axes = function(SEQS, group=TEXT)
    reference, reference_axes = tx.seqplot(SEQS, type=kind, group=TEXT)
    assert [ax.get_title() for ax in axes] == [ax.get_title() for ax in reference_axes]
    assert [_signature(ax) for ax in axes] == [_signature(ax) for ax in reference_axes]
    plt.close(fig)
    plt.close(reference)


def test_main_needs_one_title_or_one_per_level() -> None:
    with pytest.raises(ValueError, match="main"):
        tx.seqplot(SEQS, type="d", group=TEXT, main=["P", "Q", "R"])


def test_group_needs_one_value_per_sequence() -> None:
    with pytest.raises(ValueError, match="group"):
        tx.seqplot(SEQS, type="d", group=["a", "b"])


@pytest.mark.parametrize("kind, name", [("r", "seqrep"), ("rf", "seqrf")])
def test_dissimilarities_are_subset_per_panel(kind: str, name: str, monkeypatch) -> None:
    # R gives each panel diss[gindex, gindex] (seqplot.R:339, :369). The values R
    # computes once on the whole set (dmax, the MDS sort variable) are
    # R-golden tested in tests/test_r_goldens_seqplot_group_rep.py.
    rows = [list("AABB"), list("ABAB"), list("BBAA"), list("AAAA"), list("ABBB"), list("BABA")] * 2
    seqs = tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B"])
    diss = np.asarray(tx.seqdist(seqs, method="LCS"), dtype=float)
    group = ["g", "h"] * 6
    module = importlib.import_module("traminex.visualization.plots.representatives")
    original = getattr(module, name)
    received = []

    def spy(*args, **kwargs):
        received.append(np.asarray(kwargs["diss"] if "diss" in kwargs else args[1], dtype=float))
        return original(*args, **kwargs)

    monkeypatch.setattr(module, name, spy)
    extra = {"k": 2} if kind == "rf" else {}
    fig, _ = tx.seqplot(seqs, type=kind, group=group, diss=diss, **extra)
    plt.close(fig)
    assert len(received) == 2
    for matrix, level in zip(received, ["g", "h"]):
        index = [i for i, value in enumerate(group) if value == level]
        np.testing.assert_array_equal(matrix, diss[np.ix_(index, index)])


DISS5 = np.asarray(tx.seqdist(SEQS, method="LCS"), dtype=float)


# Level "a" has two sequences, so k = 2 leaves no residual degree of freedom; R's
# dissrf warns the same way.
@pytest.mark.filterwarnings(r"ignore:\s*\[!!\] Negative  or zero df:UserWarning")
@pytest.mark.parametrize(
    "function, kind, extra", [(tx.seqrplot, "r", {}), (tx.seqrfplot, "rf", {"k": 2})]
)
def test_the_representative_functions_take_group(function, kind, extra) -> None:
    # R draws the medoid panel of each level for rf (plot.seqrf.R:4; seqplot.R:342);
    # seqrfplot's two-panel default applies without group only (batch-4 review I-2).
    fig, axes = function(SEQS, diss=DISS5, group=TEXT, **extra)
    reference, reference_axes = tx.seqplot(SEQS, type=kind, group=TEXT, diss=DISS5, **extra)
    assert [ax.get_title() for ax in axes] == _r("d_auto", "factor_title")
    assert [_signature(ax) for ax in axes] == [_signature(ax) for ax in reference_axes]
    plt.close(fig)
    plt.close(reference)


def test_which_plot_both_with_group_is_refused_as_in_r() -> None:
    rep = pd.read_csv(
        Path(__file__).resolve().parent / "r_goldens" / "seqplot_group_rep.csv", dtype=str
    )
    assert rep[rep["case"] == "rf_both_with_group"]["value"].tolist() == [
        'which.plot="both" applies only when group=NULL'
    ]
    with pytest.raises(ValueError, match='which_plot="both" applies only when group is None'):
        tx.seqplot(SEQS, type="rf", group=TEXT, diss=DISS5, k=2, which_plot="both")
    with pytest.raises(ValueError, match='which_plot="both" applies only when group is None'):
        tx.seqrfplot(SEQS, DISS5, group=TEXT, k=2, which_plot="both")


@pytest.mark.parametrize(
    "function, extra",
    [
        (tx.seqfplot, {}),
        (tx.seqmtplot, {}),
        (tx.seqmsplot, {}),
        (tx.seqdHplot, {}),
        (tx.seqrplot, {"diss": DISS5}),
    ],
)
def test_the_per_type_functions_take_main(function, extra) -> None:
    # R's per-type functions take main (seqfplot.R etc.: main = "auto"), batch-4 review M-2.
    fig, axes = function(SEQS, group=TEXT, main="X", **extra)
    assert [ax.get_title() for ax in axes] == _r("d_one_main", "factor_title")
    plt.close(fig)
    fig, ax = function(SEQS, main="X", **extra)
    assert ax.get_title() == "X"
    plt.close(fig)


# Level "a" has two sequences, so k = 2 leaves no residual degree of freedom; R's
# dissrf warns the same way.
@pytest.mark.filterwarnings(r"ignore:\s*\[!!\] Negative  or zero df:UserWarning")
def test_per_sequence_weights_are_split_per_panel() -> None:
    # batch-4 review M-3: a weights vector is one value per sequence, like diss and sortv.
    weights = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    fig, axes = tx.seqdHplot(SEQS, group=TEXT, weights=weights)
    for ax, level in zip(axes, _r("d_auto", "factor_level")):
        rows = [i for i, value in enumerate(TEXT) if value == level]
        subset = tx.seqdef(
            np.array([list(STRINGS[i]) for i in rows], dtype=object), alphabet=["A", "B", "C"]
        )
        single_fig, single_ax = tx.seqdHplot(subset, weights=weights[rows])
        assert _signature(ax) == _signature(single_ax)
        plt.close(single_fig)
    plt.close(fig)
    rf_fig, _ = tx.seqrfplot(SEQS, DISS5, group=TEXT, k=2, weights=weights)
    plt.close(rf_fig)
