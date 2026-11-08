"""Batch-13 plotting items against R (``scripts/r_goldens/make_goldens_plot_items.R``).

The golden records the arguments TraMineR's plotting code passes to the graphics functions (traced in live R), so
each test compares the drawn geometry with R's: row heights, axis limits, tick positions, colours and labels.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.collections import QuadMesh  # noqa: E402

import traminex as tx  # noqa: E402

GOLDEN_PATH = Path(__file__).resolve().parent / "r_goldens" / "plot_items.csv"
with open(GOLDEN_PATH, encoding="utf-8", newline="") as handle:
    GOLDEN = {(row["case"], row["field"]): row["value"] for row in csv.DictReader(handle)}


def _numbers(case: str, field: str) -> np.ndarray:
    return np.array([np.nan if v == "NA" else float(v) for v in GOLDEN[(case, field)].split()])


def _texts(case: str, field: str) -> list[str]:
    return GOLDEN[(case, field)].split(" ;; ")


# TX-AUD-20261019-048: R's index plot draws each row with a height equal to the case's weight (plot.stslist.R:149,
# barplot(width = weights[idxs[sel]]) with space = 0); with weighted = FALSE or all weights 1 the rows are equal.
def _index_sequences(weights) -> tx.STSList:
    cells = np.array([list(row) for row in _texts("index", "cells")], dtype=object)
    return tx.seqdef(cells, alphabet=["A", "B"], weights=weights)


@pytest.mark.parametrize("case", ["index_weighted", "index_weighted_sorted"])
def test_weighted_index_rows_have_r_s_heights(case: str) -> None:
    sortv = _numbers("index", "sortv") if case.endswith("sorted") else None
    fig, ax = tx.sequence_index_plot(_index_sequences(_numbers("index", "weights")), sortv=sortv)
    meshes = [c for c in ax.collections if isinstance(c, QuadMesh)]
    assert len(meshes) == 1
    edges = meshes[0].get_coordinates()[:, 0, 1]
    assert edges[0] == 0.0
    np.testing.assert_allclose(np.diff(edges), _numbers(case, "widths"), rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        ax.get_yticks(), _numbers(case, "ylabel_positions"), rtol=0, atol=1e-12
    )
    plt.close(fig)


@pytest.mark.parametrize("case", ["index_unweighted", "index_all_one"])
def test_unweighted_index_rows_are_equal(case: str) -> None:
    assert set(_numbers(case, "widths")) == {1.0}
    weights = np.ones(6) if case == "index_all_one" else _numbers("index", "weights")
    fig, ax = tx.sequence_index_plot(_index_sequences(weights), weighted=case != "index_unweighted")
    assert not [c for c in ax.collections if isinstance(c, QuadMesh)]
    assert ax.images[0].get_array().shape[0] == 6
    plt.close(fig)


# TX-AUD-20261019-049: R's plot.subseqelistchisq gives every panel one frequency range, ylim = c(min(min(freq), 0),
# max(freq)) over all groups' frequency columns (seqecmpgroup.R:136-138); the golden's groups have different largest
# frequencies, so a per-panel range differs from R's. traminex draws the bars horizontally: R's y range is the x range.
def test_chisq_panels_share_r_s_frequency_range() -> None:
    eseq = tx.seqecreate(
        id=_numbers("chisq", "tse_id").astype(int),
        timestamp=_numbers("chisq", "tse_time"),
        event=_texts("chisq", "tse_event"),
    )
    comparison = tx.seqecmpgroup(
        tx.seqefsub(eseq, pmin_support=0.05), group=_texts("chisq", "group")
    )
    assert list(comparison.labels) == _texts("chisq", "labels")
    freq = comparison.data[[f"Freq.{g}" for g in comparison.labels]].to_numpy(dtype=float)
    np.testing.assert_allclose(
        freq.max(axis=0), _numbers("chisq", "max_freq_per_group"), rtol=0, atol=1e-12
    )
    fig, axes = tx.plot_subseqelistchisq(comparison)
    assert len(axes) == 2
    for k, ax in enumerate(axes, start=1):
        np.testing.assert_allclose(
            ax.get_xlim(), _numbers("chisq", f"ylim_panel{k}"), rtol=0, atol=1e-12
        )
    plt.close(fig)


# TX-AUD-20261019-050: R colours seqedplot's event k with cpal[k] (seqedplot.R:192-227), so a palette shorter than the
# events gives the later curves colour NA, drawn in no colour (recorded in the golden). DIVERGENCE: traminex recycles
# the palette over the events, as it does for a short state palette (docs/PLOTTING.md), so every event keeps a
# visible colour: event k takes cpal[k mod len(cpal)], in every panel and in the legend.
def _seqed_sequences():
    eseq = tx.seqecreate(
        id=_numbers("seqed", "tse_id").astype(int),
        timestamp=_numbers("seqed", "tse_time"),
        event=_texts("seqed", "tse_event"),
    )
    eseq.lengths = _numbers("seqed", "lengths")
    return eseq


@pytest.mark.parametrize("case", ["seqed_short", "seqed_short_grouped"])
def test_seqedplot_recycles_a_short_palette(case: str) -> None:
    from matplotlib.colors import to_rgba

    cpal = _texts(case, "cpal")
    events = _texts(case, "legend_text")
    r_colours = _texts(case, "curve_colours")
    npanel = len(r_colours) // len(events)
    # R's own value: the first event takes the colour, the second none.
    assert r_colours == ["red", "NA"] * npanel
    expected = [cpal[k % len(cpal)] for k in range(len(events))]
    group = _texts("seqed", "group") if case.endswith("grouped") else None
    before = plt.get_fignums()
    fig, axes = tx.seqedplot(_seqed_sequences(), cpal=cpal, group=group)
    assert len(axes) == npanel
    for ax in axes:
        assert [line.get_label() for line in ax.lines] == events
        assert [to_rgba(line.get_color()) for line in ax.lines] == [to_rgba(c) for c in expected]
        for line, r_colour in zip(ax.lines, r_colours):
            if r_colour != "NA":
                assert to_rgba(line.get_color()) == to_rgba(r_colour)
    legend = fig.legends[0]
    assert [text.get_text() for text in legend.get_texts()] == events
    # matplotlib 3.6.0, the declared floor, names them legendHandles (TX-REP-20261025-002).
    handles = legend.legend_handles if hasattr(legend, "legend_handles") else legend.legendHandles
    assert [to_rgba(h.get_color()) for h in handles] == [to_rgba(c) for c in expected]
    plt.close(fig)
    assert plt.get_fignums() == before


def test_seqedplot_refuses_an_empty_palette() -> None:
    before = plt.get_fignums()
    with pytest.raises(ValueError, match="cpal must contain at least one colour"):
        tx.seqedplot(_seqed_sequences(), cpal=[])
    assert plt.get_fignums() == before


# TX-AUD-20261019-051: a state whose spells are all censored has no event; seqsurvplot draws it as a flat line at 1
# up to its last follow-up time, as R's plot.stslist.surv draws one curve per state (TX-AUD-20261010-042). The
# docstring said such a state gets no curve.
def test_seqsurvplot_draws_a_state_without_events_flat() -> None:
    cells = np.array([list("AABCC"), list("ABBCC"), list("BBAAC")], dtype=object)
    fig, ax = tx.seqsurvplot(tx.seqdef(cells, alphabet=["A", "B", "C"]))
    curves = {line.get_label(): line for line in ax.lines}
    assert list(curves) == ["A", "B", "C"]
    np.testing.assert_array_equal(curves["C"].get_xdata(), [0.0, 2.0])
    np.testing.assert_array_equal(curves["C"].get_ydata(), [1.0, 1.0])
    plt.close(fig)
    doc = " ".join(tx.seqsurvplot.__doc__.split())
    assert "only states with at least one event get a curve" not in doc
    assert "a state without events is a flat line at 1" in doc


# TX-AUD-20261019-052: a chronogram's rows are seqstatd's states (R's seqstatd.R:3, 12): the real states, plus the
# missing state when include_missing (a zero row when the object has none, R's "phantom" nr); the void marker never.
_CHRONO_CELLS = np.array(
    [["A", "B", None, "C"], ["A", None, "B", "B"], ["C", "C", "A", None]], dtype=object
)


def _chrono_expected(x: tx.STSList, include_missing: bool) -> tuple[list[str], np.ndarray]:
    from traminex.sequences import _alphabet_codes

    freq = np.asarray(tx.seqstatd(x, include_missing=include_missing)["Frequencies"], dtype=float)
    rows = [
        freq[:, code] if code is not None else np.zeros(freq.shape[0])
        for code in _alphabet_codes(x, with_missing=include_missing)
    ]
    return tx.alphabet(x, with_missing=include_missing), np.vstack(rows)


def _cluster(x: tx.STSList, rows: np.ndarray) -> tx.STSList:
    from traminex.group_compare import _subset_sequence

    return _subset_sequence(x, rows)


def _row_names(ax) -> list[str]:
    return [label.get_text() for label in ax.get_yticklabels()]


@pytest.mark.parametrize("include_missing", [False, True])
@pytest.mark.parametrize("route", ["image", "seaborn"])
def test_chronogram_rows_are_seqstatd_states(include_missing: bool, route: str) -> None:
    x = tx.seqdef(_CHRONO_CELLS, alphabet=["A", "B", "C"])
    names, values = _chrono_expected(x, include_missing)
    assert names == (["A", "B", "C", "*"] if include_missing else ["A", "B", "C"])
    fig, ax = tx.chronogram(x, include_missing=include_missing, seaborn=route == "seaborn")
    assert _row_names(ax) == names
    drawn = ax.images[0].get_array() if route == "image" else ax.collections[0].get_array()
    np.testing.assert_allclose(
        np.asarray(drawn, dtype=float).reshape(values.shape), values, rtol=0, atol=1e-12
    )
    plt.close(fig)


def test_chronogram_includes_r_s_phantom_missing_row() -> None:
    x = tx.seqdef(np.array([list("ABC"), list("BCA")], dtype=object), alphabet=["A", "B", "C"])
    fig, ax = tx.chronogram(x, include_missing=True)
    assert _row_names(ax) == ["A", "B", "C", "*"]
    np.testing.assert_array_equal(ax.images[0].get_array()[3], np.zeros(3))
    plt.close(fig)


@pytest.mark.parametrize("include_missing", [False, True])
def test_cluster_chronograms_rows_are_seqstatd_states(include_missing: bool) -> None:
    x = tx.seqdef(_CHRONO_CELLS, alphabet=["A", "B", "C"])
    labels = np.array([1, 1, 2])
    fig, axes = tx.cluster_chronograms(x, labels, include_missing=include_missing)
    for ax, cluster in zip(axes.ravel(), (1, 2)):
        names, values = _chrono_expected(
            _cluster(x, np.flatnonzero(labels == cluster)), include_missing
        )
        assert _row_names(ax) == names
        np.testing.assert_allclose(ax.images[0].get_array(), values, rtol=0, atol=1e-12)
    plt.close(fig)


@pytest.mark.parametrize("include_missing", [False, True])
def test_cluster_chronograms_with_representatives_rows_and_markers(include_missing: bool) -> None:
    x = tx.seqdef(_CHRONO_CELLS, alphabet=["A", "B", "C"])
    labels = np.array([1, 1, 2])
    D = np.asarray(tx.seqdist(x, method="LCS", with_missing=True), dtype=float)
    fig, axes = tx.cluster_chronograms_with_representatives(
        x, labels, D, include_missing=include_missing
    )
    for ax, cluster in zip(axes.ravel(), (1, 2)):
        rows = np.flatnonzero(labels == cluster)
        names, values = _chrono_expected(_cluster(x, rows), include_missing)
        assert _row_names(ax) == names
        np.testing.assert_allclose(ax.images[0].get_array(), values, rtol=0, atol=1e-12)
        # The medoid's states are marked on their rows; a state without a row is not marked.
        medoid = int(rows[tx.medoid_index(D[np.ix_(rows, rows)])])
        states = [x.alphabet[code] for code in x.data[medoid]]
        expected_y = [names.index(s) if s in names else np.nan for s in states]
        np.testing.assert_array_equal(ax.lines[0].get_ydata(), expected_y)
    plt.close(fig)


# TX-AUD-20261011-029: the grouped helpers (Python additions) refuse labels that are not one per sequence and an
# order that is not a permutation of the rows, with a ValueError and no figure left open.
_GROUPED = np.array([list("ABAB"), list("BBAA"), list("AABB"), list("BABA")], dtype=object)


def _refused(call, match: str) -> None:
    before = plt.get_fignums()
    with pytest.raises(ValueError, match=match):
        call()
    assert plt.get_fignums() == before


@pytest.mark.parametrize("labels", [[1, 1, 2], [1, 1, 2, 2, 3]], ids=["shorter", "longer"])
@pytest.mark.parametrize(
    "helper", ["cluster_chronograms", "cluster_chronograms_with_representatives"]
)
def test_cluster_chronograms_refuse_labels_of_the_wrong_length(helper: str, labels) -> None:
    x = tx.seqdef(_GROUPED, alphabet=["A", "B"])
    args = (x, np.array(labels))
    if helper.endswith("representatives"):
        args += (np.zeros((4, 4)),)
    _refused(lambda: getattr(tx, helper)(*args), "labels length must equal number of sequences")


@pytest.mark.parametrize(
    "order",
    [[0, 0, 1, 2], [0, 1, 2], [0, 1, 2, 3, 0], [1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0]],
    ids=["repeated", "shorter", "longer", "out-of-range", "not-integer"],
)
def test_sequence_index_plot_grouped_refuses_an_order_that_is_not_a_permutation(order) -> None:
    x = tx.seqdef(_GROUPED, alphabet=["A", "B"])
    _refused(
        lambda: tx.sequence_index_plot_grouped(x, np.array([1, 2, 1, 2]), order=order),
        "order must be a permutation of the sequence rows",
    )


def test_sequence_index_plot_grouped_accepts_a_permutation() -> None:
    x = tx.seqdef(_GROUPED, alphabet=["A", "B"])
    fig, ax, order = tx.sequence_index_plot_grouped(x, np.array([1, 2, 1, 2]), order=[3, 1, 2, 0])
    np.testing.assert_array_equal(order, [3, 1, 2, 0])
    plt.close(fig)


# TX-AUD-20261011-031: R always shows a deprecation warning; Python shows a DeprecationWarning by default only when it
# is attributed to the caller's code (__main__ or a test file), not to a line inside traminex. Each deprecated
# argument or function, called from this file, is attributed to this file (traminex/_deprecated.py
# _caller_stacklevel), also when reached through another traminex function (seqtabstocc's tlim, forwarded to seqtab).
def _deprecation_calls():
    x = tx.seqdef(
        np.array([list("AABB"), list("ABAB"), list("BBAA"), list("BABA")], dtype=object),
        alphabet=["A", "B"],
    )
    group = np.array([1, 1, 2, 2])

    def leaf_summary():
        result = tx.seqtree_fit(x, group, max_depth=2)
        return tx.seqtree_leaf_summary(x, result)

    def eordplot():
        events = tx.seqecreate(
            id=[1, 1, 2, 2, 3], timestamp=[0, 2, 0, 1, 3], event=["A", "B", "B", "A", "C"]
        )
        fig, _ = tx.seqeordplot(events)
        plt.close(fig)

    return {
        "disscenter_double_center": (
            lambda: tx.disscenter_double_center(np.array([[0.0, 1.0], [1.0, 0.0]])),
            "deprecated alias",
        ),
        "seqdist_method": (lambda: tx.seqdist(x, method="LCSopt"), "'LCSopt' is deprecated"),
        "seqdist_norm": (lambda: tx.seqdist(x, method="LCS", norm=True), "norm=True"),
        "seqdiff_seqdist_arg": (
            lambda: tx.seqdiff(x, group, seqdist_arg={"method": "LCS"}, weighted=False),
            "seqdist_arg is deprecated",
        ),
        "seqformat_compressed": (
            lambda: tx.seqformat(
                np.array([["A", "B"]], dtype=object), from_="STS", to="STS", compressed=False
            ),
            "compressed is deprecated",
        ),
        "seqformat_nr": (
            lambda: tx.seqformat(
                np.array([["A", "B"]], dtype=object), from_="STS", to="STS", nr="*"
            ),
            "nr is deprecated",
        ),
        "seqtree_leaf_summary": (leaf_summary, "pass only result"),
        "seqtab_tlim": (lambda: tx.seqtab(x, tlim=1), "tlim is deprecated"),
        "seqtabstocc_tlim": (lambda: tx.seqtabstocc(x, tlim=1), "tlim is deprecated"),
        "seqipos_integration": (
            lambda: tx.seqipos(x, pos_states="A", index="integration"),
            "'integration' obsolete",
        ),
        "seqtrate_statl": (lambda: tx.seqtrate(x, statl=["A", "B"]), "statl is deprecated"),
        "seqeordplot": (eordplot, "superseded by seqpcplot"),
    }


@pytest.mark.parametrize("site", sorted(_deprecation_calls()))
def test_deprecation_warnings_name_the_caller(site: str) -> None:
    if site == "seqtree_leaf_summary":  # only seqtree needs the optional scikit-learn
        pytest.importorskip("sklearn")
    call, match = _deprecation_calls()[site]
    with pytest.warns(DeprecationWarning, match=match) as records:
        call()
    matching = [
        r for r in records if issubclass(r.category, DeprecationWarning) and match in str(r.message)
    ]
    assert matching
    assert {Path(r.filename).resolve() for r in matching} == {Path(__file__).resolve()}


# TX-AUD-20261011-067 (a): R's frequency plot labels its y axis 0% at the bottom of the first bar and the cumulative
# percentage of the shown sequences at the top of the last, titled "Cum. % freq. (n=N)" or "(weighted n=N)"
# (plot.stslist.freq.R:56-61, 102-110). traminex stacks the bands without gaps from 0 to the shown proportion.
@pytest.mark.parametrize("case", ["freq_unweighted", "freq_weighted", "freq_top2"])
def test_seqfplot_y_axis_is_r_s(case: str) -> None:
    cells = np.array([list(row) for row in _texts("freq", "cells")], dtype=object)
    weights = _numbers("freq", "weights") if case == "freq_weighted" else None
    x = tx.seqdef(cells, alphabet=["A", "B", "C"], weights=weights)
    widths = _numbers(case, "bar_widths")
    centres = _numbers(case, "bar_centres")
    # In R the two labels sit on the outer edges of the bars.
    np.testing.assert_allclose(
        _numbers(case, "ylabel_positions"),
        [centres[0] - widths[0] / 2, centres[-1] + widths[-1] / 2],
        rtol=0,
        atol=1e-9,
    )
    fig, ax = tx.seqfplot(x, top=len(_numbers(case, "idxs")))
    np.testing.assert_allclose(ax.get_yticks(), [0.0, widths.sum() / 100], rtol=0, atol=1e-12)
    assert [label.get_text() for label in ax.get_yticklabels()] == _texts(case, "ylabels")
    assert ax.get_ylabel() == GOLDEN[(case, "ylab")]
    plt.close(fig)


# TX-AUD-20261011-067 (b): R puts a time tick at seq(1, L, by = xtstep), plus L when tick.last (plot.stslist.R:160-161),
# each given its column name, and draws the labels left to right, omitting one that would come closer than an "m" to
# the last drawn (?axis; review I-1). traminex keeps every tick and applies that rule on the drawn figure. An
# index-plot cell t (0-based) is centred on t, R's on t + 0.5.
@pytest.mark.parametrize(
    "case", ["ticks_step6_72", "ticks_default_72", "ticks_step5_last_72", "ticks_step3_last_5"]
)
def test_time_ticks_follow_r_s_xtstep_grid(case: str) -> None:
    length = int(_numbers(case, "length")[0])
    xtstep = int(_numbers(case, "xtstep")[0])
    cells = np.array(
        [[("A", "B")[(row + t) % 2] for t in range(length)] for row in range(3)], dtype=object
    )
    x = tx.seqdef(
        cells,
        alphabet=["A", "B"],
        cnames=[f"c{t}" for t in range(1, length + 1)],
        xtstep=xtstep,
        tick_last=GOLDEN[(case, "tick_last")] == "TRUE",
    )
    fig, ax = tx.sequence_index_plot(x)
    positions = _numbers(case, "positions")
    np.testing.assert_allclose(ax.get_xticks(), positions - 0.5, rtol=0, atol=1e-12)
    r_labels = dict(zip(positions - 0.5, _texts(case, "labels")))
    shown = _shown_time_labels(fig, ax)
    assert shown[0][1] == _texts(case, "labels")[0]
    assert all(text == r_labels[position] for _, text, position in shown)
    gap = _m_width(fig, ax)
    for (left, _, _), (right, _, _) in zip(shown, shown[1:]):
        assert right.x0 - left.x1 >= gap - 1.0
    if case in ("ticks_step6_72", "ticks_step3_last_5"):  # these labels fit at the default size
        assert [text for _, text, _ in shown] == _texts(case, "labels")
    plt.close(fig)


# TX-AUD-20261011-067 (h): R puts the tick of position t at t - 0.5, the centre of a cell drawn on [t - 1, t]. These
# four plots draw position t (0-based) on [t, t + 1], so the tick of R's t lands on R's own coordinate, t - 0.5.
def _centre_axes(case: str):
    cells = np.array([list(row) for row in _texts("freq", "cells")], dtype=object)
    x = tx.seqdef(cells, alphabet=["A", "B", "C"])
    if case == "centres_seqfplot":
        fig, ax = tx.seqfplot(x)
    elif case == "centres_seqmsplot":
        fig, ax = tx.seqmsplot(x)
    elif case == "centres_seqrplot":
        fig, ax = tx.seqrplot(x, method="HAM")
    else:
        rows = _numbers("centres", "rf_rows").astype(int) - 1
        x20 = tx.seqdef(cells[rows], alphabet=["A", "B", "C"])
        fig, axes = tx.seqrfplot(x20, tx.seqdist(x20, method="HAM"))
        ax = axes[0]
    return fig, ax


@pytest.mark.parametrize(
    "case", ["centres_seqfplot", "centres_seqmsplot", "centres_seqrplot", "centres_seqrfplot"]
)
def test_time_ticks_sit_at_cell_centres(case: str) -> None:
    from matplotlib.patches import Rectangle

    fig, ax = _centre_axes(case)
    first = min(
        (p for p in ax.patches if isinstance(p, Rectangle) and p.get_width() == 1.0),
        key=lambda p: (p.get_x(), p.get_y()),
    )
    assert (first.get_x(), first.get_x() + first.get_width()) == (0.0, 1.0)
    np.testing.assert_allclose(ax.get_xticks(), _numbers(case, "positions"), rtol=0, atol=1e-12)
    plt.close(fig)


# TX-AUD-20261011-067 (c, d): R's seqplot takes ylab, with.legend and cpal for every type (seqplot.R:6-11), drops
# cpal for Ht (seqplot.R:299-301) and draws one state legend by default, none for Ht (seqplot.R:210). The golden
# records, per type, the y label R draws, the colours it hands to barplot and the number of state legends.
_OPTION_TYPES = ["i", "I", "d", "dH", "f", "ms", "mt", "r", "rf", "Ht"]


def _options_input(type_: str):
    cells = np.array([list(row) for row in _texts("freq", "cells")], dtype=object)
    if type_ == "rf":
        cells = cells[_numbers("centres", "rf_rows").astype(int) - 1]
    x = tx.seqdef(cells, alphabet=["A", "B", "C"])
    extra: dict = {}
    if type_ == "r":
        extra["method"] = "HAM"
    if type_ == "rf":
        extra["diss"] = tx.seqdist(x, method="HAM")
    return x, extra


def _legend_count(fig) -> int:
    return len(fig.legends) + sum(1 for ax in fig.axes if ax.get_legend() is not None)


def _drawn_colours(fig) -> set:
    from matplotlib.colors import to_rgba

    colours: set = set()
    for ax in fig.axes:
        for image in ax.images:
            colours.update(tuple(c) for c in getattr(image.cmap, "colors", []))
        for collection in ax.collections:
            colours.update(tuple(c) for c in collection.get_facecolors())
        for patch in ax.patches:
            colours.add(tuple(patch.get_facecolor()))
    return {tuple(round(v, 6) for v in to_rgba(c)) for c in colours}


def _rgba_set(colours) -> set:
    from matplotlib.colors import to_rgba

    return {tuple(round(v, 6) for v in to_rgba(c)) for c in colours}


def _main_axes(type_: str, axes) -> list:
    axes = list(axes) if isinstance(axes, (list, tuple)) else [axes]
    return axes[:1] if type_ == "rf" else axes  # rf: the medoid panel (R gives the other no ylab)


@pytest.mark.parametrize("type_", _OPTION_TYPES)
def test_seqplot_takes_r_s_common_options(type_: str) -> None:
    from traminex.visualization.palette import ensure_cmap

    case = f"options_{type_}"
    x, extra = _options_input(type_)
    cpal = _texts("options", "cpal")
    if type_ != "Ht":
        extra["cpal"] = cpal
    fig, axes = tx.seqplot(x, type=type_, ylab="my y", with_legend=False, **extra)
    assert [ax.get_ylabel() for ax in _main_axes(type_, axes)] == [GOLDEN[(case, "ylab")]] * len(
        _main_axes(type_, axes)
    )
    assert _legend_count(fig) == int(_numbers(case, "state_legends")[0])
    if GOLDEN[(case, "colours")]:
        drawn = _drawn_colours(fig)
        assert drawn & _rgba_set(_texts(case, "colours"))
        assert not drawn & _rgba_set(ensure_cmap(x, None).colors[:3])
    plt.close(fig)


@pytest.mark.parametrize("type_", _OPTION_TYPES)
def test_seqplot_draws_r_s_default_state_legend(type_: str) -> None:
    x, extra = _options_input(type_)
    fig, _ = tx.seqplot(x, type=type_, **extra)
    assert _legend_count(fig) == int(_numbers(f"options_{type_}", "state_legends_default")[0])
    plt.close(fig)


@pytest.mark.parametrize("type_", ["i", "d", "f", "ms"])
def test_grouped_seqplot_takes_the_common_options_in_every_panel(type_: str) -> None:
    x, extra = _options_input(type_)
    group = np.array([1, 1, 1, 2, 2, 2])
    cpal = _texts("options", "cpal")
    fig, axes = tx.seqplot(
        x, type=type_, group=group, ylab="my y", with_legend=False, cpal=cpal, **extra
    )
    assert [ax.get_ylabel() for ax in axes] == ["my y", "my y"]
    assert _legend_count(fig) == 0
    assert _drawn_colours(fig) & _rgba_set(cpal)
    plt.close(fig)
    fig, axes = tx.seqplot(x, type=type_, group=group, **extra)
    assert _legend_count(fig) == 1
    plt.close(fig)


# TX-AUD-20261011-067 (e): R's plot.stslist.surv draws no legend; traminex's seqsurvplot legend (an addition) names
# the states as seqsplot's legend does, by their labels (seqsplot.R:91-92), the missing state included.
@pytest.mark.parametrize("with_missing", [False, True])
def test_seqsurvplot_legend_names_states_by_label(with_missing: bool) -> None:
    cells = np.array(
        [list("AABCC"), list("ABBCC"), list("BBAAC"), ["A", "A", None, "C", "C"]], dtype=object
    )
    x = tx.seqdef(cells, alphabet=["A", "B", "C"], labels=["alpha", "beta", "gamma"])
    fig_s, _ = tx.seqsplot(x, with_missing=with_missing)
    expected = [t.get_text() for legend in fig_s.legends for t in legend.get_texts()]
    expected += [
        t.get_text() for ax in fig_s.axes if ax.get_legend() for t in ax.get_legend().get_texts()
    ]
    plt.close(fig_s)
    assert expected[:3] == ["alpha", "beta", "gamma"]
    fig, ax = tx.seqsurvplot(x, with_missing=with_missing)
    assert [t.get_text() for t in ax.get_legend().get_texts()] == expected
    plt.close(fig)


# TX-AUD-20261011-067 (f): R draws missing cells in missing.color, opaque (col2rgb alpha 255), in the index and the
# distribution plots; the golden records the colour barplot receives for the missing state.
def _missing_sequences() -> tx.STSList:
    cells = np.array(
        [[None if v == "NA" else v for v in row.split()] for row in _texts("missing", "cells")],
        dtype=object,
    )
    return tx.seqdef(cells, alphabet=["A", "B"], right=None)


def _r_missing_rgba() -> tuple:
    from matplotlib.colors import to_rgba

    assert _numbers("missing", "missing_color_alpha")[0] == 255
    return to_rgba(GOLDEN[("missing", "missing_color")])


def test_index_plot_draws_missing_cells_solid() -> None:
    from matplotlib.colors import to_rgba

    x = _missing_sequences()
    assert _texts("missing_index", "colours")[-1] == GOLDEN[("missing", "missing_color")]
    fig, ax = tx.sequence_index_plot(x)
    base = ax.images[0]
    assert to_rgba(base.cmap(base.norm(x.nr_code))) == _r_missing_rgba()
    missing = x.data == x.nr_code
    assert missing.any()
    for overlay in ax.images[1:]:
        assert np.all(np.asarray(overlay.get_array())[missing][:, 3] == 0.0)
    plt.close(fig)


def test_distribution_plot_draws_missing_bars_solid() -> None:
    from matplotlib.colors import to_rgba

    x = _missing_sequences()
    assert _texts("missing_distribution", "colours")[-1] == GOLDEN[("missing", "missing_color")]
    fig, ax = tx.state_distribution_plot(x, include_missing=True)
    bars = [c for c in ax.containers if c.get_label() == "missing"]
    assert len(bars) == 1
    for patch in bars[0].patches:
        assert to_rgba(patch.get_facecolor()) == _r_missing_rgba()
    plt.close(fig)


# TX-REP-20261023-001: R's representative plot titles its y axis "N representative(s) (n=N)" ("weighted n=N"), N
# rounded to 2 decimals, and labels none of the bands: its one side-2 axis labels the statistics overlay's rows B and
# A (plot.stslist.rep.R:109-110, 132-136, 232-236), an overlay traminex does not draw.
@pytest.mark.parametrize("case", ["rep_unweighted", "rep_weighted"])
def test_seqrplot_y_axis_is_r_s(case: str) -> None:
    assert _texts(case, "side2_labels") == ["B", "A"]
    cells = np.array([list(row) for row in _texts("freq", "cells")], dtype=object)
    weights = _numbers("freq", "weights") if case == "rep_weighted" else None
    x = tx.seqdef(cells, alphabet=["A", "B", "C"], weights=weights)
    fig, ax = tx.seqrplot(x, method="HAM")
    assert ax.get_ylabel() == GOLDEN[(case, "ylab")]
    assert list(ax.get_yticks()) == []
    plt.close(fig)


# Batch-13 review I-1: R draws the time labels left to right and omits one that would come closer than the width of
# an "m" to the last label drawn (?axis, gap.axis = 1 for labels parallel to the axis), so many ticks never give
# overlapping labels. Checked on the drawn figure at the default size.
def _shown_time_labels(fig, ax) -> list:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    shown = [t for t in ax.get_xticklabels() if t.get_visible() and t.get_text()]
    return sorted(
        ((t.get_window_extent(renderer), t.get_text(), t.get_position()[0]) for t in shown),
        key=lambda item: item[0].x0,
    )


def _m_width(fig, ax) -> float:
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import text_to_path

    size = ax.get_xticklabels()[0].get_fontsize()
    width = text_to_path.get_text_width_height_descent("m", FontProperties(size=size), False)[0]
    return width * fig.dpi / 72.0


@pytest.mark.parametrize("length,xtstep", [(72, 2), (200, 2), (200, 5), (40, 2), (72, 1)])
def test_time_labels_never_overlap(length: int, xtstep: int) -> None:
    cells = np.array(
        [[("A", "B")[(row + t) % 2] for t in range(length)] for row in range(3)], dtype=object
    )
    names = [f"c{t}" for t in range(1, length + 1)]
    x = tx.seqdef(cells, alphabet=["A", "B"], cnames=names, xtstep=xtstep)
    fig, ax = tx.sequence_index_plot(x)
    shown = _shown_time_labels(fig, ax)
    gap = _m_width(fig, ax)
    assert shown[0][1] == "c1"
    for extent, text, position in shown:
        assert text == names[int(round(position))]
    for (left, _, _), (right, _, _) in zip(shown, shown[1:]):
        assert right.x0 - left.x1 >= gap - 1.0
    plt.close(fig)


# Batch-13 review I-3: apply_time_axis is consumed by weightedcluster (docs/CONSUMED_SURFACE.md), so its behaviour
# change is recorded there (prompt section 3).
def test_consumed_surface_records_the_time_axis_change() -> None:
    doc = (Path(__file__).resolve().parents[1] / "docs" / "CONSUMED_SURFACE.md").read_text(
        encoding="utf-8"
    )
    changes = " ".join(doc.split("## Behaviour changes to listed objects", 1)[1].split())
    assert "apply_time_axis" in changes and "TX-AUD-20261011-067" in changes
    assert "fuzzy.py" in changes


# Batch-13 review I-2: seqplot's legend rule acts only on the legends its own call drew; it never removes a legend
# the caller made, and a second seqplot in the same figure gets its own state legend.
def _option_sequences() -> tx.STSList:
    cells = np.array([list(row) for row in _texts("freq", "cells")], dtype=object)
    return tx.seqdef(cells, alphabet=["A", "B", "C"])


def test_seqplot_leaves_the_caller_s_legends_alone() -> None:
    fig, (mine, target) = plt.subplots(1, 2)
    mine.plot([0, 1], [0, 1], label="mine")
    own = mine.legend()
    figure_legend = fig.legend(handles=mine.lines, labels=["mine"])
    tx.seqplot(_option_sequences(), type="d", ax=target, with_legend=False)
    assert mine.get_legend() is own and figure_legend in fig.legends
    assert target.get_legend() is None
    plt.close(fig)


def test_two_seqplots_in_one_figure_each_get_a_state_legend() -> None:
    fig, axes = plt.subplots(1, 2)
    for ax in axes:
        tx.seqplot(_option_sequences(), type="i", ax=ax)
    assert [ax.get_legend() is not None for ax in axes] == [True, True]
    plt.close(fig)


# R's plot.SAMM is seqplot(..., group = origin), which draws one state legend for the figure (seqsamm.R:56-60).
@pytest.mark.parametrize("type_", ["i", "d"])
def test_plot_samm_draws_one_state_legend(type_: str) -> None:
    from traminex.sequences import STSList

    x = STSList(
        data=np.array([[0, 1, 0, 0, 2, 2, 0], [0, 0, 1, 1, 0, 2, 2]], dtype=int),
        alphabet=["A", "B", "C"],
        nr_code=None,
        void_code=None,
    )
    samm = tx.seqsamm(x, 2)
    fig, _ = tx.plot_SAMM(samm, type=type_)
    assert _legend_count(fig) == 1
    plt.close(fig)
    fig, _ = tx.plot_SAMM(samm, type=type_, with_legend=False)
    assert _legend_count(fig) == 0
    plt.close(fig)


# Batch-13 review M-2: a renderer that hands main or group to seqplot keeps its own legend default (docs/PLOTTING.md).
@pytest.mark.parametrize(
    "renderer", ["seqfplot", "seqmsplot", "seqmtplot", "seqrplot", "seqdHplot"]
)
def test_main_and_group_do_not_change_a_renderer_s_legend(renderer: str) -> None:
    x = _option_sequences()
    extra = {"method": "HAM"} if renderer == "seqrplot" else {}
    fig, _ = getattr(tx, renderer)(x, **extra)
    direct = _legend_count(fig)
    plt.close(fig)
    for options in ({"main": "T"}, {"group": np.array([1, 1, 1, 2, 2, 2])}):
        fig, _ = getattr(tx, renderer)(x, **options, **extra)
        assert _legend_count(fig) == min(direct, 1)
        plt.close(fig)


# Batch-13 review M-1: R hands each group panel its own ylab[np] for every type (seqplot.R:173-176, 283), rf included.
def test_grouped_rf_takes_one_ylab_per_level() -> None:
    x, extra = _options_input("rf")
    group = np.repeat([1, 2], x.nseq // 2)
    fig, axes = tx.seqplot(x, type="rf", group=group, ylab=["first", "second"], **extra)
    assert [ax.get_ylabel() for ax in axes] == ["first", "second"]
    plt.close(fig)
    with pytest.raises(ValueError, match="one value per group level"):
        tx.seqplot(x, type="rf", group=group, ylab=["only one"], **extra)


# Batch-13 review M-3: R plots the first n subsequences as plot(x[1:n]), and `[` keeps those rows of data
# (subseqelist.R:28-32), so its uniform range comes from the shown rows; traminex's top=n is that subset.
def test_chisq_top_takes_the_range_of_the_shown_rows() -> None:
    eseq = tx.seqecreate(
        id=_numbers("chisq", "tse_id").astype(int),
        timestamp=_numbers("chisq", "tse_time"),
        event=_texts("chisq", "tse_event"),
    )
    comparison = tx.seqecmpgroup(
        tx.seqefsub(eseq, pmin_support=0.05), group=_texts("chisq", "group")
    )
    freq = comparison.data[[f"Freq.{g}" for g in comparison.labels]].to_numpy(dtype=float)[:3]
    expected = (min(float(freq.min()), 0.0), float(freq.max()))
    assert expected != tuple(
        _numbers("chisq", "ylim_panel1")
    )  # the shown rows' range is not all rows'
    fig, axes = tx.plot_subseqelistchisq(comparison, top=3)
    for ax in axes:
        np.testing.assert_allclose(ax.get_xlim(), expected, rtol=0, atol=1e-12)
    plt.close(fig)


# TX-REP-20261023-002 (batch-13 review M-4): R's cpal = "red" is a one-element vector; seqedplot read a string as
# its characters. One colour, recycled over the events (TX-AUD-20261019-050).
def test_seqedplot_reads_a_string_palette_as_one_colour() -> None:
    from matplotlib.colors import to_rgba

    assert GOLDEN[("seqed_short", "cpal")] == "red"
    fig, axes = tx.seqedplot(_seqed_sequences(), cpal="red")
    assert [to_rgba(line.get_color()) for line in axes[0].lines] == [to_rgba("red")] * 2
    plt.close(fig)


# Batch-13 review M-6: the weighted (mesh) index plot draws missing cells as solidly as the image route, and its
# hatching covers each cell's own row, from the cumulated weights below it to those above.
def test_weighted_index_plot_draws_missing_cells_solid_on_their_rows() -> None:
    from matplotlib.colors import to_rgba
    from matplotlib.patches import Rectangle

    x = _missing_sequences()
    x.weights = np.array([2.0, 3.0])
    fig, ax = tx.sequence_index_plot(x, apply_hatching=True)
    meshes = [c for c in ax.collections if isinstance(c, QuadMesh)]
    base, overlays = meshes[0], meshes[1:]
    assert to_rgba(base.cmap(base.norm(x.nr_code))) == _r_missing_rgba()
    assert overlays and all(overlay.get_alpha() == 0.0 for overlay in overlays)
    hatched = [p for p in ax.patches if isinstance(p, Rectangle) and p.get_hatch()]
    row, column = (int(v) for v in np.argwhere(x.data == x.nr_code)[0])
    bottom = float(np.sum(x.weights[:row]))
    assert [(p.get_x(), p.get_y(), p.get_width(), p.get_height()) for p in hatched] == [
        (column - 0.5, bottom, 1.0, float(x.weights[row]))
    ]
    plt.close(fig)
