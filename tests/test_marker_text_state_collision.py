"""A real state spelled like the nonresponse marker stays real (audit S2-09).

The SPS, SPELL and TSE routes deliberately keep a literal ``"*"`` in the data as
an ordinary state; ``seqdef`` renames its nonresponse marker (``"<NR>"``) only
when missing cells exist, so an object without them keeps ``nr="*"`` beside a
real ``"*"`` state. Three consumers then matched the marker by its text:

- ``STS_to_SPELL`` / ``seqformat(to="SPELL")`` re-ran ``seqdef`` on the decoded
  text, which read the real ``"*"`` cells as nonresponse: they were dropped
  with ``with_missing=False`` and merged into the adjacent missing spell with
  ``with_missing=True``;
- ``seqstart`` looked up its fill marker by text and so turned the real ``"*"``
  cells into missing cells;
- ``seqgen_missing`` appended a second ``"*"`` to the alphabet, two codes spelled
  the same.

Source of truth: traminex's own marker model (a marker gets a collision-free
spelling when a real state already uses its text, as ``seqdef`` and
``TSE_to_STS`` do), and, for ``seqstart``, live TraMineRextras 0.6.9, whose
fill is ``NA`` -- distinct from a literal ``"*"`` -- recorded in
``tests/r_goldens/audit20260922_s2_cases.txt``
(``scripts/r_goldens/make_goldens_audit20260922_s2.R``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex import STS_to_SPELL, STSList, alphabet, seqdef, seqformat, seqgen_missing, seqstart

FIXTURE = Path(__file__).resolve().parent / "r_goldens" / "audit20260922_s2_cases.txt"


@pytest.mark.parametrize("opt", [1, 2])
def test_seqcompare_combining_keeps_undeclared_marker_shaped_real_states(opt):
    import math

    states = ["*", "%", "A"]
    data = np.array([[0, 1, 2], [0, 0, 2], [0, 2, 2], [1, 2, 0], [1, 1, 2], [1, 0, 0]])
    data.setflags(write=False)
    first = STSList(data[:3], states)
    second = STSList(data[3:], states)
    # One balanced draw includes all three cases of each group. HAM costs do
    # not depend on corpus fitting, so both optimization paths have this exact
    # hand-counted discrepancy regardless of the seeded order within a group.
    result = tx.seqCompare(
        first,
        second,
        s=3,
        opt=opt,
        method="HAM",
        norm="none",
        squared=False,
        weighted=False,
        BFopt=1,
    )
    pair_total = sum(
        int(np.count_nonzero(data[i] != data[j])) for i in range(6) for j in range(i + 1, 6)
    )
    within = (
        sum(
            int(np.count_nonzero(data[i] != data[j]))
            for group in (range(3), range(3, 6))
            for i in group
            for j in group
            if i < j
        )
        / 3
    )
    lrt = 6 * math.log((pair_total / 6) / within)
    bic = lrt - math.log(6)
    expected = [lrt, math.erfc(math.sqrt(lrt / 2)), bic, math.exp(bic / 2)]
    np.testing.assert_allclose(result.to_numpy(), [expected], rtol=0, atol=1e-14)
    assert first.nr_code is second.nr_code is first.void_code is second.void_code is None
    np.testing.assert_array_equal(first.data, data[:3])
    np.testing.assert_array_equal(second.data, data[3:])


def _r_cells(case: str) -> list[list[str]]:
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        name, field, value = line.split("\t")
        if (name, field) == (case, "cells"):
            return [row.split("|") for row in value.split(";")]
    raise AssertionError(f"no R record {case}")


def _decoded(x: STSList) -> list[list[str]]:
    return [[x._state_for_code(int(code)) for code in row] for row in x.data]


def _star_with_missing() -> STSList:
    # Real "*" twice, one missing cell, then "A": nr is renamed to "<NR>".
    return seqdef(
        pd.DataFrame([["(*,2)", "(@,1)", "(A,1)"]]),
        informat="SPS",
        missing="@",
        right=None,
    )


def _star_without_missing() -> STSList:
    # Real "*" state, no missing cell: nr stays "*" beside the real "*".
    return seqdef(
        pd.DataFrame([["(*,2)", "(A,2)"], ["(A,3)", "(*,1)"]]),
        informat="SPS",
        missing="@",
    )


def _spells(frame: pd.DataFrame) -> list[list[object]]:
    return frame[["begin", "end", "states"]].values.tolist()


def test_sources_keep_the_literal_star_as_a_real_state():
    with_missing = _star_with_missing()
    assert alphabet(with_missing) == ["*", "A"]
    assert with_missing.nr == "<NR>"
    without_missing = _star_without_missing()
    assert alphabet(without_missing) == ["*", "A"]
    assert (without_missing.nr, without_missing.nr_code) == ("*", None)


def test_sts_to_spell_keeps_real_star_spells_apart_from_missing_ones():
    x = _star_with_missing()

    with_missing = STS_to_SPELL(x, with_missing=True)
    assert _spells(with_missing) == [[1, 2, "*"], [3, 3, "<NR>"], [4, 4, "A"]]
    assert with_missing.attrs["statl"] == ["*", "A", "<NR>"]

    without_missing = STS_to_SPELL(x, with_missing=False)
    assert _spells(without_missing) == [[1, 2, "*"], [4, 4, "A"]]
    assert without_missing.attrs["statl"] == ["*", "A"]


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*'missing'\ set\ as\ .c\('<NR>','%'\).,\ the\ 'nr'\ and\ 'void'\ code\ from\ the\ 'data'\ state\ sequence\ object:UserWarning"
)
def test_seqformat_to_spell_keeps_real_star_spells():
    spells = seqformat(_star_with_missing(), from_="STS", to="SPELL")
    assert _spells(spells) == [[1, 2, "*"], [3, 3, "<NR>"], [4, 4, "A"]]


def test_sts_to_spell_on_an_object_whose_nr_shares_the_real_state_text():
    expected = [[1, 2, "*"], [3, 4, "A"], [1, 3, "A"], [4, 4, "*"]]
    for with_missing in (True, False):
        spells = STS_to_SPELL(_star_without_missing(), with_missing=with_missing)
        assert _spells(spells) == expected
        assert spells.attrs["statl"][:2] == ["*", "A"]


def test_seqstart_fill_stays_distinct_from_a_real_star_state():
    x = _star_without_missing()

    shifted = seqstart(x, data_start=[1, 2], new_start=1)

    assert alphabet(shifted) == ["*", "A"]
    assert shifted.nr == "<NR>"
    assert shifted.nr_code is not None and shifted.alphabet[shifted.nr_code] == "<NR>"
    as_r = [["NA" if cell == "<NR>" else cell for cell in row] for row in _decoded(shifted)]
    assert as_r == _r_cells("seqstart_literal_star")


def test_seqstart_without_a_collision_keeps_the_star_marker():
    x = seqdef(pd.DataFrame([["A", "B"], ["B", "A"]]))

    shifted = seqstart(x, data_start=[1, 2], new_start=1)

    assert shifted.nr == "*"
    assert _decoded(shifted) == [["A", "B", "*"], ["*", "B", "A"]]


def test_seqgen_missing_gives_the_new_marker_its_own_spelling():
    x = _star_without_missing()

    generated = seqgen_missing(x, p_cases=1.0, p_left=0.0, p_right=0.5, random_state=0)

    assert len(set(generated.alphabet)) == len(generated.alphabet)
    assert alphabet(generated) == ["*", "A"]
    assert generated.nr == "<NR>"
    assert generated.nr_code is not None and generated.alphabet[generated.nr_code] == "<NR>"
    real_star = x.alphabet.index("*")
    assert (generated.data == real_star).sum() >= 1
    assert all(
        generated.data[row, column] == real_star or x.data[row, column] != real_star
        for row in range(x.nseq)
        for column in range(x.seqlen)
        if generated.data[row, column] != generated.nr_code
    )


@pytest.mark.parametrize("marker", ["*", "%"])
def test_literal_marker_states_keep_cost_dimensions_missing_status_and_dss(marker):
    x = STSList(
        np.array([[0, 0, 1], [1, 1, 1]]),
        [marker, "A"],
        cpal=["red", "blue"],
        labels=["literal", "ordinary"],
        weights=np.array([2.0, 1.0]),
    )
    original = x.copy()
    assert tx.seqhasmiss(x) is False
    np.testing.assert_array_equal(tx.seqhasmiss(x, per_sequence=True), [False, False])
    sm = np.array([[0.0, 2.0], [2.0, 0.0]])
    np.testing.assert_array_equal(tx.seqdist(x, "OM", subm=sm, indel=1), [[0, 4], [4, 0]])
    dhd = np.repeat(sm[None], 3, axis=0)
    np.testing.assert_array_equal(tx.seqdist(x, "DHD", subm=dhd), [[0, 4], [4, 0]])
    control = STSList(x.data.copy(), ["B", "A"], weights=x.weights.copy())
    for method in ["OM", "DHD"]:
        options = {"subm": "TRATE"} if method == "OM" else {}
        np.testing.assert_array_equal(
            tx.seqdist(x, method, **options), tx.seqdist(control, method, **options)
        )
    dss = tx.seqdss(x)
    assert alphabet(dss) == [marker, "A"]
    assert dss.nr_code is None
    assert dss.to_strings()[0] == f"{marker}-A"
    assert dss.to_strings()[1] == f"A-{dss.void}"
    assert dss.void_code not in [0, 1]
    assert dss.labels[:2] == ["literal", "ordinary"]
    assert dss.cpal[:2] == ["red", "blue"]
    np.testing.assert_array_equal(dss.weights, [2, 1])
    np.testing.assert_array_equal(x.data, original.data)
    assert x.alphabet == original.alphabet
    assert x.nr_code == original.nr_code and x.void_code == original.void_code


@pytest.mark.parametrize("marker", ["*", "%"])
def test_literal_marker_plot_artists_keep_state_colors_and_full_opacity(marker):
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    from traminex.visualization.plots.distribution import state_distribution_plot
    from traminex.visualization.plots.index import sequence_index_plot

    x = STSList(np.array([[0, 1], [1, 0]]), [marker, "A"], cpal=["red", "blue"])
    fig, ax = sequence_index_plot(x)
    try:
        fig.canvas.draw()
        assert len(ax.images) == 1
        image = ax.images[0]
        np.testing.assert_array_equal(image.cmap(image.norm(0)), to_rgba("red"))
        np.testing.assert_array_equal(image.cmap(image.norm(1)), to_rgba("blue"))
    finally:
        plt.close(fig)
    fig, ax = state_distribution_plot(x)
    try:
        fig.canvas.draw()
        assert all(p.get_alpha() == 1 for p in ax.patches)
        assert [t.get_text() for t in ax.get_legend().get_texts()] == [marker, "A"]
        np.testing.assert_array_equal(ax.patches[0].get_facecolor(), to_rgba("red"))
    finally:
        plt.close(fig)


def test_sps_literal_star_survives_distance_and_dss_consumers():
    x = _star_without_missing()
    assert tx.seqhasmiss(x) is False
    assert alphabet(tx.seqdss(x)) == ["*", "A"]
    np.testing.assert_array_equal(
        tx.seqdist(x, "OM", subm=np.array([[0, 2], [2, 0]]), indel=1), [[0, 4], [4, 0]]
    )


@pytest.mark.parametrize("fn", [tx.seqLLCS, tx.seqLLCP])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("codes", [(5, 6), (2, 6), (6, 2)])
def test_scalar_common_helpers_compare_logical_missing_identity(fn, reverse, codes):
    left = STSList(
        np.array([[0, codes[0], 9]]),
        ["A", "B"] if codes[0] > 2 else ["A", "B", "*"],
        nr_code=codes[0],
        void_code=9,
    )
    right = STSList(
        np.array([[0, codes[1], 12]]),
        ["A", "B"] if codes[1] > 2 else ["A", "B", "*"],
        nr_code=codes[1],
        void_code=12,
    )
    a, b = (right, left) if reverse else (left, right)
    assert fn(a, b) == 2
    assert fn(a, a) == 2
    assert fn(np.array([0, 5]), np.array([0, 6])) == 1
    other = STSList(np.array([[1, 6]]), ["A", "B"], nr_code=6)
    assert fn(left, other) == (1 if fn is tx.seqLLCS else 0)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ max\ substitution\ cost\ less\ or\ equal\ to\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("reverse", [False, True])
def test_external_refseq_checks_real_marker_text_in_alphabet(reverse):
    with_star = STSList(np.array([[0]], dtype=np.int32), ["*", "A"])
    without_star = STSList(np.array([[0]], dtype=np.int32), ["A"])
    source, reference = (without_star, with_star) if reverse else (with_star, without_star)
    with pytest.raises(ValueError, match="same alphabet"):
        tx.seqdist(source, method="LCS", refseq=reference)
    # A real star remains distinct from A when both codebooks agree.
    ordinary_a = STSList(np.array([[1]], dtype=np.int32), ["*", "A"])
    np.testing.assert_array_equal(tx.seqdist(with_star, method="LCS", refseq=ordinary_a), [2.0])


@pytest.mark.parametrize("alphabet_order", [["%", "A", "*", "B"], ["A", "*", "B", "%"]])
@pytest.mark.parametrize("method", ["OM", "DHD"])
@pytest.mark.parametrize("explicit_void", [False, True])
def test_distance_costs_follow_real_state_order_with_interior_markers(
    alphabet_order, method, explicit_void
):
    codes = {state: index for index, state in enumerate(alphabet_order)}
    data = np.array([[codes[state] for state in row] for row in ["A*B%", "B*A%"]], dtype=np.uint8)
    data.setflags(write=False)
    source = STSList(data, alphabet_order, nr_code=codes["*"], void_code=codes["%"])
    costs = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 2.0], [2.0, 2.0, 0.0]])
    kwargs = (
        {"subm": costs, "indel": 1.0}
        if method == "OM"
        else {"subm": np.broadcast_to(costs, (4, 3, 3))}
    )
    if explicit_void:
        kwargs["void"] = codes["%"]
    np.testing.assert_array_equal(
        tx.seqdist(source, method=method, with_missing=True, **kwargs), [[0.0, 2.0], [2.0, 0.0]]
    )
    np.testing.assert_array_equal(source.data, data)
    assert source.alphabet == alphabet_order


@pytest.mark.parametrize("start", ["beg", "end"])
@pytest.mark.parametrize("inverse", [False, True])
def test_sorting_preserves_real_marker_text(start, inverse):
    data = np.array([[2], [1], [3]], dtype=np.uint8)
    data.setflags(write=False)
    source = STSList(data, ["A", "%", "*", "B"])
    expected = np.array([1, 0, 2])
    np.testing.assert_array_equal(tx.sorti(source, start=start, sort_index=not inverse), expected)
    np.testing.assert_array_equal(tx.sortv(source, start=start), expected)
    np.testing.assert_array_equal(source.data, data)
    assert source.nr_code is None and source.void_code is None


@pytest.mark.parametrize("dtype,states", [(np.int8, 128), (np.uint8, 256)])
def test_sorting_can_create_missing_rank_beyond_narrow_storage(dtype, states):
    data = np.array([[states - 1], [0], [1]], dtype=dtype)
    data.setflags(write=False)
    source = STSList(data, ["%", *[f"S{i}" for i in range(1, states)]], void_code=0)
    np.testing.assert_array_equal(tx.sorti(source), [2, 0, 1])
    np.testing.assert_array_equal(tx.sortv(source), [1, 2, 0])
    np.testing.assert_array_equal(source.data, data)


@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("include_missing", [False, True])
def test_modal_sequence_preserves_real_marker_states_and_time_metadata(weighted, include_missing):
    data = np.array([[0, 1], [0, 1]], dtype=np.uint8)
    data.setflags(write=False)
    source = STSList(
        data,
        ["*", "%", "A"],
        labels=["Star", "Percent", "A"],
        cpal=["red", "blue", "green"],
        xtlab=["Month 1", "Month 2"],
        weights=np.array([0.25, 0.75]),
        tick_last=True,
    )
    modal = tx.seqmodst(source, weighted=weighted, include_missing=include_missing)
    assert modal.to_strings() == ["*-%"]
    assert modal.nr_code is None and modal.void_code is None
    assert modal.xtlab == source.xtlab
    assert modal.labels == source.labels and modal.cpal == source.cpal
    assert modal.tick_last
    assert modal.occurrences_ == 2
    assert modal.nbseq_ == (1.0 if weighted else 2.0)
    np.testing.assert_array_equal(modal.frequencies_, [1.0, 1.0])
    np.testing.assert_array_equal(source.data, data)
    modal.cpal[0] = "black"
    assert source.cpal[0] == "red"


def test_undefined_modal_state_gets_distinct_missing_code():
    source = STSList(np.array([[0], [1]]), ["*", "A"], weights=np.array([0.0, 0.0]))
    modal = tx.seqmodst(source)
    assert modal.nr_code is not None and modal.nr_code != 0
    assert modal.alphabet[0] == "*"
    assert modal.nr != "*"
    assert tx.seqhasmiss(modal)
    assert np.isnan(modal.frequencies_).all()
    assert modal.occurrences_ == 0


@pytest.mark.parametrize("recodes", [{}, {"Other": ["A"]}])
def test_recoding_keeps_unmapped_marker_text_states_real(recodes):
    source = STSList(
        np.array([[0, 1, 2]], dtype=np.uint8),
        ["*", "%", "A"],
        xtlab=["t1", "t2", "t3"],
        weights=np.array([0.5]),
        ids=["person"],
    )
    result = tx.seqrecode(source, recodes)
    assert result.to_strings() == ["*-%-" + ("Other" if recodes else "A")]
    assert result.nr_code is None and result.void_code is None
    assert not tx.seqhasmiss(result)
    assert tx.seqlength(result).tolist() == [3]
    assert result.xtlab == source.xtlab and result.ids == source.ids
    np.testing.assert_array_equal(result.weights, [0.5])
    np.testing.assert_array_equal(source.data, [[0, 1, 2]])


@pytest.mark.parametrize("method", ["first", "last", "first.valid", "last.valid", "mostfreq"])
def test_granularity_can_add_missing_without_duplicate_real_marker_text(method):
    data = np.array([[0, 0, 1, 1], [0, 0, 0, 0]], dtype=np.uint8)
    data.setflags(write=False)
    source = STSList(
        data,
        ["*", "%"],
        void_code=1,
        labels=["Star", "Void"],
        cpal=["red", "white"],
        weights=np.array([0.25, 0.75]),
        xtlab=["t1", "t2", "t3", "t4"],
    )
    result = tx.seqgranularity(source, tspan=2, method=method)
    assert result.to_strings() == ["*-%", "*-*"]
    assert result.void_code == 1
    assert result.nr_code is not None and result.nr_code != 0
    assert result.nr != "*"
    assert len(set(result.alphabet)) == len(result.alphabet)
    assert result.xtlab == ["t1", "t3"]
    assert result.labels[:2] == source.labels and result.cpal[:2] == source.cpal
    np.testing.assert_array_equal(result.weights, source.weights)
    np.testing.assert_array_equal(source.data, data)


@pytest.mark.parametrize("replace", [False, True])
def test_polyad_resampling_preserves_real_marker_text(replace):
    data = np.array([[0, 1, 2, 0], [1, 2, 0, 1]], dtype=np.uint8)
    data.setflags(write=False)
    roles = [STSList(data, ["*", "%", "A"]) for _ in range(2)]
    iterations, seed = 8, 123
    rng = np.random.default_rng(seed)
    sampled = rng.choice(2, size=(iterations, 2), replace=True, p=[0.5, 0.5])
    if replace:
        positions = rng.integers(0, 4, size=(iterations, 2, 4))
    else:
        positions = np.array([[rng.permutation(4) for _ in range(2)] for _ in range(iterations)])
    resampled = np.take_along_axis(data[sampled], positions, axis=2)
    expected = np.count_nonzero(resampled[:, 0] != resampled[:, 1], axis=1)
    result = tx.seqpolyads(
        roles,
        a=2,
        method="HAM",
        norm="none",
        T=iterations,
        s=seed,
        replace=replace,
        weighted=False,
    )
    np.testing.assert_array_equal(result["observed.dist"], [0.0, 0.0])
    np.testing.assert_array_equal(result["random.dist"], expected)
    np.testing.assert_array_equal(data, [[0, 1, 2, 0], [1, 2, 0, 1]])


@pytest.mark.parametrize("include_missing", [False, True])
def test_integrative_potential_addresses_real_star_by_code(include_missing):
    data = np.array([[0, 1, 0]], dtype=np.uint8)
    data.setflags(write=False)
    source = STSList(data, ["*", "A"])
    # Position weights are 1, 2, 3: '*' occupies weights 1+3 out of 6.
    np.testing.assert_allclose(
        tx.seqintegr(source, state="*", include_missing=include_missing),
        [2.0 / 3.0],
        rtol=0,
        atol=1e-15,
    )
    np.testing.assert_array_equal(source.data, data)
