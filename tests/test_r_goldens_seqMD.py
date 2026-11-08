"""R-golden tests for the CAT combined-OM multidomain seqdistmc / seqMD (D2).

Anchors traminex's faithful CAT engine to TraMineR 2.2-14
``seqdistmc`` / ``seqMD(what="diss")`` (see scripts/r_goldens/make_goldens_seqMD.R).
The CAT ("additive trick") builds MD sequences of combined domain states on the
expanded alphabet, sums the per-domain channel-weighted substitution/indel costs,
then runs one OM on the combined sequences.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import assert_matches_golden

from traminex import seqcost, seqdef, seqdistmc, seqMD

R2_CASES = Path(__file__).with_name("r_goldens") / "seqMD_r2_cases.csv"
R2_ERRORS = Path(__file__).with_name("r_goldens") / "seqMD_r2_errors.txt"

CH1 = [list("AABBCC"), list("ABBCCA"), list("CCBBAA"), list("BBBCCC")]
CH2 = [list("XXYYYZ"), list("YYZZXX"), list("ZZZXXY"), list("XYYYZZ")]


def _channels():
    ch1 = seqdef(CH1, alphabet=["A", "B", "C"])
    ch2 = seqdef(CH2, alphabet=["X", "Y", "Z"])
    return [ch1, ch2]


def test_seqdistmc_om_sum_matches_r():
    D = seqdistmc(_channels(), method="OM", sm="CONSTANT", indel=1, cweight=[1, 1])
    assert_matches_golden(D, "seqMD_OM_sum")


def test_seqdistmc_om_cweighted_matches_r():
    D = seqdistmc(_channels(), method="OM", sm="CONSTANT", indel=1, cweight=[2, 1])
    assert_matches_golden(D, "seqMD_OM_cw")


def test_seqdistmc_lcs_matches_r():
    D = seqdistmc(_channels(), method="LCS", cweight=[1, 1])
    assert_matches_golden(D, "seqMD_LCS")


def test_seqMD_what_diss_is_seqdistmc_alias():
    D = seqMD(_channels(), method="OM", sm="CONSTANT", indel=1, what="diss", cweight=[1, 1])
    assert_matches_golden(D, "seqMD_alias")


# R2 declared-state and DHD contract ------------------------------------------
@lru_cache(maxsize=1)
def _r2_rows():
    with R2_CASES.open(encoding="utf-8", newline="") as stream:
        return tuple(csv.DictReader(stream))


def _case_rows(case):
    return [row for row in _r2_rows() if row["case"] == case]


def _section(rows, name):
    return sorted(
        (row for row in rows if row["section"] == name),
        key=lambda row: (int(row["time"]), int(row["row"]), int(row["col"])),
    )


def _values(rows, name):
    return np.asarray([float(row["value"]) for row in _section(rows, name)], dtype=float)


def _texts(rows, name):
    return [row["text"] for row in _section(rows, name)]


def _matrix(rows, name):
    selected = _section(rows, name)
    nrow = max(int(row["row"]) for row in selected)
    ncol = max(int(row["col"]) for row in selected)
    out = np.empty((nrow, ncol), dtype=float)
    for row in selected:
        out[int(row["row"]) - 1, int(row["col"]) - 1] = float(row["value"])
    return out


def _cat_cost(rows):
    dims = tuple(_values(rows, "cat_dim").astype(int))
    out = np.empty(dims, dtype=float)
    for row in _section(rows, "cat_sm"):
        index = (int(row["row"]) - 1, int(row["col"]) - 1)
        if len(dims) == 3:
            index += (int(row["time"]) - 1,)
        out[index] = float(row["value"])
    return out


def _fixed_channels():
    a = seqdef([list("ABAB"), list("BABA"), list("AABB")], alphabet=["A", "B", "U"])
    b = seqdef([list("XYXY"), list("YXYX"), list("XXYY")], alphabet=["X", "Y", "V"])
    c = seqdef([list("PQQP"), list("QPPQ"), list("PPQQ")], alphabet=["P", "Q", "R"])
    return a, b, c


def _missing_fixed_channels():
    missing = seqdef(
        [["A", None, "B", "A"], list("BABA"), list("ABAB")],
        alphabet=["A", "B", "U"],
    )
    _, b, _ = _fixed_channels()
    return missing, b


def _ragged_channels():
    missing = seqdef(
        [["A", None, "B", None], list("BABA"), ["A", "B", None, None]],
        alphabet=["A", "B", "U"],
    )
    other = seqdef(
        [["X", "Y", None, None], list("YXYX"), ["X", "X", "Y", None]],
        alphabet=["X", "Y", "V"],
    )
    return missing, other


def _make_dhd_sm(k, width, scale):
    out = np.zeros((k, k, width), dtype=float)
    for time in range(width):
        for left in range(k):
            for right in range(k):
                out[left, right, time] = scale * abs(left - right) * (time + 1)
    return out


def _decode_md(md):
    return np.asarray([[md.alphabet[int(code)] for code in row] for row in md.data], dtype=object)


def _assert_mdseq(case, md):
    rows = _case_rows(case)
    assert list(md.alphabet[: len(_texts(rows, "md_alphabet"))]) == _texts(rows, "md_alphabet")
    expected_rows = _section(rows, "md_state")
    expected = np.empty(md.data.shape, dtype=object)
    for row in expected_rows:
        expected[int(row["row"]) - 1, int(row["col"]) - 1] = row["text"]
    np.testing.assert_array_equal(_decode_md(md), expected)


def _assert_channel_dims(case, channels, sm, *, time_varying):
    rows = _case_rows(case)
    expected = _section(rows, "channel_cost_dim")
    by_channel = {}
    for row in expected:
        by_channel.setdefault(int(row["row"]), []).append(int(float(row["value"])))
    for index, (channel, spec) in enumerate(zip(channels, sm), start=1):
        if isinstance(spec, str):
            actual = np.asarray(
                seqcost(
                    channel,
                    spec,
                    with_missing=bool(
                        channel.nr_code is not None and np.any(channel.data == channel.nr_code)
                    ),
                    time_varying=time_varying,
                )["sm"]
            )
        else:
            actual = np.asarray(spec)
        assert actual.shape == tuple(by_channel[index])


def _assert_cost(case, actual):
    rows = _case_rows(case)
    assert actual["sm"].shape == tuple(_values(rows, "cat_dim").astype(int))
    np.testing.assert_allclose(actual["sm"], _cat_cost(rows), rtol=0, atol=1e-12)
    assert actual["alphabet"] == _texts(rows, "md_alphabet")
    indel_rows = _section(rows, "cat_indel")
    if indel_rows[0]["text"] == "NULL":
        assert actual["indel"] is None
    else:
        expected_indel = _values(rows, "cat_indel")
        np.testing.assert_allclose(
            np.atleast_1d(actual["indel"]), expected_indel, rtol=0, atol=1e-12
        )
    np.testing.assert_allclose(actual["cweight"], _values(rows, "cweight"), rtol=0, atol=0)


def test_seqMD_ragged_mdseq_matches_observed_r_alphabet_and_states():
    with pytest.warns(UserWarning, match="different length"):
        actual = seqMD(_ragged_channels(), what="MDseq")
    _assert_mdseq("mdseq_ragged_2", actual)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ at\ least\ one\ indel\ cost\ does\ not\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize(
    "case,channel_count,method,indel,cweight,link",
    [
        ("cost_trate_declared_2", 2, "TRATE", [1, 2], [2, 1], "sum"),
        (
            "cost_indels_vector_mean",
            2,
            "EXPLICIT_INDELS",
            [[1, 2, 3], [4, 5, 6]],
            [2, 1],
            "mean",
        ),
        ("cost_indelslog_auto", 2, "INDELSLOG", "auto", [1, 1], "sum"),
        ("cost_trate_declared_3", 3, "TRATE", 1.5, [1, 2, 3], "sum"),
    ],
)
def test_seqMD_stationary_costs_match_declared_state_r_contract(
    case, channel_count, method, indel, cweight, link
):
    channels = list(_fixed_channels()[:channel_count])
    if method == "EXPLICIT_INDELS":
        sm = [np.asarray(seqcost(channel, "INDELS")["sm"]) for channel in channels]
    else:
        sm = [method] * channel_count
    _assert_channel_dims(case, channels, sm, time_varying=False)
    actual = seqMD(channels, what="cost", sm=sm, indel=indel, cweight=cweight, link=link)
    _assert_cost(case, actual)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\ ,\ B\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ At\ pos\ \d+,\ \d+\ substitution\ cost\ value\(s\)\ is\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize(
    "case,fixture,sm,cweight,link",
    [
        ("dhd_trate_cost_missing", "missing", "TRATE", [1, 2], "sum"),
        ("dhd_indels_cost_3", "three", "INDELS", [1, 0.5, 2], "sum"),
        ("dhd_indelslog_cost_mean", "two", "INDELSLOG", [2, 1], "mean"),
        ("dhd_explicit_cost_missing", "explicit_missing", None, [2, 1], "sum"),
    ],
)
def test_seqMD_dhd_cat_cost_arrays_match_every_r_time_slice(case, fixture, sm, cweight, link):
    if fixture in {"missing", "explicit_missing"}:
        channels = list(_missing_fixed_channels())
        with_missing = [True, False]
    else:
        count = 3 if fixture == "three" else 2
        channels = list(_fixed_channels()[:count])
        with_missing = None
    if fixture == "explicit_missing":
        sm_specs = [_make_dhd_sm(4, 4, 0.5), _make_dhd_sm(3, 4, 0.25)]
        indel = [1, 2]
    else:
        sm_specs = [
            np.asarray(
                seqcost(
                    channel,
                    sm,
                    with_missing=bool(
                        channel.nr_code is not None and np.any(channel.data == channel.nr_code)
                    ),
                    time_varying=True,
                )["sm"]
            )
            for channel in channels
        ]
        indel = [1, 2] if fixture == "missing" else [1] * len(channels)
    _assert_channel_dims(case, channels, sm_specs, time_varying=True)
    actual = seqMD(
        channels,
        method="DHD",
        what="cost",
        sm=sm_specs,
        indel=indel,
        with_missing=with_missing,
        cweight=cweight,
        link=link,
    )
    _assert_cost(case, actual)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\ ,\ B\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\+X\ ,\ A\+Y\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ X\ ,\ Y\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ At\ pos\ \d+,\ \d+\ substitution\ cost\ value\(s\)\ are\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ At\ pos\ \d+,\ \d+\ substitution\ cost\ value\(s\)\ is\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_seqMD_dhd_full_and_condensed_distances_match_r():
    missing_channels = list(_missing_fixed_channels())
    full = seqMD(
        missing_channels,
        method="DHD",
        what="diss",
        sm="TRATE",
        with_missing=[True, False],
        full_matrix=True,
    )
    np.testing.assert_allclose(
        full, _matrix(_case_rows("dhd_trate_diss_full"), "distance"), rtol=0, atol=1e-12
    )

    channels = list(_fixed_channels())
    explicit = [
        _make_dhd_sm(3, 4, 0.5),
        _make_dhd_sm(3, 4, 0.25),
        _make_dhd_sm(3, 4, 0.75),
    ]
    condensed = seqMD(
        channels,
        method="DHD",
        what="diss",
        sm=explicit,
        indel=1,
        cweight=[1, 2, 1],
        full_matrix=False,
    )
    np.testing.assert_allclose(
        condensed,
        _values(_case_rows("dhd_explicit_diss_condensed_3"), "distance"),
        rtol=0,
        atol=1e-12,
    )

    indels = seqMD(
        channels,
        method="DHD",
        what="diss",
        sm="INDELS",
        indel=1,
        cweight=[1, 0.5, 2],
        full_matrix=True,
    )
    np.testing.assert_allclose(
        indels,
        _matrix(_case_rows("dhd_indels_diss_full_3"), "distance"),
        rtol=0,
        atol=1e-12,
    )

    indelslog = seqMD(
        channels[:2],
        method="DHD",
        what="diss",
        sm="INDELSLOG",
        indel=1,
        cweight=[2, 1],
        link="mean",
        full_matrix=False,
    )
    np.testing.assert_allclose(
        indelslog,
        _values(_case_rows("dhd_indelslog_diss_condensed"), "distance"),
        rtol=0,
        atol=1e-12,
    )


def test_seqMD_warns_for_channel_specific_missing_declaration_mismatches():
    expected = " ".join(_texts(_case_rows("warning_with_missing"), "warning"))
    with pytest.warns(UserWarning) as caught:
        md = seqMD(_missing_fixed_channels(), what="MDseq", with_missing=[False, True])
    messages = " ".join(str(item.message) for item in caught)
    assert "domain 1" in expected and "domain 1" in messages
    assert "domain 2" in expected and "domain 2" in messages
    _assert_mdseq("warning_with_missing", md)


@lru_cache(maxsize=1)
def _r2_errors():
    rows = R2_ERRORS.read_text(encoding="utf-8").splitlines()[1:]
    return dict(row.split("=", 1) for row in rows)


@pytest.mark.parametrize(
    "case",
    ["dhd_state_width", "dhd_time_width", "stationary_state_width", "indel_state_width"],
)
def test_seqMD_source_width_errors_match_r_categories(case):
    a, b, _ = _fixed_channels()
    if case == "dhd_state_width":
        channels = list(_missing_fixed_channels())
        kwargs = dict(
            method="DHD",
            sm=[_make_dhd_sm(3, 4, 1), _make_dhd_sm(3, 4, 1)],
            with_missing=[True, False],
        )
    elif case == "dhd_time_width":
        channels = [a, b]
        kwargs = dict(method="DHD", sm=[_make_dhd_sm(3, 3, 1), _make_dhd_sm(3, 3, 1)])
    elif case == "stationary_state_width":
        channels = [a, b]
        kwargs = dict(sm=[np.ones((2, 2)), np.ones((3, 3))], indel=1)
    else:
        channels = [a, b]
        kwargs = dict(sm="TRATE", indel=[[1, 2], [1, 2, 3]])
    expected = _r2_errors()[case].strip().removeprefix("[!]").strip()
    with pytest.raises(ValueError) as caught:
        seqMD(channels, what="cost", **kwargs)
    assert expected in str(caught.value)
