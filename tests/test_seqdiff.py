"""Behavioural tests for the R-faithful seqdiff (RP4 T16).

The exact numerics are pinned against live TraMineR 2.2-14 in
``tests/test_r_goldens_seqdiff.py``; this module covers the contract / structure
of the sliding-window discrepancy analysis.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from traminex import STSList, seqdef, seqdiff


def _toy():
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B", "B", "C", "A", "B", "C"],
            "t2": ["A", "B", "B", "C", "C", "A", "C", "A"],
            "t3": ["B", "B", "C", "C", "A", "B", "A", "B"],
            "t4": ["B", "C", "C", "A", "A", "B", "B", "C"],
        }
    )
    return seqdef(df), np.array([0, 0, 1, 1, 0, 1, 0, 1])


def test_seqdiff_structure():
    x, group = _toy()
    out = seqdiff(x, group)
    assert set(out.keys()) == {"stat", "discrepancy"}
    assert list(out["stat"].columns) == [
        "Pseudo F",
        "Pseudo Fbf",
        "Pseudo R2",
        "Bartlett",
        "Levene",
    ]
    assert list(out["discrepancy"].columns) == ["0", "1", "Total"]
    # cmprange (0, 1) over 4 columns -> 3 window positions.
    assert out["stat"].shape == (3, 5)
    assert out["discrepancy"].shape == (3, 3)


def test_seqdiff_pseudo_r2_in_unit_interval():
    x, group = _toy()
    out = seqdiff(x, group)
    r2 = out["stat"]["Pseudo R2"].to_numpy()
    assert np.all((r2 >= -1e-12) & (r2 <= 1.0 + 1e-12))


def test_seqdiff_requires_stslist():
    with pytest.raises(TypeError):
        seqdiff(np.zeros((4, 3), dtype=int), np.array([0, 1, 0, 1]))


def test_seqdiff_group_length_checked():
    x, _ = _toy()
    with pytest.raises(ValueError):
        seqdiff(x, np.array([0, 1]))


def test_seqdiff_cmprange_controls_number_of_positions():
    x, group = _toy()  # 4 columns
    # cmprange (0, 1): positions 1..(L-1) = 3.
    assert seqdiff(x, group, cmprange=(0, 1))["stat"].shape[0] == 3
    # cmprange (0, 2): positions 1..(L-2) = 2.
    assert seqdiff(x, group, cmprange=(0, 2))["stat"].shape[0] == 2


def test_seqdiff_weighted_differs_from_unweighted():
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B", "B", "C", "A", "B", "C"],
            "t2": ["A", "B", "B", "C", "C", "A", "C", "A"],
            "t3": ["B", "B", "C", "C", "A", "B", "A", "B"],
        }
    )
    w = np.array([3.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 3.0])
    x = seqdef(df, weights=w)
    group = np.array([0, 0, 1, 1, 0, 1, 0, 1])
    uw = seqdiff(x, group, weighted=False)["stat"].to_numpy()
    wt = seqdiff(x, group, weighted=True)["stat"].to_numpy()
    assert not np.allclose(uw, wt)


def test_seqdiff_default_is_lcs_auto_and_complete_map_requires_method():
    x, group = _toy()

    default = seqdiff(x, group)
    mapped = seqdiff(x, group, seqdist_args={"method": "LCS", "norm": "auto"})
    unnormalized = seqdiff(x, group, method="LCS", norm=None)

    pd.testing.assert_frame_equal(default["stat"], mapped["stat"])
    pd.testing.assert_frame_equal(default["discrepancy"], mapped["discrepancy"])
    assert not np.allclose(
        default["discrepancy"].to_numpy(),
        unnormalized["discrepancy"].to_numpy(),
    )
    with pytest.raises(TypeError, match="method"):
        seqdiff(x, group, seqdist_args={})


def test_seqdiff_map_forwards_custom_om_and_translates_sm_without_mutation():
    x, group = _toy()
    subm = np.array(
        [
            [0.0, 0.5, 1.5],
            [0.5, 0.0, 2.0],
            [1.5, 2.0, 0.0],
        ]
    )
    options = {"method": "OM", "norm": "none", "sm": subm, "indel": 0.75}
    original = dict(options)
    original_subm = subm.copy()

    via_sm = seqdiff(x, group, seqdist_args=options, weighted=False)
    via_subm = seqdiff(
        x,
        group,
        seqdist_args={
            "method": "OM",
            "norm": "none",
            "subm": subm,
            "indel": 0.75,
        },
        weighted=False,
    )

    pd.testing.assert_frame_equal(via_sm["stat"], via_subm["stat"])
    pd.testing.assert_frame_equal(via_sm["discrepancy"], via_subm["discrepancy"])
    assert options == original
    assert options["sm"] is subm
    np.testing.assert_array_equal(subm, original_subm)


def test_seqdiff_map_alias_conflicts_and_reserved_options():
    x, group = _toy()

    with pytest.deprecated_call(match="seqdist_arg"):
        alias = seqdiff(
            x,
            group,
            seqdist_arg={"method": "LCS", "norm": "none"},
            weighted=False,
        )
    modern = seqdiff(
        x,
        group,
        seqdist_args={"method": "LCS", "norm": "none"},
        weighted=False,
    )
    pd.testing.assert_frame_equal(alias["stat"], modern["stat"])

    with pytest.raises(ValueError, match="seqdist_args.*seqdist_arg"):
        seqdiff(x, group, seqdist_args={}, seqdist_arg={})
    with pytest.raises(TypeError, match="mapping"):
        seqdiff(x, group, seqdist_args=["method", "LCS"])
    with pytest.raises(ValueError, match="put method and norm in the map"):
        seqdiff(x, group, method="OM", seqdist_args={"method": "LCS"})
    with pytest.raises(ValueError, match="sm.*subm"):
        seqdiff(x, group, seqdist_args={"method": "OM", "sm": 2.0, "subm": 2.0})
    with pytest.raises(ValueError, match="window owns"):
        seqdiff(x, group, seqdist_args={"method": "LCS", "seqdata": x})
    with pytest.raises(ValueError, match="window owns"):
        seqdiff(x, group, seqdist_args={"method": "LCS", "x": x})
    with pytest.raises(ValueError, match="void"):
        seqdiff(x, group, seqdist_args={"method": "LCS", "void": 999})
    with pytest.raises(TypeError, match="Unsupported kwargs"):
        seqdiff(x, group, seqdist_args={"method": "LCS", "unsupported": True})
    for option in (
        {"refseq": 1},
        {"ref_sets": ([0], [1])},
        {"condensed": True},
        {"full_matrix": False},
    ):
        with pytest.raises(ValueError, match="square distance"):
            seqdiff(x, group, seqdist_args={"method": "LCS", **option})


def test_seqdiff_map_direct_missing_override_and_distance_weight_are_separate():
    mat = np.array(
        [
            ["A", "*", "B", "C"],
            ["A", "B", "B", "%"],
            ["B", "B", "C", "C"],
            ["B", "A", "C", "%"],
            ["C", "C", "A", "A"],
            ["C", "A", "A", "B"],
        ],
        dtype=object,
    )
    x = seqdef(
        mat,
        alphabet=["A", "B", "C"],
        nr="*",
        void="%",
        weights=[4.0, 1.0, 2.0, 1.0, 3.0, 1.0],
    )
    group = np.array([0, 0, 1, 1, 0, 1])

    overridden = seqdiff(
        x,
        group,
        seqdist_args={
            "method": "LCS",
            "norm": "none",
            "with_missing": False,
        },
        with_missing=True,
    )
    direct = seqdiff(x, group, method="LCS", norm=None, with_missing=True)
    pd.testing.assert_frame_equal(overridden["stat"], direct["stat"])

    chi_weighted = seqdiff(
        x,
        group,
        seqdist_args={"method": "CHI2", "step": 2, "overlap": True, "weighted": True},
        weighted=False,
        with_missing=True,
    )
    chi_unweighted = seqdiff(
        x,
        group,
        seqdist_args={"method": "CHI2", "step": 2, "overlap": True, "weighted": False},
        weighted=False,
        with_missing=True,
    )
    assert not np.allclose(
        chi_weighted["discrepancy"].to_numpy(),
        chi_unweighted["discrepancy"].to_numpy(),
    )


def test_seqdiff_sts_group_uses_position_states_and_declared_alphabet_order():
    x, _ = _toy()
    group = seqdef(
        np.array(
            [
                ["z", "a", "m", "z"],
                ["a", "a", "m", "m"],
                ["m", "z", "z", "a"],
                ["z", "m", "a", "a"],
                ["a", "z", "m", "z"],
                ["m", "a", "z", "m"],
                ["z", "m", "a", "z"],
                ["a", "z", "m", "a"],
            ],
            dtype=object,
        ),
        alphabet=["z", "a", "m"],
    )

    out = seqdiff(x, group, method="LCS", norm=None)

    assert list(out["discrepancy"].columns) == ["z", "a", "m", "Total"]
    # R labels each comparison window with the sequence's time label
    # (seqdiff.R:29,32 over colnames(seqdata)); an object built without explicit
    # column names carries R's default T1..TL (seqdef.R:254). This expectation
    # was previously the Python-authored 1..n index (cycle-6 §9).
    assert out["stat"].index.tolist() == ["T1", "T2", "T3"]
    assert np.isfinite(out["discrepancy"].to_numpy()).all()

    short_group = STSList(data=group.data[:, :2].copy(), alphabet=list(group.alphabet))
    with pytest.raises(ValueError, match="enough positions"):
        seqdiff(x, short_group)


def test_seqdiff_squared_is_forwarded_to_discrepancy_analysis():
    x, group = _toy()
    unsquared = seqdiff(x, group, squared=False)["discrepancy"].to_numpy()
    squared = seqdiff(x, group, squared=True)["discrepancy"].to_numpy()
    assert not np.allclose(unsquared, squared)


@pytest.mark.parametrize(
    "cmprange,error",
    [
        ((0,), ValueError),
        ((0, 1, 2), ValueError),
        ((0.0, 1), TypeError),
        ((False, 1), TypeError),
        (None, TypeError),
        ((1, 0), ValueError),
        ((0, 5), ValueError),
    ],
)
def test_seqdiff_validates_cmprange(cmprange, error):
    x, group = _toy()
    with pytest.raises(error, match="cmprange"):
        seqdiff(x, group, cmprange=cmprange)
