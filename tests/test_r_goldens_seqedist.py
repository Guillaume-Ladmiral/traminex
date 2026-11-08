"""R-golden tests for the public ``TraMineRextras::seqedist`` wrapper.

Public string/default fixtures are pinned to TraMineRextras 0.6.9. Separate
internal-engine fixtures pin traminex's backward-compatible numeric/boolean
option extensions against TraMineR 2.2-14.

Generator: scripts/r_goldens/make_goldens_seqedist.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqecreate, seqedist
from traminex.events import _mining

_IDCOST = [1, 2, 3]  # A, B, C
_V = 0.5


def _eseq():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s2", "s2", "s3", "s3", "s3", "s4", "s4"],
            "timestamp": [0, 2, 2, 0, 3, 0, 1, 4, 0, 2],
            "event": ["A", "B", "C", "A", "C", "B", "A", "C", "A", "B"],
        }
    )
    return seqecreate(tse)


def _assert_dist(interval, norm, tag):
    D = seqedist(_eseq(), idcost=_IDCOST, vparam=_V, interval=interval, norm=norm)
    golden = load_golden("seqedist_" + tag)
    np.testing.assert_allclose(D, golden, rtol=0, atol=R_GOLDEN_ATOL)
    assert np.allclose(D, D.T, rtol=0, atol=R_GOLDEN_ATOL)
    assert np.all(np.diag(D) == 0)
    return D


@pytest.mark.parametrize(
    ("interval", "norm", "tag"),
    [
        ("No", "none", "public_no_none"),
        ("No", "max", "public_no_max"),
        ("No", "YujianBo", "public_no_yujianbo"),
        ("previous", "none", "public_previous_none"),
        ("next", "none", "public_next_none"),
    ],
)
def test_seqedist_public_strings_match_r(interval, norm, tag):
    _assert_dist(interval, norm, tag)


def test_seqedist_public_defaults_match_r_and_explicit_strings():
    default = seqedist(_eseq(), idcost=_IDCOST, vparam=_V)
    explicit = seqedist(
        _eseq(),
        idcost=_IDCOST,
        vparam=_V,
        interval="No",
        norm="YujianBo",
    )
    golden = load_golden("seqedist_public_default")

    np.testing.assert_allclose(default, golden, rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(default, explicit, rtol=0, atol=R_GOLDEN_ATOL)
    assert np.allclose(default, default.T, rtol=0, atol=R_GOLDEN_ATOL)
    assert np.all(np.diag(default) == 0)


@pytest.mark.parametrize(
    ("interval", "norm", "canonical_interval", "canonical_norm"),
    [
        ("p", "n", "previous", "none"),
        ("ne", "m", "next", "max"),
        ("NO", "Y", "No", "YujianBo"),
    ],
)
def test_seqedist_public_options_are_case_insensitive_unique_partials(
    interval, norm, canonical_interval, canonical_norm
):
    partial = seqedist(_eseq(), _IDCOST, _V, interval=interval, norm=norm)
    canonical = seqedist(
        _eseq(),
        _IDCOST,
        _V,
        interval=canonical_interval,
        norm=canonical_norm,
    )
    np.testing.assert_allclose(partial, canonical, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("interval", "n"),
        ("interval", ""),
        ("interval", "later"),
        ("norm", ""),
        ("norm", "unknown"),
    ],
)
def test_seqedist_rejects_ambiguous_empty_and_unknown_public_options(option, value):
    kwargs = {"interval": "No", "norm": "YujianBo", option: value}
    with pytest.raises(ValueError, match=option):
        seqedist(_eseq(), _IDCOST, _V, **kwargs)


def test_seqedist_rejects_idcost_length_mismatch_before_engine(monkeypatch):
    def forbidden_engine_call(*args, **kwargs):
        raise AssertionError("distance engine must not run")

    monkeypatch.setattr(_mining._engine, "tmrseqedist", forbidden_engine_call)
    with pytest.raises(ValueError, match="one indel cost per event"):
        seqedist(_eseq(), [1, 2], _V)


@pytest.mark.parametrize(
    ("interval", "norm", "tag"),
    [
        (1, 0, "internal_i1_n0"),
        (1, 1, "internal_i1_n1"),
        (1, 2, "internal_i1_n2"),
        (0, 0, "public_no_none"),
        (2, 0, "public_previous_none"),
        (3, 0, "public_next_none"),
    ],
)
def test_seqedist_integer_codes_remain_supported(interval, norm, tag):
    _assert_dist(interval, norm, tag)


@pytest.mark.parametrize(("interval", "norm"), [(False, False), (True, True)])
def test_seqedist_boolean_codes_remain_supported(interval, norm):
    boolean_result = seqedist(_eseq(), _IDCOST, _V, interval=interval, norm=norm)
    integer_result = seqedist(
        _eseq(),
        _IDCOST,
        _V,
        interval=int(interval),
        norm=int(norm),
    )
    np.testing.assert_allclose(boolean_result, integer_result, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("interval", -1),
        ("interval", 4),
        ("norm", -1),
        ("norm", 3),
        ("interval", 1.5),
        ("norm", 1.0),
    ],
)
def test_seqedist_rejects_invalid_numeric_codes(option, value):
    kwargs = {"interval": 0, "norm": 0, option: value}
    with pytest.raises((TypeError, ValueError), match=option):
        seqedist(_eseq(), _IDCOST, _V, **kwargs)


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("interval", ["No"]),
        ("interval", np.array(["No"])),
        ("norm", ["none"]),
        ("norm", np.array(["none"])),
    ],
)
def test_seqedist_rejects_non_scalar_options(option, value):
    kwargs = {"interval": 0, "norm": 0, option: value}
    with pytest.raises((TypeError, ValueError), match=option):
        seqedist(_eseq(), _IDCOST, _V, **kwargs)
