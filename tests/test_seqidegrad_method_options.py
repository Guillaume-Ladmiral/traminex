"""The degradation family honours TraMineR's ``method`` and ``penalized`` options.

``seqidegrad.private`` (``seqidegrad.R:74-119``, ``:130-131``, ``:334-337``)
maps a logical ``penalized`` to ``"BOTH"``/``"NO"``; refuses a ``method``
outside ``FREQ, TRATE, TRATEDSS, RANK, FREQ+, TRATE+, TRATEDSS+, RANK+, ONE``;
for the transition-rate methods refuses a ``weight.type`` outside
``ADD, INV, LOGINV`` and a ``border.effect`` not strictly above one; refuses a
``penalized`` outside ``NEG, POS, BOTH, NO``; gives every transition weight 1
for ``"ONE"``; and for the ``"+"`` methods multiplies the result by the mean
transition weight per spell.  traminex computed ``TRATE`` for every method
other than ``RANK``/``FREQ``/``TRATEDSS`` (including unknown ones), returned
zeros for ``penalized=True`` or an unknown ``penalized``, and accepted any
``weight_type`` and ``border_effect`` (2026-09-22 audit S3-02).  The private
term is shared by ``seqidegrad``, ``seqprecarity``, ``seqinsecurity`` and the
ranked ``seqindic`` columns, so each route is checked.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-02), live TraMineR 2.2-14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens" / "audit20260922_s3a.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDENS.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, *_ = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _r_values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", R[case]
    return np.array(payload.split(" "), dtype=float)


def _dg():
    rows = [
        ["A", "A", "B", "B", "C", "C"],
        ["A", "B", "B", "C", "C", "C"],
        ["B", "B", "B", "A", "A", "C"],
        ["C", "C", "A", "A", "B", "B"],
        ["A", "A", "A", "A", "A", "A"],
        ["B", "C", "B", "C", "B", "C"],
    ]
    return tx.seqdef(
        np.array(rows, dtype=object), alphabet=["A", "B", "C"], weights=[1, 2, 0.5, 1.5, 3, 1]
    )


def _dgm():
    rows = [["A", "B", None, "C"], [None, None, None, None], ["C", "C", "A", "B"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"], right=None)


def _check(actual, case: str) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float).ravel(), _r_values(case), rtol=0, atol=1e-12
    )


@pytest.mark.parametrize("method", ["RANK+", "TRATE+", "TRATEDSS+", "FREQ+", "ONE"])
def test_seqidegrad_one_and_plus_methods_match_r(method):
    _check(tx.seqidegrad(_dg(), method=method), f"s3_02_degrad_{method}")
    _check(
        tx.seqidegrad(_dg(), method=method, spell_integr=False),
        f"s3_02_degrad_{method}_nointegr",
    )


@pytest.mark.parametrize(
    "case,kwargs",
    [
        ("s3_02_degrad_penalized_true", {"penalized": True}),
        ("s3_02_degrad_penalized_true", {"penalized": np.bool_(True)}),
        ("s3_02_degrad_one_neg", {"method": "ONE", "penalized": "NEG"}),
        ("s3_02_degrad_one_partial", {"method": "ONE", "state_order": ["A", "B"]}),
        ("s3_02_degrad_one_equiv", {"method": "ONE", "state_equiv": [["B", "C"]]}),
        # weight_type and border_effect are validated only for the rate methods.
        ("s3_02_degrad_rank_bogus_weight", {"method": "RANK", "weight_type": "BOGUS"}),
        ("s3_02_degrad_rank_border_half", {"method": "RANK", "border_effect": 0.5}),
        ("s3_02_degrad_one_bogus_weight", {"method": "ONE", "weight_type": "BOGUS"}),
    ],
)
def test_seqidegrad_options_match_r(case, kwargs):
    _check(tx.seqidegrad(_dg(), **kwargs), case)


def test_seqidegrad_logical_false_penalized_is_the_zero_correction():
    """``penalized=FALSE`` is ``"NO"`` in R; the private term returns 0."""
    np.testing.assert_array_equal(tx.seqidegrad(_dg(), penalized=False), np.zeros(6))


@pytest.mark.parametrize(
    "case,kwargs,match",
    [
        ("s3_02_degrad_bogus_method", {"method": "BOGUS"}, "'method' must be one of"),
        ("s3_02_degrad_lowercase_method", {"method": "rank"}, "'method' must be one of"),
        ("s3_02_degrad_bogus_penalized", {"penalized": "BOGUS"}, "'penalized' must be one of"),
        (
            "s3_02_degrad_trateplus_bogus_weight",
            {"method": "TRATE+", "weight_type": "BOGUS"},
            "'weight_type' must be one of",
        ),
        (
            "s3_02_degrad_freq_border_one",
            {"method": "FREQ", "border_effect": 1},
            "border_effect should be strictly greater than one",
        ),
    ],
)
def test_seqidegrad_refuses_what_r_refuses(case, kwargs, match):
    assert R[case][0] == "error"
    with pytest.raises(ValueError, match=match):
        tx.seqidegrad(_dg(), **kwargs)


@pytest.mark.parametrize(
    "case,kwargs",
    [
        ("s3_02_prec_one", {"method": "ONE"}),
        ("s3_02_prec_tratedssplus", {"method": "TRATEDSS+"}),
        ("s3_02_prec_penalized_true", {"penalized": True}),
        ("s3_02_prec_penalized_false", {"penalized": False}),
        # "NO" returns before the method is validated, as in R.
        ("s3_02_prec_no_bogus_method", {"penalized": "NO", "method": "BOGUS"}),
    ],
)
def test_seqprecarity_routes_the_options_as_r(case, kwargs):
    _check(tx.seqprecarity(_dg(), **kwargs), case)


@pytest.mark.parametrize(
    "case,method", [("s3_02_insec_one", "ONE"), ("s3_02_insec_rankplus", "RANK+")]
)
def test_seqinsecurity_routes_the_methods_as_r(case, method):
    _check(tx.seqinsecurity(_dg(), method=method), case)


@pytest.mark.parametrize(
    "case,prec_args",
    [
        ("s3_02_indic_one", {"method": "ONE"}),
        ("s3_02_indic_trateplus_logical", {"method": "TRATE+", "penalized": True}),
    ],
)
def test_ranked_seqindic_columns_route_the_options_as_r(case, prec_args):
    table = tx.seqindic(_dg(), indic=["degrad", "prec", "insec"], prec_args=prec_args)
    assert list(table.columns) == ["Degrad", "Prec", "Insec"]
    _check(table.to_numpy().ravel(order="F"), case)


def test_plus_method_on_an_empty_sequence_is_nan_as_in_r():
    """R divides the zero transition weight by a zero spell count (0/0)."""
    _check(tx.seqidegrad(_dgm(), method="RANK+"), "s3_02_degrad_missing_RANK+")
    _check(
        tx.seqidegrad(_dgm(), method="RANK+", with_missing=True),
        "s3_02_degrad_missing_RANK+_wm",
    )
