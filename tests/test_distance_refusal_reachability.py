"""Every uncovered refusal in the two distance modules, adjudicated (W-10).

``traminex/distances/__init__.py`` and ``python_kernels.py`` are the
lowest-covered core modules. Cycle 6 raised each by about a point and left the
rest, having judged that several are defence-in-depth behind validation that
fires first. This module confirms that judgement one branch at a time rather
than inheriting it: each refusal is either **reached through a public path**,
and then pinned here, or **pre-empted by an earlier guard**, and then recorded
with the guard that pre-empts it.

Two things this deliberately does not do. It does not chase a percentage: the
``@njit`` bodies in ``python_kernels.py`` account for 157 of that module's 285
uncovered lines and never report Python line coverage whether or not they run,
so they are neither counted as dead nor "covered" by a test that does not
execute them. And it adds no assertion-free tests: every case below asserts the
message, because the message is the contract a user meets.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex import distances as D

MAT = np.array([[0, 1, 2, 1], [1, 0, 2, 2], [2, 2, 0, 1]], dtype=int)


def _sequences():
    return tx.seqdef(pd.DataFrame([list("ABCB"), list("BACC"), list("CCAB")]))


def _ragged():
    return tx.seqdef(pd.DataFrame([["A", "B", "C", "B"], ["B", "A", None, None]]))


# --------------------------------------------------------------------------- #
# 1. Reachable through the public dispatcher: seqdist / pairwise.
# --------------------------------------------------------------------------- #
_VOID_RAGGED = [list("ABAB"), ["A", "B", "%", "%"], ["B", "A", "B", "%"]]


@pytest.mark.parametrize(
    "extension",
    [
        {"subm_by_pos": np.stack([np.array([[0.0, 2.0], [2.0, 0.0]])] * 4)},
        {"indel_by_pos": np.ones(4)},
    ],
    ids=["subm_by_pos", "indel_by_pos"],
)
def test_position_dependent_om_costs_are_refused_for_ragged_sequences(extension):
    """The traminex-only position costs have no ragged definition (docs/OMLOC.md).

    Reachable from public ``seqdist`` but untested: removing the refusal left
    the whole suite green (2026-09-22 audit S10-09).
    """
    ragged = tx.seqdef(np.array(_VOID_RAGGED, dtype=object), alphabet=["A", "B"])
    with pytest.raises(NotImplementedError, match="not supported with ragged sequences"):
        tx.seqdist(ragged, method="OM", subm=2.0, indel=1.0, **extension)


@pytest.mark.parametrize(
    "legacy", [{"variant": "mean"}, {"smooth_window": 1}], ids=["variant", "smooth_window"]
)
def test_legacy_dhd_refuses_a_short_external_reference(legacy):
    """Only R-faithful DHD handles an external reference shorter than the data.

    The legacy variants are traminex-only; reachable from public ``seqdist``
    but untested (2026-09-22 audit S10-09).
    """
    full = tx.seqdef(
        np.array([list("ABAB"), list("BBAA"), list("BABA")], dtype=object), alphabet=["A", "B"]
    )
    short = tx.seqdef(np.array([["A", "B", "%", "%"]], dtype=object), alphabet=["A", "B"])
    with pytest.raises(NotImplementedError, match="supported only for R-faithful DHD"):
        tx.seqdist(full, method="DHD", refseq=short, **legacy)


@pytest.mark.parametrize(
    "kwargs,message",
    [
        ({"method": "DHD", "dhd_profile": "bogus"}, "dhd_profile must be"),
        (
            {"method": "OM", "subm": "TRATE", "indel": 1.0, "trate_profile": "bogus"},
            "trate_profile must be",
        ),
        ({"method": "OM", "indel": 1.0, "subm": float("nan")}, "must be finite"),
        # a valid 2x2 cost matrix, so the size refusal is what fires rather than
        # TraMineR's triangle-inequality diagnostic
        (
            {"method": "OM", "indel": 1.0, "subm": np.array([[0.0, 2.0], [2.0, 0.0]])},
            # a sequence object's alphabet fixes the width exactly, as R's
            # checkcost does (2026-09-22 audit S1b-01)
            "subm must be 3x3",
        ),
        ({"method": "OM", "indel": 0.0, "subm": 2.0}, "indel cost must be strictly positive"),
        (
            {"method": "OM", "indel": [1.0, 0.0, 1.0], "subm": 2.0},
            "indel cost must be strictly positive",
        ),
        ({"method": "HAM", "subm": 2.0}, "subm must be square"),
        (
            {"method": "TWED", "subm": 2.0, "nu": 0.5, "h": 0.5, "indel": [1.0, 1.0]},
            # R checks a vector indel's length for every method (seqdist.R:137-147,
            # TX-AUD-20261010-071), before TWED's own rule.
            "'indel' must contain a cost for each state",
        ),
        (
            {"method": "TWED", "subm": 2.0, "nu": 0.5, "h": 0.5, "indel": 0.0},
            "TWED indel cost must be strictly positive",
        ),
        ({"method": "NMS", "prox": np.zeros((2, 2))}, "prox must be of size"),
        ({"method": "DHD", "subm": np.zeros((2, 2, 2))}, "time-varying array"),
        ({"method": "CHI2", "step": 1, "weights": [1.0]}, "weights length must equal"),
    ],
)
def test_seqdist_refusals_are_reachable_and_name_the_problem(kwargs, message):
    with pytest.raises((ValueError, TypeError), match=message):
        tx.seqdist(_sequences(), **kwargs)


def test_hamming_refuses_ragged_sequences():
    """R raises the same way; the message is quoted from TraMineR."""
    with pytest.raises(ValueError, match="HAM is not defined for sequences of different length"):
        tx.seqdist(_ragged(), method="HAM")


@pytest.mark.parametrize(
    "kwargs,message",
    [
        ({"method": "OM", "indel": 1.0, "subm": np.zeros((3, 2))}, "subm must be square"),
        ({"method": "OM", "subm": "CONSTANT", "indel": 1.0}, "subm accepts float"),
        ({"method": "CHI2", "weights": [1.0]}, "weights length must equal"),
    ],
)
def test_low_level_pairwise_refusals_are_reachable(kwargs, message):
    with pytest.raises((ValueError, TypeError), match=message):
        D.pairwise(MAT, **kwargs)


@pytest.mark.parametrize(
    "method,extra",
    [("OM", {}), ("TWED", {"nu": 0.5, "h": 0.5})],
)
def test_a_string_indel_is_refused_before_the_per_method_guard_sees_it(method, extra):
    """Records a pre-emption rather than claiming to cover the later guard.

    ``_state_indel_spec`` and ``_twed_indel_value`` in ``python_kernels.py``
    both refuse a string ``indel``, but
    ``validate_indel_value`` in ``_validation.py`` refuses it first, with a
    message that also says where ``'auto'`` *is* resolved. So the two later
    guards are defence-in-depth: no public input reaches them.
    """
    with pytest.raises(ValueError, match="is not a number; pairwise accepts") as excinfo:
        D.pairwise(MAT, method=method, indel="auto", subm=2.0, **extra)

    # the message is the pre-empting guard's, not the per-method one's
    assert "indel must be a scalar or state-cost vector" not in str(excinfo.value)


# --------------------------------------------------------------------------- #
# 2. Reference-set and reference-sequence index validation.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "refseq,message",
    [
        ([[1], [99]], "outside the sequence data"),
        ([[0], [1]], "outside the sequence data"),
        ([["a"], [1]], "must contain positive integer indices"),
        ([[], [1]], "must contain positive integer indices"),
        ([[1.5], [2]], "must contain positive integer indices"),
    ],
)
def test_reference_index_sets_are_validated(refseq, message):
    with pytest.raises(ValueError, match=message):
        tx.seqdist(_sequences(), method="LCS", refseq=refseq)


def test_an_external_reference_must_share_the_missing_symbol():
    external = tx.seqdef(pd.DataFrame([list("ABCB")]), nr="#")

    with pytest.raises(ValueError, match="same missing-value symbol"):
        tx.seqdist(_sequences(), method="LCS", refseq=external)


@pytest.mark.parametrize(
    "value,message",
    [
        (3.5, "must be a 2D matrix or iterable of rows"),
        ([[[0, 1]], [[1, 0]]], "must contain one-dimensional sequence rows"),
    ],
)
def test_between_sets_row_coercion_refuses_malformed_input(value, message):
    other = np.array([[0, 1], [1, 0]], dtype=int)

    with pytest.raises(ValueError, match=message):
        D.pairwise_between_sets(value, other, method="LCS")


# --------------------------------------------------------------------------- #
# 3. Matrix-shape helpers.
# --------------------------------------------------------------------------- #
def test_sanitize_distance_matrix_requires_a_square_matrix():
    with pytest.raises(ValueError, match="D must be a square matrix"):
        D.sanitize_distance_matrix(np.zeros((2, 3)))


def test_to_condensed_requires_a_square_matrix():
    with pytest.raises(ValueError, match="D must be a square matrix"):
        D.to_condensed(np.zeros((2, 3)))


def test_from_condensed_rejects_an_impossible_length():
    with pytest.raises(ValueError, match="not a valid condensed size"):
        D.from_condensed(np.zeros(5))


# --------------------------------------------------------------------------- #
# 4. seqsubm's per-method keyword vocabularies.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "method", ["TRATE", "CONSTANT", "INDELS", "INDELSLOG", "FEATURES", "FUTURE"]
)
def test_every_seqsubm_method_refuses_an_unknown_keyword(method):
    kwargs: dict[str, object] = {"bogus": 1}
    if method == "FEATURES":
        kwargs["state_features"] = np.zeros((3, 2))

    with pytest.raises(TypeError, match=f"Unsupported kwargs for seqsubm\\(method='{method}'\\)"):
        D.seqsubm(MAT, method=method, **kwargs)


def test_seqsubm_names_the_replacement_for_the_removed_constant_keyword():
    with pytest.raises(TypeError, match="constant is not supported; use cval="):
        D.seqsubm(MAT, method="CONSTANT", constant=2.0)


def test_seqsubm_features_requires_the_feature_table():
    with pytest.raises(ValueError, match="requires state_features"):
        D.seqsubm(MAT, method="FEATURES")


def test_seqsubm_refuses_an_unknown_method_and_lists_the_supported_ones():
    with pytest.raises(ValueError, match="Unknown seqsubm method"):
        D.seqsubm(MAT, method="NOPE")


# --------------------------------------------------------------------------- #
# 5. Unreachable through any public path -- recorded, not covered.
#
# Each of these is defence-in-depth: an earlier guard fires first, so no input
# a user can supply reaches it. Adding a test that called the private helper
# directly would raise the coverage number without testing anything a user can
# meet, which is what this section exists to avoid.
# --------------------------------------------------------------------------- #
def test_map_norms_trailing_guard_cannot_fire_from_any_caller():
    """``_map_norm``'s last ``raise ValueError("norm is not supported for method")``.

    Its guard is ``m not in {"OM", "LCS", "LCP", "RLCP"}``, but ``_map_norm``
    has two guarded call sites, in ``pairwise`` and ``_low_level_rows``.
    Both require the classic method family before calling it. The condition
    is the negation of its caller's, so it is unreachable by construction.
    """
    import inspect

    source = inspect.getsource(D)
    assert source.count("_map_norm(") == 3  # definition plus two guarded calls
    call_context = source.split("mapped_norm = _map_norm(")[0]
    assert call_context.rstrip().endswith("if method.upper() in _CLASSIC_NORM_METHODS:")
    row_context = source.split("norm = _map_norm(m, norm)")[0]
    assert row_context.rstrip().endswith("if m in _CLASSIC_NORM_METHODS:")
    assert D._CLASSIC_NORM_METHODS == {"OM", "LCS", "LCP", "RLCP"}

    # Every other method carries norm straight through, so none of them reaches it.
    for method in ("DHD", "TWED", "OMSPELL"):
        options = {"subm": 2.0, "indel": 1.0}
        if method == "TWED":
            options.update(nu=0.5, h=0.5)
        if method == "DHD":
            options = {}
        D.pairwise(MAT, method=method, norm="gmean", **options)


@pytest.mark.parametrize("method,transport", [("DHD", "_n_states"), ("CHI2", "_n_states")])
def test_the_declared_state_count_guards_are_behind_kwarg_validation(method, transport):
    """``seqdist``'s DHD "declared state count" guard and ``_chi2_euclid_pairwise``'s.

    Both compare a declared state count against the observed one, and both are
    reached only through the private ``_n_states`` transport key that
    ``seqdist`` sets for itself. ``validate_pairwise_kwargs`` rejects it on the
    public path, so a user cannot get far enough to be told the count is wrong.
    """
    with pytest.raises(TypeError, match=f"Unsupported kwargs for method {method}"):
        D.pairwise(MAT, method=method, **{transport: 1})
