"""Refusal branches of the distance dispatcher and its kernels (cycle-6 §8.2).

``traminex/distances/__init__.py`` (81%) and
``traminex/distances/python_kernels.py`` (79%) are the two lowest-covered core
modules and hold the two highest-complexity functions. The cycle-6 audit showed
no *defect* in any uncovered branch, so this is hardening, not a fix: the
branches exercised below are the ones with analytical consequence -- guards
whose silent failure would return a wrong matrix or accept an input the R
contract refuses. Each asserts the message, not merely that something raised,
so a guard cannot be weakened into a different refusal unnoticed.

Numba ``@njit`` bodies never report Python line coverage, so nothing here tries
to "cover" them.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise, seqsubm
from traminex.distances import python_kernels as kernels

_RAW = np.array(
    [["A", "B", "A", "B"], ["B", "A", "B", "A"], ["A", "A", "B", "B"]],
    dtype=object,
)


def _sequences():
    return tx.seqdef(_RAW.copy(), alphabet=["A", "B"])


def _matrix():
    return np.asarray(_sequences().data, dtype=np.int32)


# --- dispatcher-layer refusals ------------------------------------------------


@pytest.mark.parametrize("method", ["NMS", "NMSMST", "SVRSPELL"])
def test_pairwise_refuses_norm_for_the_methods_traminer_refuses_it_for(method):
    # seqdist.R:205 refuses norm for these three; the low-level layer must too,
    # or a caller could silently get an unnormalized matrix it believes is scaled.
    with pytest.raises(ValueError, match="norm is not"):
        pairwise(_matrix(), method=method, norm="length")


def test_seqdist_refuses_an_unknown_dhd_profile():
    with pytest.raises(ValueError, match="dhd_profile must be 'default' or 'calibrated'"):
        pairwise(_matrix(), method="DHD", dhd_profile="bogus")


def test_seqdist_refuses_an_unknown_trate_profile():
    with pytest.raises(ValueError, match="trate_profile must be 'default' or 'calibrated'"):
        tx.seqdist(_sequences(), method="OM", subm="TRATE", indel=1.0, trate_profile="bogus")


def test_conflicting_output_forms_are_refused():
    # full_matrix=False means condensed=True, so these two agree and are accepted;
    # only a genuine contradiction is refused.
    assert np.asarray(tx.seqdist(_sequences(), method="LCS", full_matrix=False)).ndim == 1
    with pytest.raises(ValueError, match="conflicting output forms"):
        tx.seqdist(_sequences(), method="LCS", full_matrix=True, condensed=True)


@pytest.mark.parametrize(
    ("refseq", "expected"),
    [
        (2.5, "refseq must be a non-negative integer index"),
        (99, "outside the sequence data"),
        # R calls a negative index an invalid value, not an out-of-range one
        # (is.a.positive.integer; TX-AUD-20261011-044).
        (-1, "refseq must be a non-negative integer index"),
        # R uses the first two of a longer list (2026-09-22 audit S1a-12), so the
        # refusal is exercised with a single set.
        ([[1]], "two 1-based index sets"),
        ("nonsense", "refseq must be an integer, two index sets, or one STS sequence"),
    ],
)
def test_refseq_refusals_name_the_shape_they_expect(refseq, expected):
    with pytest.raises((ValueError, IndexError, TypeError), match=expected):
        tx.seqdist(_sequences(), method="LCS", refseq=refseq)


def test_seqllcs_and_seqllcp_refuse_mismatched_alphabets():
    left = tx.seqdef(np.array([["A", "B", "A"]], dtype=object), alphabet=["A", "B"])
    right = tx.seqdef(np.array([["X", "Y", "X"]], dtype=object), alphabet=["X", "Y"])
    for function in (tx.seqLLCS, tx.seqLLCP):
        with pytest.raises(ValueError, match="alphabet of both sequences"):
            function(left, right)


# --- cost-construction refusals ----------------------------------------------


def test_seqsubm_features_requires_one_row_per_state():
    with pytest.raises(ValueError, match="state_features"):
        seqsubm(_sequences(), method="FEATURES")
    with pytest.raises(ValueError, match="exactly one row per state"):
        seqsubm(_sequences(), method="FEATURES", state_features=np.zeros((5, 2)))


def test_seqdist_rejects_the_r_constant_spelling_with_a_pointer_to_cval():
    with pytest.raises(TypeError, match=r"use seqsubm\(\.\.\., method='CONSTANT', cval="):
        tx.seqdist(_sequences(), method="OM", subm="CONSTANT", constant=2.0)


# --- kernel-layer refusals ----------------------------------------------------


def test_numba_switch_rejects_an_unknown_value(monkeypatch):
    monkeypatch.setenv("TRAMINEX_NUMBA", "sometimes")
    with pytest.raises(ValueError, match="TRAMINEX_NUMBA must be one of"):
        kernels._resolve_numba_mode()


def test_om_subm_must_be_square():
    with pytest.raises(ValueError, match="subm must be square"):
        pairwise(_matrix(), method="OM", indel=1.0, subm=np.zeros((2, 3)))


def test_om_subm_must_be_finite():
    subm = np.array([[0.0, np.inf], [np.inf, 0.0]])
    with pytest.raises(ValueError, match="substitution costs must be finite"):
        pairwise(_matrix(), method="OM", indel=1.0, subm=subm)


@pytest.mark.parametrize("indel", [0.0, -1.0])
def test_om_indel_must_be_strictly_positive(indel):
    with pytest.raises(ValueError, match="indel cost must be strictly positive"):
        pairwise(_matrix(), method="OM", indel=indel, subm=2.0)


@pytest.mark.parametrize("method", ["OMSLEN", "OMSPELL"])
def test_state_cost_indel_vector_length_is_checked(method):
    # These methods accept a per-state indel vector; a wrong length must be named
    # rather than broadcast into a wrong distance.
    with pytest.raises(ValueError, match="indel vector length must equal the active state count"):
        pairwise(_matrix(), method=method, indel=np.array([1.0, 1.0, 1.0]), subm=2.0)


@pytest.mark.parametrize("method", ["OM", "LCP", "RLCP"])
def test_norm_domain_at_the_kernel_layer_names_its_own_set(method):
    options = {"indel": 1.0, "subm": 2.0} if method == "OM" else {}
    with pytest.raises(ValueError, match="low-level pairwise"):
        pairwise(_matrix(), method=method, norm="bogus", **options)


def test_omslen_link_domain():
    with pytest.raises(ValueError, match="OMslen link must be 'mean' or 'gmean'"):
        pairwise(_matrix(), method="OMSLEN", indel=1.0, subm=2.0, link="median")


def test_omslen_h_domain():
    with pytest.raises(ValueError, match="OMslen 'h' must be greater than or equal to 0"):
        pairwise(_matrix(), method="OMSLEN", indel=1.0, subm=2.0, h=-1.0)


@pytest.mark.parametrize(
    ("options", "expected"),
    [
        ({"nu": 0.0, "h": 1.0}, "TWED 'nu' must be strictly greater than 0"),
        ({"nu": 0.5, "h": -1.0}, "TWED 'h' must be greater than or equal to 0"),
    ],
)
def test_twed_parameter_domains(options, expected):
    with pytest.raises(ValueError, match=expected):
        pairwise(_matrix(), method="TWED", subm=2.0, indel=1.0, **options)


@pytest.mark.parametrize("otto", [-0.5, 1.5])
def test_omstran_otto_domain(otto):
    with pytest.raises(ValueError, match=r"OMstran 'otto' must be in \[0, 1\]"):
        pairwise(_matrix(), method="OMSTRAN", indel=1.0, subm=2.0, otto=otto)


def test_omloc_refuses_negative_context_parameters():
    with pytest.raises(ValueError, match="expcost and context must be non-negative"):
        pairwise(_matrix(), method="OMLOC", subm=2.0, expcost=-1.0, context=0.5)


def test_dhd_refuses_a_three_dimensional_subm_with_the_wrong_time_axis():
    length = _matrix().shape[1]
    bad = np.zeros((length + 3, 2, 2), dtype=float)
    with pytest.raises(ValueError, match=f"first axis L={length}"):
        pairwise(_matrix(), method="DHD", subm=bad)


def test_dhd_refuses_a_subm_of_the_wrong_rank():
    with pytest.raises(ValueError, match="2D \\(K,K\\) or 3D \\(L,K,K\\)"):
        pairwise(_matrix(), method="DHD", subm=np.zeros((2, 2, 2, 2)))


def test_dhd_refuses_an_unknown_cost_string_at_the_seqdist_layer():
    # R's text (seqdist.R:253; TX-AUD-20261019-079): CONSTANT is the one R cost name DHD
    # cannot use, since seqdist refuses names outside R's list first.
    with pytest.raises(ValueError, match="is not relevant for DHD, consider HAM instead"):
        tx.seqdist(_sequences(), method="DHD", subm="CONSTANT")


@pytest.mark.parametrize("method", ["HAM", "DHD"])
def test_equal_length_methods_refuse_unequal_rows(method):
    ragged = tx.seqdef(
        np.array([["A", "B", "A", "B"], ["B", "A", "%", "%"]], dtype=object),
        alphabet=["A", "B"],
    )
    with pytest.raises(ValueError, match="different length"):
        tx.seqdist(ragged, method=method)


def test_declared_state_count_cannot_hide_an_observed_state():
    # n_states below the largest observed code would index outside the cost table
    # and silently return a wrong matrix, so the DHD cost builders refuse it.
    with pytest.raises(ValueError, match="cannot be smaller than an observed state code"):
        kernels._dhd_costs_traminer(_matrix(), cval=4.0, n_states=1)
    with pytest.raises(ValueError, match="cannot be smaller than an observed state code"):
        kernels._dhd_indels_costs(_matrix(), method="INDELS", n_states=1)


def test_wham_weights_length_is_checked_against_the_sequence_length():
    with pytest.raises(ValueError, match="weights must be one-dimensional"):
        pairwise(_matrix(), method="WHAM", weights=np.ones(2))
