"""Consumed-component gaps, part 4: distances.

``seqdist``'s rarely used routes, the low-level ``pairwise`` and the Python kernels' refusals.
See the supervision repository's ``docs/TRAMINEX_COMPONENT_GAPS.md`` (part 4, re-measured on
2026-11-08 at ``2ed89d77``). A rare route is checked against an equivalent common one, which
the R-golden tests cover, against TraMineR's C source, or against R's values from
``scripts/r_goldens/make_goldens_gap_part4.R``.
"""

from __future__ import annotations

import csv
import types
import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

import traminex as tx
from traminex import distances as dist
from traminex.distances import _batched, _normalization
from traminex.distances import python_kernels as kernels

FRAME = pd.DataFrame(
    {"a": ["A", "B", "A", "C"], "b": ["B", "B", "A", "C"], "c": ["C", "A", "B", "A"]}
)
ALTERNATING = pd.DataFrame(
    {"a": list("ABAB"), "b": list("BABA"), "c": list("ABAB"), "d": list("BABA")}
)


def _read_golden(name: str) -> dict[tuple[str, str], str]:
    with (R_GOLDENS_DIR / name).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(line for line in handle if not line.startswith("#")))
    return {(case, key): value for case, key, value in rows[1:]}


#: R's values, from scripts/r_goldens/make_goldens_gap_part4.R.
R = _read_golden("gap_part4.csv")


def _numbers(case: str) -> np.ndarray:
    return np.array([float(value) for value in R[case, "distances"].split()])


@pytest.fixture(scope="module")
def seqs():
    return tx.seqdef(FRAME)


@pytest.fixture(scope="module")
def codes(seqs):
    return np.asarray(seqs.data)


@pytest.fixture(autouse=True)
def _quiet_cost_diagnostics():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


# --- seqdist ------------------------------------------------------------------------------------


def test_trate_with_a_zero_indel_computes_for_reference_sets_as_for_the_full_matrix():
    alternating = tx.seqdef(ALTERNATING)
    full = tx.seqdist(alternating, method="OM", subm="TRATE")
    block = tx.seqdist(alternating, method="OM", subm="TRATE", refseq=[[1, 2], [3, 4]])
    assert np.allclose(full.ravel(order="F"), _numbers("trate_zero_full"))
    assert np.allclose(block.ravel(order="F"), _numbers("trate_zero_refsets"))
    assert np.array_equal(block, full[np.ix_([0, 1], [2, 3])])
    with pytest.raises(ValueError, match="OM indel cost must be strictly positive"):
        tx.seqdist(
            alternating, method="OM", subm=np.zeros((2, 2)), indel=0, refseq=[[1, 2], [3, 4]]
        )


@pytest.mark.parametrize("subm", ["TRATE", "zero-matrix"])
def test_pairwise_reference_sets_carry_seqdists_zero_indel_allowance(subm):
    rows = np.asarray(tx.seqdef(ALTERNATING).data)
    options = {"subm": "TRATE"} if subm == "TRATE" else {"subm": np.zeros((2, 2)), "indel": 0.0}
    block = dist.pairwise(
        rows, method="OM", _allow_zero_indel=True, ref_sets=([0, 1], [2, 3]), **options
    )
    assert np.array_equal(np.asarray(block), np.zeros((2, 2)))


def test_the_modal_reference_of_a_plain_matrix_equals_that_of_the_object(seqs, codes):
    on_object = tx.seqdist(seqs, method="OM", subm="CONSTANT", indel=1, refseq=0)
    on_matrix = tx.seqdist(codes, method="OM", subm="CONSTANT", indel=1, refseq=0)
    assert np.array_equal(on_object, on_matrix)


def test_the_modal_reference_of_all_missing_sequences_is_taken_from_all_rows():
    missing = tx.seqdef(
        pd.DataFrame({"a": [None, None], "b": [None, None]}), right=None, left=None, gaps=None
    )
    result = tx.seqdist(missing, method="OM", subm="CONSTANT", indel=1, refseq=0, with_missing=True)
    assert result.shape == (2,) or result.shape == (2, 1)
    assert np.all(np.asarray(result) == 0)


def test_chi2_against_an_external_reference_reads_a_weights_keyword_as_the_objects_weights():
    external = tx.seqdef(
        pd.DataFrame({"a": ["A"], "b": ["B"], "c": ["C"]}), alphabet=["A", "B", "C"]
    )
    weights = [1.0, 2.0, 1.0, 1.0]
    by_keyword = tx.seqdist(tx.seqdef(FRAME), method="CHI2", refseq=external, weights=weights)
    by_object = tx.seqdist(
        tx.seqdef(FRAME, weights=weights), method="CHI2", refseq=external, weighted=True
    )
    assert np.allclose(by_keyword, by_object)


def test_dhd_on_a_plain_matrix_equals_dhd_on_the_object(seqs, codes):
    assert np.allclose(tx.seqdist(codes, method="DHD"), tx.seqdist(seqs, method="DHD"))


def test_generated_costs_for_an_sts_like_object_without_an_alphabet(seqs, codes):
    # seqdist accepts any non-ndarray object exposing .data (the historical STS-like contract);
    # its generated costs take the codes as the alphabet and the object's weights.
    weights = np.array([1.0, 2.0, 1.0, 1.0])
    like = types.SimpleNamespace(data=codes.copy(), weights=weights.copy())
    weighted = tx.seqdef(FRAME, weights=weights)
    assert np.allclose(
        tx.seqdist(like, method="OM", subm="INDELS"),
        tx.seqdist(weighted, method="OM", subm="INDELS", weighted=True),
    )


def test_seqdist_refuses_an_indel_that_is_neither_a_number_nor_a_vector(seqs):
    with pytest.raises(ValueError, match="'indel' invalid type or unknown value"):
        tx.seqdist(seqs, method="OM", subm="CONSTANT", indel=object())


# --- normalization and batched kernels ----------------------------------------------------------


def test_normalize_value_follows_traminex_distancecalculator():
    # DistanceCalculator::normalizeDistance (distancecalculator.h:45-70): maxlength with two empty
    # sequences gives 0; maxdist and YujianBo with a zero bound give 1; an unknown code leaves the
    # raw distance.
    assert _normalization.normalize_value(1.0, 5.0, 0, 0, 1) == 0.0
    assert _normalization.normalize_value(1.0, 0.0, 2, 2, 3) == 1.0
    assert _normalization.normalize_value(1.0, 0.0, 2, 2, 4) == 1.0
    assert _normalization.normalize_value(1.5, 2.0, 2, 2, 9) == 1.5


def test_batched_om_with_an_infinite_indel_separates_only_different_rows():
    rows = [np.array([0, 1]), np.array([0, 1]), np.array([1, 1])]
    raw = _batched.om_raw(
        rows, np.array([0, 0]), np.array([1, 2]), indel=np.inf, subm=np.ones((2, 2))
    )
    assert raw[0] == 0.0 and np.isposinf(raw[1])


# --- Python kernels: refusals and empty inputs --------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        (
            {"method": "OM", "subm": np.inf, "indel": 1.0},
            ValueError,
            "OM substitution costs must be finite",
        ),
        (
            {"method": "OM", "subm": np.zeros((2, 2)), "indel": 1.0},
            ValueError,
            r"need at least 3x3",
        ),
        (
            {"method": "OM", "subm": "x", "indel": 1.0},
            TypeError,
            "OM subm accepts float or ndarray",
        ),
        (
            {"method": "OM", "subm": 2.0, "indel": 1.0, "norm": "gmean"},
            ValueError,
            r"OM supports norm=None\|'length'\|'max'",
        ),
        (
            {"method": "OM", "subm": 2.0, "indel": -1.0},
            ValueError,
            "OM indel cost must be strictly positive",
        ),
        ({"method": "HAM", "norm": "gmean"}, ValueError, "HAM supports norm in"),
        ({"method": "LCS", "norm": "max"}, ValueError, r"LCS supports norm=None\|'length'"),
        ({"method": "LCP", "norm": "max"}, ValueError, r"LCP supports norm=None\|'length'"),
        ({"method": "RLCP", "norm": "max"}, ValueError, r"RLCP supports norm=None\|'length'"),
        (
            {"method": "DHD", "subm": np.zeros((3, 3, 2))},
            ValueError,
            "DHD subm must be square in its two state axes",
        ),
        (
            {"method": "DHD", "subm": np.array([object()]), "variant": "sum"},
            ValueError,
            r"DHD subm must be a 2D \(K,K\) or 3D",
        ),
    ],
)
def test_python_kernel_refusals(codes, kwargs, error, message):
    with pytest.raises(error, match=message):
        kernels.pairwise(codes, **kwargs)


def test_om_without_substitution_costs_uses_two(codes):
    assert np.array_equal(
        kernels.pairwise(codes, method="OM", subm=None, indel=1.0),
        kernels.pairwise(codes, method="OM", subm=2.0, indel=1.0),
    )


@pytest.mark.parametrize(
    ("indel", "error", "message"),
    [
        ("x", TypeError, "OMspell indel must be a scalar or state-cost vector"),
        (0, ValueError, "OMspell indel cost must be strictly positive"),
        ([1, 0, 1], ValueError, "OMspell indel cost must be strictly positive"),
    ],
)
def test_state_indel_refusals(indel, error, message):
    with pytest.raises(error, match=message):
        kernels._state_indel_spec(indel, n_states=3, method="OMspell")


def test_fixed_cost_refusals():
    with pytest.raises(ValueError, match="OMspell substitution costs cannot be negative"):
        kernels._refuse_invalid_fixed_costs("OMspell", -1.0, [np.array([0, 1])])
    with pytest.raises(ValueError, match=r"need at least 2x2"):
        kernels._refuse_invalid_fixed_costs("OMspell", np.zeros((1, 1)), [np.array([0, 1])])


def test_kernel_helpers_on_empty_inputs():
    assert kernels._subm_max_bound(np.zeros((0, 0))) == 0.0
    assert kernels._asymmetric_slice_flags(np.zeros(3)).size == 0
    assert kernels._dhd_max_cost_sum(np.zeros((0, 3, 3))) == 0.0
    empty = np.array([], dtype=int)
    assert kernels._in_batches(lambda a, b: np.zeros(len(a)), empty, empty, 4).size == 0


def test_explicit_dhd_costs_broadcast_a_matrix_and_refuse_other_shapes(codes):
    costs = np.ones((3, 3)) - np.eye(3)
    stationary = kernels._pairwise_dhd_explicit_subm(codes, costs)
    per_position = kernels._pairwise_dhd_explicit_subm(
        codes, np.broadcast_to(costs, (3, 3, 3)).copy()
    )
    assert np.array_equal(stationary, per_position)
    with pytest.raises(ValueError, match="time-varying DHD subm must have first axis L=3, got 2"):
        kernels._pairwise_dhd_explicit_subm(codes, np.zeros((2, 3, 3)))
    with pytest.raises(ValueError, match=r"DHD subm must be a 2D \(K,K\) or 3D"):
        kernels._pairwise_dhd_explicit_subm(codes, np.zeros((1, 3, 3, 3)))


def test_omspell_rows_infer_the_state_count_when_none_is_given():
    rows = [np.array([0, 1, 1]), np.array([1, 0])]
    options = {"indel": 1.0, "subm": 2.0}
    assert np.array_equal(
        kernels._pairwise_rows_omspell(rows, dict(options), n_states=None),
        kernels._pairwise_rows_omspell(rows, dict(options), n_states=2),
    )


def test_indels_costs_count_the_states_of_an_external_reference_with_a_void():
    base = pd.DataFrame({"a": ["A", "B", "A"], "b": ["B", "B", None], "c": ["C", None, None]})
    external = tx.seqdef(
        pd.DataFrame({"a": ["A"], "b": ["B"], "c": [None]}), alphabet=["A", "B", "C"]
    )
    distances = np.asarray(
        tx.seqdist(tx.seqdef(base), method="OM", subm="INDELS", refseq=external)
    ).ravel()
    assert np.allclose(distances, _numbers("indels_external_void"))
