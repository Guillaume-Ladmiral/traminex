"""Parity tests for OMstran (transition-augmented OM) vs TraMineR seqdist-OMstran.R.

Key invariant: with otto=1 the transition penalty vanishes (transweight=0) and tokens are
compared only on their origin state, so OMstran reduces to OM on the original sequence with
costs rescaled by indelrate = max(indel)/max(sm). For indel=1, subm=2 that factor is 0.5.
"""

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_c
from r_goldens_utils import R_GOLDENS_DIR

import traminex as tx
from traminex.distances._omstran import omstran_transform


def test_omstran_otto1_reduces_to_scaled_om():
    seqmat = np.array(
        [[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0], [0, 1, 0, 1, 2]], dtype=np.int32
    )
    D_om = tx.seqdist(seqmat, method="OM", indel=1.0, subm=2.0)
    D_stran = tx.seqdist(seqmat, method="OMstran", otto=1.0, indel=1.0, subm=2.0)
    assert np.allclose(D_stran, 0.5 * D_om)


def test_omstran_transform_token_count():
    seqmat = np.array([[0, 1, 2], [2, 1, 0]], dtype=np.int32)
    new_seqmat, indels, new_sm = omstran_transform(seqmat, otto=0.5)
    # add_column=True keeps L tokens per row
    assert new_seqmat.shape == (2, 3)
    assert new_sm.shape[0] == new_sm.shape[1] == len(indels)
    assert np.allclose(new_sm, new_sm.T)
    assert np.allclose(np.diag(new_sm), 0.0)


@pytest.mark.parametrize("transindel", ["constant", "prob", "subcost"])
def test_omstran_properties(transindel):
    seqmat = np.array(
        [[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0], [0, 0, 0, 1, 2]], dtype=np.int32
    )
    D = tx.seqdist(seqmat, method="OMstran", otto=0.3, transindel=transindel, indel=1.0, subm=2.0)
    assert D.shape == (4, 4)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert np.all(D >= -1e-9)


def test_omstran_previous_variant_runs():
    seqmat = np.array([[0, 1, 2, 1], [2, 1, 0, 1], [0, 0, 1, 2]], dtype=np.int32)
    D = tx.seqdist(seqmat, method="OMstran", otto=0.5, previous=True, indel=1.0, subm=2.0)
    assert D.shape == (3, 3)
    assert np.allclose(D, D.T)
    assert np.allclose(np.diag(D), 0.0)


def test_omstran_requires_otto():
    seqmat = np.array([[0, 1, 2], [2, 1, 0]], dtype=np.int32)
    with pytest.raises((ValueError, TypeError)):
        tx.seqdist(seqmat, method="OMstran", indel=1.0, subm=2.0)


#: R's d(1, 2) for these cases, recorded by scripts/r_goldens/make_goldens_distances_void.R
#: (TX-AUD-20261019-088: the values were copied here by hand, outside the provenance gate).
_R1_CASES = pd.read_csv(R_GOLDENS_DIR / "distances_r1_cases.csv")


def _r_first_pair(case: str) -> float:
    rows = _R1_CASES[
        (_R1_CASES["case"] == case) & (_R1_CASES["row"] == 0) & (_R1_CASES["col"] == 1)
    ]
    assert len(rows) == 1, case
    return float(rows["value"].iloc[0])


def _option_case(transindel: str, previous: bool, add_column: bool) -> str:
    if (transindel, previous, add_column) == ("constant", False, True):
        return "ragged_omstran"  # the generator's default-option case
    return (
        f"ragged_omstran_{transindel}_previous_{str(previous).lower()}"
        f"_add_column_{str(add_column).lower()}"
    )


@pytest.mark.parametrize("transindel", ["constant", "prob", "subcost"])
@pytest.mark.parametrize("previous", [False, True])
@pytest.mark.parametrize("add_column", [False, True])
def test_ragged_omstran_option_matrix_matches_r(transindel, previous, add_column):
    expected = _r_first_pair(_option_case(transindel, previous, add_column))
    actual = tx.seqdist(
        fixture_c(),
        method="OMstran",
        with_missing=True,
        norm="none",
        subm=2.0,
        indel=1.0,
        otto=0.5,
        transindel=transindel,
        previous=previous,
        add_column=add_column,
    )

    assert actual[0, 1] == pytest.approx(expected, rel=0.0, abs=1e-12)


def test_ragged_omstran_probability_rates_use_attached_weights():
    actual = tx.seqdist(
        fixture_c(weighted=True),
        method="OMstran",
        with_missing=True,
        norm="none",
        subm=2.0,
        indel=1.0,
        otto=0.5,
        transindel="prob",
    )

    expected = _r_first_pair("ragged_omstran_prob_attached_weighted")
    assert actual[0, 1] == pytest.approx(expected, rel=0.0, abs=1e-12)
