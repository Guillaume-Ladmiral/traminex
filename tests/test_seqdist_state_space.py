"""CHI2/EUCLID and the NMS family are sized on the alphabet, not the observed codes.

R sizes ``global.pdotj`` and ``prox`` by the alphabet (plus the missing state
under ``with.missing``). traminex sized them by the codes the corpus happened to
use on its dense path, and resolved ``global_pdotj="obs"`` from ``seqmeant``
rows that include the void marker: an object whose last state is unobserved, a
ragged object, and a plain integer matrix all refused ``"obs"``, while a vector
or ``prox`` sized to the observed codes only -- which R refuses -- was accepted
(2026-09-22 audit S1a-01, S1a-04, S1b-04, S1b-05). A plain matrix's alphabet is
``0..max``, so its analogue is ``seqdef(X, alphabet = 0:max)``.

Expected values: live TraMineR 2.2-14 via
``scripts/r_goldens/make_goldens_seqdist_state_space.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens"
X = np.array([[0, 1, 2, 1], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 2, 0], [1, 0, 0, 1]])
P4 = np.array([[1, 0.5, 0, 0.2], [0.5, 1, 0.3, 0.2], [0, 0.3, 1, 0.2], [0.2, 0.2, 0.2, 1.0]])


def _matrices() -> dict[str, np.ndarray]:
    cells: dict[str, dict] = defaultdict(dict)
    with (GOLDENS / "seqdist_state_space_cases.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            cells[row["case"]][(int(row["row"]), int(row["col"]))] = float(row["value"])
    out = {}
    for case, values in cells.items():
        size = max(i for i, _ in values) + 1
        matrix = np.empty((size, size))
        for (i, j), value in values.items():
            matrix[i, j] = value
        out[case] = matrix
    return out


R = _matrices()
REFUSALS = dict(
    line.split("\t", 1)
    for line in (GOLDENS / "seqdist_state_space_refusals.txt").read_text().splitlines()
    if not line.startswith("#")
)


def _unobserved_last():
    return tx.seqdef(X, alphabet=[0, 1, 2, 3])


def _ragged():
    rows = X.astype(object)
    rows[0, 3] = None
    rows[2, 2:] = None
    return tx.seqdef(rows, alphabet=[0, 1, 2])


def _with_missing():
    rows = X.astype(object)
    rows[1, 2] = None
    return tx.seqdef(rows, alphabet=[0, 1, 2], right=None)


def _letters_dense():
    data = [["A", "A", "B", "C"], ["A", "B", "B", "C"], ["C", "C", "A", "B"]]
    return tx.seqdef(np.array(data, dtype=object), alphabet=["A", "B", "C", "D"])


CASES = {
    "chi2_obs_plain": lambda: tx.seqdist(X, method="CHI2", global_pdotj="obs"),
    "chi2_obs_plain_auto": lambda: tx.seqdist(X, method="CHI2", global_pdotj="obs", norm="auto"),
    "euclid_obs_plain": lambda: tx.seqdist(X, method="EUCLID", global_pdotj="obs"),
    "chi2_obs_unobserved_last": lambda: tx.seqdist(
        _unobserved_last(), method="CHI2", global_pdotj="obs"
    ),
    "chi2_vec4_unobserved_last": lambda: tx.seqdist(
        _unobserved_last(), method="CHI2", global_pdotj=[0.2, 0.3, 0.4, 0.1]
    ),
    "euclid_obs_unobserved_last": lambda: tx.seqdist(
        _unobserved_last(), method="EUCLID", global_pdotj="obs"
    ),
    "chi2_obs_ragged": lambda: tx.seqdist(_ragged(), method="CHI2", global_pdotj="obs"),
    "euclid_obs_ragged": lambda: tx.seqdist(_ragged(), method="EUCLID", global_pdotj="obs"),
    "chi2_obs_with_missing": lambda: tx.seqdist(
        _with_missing(), method="CHI2", global_pdotj="obs", with_missing=True
    ),
    "nms_prox4_unobserved_last": lambda: tx.seqdist(_letters_dense(), method="NMS", prox=P4),
    "svrspell_prox4_unobserved_last": lambda: tx.seqdist(
        _letters_dense(), method="SVRSPELL", prox=P4
    ),
}


def test_the_fixture_covers_every_case():
    assert set(R) == set(CASES)


@pytest.mark.parametrize("case", sorted(CASES))
def test_the_state_space_is_the_alphabet(case):
    np.testing.assert_allclose(
        np.asarray(CASES[case](), dtype=float), R[case], rtol=1e-12, atol=1e-12
    )


@pytest.mark.parametrize(
    "case,call,message",
    [
        (
            "chi2_vec3_unobserved_last",
            lambda: tx.seqdist(_unobserved_last(), method="CHI2", global_pdotj=[0.2, 0.3, 0.5]),
            "global_pdotj length must equal",
        ),
        (
            "nms_prox3_unobserved_last",
            lambda: tx.seqdist(_letters_dense(), method="NMS", prox=P4[:3, :3]),
            "prox",
        ),
        (
            "svrspell_prox3_unobserved_last",
            lambda: tx.seqdist(_letters_dense(), method="SVRSPELL", prox=P4[:3, :3]),
            "prox",
        ),
    ],
)
def test_a_vector_or_matrix_sized_to_the_observed_codes_is_refused(case, call, message):
    assert case in REFUSALS, "R must refuse this input"
    with pytest.raises(ValueError, match=message):
        call()
