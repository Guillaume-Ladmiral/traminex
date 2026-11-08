"""Caller-supplied costs are checked on every OM-family and TWED path, as R does.

``checkcost`` (``TraMineR-checkcost.R``) stops on a wrong-sized substitution
matrix, a negative cost and a non-zero diagonal for OM, OMloc, OMslen, OMspell,
OMstran and TWED, on regular and ragged data alike. traminex checked only its
dense OM path: every other method, and OM on ragged data, computed -- negative
distances, and non-zero distances between identical sequences (2026-09-22 audit
S1b-01). R's outcome for each case is in ``tests/r_goldens/seqdist_cost_checks.txt``
(``scripts/r_goldens/make_goldens_seqdist_cost_checks.R``).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise

GOLDENS = Path(__file__).parent / "r_goldens"
X = np.array([[0, 1, 2, 1], [1, 1, 0, 2], [2, 0, 1, 1]])


def _data(shape: str) -> tx.STSList:
    if shape == "regular":
        return tx.seqdef(X, alphabet=[0, 1, 2])
    rows = X.astype(object)
    rows[2, 3] = None
    return tx.seqdef(rows, alphabet=[0, 1, 2])


def _defect(name: str) -> np.ndarray:
    base = np.full((3, 3), 2.0)
    np.fill_diagonal(base, 0.0)
    if name == "diagonal":
        base[1, 1] = 1.0
    elif name == "negative":
        base[0, 2] = base[2, 0] = -1.0
    elif name == "oversized":
        base = np.full((4, 4), 2.0)
        np.fill_diagonal(base, 0.0)
    elif name == "undersized":
        base = np.full((2, 2), 2.0)
        np.fill_diagonal(base, 0.0)
    return base


METHOD_ARGS = {
    "OM": ("OM", {"indel": 1.0}),
    "OMloc": ("OMLOC", {}),
    "OMslen": ("OMSLEN", {"indel": 1.0}),
    "OMspell": ("OMSPELL", {"indel": 1.0}),
    "OMstran": ("OMSTRAN", {"indel": 1.0, "otto": 0.5}),
    "TWED": ("TWED", {"indel": 1.0, "nu": 0.5, "h": 0.5}),
}

R_OUTCOMES = [
    tuple(line.split("\t"))
    for line in (GOLDENS / "seqdist_cost_checks.txt").read_text().splitlines()
    if not line.startswith("#")
]


def test_r_refuses_every_defective_matrix():
    assert len(R_OUTCOMES) == 2 * len(METHOD_ARGS) * 4
    assert all(outcome.startswith("error:") for *_, outcome in R_OUTCOMES)


@pytest.mark.parametrize("shape,r_method,defect,r_outcome", R_OUTCOMES, ids=lambda v: str(v)[:40])
def test_a_defective_substitution_matrix_is_refused(shape, r_method, defect, r_outcome):
    method, extra = METHOD_ARGS[r_method]
    # checkcost stops before its triangle-inequality diagnostic, so no cost
    # message may precede the refusal (warnings raised as errors here).
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match="subm|substitution"):
            tx.seqdist(_data(shape), method=method, subm=_defect(defect), **extra)


def _r_indel_results() -> dict[tuple[str, float], np.ndarray]:
    out = {}
    for line in (GOLDENS / "seqdist_cost_checks_indel.txt").read_text().splitlines():
        if line.startswith("#"):
            continue
        shape, indel, values = line.split("\t")
        out[(shape, float(indel))] = np.array(values.split(","), dtype=float).reshape(3, 3).T
    return out


@pytest.mark.parametrize("shape", ["regular", "ragged"])
@pytest.mark.parametrize("indel", [-1.0, 0.0])
def test_a_non_positive_indel_is_refused_on_every_path(shape, indel):
    """R never checks the indel for a cost-method string and returns negative
    distances or an all-zero matrix; traminex refuses on both paths (the dense
    path always did; the ragged one returned -8 between identical sequences)."""
    r_result = _r_indel_results()[(shape, indel)]
    assert (r_result.min() < 0) if indel < 0 else np.all(r_result == 0)

    with pytest.raises(ValueError, match="indel cost must be strictly positive"):
        tx.seqdist(_data(shape), method="OM", subm="CONSTANT", indel=indel)


@pytest.mark.parametrize("method", ["OM", "OMLOC", "OMSLEN", "OMSPELL", "TWED"])
def test_the_low_level_pairwise_refuses_the_same_matrices(method):
    """``pairwise`` reaches the row engines (void rows) and the dense branches
    without ``seqdist``; both refuse a negative cost and a non-zero diagonal."""
    extra = {"nu": 0.5, "h": 0.5} if method == "TWED" else {}
    if method != "OMLOC":
        extra["indel"] = 1.0
    ragged = np.array([[0, 1, 2, 1], [1, 1, 0, 2], [2, 0, 1, 9]])
    for defect in ("negative", "diagonal"):
        with pytest.raises(ValueError, match="substitution"):
            pairwise(ragged, method=method, subm=_defect(defect), void=9, **extra)
        if method != "OM":
            with pytest.raises(ValueError, match="substitution"):
                pairwise(X, method=method, subm=_defect(defect), **extra)
