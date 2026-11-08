"""Reference distances with asymmetric costs are d(sequence, reference) (TX-REP-20261017-001).

With an asymmetric substitution matrix d(x, y) and d(y, x) differ. R gives a
reference distance with the sequence first, and index sets as d(set-1 element,
set-2 element). traminex read reference distances off its mirrored full matrix,
so for every sequence after the reference (or a first set after the second) it
gave d(reference, sequence) instead. A normalized distance divides each
oriented cell by R's bound; the normalizers used to read one triangle and
returned 0 after the reference (batch-6 review C-1). Fixture
``scripts/r_goldens/make_goldens_refseq_asymmetric.R``.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_DIFFERENCE_ATOL

import traminex as tx
from traminex.distances import pairwise, pairwise_between_sets
from traminex.sequences import STSList

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "refseq_asymmetric.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
ASYM = np.array([[0, 1, 3], [2, 0, 1], [1, 4, 0]], dtype=float)
# R's array(asym, c(3, 3, 3)): the same costs at every position.
ARR = np.broadcast_to(ASYM, (3, 3, 3)).copy()
D4 = tx.seqdef(
    np.array([list("ABC"), list("BCA"), list("CAB"), list("AAB")], dtype=object),
    alphabet=["A", "B", "C"],
)
D5 = tx.seqdef(
    np.array([list("ABC"), list("BCA"), list("ABC"), list("CAB"), list("BCA")], dtype=object),
    alphabet=["A", "B", "C"],
)


def _row(seqs: STSList, index: int) -> STSList:
    return STSList(
        data=seqs.data[index : index + 1].copy(),
        alphabet=list(seqs.alphabet),
        nr_code=seqs.nr_code,
        void_code=seqs.void_code,
    )


def _om(seqs, **kw):
    return tx.seqdist(seqs, method="OM", subm=ASYM, indel=2, **kw)


def _ham(seqs, **kw):
    return tx.seqdist(seqs, method="HAM", subm=ASYM, **kw)


def _dhd(seqs, **kw):
    return tx.seqdist(seqs, method="DHD", subm=ARR, **kw)


SM3 = np.array([[0, 1, 5], [3, 0, 1], [1, 4, 0]], dtype=float)
DUP = tx.seqdef([list("AAB"), list("CCC"), list("AAB")], alphabet=["A", "B", "C"])
VD = tx.seqdef([list("AAB"), list("BCC"), list("CCCC")], alphabet=["A", "B", "C"], right="DEL")
DUP_CODES = np.asarray(DUP.data, dtype=np.int32)
VD_CODES = np.asarray(VD.data, dtype=np.int32)
VD_ROWS = [row[row != VD.void_code] for row in VD_CODES]


def _om3(seqs, indel=10.0, **kw):
    return tx.seqdist(seqs, method="OM", subm=SM3, indel=indel, **kw)


#: Cases whose R value is an R defect corrected by traminex, with the case of their raw distances.
R_DEFECT_CORRECTED = {
    "d4_om_gmean_refseq_1": "d4_om_refseq_1",
    "d4_om_maxdist_refseq_3": "d4_om_refseq_3",
}

CASES = {
    "d4_om_full": lambda: _om(D4),
    **{f"d4_om_refseq_{k}": (lambda k=k: _om(D4, refseq=k)) for k in (1, 2, 3, 4, 0)},
    "d4_om_sets_12_34": lambda: _om(D4, refseq=([1, 2], [3, 4])),
    "d4_om_sets_34_12": lambda: _om(D4, refseq=([3, 4], [1, 2])),
    "d4_om_external_2": lambda: _om(D4, refseq=_row(D4, 1)),
    "d4_ham_refseq_1": lambda: _ham(D4, refseq=1),
    "d4_ham_refseq_3": lambda: _ham(D4, refseq=3),
    "d4_ham_sets_34_12": lambda: _ham(D4, refseq=([3, 4], [1, 2])),
    "d4_ham_external_2": lambda: _ham(D4, refseq=_row(D4, 1)),
    "d5_om_full": lambda: _om(D5),
    **{f"d5_om_refseq_{k}": (lambda k=k: _om(D5, refseq=k)) for k in (1, 2, 4, 5, 0)},
    "d5_om_external_2": lambda: _om(D5, refseq=_row(D5, 1)),
    "d5_om_sets_45_12": lambda: _om(D5, refseq=([4, 5], [1, 2])),
    **{f"d4_om_auto_refseq_{k}": (lambda k=k: _om(D4, norm="auto", refseq=k)) for k in (1, 3, 0)},
    "d4_om_auto_sets_34_12": lambda: _om(D4, norm="auto", refseq=([3, 4], [1, 2])),
    "d4_om_auto_external_2": lambda: _om(D4, norm="auto", refseq=_row(D4, 1)),
    "d4_om_maxlength_refseq_1": lambda: _om(D4, norm="maxlength", refseq=1),
    "d4_om_gmean_refseq_1": lambda: _om(D4, norm="gmean", refseq=1),
    "d4_om_maxdist_refseq_3": lambda: _om(D4, norm="maxdist", refseq=3),
    "d4_om_yujianbo_refseq_1": lambda: _om(D4, norm="YujianBo", refseq=1),
    "d4_ham_auto_refseq_3": lambda: _ham(D4, norm="auto", refseq=3),
    "d4_ham_auto_sets_34_12": lambda: _ham(D4, norm="auto", refseq=([3, 4], [1, 2])),
    "d4_dhd_refseq_3": lambda: _dhd(D4, refseq=3),
    "d4_dhd_auto_refseq_3": lambda: _dhd(D4, norm="auto", refseq=3),
    "d4_dhd_auto_refseq_1": lambda: _dhd(D4, norm="auto", refseq=1),
    "d5_om_auto_refseq_2": lambda: _om(D5, norm="auto", refseq=2),
    "d5_om_auto_refseq_0": lambda: _om(D5, norm="auto", refseq=0),
    "dup_om_sets_12_3": lambda: _om3(DUP, refseq=([1, 2], [3])),
    "dup_om_sets_3_2": lambda: _om3(DUP, refseq=([3], [2])),
    "dup_om_sets_12_23": lambda: _om3(DUP, refseq=([1, 2], [2, 3])),
    "void_om_sets_3_1": lambda: _om3(VD, refseq=([3], [1])),
    "void_om_ivec_sets_3_1": lambda: _om3(VD, indel=np.array([10.0, 8.0, 12.0]), refseq=([3], [1])),
    "void_om_sets_12_23": lambda: _om3(VD, refseq=([1, 2], [2, 3])),
}


def test_every_r_case_is_rebuilt() -> None:
    assert set(CASES) == set(GOLDEN.index)


REFERENCE_CASES = {
    "d4_om_refseq_3": lambda **kw: _om(D4, refseq=3, **kw),
    "d4_om_refseq_0": lambda **kw: _om(D4, refseq=0, **kw),
    "d4_om_sets_34_12": lambda **kw: _om(D4, refseq=([3, 4], [1, 2]), **kw),
    "d4_om_external_2": lambda **kw: _om(D4, refseq=_row(D4, 1), **kw),
    "d4_om_auto_refseq_3": lambda **kw: _om(D4, norm="auto", refseq=3, **kw),
    "d4_ham_refseq_3": lambda **kw: _ham(D4, refseq=3, **kw),
    "d5_om_refseq_2": lambda **kw: _om(D5, refseq=2, **kw),
    "d5_om_refseq_4": lambda **kw: _om(D5, refseq=4, **kw),
    "d5_om_refseq_0": lambda **kw: _om(D5, refseq=0, **kw),
    "d5_om_sets_45_12": lambda **kw: _om(D5, refseq=([4, 5], [1, 2]), **kw),
}


@pytest.mark.parametrize("case", sorted(REFERENCE_CASES))
def test_deduplicate_keeps_r_direction(case: str) -> None:
    # deduplicate is a traminex computation option and must not change a value;
    # the index routes passed it to their second call, which then computed over
    # the distinct rows in first-occurrence order (batch-6 review I-2).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(REFERENCE_CASES[case](deduplicate=True), dtype=float)
    expected = np.array([float(v) for v in GOLDEN.loc[case, "value"].split()])
    np.testing.assert_allclose(actual.ravel(order="F"), expected, atol=1e-12, rtol=0)


@pytest.mark.parametrize("case", sorted(CASES))
def test_reference_distances_match_r(case: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # R's asymmetry and triangle messages
        actual = np.asarray(CASES[case](), dtype=float)
    expected = np.array([float(v) for v in GOLDEN.loc[case, "value"].split()])
    raw_case = R_DEFECT_CORRECTED.get(case)
    if raw_case is not None:
        # The bound takes ASYM's whole off-diagonal maximum, min(4, 2 * 2) = 4, so it is 4 * 3 = 12
        # (R's upper-triangle 3 gives 9); with p = q = 2 * 3, 2 * sqrt(p * q) = 12 too. From R's raw
        # distances d: maxdist = d / 12 and gmean = 1 - (12 - d) / 12 = d / 12
        # (TX-AUD-20261019-008, an R defect corrected; R's values stay recorded).
        raw = np.array([float(v) for v in GOLDEN.loc[raw_case, "value"].split()])
        np.testing.assert_allclose(actual.ravel(order="F"), raw / 12.0, atol=1e-12, rtol=0)
        assert not np.allclose(actual.ravel(order="F"), expected, rtol=0, atol=R_DIFFERENCE_ATOL)
        return
    np.testing.assert_allclose(actual.ravel(order="F"), expected, atol=1e-12, rtol=0)


_OM10 = {"subm": SM3, "indel": 10.0}
#: The low-level routes of each R case: pairwise(ref_sets=) (zero-based) and
#: pairwise_between_sets, which compute d(set-1 sequence, set-2 sequence) as R does
#: (TX-AUD-20261019-009; ruling R6-j covered seqdist only).
LOW_LEVEL = {
    "dup_om_sets_12_3": [
        lambda **kw: pairwise(DUP_CODES, "OM", ref_sets=([0, 1], [2]), **_OM10, **kw),
        lambda **kw: pairwise_between_sets(
            DUP_CODES[[0, 1]], DUP_CODES[[2]], method="OM", **_OM10, **kw
        ),
    ],
    "dup_om_sets_3_2": [
        lambda **kw: pairwise(DUP_CODES, "OM", ref_sets=([2], [1]), **_OM10, **kw),
        lambda **kw: pairwise_between_sets(
            DUP_CODES[[2]], DUP_CODES[[1]], method="OM", **_OM10, **kw
        ),
    ],
    "dup_om_sets_12_23": [
        lambda **kw: pairwise(DUP_CODES, "OM", ref_sets=([0, 1], [1, 2]), **_OM10, **kw),
        lambda **kw: pairwise_between_sets(
            DUP_CODES[[0, 1]], DUP_CODES[[1, 2]], method="OM", **_OM10, **kw
        ),
    ],
    "void_om_sets_3_1": [
        lambda **kw: pairwise(
            VD_CODES, "OM", void=VD.void_code, ref_sets=([2], [0]), **_OM10, **kw
        ),
        lambda **kw: pairwise_between_sets(VD_ROWS[2:3], VD_ROWS[0:1], method="OM", **_OM10, **kw),
    ],
    "void_om_ivec_sets_3_1": [
        lambda **kw: pairwise(
            VD_CODES,
            "OM",
            subm=SM3,
            indel=np.array([10.0, 8.0, 12.0]),
            void=VD.void_code,
            ref_sets=([2], [0]),
            **kw,
        ),
    ],
    "void_om_sets_12_23": [
        lambda **kw: pairwise(
            VD_CODES, "OM", void=VD.void_code, ref_sets=([0, 1], [1, 2]), **_OM10, **kw
        ),
        lambda **kw: pairwise_between_sets(VD_ROWS[0:2], VD_ROWS[1:3], method="OM", **_OM10, **kw),
    ],
}


@pytest.mark.parametrize("deduplicate", [False, True])
@pytest.mark.parametrize(
    "case, route",
    [(case, k) for case, calls in sorted(LOW_LEVEL.items()) for k in range(len(calls))],
)
def test_low_level_reference_routes_match_r(case: str, route: int, deduplicate: bool) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(LOW_LEVEL[case][route](deduplicate=deduplicate), dtype=float)
    expected = np.array([float(v) for v in GOLDEN.loc[case, "value"].split()])
    np.testing.assert_allclose(actual.ravel(order="F"), expected, atol=1e-12, rtol=0)
