"""Subsetting a SubSeqElist gives a SubSeqElist (TX-AUD-20261010-038).

R's ``[.subseqelist`` (subseqelist.R:28-45) returns a subseqelist of the same
class -- subseq and data rows subset, labels and bonferroni kept -- so
``seqecmpgroup(fs[1:2], ...)`` works. traminex returned the label strings, so
``seqecmpgroup(fs[0:2], ...)`` raised "subseq should be a subseqelist", and a
list of positions or a boolean mask raised a raw TypeError. An integer index
still returns the label: Python's scalar-index convention (ledger
``[.subseqelist`` row). Fixture:
``scripts/r_goldens/make_goldens_subseqelist_subset.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx
from traminex.events import SubSeqElist, seqecmpgroup, seqefsub

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "subseqelist_subset.csv"


def _golden() -> dict[tuple[str, str], list[str]]:
    out: dict[tuple[str, str], list[str]] = defaultdict(list)
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            out[(row["case"], row["kind"])].append(row["value"])
    return out


R = _golden()
EVENTS = tx.seqecreate(
    id=[1, 1, 2, 2, 3, 3, 4, 4],
    timestamp=[1, 2, 1, 3, 2, 4, 1, 2],
    event=["A", "B", "A", "C", "B", "C", "A", "B"],
)
FS = seqefsub(EVENTS, pmin_support=0.25)
SUBSETS = {
    "slice_1_2": slice(0, 2),
    "index_1_3": [0, 2],
    "single_2": [1],
    "mask_1_3_5": [True, False, True, False, True, False],
}


def test_full_list_matches_r() -> None:
    assert FS.subseq == R[("full", "subseq")]


@pytest.mark.parametrize("case", sorted(SUBSETS))
def test_subset_is_a_subseqelist_matching_r(case: str) -> None:
    sub = FS[SUBSETS[case]]
    assert type(sub) is SubSeqElist
    assert sub.subseq == R[(case, "subseq")]
    np.testing.assert_allclose(
        sub.data["Support"].to_numpy(float),
        [float(v) for v in R[(case, "Support")]],
        atol=1e-12,
        rtol=0,
    )
    np.testing.assert_allclose(
        sub.data["Count"].to_numpy(float),
        [float(v) for v in R[(case, "Count")]],
        atol=1e-12,
        rtol=0,
    )
    assert sub.eseq is FS.eseq and sub.constraint is FS.constraint and sub.type == FS.type


def test_numpy_index_array_is_accepted() -> None:
    assert FS[np.array([0, 2])].subseq == R[("index_1_3", "subseq")]


@pytest.mark.parametrize(
    ("case", "index", "group"),
    [
        ("cmp_slice_1_2", slice(0, 2), ["g1", "g1", "g2", "g2"]),
        ("cmp_index_1_3", [0, 2], ["g1", "g2", "g1", "g2"]),
    ],
)
def test_seqecmpgroup_on_a_subset_matches_r(case: str, index, group: list[str]) -> None:
    result = seqecmpgroup(FS[index], group=group)
    assert result.subseq == R[(case, "subseq")]
    for column in ("Support", "p.value", "statistic", "index", "Freq.g1", "Freq.g2"):
        np.testing.assert_allclose(
            result.data[column].to_numpy(float),
            [float(v) for v in R[(case, column)]],
            rtol=0,
            atol=R_GOLDEN_ATOL,
            err_msg=column,
        )
    for column in ("Resid.g1", "Resid.g2"):
        np.testing.assert_allclose(
            result.data[column].to_numpy(float),
            [float(v) for v in R[(case, column)]],
            rtol=0,
            atol=1e-12,
            err_msg=column,
        )


def test_chisq_subset_keeps_its_class_labels_and_bonferroni() -> None:
    whole = seqecmpgroup(FS, group=["g1", "g1", "g2", "g2"])
    part = whole[0:2]
    assert type(part) is type(whole)
    assert R[("chisq_slice_1_2", "class")] == ["subseqelistchisq", "subseqelist"]
    assert part.subseq == R[("chisq_slice_1_2", "subseq")]
    assert part.labels == R[("chisq_slice_1_2", "labels")] == ["g1", "g2"]
    assert part.bonferroni == whole.bonferroni == {"used": False, "ntest": 6}
    assert R[("chisq_slice_1_2", "bonferroni")] == ["0", "6"]


def test_integer_index_returns_the_label() -> None:
    assert FS[1] == "(B)"
    assert list(FS) == FS.subseq
