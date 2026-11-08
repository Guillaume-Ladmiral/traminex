"""``seqdist`` must apply R's ``with.missing`` reset to plain integer arrays too.

``seqdist.R:111-118`` resets ``with.missing`` to ``FALSE`` when neither
``seqdata`` nor a sequence-type ``refseq`` contains a missing value, and stops
when either does without ``with.missing``.  traminex applied that rule only on
its sequence-object path: a plain integer array kept ``with_missing=True``, so
``seqcost`` (which since TraMineR 2.2-14 no longer downgrades the flag itself)
widened the generated costs by a missing state nothing could reach.  The OM
family then either refused the resulting indel vector (``ValueError: OM indel
vector length must equal the active state count``) or silently returned
different distances (OMslen, OMspell, OMloc, OMstran, TWED).  See
``ISSUE_SEQDIST_WITH_MISSING_RAW_ARRAY_20260921.md``.

A plain array carries no missing-state code, so R's nearest analogue is
``seqdef(X, alphabet = 0:2)`` with no missing cell.  Every expected value is
live TraMineR 2.2-14, written by
``scripts/r_goldens/make_goldens_seqdist_with_missing_reset.R``.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens"

X = np.array([[0, 1, 2, 1], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 2, 0], [1, 0, 0, 1]])

#: traminex spelling of each R method name, plus the arguments R requires.
METHODS = {
    "OM": ("OM", {}),
    "OMslen": ("OMSLEN", {}),
    "OMspell": ("OMSPELL", {}),
    "OMloc": ("OMLOC", {}),
    "OMstran": ("OMSTRAN", {"otto": 0.5}),
    "TWED": ("TWED", {"nu": 0.5, "h": 0.5}),
    "HAM": ("HAM", {}),
}

#: The recorded values carry 17 significant digits; the two implementations
#: accumulate the same sums in a different order, which moves the last bits.
TOLERANCE = {"rtol": 1e-12, "atol": 1e-12}


def _load_matrices() -> dict[tuple[str, str, str, bool], np.ndarray]:
    cells: dict[tuple[str, str, str, bool], dict[tuple[int, int], float]] = defaultdict(dict)
    with (GOLDENS / "seqdist_with_missing_reset.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["method"], row["sm"], row["indel"], row["with_missing"] == "TRUE")
            cells[key][(int(row["row"]), int(row["col"]))] = float(row["value"])
    matrices = {}
    for key, values in cells.items():
        size = max(i for i, _ in values) + 1
        matrix = np.empty((size, size))
        for (i, j), value in values.items():
            matrix[i, j] = value
        matrices[key] = matrix
    return matrices


R_MATRICES = _load_matrices()


def _sequence_object() -> tx.STSList:
    sequences = tx.seqdef(X, alphabet=[0, 1, 2])
    assert np.array_equal(np.asarray(sequences.data), X), "codes must stay 0, 1, 2"
    return sequences


def _inputs():
    return {"plain array": X.copy(), "sequence object": _sequence_object()}


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ seqdist.\ 'with\.missing'\ set\ as\ FALSE\ as\ 'seqdata'\ has\ no\ non\-void\ missing\ values:UserWarning"
)
def test_the_issue_reproduction_returns_the_r_matrix():
    """The exact call from the issue note: a plain array, INDELSLOG, with_missing=True."""
    expected = R_MATRICES[("OM", "INDELSLOG", "auto", True)]

    result = tx.seqdist(X, method="OM", subm="INDELSLOG", with_missing=True)

    np.testing.assert_allclose(np.asarray(result, dtype=float), expected, rtol=0, atol=1e-12)


def test_r_itself_resets_with_missing_on_missing_free_data():
    """Guard the fixture: R's TRUE and FALSE results coincide in every configuration."""
    for (method, sm, indel, with_missing), matrix in R_MATRICES.items():
        if with_missing:
            np.testing.assert_array_equal(matrix, R_MATRICES[(method, sm, indel, False)])


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ seqdist.\ 'with\.missing'\ set\ as\ FALSE\ as\ 'seqdata'\ has\ no\ non\-void\ missing\ values:UserWarning"
)
@pytest.mark.parametrize("input_kind", ["plain array", "sequence object"])
@pytest.mark.parametrize(
    "r_method,sm,indel,with_missing",
    sorted(R_MATRICES),
    ids=lambda value: str(value),
)
def test_every_generated_cost_method_matches_r(r_method, sm, indel, with_missing, input_kind):
    method, extra = METHODS[r_method]
    kwargs = dict(extra, subm=sm, with_missing=with_missing)
    if indel != "auto":
        kwargs["indel"] = [float(value) for value in indel.split(";")]

    if r_method == "TWED" and indel != "auto":
        # seqdist.R:301-305 with msg.warn's " [!!]" prefix (TX-REP-20261023-008).
        twed_vector = (
            " [!!] indel vector not supported by the chosen method, max(indel) used instead!"
        )
        with pytest.warns(UserWarning, match=re.escape(twed_vector)):
            result = tx.seqdist(_inputs()[input_kind], method=method, **kwargs)
    else:
        result = tx.seqdist(_inputs()[input_kind], method=method, **kwargs)

    np.testing.assert_allclose(
        np.asarray(result, dtype=float),
        R_MATRICES[(r_method, sm, indel, with_missing)],
        **TOLERANCE,
    )


def _external_reference_with_missing() -> tx.STSList:
    reference = tx.seqdef(np.array([[0, 1, None, 2]], dtype=object), alphabet=[0, 1, 2])
    assert reference.nr_code is not None and np.any(reference.data == reference.nr_code)
    return reference


def _load_reference_distances() -> dict[str, np.ndarray]:
    values: dict[str, list[tuple[int, float]]] = defaultdict(list)
    with (GOLDENS / "seqdist_with_missing_reset_refseq.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["with_missing"] == "TRUE":
                values[row["sm"]].append((int(row["row"]), float(row["value"])))
            else:
                assert row["row"] == "error", "R must refuse with.missing=FALSE here"
    return {sm: np.array([value for _, value in sorted(pairs)]) for sm, pairs in values.items()}


R_REFERENCE_DISTANCES = _load_reference_distances()


@pytest.mark.parametrize("input_kind", ["plain array", "sequence object"])
@pytest.mark.parametrize("sm", sorted(R_REFERENCE_DISTANCES))
def test_a_missing_bearing_reference_widens_the_costs(sm, input_kind):
    """R's reset also looks at the reference: its missing cell keeps with.missing on."""
    result = tx.seqdist(
        _inputs()[input_kind],
        method="OM",
        subm=sm,
        refseq=_external_reference_with_missing(),
        with_missing=True,
    )

    np.testing.assert_allclose(
        np.asarray(result, dtype=float), R_REFERENCE_DISTANCES[sm], **TOLERANCE
    )


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("input_kind", ["plain array", "sequence object"])
@pytest.mark.parametrize("sm", ["CONSTANT", "TRATE", "INDELSLOG"])
def test_a_reference_whose_missing_code_is_a_source_state_is_refused(sm, input_kind):
    """A code must not mean a real state in the source and ``nr`` in the reference.

    The reference's alphabet is ``{0, 1}``, so its missing cell is code 2, which
    the source uses as an ordinary state.  R refuses the pair ("'refseq' and
    'seqdata' must have the same alphabet").  A plain array declares no alphabet,
    so traminex used to skip the check for it and silently scored every source
    ``2`` as missing.
    """
    r_message = (GOLDENS / "seqdist_with_missing_reset_alphabet_refusal.txt").read_text().strip()
    assert "same alphabet" in r_message
    codes = np.array([[0, 1, 2, 0], [1, 2, 0, 0], [0, 0, 1, 2]])
    source = codes if input_kind == "plain array" else tx.seqdef(codes, alphabet=[0, 1, 2])
    reference = tx.seqdef(np.array([[0, 1, None, 1]], dtype=object), alphabet=[0, 1])
    assert reference.nr_code == 2, "the collision this test is about"

    with pytest.raises(ValueError, match="external refseq must use the same alphabet"):
        tx.seqdist(source, method="OM", subm=sm, refseq=reference, with_missing=True)


@pytest.mark.parametrize("input_kind", ["plain array", "sequence object"])
@pytest.mark.parametrize("sm", sorted(R_REFERENCE_DISTANCES))
def test_a_missing_bearing_reference_requires_with_missing(sm, input_kind):
    """R stops with "'with.missing' must be TRUE when 'seqdata' or 'refseq' contain missing values"."""
    r_message = (GOLDENS / "seqdist_with_missing_reset_refusal.txt").read_text().strip()
    assert "must be TRUE" in r_message and "refseq" in r_message

    with pytest.raises(ValueError, match="with_missing must be true when sequences contain"):
        tx.seqdist(
            _inputs()[input_kind],
            method="OM",
            subm=sm,
            refseq=_external_reference_with_missing(),
            with_missing=False,
        )
