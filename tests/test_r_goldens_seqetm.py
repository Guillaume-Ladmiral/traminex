"""R-golden tests for ``seqetm`` -- the state-alphabet transition/event matrix.

``TraMineR::seqetm`` builds an ``nbstat x nbstat`` character matrix from a
state-sequence object's alphabet (and optional long labels). We pin the three
encoding methods (``transition`` / ``period`` / ``state``), with and without
labels, against live TraMineR 2.2-14.

Generator: scripts/r_goldens/make_goldens_seqetm.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

from traminex import STSList, seqetm


def _load_tm(name: str) -> pd.DataFrame:
    path = R_GOLDENS_DIR / f"{name}.tsv"
    lines = [
        ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")
    ]
    cols = lines[0].split("\t")
    idx: list[str] = []
    rows: list[list[str]] = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        idx.append(parts[0])
        rows.append(parts[1:])
    return pd.DataFrame(rows, index=idx, columns=cols)


def _make_sd() -> STSList:
    # Mirrors the seqdef() in scripts/r_goldens/make_goldens_seqetm.R.
    mat = [
        ["A", "A", "B", "C"],
        ["B", "B", "B", "A"],
        ["C", "A", "A", "B"],
        ["A", "B", "C", "C"],
    ]
    alphabet = ["A", "B", "C"]
    code = {s: i for i, s in enumerate(alphabet)}
    data = np.array([[code[s] for s in row] for row in mat], dtype=int)
    return STSList(data=data, alphabet=alphabet, labels=["Apple", "Banana", "Cherry"])


def test_seqetm_transition_matches_r():
    pd.testing.assert_frame_equal(
        seqetm(_make_sd(), method="transition"),
        _load_tm("seqetm_transition"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqetm_period_matches_r():
    pd.testing.assert_frame_equal(
        seqetm(_make_sd(), method="period"),
        _load_tm("seqetm_period"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqetm_state_matches_r():
    pd.testing.assert_frame_equal(
        seqetm(_make_sd(), method="state"),
        _load_tm("seqetm_state"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqetm_transition_no_labels_matches_r():
    got = seqetm(_make_sd(), method="transition", use_labels=False)
    pd.testing.assert_frame_equal(
        got, _load_tm("seqetm_transition_nolab"), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


def _make_sdm() -> STSList:
    # Mirrors seqdef(matm, alphabet=c("A","B","C"), missing="*") in the generator:
    # interior '*' stays nr, the trailing '*' of row 2 becomes void '%'. R's seqetm
    # appends nr and void ONCE each -> alphabet A,B,C,*,% (N+1).
    alphabet = ["A", "B", "C", "*", "%"]  # codes 0..4; nr='*'(3), void='%'(4)
    data = np.array(
        [
            [0, 0, 3, 1],  # A A * B   (interior nr)
            [1, 1, 0, 4],  # B B A %   (trailing void)
            [2, 0, 0, 1],  # C A A B
        ],
        dtype=int,
    )
    return STSList(
        data=data, alphabet=alphabet, labels=None, nr="*", void="%", nr_code=3, void_code=4
    )


def test_seqetm_missing_void_no_double_append_matches_r():
    # Regression for the nr/void double-append bug: traminex's x.alphabet already
    # includes nr/void, so seqetm must build statl from the REAL states only.
    got = seqetm(_make_sdm(), method="transition", use_labels=False)
    pd.testing.assert_frame_equal(
        got, _load_tm("seqetm_transition_miss"), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqetm_missing_void_state_method_matches_r():
    got = seqetm(_make_sdm(), method="state", use_labels=False)
    pd.testing.assert_frame_equal(
        got, _load_tm("seqetm_state_miss"), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )
