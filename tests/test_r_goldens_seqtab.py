"""Public TraMineR 2.2-14 evidence for ``seqtab`` and ``seqtabstocc``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_a, fixture_b
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import STSList, seqconc
from traminex.format import seqdef
from traminex.stats import seqtab, seqtabstocc

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"

_S = np.array(
    [
        ["A", "A", "B", "C"],
        ["A", "A", "B", "C"],
        ["A", "B", "B", "C"],
        ["A", "A", "B", "C"],
        ["A", "B", "B", "C"],
    ],
    dtype=object,
)
_WS = np.array([1.0, 2.0, 1.0, 3.0, 1.0])


def _fixture_s(weighted: bool = False):
    return seqdef(_S, alphabet=["A", "B", "C"], weights=_WS if weighted else None)


_T = np.array(
    [
        ["A", "A", "A", "A"],
        ["A", "A", "A", "B"],
        ["A", "A", "B", "B"],
        ["A", "B", "B", "B"],
        ["B", "B", "B", "B"],
        ["B", "B", "B", "C"],
        ["B", "B", "C", "C"],
        ["B", "C", "C", "C"],
        ["C", "C", "C", "C"],
        ["C", "C", "C", "A"],
        ["A", "*", "B", "%"],
        ["C", "A", "B", "%"],
        ["A", "A", "A", "A"],
        ["A", "A", "A", "B"],
    ],
    dtype=object,
)
_WT = np.array([2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 2], dtype=float)


def _fixture_t(weighted: bool = False) -> STSList:
    x = seqdef(
        _T,
        alphabet=["A", "B", "C"],
        nr="*",
        weights=_WT if weighted else None,
        id=[f"case-{index}" for index in range(1, len(_T) + 1)],
        cpal=["#111111", "#222222", "#333333"],
        cnames=["t1", "t2", "t3", "t4"],
        start=5,
    )
    x.labels = ["Alpha", "Beta", "Gamma", "Missing", "Void"]
    return x


def _fixture_all_nr() -> STSList:
    return seqdef(
        np.array(
            [["*", "*", "*"], ["A", "A", "A"], ["A", "A", "A"]],
            dtype=object,
        ),
        alphabet=["A", "B"],
        nr="*",
        id=["empty", "first", "second"],
    )


def _fixture_locale_ties() -> STSList:
    states = ["a", "B", "A", "b", "á", "ä", "Z"]
    return seqdef(
        np.array([[state, state] for state in states], dtype=object),
        alphabet=states,
    )


PUBLIC_CASES = {
    "default": (_fixture_t, {}),
    "all": (_fixture_t, {"idxs": 0}),
    "scalar_rank": (_fixture_t, {"idxs": 2}),
    "reversed_ranks": (_fixture_t, {"idxs": [2, 1]}),
    "duplicate_ranks": (_fixture_t, {"idxs": [2, 1, 2]}),
    "out_of_range": (_fixture_t, {"idxs": 20}),
    "weighted_all": (lambda: _fixture_t(weighted=True), {"idxs": 0, "weighted": True}),
    "weighted_ignored": (
        lambda: _fixture_t(weighted=True),
        {"idxs": 0, "weighted": False},
    ),
    "sts_all": (_fixture_t, {"idxs": 0, "format": "STS"}),
    "locale_ties": (
        _fixture_locale_ties,
        {"idxs": 0, "weighted": False, "format": "STS"},
    ),
    "all_nr_removed": (_fixture_all_nr, {"idxs": 0}),
    "tlim_only": (_fixture_t, {"tlim": [2, 1]}),
    "tlim_with_default": (
        _fixture_t,
        {"idxs": list(range(1, 11)), "tlim": [2, 1]},
    ),
}


def _seqtab_numbers(df):
    """count, prop*100 sorted descending by count (to align with R's sorted Freq)."""
    arr = np.column_stack([df["count"].to_numpy(), df["prop"].to_numpy() * 100.0])
    order = np.argsort(-arr[:, 0], kind="stable")
    return arr[order]


def _public_warning(case: str) -> str:
    records = pd.read_csv(
        GOLDEN_DIR / "seqtab_public_warnings.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    return records.loc[records["case"] == case, "message"].item()


def test_seqtab_unweighted():
    got = _seqtab_numbers(seqtab(_fixture_s()))
    np.testing.assert_allclose(got, load_golden("seqtab_S"), rtol=0, atol=R_GOLDEN_ATOL)


def test_seqtab_weighted():
    got = _seqtab_numbers(seqtab(_fixture_s(weighted=True)))
    np.testing.assert_allclose(got, load_golden("seqtab_Sw"), rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("case", list(PUBLIC_CASES))
def test_seqtab_complete_public_result_matches_r(case: str):
    builder, kwargs = PUBLIC_CASES[case]
    expected = pd.read_csv(
        GOLDEN_DIR / "seqtab_public_cases.csv",
        dtype={
            "case": str,
            "call": str,
            "sequence": str,
            "weighted": str,
            "representative_sts": str,
        },
        keep_default_na=False,
    )
    expected = expected.loc[expected["case"] == case].reset_index(drop=True)
    assert not expected.empty

    if case.startswith("tlim"):
        assert "tlim is deprecated" in _public_warning(case)
        with pytest.warns(DeprecationWarning, match="tlim.*deprecated"):
            actual = seqtab(builder(), **kwargs)
    else:
        actual = seqtab(builder(), **kwargs)

    assert list(actual.columns) == ["sequence", "count", "prop", "cumprop"]
    assert actual["sequence"].tolist() == expected["sequence"].tolist()
    np.testing.assert_allclose(actual["count"], expected["Freq"], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(
        actual["prop"],
        expected["Percent"] / 100.0,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        actual["cumprop"],
        np.cumsum(expected["Percent"].to_numpy(dtype=float) / 100.0),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )

    freq = actual.attrs["freq"]
    assert isinstance(freq, pd.DataFrame)
    assert list(freq.columns) == ["Freq", "Percent"]
    assert freq.index.tolist() == actual["sequence"].tolist()
    np.testing.assert_allclose(freq["Freq"], expected["Freq"], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(freq["Percent"], expected["Percent"], rtol=0, atol=R_GOLDEN_ATOL)
    assert actual.attrs["nbseq"] == pytest.approx(
        float(expected["nbseq"].iloc[0]), rel=0, abs=R_GOLDEN_ATOL
    )
    assert actual.attrs["weighted"] is (expected["weighted"].iloc[0] == "TRUE")
    assert actual.attrs["idxs"] == expected["idxs"].astype(int).tolist()
    # R's idxf is one-based (TX-AUD-20261019-034); source_indices are zero-based rows.
    assert actual.attrs["idxf"] == expected["idxf"].astype(int).tolist()
    assert actual.attrs["source_indices"] == (expected["source_index"].astype(int) - 1).tolist()
    assert actual.attrs["format"] == expected["format"].iloc[0]

    representatives = actual.attrs["representatives"]
    assert isinstance(representatives, STSList)
    assert representatives is not builder()
    assert representatives.alphabet == builder().alphabet
    assert representatives.labels == builder().labels
    assert representatives.cpal == builder().cpal
    assert representatives.xtlab == builder().xtlab
    assert representatives.start == builder().start
    assert seqconc(representatives).tolist() == expected["representative_sts"].tolist()
    np.testing.assert_allclose(
        representatives.weights,
        expected["representative_weight"],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    assert not np.shares_memory(representatives.data, builder().data)


def test_seqtab_rejects_invalid_public_rank_and_format_inputs():
    x = _fixture_t()
    errors = pd.read_csv(
        GOLDEN_DIR / "seqtab_public_errors.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    assert set(errors["case"]) == {
        "idxs_null",
        "idxs_negative",
        "idxs_mixed_zero",
        "idxs_tlim_conflict",
        "invalid_format",
    }

    # [0] is not here: R's c(0) is the scalar 0, i.e. every rank (seqtab.R:12,
    # 50; live R in tests/r_goldens/audit20260922_s3c.txt; audit S3-21).
    for value in (None, -1, [-1], [0, 1]):
        with pytest.raises((TypeError, ValueError), match="idxs"):
            seqtab(x, idxs=value)
    for value in (True, 1.5, [True], [1.5], [1, 2.5], [], [[1, 2]], "1"):
        with pytest.raises((TypeError, ValueError), match="idxs"):
            seqtab(x, idxs=value)
    with pytest.raises(ValueError, match="format"):
        seqtab(x, format="sps")
    with pytest.raises(ValueError, match="cannot be specified together"):
        seqtab(x, idxs=[2, 1], tlim=1)


def test_seqtab_python_extensions_remain_explicit_and_nonmutating():
    x = _fixture_t(weighted=True)
    original = x.copy()
    explicit = np.arange(1, x.nseq + 1, dtype=float)
    explicit_before = explicit.copy()

    sts = seqtab(x, idxs=0, format="STS", sep="|", weights=explicit)
    assert "|" in sts["sequence"].iloc[0]
    assert sts.attrs["nbseq"] == pytest.approx(explicit.sum(), rel=0, abs=R_GOLDEN_ATOL)

    ignored = seqtab(x, idxs=0, weighted=False, weights=explicit)
    assert ignored.attrs["nbseq"] == pytest.approx(x.nseq, rel=0, abs=R_GOLDEN_ATOL)
    assert ignored.attrs["weighted"] is False

    unsorted = seqtab(x, idxs=0, sort=False)
    assert unsorted["sequence"].tolist() == sorted(unsorted["sequence"].tolist())
    unnormalized = seqtab(x, idxs=0, normalize=False)
    np.testing.assert_array_equal(unnormalized[["prop", "cumprop"]], 0.0)

    without_missing = seqtab(x, idxs=0, include_missing=False)
    assert all("*" not in label for label in without_missing["sequence"])
    assert all("%" not in label for label in without_missing["sequence"])

    np.testing.assert_array_equal(x.data, original.data)
    np.testing.assert_array_equal(x.weights, original.weights)
    np.testing.assert_array_equal(explicit, explicit_before)
    assert x.ids == original.ids


def test_seqtab_validates_explicit_weights_before_filtering_and_warns_for_dash_state():
    x = _fixture_all_nr()
    with pytest.raises(ValueError, match="weights.*length"):
        seqtab(x, idxs=0, weights=np.ones(2))

    dash = seqdef(
        np.array([["A-B", "C"], ["C", "A-B"]], dtype=object),
        alphabet=["A-B", "C"],
    )
    with pytest.warns(UserWarning, match="'-'|character"):
        seqtab(dash, idxs=0)


def _stocc_counts(df):
    """Multiset of co-occurrence counts (order-independent: labels differ R vs Py)."""
    return np.sort(df["count"].to_numpy())[::-1]


def test_seqtabstocc_a():
    got = _stocc_counts(seqtabstocc(fixture_a()))
    golden = load_golden("seqtabstocc_A")
    np.testing.assert_allclose(got, np.sort(golden[:, 0])[::-1], rtol=0, atol=R_GOLDEN_ATOL)


def test_seqtabstocc_b():
    got = _stocc_counts(seqtabstocc(fixture_b()))
    golden = load_golden("seqtabstocc_B")
    np.testing.assert_allclose(got, np.sort(golden[:, 0])[::-1], rtol=0, atol=R_GOLDEN_ATOL)
