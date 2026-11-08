"""R-golden coverage for the complete ``TraMineR::seqlength`` contract.

Generator: ``scripts/r_goldens/make_goldens_seqlength_full.R`` against
TraMineR 2.2-14.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden

from traminex import STSList, seqdef, seqlength

ALPHABET = ["A", "B", "C"]
SEQUENCES_B = [
    ["A", "A", "*", "B", "B", "C"],
    ["*", "A", "A", "A", "*", "*"],
    ["A", "*", "A", "B", "B", "B"],
    ["C", "C", "C", "C", "C", "C"],
]
SEQUENCES_C = [
    ["A", "A", "B", "C", "C", "C"],
    ["A", "B", "B", "%", "%", "%"],
    ["A", "*", "A", "B", "%", "%"],
    ["C", "C", "%", "%", "%", "%"],
]


def _fixture_b() -> STSList:
    return seqdef(SEQUENCES_B, alphabet=ALPHABET)


def _fixture_c() -> STSList:
    return seqdef(SEQUENCES_C, alphabet=ALPHABET)


def _golden(name: str) -> np.ndarray:
    return load_golden(name).astype(int).ravel()


def _golden_row_names(name: str) -> list[str]:
    path = R_GOLDENS_DIR / f"{name}.csv"
    prefix = "# rownames: "
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix):
            return line.removeprefix(prefix).split(",")
    raise AssertionError(f"{path} has no row-name provenance")


def test_seqlength_interior_missing_matches_r_with_both_settings():
    x = _fixture_b()

    np.testing.assert_array_equal(seqlength(x, with_missing=False), _golden("seqlength_B_wm0"))
    np.testing.assert_array_equal(seqlength(x, with_missing=True), _golden("seqlength_B_wm1"))
    np.testing.assert_array_equal(seqlength(x), _golden("seqlength_B_wm1"))
    assert _golden_row_names("seqlength_B_wm0") == ["B1", "B2", "B3", "B4"]


def test_seqlength_ragged_void_and_missing_matches_r():
    x = _fixture_c()

    without_missing = seqlength(x, with_missing=False)
    with_missing = seqlength(x, with_missing=True)
    np.testing.assert_array_equal(without_missing, _golden("seqlength_C_wm0"))
    np.testing.assert_array_equal(with_missing, _golden("seqlength_C_wm1"))
    assert without_missing.shape == (4,)
    assert np.issubdtype(without_missing.dtype, np.integer)
    assert _golden_row_names("seqlength_C_wm1") == ["C1", "C2", "C3", "C4"]


def test_seqlength_all_void_and_empty_width_are_zero():
    all_void = seqdef([["%", "%", "%"], ["%", "%", "%"]], alphabet=[])
    np.testing.assert_array_equal(seqlength(all_void), np.array([0, 0], dtype=int))
    np.testing.assert_array_equal(
        seqlength(all_void, with_missing=False), np.array([0, 0], dtype=int)
    )

    empty_width = STSList(np.empty((3, 0), dtype=int), alphabet=[])
    np.testing.assert_array_equal(seqlength(empty_width), np.zeros(3, dtype=int))
    np.testing.assert_array_equal(
        seqlength(empty_width, with_missing=False), np.zeros(3, dtype=int)
    )


def test_seqlength_dataframe_extension_uses_fallback_then_metadata_ids():
    x = _fixture_c()
    expected_lengths = _golden("seqlength_C_wm1")

    fallback = seqlength(x, as_dataframe=True)
    pd.testing.assert_frame_equal(
        fallback,
        pd.DataFrame({"id": [1, 2, 3, 4], "seqlength": expected_lengths}),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )

    x.ids = ["C1", "C2", "C3", "C4"]  # type: ignore[attr-defined]
    with_ids = seqlength(x, as_dataframe=True)
    pd.testing.assert_frame_equal(
        with_ids,
        pd.DataFrame({"id": x.ids, "seqlength": expected_lengths}),  # type: ignore[attr-defined]
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqlength_rejects_non_sequence_input():
    with pytest.raises(TypeError, match="state sequence object"):
        seqlength([["A", "B"]])
