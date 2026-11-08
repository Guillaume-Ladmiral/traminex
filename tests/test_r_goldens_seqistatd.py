"""R-golden parity for seqistatd (T8) vs TraMineR 2.2-14.

Per-sequence state-duration counts (n x nbstat). R's columns are statl (real
alphabet, + nr when with.missing); traminex carries a trailing void column that
R never produces, so the test selects the matching state columns and excludes
void. Goldens: scripts/r_goldens/make_goldens_seqistatd.R.
"""

from __future__ import annotations

import csv

import numpy as np
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden, load_golden_lines

from traminex.stats import seqistatd


def _state_cols(x, include_missing: bool):
    return [
        i
        for i in range(len(x.alphabet))
        if i != x.void_code and (include_missing or i != x.nr_code)
    ]


def _check(x, tag: str, include_missing: bool):
    out = np.asarray(seqistatd(x, include_missing=include_missing), dtype=float)
    cols = _state_cols(x, include_missing)
    np.testing.assert_allclose(
        out[:, cols], load_golden(f"seqistatd_{tag}"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqistatd_a():
    _check(fixture_a(), "A", False)


def test_seqistatd_b_with_and_without_missing():
    _check(fixture_b(), "B_wm0", False)
    _check(fixture_b(), "B_wm1", True)


def test_seqistatd_c_ragged_with_and_without_missing():
    _check(fixture_c(), "C_wm0", False)
    _check(fixture_c(), "C_wm1", True)


def _with_phantom_missing(x):
    x.alphabet.append(x.nr)
    x.nr_code = len(x.alphabet) - 1
    if x.labels is not None:
        x.labels.append(x.nr)
    return x


def _all_void():
    x = fixture_a()
    x.alphabet.append(x.void)
    x.void_code = len(x.alphabet) - 1
    if x.labels is not None:
        x.labels.append(x.void)
    x.data = np.full((2, 3), x.void_code, dtype=int)
    return x


def _r2_cases():
    return {
        "dense_wm0": (fixture_a(), False),
        "dense_wm1": (_with_phantom_missing(fixture_a()), True),
        "missing_wm0": (fixture_b(), False),
        "missing_wm1": (fixture_b(), True),
        "ragged_wm0": (fixture_c(), False),
        "ragged_wm1": (fixture_c(), True),
        "all_void_wm0": (_all_void(), False),
    }


def _r2_rows():
    with (R_GOLDENS_DIR / "seqistatd_r2_prop.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_seqistatd_r2_counts_and_proportions_match_r():
    rows = _r2_rows()
    assert {row["case"] for row in rows} == set(_r2_cases())

    for case, (x, include_missing) in _r2_cases().items():
        expected = [row for row in rows if row["case"] == case]
        counts = seqistatd(x, include_missing=include_missing, prop=False)
        proportions = seqistatd(x, include_missing=include_missing, prop=True)
        for row in expected:
            i = int(row["row"]) - 1
            j = x.alphabet.index(row["state"])
            np.testing.assert_allclose(
                counts[i, j], float(row["count"]), rtol=0, atol=R_GOLDEN_ATOL
            )
            np.testing.assert_allclose(
                proportions[i, j],
                float(row["proportion"]),
                rtol=0,
                atol=R_GOLDEN_ATOL,
                equal_nan=True,
            )


def test_seqistatd_is_quiet_despite_the_recorded_r_message(capsys):
    assert load_golden_lines("seqistatd_r2_messages")[0] == (
        "progress= [>] computing state distribution for 4 sequences ..."
    )
    seqistatd(fixture_a())
    assert capsys.readouterr() == ("", "")


def test_seqistatd_zero_width_keeps_safe_python_extension():
    assert load_golden_lines("seqistatd_r2_messages")[1:] == [
        "zero_width_count=ERROR:subscript out of bounds",
        "zero_width_prop=ERROR:subscript out of bounds",
    ]
    x = fixture_a()
    x.data = x.data[:, :0]
    counts = seqistatd(x, prop=False)
    proportions = seqistatd(x, prop=True)
    np.testing.assert_array_equal(counts, np.zeros((x.nseq, len(x.alphabet))))
    assert np.isnan(proportions).all()
