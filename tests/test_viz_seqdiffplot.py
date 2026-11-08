"""RP7a: seqdiffplot (plot.seqdiff) structure test.

A single-stat line == the per-position statistic (default "Pseudo R2") and the
multi-line discrepancy mode == one line per group + Total. Asserted against the
seqdiff_uw_stat / seqdiff_uw_disc goldens (rows = positions).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.format import seqdef
from traminex.visualization.plots.discrepancy import seqdiffplot

_ALPH = ["1", "2", "3", "4"]


def _fixture():
    codes = load_golden("seqdiff_input_codes").astype(int)
    w = load_golden("seqdiff_input_weights").ravel()
    chars = np.array([[_ALPH[c] for c in row] for row in codes], dtype=object)
    return seqdef(chars, alphabet=_ALPH, weights=w)


def test_seqdiffplot_pseudo_r2_line_matches_golden():
    x = _fixture()
    g = load_golden("seqdiff_input_group").ravel().astype(int)
    fig, ax = seqdiffplot(x, g, method="LCS", norm=None, weighted=False, stat="Pseudo R2")
    (line,) = ax.lines
    # stat column order is [PseudoF, PseudoFbf, PseudoR2, Bartlett, Levene] -> R2 is col 2
    expected = load_golden("seqdiff_uw_stat")[:, 2]
    np.testing.assert_allclose(line.get_ydata(), expected, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqdiffplot_discrepancy_lines_match_golden():
    x = _fixture()
    g = load_golden("seqdiff_input_group").ravel().astype(int)
    fig, ax = seqdiffplot(x, g, method="LCS", norm=None, weighted=False, stat="discrepancy")
    disc = load_golden("seqdiff_uw_disc")  # columns = groups + Total
    assert len(ax.lines) == disc.shape[1]
    got = np.column_stack([line.get_ydata() for line in ax.lines])
    np.testing.assert_allclose(got, disc, rtol=0, atol=R_GOLDEN_ATOL)


# R's plot.seqdiff draws the object seqdiff returned, and seqdiff's default is
# seqdist.args = list(method = "LCS", norm = "auto"); seqdiffplot must plot the
# same series by default (TX-AUD-20261010-007). R values: the public cases of
# scripts/r_goldens/make_goldens_seqdiff.R.
def _public_case(case: str, kind: str) -> dict[str, np.ndarray]:
    import csv
    from collections import defaultdict

    from r_goldens_utils import R_GOLDENS_DIR

    columns: dict[str, dict[int, float]] = defaultdict(dict)
    with (R_GOLDENS_DIR / "seqdiff_public_cases.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["case"] == case and row["kind"] == kind:
                columns[row["column"]][int(row["row_position"])] = float(row["value"])
    return {name: np.array([values[i] for i in sorted(values)]) for name, values in columns.items()}


def test_seqdiffplot_default_series_match_r_seqdiff_default():
    x = _fixture()
    g = load_golden("seqdiff_input_group").ravel().astype(int)
    fig, ax = seqdiffplot(x, g)
    (line,) = ax.lines
    np.testing.assert_allclose(
        line.get_ydata(),
        _public_case("default_lcs_auto", "stat")["Pseudo R2"],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    fig, ax = seqdiffplot(x, g, stat="discrepancy")
    expected = _public_case("default_lcs_auto", "discrepancy")
    assert [line.get_label() for line in ax.lines] == list(expected)
    for line in ax.lines:
        np.testing.assert_allclose(
            line.get_ydata(), expected[line.get_label()], rtol=0, atol=R_GOLDEN_ATOL
        )


def test_seqdiffplot_accepts_a_seqdist_args_map():
    x = _fixture()
    g = load_golden("seqdiff_input_group").ravel().astype(int)
    fig, ax = seqdiffplot(
        x, g, seqdist_args={"method": "LCS", "norm": "none"}, weighted=False, stat="discrepancy"
    )
    expected = _public_case("lcs_none_map", "discrepancy")
    assert [line.get_label() for line in ax.lines] == list(expected)
    for line in ax.lines:
        np.testing.assert_allclose(
            line.get_ydata(), expected[line.get_label()], rtol=0, atol=R_GOLDEN_ATOL
        )
