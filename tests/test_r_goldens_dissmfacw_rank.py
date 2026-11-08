"""dissmfacw projects onto the column space of a rank-deficient design (TX-AUD-20261010-040).

R's dissmfacw builds its hat matrices from all ``ncol(X)`` columns of
``qr.Q(qr(X))`` and counts degrees of freedom in columns (dissmfacw.R:67-69),
so a design whose rank is lower than its column count gains one arbitrary
direction; traminex copied the construction with NumPy's QR, whose arbitrary
direction differs. Both are artefacts. The corrected statistics project onto the
column space and count ranks: the model's values are those of the same design
without its redundant column, where R is right, and a term that adds no rank
given the others has no test (NaN, with a warning). Fixture:
``scripts/r_goldens/make_goldens_dissmfacw_rank.R`` (mvad 1:40, LCS).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_DIFFERENCE_ATOL

import traminex as tx

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(HERE / "dissmfacw_rank.csv", dtype={"PseudoF": str, "PseudoR2": str})
DATA = pd.read_csv(HERE / "dissmfacw_rank_data.csv", dtype={"weight": float})
DISS = np.loadtxt(HERE / "dissmfacw_rank_diss.csv", delimiter=",")
REGIONS = ["Belfast", "N.Eastern", "Southern", "S.Eastern", "Western"]


def _r(case: str) -> pd.DataFrame:
    rows = GOLDEN[GOLDEN["case"] == case].set_index("Variable")
    return rows[["PseudoF", "PseudoR2"]].astype(float)


def _ours(columns: dict) -> pd.DataFrame:
    result = tx.dissmfacw(DISS, columns, nperm=1)
    return result.mfac.set_index("Variable")[["PseudoF", "PseudoR2"]].astype(float)


def _close(ours: pd.Series, r: pd.Series) -> None:
    np.testing.assert_allclose(ours.to_numpy(), r.to_numpy(), rtol=0, atol=1e-10)


def test_rank_deficient_regions_give_the_reduced_model() -> None:
    columns = {"male": DATA["male"].tolist(), **{name: DATA[name].tolist() for name in REGIONS}}
    with pytest.warns(UserWarning, match="rank"):
        ours = _ours(columns)
    reduced = _r("regions_reduced")
    _close(ours.loc["Total"], reduced.loc["Total"])
    _close(ours.loc["male"], reduced.loc["male"])
    assert ours.loc[REGIONS, "PseudoF"].isna().all()
    assert ours.loc[REGIONS, "PseudoR2"].isna().all()


def test_r_rank_deficient_values_are_not_the_reduced_model() -> None:
    full, reduced = _r("regions_full"), _r("regions_reduced")
    assert abs(full.loc["Total", "PseudoF"] - reduced.loc["Total", "PseudoF"]) > 0.1
    assert abs(full.loc["male", "PseudoF"] - reduced.loc["male", "PseudoF"]) > 0.1


def test_a_collinear_numeric_copy_gives_the_single_predictor_model() -> None:
    weight = DATA["weight"].to_numpy()
    with pytest.warns(UserWarning, match="rank"):
        ours = _ours({"weight": weight, "weight2": 2 * weight})
    alone = _r("numeric_alone")
    _close(ours.loc["Total"], alone.loc["Total"])
    assert ours.loc[["weight", "weight2"], "PseudoF"].isna().all()
    assert _r("numeric_collinear").loc["Total", "PseudoF"] != pytest.approx(
        alone.loc["Total", "PseudoF"], rel=0, abs=R_DIFFERENCE_ATOL
    )


def test_a_full_rank_design_keeps_r_values() -> None:
    columns = {"male": DATA["male"].tolist(), **{name: DATA[name].tolist() for name in REGIONS[:4]}}
    ours = _ours(columns)
    reduced = _r("regions_reduced")
    for variable in ["male", *REGIONS[:4], "Total"]:
        _close(ours.loc[variable], reduced.loc[variable])


def test_a_constant_predictor_has_no_test() -> None:
    constant = np.ones(len(DATA))
    with pytest.warns(UserWarning, match="rank 1 for 2 columns"):
        ours = _ours({"constant": constant})
    assert ours["PseudoF"].isna().all()
    assert np.isfinite(_r("constant")["PseudoF"]).all()  # R's artefact of its column count
