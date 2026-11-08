"""R-golden tests for the faithful seqerules -- association rules from frequent
event subsequences with the full TraMineR statistics (Support / Conf / Lift /
Standardlift / JMeasure / ImplicStat / p.value / p.valueB1 / p.valueB2), including
R's exact NaN / +-Inf results for degenerate rules.

Generator: scripts/r_goldens/make_goldens_seqerules.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden, load_golden_lines

from traminex import seqecreate, seqefsub, seqerules

_NUMCOLS = [
    "Support",
    "Conf",
    "Lift",
    "Standardlift",
    "JMeasure",
    "ImplicStat",
    "p.value",
    "p.valueB1",
    "p.valueB2",
]


def _actcal_eseq():
    df = pd.read_csv(R_GOLDENS_DIR / "actcal_tse.csv")
    levels = load_golden_lines("eseq_actcal_levels")
    df["event"] = pd.Categorical(df["event"].astype(str), categories=levels)
    return seqecreate(df, time_col="time")


def _syn_eseq():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"],
            "timestamp": [0, 1, 1, 3, 0, 2, 4, 0, 1, 2],
            "event": ["A", "B", "C", "A", "A", "B", "C", "A", "B", "C"],
        }
    )
    return seqecreate(tse)


def _assert_rules(rules: pd.DataFrame, tag: str):
    exp_rules = load_golden_lines("seqerules_" + tag + "_rules")
    golden = load_golden("seqerules_" + tag)
    assert list(rules["Rules"]) == exp_rules
    actual = rules[_NUMCOLS].to_numpy(dtype=float)
    np.testing.assert_allclose(actual, golden, rtol=0, atol=R_GOLDEN_ATOL, equal_nan=True)


def test_seqerules_syn():
    rules = seqerules(seqefsub(_syn_eseq(), min_support=1, max_k=3))
    _assert_rules(rules, "syn")


def test_seqerules_actcal():
    rules = seqerules(seqefsub(_actcal_eseq(), pmin_support=0.02, max_k=3))
    _assert_rules(rules, "actcal")


def test_seqerules_columns_and_type():
    rules = seqerules(seqefsub(_syn_eseq(), min_support=1, max_k=3))
    assert isinstance(rules, pd.DataFrame)
    assert list(rules.columns) == ["Rules"] + _NUMCOLS


def test_seqerules_requires_subseqelist():
    import pytest

    with pytest.raises(TypeError):
        seqerules(_syn_eseq())  # an eseq, not a subseqelist


def test_emitted_column_names_match_the_r_data_frame():
    """cycle-6 §9 (S7): the emitted columns must be R's, not a near-miss.

    R's ``seqerules`` data frame names its columns
    ``Rules, Support, Conf, Lift, Standardlift, JMeasure, ImplicStat, p.value,
    p.valueB1, p.valueB2`` -- as the golden CSVs' own header comments record.
    traminex emitted ``p-value`` for the seventh statistic while keeping the dots
    in its two siblings ``p.valueB1``/``p.valueB2``, and while the neighbouring
    ``seqecmpgroup`` in the same module emits ``p.value``. Porting
    ``rules$p.value`` therefore raised a KeyError on one column of nine, and the
    parity ledger recorded the row as EXACT with "the emitted columns match
    TraMineR".

    The expected names are read from the golden file's recorded header rather
    than restated here, so they stay anchored to what R actually produced.
    """
    header = (R_GOLDENS_DIR / "seqerules_syn.csv").read_text(encoding="utf-8").splitlines()[0]
    assert header.startswith("# columns: ")
    r_numeric_columns = header.removeprefix("# columns: ").split(",")

    rules = seqerules(seqefsub(_syn_eseq(), min_support=1, max_k=3))
    assert list(rules.columns) == ["Rules"] + r_numeric_columns
