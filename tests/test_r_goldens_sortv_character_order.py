"""A text sortv is ordered as R's order() orders it (TX-AUD-20261010-060).

plot.stslist sorts rows with order(sortv): text in R's collation (ICU), NA
last. traminex's index-plot resolver used NumPy's code-point argsort, so
['b', 'A', 'a', 'B', '_x', 'c'] came out in a different order, and a None
raised TypeError where R puts NA last. A numeric or logical sortv holding
None was then ordered as the inverse of R's order (batch-2 review I-1).
Fixture: ``scripts/r_goldens/make_goldens_sortv_character_order.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.visualization.series import resolve_sort_order

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "sortv_character_order.csv"


def _golden() -> dict[str, dict[str, list[str]]]:
    out: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            out[row["case"]][row["kind"]].append(row["value"])
    return out


R = _golden()


def _seqs(n: int) -> tx.STSList:
    return tx.seqdef(np.array([["A"]] * n, dtype=object), alphabet=["A"])


def _sortv(case: str) -> list[object]:
    values: list[object] = []
    for text in R[case]["value"]:
        if text == "<NA>":
            values.append(None)
        elif case.startswith("numeric"):
            values.append(float(text))
        elif case.startswith("logical"):
            values.append({"TRUE": True, "FALSE": False}[text])
        else:
            values.append(text)
    return values


@pytest.mark.parametrize("case", ["text_case_punctuation", "text_missing_accent"])
def test_text_sortv_order_matches_r(case: str) -> None:
    sortv = _sortv(case)
    order = resolve_sort_order(_seqs(len(sortv)), sortv)
    assert (np.asarray(order) + 1).tolist() == [int(v) for v in R[case]["order"]]


NUMERIC_CASES = [
    "numeric_missing",
    "numeric_missing_five",
    "numeric_infinite_missing",
    "logical_missing",
]


@pytest.mark.parametrize("case", NUMERIC_CASES)
@pytest.mark.parametrize("missing", [None, float("nan")], ids=["None", "NaN"])
def test_numeric_sortv_with_missing_matches_r(case: str, missing: object) -> None:
    sortv = [missing if value is None else value for value in _sortv(case)]
    order = resolve_sort_order(_seqs(len(sortv)), sortv)
    assert (np.asarray(order) + 1).tolist() == [int(v) for v in R[case]["order"]]


def test_mixed_text_and_numbers_is_refused() -> None:
    with pytest.raises(TypeError, match="mix"):
        resolve_sort_order(_seqs(3), ["a", 1, "b"])


def _factor_spellings():
    import pandas as pd

    values = [None if v == "<NA>" else v for v in R["factor_levels_z_a"]["value"]]
    levels = R["factor_levels_z_a"]["levels"]
    return {
        "Categorical": pd.Categorical(values, categories=levels),
        "Series": pd.Series(values, dtype=pd.CategoricalDtype(levels)),
    }


@pytest.mark.parametrize("spelling", ["Categorical", "Series"])
def test_categorical_sortv_follows_its_levels_as_r_does(spelling: str) -> None:
    # R orders a factor by its level codes, NA last (TX-AUD-20261019-019): a category Series
    # was ordered as text.
    sortv = _factor_spellings()[spelling]
    order = resolve_sort_order(_seqs(len(sortv)), sortv)
    assert (np.asarray(order) + 1).tolist() == [int(v) for v in R["factor_levels_z_a"]["order"]]


@pytest.mark.parametrize("spelling", ["Categorical", "Series"])
@pytest.mark.parametrize("grouped", [False, True])
def test_index_plot_rows_follow_the_levels_on_every_route(spelling: str, grouped: bool) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sortv = _factor_spellings()[spelling]
    ids = [f"r{i}" for i in range(1, len(sortv) + 1)]
    seqs = tx.seqdef(np.array([["A"]] * len(sortv), dtype=object), alphabet=["A"], id=ids)
    kwargs = {"group": ["g"] * len(sortv)} if grouped else {}
    fig, axes = tx.seqplot(seqs, type="i", sortv=sortv, **kwargs)
    ax = axes[0] if grouped else axes
    labels = [text.get_text() for text in ax.get_yticklabels()]
    plt.close(fig)
    assert labels == [ids[int(v) - 1] for v in R["factor_levels_z_a"]["order"]]
