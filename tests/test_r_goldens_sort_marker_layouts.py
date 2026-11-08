"""Sorting ranks missing and void after the states in any codebook layout.

TX-AUD-20261011-006 and TX-AUD-20261011-014. R orders cells over the factor
levels ``c(alphabet, nr, void)``: plot.stslist and seqrf sort with
``do.call(order, ...)`` over the state columns, and TraMineRextras' sorti/sortv
first replace void by nr. traminex ranked cells by their raw integer codes, so
an STSList whose markers are stored inside or before the real states sorted
differently from the same cells stored as seqdef stores them (missing last).
Fixtures: ``scripts/r_goldens/make_goldens_sort_marker_layouts.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx
from traminex import STSList
from traminex.visualization.series import resolve_sort_order

GOLDENS = Path(__file__).resolve().parent / "r_goldens"


def _read(name: str, *keys: str) -> dict[tuple[str, ...], list[str]]:
    out: dict[tuple[str, ...], list[str]] = defaultdict(list)
    with (GOLDENS / name).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            out[tuple(row[key] for key in keys)].append(row["value"])
    return out


R = {case: values for (case,), values in _read("sort_marker_layouts.csv", "case").items()}
R_RF = _read("sort_marker_layouts_seqrf.csv", "case", "kind")
CELLS = [R["cells"][3 * row : 3 * row + 3] for row in range(6)]  # A/B/C real, '*' missing, '%' void


def _layout(alphabet: list[str], nr_code: int, void_code: int) -> STSList:
    position = {state: alphabet.index(state) for state in ("A", "B", "C")}
    position.update({"*": nr_code, "%": void_code})
    data = np.array([[position[c] for c in row] for row in CELLS], dtype=np.uint8)
    return STSList(data, alphabet, nr_code=nr_code, void_code=void_code)


LAYOUTS = {
    "seqdef": tx.seqdef(
        [[None if c in ("*", "%") else c for c in row] for row in CELLS], alphabet=["A", "B", "C"]
    ),
    "trailing": _layout(["A", "B", "C", "*", "%"], 3, 4),
    "interior": _layout(["A", "*", "B", "%", "C"], 1, 3),
    "leading": _layout(["*", "%", "A", "B", "C"], 0, 1),
    "external": _layout(["A", "B", "C"], 7, 8),
}
ORDERS = {
    "sorti_beg": lambda x: tx.sorti(x, start="beg"),
    "sorti_end": lambda x: tx.sorti(x, start="end"),
    "sortv_beg": lambda x: tx.sortv(x, start="beg"),
    "sortv_end": lambda x: tx.sortv(x, start="end"),
    "plot_from_start": lambda x: resolve_sort_order(x, "from.start"),
    "plot_from_end": lambda x: resolve_sort_order(x, "from.end"),
}


@pytest.mark.parametrize("layout", sorted(LAYOUTS))
def test_layouts_hold_the_same_cells(layout: str) -> None:
    x = LAYOUTS[layout]
    decoded = [[x._state_for_code(int(code)) for code in row] for row in np.asarray(x.data)]
    assert decoded == CELLS


@pytest.mark.parametrize("order", sorted(ORDERS))
@pytest.mark.parametrize("layout", sorted(LAYOUTS))
def test_sort_orders_match_r_in_every_layout(layout: str, order: str) -> None:
    got = (np.asarray(ORDERS[order](LAYOUTS[layout])) + 1).tolist()
    assert got == [int(v) for v in R[order]]


def _mvad_layouts() -> dict[str, STSList]:
    cells = pd.read_csv(
        GOLDENS / "sort_marker_layouts_seqrf_cells.csv", dtype=str, keep_default_na=False
    )
    rows = [[None if c == "" else c for c in row] for row in cells.to_numpy(dtype=object).tolist()]
    alphabet = ["employment", "FE", "HE", "joblessness", "school", "training"]
    seqdef_layout = tx.seqdef(rows, alphabet=alphabet, right=None)
    nr = seqdef_layout.nr_code
    order = [nr] + [code for code in range(len(seqdef_layout.alphabet)) if code != nr]
    remap = {old: new for new, old in enumerate(order)}
    data = np.vectorize(remap.get)(np.asarray(seqdef_layout.data)).astype(int)
    missing_first = STSList(
        data=data,
        alphabet=[seqdef_layout.alphabet[code] for code in order],
        nr=seqdef_layout.nr,
        void=seqdef_layout.void,
        nr_code=remap[nr],
        void_code=None,
    )
    return {"seqdef": seqdef_layout, "missing_first": missing_first}


MVAD = _mvad_layouts()
MVAD_DISS = pd.read_csv(GOLDENS / "sort_marker_layouts_seqrf_diss.csv").to_numpy(float)


@pytest.mark.parametrize("layout", ["seqdef", "missing_first"])
@pytest.mark.parametrize("sortv", ["from.start", "from.end"])
def test_seqrf_order_is_layout_independent(sortv: str, layout: str) -> None:
    rf = tx.seqrf(MVAD[layout], MVAD_DISS, k=10, sortv=sortv, weighted=False)["rf"]
    np.testing.assert_allclose(
        float(rf["R2"]), float(R_RF[(sortv, "R2")][0]), rtol=0, atol=R_GOLDEN_ATOL
    )
    assert (np.asarray(rf["medoids"]) + 1).tolist() == [int(v) for v in R_RF[(sortv, "medoid")]]
