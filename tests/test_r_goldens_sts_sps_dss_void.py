"""to_SPS and to_DSS leave void cells out, as R's SPS and DSS forms do (TX-AUD-20261010-033).

R's seqformat(to="SPS"/"DSS") turns void cells into NA, so a shorter sequence
simply ends, while a missing cell stays a state. traminex's string helpers
printed the void marker as a state ("A/1-B/1-%/2", "A-B-%"). Fixture:
``scripts/r_goldens/make_goldens_sts_sps_dss_void.R``.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

import traminex as tx
from traminex.seqformat import to_DSS, to_SPS

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "sts_sps_dss_void.csv", keep_default_na=False
)
RAGGED = tx.seqdef(
    np.array([["A", "B", None, None], ["A", None, "B", "B"], ["B", "B", None, None]], dtype=object),
    alphabet=["A", "B"],
)


def _r_strings(form: str, *, sps: bool) -> list[str]:
    out = []
    rows = GOLDEN[GOLDEN["form"] == form]
    for _, tokens in rows.groupby("sequence", sort=True):
        kept = [token for token in tokens.sort_values("position")["token"] if token != "<NA>"]
        if sps:
            kept = [re.sub(r"^\((.*),(\d+)\)$", r"\1/\2", token) for token in kept]
        out.append("-".join(kept))
    return out


def test_to_sps_leaves_void_out_as_r_sps_does() -> None:
    assert to_SPS(RAGGED) == _r_strings("SPS", sps=True) == ["A/1-B/1", "A/1-*/1-B/2", "B/2"]


def test_to_dss_leaves_void_out_as_r_dss_does() -> None:
    assert to_DSS(RAGGED) == _r_strings("DSS", sps=False) == ["A-B", "A-*-B", "B"]
