"""seqMD details against R (TX-AUD-20261010-074; ``scripts/r_goldens/make_goldens_seqmd_details.R``).

- The deprecated ``what="sm"``/``"seqmc"`` warn with R's text, as a DeprecationWarning (the package's warning for
  R's deprecations); it was a UserWarning with other words.
- The combined sequence object takes seqdef's default markers, ``*`` and ``%``, as R's ``seqdef(newseqdata,
  cnames=md.cnames)`` does (seqMD.R:197); it took channel 1's.
- R also takes seqdef's default start, 1, but copies the channels' column names; traminex labels time columns from
  ``start``, so it keeps the channels' start, and its labels equal R's names (no change).
- ``cweight=[0, 0]``: with ``link="mean"`` R divides by zero and returns NaN costs (its distance route then refuses);
  traminex refuses on both routes (a divergence: a mean over zero total weight is undefined). With ``link="sum"``
  the costs are 0, as in R, and the distances are refused by traminex's zero-indel rule where R returns zeros.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqmd_details.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")["value"]
D1 = tx.seqdef(np.array([list("AABB"), list("ABBB"), list("BBAA")], dtype=object))
D2 = tx.seqdef(np.array([list("XXYY"), list("XYYY"), list("YYXX")], dtype=object))


@pytest.mark.parametrize("what", ["sm", "seqmc"])
def test_deprecated_what_values_warn_with_r_text(what: str) -> None:
    with pytest.warns(DeprecationWarning, match=re.escape(GOLDEN[f"deprecated_{what}"])):
        tx.seqMD([D1, D2], sm=["CONSTANT", "CONSTANT"], what=what)


def _combined():
    e1 = tx.seqdef(
        np.array([["A", None, "B", "B"], ["A", "B", "B", None]], dtype=object),
        nr="?",
        void="~",
        right=None,
        start=5,
    )
    e2 = tx.seqdef(np.array([list("XXYY"), list("XYYY")], dtype=object), nr="?", void="~", start=5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return tx.seqMD([e1, e2], what="MDseq", with_missing=True)


def test_the_combined_object_takes_seqdef_default_markers() -> None:
    combined = _combined()
    assert (combined.nr, combined.void) == (GOLDEN["mdseq_nr"], GOLDEN["mdseq_void"])


def test_the_combined_object_labels_time_as_r_names_its_columns() -> None:
    combined = _combined()
    labels = (
        list(combined.xtlab)
        if combined.xtlab is not None
        else [f"T{combined.start + offset}" for offset in range(combined.seqlen)]
    )
    assert labels == GOLDEN["mdseq_colnames"].split()
    assert GOLDEN["mdseq_start"] == "1"  # R's attribute, which its own names contradict


@pytest.mark.parametrize("what", ["cost", "diss"])
def test_zero_channel_weights_with_a_mean_link_are_refused(what: str) -> None:
    assert GOLDEN["cweight0_mean_cost"].split()[0] == "NaN"
    assert GOLDEN["cweight0_mean_diss"].startswith("ERROR")
    with pytest.raises(ValueError, match="non-zero sum"):
        tx.seqMD(
            [D1, D2],
            method="OM" if what == "diss" else None,
            sm=["CONSTANT", "CONSTANT"],
            cweight=[0, 0],
            link="mean",
            what=what,
        )


def test_zero_channel_weights_with_a_sum_link() -> None:
    r_costs = GOLDEN["cweight0_sum_cost"].split()
    costs = tx.seqMD([D1, D2], sm=["CONSTANT", "CONSTANT"], cweight=[0, 0], link="sum", what="cost")
    np.testing.assert_array_equal(
        np.asarray(costs["sm"], float).ravel(), [float(v) for v in r_costs[:4]]
    )
    assert float(costs["indel"]) == float(r_costs[-1])
    assert set(GOLDEN["cweight0_sum_diss"].split()) == {"0"}
    with pytest.raises(ValueError, match="indel cost must be strictly positive"):
        tx.seqMD(
            [D1, D2],
            method="OM",
            sm=["CONSTANT", "CONSTANT"],
            cweight=[0, 0],
            link="sum",
            what="diss",
        )
