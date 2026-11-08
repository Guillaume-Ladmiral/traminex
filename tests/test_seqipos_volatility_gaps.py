"""``seqipos(index="volatility")`` across an internal gap: a recorded divergence.

TraMineR's ``seqipos.R`` (line 51) recodes every state outside ``pos.states``
and ``neg.states`` to void; with ``with.missing = FALSE`` that includes the
missing state, and when ``neg.states`` is given it can include ordinary
states. The void can then sit *inside* a sequence, but ``seqdss.R`` (lines
12-50) walks only the first ``seqlength`` cells, as if void were always
trailing padding, so the spells after the gap are lost: on ``A B * * C A``
with ``pos.states = "A"`` R reads ``p n`` instead of ``p n n p`` and returns
0.667 (audit S3-12). ``seqindic``'s ``Vpos`` inherits the same value.

traminex strips the void wherever it lies and counts every observed cell,
which is R's own ``seqivolatility`` applied to the gap-free binary sequence
``p n n p`` (0.833). The departure is deliberate and recorded in
``docs/API_PARITY.md`` and the ``seqipos`` parity row. Both R values come from
``scripts/r_goldens/make_goldens_audit20260922_s3b.R``; the test fails if a
future oracle stops losing the spells, so the record can be revisited.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s3b.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, _warnings = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _value(case: str) -> float:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return float(payload)


EXPECTED = _value("s312_seqivolatility_without_gap")


def _gap_missing():
    rows = [["A", "B", None, None, "C", "A"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])


def _gap_unlisted():
    rows = [["A", "B", "C", "C", "B", "A"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])


@pytest.mark.parametrize(
    ("case", "call"),
    [
        (
            "s312_r_seqipos_missing_gap",
            lambda: tx.seqipos(_gap_missing(), pos_states="A", index="volatility"),
        ),
        (
            "s312_r_seqipos_unlisted_gap",
            lambda: tx.seqipos(_gap_unlisted(), pos_states="A", neg_states="B", index="volatility"),
        ),
        (
            "s312_r_seqindic_vpos_missing_gap",
            lambda: tx.seqindic(_gap_missing(), indic=["vpos"], ipos_args={"pos_states": "A"})[
                "Vpos"
            ],
        ),
    ],
)
def test_volatility_counts_the_spells_after_an_internal_gap(case, call):
    result = np.asarray(call(), dtype=float)
    np.testing.assert_allclose(result, [EXPECTED], rtol=0, atol=1e-12)
    # The recorded R value loses the spells after the gap; if a future oracle
    # agrees with traminex, the divergence record should be removed.
    assert abs(_value(case) - EXPECTED) > 0.1
