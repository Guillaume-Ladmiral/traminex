"""One ``seqdist(method="DHD")`` call builds DHD's generated costs once (batch-6 review M-6).

The costs checked for seqcost's messages were built a second time inside the
kernel, doubling the time-varying transition-rate work on every DHD call with
generated costs (TX-REP-20261019-007).
"""

from __future__ import annotations

import pytest

import traminex as tx
from traminex.distances import python_kernels as pk

ROWS = [
    ["A", "A", "B", "C"],
    ["A", "B", "B", "C"],
    ["C", "C", "B", "A"],
    ["B", "C", "C", "A"],
    ["A", "A", "C", "C"],
]


@pytest.mark.parametrize("refseq", [None, 2])
@pytest.mark.parametrize("subm", [None, "TRATE"])
def test_dhd_generated_costs_built_once(subm, refseq, monkeypatch) -> None:
    seqs = tx.seqdef(ROWS)
    expected = tx.seqdist(seqs, method="DHD", subm=subm, refseq=refseq)
    calls: list[int] = []
    original = pk._dhd_costs_traminer

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(pk, "_dhd_costs_traminer", counting)
    result = tx.seqdist(seqs, method="DHD", subm=subm, refseq=refseq)
    assert len(calls) == 1
    assert (result == expected).all()
