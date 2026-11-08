"""``seqefsub`` warns, as R does, that count method 6 is for internal use.

``seqefsub.R:26-28`` warns " [!] count method 6 is only for internal use."
when the constraint's ``count.method`` is 6 (the age of first occurrence,
which ``seqeapplysub(method="age")`` uses internally), then mines anyway.
traminex mined silently (2026-09-22 audit S5-09).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_event_objects.txt"
R = {
    fields[1]: fields[2]
    for fields in (line.split("\t") for line in GOLDEN.read_text().splitlines())
    if fields[0] == "count_method_6"
}


def test_count_method_six_warns_and_still_mines():
    eseq = tx.seqecreate(
        id=[1, 1, 2, 2, 2, 3, 3], timestamp=[0, 2, 0, 1, 3, 1, 4], event=list("abacbbc")
    )
    constraint = tx.seqeconstraint(count_method=6)
    with pytest.warns(UserWarning, match=re.escape(R["warning"].strip())) as caught:
        result = tx.seqefsub(eseq, min_support=1, constraint=constraint)
    assert [str(w.message) for w in caught] == [R["warning"]]
    assert result.subseq == R["subseq"].split(" ")
    np.testing.assert_array_equal(
        result.data["Count"].to_numpy(), np.array(R["count"].split(" "), dtype=float)
    )
