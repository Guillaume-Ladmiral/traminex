"""``seqsurv`` orders and labels groups as R's ``factor()`` does (audit S3-22).

TraMineRextras ``seqsurv.R`` (line 38) turns ``groups`` into
``factor(groups)``, and the per-state curves follow its levels. Character
levels come in R's collation (``a``, ``A``, ``b``; ``ab``, ``aB``, ``Ab``).
traminex sorted them by code point (``A``, ``a``, ``b``), so the curves, their
labels and the default palette came in another order. The recorded plot-panel
divergence in ``docs/API_PARITY.md`` names ``seqedplot``, ``plot_SAMM`` and the
grouped transversal-entropy helper, not ``seqsurv``, whose parity row claims
that the group routes match R; analytical functions in the package
(``seqimplic``, ``seqdiff``, ``rowmode``, ``seqtab``) use R's collation through
``r_character_order``, so ``seqsurv`` does too.

Found while verifying: numeric and logical groups were labelled with Python's
``str`` (``2.0``, ``10.0``, ``True``) where R's ``as.character`` gives ``2``,
``10``, ``TRUE``. Expected values:
``scripts/r_goldens/make_goldens_audit20260922_s3c.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s3c.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, _warnings = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _text(case: str) -> list[str]:
    kind, payload = R[case]
    assert kind == "text", (case, kind, payload)
    return payload.split("|")


def _rows(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return np.array([float(v) for v in payload.split(" ")]).reshape(-1, 4)


ROWS = [
    ["A", "A", "B", "B"],
    ["B", "B", "A", "A"],
    ["A", "B", "B", "B"],
    ["B", "A", "A", "A"],
    ["A", "A", "A", "B"],
    ["B", "B", "B", "A"],
]
GROUPS = {
    "letter_case": ["b", "A", "a", "b", "A", "a"],
    "mixed_case_words": ["ab", "Ab", "aB", "ab", "Ab", "aB"],
    "fractional": [1.5, 2.0, 10.0, 2.0, 1.5, 10.0],
    "logical": [True, False, True, False, True, False],
}


@pytest.mark.parametrize("name", sorted(GROUPS))
def test_seqsurv_group_levels_labels_and_curves_match_r(name):
    x = tx.seqdef(np.array(ROWS, dtype=object), alphabet=["A", "B"])
    result = tx.seqsurv(x, groups=GROUPS[name], per_state=True)
    strata = _text(f"s322_{name}_strata")
    assert result.groups == strata
    assert result.ltext == _text(f"s322_{name}_ltext")
    assert list(result.curves) == strata
    for position, key in enumerate(strata, start=1):
        np.testing.assert_allclose(
            result.curves[key][:, :4],
            _rows(f"s322_{name}_curve{position}"),
            rtol=0,
            atol=1e-12,
            err_msg=f"curve {key}",
        )
