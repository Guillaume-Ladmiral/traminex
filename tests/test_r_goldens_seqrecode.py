"""R-golden parity for seqrecode (T17 VERIFY) vs TraMineR 2.2-14.

seqrecode(seqdata, recodes, otherwise) collapses the alphabet. New alphabet =
names(recodes) then, if otherwise is NULL, the unlisted original states in their
original order (else `otherwise`). Goldens (scripts/r_goldens/make_goldens_seqrecode.R)
pin the new alphabet and the recoded sequence codes (seqasnum, 0-based).
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import seqrecode
from traminex.stats import seqasnum


def _codes(res):
    out = np.asarray(seqasnum(res), dtype=float).copy()
    out[np.isnan(out)] = -1
    return out


def test_seqrecode_keep_unlisted():
    res = seqrecode(fixture_a(), {"Work": ["A", "B"], "Leisure": ["C"]})
    assert list(res.alphabet) == load_golden_lines("seqrecode_keep_alph")[0].split()
    np.testing.assert_allclose(
        _codes(res), load_golden("seqrecode_keep_codes"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqrecode_otherwise():
    res = seqrecode(fixture_a(), {"Active": ["A", "B", "C"]}, otherwise="Inactive")
    assert list(res.alphabet) == load_golden_lines("seqrecode_otherwise_alph")[0].split()
    np.testing.assert_allclose(
        _codes(res), load_golden("seqrecode_otherwise_codes"), rtol=0, atol=R_GOLDEN_ATOL
    )
