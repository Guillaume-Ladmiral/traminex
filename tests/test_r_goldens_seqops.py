"""R-golden parity for the RP6a sequence-operation indicators vs TraMineR 2.2-14.

seqconc / seqdecomp / seqdim / seqmpos / seqcomp / seqsubsn / seqpm. Goldens:
scripts/r_goldens/make_goldens_seqops.R (over the corestats A/B/C fixtures).
"""

from __future__ import annotations

import numpy as np
import pytest
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import seqcomp, seqconc, seqdecomp, seqdim, seqmpos, seqpm, seqsubsn


def _data_lines(name):
    """Text-golden lines with the leading ``# R: ...`` provenance comment dropped."""
    return [ln for ln in load_golden_lines(name) if not ln.startswith("#")]


def _row(x, i):
    r = x.copy()
    r.data = x.data[i : i + 1]
    return r


def _pairwise(fn, x, **kw):
    n = x.nseq
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            out[i, j] = fn(_row(x, i), _row(x, j), **kw)
    return out


# ---- seqconc ----
def test_seqconc():
    for tag, fx in (("A", fixture_a), ("B", fixture_b), ("C", fixture_c)):
        got = list(seqconc(fx(), sep="-"))
        exp = _data_lines(f"seqconc_{tag}")
        assert got == exp, tag


# ---- seqdecomp ----
def _decomp_lines(m):
    return [",".join("NA" if c is None else str(c) for c in row) for row in m]


def test_seqdecomp():
    for tag, fx in (("A", fixture_a), ("C", fixture_c)):
        strs = seqconc(fx(), sep="-")
        got = _decomp_lines(seqdecomp(strs, sep="-"))
        assert got == _data_lines(f"seqdecomp_{tag}"), tag


# ---- seqdim ----
def test_seqdim():
    for tag, fx in (("A", fixture_a), ("C", fixture_c)):
        strs = seqconc(fx(), sep="-")
        np.testing.assert_array_equal(seqdim(strs), load_golden(f"seqdim_{tag}").ravel())


def test_seqdim_on_stslist():
    a = fixture_a()
    np.testing.assert_array_equal(seqdim(a), [a.nseq, a.seqlen])


# ---- seqmpos ----
def test_seqmpos():
    np.testing.assert_allclose(
        _pairwise(seqmpos, fixture_a()), load_golden("seqmpos_A"), rtol=0, atol=R_GOLDEN_ATOL
    )
    for tag, fx, wm in (
        ("B_wm0", fixture_b, False),
        ("B_wm1", fixture_b, True),
        ("C_wm0", fixture_c, False),
        ("C_wm1", fixture_c, True),
    ):
        np.testing.assert_allclose(
            _pairwise(seqmpos, fx(), with_missing=wm),
            load_golden(f"seqmpos_{tag}"),
            rtol=0,
            atol=R_GOLDEN_ATOL,
        )


def test_seqmpos_counts_every_row_and_requires_sequence_objects():
    x = fixture_c()
    expected = int(np.sum((x.data != x.void_code) & (x.data != x.nr_code)))

    assert seqmpos(x, x) == expected
    with pytest.raises(TypeError, match="sequences must be sequence objects"):
        seqmpos(x.data[0], x.data[1])


# ---- seqcomp ----
def test_seqcomp():
    for tag, fx in (("A", fixture_a), ("C", fixture_c)):
        np.testing.assert_allclose(
            _pairwise(seqcomp, fx()), load_golden(f"seqcomp_{tag}"), rtol=0, atol=R_GOLDEN_ATOL
        )


def test_seqcomp_requires_one_row_sequence_objects():
    x = fixture_a()

    with pytest.raises(ValueError, match="exactly one sequence"):
        seqcomp(x, x)
    with pytest.raises(TypeError, match="sequences must be sequence objects"):
        seqcomp(x.data[0], x.data[1])


# ---- seqsubsn ----
def test_seqsubsn():
    # R returns an n x 1 matrix; traminex's per-sequence indicators are 1-D
    # (docs/API_CONSISTENCY.md section 1c; TX-AUD-20261010-029).
    np.testing.assert_allclose(
        seqsubsn(fixture_a(), dss=True),
        load_golden("seqsubsn_A_dss").reshape(-1),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        seqsubsn(fixture_a(), dss=False),
        load_golden("seqsubsn_A_full").reshape(-1),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        seqsubsn(fixture_c(), dss=True, with_missing=False),
        load_golden("seqsubsn_C_dss_wm0").reshape(-1),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        seqsubsn(fixture_c(), dss=True, with_missing=True),
        load_golden("seqsubsn_C_dss_wm1").reshape(-1),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


# ---- seqpm ----
def _check_pm(tag, pattern):
    lines = _data_lines(f"seqpm_{tag}")
    exp_nbocc = int(lines[0])
    exp_idx = [int(t) for t in lines[1].split()] if len(lines) > 1 and lines[1].strip() else []
    res = seqpm(fixture_a(), pattern, sep="")
    assert int(res["MTab"]["nbocc"].iloc[0]) == exp_nbocc, tag
    assert list(res["MIndex"]) == exp_idx, tag


def test_seqpm():
    _check_pm("A_BC", "BC")
    _check_pm("A_DD", "DD")
    _check_pm("A_A.C", "A.C")
