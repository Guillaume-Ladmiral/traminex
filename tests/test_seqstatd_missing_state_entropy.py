"""``seqstatd(include_missing=True)`` counts the missing state even when absent.

``seqstatd.R:12-17`` appends the ``nr`` state (and ``missing.color``) to
``statl`` whenever ``with.missing=TRUE``, so the normalised entropy is divided
by ``log(K + 1)`` and the palette gains the missing colour even when no cell is
missing.  traminex counted only the codes present in the object, so on
missing-free data it divided by ``log(K)`` (``0.6309`` where R gives ``0.5`` on
a three-state alphabet; ``NaN`` where R gives ``0`` on a one-state alphabet)
and omitted the missing colour (2026-09-22 audit S3-04).  ``seqient`` already
handled this case.  The frequency matrix keeps its documented marker-padded
shape (``docs/API_PARITY_DECISIONS.csv``, seqstatd DIVERGENCE row): R's extra
``*`` row is all zeros and traminex has no column for an nr code the object
does not carry.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-04), live TraMineR 2.2-14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.sequences import STSList

GOLDENS = Path(__file__).parent / "r_goldens" / "audit20260922_s3a.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDENS.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, *_ = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _r_payload(case: str) -> list[str]:
    kind, payload = R[case]
    assert kind == "value", R[case]
    return payload.split(" ")


def _r_values(case: str) -> np.ndarray:
    return np.array(_r_payload(case), dtype=float)


def _sts(rows, **kwargs) -> STSList:
    return tx.seqdef(np.array(rows, dtype=object), **kwargs)


def _sf() -> STSList:
    return _sts([["A", "A", "B"], ["B", "C", "C"]], alphabet=["A", "B", "C"])


def _sk1() -> STSList:
    return _sts([["A", "A", "A"], ["A", "A", "A"]], alphabet=["A"])


def _close(actual, case: str) -> None:
    np.testing.assert_allclose(np.asarray(actual, dtype=float), _r_values(case), rtol=0, atol=1e-15)


@pytest.mark.parametrize(
    "case,make,kwargs",
    [
        ("s3_04_statd_entropy_wm", _sf, {}),
        ("s3_04_statd_entropy_wm_raw", _sf, {"norm": False}),
        ("s3_04_statd_entropy_k1_wm", _sk1, {}),
    ],
)
def test_entropy_counts_the_missing_state_as_r(case, make, kwargs):
    x = make()
    assert x.nr_code is None
    _close(tx.seqstatd(x, include_missing=True, **kwargs)["Entropy"], case)


def test_weighted_entropy_counts_the_missing_state_as_r():
    x = _sts(
        [["A", "A", "B"], ["B", "C", "C"], ["C", "A", "A"]],
        alphabet=["A", "B", "C"],
        weights=[1, 2.5, 0.5],
    )
    _close(tx.seqstatd(x, include_missing=True)["Entropy"], "s3_04_statd_entropy_wm_weighted")


def test_frequencies_and_valid_states_match_r_on_the_real_states():
    result = tx.seqstatd(_sf(), include_missing=True)
    r_freq = _r_values("s3_04_statd_freq_wm").reshape(3, 4)  # time x (A, B, C, *)
    np.testing.assert_allclose(result["Frequencies"], r_freq[:, :3], rtol=0, atol=1e-15)
    np.testing.assert_array_equal(r_freq[:, 3], 0.0)
    _close(result["ValidStates"], "s3_04_statd_valid_wm")


@pytest.mark.parametrize(
    "case,make", [("s3_04_statd_cpal_wm", _sf), ("s3_04_statd_cpal_k1_wm", _sk1)]
)
def test_palette_gains_the_missing_colour_as_r(case, make):
    assert tx.seqstatd(make(), include_missing=True).cpal == _r_payload(case)


def test_missing_free_and_declared_missing_objects_agree():
    """An object whose nr code is declared (in or out of the alphabet) but never
    used gives the same entropy and palette as one with no nr code at all."""
    plain = _sf()
    declared_inside = STSList(
        plain.data.copy(),
        alphabet=["A", "B", "C", "*"],
        nr="*",
        nr_code=3,
        cpal=list(plain.cpal) + [plain.missing_color],
    )
    declared_outside = STSList(
        plain.data.copy(), alphabet=["A", "B", "C"], nr="*", nr_code=3, cpal=list(plain.cpal)
    )
    expected = tx.seqstatd(plain, include_missing=True)
    for other in (declared_inside, declared_outside):
        result = tx.seqstatd(other, include_missing=True)
        np.testing.assert_allclose(result["Entropy"], expected["Entropy"], rtol=0, atol=1e-15)
        assert result.cpal == expected.cpal


def test_excluding_missing_is_unchanged():
    result = tx.seqstatd(_sf(), include_missing=False)
    np.testing.assert_allclose(result["Entropy"], np.log(2) / np.log(3), rtol=0, atol=1e-15)
    assert result.cpal == list(_sf().cpal)
