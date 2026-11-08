"""Argument checks of four statistics functions follow TraMineR (audit S3-17).

* ``seqST`` passes ``type`` to ``seqivardur``, which stops unless it is 1 or 2
  (``seqivardur.R`` line 5); traminex computed type-1 values for any other
  type.
* ``seqivolatility`` stops on a non-logical ``adjust`` (``seqivolatility.R``
  line 5); traminex treated any truthy value as True.
* ``seqintegr`` stops on more than one state or on a state outside
  ``alphabet(seqdata, with.missing)`` (``seqintegr.R`` lines 9-20). traminex
  looked the state up in its internal codebook, so the missing or void marker
  was accepted without ``include_missing`` (a "proportion" of 1.5), and the
  missing state was refused with ``include_missing=True`` when no cell was
  missing (R returns zeros).
* ``seqipos`` maps the obsolete ``index="integration"`` to ``"integr"``
  (``seqipos.R`` lines 11-14, with a message); traminex refused it. It is
  accepted with a ``DeprecationWarning``, the package's form for R's notices
  about obsolete arguments.

Expected values and messages:
``scripts/r_goldens/make_goldens_audit20260922_s3b.R``.
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


def _values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return np.array([float(v) for v in payload.split(" ")])


def _error(case: str) -> str:
    kind, payload = R[case]
    assert kind == "error", (case, kind, payload)
    return payload


N = None


def _plain():
    rows = [["A", "A", "B", "C"], ["B", "B", "C", "C"], ["A", "B", "A", "B"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])


def _missing():
    rows = [["A", N, "B", "C"], ["B", "B", N, "C"], ["A", "B", "A", "B"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"], right=None)


def _ragged():
    rows = [["A", "B", "C", "C"], ["B", "B", N, N]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])


@pytest.mark.parametrize(
    ("case", "kwargs"),
    [("s317_seqST_type3", {"type": 3}), ("s317_seqST_type0_norm", {"type": 0, "norm": True})],
)
def test_seqst_refuses_a_type_other_than_1_or_2(case, kwargs):
    assert "type should be 1 or 2" in _error(case)
    with pytest.raises(ValueError, match="type should be 1 or 2"):
        tx.seqST(_plain(), **kwargs)


@pytest.mark.parametrize(
    ("case", "adjust"),
    [("s317_seqivolatility_adjust_text", "yes"), ("s317_seqivolatility_adjust_number", 1)],
)
def test_seqivolatility_refuses_a_non_logical_adjust(case, adjust):
    assert "adjust should be logical" in _error(case)
    with pytest.raises(TypeError, match="adjust should be logical"):
        tx.seqivolatility(_plain(), adjust=adjust)


@pytest.mark.parametrize(
    ("case", "make", "state", "include_missing"),
    [
        ("s317_seqintegr_nr_without_missing_flag", _missing, "*", False),
        ("s317_seqintegr_void", _ragged, "%", False),
        ("s317_seqintegr_unknown", _plain, "Z", False),
    ],
)
def test_seqintegr_refuses_a_state_outside_the_alphabet(case, make, state, include_missing):
    assert "not in the alphabet" in _error(case)
    with pytest.raises(ValueError, match="not in the alphabet"):
        tx.seqintegr(make(), state=state, include_missing=include_missing)


def test_seqintegr_refuses_more_than_one_state():
    message = _error("s317_seqintegr_two_states")
    with pytest.raises(ValueError, match=message):
        tx.seqintegr(_plain(), state=["A", "B"])


@pytest.mark.parametrize(
    ("case", "make", "state", "include_missing"),
    [
        ("s317_seqintegr_nr_with_missing_flag", _missing, "*", True),
        ("s317_seqintegr_nr_absent_from_data", _plain, "*", True),
        ("s317_seqintegr_valid", _missing, "A", False),
    ],
)
def test_seqintegr_accepts_the_states_r_accepts(case, make, state, include_missing):
    result = tx.seqintegr(make(), state=state, include_missing=include_missing)
    np.testing.assert_allclose(result, _values(case), rtol=0, atol=1e-12)


def test_seqipos_maps_the_obsolete_integration_index():
    np.testing.assert_allclose(
        _values("s317_seqipos_integration"), _values("s317_seqipos_integr"), rtol=0, atol=0
    )
    with pytest.warns(DeprecationWarning, match="index='integration' obsolete"):
        result = tx.seqipos(_plain(), pos_states="A", index="integration")
    np.testing.assert_allclose(result, _values("s317_seqipos_integration"), rtol=0, atol=1e-12)
