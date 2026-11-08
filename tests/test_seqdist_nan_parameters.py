"""NaN method parameters are refused, as R refuses them (TX-AUD-20261010-013).

R stops on a NaN ``nu``/``h`` (TWED), ``expcost`` (OMspell, OMloc), ``context``
(OMloc), ``h`` (OMslen) and ``otto`` (OMstran) -- only through its internal
"missing value where TRUE/FALSE needed", because its sign checks meet an NA.
traminex's ``x < 0`` checks were simply False for NaN, so the call returned a
matrix of NaN. Fixture: ``scripts/r_goldens/make_goldens_seqdist_nan_parameters.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "seqdist_nan_parameters.txt"
R = {
    case: (kind, payload)
    for case, kind, payload in (
        line.split("|", 2) for line in GOLDEN.read_text(encoding="utf-8").splitlines()
    )
}
SEQS = tx.seqdef(
    np.array([list(s) for s in ("AABBC", "ABBCC", "CCBAA", "BBCCA")], dtype=object),
    alphabet=["A", "B", "C"],
)
CODES = np.asarray(SEQS.data, dtype=np.int32)
NAN = float("nan")

#: (golden case, method, seqdist options, low-level method and options)
REFUSED = {
    "twed_nu": ("TWED", dict(subm="CONSTANT", indel=1, nu=NAN, h=0.5), dict(nu=NAN, h=0.5)),
    "twed_h": ("TWED", dict(subm="CONSTANT", indel=1, nu=0.5, h=NAN), dict(nu=0.5, h=NAN)),
    "omspell_expcost": ("OMspell", dict(subm="CONSTANT", indel=1, expcost=NAN), dict(expcost=NAN)),
    "omloc_expcost": ("OMloc", dict(subm="CONSTANT", expcost=NAN), dict(expcost=NAN)),
    "omloc_context": (
        "OMloc",
        dict(subm="CONSTANT", expcost=0.25, context=NAN),
        dict(expcost=0.25, context=NAN),
    ),
    "omslen_h": ("OMslen", dict(subm="CONSTANT", indel=1, h=NAN), dict(h=NAN)),
}


@pytest.mark.parametrize("case", sorted(REFUSED))
def test_nan_parameter_is_refused(case: str) -> None:
    method, options, low_level = REFUSED[case]
    with pytest.raises(ValueError, match="must not be NaN"):
        tx.seqdist(SEQS, method=method, **options)
    with pytest.raises(ValueError, match="must not be NaN"):
        pairwise(CODES, method=method.upper(), **low_level)


def test_r_refuses_the_same_nan_parameters() -> None:
    for case in [*REFUSED, "omstran_otto"]:
        assert R[case] == ("error", "missing value where TRUE/FALSE needed"), case


def test_nan_otto_is_still_refused() -> None:
    with pytest.raises(ValueError, match="otto"):
        tx.seqdist(SEQS, method="OMstran", subm="CONSTANT", indel=1, otto=NAN)


@pytest.mark.parametrize(
    ("case", "options"),
    [
        ("omspell_tpow", dict(method="OMspell", subm="CONSTANT", indel=1, tpow=NAN)),
        ("twed_nu_inf", dict(method="TWED", subm="CONSTANT", indel=1, nu=float("inf"), h=0.5)),
    ],
)
def test_r_computations_are_kept(case: str, options: dict) -> None:
    kind, payload = R[case]
    assert kind == "values"
    expected = [float(token) for token in payload.split()]
    np.testing.assert_array_equal(tx.seqdist(SEQS, **options)[0], expected)
