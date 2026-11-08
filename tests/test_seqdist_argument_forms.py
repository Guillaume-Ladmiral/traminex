"""Argument forms ``seqdist`` accepts, maps or refuses, as TraMineR does.

``seqdist.R`` maps the deprecated methods ``OMopt``/``LCSopt`` and a logical
``norm`` (``TRUE`` -> ``"auto"``, ``FALSE`` -> ``"none"``), uses the first two
elements of a longer ``refseq`` list, and refuses ``refseq = TRUE``. traminex
refused the first three and silently read ``True`` as row 1 (2026-09-22 audit
S1a-12). R's notices are ``message()`` calls; traminex raises a
``DeprecationWarning`` or ``UserWarning``. One difference is deliberate: R
silently drops a ``0`` inside ``refseq`` index sets (its helper counts 0 as
non-negative and R indexing skips it); traminex refuses it by name.
Expected values: ``scripts/r_goldens/make_goldens_seqdist_argument_forms.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens"
X = np.array([[0, 1, 2, 1], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 2, 0]])


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in (GOLDENS / "seqdist_argument_forms.txt").read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, *_ = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value"
    flat = np.array(payload.split(" "), dtype=float)
    size = int(round(np.sqrt(flat.size)))
    return flat.reshape(size, size).T if size * size == flat.size and size > 1 else flat


def _s():
    return tx.seqdef(X, alphabet=[0, 1, 2])


@pytest.mark.parametrize(
    "case,method,extra",
    [("omopt", "OMopt", {"subm": "CONSTANT"}), ("lcsopt", "LCSopt", {})],
)
def test_deprecated_method_names_map_to_their_equivalents(case, method, extra):
    with pytest.warns(DeprecationWarning, match="deprecated"):
        result = tx.seqdist(_s(), method=method, **extra)
    np.testing.assert_allclose(result, _values(case), atol=1e-12)


@pytest.mark.parametrize("case,norm", [("norm_true", True), ("norm_false", False)])
def test_a_logical_norm_maps_to_auto_or_none(case, norm):
    with pytest.warns(DeprecationWarning, match="norm"):
        result = tx.seqdist(_s(), method="OM", subm="CONSTANT", norm=norm)
    np.testing.assert_allclose(result, _values(case), atol=1e-12)


def test_a_longer_refseq_list_uses_its_first_two_sets():
    with pytest.warns(UserWarning, match="first two"):
        result = tx.seqdist(_s(), method="OM", subm="CONSTANT", refseq=[[1, 2], [3], [4]])
    np.testing.assert_allclose(np.ravel(result), _values("refseq_three_sets"), atol=1e-12)


def test_refseq_true_is_refused():
    assert R["refseq_true"][0] == "error"
    with pytest.raises(ValueError, match="refseq"):
        tx.seqdist(_s(), method="OM", subm="CONSTANT", refseq=True)


def test_a_zero_in_refseq_index_sets_is_refused_where_r_drops_it():
    """Deliberate: R returns the distance for the remaining index only."""
    assert R["refseq_set_with_zero"] == ("value", "2")
    with pytest.raises(ValueError, match="refseq"):
        tx.seqdist(_s(), method="OM", subm="CONSTANT", refseq=[[0, 1], [3]])
