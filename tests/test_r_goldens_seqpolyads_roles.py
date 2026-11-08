"""seqpolyads compares its roles as R does (audit finding TX-AUD-20261019-004).

R golden: ``scripts/r_goldens/make_goldens_seqpolyads_roles.R`` (``seqpolyads_roles.csv``).
R requires the same size and identical ``alphabet(seqdata, with.missing)`` (seqpolyads.R:52-57),
and ``alphabet()`` appends the nr code whenever ``with.missing`` is TRUE (alphabet.R:9-11).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex as tx

GOLDEN = pd.read_csv(R_GOLDENS_DIR / "seqpolyads_roles.csv", dtype={"value": str})

R1 = [["a", "b", "c"], ["b", "b", "a"], ["c", "a", "a"]]
R2 = [["b", "a", "c"], ["a", "b", "b"], ["c", "c", "a"]]
R3 = [["a", "b", "c"], ["b", None, "a"], ["c", "a", "a"]]


def _r(case: str, kind: str) -> list[str]:
    rows = GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["kind"] == kind)]
    assert not rows.empty, (case, kind)
    return rows.sort_values("position")["value"].tolist()


def _role(rows, **kwargs) -> tx.STSList:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return tx.seqdef(rows, **kwargs)


def _polyads(roles, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return tx.seqpolyads(roles, T=10, method="HAM", **kwargs)


def test_default_roles_whose_states_first_appear_in_different_orders_are_computed() -> None:
    result = _polyads([_role(R1), _role(R2)])
    expected = np.asarray(_r("default_roles_first_appearance_differs", "observed.dist"), float)
    np.testing.assert_allclose(result["observed.dist"], expected, rtol=0, atol=R_GOLDEN_ATOL)


def test_with_missing_true_accepts_a_role_without_missing_cells() -> None:
    result = _polyads([_role(R1), _role(R3)], with_missing=True)
    expected = np.asarray(_r("with_missing_true_one_role_missing", "observed.dist"), float)
    np.testing.assert_allclose(result["observed.dist"], expected, rtol=0, atol=R_GOLDEN_ATOL)


def test_with_missing_false_and_a_missing_cell_is_refused() -> None:
    # R stops in seqdist; its text is recorded in the fixture. The wording of seqdist's
    # refusals is finding TX-AUD-20261019-079's scope.
    assert _r("with_missing_false_one_role_missing", "error") == [
        "'with.missing' must be TRUE when 'seqdata' or 'refseq' contain missing values"
    ]
    with pytest.raises(ValueError, match="with_missing must be true"):
        _polyads([_role(R1), _role(R3)])


def test_roles_with_the_same_states_in_different_explicit_orders_are_refused() -> None:
    message = _r("explicit_orders_differ", "error")[0]
    assert message == "All stslist objects must have same size and alphabet"
    with pytest.raises(ValueError, match=message):
        _polyads([_role(R1, alphabet=["a", "b", "c"]), _role(R2, alphabet=["b", "a", "c"])])
