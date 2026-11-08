"""``r_character_order`` must reproduce R's C.UTF-8 collation (BLOCKED-3).

The cycle-4b audit recorded this path as blocked because "PyICU is genuinely
absent". That premise is wrong: ``traminex/_categorical_order.py`` never imports
PyICU. It loads the **system ICU C library** through
``ctypes.util.find_library("icui18n")`` and calls the stable C API directly, so
the Python binding is irrelevant to it. On this machine the library resolves
(``libicui18n.so.78``) and the ICU branch is the *active* one, not the unexercised
one.

Expected order: live R 4.6.1 under ``LC_COLLATE=C.UTF-8``, recorded by
``scripts/r_goldens/make_goldens_audit20260922_s10.R``:

    sort(c("b","A","a","B","Z","e-acute","e","_x","1","10","2"))
      -> _x 1 10 2 a A b B e e-acute Z

which is materially different from Python's code-point ``sorted()``
(``1 10 2 A B Z _x a b e e-acute``). Getting this wrong permutes every
collation-dependent result -- ``seqstatl``'s state list, ``rowmode``, and any
alphabet derived from observed states.

The test skips only if the ICU library is genuinely unavailable, and says so.
"""

from __future__ import annotations

import pytest
from r_goldens_utils import r_anchor

from traminex._categorical_order import r_character_order

try:
    from traminex._categorical_order import _icu_root_collation_api

    _icu_root_collation_api()
    ICU_AVAILABLE = True
    ICU_REASON = ""
except Exception as exc:  # noqa: BLE001 - availability probe
    ICU_AVAILABLE = False
    ICU_REASON = f"system ICU C library unavailable: {type(exc).__name__}: {exc}"

pytestmark = pytest.mark.skipif(not ICU_AVAILABLE, reason=ICU_REASON)

# R 4.6.1, LC_COLLATE=C.UTF-8: sort(unique(INPUT)), from tests/r_goldens/inline_r_anchors.txt.
R_SORTED = r_anchor("sort_c_utf8", "strings")
INPUT = ["b", "A", "a", "B", "Z", "é", "e", "_x", "1", "10", "2"]


def test_matches_r_c_utf8_collation():
    assert r_character_order(INPUT) == R_SORTED


def test_it_is_not_merely_python_code_point_order():
    """Guard against a silent fallback to ``sorted()``."""
    assert r_character_order(INPUT) != sorted(INPUT)


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        (["B", "a"], ["a", "B"]),
        (["b", "A"], ["A", "b"]),
        (["Z", "a"], ["a", "Z"]),
        # Case is a tertiary difference: lower case sorts before upper case.
        (["A", "a"], ["a", "A"]),
    ],
)
def test_case_ordering_follows_icu_not_ascii(values, expected):
    assert r_character_order(values) == expected


def test_short_inputs_are_returned_unchanged():
    assert r_character_order([]) == []
    assert r_character_order(["only"]) == ["only"]
    assert r_character_order(["dup", "dup"]) == ["dup"]
