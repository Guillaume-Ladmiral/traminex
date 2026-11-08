"""Private deterministic ordering for factor-like Python categories."""

from __future__ import annotations

import ctypes
import ctypes.util
import math
import os
import re
from collections.abc import Sequence
from functools import lru_cache
from numbers import Real
from typing import Any, cast

import numpy as np
import pandas as pd


def stable_sorted_unique(values: Sequence[object]) -> list[object]:
    """Return nonmissing unique values in R's factor order.

    Numbers sort by value; text sorts in R's collation (:func:`r_character_order`,
    which needs ICU), as ``levels(factor(x))`` does; a mix sorts by the text of
    each value in that collation. This replaced a portable case-folded key that
    differed from R for punctuation and accents (TX-AUD-20261011-015, -017).
    """
    observed: list[object] = []
    for value in values:
        missing = pd.isna(value)
        if isinstance(missing, (bool, np.bool_)) and bool(missing):
            continue
        if value not in observed:
            observed.append(value)

    if observed and all(isinstance(value, Real) for value in observed):
        return list(sorted(cast(list[Real], observed)))
    rank = r_character_rank([str(value) for value in observed])
    return sorted(observed, key=lambda value: rank[str(value)])


#: The ICU major version in the C-API symbol suffix, from the library's file name:
#: Linux libicui18n.so.N[.M], macOS libicui18n.N[.M].dylib.
_ICU_LIBRARY_VERSION = (
    re.compile(r"\.so\.(\d+)(?:\D|$)"),
    re.compile(r"libicui18n\.(\d+)(?:\.\d+)*\.dylib$"),
)


def _icu_version(library_name: str) -> str:
    """The ICU major version of the library ``library_name`` (TX-AUD-20261010-030).

    ``ctypes.util.find_library("icui18n")`` gives a soname on Linux
    (``libicui18n.so.N``) and, on macOS, the path dyld finds, often the
    unversioned ``libicui18n.dylib`` symlink, so an absolute path is read through
    its symlinks first (TX-REP-20261019-005). On Windows it can give only an
    ``icui18n[.dll]`` name, which holds no version, and is refused. Only the
    Linux form has been checked on a real host.
    """
    resolved = os.path.realpath(library_name) if os.path.isabs(library_name) else library_name
    for pattern in _ICU_LIBRARY_VERSION:
        match = pattern.search(resolved)
        if match is not None:
            return match.group(1)
    shown = repr(library_name) if resolved == library_name else f"{library_name!r} ({resolved!r})"
    raise RuntimeError(
        f"cannot determine the ICU version from the library name {shown} "
        "(expected libicui18n.so.N on Linux, or a path resolving to libicui18n.N.dylib)"
    )


@lru_cache(maxsize=1)
def _icu_root_collation_api() -> tuple[Any, Any, Any, Any]:
    """Load the versioned ICU collation functions used by the R oracle."""
    library_name = ctypes.util.find_library("icui18n")
    if library_name is None:
        raise RuntimeError(
            "R-compatible character collation requires ICU, which could not be "
            "loaded (looked for 'icui18n'). pip cannot install it: on Ubuntu "
            "24.04 run `sudo apt install libicu74`. traminex deliberately has "
            "no fallback ordering: no portable key reproduces R's collation for "
            "case, accents and punctuation, so a fallback would silently permute "
            "the state alphabet. See docs/DEV_BUILD.md."
        )
    version = _icu_version(library_name)
    library = ctypes.CDLL(library_name)
    open_collator = getattr(library, f"ucol_open_{version}")
    get_sort_key = getattr(library, f"ucol_getSortKey_{version}")
    close_collator = getattr(library, f"ucol_close_{version}")

    open_collator.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int32)]
    open_collator.restype = ctypes.c_void_p
    get_sort_key.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_uint16),
        ctypes.c_int32,
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.c_int32,
    ]
    get_sort_key.restype = ctypes.c_int32
    close_collator.argtypes = [ctypes.c_void_p]
    close_collator.restype = None
    return library, open_collator, get_sort_key, close_collator


def r_character_order(values: Sequence[str]) -> list[str]:
    """Sort unique strings with the ICU root collation used by R under C.UTF-8.

    TraMineR's factor and named-vector ordering is locale-aware. The supported
    Linux oracle resolves C.UTF-8 through ICU root collation, which differs
    materially from Python code-point ordering for case, accents, and
    punctuation. ICU is loaded through its stable C API; R is not invoked.
    """
    observed = list(dict.fromkeys(str(value) for value in values))
    if len(observed) < 2:
        return observed

    library, open_collator, get_sort_key, close_collator = _icu_root_collation_api()
    del library  # The function pointers retain the cached library handle.
    status = ctypes.c_int32(0)
    collator = open_collator(b"root", ctypes.byref(status))
    if not collator or status.value > 0:
        raise RuntimeError(f"failed to create the ICU root collator (status={status.value})")

    try:
        keyed: list[tuple[bytes, int, str]] = []
        for position, value in enumerate(observed):
            encoded = value.encode("utf-16-le")
            code_unit_count = len(encoded) // 2
            source = (
                (ctypes.c_uint16 * code_unit_count).from_buffer_copy(encoded)
                if code_unit_count
                else None
            )
            required = int(get_sort_key(collator, source, code_unit_count, None, 0))
            if required <= 0:
                raise RuntimeError("ICU returned an invalid character sort key")
            buffer = (ctypes.c_uint8 * required)()
            written = int(get_sort_key(collator, source, code_unit_count, buffer, required))
            if written != required:
                raise RuntimeError("ICU returned an inconsistent character sort key")
            keyed.append((bytes(buffer), position, value))
    finally:
        close_collator(collator)

    return [value for _, _, value in sorted(keyed)]


def r_character_rank(values: Sequence[object]) -> dict[str, int]:
    """Map each distinct string to its rank in :func:`r_character_order`.

    A sort key for rows that R orders with ``order()`` on a character column,
    which follows the session collation rather than code points.
    """
    ordered = r_character_order([str(value) for value in values])
    return {value: rank for rank, value in enumerate(ordered)}


_R_INTEGER_MAX = 2**31 - 1


def _as_r_value(value: object) -> object:
    """A whole number as R holds it: beyond R's integer range (+-2147483647) it is a double.

    So two 64-bit integers equal as doubles are one level, named as R writes the double
    (batch-12 review M-5).
    """
    if (
        isinstance(value, (int, np.integer))
        and not isinstance(value, (bool, np.bool_))
        and abs(int(value)) > _R_INTEGER_MAX
    ):
        return float(value)
    return value


def _type_mix(values: Sequence[object]) -> str:
    """The kinds of value a mixed vector holds, named in a refusal: "text and numbers",
    "numbers and booleans", "text, numbers and booleans" (batch-15 review M-4)."""
    kinds = []
    if any(isinstance(value, str) for value in values):
        kinds.append("text")
    if any(isinstance(value, Real) and not isinstance(value, (bool, np.bool_)) for value in values):
        kinds.append("numbers")
    if any(isinstance(value, (bool, np.bool_)) for value in values):
        kinds.append("booleans")
    if any(not isinstance(value, (str, Real, bool, np.bool_)) for value in values):
        kinds.append("other values")
    if len(kinds) < 2:
        return "values that are not text, numbers or booleans"
    return ", ".join(kinds[:-1]) + " and " + kinds[-1]


def r_factor_levels(
    values: Sequence[object] | np.ndarray, declared: Sequence[object] | None = None
) -> list[object]:
    """The levels R's ``factor()`` gives a grouping variable, in R's order.

    Missing values (``pd.isna``) are dropped. A declared category order (a
    pandas categorical) is kept, unused categories dropped. Otherwise booleans
    sort ``False`` before ``True``, real numbers by value, and text in R's
    collation (:func:`r_character_order`, which needs ICU), so the discrepancy
    and tree functions list text groups as R does (TX-AUD-20261010-039). Text
    mixed with numbers is refused: an R vector cannot hold both.
    """
    observed: list[object] = []
    seen: set[object] = set()
    for value in np.asarray(values, dtype=object).ravel().tolist():
        missing = pd.isna(value)
        if isinstance(missing, (bool, np.bool_)) and bool(missing):
            continue
        value = _as_r_value(value)
        if value not in seen:
            seen.add(value)
            observed.append(value)
    if declared is not None:
        return [_as_r_value(level) for level in declared if _as_r_value(level) in seen]
    if all(isinstance(value, (bool, np.bool_)) for value in observed):
        return sorted(observed, key=bool)
    if all(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)) for value in observed
    ):
        return sorted(observed, key=lambda value: cast(float, value))
    if all(isinstance(value, str) for value in observed):
        rank = r_character_rank(observed)
        return sorted(observed, key=lambda value: rank[cast(str, value)])
    # Named from the values as given: the set above merges True with 1.
    given = [
        v
        for v in np.asarray(values, dtype=object).ravel().tolist()
        if not (isinstance(pd.isna(v), (bool, np.bool_)) and bool(pd.isna(v)))
    ]
    raise TypeError(f"group levels mix {_type_mix(given)}; R's factor() takes one type")


def r_factor(
    values: Sequence[object] | np.ndarray, declared: Sequence[object] | None = None
) -> tuple[list[str], np.ndarray, list[object]]:
    """R's ``factor(values)``: its level labels, each value's level code and each level's first value.

    R takes ``levels <- unique(as.character(y)[order(y)])`` for the unique values ``y`` and
    matches ``as.character(values)`` against them, so a level is named as R prints it (``1``,
    ``1e+05``) and numbers equal to 15 significant digits share one level (``0.3`` and
    ``0.1 + 0.2``; TX-AUD-20261019-046, TX-AUD-20261011-055). The order is :func:`r_factor_levels`'s.
    Codes are 0-based, ``-1`` for a missing value; the third list holds, per level, the first
    value of :func:`r_factor_levels` that carries its label.
    """
    from ._rformat import r_as_character

    labels: list[str] = []
    firsts: list[object] = []
    for level in r_factor_levels(values, declared):
        label = cast(str, r_as_character(level))
        if label not in labels:
            labels.append(label)
            firsts.append(level)
    position = {label: index for index, label in enumerate(labels)}
    flat = np.asarray(values, dtype=object).ravel().tolist()
    codes = np.full(len(flat), -1, dtype=np.int64)
    for index, value in enumerate(flat):
        missing = pd.isna(value)
        if isinstance(missing, (bool, np.bool_)) and bool(missing):
            continue
        codes[index] = position.get(cast(str, r_as_character(_as_r_value(value))), -1)
    return labels, codes, firsts


_R_SPACE = " \t\n\v\f\r"


def r_order_keys(values: Sequence[object]) -> np.ndarray | None:
    """Sort keys reproducing R's ``order()`` on a plain vector, or None for numbers.

    Text is ranked in R's collation (ICU), missing values (None, NaN) last, as
    ``order(sortv)`` does in ``plot.stslist`` (TX-AUD-20261010-060) and in
    ``dissrf`` (dissrf.R:94; TX-REP-20261021-001). A purely
    numeric or boolean vector returns None: NumPy's stable argsort already
    orders it as R does, NaN last. An R vector holds one type, so text mixed
    with numbers is refused.
    """
    present = [value for value in values if not _is_missing(value)]
    if not any(isinstance(value, str) for value in present):
        if any(value is None for value in values):
            numbers = [math.inf if _is_missing(v) else float(cast(float, v)) for v in values]
            missing = [_is_missing(v) for v in values]
            # lexsort gives the order; the caller argsorts keys, so return ranks.
            order = np.lexsort((np.asarray(numbers), np.asarray(missing)))
            ranks = np.empty(len(values), dtype=np.int64)
            ranks[order] = np.arange(len(values))
            return ranks
        return None
    if not all(isinstance(value, str) for value in present):
        raise TypeError(f"sortv values mix {_type_mix(present)}; R's order() takes one type")
    rank = r_character_rank(present)
    last = len(rank)
    return np.asarray([last if _is_missing(v) else rank[cast(str, v)] for v in values])


def _is_missing(value: object) -> bool:
    missing = pd.isna(value)
    return isinstance(missing, (bool, np.bool_)) and bool(missing)


_R_HEX = re.compile(r"0[xX]([0-9a-fA-F]*)(?:\.([0-9a-fA-F]*))?(?:[pP]([+-]?)([0-9]+))?")
_R_DECIMAL = re.compile(r"([0-9]*)(?:\.([0-9]*))?(?:[eE]([+-]?)([0-9]+))?")


def r_string_to_double(text: str) -> float | None:
    """R's ``as.numeric`` for one string; ``None`` stands for R's ``NA``.

    Follows R's string-to-double reader (``R_strtod``): optional leading
    white space and sign, then ``NaN``, ``Inf``/``Infinity`` (any case), a
    hexadecimal number (``0x`` digits, optional fraction and binary ``p``
    exponent) or a decimal number with at least one digit and an optional
    exponent with at least one digit; only white space may follow. Unlike
    Python's ``float``, digit separators (``1_0``) are not numbers. ``"NaN"``
    gives ``math.nan``, which R's ``is.na`` also reports as missing.
    """
    body = text.lstrip(_R_SPACE)
    if not body.strip(_R_SPACE):
        return None
    sign = 1.0
    if body[0] in "+-":
        sign = -1.0 if body[0] == "-" else 1.0
        body = body[1:]

    def finished(rest: str) -> bool:
        return not rest.strip(_R_SPACE)

    folded = body[:8].casefold()
    if folded.startswith("nan"):
        return math.nan if finished(body[3:]) else None
    if folded.startswith("inf"):
        rest = body[8:] if folded == "infinity" else body[3:]
        return sign * math.inf if finished(rest) else None

    if len(body) > 2 and body[:2] in ("0x", "0X"):
        match = _R_HEX.match(body)
        assert match is not None  # the prefix alone matches
        if not finished(body[match.end() :]):
            return None
        whole, fraction, exponent_sign, exponent = match.groups()
        digits = whole + (fraction or "")
        power = (-1 if exponent_sign == "-" else 1) * int(exponent or "0")
        try:
            return sign * math.ldexp(float(int(digits or "0", 16)), power - 4 * len(fraction or ""))
        except OverflowError:
            return sign * math.inf

    match = _R_DECIMAL.match(body)
    assert match is not None  # the pattern matches the empty string
    whole, fraction, exponent_sign, exponent = match.groups()
    if not (whole or fraction) or not finished(body[match.end() :]):
        return None
    return sign * float(f"{whole or '0'}.{fraction or '0'}e{exponent_sign or ''}{exponent or '0'}")


def r_seqstatl_order(states: Sequence[str], seen: Sequence[str] = ()) -> list[str]:
    """Distinct states in the order TraMineR's ``seqstatl`` returns them.

    ``seqstatl`` sorts the distinct values and, when every one converts with
    ``as.numeric`` (``NaN`` counts as missing there), sorts them as numbers
    instead (``seqstatl.R:21-30``). ``seen`` holds further values ``seqstatl``
    also sees but that are not returned as states, such as a ``missing`` code
    that ``seqdef`` drops afterwards (``seqdef.R:59-62``): they take part in
    the numeric-or-character decision only. Equal numbers spelled differently
    keep R's collation order. ICU collation is needed only for character order
    or for such ties.
    """
    distinct = list(dict.fromkeys(str(state) for state in states))
    values = {state: r_string_to_double(state) for state in distinct}
    numbers = [value for value in values.values() if value is not None and not math.isnan(value)]
    others = (r_string_to_double(str(value)) for value in dict.fromkeys(seen))
    if len(numbers) < len(distinct) or any(value is None or math.isnan(value) for value in others):
        return r_character_order(distinct)
    tie_rank = (
        r_character_rank(distinct)
        if len(set(numbers)) < len(numbers)
        else dict.fromkeys(distinct, 0)
    )
    return sorted(distinct, key=lambda state: (cast(float, values[state]), tie_rank[state]))
