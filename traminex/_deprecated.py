"""Private helpers for source-declared deprecated argument aliases."""

from __future__ import annotations

import os
import sys
import warnings
from collections.abc import Mapping, Sequence
from types import FrameType
from typing import Any, Final

import numpy as np


class _Missing:
    def __repr__(self) -> str:
        return "..."


MISSING: Final = _Missing()

_PACKAGE_PREFIX: Final = os.path.dirname(os.path.abspath(__file__)) + os.sep


def _caller_stacklevel() -> int:
    """Return the ``warnings.warn`` stack level of the first frame outside traminex.

    Called from :func:`resolve_deprecated_alias`, where level 1 is that function
    and level 2 its caller. Some public functions are wrapped (``seqrep`` by
    ``sequences._require_stslist``), so a fixed level can land inside the
    package, where Python's default filters hide a ``DeprecationWarning``; R
    always shows the equivalent warning (audit S7-07).
    """
    level = 1
    frame: FrameType | None = sys._getframe(1)
    while frame is not None and os.path.abspath(frame.f_code.co_filename).startswith(
        _PACKAGE_PREFIX
    ):
        frame = frame.f_back
        level += 1
    return level


def _same_value(left: Any, right: Any) -> bool:
    if left is right:
        return True
    if type(left) is not type(right):
        return False
    if isinstance(left, Mapping):
        return left.keys() == right.keys() and all(
            _same_value(left[key], right[key]) for key in left
        )
    if isinstance(left, np.ndarray):
        if left.shape != right.shape or left.dtype != right.dtype:
            return False
        if left.dtype == object:
            return all(
                _same_value(left_value, right_value)
                for left_value, right_value in zip(left.flat, right.flat)
            )
        try:
            return bool(np.array_equal(left, right, equal_nan=True))
        except TypeError:
            return bool(np.array_equal(left, right))
    if isinstance(left, Sequence) and not isinstance(left, (str, bytes, bytearray)):
        return len(left) == len(right) and all(
            _same_value(left_value, right_value) for left_value, right_value in zip(left, right)
        )
    try:
        if bool(np.isscalar(left)) and bool(np.isscalar(right)):
            left_missing = bool(np.asarray(np.not_equal(left, left)).item())
            right_missing = bool(np.asarray(np.not_equal(right, right)).item())
            if left_missing or right_missing:
                return left_missing and right_missing
    except (TypeError, ValueError):
        pass
    try:
        equal = left == right
    except Exception:
        return False
    if isinstance(equal, (bool, np.bool_)):
        return bool(equal)
    try:
        return bool(np.all(equal))
    except (TypeError, ValueError):
        return False


def resolve_deprecated_alias(
    canonical_value: Any,
    alias_value: Any,
    *,
    canonical_name: str,
    alias_name: str,
    function_name: str,
    default: Any = MISSING,
) -> Any:
    """Resolve one R-style deprecated alias without silently discarding it."""
    canonical_missing = canonical_value is MISSING
    alias_missing = alias_value is MISSING
    if alias_missing:
        if canonical_missing:
            raise TypeError(
                f"{function_name}() missing required argument: {canonical_name!r} "
                f"(deprecated alias {alias_name!r} is also accepted)"
            )
        return canonical_value

    if not canonical_missing and (default is MISSING or not _same_value(canonical_value, default)):
        raise ValueError(
            f"{function_name}() cannot receive both {canonical_name!r} and its "
            f"deprecated alias {alias_name!r}"
        )
    warnings.warn(
        f"{alias_name!r} is deprecated for {function_name}(); use {canonical_name!r} instead",
        DeprecationWarning,
        stacklevel=_caller_stacklevel(),
    )
    return alias_value
