"""R's text rendering of numbers, for labels that must read as TraMineR's do.

R chooses between fixed and scientific notation per number: with ``digits``
significant digits it keeps the fewest digits that still represent the value,
then prints fixed notation unless that is wider than scientific notation
(``getOption("scipen") = 0``).  So ``as.character(100000)`` is ``"1e+05"`` but
``as.character(123456)`` is ``"123456"``, and ``format(1234.5, digits = 5)`` is
``"1234.5"`` where C's ``%.5g`` gives ``"1234.5"`` but ``%.3g`` gives
``"1.23e+03"``.
"""

from __future__ import annotations

import datetime
import math


def r_number_text(value: float, digits: int = 15) -> str:
    """Render one number as R's ``as.character`` / ``format(digits=)`` would.

    ``digits=15`` is ``as.character`` for doubles; smaller values match
    ``format(x, digits = digits)`` for a single number.
    """
    x = float(value)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0.0:
        return "0"
    # Round to ``digits`` significant digits, then drop trailing zeros: what is
    # left is the fewest digits that represent the value at that precision.
    mantissa, exponent_text = f"{x:.{digits - 1}e}".split("e")
    if "." in mantissa:
        mantissa = mantissa.rstrip("0").rstrip(".")
    exponent = int(exponent_text)
    significant = sum(char.isdigit() for char in mantissa)
    scientific = f"{mantissa}e{'+' if exponent >= 0 else '-'}{abs(exponent):02d}"
    # Fixed notation shows every integer digit and only the decimals the
    # significant digits need (R: sprintf("%.*f", rgt, x)).
    fixed = f"{x:.{max(0, significant - 1 - exponent)}f}"
    return fixed if len(fixed) <= len(scientific) else scientific


# Constants of R's nmath/fround.c: M_LOG10_2, DBL_DIG and DBL_MAX_10_EXP.
_LOG10_2 = 0.301029995663981195213738894724493027
_DBL_DIG = 15
_MAX10E = 308


def _r_pow_di(x: float, n: int) -> float:
    """R's ``R_pow_di`` (arithmetic.c): ``x**n`` by repeated squaring of ``|n|``, inverted for ``n < 0``."""
    if math.isnan(x):
        return x
    xn = 1.0
    if n != 0:
        if not math.isfinite(x):
            return math.pow(x, n)
        is_neg = n < 0
        n = abs(n)
        while True:
            if n & 1:
                xn *= x
            n >>= 1
            if not n:
                break
            x *= x
        if is_neg:
            xn = 1.0 / xn
    return xn


def r_round(value: float, digits: float = 0.0) -> float:
    """R's ``round(x, digits)`` for one number: a port of R 4.x's ``nmath/fround.c``.

    Of ``floor(x * 10^d) / 10^d`` and ``ceiling(x * 10^d) / 10^d`` it returns the
    nearer to ``x``, and on an exact tie the one whose ``floor`` is even; so
    ``round(6.165, 2)`` is 6.16 and ``round(5e-05, 4)`` is 0, where Python's
    ``round`` gives 6.17 and 1e-04 (batch-11 review M-4).
    """
    x = float(value)
    if math.isnan(x) or math.isnan(digits):
        return x + digits
    if not math.isfinite(x) or digits > _MAX10E + _DBL_DIG or x == 0.0:
        return x
    if digits < -_MAX10E:
        return 0.0
    if digits == 0.0:
        return float(round(x))  # nearbyint: half to even
    dig = int(math.floor(digits + 0.5))
    sgn = 1.0
    if x < 0.0:
        sgn, x = -1.0, -x
    if _LOG10_2 * (0.5 + (math.frexp(x)[1] - 1)) + dig > _DBL_DIG:  # logb(x) = frexp's exponent - 1
        return sgn * x
    if dig <= _MAX10E:
        pow10 = _r_pow_di(10.0, dig)
        x10 = x * pow10
        i10 = float(math.floor(x10))
        xd = i10 / pow10
        xu = float(math.ceil(x10)) / pow10
    else:
        p10 = _r_pow_di(10.0, dig - _MAX10E)
        pow10 = _r_pow_di(10.0, _MAX10E)
        x10 = (x * pow10) * p10
        i10 = float(math.floor(x10))
        xd = i10 / pow10 / p10
        xu = float(math.ceil(x10)) / pow10 / p10
    du = xu - x
    dd = x - xd
    return sgn * (xu if (du < dd or (math.fmod(i10, 2.0) == 1.0 and du == dd)) else xd)


def r_as_character(value: object) -> str | None:
    """R's ``as.character`` for one cell; ``None`` stands for R's ``NA``."""
    if value is None:
        return None
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    try:
        import numpy as np

        if isinstance(value, np.bool_):
            return "TRUE" if bool(value) else "FALSE"
        if isinstance(value, np.integer):
            return str(int(value))
        if isinstance(value, np.floating):
            value = float(value)
    except ImportError:  # pragma: no cover - numpy is a hard dependency
        pass
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return None if math.isnan(value) else r_number_text(value)
    try:
        import pandas as pd

        if value is pd.NA or value is pd.NaT:
            return None
    except ImportError:  # pragma: no cover - pandas is a hard dependency
        pass
    if isinstance(value, datetime.datetime):
        return _r_datetime_text(value)
    return str(value)


def _r_datetime_text(value: datetime.datetime) -> str:
    """R's ``as.character`` of one date-time: its wall-clock time, without the zone.

    A midnight value reads as its date alone; otherwise the clock time follows,
    with fractional seconds but no trailing zeros (``"2005-01-01 12:30:00.5"``).
    """
    fraction = value.microsecond * 1000 + int(getattr(value, "nanosecond", 0))
    if not (value.hour or value.minute or value.second or fraction):
        return value.strftime("%Y-%m-%d")
    text = value.strftime("%Y-%m-%d %H:%M:%S")
    return text + f".{fraction:09d}".rstrip("0") if fraction else text


def r_match_arg(value: object, choices: tuple[str, ...]) -> str | None:
    """R's ``match.arg(arg, choices)`` for one value: an exact choice, else the one choice
    it begins, case-sensitive; None where R stops (no match, several, or not text).
    Used by ``seqpcfilter`` (seqpcplot.R:1100) and ``select_index_rows``
    (plot.stslist.R:16) (TX-AUD-20261011-002, TX-REP-20261023-007)."""
    if not isinstance(value, str):
        return None
    if value in choices:
        return value
    matches = [choice for choice in choices if value and choice.startswith(value)]
    return matches[0] if len(matches) == 1 else None
