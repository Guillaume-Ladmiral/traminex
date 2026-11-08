from __future__ import annotations

from collections.abc import Mapping
from numbers import Integral, Real

import numpy as np

BASE_ALLOWED_PAIRWISE_KWARGS = {
    "OM": {"indel", "subm", "norm", "indel_by_pos", "subm_by_pos"},
    "OMLOC": {"subm", "norm", "expcost", "context"},
    "HAM": {"norm", "subm"},
    "WHAM": {"weights"},
    "LCS": {"norm"},
    "LCP": {"norm"},
    "RLCP": {"norm"},
    "DUR": {"metric"},
    "JSD_TRANS": {"eps"},
    "DHD": {
        "subm",
        "norm",
        "variant",
        "smooth_window",
        "dhd_profile",
        "calibrated_params_path",
        "calibrated_params_paths",
    },
    "CHI2": {"breaks", "step", "overlap", "norm", "weights", "global_pdotj"},
    "EUCLID": {"breaks", "step", "overlap", "norm", "weights", "global_pdotj"},
    "TWED": {"nu", "h", "subm", "indel", "norm"},
    "OMSLEN": {"indel", "subm", "link", "h", "norm"},
    "OMSPELL": {"indel", "subm", "expcost", "tpow", "norm"},
    "OMSTRAN": {"indel", "subm", "otto", "transindel", "previous", "add_column", "weights", "norm"},
    "NMS": {"kweights", "prox"},
    "NMSMST": {"kweights", "tpow"},
    "SVRSPELL": {"kweights", "tpow", "prox"},
}

#: Keywords a method once accepted and now always refuses, with the reason the
#: refusal gives.  They are deliberately kept OUT of
#: ``BASE_ALLOWED_PAIRWISE_KWARGS`` so the generic "Allowed kwargs" line cannot
#: advertise a parameter the next call rejects, and they are re-admitted just
#: long enough to raise this targeted message instead of the generic one.
#: ``DUR``'s ``eps`` was removed as unused (``DURConfig`` rejects it too).
ALWAYS_REFUSED_PAIRWISE_KWARGS = {
    "DUR": {"eps": "eps is not supported for method DUR"},
}

# TraMineR seqdist refuses 'norm' for these (seqdist.R:206-208: msg.stop.impl("norm", ...)).
NORM_REFUSED_METHODS = {"NMS", "NMSMST", "SVRSPELL"}

TRATE_OM_KWARGS = {
    "trate_profile",
    "trate_smoothing",
    "trate_mode",
    "trate_normalize",
    "trate_weights",
    "calibrated_params_path",
    "calibrated_params_paths",
}


#: The only ``indel`` string ``seqdist`` accepts. ``seqdist.R:359-365`` resolves an
#: automatic indel as ``max(sm)/2`` over the whole substitution matrix; the
#: low-level ``pairwise`` is handed a resolved cost and has no matrix to resolve
#: against, so its vocabulary is deliberately narrower -- the same two-layer split
#: the ``subm`` diagnostics already document above.
SEQDIST_INDEL_STRINGS = ("auto",)


def validate_indel_value(value: object, *, layer: str) -> None:
    """Reject an ``indel`` that is neither a number nor this layer's string.

    Left to itself the value reached ``float(np.asarray(indel).item())`` and
    produced a bare ``could not convert string to float: 'atuo'`` -- naming
    neither the parameter nor the accepted values, at either layer (cycle-6
    F-009). Numeric values are not judged here: their positivity and finiteness
    are per-method contracts checked at the point of use, and this message
    deliberately avoids that vocabulary so the two refusals stay distinguishable.
    """
    if value is None or not isinstance(value, str):
        return
    if layer == "seqdist":
        if value.lower() in SEQDIST_INDEL_STRINGS:
            return
        accepted = ", ".join(repr(name) for name in SEQDIST_INDEL_STRINGS)
        raise ValueError(
            f"indel={value!r} is not a recognised spelling; seqdist accepts a "
            f"number, a per-state array of numbers, None, or {accepted} "
            "(resolved as max(subm)/2)"
        )
    raise ValueError(
        f"indel={value!r} is not a number; pairwise accepts a number or a "
        "per-state array of numbers. The string 'auto' is resolved by seqdist, "
        "which has the substitution matrix to resolve it against."
    )


#: R's 16 seqdist methods (seqdist.R:58-60), spelled as traminex dispatches them.
R_SEQDIST_METHODS = frozenset(
    {"OM", "OMLOC", "OMSLEN", "OMSPELL", "OMSTRAN", "HAM", "DHD", "CHI2", "EUCLID", "LCS", "LCP"}
    | {"RLCP", "NMS", "NMSMST", "SVRSPELL", "TWED"}
)

#: R seqdist arguments (traminex spelling) that only some methods use, with those
#: methods (seqdist.R:206-281, 604-640, 700-815). R accepts each of them for every
#: method and ignores it where it is unused (TX-AUD-20261010-071). LCS replaces
#: ``sm`` by constant costs (seqdist.R:370-374); NMS and SVRspell stop on it.
R_SEQDIST_ARGUMENT_USERS: dict[str, frozenset[str]] = {
    "indel": frozenset({"OM", "OMSLEN", "OMSPELL", "OMSTRAN", "TWED"}),
    "subm": frozenset({"OM", "OMLOC", "OMSLEN", "OMSPELL", "OMSTRAN", "HAM", "DHD", "TWED"}),
    "kweights": frozenset({"NMS", "NMSMST", "SVRSPELL"}),
    "tpow": frozenset({"OMSPELL", "NMSMST", "SVRSPELL"}),
    "expcost": frozenset({"OMLOC", "OMSPELL"}),
    "context": frozenset({"OMLOC"}),
    "link": frozenset({"OMSLEN"}),
    "h": frozenset({"OMSLEN", "TWED"}),
    "nu": frozenset({"TWED"}),
    "transindel": frozenset({"OMSTRAN"}),
    "otto": frozenset({"OMSTRAN"}),
    "previous": frozenset({"OMSTRAN"}),
    "add_column": frozenset({"OMSTRAN"}),
    "breaks": frozenset({"CHI2", "EUCLID"}),
    "step": frozenset({"CHI2", "EUCLID"}),
    "overlap": frozenset({"CHI2", "EUCLID"}),
    "global_pdotj": frozenset({"CHI2", "EUCLID"}),
    "prox": frozenset({"NMS", "SVRSPELL"}),
}

R_SEQDIST_LINKS = ("none", "mean", "gmean")
R_SEQDIST_COST_METHODS = ("TRATE", "CONSTANT", "INDELS", "INDELSLOG")


def is_r_positive_integer(value: object) -> bool:
    """TraMineR's ``is.a.positive.integer``: one number, >= 0, with no fraction.

    ``TraMineR-is_helpers.R:35-37``; a boolean is not numeric in R.
    """
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (Real, np.number)):
        return False
    number = float(value)
    return np.isfinite(number) and number >= 0 and number % 1 == 0


#: seqdist's norm spellings: R's six names and traminex's aliases, case-insensitive
#: (docs/DISTANCES_NORMALIZATION.md); whether the aliases and case folding stay is
#: with the owner (batch-9 triage R9-g).
SEQDIST_NORM_SPELLINGS = frozenset(
    {"none", "maxlength", "gmean", "maxdist", "yujianbo", "auto", "length", "max"}
)
SEQDIST_NORM_MESSAGE = (
    "norm must be one of: None, 'length', 'max', 'none', 'maxlength', 'gmean', "
    "'maxdist', 'YujianBo', 'auto' (case-insensitive)"
)


def resolve_r_seqdist_arguments(method: str, kwargs: dict, *, n_states: int | None) -> None:
    """Apply R's seqdist argument rules for one of R's methods, in place.

    R checks some arguments for every method (seqdist.R:130-190) and ignores an
    argument the method does not use (TX-AUD-20261010-071). This runs those
    checks, R's stops of NMS/SVRspell on ``sm``, and the ``check.costs`` and
    ``opt.args`` rules, then removes every argument the method does not use, so
    the per-method keyword check that follows sees only the method's own.
    """
    m = method.upper()
    if m not in R_SEQDIST_METHODS:
        return
    # R checks norm, indel, sm, prox, link and step in this order (seqdist.R:127-190), so
    # with two invalid arguments it names the earlier one (TX-REP-20261019-010).
    norm = kwargs.get("norm", None)
    # R compares the spelling as given: no spaces are trimmed (batch-9 review M-4).
    if norm is not None and str(norm).lower() not in SEQDIST_NORM_SPELLINGS:
        raise ValueError(SEQDIST_NORM_MESSAGE)
    indel = kwargs.get("indel", None)
    if indel is not None:
        if isinstance(indel, str):
            validate_indel_value(indel, layer="seqdist")
        elif isinstance(indel, (bool, np.bool_)) or not (
            isinstance(indel, (Real, np.number)) or np.ndim(indel) == 1
        ):
            raise ValueError("'indel' invalid type or unknown value")
        elif np.ndim(indel) == 1 and n_states is not None and len(indel) != n_states:
            raise ValueError("when a vector, 'indel' must contain a cost for each state")
    subm = kwargs.get("subm", None)
    if is_refused_cost_scalar(subm):
        raise ValueError("'sm' invalid type or unknown value")
    if isinstance(subm, str) and subm.upper() not in R_SEQDIST_COST_METHODS:
        raise ValueError(f"'sm' must be one of: {', '.join(R_SEQDIST_COST_METHODS)}")
    if subm is None:
        # R's sm = NULL is its default: the methods that need costs stop or
        # generate them exactly as when it is omitted.
        kwargs.pop("subm", None)
    elif m in {"NMS", "SVRSPELL"}:
        raise ValueError("use 'prox' instead of 'sm'")
    elif (
        m not in R_SEQDIST_ARGUMENT_USERS["subm"]
        and not isinstance(subm, str)
        and np.ndim(subm) < 2
    ):
        # R's sm is NULL, a matrix, an array or a method name for every method
        # (seqdist.R:153-168). A number is traminex's extension for the methods
        # that take costs; where none does, R's stop applies (batch-6 review M-3).
        raise ValueError("'sm' invalid type or unknown value")
    prox = kwargs.get("prox", None)
    if prox is not None and np.ndim(prox) != 2:
        raise ValueError("'prox' invalid type or unknown value: it must be a matrix")
    if "link" in kwargs and kwargs["link"] not in R_SEQDIST_LINKS:
        raise ValueError(f"'link' must be one of: {', '.join(R_SEQDIST_LINKS)}")
    if "step" in kwargs and not is_r_positive_integer(kwargs["step"]):
        raise ValueError("'step' must be a positive integer")
    if "check_costs" in kwargs:
        check_costs_requested(kwargs.pop("check_costs"))
    kwargs.pop("check_max_size", None)  # R's integer-index limit; traminex has none
    opt_args = kwargs.pop("opt_args", None)
    # R reads opt.args only when it is a list, and only its tokdep.coeff entry
    # (seqdist.R:690-691); anything else is ignored (batch-6 review M-2).
    if (
        m == "OMSPELL"
        and isinstance(opt_args, Mapping)
        and any(opt_args.get(key) is not None for key in ("tokdep.coeff", "tokdep_coeff"))
    ):
        raise ValueError(
            "opt_args with tokdep.coeff selects TraMineR's undocumented development "
            "method OMtspell (seqdist.Rd: 'non-documented arguments for development "
            "usage'), which traminex does not implement"
        )
    for name, users in R_SEQDIST_ARGUMENT_USERS.items():
        if name in kwargs and m not in users:
            kwargs.pop(name)


#: Why traminex keeps the cost checks R's ``check.costs = FALSE`` skips
#: (TX-AUD-20261011-041): R then computes OM with a negative indel or negative
#: substitution costs and returns negative "distances".
CHECK_COSTS_FALSE_MESSAGE = (
    "check_costs=False is not supported: traminex always checks the costs, because "
    "TraMineR with check.costs = FALSE returns negative distances for negative costs; "
    "the non-fatal cost messages can be silenced with a warnings filter on "
    "CostDiagnosticMessage"
)


#: The strings R's as.logical reads as TRUE and FALSE.
_R_TRUE_STRINGS = frozenset({"TRUE", "true", "True", "T"})
_R_FALSE_STRINGS = frozenset({"FALSE", "false", "False", "F"})


def is_refused_cost_scalar(value: object) -> bool:
    """A logical or complex scalar (a 0-d array included) given as ``sm``.

    R's sm is a matrix, an array or a method name (seqdist.R:153-168) and refuses
    anything else; traminex adds real numbers, which neither is (batch-9 review M-1, M-5).
    """
    if isinstance(value, np.ndarray):
        return value.ndim == 0 and (value.dtype == bool or np.iscomplexobj(value))
    return isinstance(value, (bool, np.bool_, complex, np.complexfloating))


def check_costs_requested(value: object) -> bool:
    """Read ``check_costs`` as R's ``if (check.costs)`` does, then refuse R's FALSE.

    TRUE and FALSE, a number (0 is FALSE) and the strings as.logical reads; any other
    string, NA (a float NaN) and NULL (None) stop with R's messages (live R:
    ``scripts/r_goldens/make_goldens_seqdist_check_order.R``; TX-REP-20261019-009). FALSE
    is refused because traminex always checks the costs (TX-AUD-20261011-041). A 0-d array
    is read as the scalar it holds (batch-9 review I-2).

    R reads ``check.costs`` only where it checks costs; traminex reads it on every route,
    as it checks the costs and refuses FALSE on every route (the seqdist ledger row lists
    the routes R ignores it on).
    """
    if isinstance(value, np.ndarray) and value.ndim == 0:
        value = value.item()
    if value is None:
        raise ValueError("check_costs: argument is of length zero")
    if isinstance(value, (bool, np.bool_)):
        reading = bool(value)
    elif isinstance(value, str):
        if value in _R_TRUE_STRINGS:
            reading = True
        elif value in _R_FALSE_STRINGS:
            reading = False
        else:
            raise ValueError("check_costs: argument is not interpretable as logical")
    elif isinstance(value, (Real, np.number)):
        if np.isnan(float(value)):
            raise ValueError("check_costs: missing value where TRUE/FALSE needed")
        reading = float(value) != 0
    else:
        raise ValueError("check_costs: argument is not interpretable as logical")
    if not reading:
        raise ValueError(CHECK_COSTS_FALSE_MESSAGE)
    return True


def normalize_dhd_variant(variant: str | None) -> str:
    variant_l = "mean" if variant is None else str(variant).lower()
    if variant_l not in {"mean", "sum"}:
        raise ValueError("variant must be 'mean' or 'sum'")
    return variant_l


def validate_smooth_window(value: object) -> int | None:
    """Validate a positive integer moving-average window."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("smooth_window must be an integer")
    window = int(value)
    if window < 1:
        raise ValueError(f"smooth_window must be positive, got {window}")
    return window


def _allowed_pairwise_kwargs(*, dispatcher: bool) -> dict[str, set[str]]:
    allowed = {method: set(values) for method, values in BASE_ALLOWED_PAIRWISE_KWARGS.items()}
    if dispatcher:
        allowed["OM"].update(TRATE_OM_KWARGS)
    return allowed


def refuse_unread_calibration_path(kwargs: dict) -> None:
    """Refuse a calibration path that no calibrated profile would read (TX-AUD-20261019-026)."""
    if any(
        kwargs.get(name) is not None
        for name in ("calibrated_params_path", "calibrated_params_paths")
    ):
        raise TypeError(
            "calibrated_params_path is read only with trate_profile='calibrated' or "
            "dhd_profile='calibrated'; without one it would be ignored"
        )


def validate_pairwise_kwargs(
    method: str, kwargs: dict, *, dispatcher: bool, layer: str = "pairwise"
) -> None:
    """Validate dispatcher keyword arguments.

    ``layer`` selects the vocabulary the diagnostics quote.  The two layers
    genuinely differ: ``seqdist`` accepts all four TraMineR cost-method strings
    (``SOURCES/TraMineR/R/seqdist.R:151``), while the low-level ``pairwise``
    accepts only ``'TRATE'``.  A single message cannot be accurate for both, and
    quoting the narrower one at the ``seqdist`` layer told a user porting
    ``sm="CONSTANT"`` that ``CONSTANT`` was unsupported when it is not.
    """
    if "normalize" in kwargs:
        raise TypeError("normalize is no longer supported; use norm='length' or norm='max'")
    if "sm" in kwargs:
        if layer == "seqdist":
            raise TypeError(
                "sm is not supported; use subm=<float|matrix|str>, where str is one of: "
                "TRATE, CONSTANT, INDELS, INDELSLOG"
            )
        raise TypeError("sm is not supported; use subm=<float|matrix|'TRATE'>")
    if "constant" in kwargs:
        raise TypeError(
            "constant is not supported; use seqsubm(..., method='CONSTANT', cval=...) and pass subm=<matrix>"
        )

    if not any(
        str(kwargs.get(profile, "default")).lower() == "calibrated"
        for profile in ("trate_profile", "dhd_profile")
    ):
        refuse_unread_calibration_path(kwargs)
    m = method.upper()
    if "norm" in kwargs and m in NORM_REFUSED_METHODS:
        raise ValueError(
            f"norm is not implemented for method {m}; TraMineR's seqdist() refuses 'norm' "
            f"for {m} (seqdist.R:206-208). Only the final sqrt transform is applied."
        )
    # Raise the specific explanation before the generic unknown-kwarg message,
    # while keeping the refused name out of the advertised set above.
    for name, reason in ALWAYS_REFUSED_PAIRWISE_KWARGS.get(m, {}).items():
        if name in kwargs:
            raise ValueError(reason)

    allowed = _allowed_pairwise_kwargs(dispatcher=dispatcher).get(m)
    if allowed is None:
        return
    unknown = sorted(set(kwargs) - allowed)
    if unknown:
        raise TypeError(
            f"Unsupported kwargs for method {m}: {unknown}. Allowed kwargs: {sorted(allowed)}"
        )

    if "indel" in kwargs:
        validate_indel_value(kwargs["indel"], layer=layer)
    if m == "DHD" and "variant" in kwargs:
        normalize_dhd_variant(kwargs["variant"])
    if m == "DHD" and "smooth_window" in kwargs:
        validate_smooth_window(kwargs["smooth_window"])
