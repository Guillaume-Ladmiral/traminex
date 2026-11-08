"""CAT multidomain sequences and costs derived from TraMineR 2.2-14 ``seqMD``."""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Optional, Sequence, cast

import numpy as np

from ._categorical_order import r_seqstatl_order
from ._deprecated import _caller_stacklevel
from .distances import seqcost, seqdist
from .distances._validation import check_costs_requested
from .distances.python_kernels import CostDiagnosticMessage, _warn_cost_quality
from .format import seqdef
from .sequences import STSList

_METLIST = ("OM", "LCS", "DHD", "HAM")
_SM_METHODS = ("CONSTANT", "TRATE", "INDELS", "INDELSLOG")
_WHATLIST = ("MDseq", "cost", "diss")


@dataclass(frozen=True)
class _ChannelDecoder:
    labels: tuple[str, ...]
    kinds: tuple[str, ...]
    cost_states: tuple[str, ...]
    state_lookup: dict[str, int]
    has_missing: bool


@dataclass(frozen=True)
class _CombinedRegistry:
    sequence: STSList
    alphabet: tuple[str, ...]
    components: tuple[tuple[str, ...], ...]


def _channel_decoder(channel: STSList) -> _ChannelDecoder:
    labels: list[str] = []
    kinds: list[str] = []
    real_states: list[str] = []
    for code, state in enumerate(channel.alphabet):
        label = str(state)
        if channel.nr_code is not None and code == int(channel.nr_code):
            kind = "missing"
            label = str(channel.nr)
        elif channel.void_code is not None and code == int(channel.void_code):
            kind = "void"
            label = str(channel.void)
        else:
            kind = "state"
            real_states.append(label)
        labels.append(label)
        kinds.append(kind)

    has_missing = bool(channel.nr_code is not None and np.any(channel.data == int(channel.nr_code)))
    cost_states = real_states + ([str(channel.nr)] if has_missing else [])
    return _ChannelDecoder(
        labels=tuple(labels),
        kinds=tuple(kinds),
        cost_states=tuple(cost_states),
        state_lookup={state: index for index, state in enumerate(cost_states)},
        has_missing=has_missing,
    )


def _effective_lengths(channel: STSList) -> np.ndarray:
    if channel.void_code is None:
        return np.full(channel.nseq, channel.seqlen, dtype=int)
    return np.sum(channel.data != int(channel.void_code), axis=1, dtype=int)


def _combined_registry(
    channels: Sequence[STSList],
    decoders: Sequence[_ChannelDecoder],
    *,
    ch_sep: str,
    fill_with_miss: bool,
) -> _CombinedRegistry:
    ncase = channels[0].nseq
    width = max(channel.seqlen for channel in channels)
    tuples = np.empty((ncase, width), dtype=object)
    labels = np.empty((ncase, width), dtype=object)

    for row in range(ncase):
        for time in range(width):
            components: list[str] = []
            all_void = True
            for channel, decoder in zip(channels, decoders):
                if time >= channel.seqlen:
                    # seqMD.R:169-172 pads a shorter channel with void, and
                    # fill.with.miss (lines 175-178) touches only the
                    # channel's own columns (audit S2-11).
                    component = str(channel.void)
                    kind = "void"
                else:
                    code = int(channel.data[row, time])
                    if code < 0 or code >= len(decoder.labels):
                        raise ValueError("channel contains a state code outside its alphabet")
                    component = decoder.labels[code]
                    kind = decoder.kinds[code]
                    if fill_with_miss and decoder.has_missing and kind == "void":
                        component = str(channel.nr)
                        kind = "missing"
                all_void = all_void and kind == "void"
                components.append(component)
            if all_void:
                tuples[row, time] = None
                labels[row, time] = None
            else:
                key = tuple(components)
                tuples[row, time] = key
                labels[row, time] = ch_sep.join(key)

    # seqMD.R:197 takes seqdef's default alphabet: seqstatl's order, numbers by
    # value when every combined label is an R number (ch_sep "."), R's collation
    # otherwise (audit S2-08; TX-AUD-20261019-038).
    observed = r_seqstatl_order(sorted({str(value) for value in labels.flat if value is not None}))
    component_by_label = {
        ch_sep.join(value): value for value in tuples.flat if isinstance(value, tuple)
    }
    longest = max(channels, key=lambda channel: channel.seqlen)
    combined = seqdef(
        labels,
        alphabet=observed,
        id=None if channels[0].ids is None else list(channels[0].ids),
        start=longest.start,
        cnames=None if longest.xtlab is None else list(longest.xtlab),
        # seqdef's default markers, as R's seqdef(newseqdata, cnames=md.cnames)
        # (seqMD.R:197) takes them; the start stays the channels', from which
        # traminex labels the columns R names after the channels (TX-AUD-20261010-074).
    )
    return _CombinedRegistry(
        sequence=combined,
        alphabet=tuple(observed),
        components=tuple(component_by_label[label] for label in observed),
    )


def _expand_per_channel(value: object, nchannels: int, *, name: str) -> list[object]:
    if isinstance(value, str) or np.isscalar(value):
        return [value] * nchannels
    if isinstance(value, np.ndarray) and value.ndim >= 2:
        raise ValueError(f"{name} must provide one value per domain")
    if not isinstance(value, Iterable):
        raise ValueError(f"{name} must provide one value per domain")
    expanded = list(value)
    if len(expanded) != nchannels:
        raise ValueError(f"{name} must provide one value per domain")
    return expanded


def _resolve_missing_declarations(
    with_missing: object,
    decoders: Sequence[_ChannelDecoder],
) -> list[bool] | None:
    """Warn about declarations that disagree with the data; return them."""
    if with_missing is None:
        return None
    if isinstance(with_missing, (bool, np.bool_)):
        declared = [bool(with_missing)] * len(decoders)
    else:
        if not isinstance(with_missing, Iterable):
            raise ValueError("with_missing must provide one value per domain")
        declared = [bool(value) for value in with_missing]
        if len(declared) != len(decoders):
            raise ValueError("with_missing must provide one value per domain")
    for index, (value, decoder) in enumerate(zip(declared, decoders), start=1):
        if value != decoder.has_missing:
            warnings.warn(
                f"Bad with.missing value for domain {index}. " f"I set it as {decoder.has_missing}",
                UserWarning,
                stacklevel=3,
            )
    return declared


def _channel_cost(
    channel: STSList,
    decoder: _ChannelDecoder,
    sm_spec: object,
    indel_spec: object,
    *,
    cval: float,
    miss_cost: float,
    time_varying: bool,
    width: int,
) -> tuple[np.ndarray, float | np.ndarray]:
    nstates = len(decoder.cost_states)
    if isinstance(sm_spec, str):
        if sm_spec not in _SM_METHODS:
            raise ValueError(f"unknown substitution-cost method {sm_spec!r}")
        # R calls suppressMessages(seqcost(..., check.costs = FALSE)) here and
        # checks the combined costs itself (seqMD.R:241-243, 266).
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", CostDiagnosticMessage)
            generated = seqcost(
                channel,
                sm_spec,
                cval=cval,
                with_missing=decoder.has_missing,
                miss_cost=miss_cost,
                time_varying=time_varying,
            )
        matrix = np.asarray(generated["sm"], dtype=float)
        generated_indel = generated["indel"]
    else:
        matrix = np.asarray(sm_spec, dtype=float)
        generated_indel = 0.5 * float(np.max(matrix))

    expected = (nstates, nstates, width) if time_varying else (nstates, nstates)
    if matrix.shape != expected:
        shape = "x".join(str(size) for size in expected)
        raise ValueError(f"[!] size of substitution cost matrix must be {shape}")

    if isinstance(indel_spec, str):
        if indel_spec != "auto":
            raise ValueError("indel string must be 'auto'")
        chosen = generated_indel
    else:
        chosen = indel_spec
    if np.isscalar(chosen):
        return matrix, float(cast(Any, chosen))
    values = np.asarray(chosen, dtype=float)
    if values.ndim == 2 and time_varying and values.shape == (nstates, width):
        return matrix, values
    if values.ndim != 1 or values.size != nstates:
        raise ValueError(
            "[!] indel length does not match size of alphabet for at least one channel"
        )
    return matrix, values


def _cat_costs(
    registry: _CombinedRegistry,
    decoders: Sequence[_ChannelDecoder],
    channel_costs: Sequence[np.ndarray],
    channel_indels: Sequence[float | np.ndarray],
    channel_weights: np.ndarray,
    *,
    time_varying: bool,
    link: str,
) -> dict[str, object]:
    # A component with no cost entry -- a void component, from a shorter
    # channel's padding or a channel without missing values -- is R's
    # match() NA: every cost pairing that state with another is NA, while its
    # diagonal keeps R's initial 0 (seqMD.R:289-335; audit S2-19). The
    # position one past each channel's cost states holds that NaN.
    indexes = np.empty((len(registry.alphabet), len(decoders)), dtype=int)
    for state, components in enumerate(registry.components):
        for channel, component in enumerate(components):
            lookup = decoders[channel].state_lookup
            indexes[state, channel] = lookup.get(component, len(lookup))
    undefined = np.zeros(len(registry.alphabet), dtype=bool)
    for channel, decoder in enumerate(decoders):
        undefined |= indexes[:, channel] == len(decoder.state_lookup)

    ncombined = len(registry.alphabet)
    if time_varying:
        width = registry.sequence.seqlen
        combined_sm = np.zeros((ncombined, ncombined, width), dtype=float)
    else:
        combined_sm = np.zeros((ncombined, ncombined), dtype=float)
    for channel, (matrix, weight) in enumerate(zip(channel_costs, channel_weights)):
        size = matrix.shape[0]
        padded = np.full((size + 1, size + 1, *matrix.shape[2:]), np.nan, dtype=float)
        padded[:size, :size] = matrix
        state_indexes = indexes[:, channel]
        if time_varying:
            contribution = padded[state_indexes[:, None], state_indexes[None, :], :]
        else:
            contribution = padded[np.ix_(state_indexes, state_indexes)]
        combined_sm += float(weight) * contribution
    # seqMD.R:307-316 computes the cost of each pair i < j once and mirrors it,
    # so an asymmetric channel matrix still gives symmetric CAT costs (audit
    # S10-13). The same holds at every time position.
    upper_rows, upper_cols = np.triu_indices(ncombined, 1)
    combined_sm[upper_cols, upper_rows] = combined_sm[upper_rows, upper_cols]
    undefined_states = np.flatnonzero(undefined)
    combined_sm[undefined_states, undefined_states] = 0.0

    combined_indel: float | np.ndarray | None
    if any(not np.isscalar(value) for value in channel_indels):
        if time_varying:
            combined_indel = None
        else:
            combined_indel = np.zeros(ncombined, dtype=float)
            for channel, (value, weight) in enumerate(zip(channel_indels, channel_weights)):
                resolved = (
                    np.full(
                        len(decoders[channel].cost_states),
                        float(cast(Any, value)),
                        dtype=float,
                    )
                    if np.isscalar(value)
                    else np.asarray(value, dtype=float)
                )
                resolved = np.append(resolved, np.nan)
                combined_indel += float(weight) * resolved[indexes[:, channel]]
    else:
        combined_indel = float(
            sum(
                float(weight) * float(value)
                for value, weight in zip(channel_indels, channel_weights)
            )
        )

    if link == "mean":
        divisor = float(np.sum(channel_weights))
        if divisor == 0:
            raise ValueError("channel weights must have a non-zero sum for link='mean'")
        combined_sm /= divisor
        if combined_indel is not None:
            combined_indel = combined_indel / divisor
    return {
        "sm": combined_sm,
        "indel": combined_indel,
        "alphabet": list(registry.alphabet),
        "cweight": channel_weights.copy(),
    }


def seqMD(
    channels: Sequence[STSList],
    method: Optional[str] = None,
    norm: str = "none",
    indel="auto",
    sm=None,
    with_missing=None,
    full_matrix: bool = True,
    link: str = "sum",
    cval: float = 2.0,
    miss_cost: float = 2.0,
    cweight: Optional[Sequence[float]] = None,
    what: str = "MDseq",
    ch_sep: str = "+",
    fill_with_miss: bool = True,
    check_costs: bool = True,
):
    """Build combined-domain sequences, CAT costs, or CAT-based distances.

    ``check_costs`` is R's ``check.costs``: ``True`` (R's default) is accepted,
    and the costs are always checked; ``False`` is refused, because R then
    returns negative distances for negative costs (TX-AUD-20261011-041).

    The combined alphabet is ``seqdef``'s default for the combined states, as in
    ``seqMD.R``: sorted in R's collation order (``a+x, B+x, B+Y``), which also
    orders the rows and columns of ``what="cost"``. A channel shorter than the
    longest is padded with void; ``fill_with_miss`` turns a channel's own void
    cells into its nonresponse marker but not that padding, as in R.
    With ``what="cost"`` or ``"diss"``, a channel declared
    ``with_missing=False`` that has missing values is refused, as in R;
    ``what="MDseq"`` only warns.

    A combined state with a void component (the padding of a shorter channel,
    or a trailing void in a channel without missing values) has no CAT cost:
    as in R, ``what="cost"`` gives NaN (R's NA) for every cost pairing it with
    another state and for its entry of a per-state indel vector, keeps 0 on
    its diagonal, and ``what="diss"`` refuses with R's message.

    With ``what="cost"`` or ``"diss"`` each channel's costs, explicit or
    generated from a method string, get TraMineR's ``checkcost`` diagnostics
    in R's order, as ``seqMD.R`` checks them (docs/API_PARITY.md). With
    ``what="cost"``, ``method`` does not make the costs time-varying; only the
    dimensions of ``sm`` do, as in R.
    """
    check_costs_requested(check_costs)
    channel_list = list(channels)
    if len(channel_list) < 2:
        raise ValueError("[!] please specify at least two domains")
    if any(not isinstance(channel, STSList) for channel in channel_list):
        raise TypeError("each channel must be a state sequence object")
    # Both marker representations are equivalent inputs; fold declared
    # external nr/void codes into each channel's alphabet before decoding
    # (identity for in-alphabet channels).
    channel_list = [channel._canonicalized() for channel in channel_list]
    if not isinstance(ch_sep, str) or not ch_sep:
        raise ValueError("ch_sep must be a non-empty string")
    if link not in {"sum", "mean"}:
        raise ValueError("link must be one of 'sum' or 'mean'")

    if what == "sm":
        # R's deprecation (seqMD.R:17-24) as the package warns R's deprecations
        # (TX-AUD-20261010-074).
        warnings.warn(
            " [!!] what='sm' deprecated! Use what='cost' instead.",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        what = "cost"
    elif what == "seqmc":
        warnings.warn(
            " [!!] what='seqmc' deprecated! Use what='MDseq' instead.",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        what = "MDseq"
    if what not in _WHATLIST:
        raise ValueError("what should be one of 'MDseq','cost','diss'")
    if what == "diss" and method is None:
        raise ValueError("a valid method must be provided when what = 'diss'")
    if what == "cost" and sm is None:
        raise ValueError("sm cannot be NULL when what = 'cost'")

    ncase = channel_list[0].nseq
    if any(channel.nseq != ncase for channel in channel_list):
        raise ValueError("[!] sequence objects have different numbers of rows")
    decoders = [_channel_decoder(channel) for channel in channel_list]
    for channel, decoder in zip(channel_list, decoders):
        symbols = set(decoder.labels) | {str(channel.nr), str(channel.void)}
        if any(ch_sep in symbol for symbol in symbols):
            raise ValueError("[!] ch_sep symbol occurs in alphabet of at least one channel")
    declared_missing = _resolve_missing_declarations(with_missing, decoders)

    lengths = np.vstack([_effective_lengths(channel) for channel in channel_list])
    if np.any(lengths != lengths[0]):
        warnings.warn(
            "Cases with sequences of different length across domains!",
            UserWarning,
            stacklevel=2,
        )
    registry = _combined_registry(
        channel_list,
        decoders,
        ch_sep=ch_sep,
        fill_with_miss=bool(fill_with_miss),
    )
    if what == "MDseq":
        return registry.sequence

    assert what in {"cost", "diss"}
    resolved_method = None if method is None else str(method).upper()
    if what == "diss":
        if resolved_method not in _METLIST:
            raise ValueError(f"[!] method must be one of: {' '.join(_METLIST)}")
        if resolved_method == "LCS":
            resolved_method = "OM"
            sm = "CONSTANT"
            indel = 1
            cval = 2.0
            miss_cost = 2.0
        elif sm is None:
            sm = "TRATE" if resolved_method == "DHD" else "CONSTANT"

    # seqMD.R:216-223: costs and distances stop on a channel declared without
    # missing values that has some; only the sequences are built (audit S2-12).
    if declared_missing is not None:
        for index, (declared, decoder) in enumerate(zip(declared_missing, decoders), start=1):
            if not declared and decoder.has_missing:
                raise ValueError(
                    f"[!] found missing values in channel {index}, "
                    "set with_missing as True for that channel"
                )

    nchannels = len(channel_list)
    sm_specs = _expand_per_channel(sm, nchannels, name="sm")
    indel_specs = _expand_per_channel(indel, nchannels, name="indel")
    if any(isinstance(value, str) and value == "auto" for value in indel_specs) and not all(
        isinstance(value, str) and value == "auto" for value in indel_specs
    ):
        raise ValueError("'auto' not allowed in vector or list indel")

    if cweight is None:
        channel_weights = np.ones(nchannels, dtype=float)
    else:
        channel_weights = np.asarray(cweight, dtype=float)
        if channel_weights.ndim != 1 or channel_weights.size != nchannels:
            raise ValueError("cweight must provide one weight per domain")
        if not np.all(np.isfinite(channel_weights)):
            raise ValueError("cweight must contain only finite values")

    dimensions = [np.asarray(value).ndim for value in sm_specs if not isinstance(value, str)]
    if any(dimension == 3 for dimension in dimensions) and not all(
        dimension == 3 for dimension in dimensions
    ):
        raise ValueError("One sm is 3-dimensional and some are not!")
    # seqMD.R:37 decides from sm's dimensions; only what="diss" turns DHD
    # time-varying (seqMD.R:93), so what="cost" ignores method (audit S10-12).
    time_varying = (what == "diss" and resolved_method == "DHD") or bool(
        dimensions and all(d == 3 for d in dimensions)
    )
    width = registry.sequence.seqlen
    channel_costs: list[np.ndarray] = []
    channel_indels: list[float | np.ndarray] = []
    for channel, decoder, sm_spec, indel_spec in zip(
        channel_list,
        decoders,
        sm_specs,
        indel_specs,
    ):
        matrix, channel_indel = _channel_cost(
            channel,
            decoder,
            sm_spec,
            indel_spec,
            cval=float(cval),
            miss_cost=float(miss_cost),
            time_varying=time_varying,
            width=width,
        )
        # seqMD.R:266-267 checks every channel's costs before weighting them,
        # explicit or generated from a method string (seqMD.R:241-243 builds
        # those with check.costs = FALSE because they are checked here); the
        # call is seqMD's own, in 2.2-13 too (audit S10-11).
        _warn_cost_quality(
            np.moveaxis(matrix, 2, 0) if time_varying else matrix,
            indel=channel_indel,
            labels=decoder.cost_states,
            symmetry=True,
            stacklevel=3,
        )
        channel_costs.append(matrix)
        channel_indels.append(channel_indel)

    costs = _cat_costs(
        registry,
        decoders,
        channel_costs,
        channel_indels,
        channel_weights,
        time_varying=time_varying,
        link=link,
    )
    if what == "cost":
        return costs

    combined_sm = np.asarray(costs["sm"], dtype=float)
    combined_indel = costs["indel"]
    if np.any(np.isnan(combined_sm)) or (
        combined_indel is not None and np.any(np.isnan(np.asarray(combined_indel, dtype=float)))
    ):
        raise ValueError("NA indel and/or substitution CAT costs, cannot compute MD distances!")
    backend_sm = np.moveaxis(combined_sm, 2, 0) if time_varying else combined_sm
    distance_options: dict[str, object] = {
        "norm": norm,
        "subm": backend_sm,
        "with_missing": False,
        "full_matrix": full_matrix,
    }
    if resolved_method == "OM":
        distance_options["indel"] = combined_indel
    assert resolved_method is not None
    return seqdist(registry.sequence, method=resolved_method, **distance_options)


def seqdistmc(
    channels: Sequence[STSList],
    what: str = "diss",
    ch_sep: str = "@@@@TraMineRSep@@@@",
    **kwargs,
):
    """Backward-compatible alias for :func:`seqMD`."""
    return seqMD(channels, what=what, ch_sep=ch_sep, **kwargs)
