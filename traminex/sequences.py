from __future__ import annotations

import inspect
from collections.abc import Iterable
from dataclasses import dataclass
from functools import wraps
from numbers import Integral, Real
from typing import Callable, List, Literal, Optional, ParamSpec, TypeVar, cast

import numpy as np


def _validate_sequence_weights(
    weights: object,
    n_sequences: int,
    *,
    name: str = "weights",
) -> np.ndarray:
    """Return a finite, non-negative case-weight vector of the required length.

    The result is always a new array: a sequence object must not share its
    weights with the caller's array (audit S2-13).
    """
    values = np.array(weights, dtype=float)
    if values.ndim != 1 or values.shape[0] != n_sequences:
        raise ValueError(f"{name} must have length {n_sequences}")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite numeric values")
    if np.any(values < 0):
        raise ValueError(f"{name} must be non-negative")
    return values


def _validate_analytical_weights(
    values: object,
    *,
    name: str,
    n_channels: int,
) -> List[float]:
    """Return a valid non-negative analytical channel-weight vector."""
    raw = np.asarray(values, dtype=object)
    if raw.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if raw.shape[0] != n_channels:
        raise ValueError(f"{name} must have {n_channels} elements")
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) for value in raw):
        raise TypeError(f"{name} must contain numeric values")
    weights = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(weights)):
        raise ValueError(f"{name} must contain finite numbers")
    if np.any(weights < 0):
        raise ValueError(f"{name} must be non-negative")
    if not np.any(weights > 0):
        raise ValueError(f"{name} must have positive sum")
    return weights.tolist()


@dataclass
class STSList:
    """State sequence list.

    Stores a 2D numpy array of integer state codes with attributes similar to TraMineR.
    - data: shape (n_sequences, n_positions), dtype=int
    - alphabet: list of state codes (strings)
    - labels: optional long labels, aligned with alphabet
    - cpal: optional color palette (hex strings)
    - xtlab/xtstep/tick_last: axis info for plotting
    - nr/void: symbols for missing and void and their codes, if available;
      declared marker codes may lie outside the ordinary alphabet
    - weights: optional per-sequence weights
    - start: starting time (int)
    - ids: optional row identities, aligned with sequences
    - missing_color: display color for missing states
    - issues: optional issue categories mapped to 1-based sequence indexes
    """

    data: np.ndarray
    alphabet: List[str]
    labels: Optional[List[str]] = None
    cpal: Optional[List[str]] = None
    xtlab: Optional[List[str]] = None
    xtstep: int = 1
    tick_last: bool = False
    nr: str = "*"
    void: str = "%"
    nr_code: Optional[int] = None
    void_code: Optional[int] = None
    weights: Optional[np.ndarray] = None
    start: int = 1
    ids: Optional[List[object]] = None
    missing_color: str = "darkgrey"
    issues: Optional[dict[str, list[int]]] = None

    def __setattr__(self, name: str, value: object) -> None:
        data = self.__dict__.get("data")
        if (
            name == "weights"
            and value is not None
            and isinstance(data, np.ndarray)
            and data.ndim == 2
        ):
            value = _validate_sequence_weights(
                value,
                data.shape[0],
                name="weights",
            )
        super().__setattr__(name, value)

    def __post_init__(self) -> None:
        if not isinstance(self.data, np.ndarray):
            raise TypeError("data must be a numpy array")
        if self.data.ndim != 2:
            raise ValueError("data must be 2D (n_sequences, n_positions)")
        if np.issubdtype(self.data.dtype, np.bool_) or not np.issubdtype(
            self.data.dtype,
            np.integer,
        ):
            raise TypeError("data must have an integral dtype")

        self.alphabet = list(self.alphabet)
        original_alphabet_size = len(self.alphabet)
        # Two codes spelled the same could not be told apart, and time labels
        # that do not line up with the positions would be silently shifted
        # (audit S2-15; seqdef refuses both). A cpal for the real states only
        # stays accepted: visualization.palette.ensure_cmap pads it.
        duplicated = [
            state
            for position, state in enumerate(self.alphabet)
            if state in self.alphabet[:position]
        ]
        if duplicated:
            raise ValueError(f"duplicated element(s) in alphabet: {duplicated[0]}")
        if self.labels is not None and len(self.labels) != original_alphabet_size:
            raise ValueError("labels must match alphabet length")
        if self.labels is not None:
            self.labels = list(self.labels)
        if self.xtlab is not None and len(self.xtlab) != self.data.shape[1]:
            raise ValueError(f"xtlab must have length {self.data.shape[1]}")
        if not isinstance(self.missing_color, str):
            raise TypeError("missing_color must be a string")

        if self.nr is not None and self.nr == self.void:
            # seqdef's rule: one text cannot be both missing and void (TX-AUD-20261019-042).
            raise ValueError("'nr' and 'void' must be distinct character markers")
        marker_values: dict[int, tuple[str, str]] = {}
        for attribute, marker in (("nr_code", self.nr), ("void_code", self.void)):
            value = getattr(self, attribute)
            if value is None:
                continue
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
                raise TypeError(f"{attribute} must be an integer")
            code = int(value)
            if code < 0:
                raise ValueError(f"{attribute} must be non-negative")
            setattr(self, attribute, code)
            if code in marker_values:
                raise ValueError("nr_code and void_code must be distinct")
            marker_values[code] = (attribute, marker)

        if self.data.size:
            if np.any(self.data < 0):
                raise ValueError("data state codes must be non-negative")
            invalid_mask = self.data >= len(self.alphabet)
            for marker_code in marker_values:
                invalid_mask &= self.data != marker_code
            invalid_codes = self.data[invalid_mask]
            if invalid_codes.size:
                raise ValueError(
                    f"data contains undeclared state code {int(np.min(invalid_codes))}"
                )
        if self.weights is not None:
            self.weights = _validate_sequence_weights(
                self.weights,
                self.data.shape[0],
                name="weights",
            )
        if self.ids is not None:
            if len(self.ids) != self.data.shape[0]:
                raise ValueError(f"ids must have length {self.data.shape[0]}")
            self.ids = list(self.ids)
        if self.issues is not None:
            if not isinstance(self.issues, dict):
                raise TypeError("issues must be a dictionary")
            normalized_issues: dict[str, list[int]] = {}
            for category, indexes in self.issues.items():
                if not isinstance(category, str):
                    raise TypeError("issues keys must be strings")
                if not isinstance(indexes, list):
                    raise TypeError("issues values must be lists of sequence indexes")
                normalized_indexes: list[int] = []
                for index in indexes:
                    if isinstance(index, bool) or not isinstance(index, Integral):
                        raise TypeError("issue indexes must be integers")
                    normalized_index = int(index)
                    if not 1 <= normalized_index <= self.data.shape[0]:
                        raise ValueError(
                            f"issue indexes must be between 1 and {self.data.shape[0]}"
                        )
                    normalized_indexes.append(normalized_index)
                normalized_issues[category] = normalized_indexes
            self.issues = normalized_issues

    @property
    def nseq(self) -> int:
        return self.data.shape[0]

    @property
    def seqlen(self) -> int:
        return self.data.shape[1]

    def copy(self) -> "STSList":
        return STSList(
            data=self.data.copy(),
            alphabet=list(self.alphabet),
            labels=None if self.labels is None else list(self.labels),
            cpal=None if self.cpal is None else list(self.cpal),
            xtlab=None if self.xtlab is None else list(self.xtlab),
            xtstep=self.xtstep,
            tick_last=self.tick_last,
            nr=self.nr,
            void=self.void,
            nr_code=self.nr_code,
            void_code=self.void_code,
            weights=None if self.weights is None else self.weights.copy(),
            start=self.start,
            ids=None if self.ids is None else list(self.ids),
            missing_color=self.missing_color,
            issues=self.issues,
        )

    def to_strings(self, sep: str = "-") -> List[str]:
        """Return compressed SPS-style strings using current alphabet order.

        Note: This is for display/testing, not a full SPS writer.
        """
        out: List[str] = []
        for row in self.data:
            out.append(sep.join(self._state_for_code(int(code)) for code in row))
        return out

    def _state_for_code(self, code: int) -> str:
        """Decode an already-validated state or marker code."""
        if self.nr_code is not None and code == self.nr_code:
            return self.nr
        if self.void_code is not None and code == self.void_code:
            return self.void
        return self.alphabet[code]

    def _codebook_codes(
        self,
        *,
        include_external_nr: bool = True,
        include_external_void: bool = True,
    ) -> List[int]:
        """Return ordinary positions followed by requested external marker codes."""
        codes = list(range(len(self.alphabet)))
        marker_codes = (
            self.nr_code if include_external_nr else None,
            self.void_code if include_external_void else None,
        )
        for marker_code in marker_codes:
            if marker_code is not None and marker_code not in codes:
                codes.append(marker_code)
        return codes

    def _code_for_state(self, state: str) -> int:
        """Encode a state or marker text (inverse of :meth:`_state_for_code`)."""
        if self.nr_code is not None and state == self.nr:
            return int(self.nr_code)
        if self.void_code is not None and state == self.void:
            return int(self.void_code)
        try:
            return self.alphabet.index(state)
        except ValueError:
            raise ValueError(f"state {state!r} is not in the codebook") from None

    def _canonical_code_ranks(self, *, void_as_missing: bool) -> np.ndarray:
        """Rank of every code as R orders cells: real states, then nr, then void.

        R sorts over the factor levels ``c(alphabet, nr, void)`` (seqdef.R:158),
        whatever integer codes traminex stores the markers under; sorti/sortv
        rank void as nr (TraMineRextras sortv.R:5), hence ``void_as_missing``.
        Real states keep their stored alphabet order (TX-AUD-20261011-006).
        The result is indexed by code and covers external marker codes too.
        """
        real = _real_state_codes(self)
        data = np.asarray(self.data)
        largest = max(
            [len(self.alphabet) - 1]
            + [int(code) for code in (self.nr_code, self.void_code) if code is not None]
            + ([int(data.max())] if data.size else [])
        )
        ranks = np.full(largest + 1, len(real) + 2, dtype=np.int64)
        ranks[real] = np.arange(len(real))
        if self.nr_code is not None:
            ranks[int(self.nr_code)] = len(real)
        if self.void_code is not None:
            ranks[int(self.void_code)] = len(real) + (0 if void_as_missing else 1)
        return ranks

    def _canonicalized(self) -> "STSList":
        """Return an equivalent object whose declared markers are alphabet members.

        The identity for objects already in the in-alphabet representation.
        For a declared external ``nr_code``/``void_code`` (a code outside
        ``range(len(alphabet))``) this returns a copy with the marker folded
        into the alphabet: when the marker text already names an alphabet
        member that no cell uses (the combined objects ``group_compare``
        builds from an in-alphabet source leave such a slot behind), the
        external code is merged into that existing slot; when a real state
        uses that text, the marker is appended under a collision-free spelling
        (``<NR>``, ``<VOID>``) and ``nr``/``void`` take it; otherwise the
        marker text is appended
        (nr first, then void, matching :meth:`_codebook_codes` order), with
        ``labels`` extended with the marker text and ``cpal`` extended with
        ``missing_color`` for nr and ``"white"`` for void.  The data codes
        are remapped accordingly, so consumers that canonicalize at entry
        treat both marker representations of the same data identically.
        """
        original_size = len(self.alphabet)
        external = [
            (attribute, marker, int(code))
            for attribute, marker, code in (
                ("nr_code", self.nr, self.nr_code),
                ("void_code", self.void, self.void_code),
            )
            if code is not None and int(code) >= original_size
        ]
        if not external:
            return self
        result = self.copy()
        # Collect the masks before any remapping so a swap between an old
        # external code and a freshly assigned one cannot corrupt the data.
        masks = {code: result.data == code for _, _, code in external}
        for attribute, marker, code in external:
            slot = result.alphabet.index(marker) if marker in result.alphabet else None
            if slot is not None and slot < original_size and not np.any(self.data == slot):
                # The marker text names an alphabet slot no cell uses (the
                # combined objects group_compare builds leave one behind):
                # merge into it rather than append a duplicate symbol.
                new_code = slot
            else:
                # A real state may already be spelled like the marker; then
                # the marker gets a collision-free spelling, as seqdef gives
                # it (TX-AUD-20261010-034).
                role = "NR" if attribute == "nr_code" else "VOID"
                text = marker
                if text is None or text in result.alphabet:
                    text = f"<{role}>"
                    suffix = 1
                    while text in result.alphabet:
                        text = f"<{role}:{suffix}>"
                        suffix += 1
                    setattr(result, "nr" if attribute == "nr_code" else "void", text)
                new_code = len(result.alphabet)
                result.alphabet.append(text)
                if result.labels is not None:
                    result.labels.append(text)
                if result.cpal is not None:
                    result.cpal.append(self.missing_color if attribute == "nr_code" else "white")
            result.data[masks[code]] = new_code
            setattr(result, attribute, new_code)
        return result


_P = ParamSpec("_P")
_R = TypeVar("_R")


def _require_stslist(function: Callable[_P, _R]) -> Callable[_P, _R]:
    """Validate the first supplied argument of an STSList public consumer."""
    first_parameter = next(iter(inspect.signature(function).parameters))

    @wraps(function)
    def checked(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        if args:
            value = args[0]
        elif first_parameter in kwargs:
            value = kwargs[first_parameter]
        else:
            return function(*args, **kwargs)
        if not isinstance(value, STSList):
            raise TypeError(f"{first_parameter} must be an STSList state sequence object")
        return function(*args, **kwargs)

    return cast(Callable[_P, _R], checked)


def _remap_issues(
    issues: dict[str, list[int]] | None,
    row_indices: np.ndarray,
) -> dict[str, list[int]] | None:
    """Remap 1-based issue indexes for a zero-based row subset."""
    if issues is None:
        return None
    selected = [int(index) + 1 for index in np.asarray(row_indices).reshape(-1)]
    remapped: dict[str, list[int]] = {}
    for category, indexes in issues.items():
        issue_indexes = set(indexes)
        remapped[category] = [
            new_index
            for new_index, original_index in enumerate(selected, start=1)
            if original_index in issue_indexes
        ]
    return remapped


# R-compat convenience accessors


def _real_state_codes(x: STSList) -> List[int]:
    """Internal codebook positions of the real states, in codebook order.

    ``seqdef`` currently stores observed missing/void markers alongside real
    states.  TraMineR keeps those markers as separate metadata, so the marker
    positions are excluded here.  Selection is by *encoded marker position*
    rather than by order or by marker text, so a loosely built object whose
    ordinary state reuses the marker text is still classified correctly.
    """
    marker_codes = {
        code
        for code in (getattr(x, "nr_code", None), getattr(x, "void_code", None))
        if code is not None and 0 <= code < len(x.alphabet)
    }
    return [code for code in range(len(x.alphabet)) if code not in marker_codes]


def _alphabet_codes(x: STSList, with_missing: bool = False) -> List[Optional[int]]:
    """Codebook positions addressed by :func:`alphabet`, one per returned state.

    ``None`` marks a state for which the internal codebook carries no row: R's
    "phantom" ``nr``, which ``with.missing=TRUE`` appends to ``statl`` whether or
    not any cell is actually missing (``SOURCES/TraMineR/R/seqistatd.R`` builds
    ``statl <- c(alphabet(seqdata), attr(seqdata, "nr"))`` unconditionally).
    Callers that index a per-state result by this list must supply the R value
    for that phantom row themselves -- it is a zero count in every R statistic
    that reaches it.

    Kept beside :func:`alphabet` and sharing its marker rule so the two cannot
    drift apart.
    """
    state_codes = _real_state_codes(x)
    codes: List[Optional[int]] = list(state_codes)
    if with_missing and x.nr not in [x.alphabet[code] for code in state_codes]:
        nr_code = x.nr_code
        codes.append(nr_code if nr_code is not None and 0 <= nr_code < len(x.alphabet) else None)
    return codes


def alphabet(x: object, with_missing: bool = False) -> List[str]:
    """Return the state or event alphabet (R: ``alphabet``)."""
    if isinstance(x, STSList):
        # Read through the canonical codebook, as cpal() and stlab() are, so both
        # marker representations give one alphabet (TX-AUD-20261019-042).
        x = x._canonicalized()
        states = [x.alphabet[code] for code in _real_state_codes(x)]
        if with_missing and x.nr not in states:
            states.append(x.nr)
        return states

    # Import locally because the event-object module imports STSList.
    from .events._object import SeqElist

    if isinstance(x, SeqElist):
        return list(x.dictionary)
    raise TypeError("x must be a state sequence object or an event sequence object")


@_require_stslist
def cpal(x: STSList) -> Optional[List[str]]:
    """Return the color palette associated with the sequence object (R: cpal).

    One color per codebook entry, markers included, for both marker
    representations: declared external markers are read through the
    canonical codebook (audit S2-14).
    """
    x = x._canonicalized()
    return None if x.cpal is None else list(x.cpal)


@_require_stslist
def stlab(x: STSList) -> List[str]:
    """Return state labels (R: stlab). If absent, falls back to alphabet.

    One label per codebook entry, markers included, for both marker
    representations: declared external markers are read through the
    canonical codebook (audit S2-14).
    """
    x = x._canonicalized()
    return list(x.labels) if x.labels is not None else list(x.alphabet)


def _seqgen_probability(name: str, value: float) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    resolved = float(value)
    if not np.isfinite(resolved) or not 0 <= resolved <= 1:
        raise ValueError(f"{name} must be finite and between 0 and 1")
    return resolved


def _seqgen_missing_code(x: STSList, missing_type: Literal["nr", "void"]) -> int:
    attribute = f"{missing_type}_code"
    code = getattr(x, attribute)
    if code is not None:
        return int(code)

    # A real state may already be spelled like the marker (a literal "*" kept
    # by the SPS/SPELL/TSE routes): the new marker then takes a spelling no
    # state uses, as seqdef gives it, instead of a second "*" (audit S2-09).
    from .seqformat import _collision_free_marker

    default = "*" if missing_type == "nr" else "%"
    current = getattr(x, missing_type)
    marker = _collision_free_marker(
        default if current is None else str(current),
        x.alphabet,
        missing_type.upper(),
    )
    setattr(x, missing_type, marker)
    code = len(x.alphabet)
    x.alphabet.append(marker)
    if x.labels is not None:
        x.labels.append(marker)
    if x.cpal is not None:
        x.cpal.append(x.missing_color if missing_type == "nr" else "white")
    setattr(x, attribute, code)
    return code


def seqgen_missing(
    x: STSList,
    *,
    p_cases: float = 0.1,
    p_left: float = 0.2,
    p_gaps: float = 0.0,
    p_right: float = 0.3,
    mt_left: Literal["nr", "void"] = "nr",
    mt_gaps: Literal["nr", "void"] = "nr",
    mt_right: Literal["nr", "void"] = "nr",
    random_state: int | np.random.Generator | None = None,
) -> STSList:
    """Insert random nonresponse/void cells using extras' deterministic mechanics.

    Case and position counts use R's round-to-even rule. Selected rows are
    sorted before processing; gaps are written first, then left missing cells,
    then right missing cells, so later regions overwrite earlier overlaps.
    NumPy's generator is intentionally local: a seed is reproducible within
    Python but does not claim byte-for-byte compatibility with R's global RNG.
    A marker the object does not yet use is added after the states; if a real
    state is already spelled like it (a literal ``"*"``), the marker takes a
    spelling no state uses (``"<NR>"``/``"<VOID>"``), as ``seqdef`` does.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList state sequence object")

    probabilities = {
        "p_cases": _seqgen_probability("p_cases", p_cases),
        "p_left": _seqgen_probability("p_left", p_left),
        "p_gaps": _seqgen_probability("p_gaps", p_gaps),
        "p_right": _seqgen_probability("p_right", p_right),
    }
    missing_types = {
        "mt_left": mt_left,
        "mt_gaps": mt_gaps,
        "mt_right": mt_right,
    }
    for name, missing_type in missing_types.items():
        if missing_type not in {"nr", "void"}:
            raise ValueError(f'{name} must be either "nr" or "void"')

    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-09) — with an external
    # void_code == len(alphabet), the appended nr would take the same code.
    x = x._canonicalized()
    result = x.copy()
    if x.void_code is None:
        lengths = np.full(x.nseq, x.seqlen, dtype=int)
    else:
        lengths = np.sum(x.data != int(x.void_code), axis=1, dtype=int)
    max_length = int(np.max(lengths, initial=0))
    selected_count = int(np.rint(probabilities["p_cases"] * x.nseq))

    rng = (
        random_state
        if isinstance(random_state, np.random.Generator)
        else np.random.default_rng(random_state)
    )
    selected = np.sort(rng.choice(x.nseq, size=selected_count, replace=False))
    left_draws = rng.uniform(0, probabilities["p_left"], size=x.nseq)
    gap_draws = rng.uniform(0, probabilities["p_gaps"], size=x.nseq)
    right_draws = rng.uniform(0, probabilities["p_right"], size=x.nseq)

    codes: dict[str, int] = {}
    active_types = (
        (mt_left, probabilities["p_left"]),
        (mt_gaps, probabilities["p_gaps"]),
        (mt_right, probabilities["p_right"]),
    )
    if selected_count and max_length:
        needed_types = {missing_type for missing_type, probability in active_types if probability}
        for missing_type in ("nr", "void"):
            if missing_type in needed_types:
                codes[missing_type] = _seqgen_missing_code(result, missing_type)

    for row in selected:
        gap_count = int(np.rint(gap_draws[row] * max_length))
        gaps = rng.choice(max_length, size=gap_count, replace=False)
        if gap_count:
            result.data[row, gaps] = codes[mt_gaps]

        left_count = int(np.rint(left_draws[row] * max_length))
        if left_count:
            result.data[row, :left_count] = codes[mt_left]

        right_count = int(np.rint(right_draws[row] * max_length))
        if right_count:
            result.data[row, max_length - right_count : max_length] = codes[mt_right]
    return result


@dataclass
class SeqMD:
    """Lightweight multi-domain container (channels + optional names).

    Holds the per-domain :class:`STSList` channels for the additive,
    non-R-faithful weight-search utilities (:func:`seqMD_distances`,
    ``seqMD_gridsearch_weights``, ``seqMD_optimize_weights``). The faithful
    TraMineR CAT engine lives in :func:`traminex.seqMD` /
    :func:`traminex.seqdistmc` (module ``multidomain_cat``).
    """

    channels: List[STSList]
    names: Optional[List[str]] = None

    def __post_init__(self) -> None:
        self.channels, _ = _validate_seqmd_channels(self.channels)
        if self.names is not None:
            self.names = list(self.names)
            if len(self.names) != len(self.channels):
                raise ValueError(
                    f"names must give one name per channel ({len(self.channels)} channels, "
                    f"{len(self.names)} names)"
                )


def _validate_seqmd_channels(channels: object) -> tuple[List[STSList], int]:
    """Return a nonempty list of aligned state-sequence channels."""
    if not isinstance(channels, Iterable):
        raise TypeError("channels must be a collection of STSList objects") from None
    resolved = list(channels)
    if not resolved:
        raise ValueError("channels must contain at least one channel")
    if any(not isinstance(channel, STSList) for channel in resolved):
        raise TypeError("channels must contain only STSList state sequence objects")
    n_sequences = resolved[0].nseq
    if any(channel.nseq != n_sequences for channel in resolved[1:]):
        raise ValueError("channels must contain the same number of sequences")
    return resolved, n_sequences


def seqMD_distances(
    obj: SeqMD,
    *,
    method: str = "LCS",
    weights: Optional[List[float]] = None,
    per_channel_kwargs: Optional[List[dict]] = None,
) -> np.ndarray:
    """Additive combination of per-channel distances.

    Returns an n×n matrix D = sum_c w_c * D_c where D_c = seqdist(channel_c, method, **kwargs_c).

    Raises
    ------
    ValueError
        If ``obj`` carries no channels. There is no n×n matrix to return for an
        empty sum, and the previous ``None`` violated the declared return type
        and surfaced as an unrelated ``AttributeError`` at the caller.
    """
    if not isinstance(obj, SeqMD):
        raise TypeError("obj must be a SeqMD multichannel object (seqMD)")
    # Import here to avoid circular dependency at module load time
    from .distances import seqdist

    channels, n_sequences = _validate_seqmd_channels(obj.channels)
    C = len(channels)
    if weights is None:
        w = [1.0] * C
    else:
        w = _validate_analytical_weights(weights, name="weights", n_channels=C)
    if per_channel_kwargs is None:
        per_channel_kwargs = [{} for _ in range(C)]
    elif len(per_channel_kwargs) != C:
        raise ValueError("per_channel_kwargs length must equal number of channels")
    Dsum: np.ndarray | None = None
    for c in range(C):
        D = np.asarray(seqdist(channels[c], method=method, **per_channel_kwargs[c]), dtype=float)
        if D.shape != (n_sequences, n_sequences):
            raise ValueError("channel distances must have the same square shape")
        D = w[c] * D
        Dsum = D if Dsum is None else (Dsum + D)
    assert Dsum is not None  # C >= 1 was established above.
    return Dsum
