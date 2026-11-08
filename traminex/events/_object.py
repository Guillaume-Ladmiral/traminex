"""Faithful (R-parity) event-sequence object and accessors -- RP5a.

The TraMineR ``eseq``/``seqelist`` object (:class:`Eseq`/:class:`SeqElist`,
:func:`seqecreate`), its accessors (:func:`seqeid`, :func:`seqeweight`,
:func:`seqelength`, :func:`seqecontain`), the type guards
(:func:`is_eseq`/:func:`is_seqelist`/:func:`is_stslist`) and the BUILD-from-R
transition machinery (:func:`seqetm`, :func:`seqe2stm`). Every numeric item here
is pinned to live TraMineR 2.2-14 / TraMineRextras 0.6.9 goldens.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from numbers import Number
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .._categorical_order import stable_sorted_unique
from .._deprecated import MISSING, resolve_deprecated_alias
from .._rformat import r_as_character
from ..sequences import STSList, _require_stslist
from . import _engine
from ._helpers import _ReadOnlyList


@_require_stslist
def seqetm(
    x=MISSING,
    method: str = "transition",
    *,
    use_labels: bool = True,
    sep: str = ">",
    bp: str = "",
    ep: str = "end",
    seq=MISSING,
) -> pd.DataFrame:
    """Transition / event matrix from a state-sequence alphabet (R: ``seqetm``).

    Faithful port of ``TraMineR::seqetm`` (``seqetm.R``): build the
    ``nbstat x nbstat`` character matrix whose rows/columns are the state
    alphabet and whose cells label the event associated with each (from-state,
    to-state) pair.

    - the diagonal cell ``(i, i)`` is the state itself -- its long label when
      ``use_labels`` is set and labels are present, otherwise the alphabet code;
    - ``method="transition"`` -> ``"<from><sep><to>"`` (default ``sep=">"``);
    - ``method="period"``     -> ``"<ep><from>,<bp><to>"`` (defaults ``ep="end"``,
      ``bp=""``, e.g. ``"endA,B"``);
    - ``method="state"``      -> ``"<to>"``.

    Missing (``nr``) and void states actually present in the data are appended
    to the alphabet, exactly as R does. The returned matrix is suitable as the
    ``tevent`` argument of a state->event conversion.

    Parameters
    ----------
    x : STSList
        State-sequence object (R ``stslist``).
    method : {"transition", "period", "state"}, default "transition"
        Event-encoding rule for off-diagonal cells.
    use_labels : bool, default True
        Use the long state labels (``x.labels``) instead of the alphabet codes.
    sep : str, default ">"
        Separator for ``method="transition"``.
    bp, ep : str, defaults "" / "end"
        Begin/end prefixes for ``method="period"``.

    Returns
    -------
    pandas.DataFrame
        Square character matrix indexed by the (code) alphabet on both axes.
    """
    x = resolve_deprecated_alias(
        x,
        seq,
        canonical_name="x",
        alias_name="seq",
        function_name="seqetm",
    )
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList state sequence object")
    if method not in ("transition", "period", "state"):
        raise ValueError("method must be one of 'transition', 'period', 'state'")
    # Both marker representations give one matrix (TX-REP-20261023-003).
    x = x._canonicalized()

    # R's seqetm starts from alphabet(seqdata) = the REAL states only (it holds nr
    # and void as separate attributes), then appends each once if present. traminex's
    # x.alphabet already folds nr/void in, so start from the real-state codes to avoid
    # appending them a second time (which produced an N+2 matrix with duplicate cols).
    real_codes = [c for c in range(len(x.alphabet)) if c != x.nr_code and c != x.void_code]
    statl = [x.alphabet[c] for c in real_codes]
    has_nr = x.nr_code is not None and bool(np.any(x.data == x.nr_code))
    if has_nr:
        statl = statl + [x.nr]
    has_void = x.void_code is not None and bool(np.any(x.data == x.void_code))
    if has_void:
        statl = statl + [x.void]
    nbstat = len(statl)

    alpha = list(statl)
    if use_labels and x.labels is not None:
        label = [x.labels[code] for code in real_codes]
        if has_nr:
            label = label + [x.nr]
        if has_void:
            label = label + [x.void]
        if len(label) == nbstat:
            alpha = label
        elif len(label) > 0:
            warnings.warn("Length of the labels and of the alphabet are not equal")

    if any("," in a for a in alpha):
        warnings.warn(
            " [!] Alphabet and/or state labels should not contain commas ',' which are "
            "reserved for separating multiple events of a same transition!"
        )

    cells = [["" for _ in range(nbstat)] for _ in range(nbstat)]
    for i in range(nbstat):
        for j in range(nbstat):
            if i == j:
                cells[i][j] = alpha[i]
            elif method == "transition":
                cells[i][j] = f"{alpha[i]}{sep}{alpha[j]}"
            elif method == "period":
                cells[i][j] = f"{ep}{alpha[i]},{bp}{alpha[j]}"
            else:  # state
                cells[i][j] = alpha[j]
    return pd.DataFrame(cells, index=list(statl), columns=list(statl))


def seqe2stm(
    events: List[str],
    *,
    drop_matrix: Optional[pd.DataFrame] = None,
    drop_list: Optional[Dict[str, object]] = None,
    first_state: str = "None",
) -> pd.DataFrame:
    """Event -> state-machine transition table (R: ``seqe2stm``).

    Faithful port of ``TraMineRextras::seqe2stm`` (``seqe2stm.R``): build the
    deterministic automaton whose *states* are the reachable sets of events that
    have occurred (rendered as dot-joined event names in ``events`` order, or
    ``first_state`` for the empty set) and whose *columns* are the events. Cell
    ``[s, e]`` is the state reached from state ``s`` when event ``e`` occurs.

    An optional drop rule lets an event clear previously-occurred events:

    - ``drop_list`` maps an event to the event(s) it clears, e.g.
      ``{"B": "A"}`` (occurrence of ``B`` removes ``A`` from the set);
    - ``drop_matrix`` is the equivalent boolean ``events x events`` table
      (``True`` at ``[a, b]`` means occurrence of ``a`` clears ``b``).

    Only the reachable states are returned, ordered by first appearance starting
    from ``first_state`` (matching R's ``unique(c(firstState, stmat))`` row pick).

    Parameters
    ----------
    events : list of str
        Event alphabet.
    drop_matrix : pandas.DataFrame, optional
        Boolean ``events x events`` drop table (labelled rows/cols).
    drop_list : dict, optional
        ``{event: event | list_of_events}`` drop rules.
    first_state : str, default "None"
        Name of the empty (initial) state.

    Returns
    -------
    pandas.DataFrame
        Character matrix indexed by reachable state names, columns = events.
    """
    events = [str(e) for e in events]
    nevent = len(events)
    ev_idx = {e: i for i, e in enumerate(events)}

    dropM: Optional[np.ndarray] = None
    if drop_list is not None:
        dropM = np.zeros((nevent, nevent), dtype=bool)
        for k, vals in drop_list.items():
            seq = [vals] if isinstance(vals, (str, bytes)) else list(vals)  # type: ignore[call-overload]
            for v in seq:
                dropM[ev_idx[str(k)], ev_idx[str(v)]] = True
    elif drop_matrix is not None:
        dm = pd.DataFrame(drop_matrix)
        dropM = np.zeros((nevent, nevent), dtype=bool)
        for a in dm.index:
            for b in dm.columns:
                if bool(dm.loc[a, b]):
                    dropM[ev_idx[str(a)], ev_idx[str(b)]] = True

    def state_names(occ: List[bool], add: Optional[int] = None) -> str:
        occ = list(occ)
        if add is not None:
            occ[add] = True
            if dropM is not None:
                for b in range(nevent):
                    if dropM[add, b]:
                        occ[b] = False
        if not any(occ):
            return first_state
        return ".".join(events[i] for i in range(nevent) if occ[i])

    nbstate = 1 << nevent
    stmatL: List[List[bool]] = []

    def rec(i: int, occured: List[bool]) -> None:
        if i == nevent:
            stmatL.append(list(occured))
            return
        occured[i] = False
        rec(i + 1, occured)
        occured[i] = True
        rec(i + 1, occured)

    rec(0, [False] * nevent)

    stlist = [state_names(row) for row in stmatL]
    stmat = [[state_names(row, add=j) for j in range(nevent)] for row in stmatL]

    # Reachable states in R order: unique(c(firstState, stmat[column-major])).
    flat = [first_state]
    for j in range(nevent):
        for i in range(nbstate):
            flat.append(stmat[i][j])
    seen: List[str] = []
    seen_set: set[str] = set()
    for nm in flat:
        if nm not in seen_set:
            seen.append(nm)
            seen_set.add(nm)
    name_to_row: Dict[str, int] = {}
    for i, nm in enumerate(stlist):
        if nm not in name_to_row:
            name_to_row[nm] = i
    rows_out = [nm for nm in seen if nm in name_to_row]
    data = [stmat[name_to_row[nm]] for nm in rows_out]
    return pd.DataFrame(data, index=rows_out, columns=events)


@dataclass
class Eseq:
    """A single event sequence (TraMineR ``eseq``).

    Events are stored sorted by (time, dictionary code); ``codes`` are 1-based
    indices into ``dictionary``. ``weight`` defaults to 1.0 and ``obs_time``
    (the stored sequence length, R ``seqelength``) defaults to -1.0 for
    TSE-created sequences.

    TraMineR keeps an event sequence as a chain of gaps built by inserting the
    events in their input order; with fractional times that order decides how
    the gaps round. :func:`seqecreate` records it so that the printed gaps,
    event ages and subsequence counts reproduce R's chain; an ``Eseq`` built
    directly is read in its array order.
    """

    id_code: int
    times: np.ndarray
    codes: np.ndarray
    dictionary: List[str]
    weight: float = 1.0
    obs_time: float = -1.0
    # Positions into ``times``/``codes`` in the order the events were given to
    # seqecreate (R inserts them in that order); None means array order.
    _insertion_order: Optional[np.ndarray] = field(
        default=None, init=False, repr=False, compare=False
    )

    def __eq__(self, other: object) -> bool:
        # Value equality over the compared fields. The generated __eq__ compared the arrays inside a
        # tuple and raised for more than one event (TX-AUD-20261011-059); R defines no equality
        # (seqe.R:41). Defining __eq__ leaves __hash__ None: an Eseq is mutable. With fractional
        # times the input order decides how the gaps round, so the gap chain the engine reads is
        # compared too (batch-11 review M-6).
        if other.__class__ is not self.__class__:
            return NotImplemented
        return (
            self.id_code == other.id_code
            and np.array_equal(self.times, other.times)
            and np.array_equal(self.codes, other.codes)
            and self.dictionary == other.dictionary
            and self.weight == other.weight
            and self.obs_time == other.obs_time
            and self._gap_chain() == other._gap_chain()
        )

    def _gap_chain(self) -> List[tuple]:
        """The (gap, event code) pairs of the engine's node chain, in chain order."""
        chain = []
        node = _engine.eseq_to_engine(self).event
        while node is not None:
            chain.append((node.gap, node.type))
            node = node.next
        return chain

    def as_character(self) -> str:
        """The TraMineR event-label string (R ``as.character.eseq``)."""
        return _engine.eseq_to_engine(self).sprint()

    def __str__(self) -> str:
        return self.as_character()

    def __repr__(self) -> str:
        return f"Eseq({self.as_character()!r})"


def _copy_eseq(source: Eseq) -> Eseq:
    """An independent copy of ``source``, including its recorded insertion order."""
    duplicate = Eseq(
        id_code=source.id_code,
        times=source.times.copy(),
        codes=source.codes.copy(),
        dictionary=list(source.dictionary),
        weight=source.weight,
        obs_time=source.obs_time,
    )
    if source._insertion_order is not None:
        duplicate._insertion_order = source._insertion_order.copy()
    return duplicate


class SeqElist:
    """A list of event sequences (TraMineR ``seqelist``), returned by :func:`seqecreate`.

    Iterable / indexable over its :class:`Eseq` elements: an integer gives one
    :class:`Eseq`; a slice, a list or array of positions, or a boolean mask gives
    a ``SeqElist`` of copies, with NumPy's 0-based rules (a mask of another
    length raises ``IndexError``). ``dictionary`` is the
    event alphabet (R ``levels``). For backward compatibility with the non-R
    (approximate) event miner and helpers in this module, the flat per-event
    views ``id`` / ``time`` / ``label`` are preserved in input order. They are
    read-only: they must describe the same events as the sequences, so an
    in-place write raises instead of changing one side only.
    """

    def __init__(
        self,
        sequences: List["Eseq"],
        dictionary: List[str],
        *,
        flat_id: List[str],
        flat_time: np.ndarray,
        flat_label: List[str],
    ) -> None:
        self.sequences = sequences
        self.dictionary = list(dictionary)
        # The flat views must describe the same events as the sequences (cycle-8
        # F-16) and have no setter, so they are stored read-only: an in-place
        # write would change what the flat-view helpers see but not what the
        # engine reads (2026-09-22 audit S5-05). Copies, so the caller's
        # vectors stay writable.
        self._flat_id = _ReadOnlyList(flat_id)
        self._flat_time = np.array(flat_time)
        self._flat_time.flags.writeable = False
        self._flat_label = _ReadOnlyList(flat_label)
        # Flat-row positions of each sequence's events, set by the constructors
        # in this package; None falls back to matching runs of equal id text.
        self._flat_groups: Optional[List[List[int]]] = None

    def _flat_rows_by_sequence(self) -> List[List[int]]:
        """The flat-row positions that belong to each sequence, in sequence order."""
        if self._flat_groups is not None and len(self._flat_groups) == len(self.sequences):
            return self._flat_groups
        runs: List[List[int]] = []
        for position, identifier in enumerate(self._flat_id):
            if not runs or self._flat_id[runs[-1][0]] != identifier:
                runs.append([])
            runs[-1].append(position)
        if not runs:
            return [[] for _ in self.sequences]
        if len(runs) != len(self.sequences):
            # An id holding only its end event has no flat row, and ids that
            # R's as.integer merges (1.2 and 1.7) make two runs: guessing would
            # hand one sequence another's events (2026-09-22 audit S5-03).
            raise ValueError(
                f"the flat id/time/label views hold {len(runs)} runs of equal ids for "
                f"{len(self.sequences)} sequences, so they cannot be matched to the "
                "sequences to slice them; build the list with seqecreate"
            )
        return runs

    @property
    def id(self) -> List[str]:
        return self._flat_id

    @property
    def time(self) -> np.ndarray:
        return self._flat_time

    @property
    def label(self) -> List[str]:
        return self._flat_label

    @property
    def weights(self) -> np.ndarray:
        """Per-sequence weights (R ``seqeweight``); assignable (R ``seqeweight<-``)."""
        return np.array([s.weight for s in self.sequences], dtype=float)

    @weights.setter
    def weights(self, value) -> None:
        # R compares length(value): a single number is a length-1 vector (S5-10).
        arr = np.atleast_1d(np.asarray(value, dtype=float))
        if len(arr) != len(self.sequences):
            raise ValueError("eseq and weights should be of the same size.")
        for s, w in zip(self.sequences, arr):
            s.weight = float(w)

    @property
    def lengths(self) -> np.ndarray:
        """Per-sequence stored lengths (R ``seqelength``); assignable (R ``seqelength<-``)."""
        return np.array([s.obs_time for s in self.sequences], dtype=float)

    @lengths.setter
    def lengths(self, value) -> None:
        # R compares length(value): a single number is a length-1 vector (S5-10).
        arr = np.atleast_1d(np.asarray(value, dtype=float))
        if len(arr) != len(self.sequences):
            raise ValueError("eseq and len should be of the same size.")
        for s, v in zip(self.sequences, arr):
            s.obs_time = float(v)

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            selected_indices = list(range(len(self.sequences)))[idx]
        elif np.ndim(idx) == 0:
            return self.sequences[idx]
        else:
            # R's "[.seqelist" takes any index vector (seqe.R:22-32); lists, integer arrays and
            # boolean masks select with NumPy's 0-based rules (TX-AUD-20261011-058).
            positions = np.arange(len(self.sequences))[idx]
            if positions.ndim != 1:
                raise IndexError("a SeqElist index must select along one axis")
            selected_indices = positions.tolist()
        selected_sequences = [_copy_eseq(self.sequences[index]) for index in selected_indices]

        rows_by_sequence = self._flat_rows_by_sequence()
        selected_positions: List[int] = []
        selected_groups: List[List[int]] = []
        for index in selected_indices:
            rows = rows_by_sequence[index]
            start = len(selected_positions)
            selected_groups.append(list(range(start, start + len(rows))))
            selected_positions.extend(rows)
        result = SeqElist(
            selected_sequences,
            self.dictionary,
            flat_id=[self._flat_id[position] for position in selected_positions],
            flat_time=np.asarray(self._flat_time)[selected_positions].copy(),
            flat_label=[self._flat_label[position] for position in selected_positions],
        )
        result._flat_groups = selected_groups
        return result

    def __iter__(self):
        return iter(self.sequences)

    def as_character(self) -> List[str]:
        """The TraMineR event-label strings, one per sequence (R ``as.character``)."""
        return [s.as_character() for s in self.sequences]

    def __str__(self) -> str:
        return "\n".join(self.as_character())

    def __repr__(self) -> str:
        return f"SeqElist(n={len(self.sequences)}, events={len(self.dictionary)})"


def _event_label(value: object) -> str:
    """An event value as R's ``factor(event)`` level (``as.character``: 1e+05, TRUE)."""
    text = r_as_character(value)
    return "NA" if text is None else text


def seqecreate(
    data: pd.DataFrame | STSList | None = None,
    *,
    id: Sequence[object] | None = None,
    timestamp: Sequence[float] | None = None,
    event: Sequence[object] | None = None,
    end_event: str | float | None = None,
    tevent: str | np.ndarray = "transition",
    use_labels: bool = True,
    weighted: bool = True,
    id_col: str = "id",
    time_col: str = "timestamp",
    event_col: str = "event",
    endEvent=MISSING,
) -> "SeqElist":
    """Create an event-sequence object from TSE vectors/frames or an STSList.

    This completes TraMineR's paired-vector and state-to-event routes while
    retaining the established custom-column DataFrame call. STS conversion is
    delegated to :func:`seqetm` and :func:`traminex.seqformat.seqformat`, so the
    package has one state-to-event implementation.

    Numeric and logical event values are labelled as R's ``factor(event)``
    labels them (``as.character``: ``1e+05``, ``TRUE``), and ``end_event`` is
    matched after the same rendering, so ``end_event=99`` ends observation at
    the event ``99``. An ``end_event`` that matches no event is ignored, as in
    R.
    """
    end_event = resolve_deprecated_alias(
        end_event,
        endEvent,
        canonical_name="end_event",
        alias_name="endEvent",
        function_name="seqecreate",
        default=None,
    )
    state_source: STSList | None = None
    if isinstance(data, STSList):
        from ..seqformat import seqformat

        state_source = data
        if isinstance(tevent, str):
            transition_matrix = seqetm(data, method=tevent, use_labels=use_labels)
        elif isinstance(tevent, pd.DataFrame):
            transition_matrix = tevent.copy(deep=True)
        else:
            values = np.asarray(tevent, dtype=object)
            template = seqetm(data, method="transition", use_labels=False)
            if values.ndim != 2 or values.shape != template.shape:
                raise ValueError("tevent matrix must match the state transition dimensions")
            transition_matrix = pd.DataFrame(
                values.copy(), index=template.index.copy(), columns=template.columns.copy()
            )
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="'missing' set as .*state sequence object",
                category=UserWarning,
            )
            converted = seqformat(
                data,
                from_="STS",
                to="TSE",
                tevent=transition_matrix,
                id=list(range(1, data.nseq + 1)),
            )
        assert isinstance(converted, pd.DataFrame)
        data = converted
    elif data is not None and not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a DataFrame, STSList, or None")

    if isinstance(data, pd.DataFrame):
        resolved_time_col = time_col
        if time_col == "timestamp" and time_col not in data.columns and "time" in data.columns:
            resolved_time_col = "time"
        if {id_col, resolved_time_col, event_col}.issubset(data.columns):
            id = data[id_col].copy(deep=True)
            timestamp = data[resolved_time_col].copy(deep=True)
            event = data[event_col].copy(deep=True)

    if id is None:
        raise ValueError(" [!] Could not find an id argument")
    if timestamp is None:
        raise ValueError(" [!] Could not find a timestamp argument")
    if event is None:
        raise ValueError(" [!] Could not find an event argument")

    id_series = id.copy(deep=True) if isinstance(id, pd.Series) else pd.Series(id)
    time_series = (
        timestamp.copy(deep=True) if isinstance(timestamp, pd.Series) else pd.Series(timestamp)
    )
    if isinstance(event, pd.Series):
        event_series = event.copy(deep=True)
    elif isinstance(event, pd.Categorical):
        event_series = pd.Series(event.copy())
    else:
        event_series = pd.Series(event)
    if not (len(id_series) == len(time_series) == len(event_series)):
        raise ValueError("Time ,idpers and event vector should have the same size")
    if id_series.isna().any() or time_series.isna().any() or event_series.isna().any():
        raise ValueError(" [!] Missing values not supported")

    try:
        times_num = pd.to_numeric(time_series, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError):
        raise ValueError("timestamp values must be numeric") from None

    if isinstance(event_series.dtype, pd.CategoricalDtype):
        event_factor = event_series.array
        dictionary = [_event_label(value) for value in event_factor.categories]
        event_codes = np.asarray(event_factor.codes, dtype=int) + 1
    elif all(isinstance(value, str) for value in event_series):
        raw_labels = event_series.astype(str).tolist()
        dictionary = [str(value) for value in stable_sorted_unique(raw_labels)]
        code_of = {label: code for code, label in enumerate(dictionary, start=1)}
        event_codes = np.asarray([code_of[label] for label in raw_labels], dtype=int)
    else:
        event_factor = pd.Categorical(
            event_series,
            categories=stable_sorted_unique(event_series.tolist()),
        )
        dictionary = [_event_label(value) for value in event_factor.categories]
        event_codes = np.asarray(event_factor.codes, dtype=int) + 1
    labels = [dictionary[code - 1] for code in event_codes]
    if any(character in label for label in dictionary for character in "(),"):
        warnings.warn(
            " [!] some of your events contain '(', ')' or ',' characters. "
            "The search of specific subsequences may not work properly.",
            UserWarning,
            stacklevel=2,
        )

    if isinstance(id_series.dtype, pd.CategoricalDtype):
        id_codes = id_series.cat.codes.to_numpy(dtype=int) + 1
    elif pd.api.types.is_numeric_dtype(id_series.dtype) or pd.api.types.is_bool_dtype(
        id_series.dtype
    ):
        id_codes = id_series.to_numpy(dtype=float).astype(int)
    else:
        id_codes = (
            np.asarray(
                pd.Categorical(
                    id_series,
                    categories=stable_sorted_unique(id_series.tolist()),
                ).codes,
                dtype=int,
            )
            + 1
        )
    flat_ids = [str(value) for value in id_series]

    groups: list[list[int]] = []
    seen_ids: set[int] = set()
    previous_id: int | None = None
    for position, id_code in enumerate(id_codes):
        numeric_id = int(id_code)
        if numeric_id != previous_id:
            if numeric_id in seen_ids:
                raise ValueError(" [!] Events not grouped by id! See seqecreate help page.")
            seen_ids.add(numeric_id)
            groups.append([])
            previous_id = numeric_id
        groups[-1].append(position)

    # R compares ``levels(event) == end.event``, which renders a number or a
    # logical with as.character; an end event matching no level is ignored
    # (seqecreate.R:63-73, where the not-found stop is unreachable).
    end_code = None
    end_label = None if end_event is None else _event_label(end_event)
    if end_label is not None and end_label in dictionary:
        end_code = dictionary.index(end_label) + 1
    sequences: List[Eseq] = []
    for positions in groups:
        obs_time = -1.0
        items: list[tuple[float, int]] = []
        for position in positions:
            item = (float(times_num[position]), int(event_codes[position]))
            if end_code is not None and item[1] == end_code:
                obs_time = item[0]
            else:
                items.append(item)
        # Stored sorted by (time, code); R builds its gap chain by inserting the
        # events in input order, which decides how fractional gaps round, so
        # that order is kept alongside (see _engine.eseq_to_engine).
        ranked = sorted(range(len(items)), key=lambda rank: items[rank])
        insertion_order = np.empty(len(items), dtype=int)
        insertion_order[np.asarray(ranked, dtype=int)] = np.arange(len(items))
        sequence = Eseq(
            id_code=int(id_codes[positions[0]]),
            times=np.asarray([items[rank][0] for rank in ranked], dtype=float),
            codes=np.asarray([items[rank][1] for rank in ranked], dtype=int),
            dictionary=dictionary,
            obs_time=obs_time,
        )
        sequence._insertion_order = insertion_order
        sequences.append(sequence)

    if end_code is not None:
        # The flat id/time/label views must describe the same events as the
        # per-sequence objects: end_event rows terminate observation and are
        # excluded there, so counting helpers reading the flat views would
        # otherwise silently include them (cycle-8 audit F-16).
        kept = [
            position for position in range(len(labels)) if int(event_codes[position]) != end_code
        ]
        flat_ids = [flat_ids[position] for position in kept]
        flat_times = times_num[kept].copy()
        flat_labels = [labels[position] for position in kept]
    else:
        kept = list(range(len(labels)))
        flat_times = times_num.copy()
        flat_labels = labels
    flat_row = {position: row for row, position in enumerate(kept)}
    result = SeqElist(
        sequences,
        dictionary,
        flat_id=flat_ids,
        flat_time=flat_times,
        flat_label=flat_labels,
    )
    result._flat_groups = [
        [flat_row[position] for position in positions if position in flat_row]
        for positions in groups
    ]
    if state_source is not None:
        from ..stats_basic import seqlength

        result.lengths = np.asarray(seqlength(state_source), dtype=float)
        if weighted and state_source.weights is not None:
            result.weights = np.asarray(state_source.weights, dtype=float)
    return result


def seqeid(eseq=MISSING, *, s=MISSING):
    """Per-sequence integer id (R ``seqeid``) -- the id factor code.

    Returns a 1-D int array for a :class:`SeqElist`, or a scalar for a single
    :class:`Eseq`.
    """
    eseq = resolve_deprecated_alias(
        eseq,
        s,
        canonical_name="eseq",
        alias_name="s",
        function_name="seqeid",
    )
    if isinstance(eseq, SeqElist):
        return np.array([s.id_code for s in eseq.sequences], dtype=int)
    if isinstance(eseq, Eseq):
        return int(eseq.id_code)
    raise TypeError("eseq should be a SeqElist or Eseq (see seqecreate).")


def seqeweight(eseq=MISSING, *, s=MISSING):
    """Per-sequence weight (R ``seqeweight``); defaults to 1.0.

    Returns a 1-D float array for a :class:`SeqElist`, or a scalar for a single
    :class:`Eseq`. The replacement form ``seqeweight(eseq) <- value`` is the
    assignable :attr:`SeqElist.weights` property (``eseq.weights = value``).
    """
    eseq = resolve_deprecated_alias(
        eseq,
        s,
        canonical_name="eseq",
        alias_name="s",
        function_name="seqeweight",
    )
    if isinstance(eseq, SeqElist):
        return np.array([s.weight for s in eseq.sequences], dtype=float)
    if isinstance(eseq, Eseq):
        return float(eseq.weight)
    raise TypeError("eseq should be a SeqElist or Eseq (see seqecreate).")


def seqecontain(
    eseq=MISSING,
    event_list=MISSING,
    *,
    unknown_exclude: bool = False,
    seq=MISSING,
    eventList=MISSING,
    exclude=MISSING,
):
    """Whether each sequence contains the given events (R ``seqecontain``).

    Port of ``TraMineR::seqecontain`` (``seqecontain.R`` /
    ``Sequence::contain``): ``event_list`` is a flat list of event labels.

    - ``unknown_exclude=False`` (default): a sequence matches if it contains
      **any** event in ``event_list``;
    - ``unknown_exclude=True``: a sequence matches only if **all** of its events
      are in ``event_list``.

    Labels in ``event_list`` that are not in the event dictionary are ignored
    (as in R's ``factor(event.list, levels=dict)``); each value is read as R's
    ``as.character`` reads it, so ``1.0`` names the event ``1``, and a single
    value may be given without a list. Returns a 1-D bool array for
    a :class:`SeqElist`, or a scalar for a single :class:`Eseq`, which R refuses
    ("eseq should be a seqelist"; a recorded divergence, TX-AUD-20261011-005). A
    :class:`~traminex.events.SubSeqElist` (a :func:`seqefsub` result) is
    accepted as in R (``if (is.subseqelist(eseq)) eseq <- eseq$subseq``): each
    of its subsequences is tested, giving one value per subsequence.
    """
    eseq = resolve_deprecated_alias(
        eseq,
        seq,
        canonical_name="eseq",
        alias_name="seq",
        function_name="seqecontain",
    )
    event_list = resolve_deprecated_alias(
        event_list,
        eventList,
        canonical_name="event_list",
        alias_name="eventList",
        function_name="seqecontain",
    )
    unknown_exclude = resolve_deprecated_alias(
        unknown_exclude,
        exclude,
        canonical_name="unknown_exclude",
        alias_name="exclude",
        function_name="seqecontain",
        default=False,
    )
    from ._mining import SubSeqElist  # _mining imports this module

    if isinstance(eseq, SubSeqElist):
        # seqecontain.R:6 -- a subsequence list is tested subsequence by
        # subsequence (2026-09-22 audit S5-06).
        event_codes = [_chain_codes(subsequence) for subsequence in eseq._subseq_objs]
        dictionary = eseq.eseq.dictionary
        scalar = False
    elif isinstance(eseq, SeqElist):
        event_codes = [list(sequence.codes) for sequence in eseq.sequences]
        dictionary = eseq.dictionary
        scalar = False
    elif isinstance(eseq, Eseq):
        event_codes = [list(eseq.codes)]
        dictionary = eseq.dictionary
        scalar = True
    else:
        raise TypeError(
            "eseq should be a SeqElist, Eseq or SubSeqElist (see seqecreate and seqefsub)."
        )

    # R's factor(event.list, levels = dict) (seqecontain.R:10): a scalar is a one-value
    # list, and each value is read as R's as.character reads it, so 1.0 names the event 1
    # (TX-AUD-20261019-044). Only a single value is wrapped: np.ndim is 0 for a set, a
    # generator or a dict view too, which must be read as collections (batch-11 review M-2).
    if isinstance(event_list, np.ndarray) and event_list.ndim == 0:
        event_list = [event_list.item()]
    elif isinstance(event_list, (str, bytes, Number, np.generic)):
        event_list = [event_list]
    code_of = {lab: i + 1 for i, lab in enumerate(dictionary)}
    labels = [_event_label(value) for value in event_list]
    wanted = {code_of[label] for label in labels if label in code_of}

    def _contain(codes: List[int]) -> bool:
        if len(codes) == 0:
            return False
        for c in codes:
            if int(c) in wanted:
                if not unknown_exclude:
                    return True
            elif unknown_exclude:
                return False
        return unknown_exclude

    out = np.array([_contain(codes) for codes in event_codes], dtype=bool)
    return bool(out[0]) if scalar else out


def _chain_codes(sequence: "_engine._Sequence") -> List[int]:
    """The event codes of an engine sequence, in chain order."""
    codes: List[int] = []
    node = sequence.event
    while node is not None:
        codes.append(node.type)
        node = node.next
    return codes


def is_eseq(x=MISSING, *, s=MISSING) -> bool:
    """Whether ``x`` is a single event sequence (R ``is.eseq``)."""
    x = resolve_deprecated_alias(
        x,
        s,
        canonical_name="x",
        alias_name="s",
        function_name="is_eseq",
    )
    return isinstance(x, Eseq)


def is_seqelist(x=MISSING, *, s=MISSING) -> bool:
    """Whether ``x`` is an event-sequence list (R ``is.seqelist``)."""
    x = resolve_deprecated_alias(
        x,
        s,
        canonical_name="x",
        alias_name="s",
        function_name="is_seqelist",
    )
    return isinstance(x, SeqElist)


def is_stslist(x: object) -> bool:
    """Whether ``x`` is a state-sequence object (R ``is.stslist``).

    True for every :class:`STSList`. R's ``is.stslist`` is FALSE for an object its own
    ``seqdef`` built with ``id=`` and ``weights=`` (it requires the weights' names to equal
    the row names, which seqdef does not give them): a recorded divergence
    (TX-AUD-20261011-005).
    """
    return isinstance(x, STSList)


def seqelength(eseq=MISSING, *, s=MISSING):
    """Per-sequence stored length / observation time (R ``seqelength``).

    Faithful port of ``TraMineR::seqelength`` (``getObsTime``): for an event
    sequence created from time-stamped events this is the *stored* length, which
    defaults to **-1** -- it is only set when the source is a state-sequence
    object or via the assignable :attr:`SeqElist.lengths` property (R
    ``seqelength<-``). A set length appears in ``as.character`` as a trailing gap
    to the end of observation (e.g. ``"(A)-2-(B)-3"`` for events at 0, 2 with
    length 5).

    NOTE: this is **not** the number of events per sequence (that prior,
    non-R behaviour is available as ``len(eseq[i].codes)``). Returns a 1-D float
    array for a :class:`SeqElist`, or a scalar for a single :class:`Eseq`.
    """
    eseq = resolve_deprecated_alias(
        eseq,
        s,
        canonical_name="eseq",
        alias_name="s",
        function_name="seqelength",
    )
    if isinstance(eseq, SeqElist):
        return np.array([s.obs_time for s in eseq.sequences], dtype=float)
    if isinstance(eseq, Eseq):
        return float(eseq.obs_time)
    raise TypeError("eseq should be a SeqElist or Eseq (see seqecreate).")


def seqeage(eseq, event_list) -> np.ndarray:
    """Event-age table for one event (R unexported ``seqeage`` / C ``tmreventinseq``).

    Faithful pure-Python port of ``TraMineR:::seqeage`` (``seqedist.R:7-9``), which
    calls the C routine ``tmreventinseq`` (``tmrsequence.cpp:430-484``). For a single
    event code it returns an ``(n_sequences, max_occurrences)`` float matrix whose
    cell ``[i, k]`` is the **age** (cumulative inter-event gap) of the ``(k+1)``-th
    occurrence of that event in sequence ``i``, or **-1** when the event occurs
    fewer than ``k+1`` times. The column width is ``max(1, max occurrences across
    sequences)`` (so an event that never occurs yields a single column of -1).

    Like the C routine, the age is re-added from the sequence's gap chain (built
    as TraMineR builds it, see :class:`Eseq`), so with fractional times it
    carries R's rounding -- an event at 2.0 after events at 0.1 and 0.2 has age
    ``1.9999999999999998`` -- while whole-number times give the absolute time.

    ``event_list`` is a 1-based event code or a list of them; like the C routine,
    **only the first element is used** (a multi-event ``event.list`` silently
    selects ``event_list[0]``). Returns a 2-D array for a :class:`SeqElist` or a
    single :class:`Eseq` (one row).
    """
    if isinstance(eseq, SeqElist):
        sequences = eseq.sequences
    elif isinstance(eseq, Eseq):
        sequences = [eseq]
    else:
        raise TypeError("eseq should be a SeqElist or Eseq (see seqecreate).")

    ev = np.atleast_1d(np.asarray(event_list)).ravel()
    if ev.size == 0:
        raise ValueError("event.list should have at least one element.")
    event = int(ev[0])

    per_seq: List[List[float]] = []
    for sequence in sequences:
        ages: List[float] = []
        age = 0.0
        node = _engine.eseq_to_engine(sequence).event
        while node is not None:
            age += node.gap
            if node.type == event:
                ages.append(age)
            node = node.next
        per_seq.append(ages)
    width = max(1, max((len(p) for p in per_seq), default=1))
    out = np.full((len(sequences), width), -1.0, dtype=float)
    for i, p in enumerate(per_seq):
        if len(p):
            out[i, : len(p)] = p
    return out
