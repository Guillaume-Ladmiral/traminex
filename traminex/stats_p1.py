"""Priority 1 parity functions: modal sequences, mean durations, position analysis, finding."""

from __future__ import annotations

from typing import Any, Sequence, cast

import numpy as np
import pandas as pd

from ._categorical_order import r_character_order, r_string_to_double
from .seqformat import _collision_free_marker
from .sequences import STSList, _require_stslist
from .stats_basic import seqistatd, seqstatd


@_require_stslist
def seqmodst(
    x: STSList,
    *,
    weighted: bool = True,
    include_missing: bool = False,
) -> STSList:
    """Modal state sequence (most frequent state at each position).

    Constructs a sequence of length L where each position contains the most
    frequent state at that position across all sequences. In case of ties,
    the first of the tied states in the alphabet order is selected, as R's
    ``which(freq == max(freq))[1]`` does (seqmodst.R:29); without an explicit
    ``alphabet=`` that order is R's ``seqdef`` order. A position without
    any valid cell (all excluded, or a zero total weight) has no modal state:
    R leaves an ``NA`` cell there, which traminex marks with the missing state,
    with frequency NaN; such a modal sequence matches no sequence.

    Parameters
    ----------
    x : STSList
        Input sequence object
    weighted : bool, default True
        If True, use sequence weights (if available)
    include_missing : bool, default False
        If True, include the missing state nr; structural void remains excluded

    Returns
    -------
    STSList
        A single-row STSList object representing the modal sequence.
        Attributes include:
        - frequencies_: array of shape (L,) with the frequency of the modal state at each position
        - occurrences_: int, number of sequences that exactly match the modal sequence

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqmodst
    >>> data = np.array([
    ...     ['A', 'A', 'B', 'C'],
    ...     ['A', 'B', 'B', 'C'],
    ...     ['A', 'A', 'B', 'B'],
    ... ])
    >>> sts = seqdef(data)
    >>> modal = seqmodst(sts)
    >>> modal.to_strings()[0]
    'A-A-B-C'
    """
    # seqstatd canonicalizes; index its columns with the same codebook.
    x = x._canonicalized()

    # Get position-wise state distributions. ``weighted`` must be forwarded:
    # ``seqstatd(weights=None)`` alone would fall back to the attached weights.
    freq_data = seqstatd(x, include_missing=include_missing, weighted=weighted)
    freq = freq_data["Frequencies"]  # shape (L, K)

    L = x.seqlen
    state_codes = x._codebook_codes(
        include_external_nr=include_missing,
        include_external_void=False,
    )

    # Build modal sequence: most frequent state at each position
    modal_columns = np.zeros(L, dtype=int)
    modal_freqs = np.zeros(L, dtype=float)
    valid = np.asarray(freq_data["ValidStates"], dtype=float)

    for t in range(L):
        if not valid[t] > 0:
            # No valid state here (every cell excluded, or zero total weight):
            # R's which(freq == max(freq))[1] is NA on the all-NaN column, so
            # its modal cell and frequency are NA (seqmodst.R:28-32).  STSList
            # has no NA cell; the missing marker stands for it (column -1).
            modal_columns[t] = -1
            modal_freqs[t] = np.nan
            continue
        # Ties go to the first state in the alphabet order (seqmodst.R:29).
        idx_max = int(np.argmax(freq[t, :]))
        modal_columns[t] = idx_max
        modal_freqs[t] = freq[t, idx_max]

    # Code -1 marks a position without a modal state; it matches no sequence
    # below, as R's NA cell never matches in seqfind.
    modal_source_codes = np.array(
        [state_codes[column] if column >= 0 else -1 for column in modal_columns],
        dtype=int,
    )
    # Keep logical codes and their metadata together. Re-running seqdef on
    # decoded strings would reinterpret real '*'/'%' states as markers and
    # reorder their colors. Undefined modal cells need an actual missing code.
    modal_data = modal_source_codes.reshape(1, -1).copy()
    result_alphabet = list(x.alphabet)
    result_labels = None if x.labels is None else list(x.labels)
    result_palette = None if x.cpal is None else list(x.cpal)
    nr, nr_code = x.nr, x.nr_code
    undefined = modal_data < 0
    if undefined.any():
        if nr_code is None:
            nr = _collision_free_marker(nr, [*result_alphabet, x.void], "NR")
            nr_code = len(result_alphabet)
            result_alphabet.append(nr)
            if result_labels is not None:
                result_labels.append(nr)
            if result_palette is not None:
                result_palette.extend(
                    [x.missing_color] * (len(result_alphabet) - len(result_palette))
                )
        modal_data[undefined] = nr_code
    result = STSList(
        data=modal_data,
        alphabet=result_alphabet,
        labels=result_labels,
        cpal=result_palette,
        nr=nr,
        nr_code=nr_code,
        void=x.void,
        void_code=x.void_code,
        xtlab=None if x.xtlab is None else list(x.xtlab),
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        start=x.start,
        missing_color=x.missing_color,
    )

    # Count how many sequences exactly match the modal sequence
    matches = np.all(x.data == modal_source_codes[np.newaxis, :], axis=1)
    n_occ = int(np.sum(matches))

    # Attach attributes (following R convention)
    setattr(result, "frequencies_", modal_freqs)
    setattr(result, "occurrences_", n_occ)
    # R: FALSE without weights, with weighted=FALSE, or when every weight is 1
    # (seqmodst.R:57-66; audit S3-16).
    setattr(
        result,
        "weighted_",
        bool(
            weighted
            and x.weights is not None
            and not np.all(np.asarray(x.weights, dtype=float) == 1.0)
        ),
    )
    setattr(
        result,
        "nbseq_",
        float(np.sum(x.weights)) if weighted and x.weights is not None else float(x.nseq),
    )

    return result


@_require_stslist
def seqmeant(
    x: STSList,
    *,
    weighted: bool = True,
    include_missing: bool = False,
    prop: bool = False,
    serr: bool = False,
) -> pd.DataFrame:
    """Mean durations (time spent) in each state.

    Computes the weighted mean time spent in each state across all sequences.
    Optionally returns variance, standard deviation, and standard error.

    Parameters
    ----------
    x : STSList
        Input sequence object
    weighted : bool, default True
        If True, use sequence weights (if available)
    include_missing : bool, default False
        If True, count the missing state nr. Void cells are never counted, but the
        result still has a row for each marker the object stores (see Returns).
    prop : bool, default False
        If True, return proportions (mean durations as proportion of sequence length).
        As in R, a sequence without any valid cell makes every mean NaN.
    serr : bool, default False
        If True, also compute variance, standard deviation, and standard error

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per entry of ``x.alphabet``, so also a row for a stored
        void or nr marker, which is 0 when its cells are not counted (R returns rows for
        the real states, plus nr when missing is included; docs/API_PARITY_DECISIONS.csv),
        and columns:
        - Mean: mean duration (or proportion if prop=True)
        - Var: variance (if serr=True)
        - Stdev: standard deviation (if serr=True)
        - SE: standard error (if serr=True)

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqmeant
    >>> data = np.array([
    ...     ['A', 'A', 'B', 'C'],
    ...     ['A', 'B', 'B', 'C'],
    ... ])
    >>> sts = seqdef(data)
    >>> seqmeant(sts)
              Mean
    A         1.5
    B         1.5
    C         1.0
    """
    # seqistatd canonicalizes; index its columns with the same codebook.
    x = x._canonicalized()

    # Get per-sequence state durations
    istatd = np.asarray(seqistatd(x, include_missing=include_missing), dtype=float)  # (n, K)
    state_codes = x._codebook_codes(
        include_external_nr=include_missing,
        include_external_void=False,
    )

    if prop:
        # R divides each sequence's durations by its own total over statl (the
        # real alphabet, + nr when with.missing; never void) BEFORE averaging,
        # so ragged / missing rows weight by their valid length, not the
        # rectangular L. (TraMineR seqmeant -> seqistatd(prop=TRUE).)
        included_columns = [
            column
            for column, code in enumerate(state_codes)
            if code != x.void_code and (include_missing or code != x.nr_code)
        ]
        denom = istatd[:, included_columns].sum(axis=1, keepdims=True)
        # A sequence without any valid cell is 0/0 = NaN in R's
        # seqistatd(prop=TRUE), which makes every mean NaN; the marker-padded
        # columns keep their zeros.
        with np.errstate(divide="ignore", invalid="ignore"):
            istatd[:, included_columns] = istatd[:, included_columns] / denom

    # Resolve weights
    if weighted and x.weights is not None:
        w = np.asarray(x.weights, dtype=float)
    else:
        w = np.ones(x.nseq, dtype=float)
        weighted = False
    if np.all(w == 1.0):
        # R: all-one weights are recorded as unweighted (seqmeant.R:19-20; audit S3-16).
        weighted = False

    wtot = float(w.sum())

    # Weighted mean durations
    mtime = (istatd * w[:, np.newaxis]).sum(axis=0) / wtot

    result = pd.DataFrame(
        {"Mean": mtime},
        index=[x._state_for_code(code) for code in state_codes],
    )

    if serr:
        # Compute variance, stdev, SE
        w2tot = float((w**2).sum())
        # Centered values
        vcent = istatd - mtime[np.newaxis, :]
        # Weighted variance (with finite population correction)
        var = (w[:, np.newaxis] * (vcent**2)).sum(axis=0) * wtot / (wtot**2 - w2tot)
        sd = np.sqrt(var)
        se = np.sqrt(var / wtot)

        result["Var"] = var
        result["Stdev"] = sd
        result["SE"] = se

    # Attach metadata (following R convention)
    result.attrs["nbseq"] = wtot
    result.attrs["weighted"] = weighted

    return result


@_require_stslist
def seqstatf(
    x: STSList,
    *,
    weighted: bool = True,
    include_missing: bool = False,
) -> pd.DataFrame:
    """Overall state frequency table (TraMineR ``seqstatf``).

    Computes the weighted total time spent in each state across all sequences,
    i.e. ``colSums(seqistatd(seqdata, with.missing) * weights)`` -- NOT the
    first-position distribution. Index is ``statl`` (the alphabet, plus the
    missing code ``nr`` when ``include_missing``, even if no cell is missing);
    the structural ``void`` symbol is never a state. As in R, ``Percent`` is
    NaN when the total weight is zero.

    Parameters
    ----------
    x : STSList
        Input sequence object
    weighted : bool, default True
        If True, use sequence weights (if available)
    include_missing : bool, default False
        If True, add the missing code (``nr``) to the state set (R
        ``with.missing``)

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by state with columns ``Freq`` (weighted total
        occurrences) and ``Percent`` (``Freq / sum(Freq) * 100``).

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqstatf
    >>> data = np.array([
    ...     ['A', 'A', 'B', 'C'],
    ...     ['A', 'B', 'B', 'C'],
    ...     ['B', 'A', 'B', 'B'],
    ... ])
    >>> sts = seqdef(data)
    >>> seqstatf(sts)['Freq'].tolist()
    [4.0, 6.0, 2.0]
    """
    # Both marker representations are equivalent inputs (TX-AUD-20261010-094).
    x = x._canonicalized()
    istatd = np.asarray(seqistatd(x, include_missing=include_missing), dtype=float)  # (n, K)

    if weighted and x.weights is not None:
        w = np.asarray(x.weights, dtype=float)
    else:
        w = np.ones(x.nseq, dtype=float)

    freq_all = (istatd * w[:, np.newaxis]).sum(axis=0)  # length K (full alphabet)

    # statl: real alphabet, + nr when include_missing; never void.
    state_codes = x._codebook_codes(
        include_external_nr=include_missing,
        include_external_void=False,
    )
    selected_codes = [
        code
        for code in state_codes
        if code != x.void_code and (include_missing or code != x.nr_code)
    ]
    selected_columns = [state_codes.index(code) for code in selected_codes]
    freq = freq_all[selected_columns]
    labels = [x._state_for_code(code) for code in selected_codes]
    if include_missing and (x.nr_code is None or x.nr_code not in selected_codes):
        # R's seqistatd has an nr column under with.missing even when no cell
        # is missing (seqistatd.R:12-14), so seqstatf keeps a zero nr row.
        freq = np.append(freq, 0.0)
        labels.append(x.nr)
    total = freq.sum()
    # R's Freq / sum(Freq) * 100 is 0/0 = NaN for a zero total weight.
    with np.errstate(divide="ignore", invalid="ignore"):
        percent = freq / total * 100.0

    return pd.DataFrame(
        {"Freq": freq, "Percent": percent},
        index=labels,
    )


def seqstatl(
    data,
    var: Sequence[int | str] | int | str | None = None,
    format: str = "STS",
) -> list[str] | list[float]:
    """Return the distinct observed states using TraMineR's public contract.

    Raw matrices, data frames, compressed ``"-"``/``":"`` strings, SPS
    strings, and :class:`STSList` objects are accepted. ``var`` uses one-based
    integer or column-name selection. Exact source format spellings are
    ``"STS"``, ``"SPS"``, ``"SPELL"``, and ``"DSS"``.

    Missing raw cells are omitted. Actual missing and structural-void symbols
    stored in an ``STSList`` remain observable. If every distinct state
    converts numerically, the result contains sorted floats, matching R's
    double-vector return; otherwise it contains strings ordered with the R/ICU
    character collation.
    """
    from .seqformat import (
        _cell_is_na,
        _r_state_text,
        _stslist_column_positions,
        seqformat,
    )
    from .stats_p2 import _seqconc_raw_matrix, _seqconc_var, seqdecomp

    formats = ("STS", "SPS", "SPELL", "DSS")
    if format not in formats:
        raise ValueError("Input format must be one of: STSSPSSPELLDSS")

    if isinstance(data, STSList):
        normalized_var = _seqconc_var(var)
        if normalized_var == 0:
            positions: list[int] = []
        else:
            positions = _stslist_column_positions(
                data,
                cast(Sequence[int | str] | int | str | None, normalized_var),
            )
        matrix = np.empty((data.nseq, len(positions)), dtype=object)
        for index, code in np.ndenumerate(data.data[:, positions]):
            matrix[index] = data._state_for_code(int(code))
    else:
        matrix, _ = _seqconc_raw_matrix(data, var)

    if format == "SPS":
        converted = seqformat(matrix, from_="SPS", to="STS")
        assert isinstance(converted, STSList)
        labels = np.asarray(converted.alphabet, dtype=object)
        matrix = labels[converted.data]

    nonmissing = [value for value in matrix.ravel(order="F") if not _cell_is_na(value)]
    numeric_matrix = bool(nonmissing) and all(
        isinstance(value, (int, float, np.integer, np.floating))
        and not isinstance(value, (bool, np.bool_))
        for value in nonmissing
    )
    separator = None
    if not numeric_matrix and matrix.shape[1] == 1:
        texts = [_r_state_text(value) for value in nonmissing]
        if any("-" in value for value in texts):
            separator = "-"
        elif any(":" in value for value in texts):
            separator = ":"
    if separator is not None:
        matrix = seqdecomp(matrix, sep=separator)

    distinct: list[object] = []
    for value in matrix.ravel(order="F"):
        if _cell_is_na(value) or value in distinct:
            continue
        distinct.append(value)

    numeric: list[float] = []
    for value in distinct:
        if isinstance(value, (bool, np.bool_)):
            converted_value = float(bool(value))
        elif isinstance(value, str):
            # R's as.numeric, not Python's float(): "1_1" is text, "0x10" is 16
            # (TX-AUD-20261010-050).
            parsed = r_string_to_double(value)
            if parsed is None:
                break
            converted_value = parsed
        else:
            try:
                converted_value = float(cast(Any, value))
            except (TypeError, ValueError):
                break
        if np.isnan(converted_value):
            break
        numeric.append(converted_value)
    else:
        return sorted(numeric)

    return r_character_order([_r_state_text(value) for value in distinct])


def seqfind(ref: STSList, target: STSList) -> np.ndarray:
    """Find occurrences of reference sequence(s) in target sequences.

    Matches each reference in its original order against all target rows in
    target order. Structural void is removed independently from each object;
    missing states remain part of the comparison. Duplicate references repeat
    their complete target match list, as in TraMineR.

    Parameters
    ----------
    ref : STSList
        Reference sequence(s) to search for
    target : STSList
        Target sequence set to search in

    Returns
    -------
    np.ndarray
        One-dimensional integer array of zero-based target indexes. This is the
        maintained Python adaptation of R's one-based row indexes.

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqfind
    >>> data1 = np.array([['A', 'B', 'C']])
    >>> data2 = np.array([
    ...     ['A', 'B', 'C'],
    ...     ['A', 'B', 'B'],
    ...     ['A', 'B', 'C'],
    ... ])
    >>> ref = seqdef(data1)
    >>> target = seqdef(data2)
    >>> seqfind(ref, target)
    array([0, 2])
    """
    if not isinstance(ref, STSList):
        raise TypeError("ref must be an STSList state sequence object")
    if not isinstance(target, STSList):
        raise TypeError("target must be an STSList state sequence object")

    from .stats_p2 import seqconc

    references = seqconc(ref)
    targets = seqconc(target)
    matches = [
        target_index
        for reference in references
        for target_index, candidate in enumerate(targets)
        if candidate == reference
    ]
    return np.asarray(matches, dtype=int)


def seqfpos(
    x: STSList,
    state: str | None = None,
    *,
    states: Sequence[str] | None = None,
    include_missing: bool = False,
) -> np.ndarray | pd.DataFrame:
    """First occurrence positions for one state or a vectorized state set.

    Supplying ``state`` follows TraMineR and returns a vector containing the
    raw 1-based rectangular position, with NaN when absent. Omitting ``state``
    preserves the existing DataFrame helper (0 when a selected state is absent).

    Parameters
    ----------
    x : STSList
        Input sequence object
    states : list of str, optional
        Subset of states to analyze. If None, use all states.
    include_missing : bool, default False
        If True, include the missing state (in the table mode, as a column). The
        void marker is not a state and never gets a column.

    Returns
    -------
    pd.DataFrame
        DataFrame of shape (nseq, n_states) with first occurrence positions (1-based)

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqfpos
    >>> data = np.array([
    ...     ['A', 'B', 'C', 'B'],
    ...     ['B', 'A', 'C', 'C'],
    ...     ['C', 'C', 'A', 'A'],
    ... ])
    >>> sts = seqdef(data)
    >>> seqfpos(sts)
       A  B  C
    0  1  2  3
    1  2  1  3
    2  3  0  1
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")
    x = x._canonicalized()
    if state is not None and states is not None:
        raise ValueError("state and states cannot be supplied together")

    if state is not None:
        result = np.full(x.nseq, np.nan, dtype=float)
        try:
            code = x.alphabet.index(state)
        except ValueError:
            return result
        for i, row in enumerate(x.data):
            positions = np.flatnonzero(row == code)
            if positions.size:
                result[i] = float(positions[0] + 1)
        return result

    # Compatibility DataFrame mode: determine which states to analyze.
    label_to_code = {lbl: idx for idx, lbl in enumerate(x.alphabet)}
    if states is None:
        target_codes = list(range(len(x.alphabet)))
    else:
        target_codes = [label_to_code[s] for s in states if s in label_to_code]

    # The missing state when include_missing; void never, since it is not part of
    # a sequence (TX-AUD-20261019-070).
    excluded = {int(x.void_code)} if x.void_code is not None else set()
    if not include_missing and x.nr_code is not None:
        excluded.add(int(x.nr_code))
    target_codes = [c for c in target_codes if c not in excluded]
    state_list = [x.alphabet[c] for c in target_codes]

    # Compute first occurrence positions
    result = np.zeros((x.nseq, len(target_codes)), dtype=int)

    for i, row in enumerate(x.data):
        for j, code in enumerate(target_codes):
            # Find first occurrence (1-based, 0 if not found)
            pos = np.where(row == code)[0]
            if pos.size > 0:
                result[i, j] = int(pos[0]) + 1  # 1-based
            else:
                result[i, j] = 0

    return pd.DataFrame(result, columns=state_list)


def seqfposend(
    x: STSList,
    state: str | None = None,
    *,
    states: Sequence[str] | None = None,
    include_missing: bool = False,
    lead: int = 0,
    from_seq_start: bool = True,
) -> np.ndarray | pd.DataFrame:
    """End of the first spell for one state, or compatibility DataFrame.

    Supplying ``state`` follows TraMineR's spell-duration contract and returns
    0 when that state has no spell. Omitting ``state`` preserves the existing
    vectorized DataFrame of last rectangular occurrences.

    Parameters
    ----------
    x : STSList
        Input sequence object
    states : list of str, optional
        Subset of states to analyze. If None, use all states.
    include_missing : bool, default False
        If True, include the missing state (in the table mode, as a column). The
        void marker is not a state and never gets a column.

    Returns
    -------
    pd.DataFrame
        DataFrame of shape (nseq, n_states) with last occurrence positions (1-based)

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqfposend
    >>> data = np.array([
    ...     ['A', 'B', 'C', 'B'],
    ...     ['B', 'A', 'C', 'C'],
    ...     ['C', 'C', 'A', 'A'],
    ... ])
    >>> sts = seqdef(data)
    >>> seqfposend(sts)
       A  B  C
    0  1  4  3
    1  2  1  4
    2  4  0  2
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")
    x = x._canonicalized()
    if state is not None and states is not None:
        raise ValueError("state and states cannot be supplied together")

    if state is not None:
        from .stats_complexity import seqdss, seqdur

        spells = seqdss(x, include_missing=include_missing)
        positions = seqfpos(spells, state)
        durations = seqdur(x, include_missing=include_missing)
        if from_seq_start:
            durations = np.cumsum(durations, axis=1)

        result = np.zeros(x.nseq, dtype=float)
        for i, position in enumerate(positions):
            if not np.isnan(position):
                result[i] = durations[i, int(position) - 1] + lead
        return result

    # Compatibility DataFrame mode: determine which states to analyze.
    label_to_code = {lbl: idx for idx, lbl in enumerate(x.alphabet)}
    if states is None:
        target_codes = list(range(len(x.alphabet)))
    else:
        target_codes = [label_to_code[s] for s in states if s in label_to_code]

    # The missing state when include_missing; void never, since it is not part of
    # a sequence (as seqfpos's table, TX-REP-20261023-004).
    excluded = {int(x.void_code)} if x.void_code is not None else set()
    if not include_missing and x.nr_code is not None:
        excluded.add(int(x.nr_code))
    target_codes = [c for c in target_codes if c not in excluded]
    state_list = [x.alphabet[c] for c in target_codes]

    # Compute last occurrence positions
    result = np.zeros((x.nseq, len(target_codes)), dtype=int)

    for i, row in enumerate(x.data):
        for j, code in enumerate(target_codes):
            # Find last occurrence (1-based, 0 if not found)
            pos = np.where(row == code)[0]
            if pos.size > 0:
                result[i, j] = int(pos[-1]) + 1  # 1-based
            else:
                result[i, j] = 0

    return pd.DataFrame(result, columns=state_list)
