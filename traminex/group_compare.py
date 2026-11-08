"""LRT / BIC group comparison from dissimilarities — R-faithful.

BUILD-from-R port of TraMineRextras 0.6.9 ``dissCompare`` / ``seqCompare`` /
``seqLRT`` / ``seqBIC`` (``SOURCES/TraMineRextras/R/dissCompare.R`` and
``seqCompare.R``). The shared worker ``_seqxcomp`` computes the likelihood-ratio
test and the BIC / Bayes-factor of a two-group split from a dissimilarity matrix
(Studer & Liao's procedure). ``dissCompare`` operates on a precomputed
dissimilarity matrix; ``seqCompare`` / ``seqLRT`` / ``seqBIC`` compute the
distance internally (see :mod:`traminex.group_compare` -- ``seqCompare`` lands in
the seqCompare commit).

Both the deterministic ``s=0`` route and the repeated balanced-sampling route
are implemented. NumPy generators replace R's global RNG while preserving the
draw construction and aggregation mechanics.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Callable, Literal, Optional, cast

import numpy as np
import pandas as pd

from ._categorical_order import r_factor
from .diss import disscenter
from .distances import seqdist
from .sequences import STSList


def _resolve_squared(squared) -> tuple[bool, int]:
    """Map R's ``squared`` (logical or ``"LRTonly"``) to ``(squared, LRTpow)``."""
    if isinstance(squared, str):
        if squared != "LRTonly":
            raise ValueError("squared must be a bool or 'LRTonly'")
        return False, 2
    if not isinstance(squared, (bool, np.bool_)):
        raise TypeError("squared must be a bool or 'LRTonly'")
    return bool(squared), 1


def _seqxcomp(
    r1: np.ndarray,
    r2: np.ndarray,
    diss: np.ndarray,
    weights: Optional[np.ndarray],
    *,
    is_lrt: bool,
    is_bic: bool,
    squared: bool = False,
    weighted: bool = True,
    weight_by: str = "global",
    lrt_pow: int = 2,
) -> np.ndarray:
    """Port of ``seqxcomp`` (LRT / BIC of a two-group split of ``diss``)."""
    from scipy import stats as _spstats

    r1 = np.asarray(r1, dtype=int)
    r2 = np.asarray(r2, dtype=int)
    n1, n2 = r1.size, r2.size
    n0 = n1 + n2

    weighted = weighted and weights is not None
    if weighted:
        w = np.asarray(weights, dtype=float)
        w1 = w[r1]
        w2 = w[r2]
        if weight_by == "by.group":
            w1 = n1 / w1.sum() * w1
            w2 = n2 / w2.sum() * w2
        wcat = np.concatenate([w1, w2])
    else:
        w1 = np.ones(n1)
        w2 = np.ones(n2)
        wcat = np.ones(n0)

    both = np.concatenate([r1, r2])
    dist_s = np.asarray(
        disscenter(diss[np.ix_(both, both)], weights=wcat, squared=squared), dtype=float
    )
    dist_s1 = np.asarray(disscenter(diss[np.ix_(r1, r1)], weights=w1, squared=squared), dtype=float)
    dist_s2 = np.asarray(disscenter(diss[np.ix_(r2, r2)], weights=w2, squared=squared), dtype=float)

    ss = float(np.sum(wcat * dist_s**lrt_pow))
    ss1 = float(np.sum(w1 * dist_s1**lrt_pow))
    ss2 = float(np.sum(w2 * dist_s2**lrt_pow))

    del is_lrt, is_bic  # R's worker always computes all four intermediates.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        lrt = n0 * (np.log(ss / n0) - np.log((ss1 + ss2) / n0))
        p_lrt = float(_spstats.chi2.sf(lrt, 1))
        bic = lrt - np.log(n0)
        bf = np.exp(bic / 2.0)
    return np.array([lrt, p_lrt, bic, bf], dtype=float)


def _stat_columns(stat: str) -> tuple[bool, bool, list[str]]:
    if stat == "all":
        is_lrt = is_bic = True
    else:
        is_lrt = stat == "LRT"
        is_bic = stat == "BIC"
    if not (is_lrt or is_bic):
        raise ValueError("stat must be one of 'LRT', 'BIC', or 'all'")
    cols: list[str] = []
    if is_lrt:
        cols += ["LRT", "p-value"]
    if is_bic:
        cols += ["Delta BIC", "Bayes Factor"]
    return is_lrt, is_bic, cols


def _build_sampling_plan(
    max_size: int,
    min_size: int,
    sample_size: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Build R's repeated balanced relative-index matrices (zero-based)."""
    base_repeats = sample_size // max_size
    max_repeats = max(1, base_repeats)
    if sample_size < max_size:
        max_remainder = sample_size - max_size % sample_size
    else:
        max_remainder = sample_size - base_repeats * max_size

    total = max_repeats * max_size + max_remainder
    min_repeats = total // min_size
    min_remainder = total - min_repeats * min_size

    first = rng.permutation(np.tile(np.arange(max_size, dtype=int), max_repeats))
    if max_remainder:
        first = np.concatenate([first, rng.choice(max_size, size=max_remainder, replace=False)])
    second = rng.permutation(np.tile(np.arange(min_size, dtype=int), min_repeats))
    if min_remainder:
        second = np.concatenate([second, rng.choice(min_size, size=min_remainder, replace=False)])
    return (
        first.reshape((-1, sample_size), order="F"),
        second.reshape((-1, sample_size), order="F"),
    )


@dataclass(frozen=True)
class _CompareSet:
    label: str
    first: np.ndarray
    second: np.ndarray


def _compare_factor(values: np.ndarray, valid: np.ndarray) -> tuple[list[str], np.ndarray]:
    """R's ``factor()`` of the valid values (dissCompare.R:92-105): labels and full-length codes.

    Levels as R names and groups them (numbers as R prints them, equal to 15
    significant digits one level, TX-AUD-20261019-046); codes are -1 outside ``valid``.
    """
    labels, codes, _ = r_factor(values[valid])
    full = np.full(values.shape[0], -1, dtype=np.int64)
    full[valid] = codes
    return labels, full


def _sampling_plans(
    compare_sets: Sequence[_CompareSet],
    sample_size: int,
    seed: int | np.random.Generator | None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Precompute every set's draw plan before any statistic or distance."""
    shared_rng: np.random.Generator | None
    if isinstance(seed, np.random.Generator):
        shared_rng = seed
    elif seed is None:
        shared_rng = np.random.default_rng()
    elif isinstance(seed, (int, np.integer)) and not isinstance(seed, (bool, np.bool_)):
        shared_rng = None
    else:
        raise TypeError("seed must be an integer, numpy.random.Generator, or None")

    plans: list[tuple[np.ndarray, np.ndarray]] = []
    for comparison in compare_sets:
        first_size = comparison.first.size
        second_size = comparison.second.size
        if min(first_size, second_size) == 0:
            plans.append((np.empty((0, 0), dtype=int), np.empty((0, 0), dtype=int)))
        elif sample_size == 0:
            plans.append(
                (
                    np.arange(first_size, dtype=int)[None, :],
                    (np.arange(second_size, dtype=int) + first_size)[None, :],
                )
            )
        else:
            rng = shared_rng if shared_rng is not None else np.random.default_rng(cast(int, seed))
            larger_draws, smaller_draws = _build_sampling_plan(
                max(first_size, second_size),
                min(first_size, second_size),
                sample_size,
                rng,
            )
            # The draws from 1:max(n) belong to the larger group, which seqCompare
            # makes the first (seqCompare.R:143-157). dissCompare keeps the first
            # level first (dissCompare.R:117-131) and so samples the other group's
            # cases when the first level is the smaller group (TX-AUD-20261011-022).
            if first_size >= second_size:
                plans.append((larger_draws, smaller_draws + first_size))
            else:
                plans.append((larger_draws + first_size, smaller_draws))
    return plans


def _aggregate_comparisons(
    compare_sets: Sequence[_CompareSet],
    plans: Sequence[tuple[np.ndarray, np.ndarray]],
    evaluator: Callable[[int, np.ndarray, np.ndarray], np.ndarray],
    *,
    is_lrt: bool,
    is_bic: bool,
    BFopt: Literal[1, 2] | None,
    inc_df: bool,
) -> pd.DataFrame:
    """Shared ``dissCompare`` aggregation for matrix and sequence routes."""
    raw = np.full((len(compare_sets), 4), np.nan, dtype=float)
    dfs = np.zeros(len(compare_sets), dtype=float)
    multiple_samples = False
    for index, plan in enumerate(plans):
        first_draws, second_draws = plan
        if first_draws.shape[0] == 0:
            continue
        per_draw = np.empty((first_draws.shape[0], 4), dtype=float)
        for draw in range(first_draws.shape[0]):
            per_draw[draw] = evaluator(index, first_draws[draw], second_draws[draw])
        raw[index] = np.mean(per_draw, axis=0)
        draw_count = first_draws.shape[0]
        dfs[index] = draw_count
        multiple_samples = multiple_samples or draw_count > 1
        if is_lrt and draw_count > 1:
            from scipy import stats as _spstats

            raw[index, 0] *= draw_count
            raw[index, 1] = _spstats.chi2.sf(raw[index, 0], draw_count)

    columns = ["LRT", "p-value", "Delta BIC", "Bayes Factor"]
    result_values = raw
    if is_bic and multiple_samples:
        if BFopt is None:
            with np.errstate(over="ignore", invalid="ignore"):
                from_average_bic = np.exp(raw[:, 2] / 2.0)
            result_values = np.column_stack([raw, from_average_bic])
            columns = [
                "LRT",
                "p-value",
                "Delta BIC",
                "Bayes Factor (Avg)",
                "Bayes Factor (From Avg BIC)",
            ]
        elif BFopt == 1:
            columns[3] = "Bayes Factor (Avg)"
        else:
            result_values = raw.copy()
            with np.errstate(over="ignore", invalid="ignore"):
                result_values[:, 3] = np.exp(raw[:, 2] / 2.0)
            columns[3] = "Bayes Factor (From Avg BIC)"

    if inc_df:
        result_values = np.column_stack([result_values, dfs])
        columns.append("df")

    keep: list[int] = []
    if is_lrt:
        keep.extend([0, 1])
    if is_bic:
        keep.extend(range(2, len(columns) - int(bool(inc_df))))
    if inc_df:
        keep.append(len(columns) - 1)
    return pd.DataFrame(
        result_values[:, keep],
        columns=[columns[column] for column in keep],
        index=[comparison.label for comparison in compare_sets],
    )


def dissCompare(
    diss: np.ndarray,
    group: Sequence[object],
    *,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    stat: Literal["all", "LRT", "BIC"] = "all",
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    weights: Sequence[float] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
) -> pd.DataFrame:
    """Compare two groups within one or more sets from dissimilarities.

    Nonzero ``s`` constructs the same repeated balanced samples as
    TraMineRextras and aggregates their LRT/BIC statistics. Integer seeds are
    reset for each set as in R; a supplied NumPy generator is consumed
    directly. No route touches NumPy's global random state.
    """
    D = np.asarray(diss, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("diss must be a square matrix")
    n = D.shape[0]
    group_values = np.asarray(group, dtype=object)
    if group_values.ndim != 1 or group_values.shape[0] != n:
        raise ValueError("group length not equal to nrow(diss)!")
    if set is None:
        set_values = np.ones(n, dtype=int)
    else:
        set_values = np.asarray(set, dtype=object)
        if set_values.ndim != 1 or set_values.shape[0] != n:
            raise ValueError("set length must equal nrow(diss)")

    if isinstance(s, (bool, np.bool_)) or not isinstance(s, (int, np.integer)):
        raise TypeError("s must be a nonnegative integer")
    s = int(s)
    if s < 0:
        raise ValueError("s must be a nonnegative integer")
    if BFopt not in {None, 1, 2}:
        raise ValueError("BFopt must be None, 1, or 2")
    if not isinstance(inc_df, (bool, np.bool_)):
        raise TypeError("inc_df must be a bool")

    sq, lrt_pow = _resolve_squared(squared)
    is_lrt, is_bic, _ = _stat_columns(stat)

    if isinstance(weighted, (bool, np.bool_)):
        use_weights = bool(weighted)
        weight_by = "global"
    elif weighted == "by.group":
        use_weights = True
        weight_by = "by.group"
    else:
        raise ValueError("weighted must be logical or 'by.group'")
    if weights is None or not use_weights:
        w = np.ones(n, dtype=float)
    else:
        w = np.asarray(weights, dtype=float)
        if w.ndim != 1 or w.shape[0] != n:
            raise ValueError("weights must have length equal to nrow(diss)")

    valid = ~pd.isna(group_values) & ~pd.isna(set_values)
    removed = int(n - np.count_nonzero(valid))
    if removed:
        warnings.warn(
            f"[!!] {removed} cases removed because of NA values of the grouping variable(s)",
            UserWarning,
            stacklevel=2,
        )
    group_labels, group_codes = _compare_factor(group_values, valid)
    if len(group_labels) < 2:
        raise ValueError("There is only one group among valid cases!")
    if len(group_labels) > 2:
        raise ValueError("Currently dissCompare supports only 2 groups!")
    set_labels, set_codes = _compare_factor(set_values, valid)

    compare_sets: list[_CompareSet] = []
    for set_code, set_label in enumerate(set_labels):
        in_set = set_codes == set_code
        first = np.flatnonzero(in_set & (group_codes == 0))
        second = np.flatnonzero(in_set & (group_codes == 1))
        compare_sets.append(_CompareSet(set_label, first, second))

    plans = _sampling_plans(compare_sets, s, seed)
    prepared: list[tuple[np.ndarray, np.ndarray]] = []
    for comparison in compare_sets:
        order = np.concatenate([comparison.first, comparison.second])
        diss_i = D[np.ix_(order, order)]
        weights_i = w[order]
        prepared.append((diss_i, weights_i))

    def evaluate(index: int, first: np.ndarray, second: np.ndarray) -> np.ndarray:
        diss_i, weights_i = prepared[index]
        return _seqxcomp(
            first,
            second,
            diss_i,
            weights_i,
            is_lrt=is_lrt,
            is_bic=is_bic,
            squared=sq,
            weighted=use_weights,
            weight_by=weight_by,
            lrt_pow=lrt_pow,
        )

    return _aggregate_comparisons(
        compare_sets,
        plans,
        evaluate,
        is_lrt=is_lrt,
        is_bic=is_bic,
        BFopt=BFopt,
        inc_df=bool(inc_df),
    )


def _subset_sequence(x: STSList, rows: np.ndarray) -> STSList:
    indexes = np.asarray(rows, dtype=int)
    return STSList(
        data=x.data[indexes].copy(),
        alphabet=list(x.alphabet),
        labels=None if x.labels is None else list(x.labels),
        cpal=None if x.cpal is None else list(x.cpal),
        xtlab=None if x.xtlab is None else list(x.xtlab),
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        nr=x.nr,
        void=x.void,
        nr_code=x.nr_code,
        void_code=x.void_code,
        weights=None if x.weights is None else x.weights[indexes].copy(),
        start=x.start,
        ids=None if x.ids is None else [x.ids[int(index)] for index in indexes],
        missing_color=x.missing_color,
        issues=None,
    )


def _recode_sequence_data(
    x: STSList,
    alphabet: Sequence[str],
    *,
    nr_code: int,
    void_code: int,
) -> np.ndarray:
    state_codes = {state: index for index, state in enumerate(alphabet)}
    result = np.empty_like(x.data, dtype=int)
    for index, code in np.ndenumerate(x.data):
        numeric = int(code)
        if x.nr_code is not None and numeric == x.nr_code:
            result[index] = nr_code
        elif x.void_code is not None and numeric == x.void_code:
            result[index] = void_code
        elif 0 <= numeric < len(x.alphabet):
            result[index] = state_codes[x.alphabet[numeric]]
        else:
            raise ValueError(f"sequence data contains unknown state code {numeric}")
    return result


def _combine_sequences(first: STSList, second: STSList) -> STSList:
    if first.nr != second.nr or first.void != second.void:
        raise ValueError("nr and/or void mismatch between STSList objects")
    first_states = set(first.alphabet)
    second_states = set(second.alphabet)
    if first_states.issubset(second_states):
        provider = second
    elif second_states.issubset(first_states):
        provider = first
    else:
        raise ValueError("alphabet mismatch between STSList objects")
    alphabet = list(provider.alphabet)
    from .seqformat import _collision_free_marker

    # The combined object declares external marker slots, including padding.
    # When neither input declares a marker, its default text can name a real
    # state. Give the new slot a distinct text so canonicalization cannot
    # merge that real state into missing or void.
    nr = first.nr
    if first.nr_code is None and second.nr_code is None:
        nr = _collision_free_marker(nr, [*alphabet, first.void], "NR")
    void = first.void
    if first.void_code is None and second.void_code is None:
        void = _collision_free_marker(void, [*alphabet, nr], "VOID")
    nr_code = len(alphabet)
    void_code = nr_code + 1
    first_data = _recode_sequence_data(first, alphabet, nr_code=nr_code, void_code=void_code)
    second_data = _recode_sequence_data(second, alphabet, nr_code=nr_code, void_code=void_code)
    width = max(first.seqlen, second.seqlen)

    def pad(data: np.ndarray) -> np.ndarray:
        if data.shape[1] == width:
            return data
        padding = np.full((len(data), width - data.shape[1]), void_code, dtype=int)
        return np.column_stack([data, padding])

    if (first.weights is None) != (second.weights is None):
        raise ValueError("cannot combine STSList objects with and without weights")
    weights = (
        None
        if first.weights is None
        else np.concatenate(
            [np.asarray(first.weights, dtype=float), np.asarray(second.weights, dtype=float)]
        )
    )
    return STSList(
        data=np.vstack([pad(first_data), pad(second_data)]),
        alphabet=alphabet,
        labels=None if provider.labels is None else list(provider.labels),
        cpal=None if provider.cpal is None else list(provider.cpal),
        xtlab=None,
        xtstep=first.xtstep,
        tick_last=first.tick_last,
        nr=nr,
        void=void,
        nr_code=nr_code,
        void_code=void_code,
        weights=weights,
        start=first.start,
        ids=None,
        missing_color=first.missing_color,
        issues=None,
    )


@dataclass(frozen=True)
class _SequencePair:
    label: str
    sequences: STSList
    first_size: int
    second_size: int


def _sequence_pair(label: str, first: STSList, second: STSList) -> _SequencePair:
    if first.nseq < second.nseq:
        first, second = second, first
    return _SequencePair(
        label=label,
        sequences=_combine_sequences(first, second),
        first_size=first.nseq,
        second_size=second.nseq,
    )


def _sequence_pairs(
    seqdata: STSList | Sequence[STSList],
    seqdata2: STSList | Sequence[STSList] | None,
    group: Sequence[object] | None,
    set_values: Sequence[object] | None,
) -> list[_SequencePair]:
    if isinstance(seqdata, STSList):
        if seqdata2 is not None:
            if group is not None:
                raise ValueError("group and seqdata2 cannot both be provided")
            if set_values is not None:
                raise ValueError("set and seqdata2 cannot both be provided")
            if not isinstance(seqdata2, STSList):
                raise TypeError("seqdata2 must be an STSList")
            return [_sequence_pair("1", seqdata, seqdata2)]
        if group is None:
            if set_values is not None:
                raise ValueError("set cannot be provided when group is None")
            raise ValueError("seqdata2 and group cannot both be None")

        groups = np.asarray(group, dtype=object)
        if groups.ndim != 1 or groups.shape[0] != seqdata.nseq:
            raise ValueError("group length must equal the number of sequences")
        if set_values is None:
            sets = np.ones(seqdata.nseq, dtype=int)
            include_set_labels = False
        else:
            sets = np.asarray(set_values, dtype=object)
            if sets.ndim != 1 or sets.shape[0] != seqdata.nseq:
                raise ValueError("set length must equal the number of sequences")
            include_set_labels = True
        valid = ~pd.isna(groups) & ~pd.isna(sets)
        removed = int(seqdata.nseq - np.count_nonzero(valid))
        if removed:
            warnings.warn(
                f"[!!] {removed} sequences removed because of NA values of the grouping variable(s)",
                UserWarning,
                stacklevel=3,
            )
        group_labels, group_codes = _compare_factor(groups, valid)
        if len(group_labels) < 2:
            raise ValueError("There is only one group among valid cases!")
        if len(group_labels) > 2:
            raise ValueError("Currently seqLRT supports only 2 groups!")
        set_labels, set_codes = _compare_factor(sets, valid)
        pairs: list[_SequencePair] = []
        for set_code, set_label in enumerate(set_labels):
            in_set = set_codes == set_code
            first_rows = np.flatnonzero(in_set & (group_codes == 0))
            second_rows = np.flatnonzero(in_set & (group_codes == 1))
            label = set_label if include_set_labels else "1"
            pairs.append(
                _sequence_pair(
                    label,
                    _subset_sequence(seqdata, first_rows),
                    _subset_sequence(seqdata, second_rows),
                )
            )
        return pairs

    if group is not None or set_values is not None:
        raise ValueError("group and set can only be used with a single STSList")
    if isinstance(seqdata, (str, bytes)) or not isinstance(seqdata, Sequence):
        raise TypeError("seqdata must be an STSList or a sequence of STSList objects")
    if (
        seqdata2 is None
        or isinstance(seqdata2, (STSList, str, bytes))
        or not isinstance(seqdata2, Sequence)
    ):
        raise TypeError("seqdata2 must be a sequence of STSList objects")
    first_list = list(seqdata)
    second_list = list(seqdata2)
    if len(first_list) != len(second_list):
        raise ValueError("seqdata and seqdata2 lists must have the same length")
    if not first_list:
        raise ValueError("sequence comparison lists must not be empty")
    if any(not isinstance(value, STSList) for value in [*first_list, *second_list]):
        raise TypeError("at least one list element is not an STSList")
    return [
        _sequence_pair(str(index), first, second)
        for index, (first, second) in enumerate(zip(first_list, second_list, strict=True), start=1)
    ]


def seqCompare(
    seqdata: STSList | Sequence[STSList],
    seqdata2: STSList | Sequence[STSList] | None = None,
    *,
    group: Sequence[object] | None = None,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    stat: Literal["all", "LRT", "BIC"] = "all",
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    opt: Literal[1, 2] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
    method: str = "LCS",
    **distance_kwargs: object,
) -> pd.DataFrame:
    """Compare two groups of sequences through the shared comparison engine.

    The R default sampling route and explicit ``s=0`` extension are supported.
    ``opt=1`` builds one public :func:`seqdist` matrix per balanced draw;
    ``opt=2`` builds one per pair of sets. Both feed the same aggregation used
    by :func:`dissCompare`.
    """
    if isinstance(s, (bool, np.bool_)) or not isinstance(s, (int, np.integer)):
        raise TypeError("s must be a nonnegative integer")
    sample_size = int(s)
    if sample_size < 0:
        raise ValueError("s must be a nonnegative integer")
    if opt not in {None, 1, 2}:
        raise ValueError("opt must be None, 1, or 2")
    if BFopt not in {None, 1, 2}:
        raise ValueError("BFopt must be None, 1, or 2")
    if not isinstance(inc_df, (bool, np.bool_)):
        raise TypeError("inc_df must be a bool")
    if not isinstance(method, str) or not method:
        raise ValueError("method must be a nonempty string")
    incompatible = {"condensed", "full_matrix", "refseq", "ref_sets"}.intersection(distance_kwargs)
    if incompatible:
        raise ValueError(
            "distance kwargs cannot override square-matrix construction: "
            + ", ".join(sorted(incompatible))
        )

    sq, lrt_pow = _resolve_squared(squared)
    is_lrt, is_bic, _ = _stat_columns(stat)
    if isinstance(weighted, (bool, np.bool_)):
        use_weights = bool(weighted)
        weight_by = "global"
    elif weighted == "by.group":
        use_weights = True
        weight_by = "by.group"
    else:
        raise ValueError("weighted must be logical or 'by.group'")

    pairs = _sequence_pairs(seqdata, seqdata2, group, set)
    compare_sets = [
        _CompareSet(
            pair.label,
            np.arange(pair.first_size, dtype=int),
            np.arange(pair.second_size, dtype=int) + pair.first_size,
        )
        for pair in pairs
    ]
    plans = _sampling_plans(compare_sets, sample_size, seed)
    effective_options = [
        (
            2
            if sample_size == 0
            else (opt if opt is not None else (1 if pair.sequences.nseq > 2 * sample_size else 2))
        )
        for pair in pairs
    ]

    def distance_matrix(sequences: STSList) -> np.ndarray:
        matrix = np.asarray(
            seqdist(
                sequences,
                method=method,
                weighted=use_weights,
                **distance_kwargs,
            ),
            dtype=float,
        )
        if matrix.shape != (sequences.nseq, sequences.nseq):
            raise ValueError("seqdist must return a square matrix for seqCompare")
        return matrix

    full_distances: dict[int, np.ndarray] = {}
    for index, pair in enumerate(pairs):
        if effective_options[index] == 2 and min(pair.first_size, pair.second_size) > 0:
            full_distances[index] = distance_matrix(pair.sequences)

    def evaluate(index: int, first: np.ndarray, second: np.ndarray) -> np.ndarray:
        pair = pairs[index]
        if effective_options[index] == 2:
            matrix = full_distances[index]
            weights = pair.sequences.weights
            first_rows = first
            second_rows = second
        else:
            selected = np.concatenate([first, second])
            sampled = _subset_sequence(pair.sequences, selected)
            matrix = distance_matrix(sampled)
            weights = sampled.weights
            first_rows = np.arange(len(first), dtype=int)
            second_rows = np.arange(len(second), dtype=int) + len(first)
        return _seqxcomp(
            first_rows,
            second_rows,
            matrix,
            None if weights is None else np.asarray(weights, dtype=float),
            is_lrt=is_lrt,
            is_bic=is_bic,
            squared=sq,
            weighted=use_weights,
            weight_by=weight_by,
            lrt_pow=lrt_pow,
        )

    return _aggregate_comparisons(
        compare_sets,
        plans,
        evaluate,
        is_lrt=is_lrt,
        is_bic=is_bic,
        BFopt=BFopt,
        inc_df=bool(inc_df),
    )


def seqLRT(
    seqdata: STSList | Sequence[STSList],
    seqdata2: STSList | Sequence[STSList] | None = None,
    *,
    group: Sequence[object] | None = None,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    opt: Literal[1, 2] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
    method: str = "LCS",
    **distance_kwargs: object,
) -> pd.DataFrame:
    """LRT of a two-group split of sequences — ``seqCompare(stat="LRT")``."""
    return seqCompare(
        seqdata,
        seqdata2,
        group=group,
        set=set,
        s=s,
        seed=seed,
        stat="LRT",
        squared=squared,
        weighted=weighted,
        opt=opt,
        BFopt=BFopt,
        inc_df=inc_df,
        method=method,
        **distance_kwargs,
    )


def seqBIC(
    seqdata: STSList | Sequence[STSList],
    seqdata2: STSList | Sequence[STSList] | None = None,
    *,
    group: Sequence[object] | None = None,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    opt: Literal[1, 2] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
    method: str = "LCS",
    **distance_kwargs: object,
) -> pd.DataFrame:
    """BIC / Bayes factor of a two-group split — ``seqCompare(stat="BIC")``."""
    return seqCompare(
        seqdata,
        seqdata2,
        group=group,
        set=set,
        s=s,
        seed=seed,
        stat="BIC",
        squared=squared,
        weighted=weighted,
        opt=opt,
        BFopt=BFopt,
        inc_df=inc_df,
        method=method,
        **distance_kwargs,
    )
