"""Cross-domain association for categorical state sequences."""

from __future__ import annotations

from itertools import combinations
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .sequences import STSList
from .sequences import alphabet as sequence_alphabet


def _independence_lrt_pearson(
    xtab: np.ndarray,
    struct_zero: bool,
) -> tuple[float, float, int]:
    """Return likelihood-ratio, Pearson statistic, and adjusted degrees of freedom."""
    observed = np.asarray(xtab, dtype=float)
    total = float(observed.sum())
    rows, columns = observed.shape
    base_df = max((rows - 1) * (columns - 1), 0)
    zero_count = int(np.count_nonzero(observed == 0)) if struct_zero else 0
    degrees = base_df - zero_count if zero_count < base_df else base_df
    if total <= 0 or rows == 0 or columns == 0:
        return float("nan"), float("nan"), degrees

    # Reproduce stats::loglin's numerical path as well as its model: initialize
    # admissible cells at one, then fit the one-way margins by proportional
    # adjustment.  Running this even without zeros preserves loglin's floating
    # arithmetic rather than substituting the algebraic outer-product shortcut.
    expected = (observed != 0).astype(float) if struct_zero else np.ones_like(observed, dtype=float)
    row_targets = observed.sum(axis=1)
    column_targets = observed.sum(axis=0)
    final_adjustment = False
    for _ in range(20):
        row_sums = expected.sum(axis=1)
        expected *= np.divide(
            row_targets,
            row_sums,
            out=np.zeros_like(row_targets),
            where=row_sums > 0,
        )[:, None]
        column_sums = expected.sum(axis=0)
        expected *= np.divide(
            column_targets,
            column_sums,
            out=np.zeros_like(column_targets),
            where=column_sums > 0,
        )[None, :]
        if final_adjustment:
            break
        row_deviation = np.max(np.abs(expected.sum(axis=1) - row_targets), initial=0.0)
        column_deviation = np.max(np.abs(expected.sum(axis=0) - column_targets), initial=0.0)
        # stats::loglin uses eps=0.1 and its C routine performs the next
        # margin-adjustment sweep after the convergence deviation is met.
        final_adjustment = max(row_deviation, column_deviation) < 0.1
    valid_expected = expected > 0
    pearson = float(
        np.sum(
            np.divide(
                (observed - expected) ** 2,
                expected,
                out=np.zeros_like(expected),
                where=valid_expected,
            )
        )
    )
    positive = (observed > 0) & valid_expected
    likelihood = float(
        2.0 * np.sum(observed[positive] * np.log(observed[positive] / expected[positive]))
    )
    return likelihood, pearson, degrees


def _decoded_panel(sequence: STSList, with_missing: bool) -> tuple[np.ndarray, list[str]]:
    """Decode an STSList into state labels and its cross-table level order."""
    if not isinstance(sequence, STSList):
        raise TypeError("seqdata_dom must contain STSList objects")
    sequence = sequence._canonicalized()
    labels = np.asarray(sequence.alphabet, dtype=object)
    if np.any(sequence.data < 0) or np.any(sequence.data >= labels.size):
        raise ValueError("sequence data contain codes outside the declared alphabet")
    decoded = labels[sequence.data]
    levels = list(sequence_alphabet(sequence, with_missing=False))
    # R keeps the missing row and column of every domain with with.missing,
    # missing cells or not (seqdomassoc.R:206-208; TX-AUD-20261011-013).
    if with_missing and sequence.nr not in levels:
        levels.append(sequence.nr)
    return decoded, levels


def _representative_panels(
    domains: Sequence[STSList],
    rep_method: str,
    diss_dom: Sequence[np.ndarray] | None,
    weights: np.ndarray,
    *,
    weighted: bool,
    seqrep_args: Mapping[str, object] | None,
    seqrf_args: Mapping[str, object] | None,
) -> tuple[list[np.ndarray], np.ndarray]:
    """Resolve overall, representative, or equal-frequency panels."""
    if rep_method == "overall":
        return [np.asarray(domain.data) for domain in domains], weights
    if diss_dom is None:
        raise ValueError("Except for rep_method='overall', diss_dom must be supplied")
    if len(diss_dom) != len(domains):
        raise ValueError("diss_dom should be a list of same length as seqdata_dom")

    panels: list[np.ndarray] = []
    if rep_method == "rep":
        from .cluster import seqrep

        # As R (seqdomassoc.R:79-80), a given mapping replaces the defaults; seqrep's
        # own defaults then apply to what it leaves out (TX-AUD-20261011-012).
        options: dict[str, object] = (
            {"coverage": 0.8, "pradius": 0.1} if seqrep_args is None else dict(seqrep_args)
        )
        for domain, distances in zip(domains, diss_dom):
            result = seqrep(
                domain,
                diss=np.asarray(distances, dtype=float),
                weighted=weighted,
                **options,  # type: ignore[arg-type]
            )
            distance_map = np.asarray(result["distance_matrix"], dtype=float)
            nearest = np.nanargmin(distance_map, axis=1)
            representatives = np.asarray(result["representatives"])
            panels.append(representatives[nearest])
        return panels, weights

    if rep_method == "eq.group":
        from .repseq_rf import seqrf

        options = {"k": 20} if seqrf_args is None else dict(seqrf_args)
        last_result: dict[str, Any] | None = None
        for domain, distances in zip(domains, diss_dom):
            last_result = seqrf(
                domain,
                np.asarray(distances, dtype=float),
                weighted=weighted,
                **options,  # type: ignore[arg-type]
            )
            panels.append(np.asarray(last_result["seqtoplot"]))
        assert last_result is not None
        heights = np.asarray(last_result["rf"]["heights"], dtype=float)
        return panels, heights * domains[0].nseq

    raise ValueError(f"Unknown rep_method: {rep_method}")


def _code_panels_as_labels(
    panels: Sequence[np.ndarray],
    domains: Sequence[STSList],
) -> list[np.ndarray]:
    decoded: list[np.ndarray] = []
    for panel, domain in zip(panels, domains):
        labels = np.asarray(domain.alphabet, dtype=object)
        codes = np.asarray(panel)
        if np.any(codes < 0) or np.any(codes >= labels.size):
            raise ValueError("representative sequence contains an unknown state code")
        decoded.append(labels[codes])
    return decoded


def _cross_table(
    left: np.ndarray,
    right: np.ndarray,
    left_levels: Sequence[str],
    right_levels: Sequence[str],
    weights: np.ndarray,
    wrange: tuple[int, int] | None,
) -> np.ndarray:
    """Accumulate a position-wise or windowed weighted contingency table."""
    left_lookup = {state: index for index, state in enumerate(left_levels)}
    right_lookup = {state: index for index, state in enumerate(right_levels)}
    table = np.zeros((len(left_levels), len(right_levels)), dtype=float)
    n_cases, width = left.shape
    if wrange is not None:
        # Corrected R defect (TX-AUD-20261011-009): R's seqxtabs.win windows the
        # column-major unlist of the domains in blocks of ncol/ncases, so it
        # counts states of different people, ignores the weights and uses the
        # alphabet only (seqdomassoc.R:170,210,277-301). Here the window runs
        # within each person's own sequence, weighted by the person's weight,
        # over the levels of the position-wise table.
        low, high = wrange
        for case in range(n_cases):
            case_weight = float(weights[case])
            for position in range(width):
                left_index = left_lookup.get(left[case, position])
                if left_index is None:
                    continue
                for other in range(max(0, position + low), min(width - 1, position + high) + 1):
                    right_index = right_lookup.get(right[case, other])
                    if right_index is not None:
                        table[left_index, right_index] += case_weight
        return table

    for case in range(n_cases):
        case_weight = float(weights[case])
        for position in range(width):
            left_index = left_lookup.get(left[case, position])
            if left_index is None:
                continue
            right_index = right_lookup.get(right[case, position])
            if right_index is not None:
                table[left_index, right_index] += case_weight
    return table


def seqdomassoc(
    seqdata_dom: Sequence[STSList],
    *,
    assoc: str | Sequence[str] = ("LRT", "V"),
    rep_method: str = "overall",
    diss_dom: Sequence[np.ndarray] | None = None,
    wrange: Sequence[float] | None = None,
    cross_table: bool = False,
    seqrep_args: Mapping[str, object] | None = None,
    seqrf_args: Mapping[str, object] | None = None,
    with_missing: bool = False,
    weighted: bool = True,
    struct_zero: bool = True,
    p_value: bool = True,
    dnames: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Compute pairwise LRT and Cramer's V associations between domains.

    With ``wrange = (low, high)`` the cross-table counts, for each person, the
    state at position p in the first domain against the states at p + low ..
    p + high in the second, weighted by the person's weight, over the levels of
    the position-wise table. This corrects R, whose seqxtabs.win windows the
    column-major unlist of all cases and so mixes people (TX-AUD-20261011-009).

    ``weighted=False`` also selects the representatives (``rep``) and the
    groups (``eq.group``) without weights. R means to pass unit weights to
    seqrep, but seqrep drops them, and its seqrf keeps the object weights
    (seqdomassoc.R:86, 121-122), so R's representatives stay weighted; this
    is a recorded correction (TX-AUD-20261011-011).
    """
    if not isinstance(seqdata_dom, (list, tuple)) or len(seqdata_dom) < 2:
        raise ValueError("seqdata_dom should be a list of at least two state sequence objects")
    domains = list(seqdata_dom)
    if any(not isinstance(domain, STSList) for domain in domains):
        raise TypeError("seqdata_dom must contain STSList objects")
    domains = [domain._canonicalized() for domain in domains]
    n_cases = domains[0].nseq
    if any(domain.nseq != n_cases for domain in domains):
        raise ValueError("all domains must have the same number of cases")
    if any(domain.seqlen != domains[0].seqlen for domain in domains):
        raise ValueError("all domains must have the same sequence width")

    # One measure may be given as a string, as R's character vector of length 1
    # (seqdomassoc.Rd; TX-AUD-20261010-062).
    requested = tuple(
        str(value).upper() for value in ([assoc] if isinstance(assoc, str) else assoc)
    )
    if not requested or any(value not in {"LRT", "V"} for value in requested):
        raise ValueError("bad assoc values; allowed values are LRT and V")
    if rep_method not in {"overall", "rep", "eq.group"}:
        raise ValueError(f"Unknown rep_method: {rep_method}")
    if rep_method != "overall":
        if diss_dom is None or not isinstance(diss_dom, (list, tuple)):
            raise ValueError("Except for rep_method='overall', diss_dom must be a list")
        if len(diss_dom) != len(domains):
            raise ValueError("diss_dom should be a list of same length as seqdata_dom")

    resolved_window: tuple[int, int] | None = None
    if wrange is not None and not np.all(np.asarray(wrange) == 0):
        window = np.asarray(wrange, dtype=int).reshape(-1)
        if window.size != 2 or window[0] > window[1]:
            raise ValueError("wrange[0] cannot be larger than wrange[1]")
        resolved_window = (int(window[0]), int(window[1]))

    names = list(dnames) if dnames is not None else [f"Dom{i}" for i in range(1, len(domains) + 1)]
    if len(names) != len(domains):
        raise ValueError("dnames must have one name per domain")
    raw_weights = domains[0].weights if weighted else None
    if raw_weights is None:
        case_weights = np.ones(n_cases, dtype=float)
    else:
        case_weights = np.asarray(raw_weights, dtype=float).copy()
        if case_weights.shape != (n_cases,):
            raise ValueError("sequence weights must have one value per case")

    levels: list[list[str]] = []
    for domain in domains:
        _, domain_levels = _decoded_panel(domain, with_missing)
        levels.append(domain_levels)
    code_panels, analysis_weights = _representative_panels(
        domains,
        rep_method,
        diss_dom,
        case_weights,
        weighted=weighted,
        seqrep_args=seqrep_args,
        seqrf_args=seqrf_args,
    )
    panels = _code_panels_as_labels(code_panels, domains)
    if any(panel.shape != panels[0].shape for panel in panels):
        raise ValueError(
            "all representative panels must have the same number of cases and positions"
        )
    if analysis_weights.shape != (panels[0].shape[0],):
        raise ValueError("analysis weights must have one value per representative case")

    rows: list[dict[str, float | int]] = []
    row_names: list[str] = []
    tables: dict[str, pd.DataFrame] = {}
    for left_index, right_index in combinations(range(len(domains)), 2):
        if set(names[left_index].split(" x ")) & set(names[right_index].split(" x ")):
            continue
        table_values = _cross_table(
            panels[left_index],
            panels[right_index],
            levels[left_index],
            levels[right_index],
            analysis_weights,
            resolved_window,
        )
        pair_name = f"{names[left_index]}_with_{names[right_index]}"
        table = pd.DataFrame(
            table_values,
            index=pd.Index(levels[left_index], name=names[left_index]),
            columns=pd.Index(levels[right_index], name=names[right_index]),
        )
        if cross_table:
            tables[pair_name] = table

        likelihood, pearson, degrees = _independence_lrt_pearson(table_values, struct_zero)
        result: dict[str, float | int] = {"df": int(degrees)}
        if "LRT" in requested:
            result["LRT"] = likelihood
            if p_value:
                result["p(LRT)"] = float(chi2.sf(likelihood, degrees))
        if "V" in requested:
            total = float(table_values.sum())
            dimension = min(table_values.shape) - 1
            result["v"] = (
                float(np.sqrt(pearson / (total * dimension)))
                if total > 0 and dimension > 0
                else float("nan")
            )
            if p_value:
                result["p(v)"] = float(chi2.sf(pearson, degrees))
        rows.append(result)
        row_names.append(pair_name)

    output = pd.DataFrame(rows, index=row_names)
    if "df" in output:
        output["df"] = output["df"].astype(np.int64)
    if cross_table:
        output.attrs["cross_tables"] = tables
    return output
