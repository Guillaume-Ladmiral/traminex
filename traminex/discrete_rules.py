"""Discrete person-period data for sequential association-rule models."""

from __future__ import annotations

import sys
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Real
from typing import Literal

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import norm

from ._rformat import r_as_character, r_round
from ._rlinpack import dqrdc2_pivot
from .events import SubSeqElist, seqeapplysub


@dataclass(frozen=True)
class _BinomialFit:
    coefficients: np.ndarray
    pvalues: np.ndarray


def _binomial_mean_and_derivative(
    eta: np.ndarray, link: Literal["cloglog", "logit"]
) -> tuple[np.ndarray, np.ndarray]:
    epsilon = np.finfo(float).eps
    if link == "logit":
        # R's logit_linkinv and logit_mu_eta (stats' C code): outside [-30, 30] the mean is
        # eps/(1+eps) or 1/(1+eps) and the derivative eps, which set the working weights of
        # separated rows (TX-REP-20261022-004).
        exp_eta = np.exp(np.clip(eta, -30.0, 30.0))
        tmp = np.where(eta < -30.0, epsilon, np.where(eta > 30.0, 1.0 / epsilon, exp_eta))
        opexp = 1.0 + exp_eta
        derivative = np.where(np.abs(eta) > 30.0, epsilon, exp_eta / (opexp * opexp))
        return tmp / (1.0 + tmp), derivative
    exp_eta = np.exp(np.minimum(eta, 700.0))
    mean = -np.expm1(-exp_eta)
    derivative = np.exp(np.minimum(eta - exp_eta, 700.0))
    mean = np.clip(mean, epsilon, 1.0 - epsilon)
    derivative = np.maximum(derivative, epsilon)
    return mean, derivative


def _binomial_deviance(y: np.ndarray, mean: np.ndarray) -> float:
    return float(
        -2.0
        * np.sum(
            np.where(y == 1.0, np.log(mean), np.log1p(-mean)),
            dtype=float,
        )
    )


# glm.fit's rank tolerance: Cdqrls with min(1e-7, epsilon / 1000), glm.control's epsilon 1e-8.
_GLM_RANK_TOLERANCE = 1e-11


def _binomial_irls(
    y: np.ndarray,
    design: np.ndarray,
    *,
    link: Literal["cloglog", "logit"],
) -> _BinomialFit:
    """Fit R ``glm(..., family=binomial)``-compatible models with IRLS.

    As ``glm.fit``, each step solves the weighted least squares on the columns R's
    dqrdc2 keeps at ``tol = 1e-11`` in the design weighted by that step's working
    weights; a column it sets aside gets coefficient 0 in the step, and one set
    aside in the last step is aliased: NaN (R's NA) coefficient and p-value
    (TX-AUD-20261011-010, TX-REP-20261019-002). After 25 steps without
    convergence it warns "glm.fit: algorithm did not converge" and keeps the
    last iterate, as glm.fit does (TX-AUD-20261019-085).
    """
    response = np.asarray(y, dtype=float)
    matrix = np.asarray(design, dtype=float)
    if response.ndim != 1 or matrix.ndim != 2 or len(response) != len(matrix):
        raise ValueError("y and design must have compatible one/two-dimensional shapes")
    if len(response) == 0 or matrix.shape[1] == 0:
        raise ValueError("binomial model inputs must not be empty")
    if not np.isfinite(response).all() or not np.isfinite(matrix).all():
        raise ValueError("binomial model inputs must be finite")
    if not np.isin(response, [0.0, 1.0]).all():
        raise ValueError("binomial response must contain only 0 or 1")
    if link not in ("cloglog", "logit"):
        raise ValueError("link must be 'cloglog' or 'logit'")

    mean = (response + 0.5) / 2.0
    if link == "logit":
        eta = np.log(mean / (1.0 - mean))
    else:
        eta = np.log(-np.log1p(-mean))
    old_deviance = _binomial_deviance(response, mean)
    coefficients = np.zeros(matrix.shape[1], dtype=float)
    weights = np.ones(len(response), dtype=float)

    for _ in range(25):
        mean, derivative = _binomial_mean_and_derivative(eta, link)
        variance = mean * (1.0 - mean)
        weights = derivative * derivative / variance
        working_response = eta + (response - mean) / derivative
        sqrt_weights = np.sqrt(weights)
        weighted_design = matrix * sqrt_weights[:, None]
        weighted_response = working_response * sqrt_weights
        rank, pivot = dqrdc2_pivot(weighted_design, _GLM_RANK_TOLERANCE)
        kept = pivot[:rank]
        coefficients = np.zeros(matrix.shape[1], dtype=float)
        r_matrix = np.zeros((0, 0))
        if rank:
            q_matrix, r_matrix = np.linalg.qr(weighted_design[:, kept])
            coefficients[kept] = solve_triangular(r_matrix, q_matrix.T @ weighted_response)
        eta = matrix @ coefficients
        mean, _ = _binomial_mean_and_derivative(eta, link)
        deviance = _binomial_deviance(response, mean)
        if abs(deviance - old_deviance) / (0.1 + abs(deviance)) < 1e-8:
            break
        old_deviance = deviance
    else:
        # glm.fit stops after maxit = 25 iterations with R's warning and keeps the last iterate
        # (TX-AUD-20261019-085).
        warnings.warn("glm.fit: algorithm did not converge", UserWarning, stacklevel=2)
    # glm.fit's check on the fitted means, after the convergence warning (TX-REP-20261022-005).
    tiny = 10.0 * np.finfo(float).eps
    if np.any(mean > 1.0 - tiny) or np.any(mean < tiny):
        warnings.warn(
            "glm.fit: fitted probabilities numerically 0 or 1 occurred", UserWarning, stacklevel=2
        )

    # R's ``glm.fit`` retains the weighted QR decomposition from the final
    # IRLS solve. Those weights precede the just-updated mean by one half-step;
    # recomputing them here produces subtly different standard errors. As
    # summary.glm (chol2inv of that QR's R), the covariance is inv(R) inv(R)';
    # inverting X'WX squares the condition number of a near-collinear design.
    r_inverse = solve_triangular(r_matrix, np.eye(rank))
    standard_errors = np.sqrt(np.sum(r_inverse * r_inverse, axis=1))
    estimates = np.full(matrix.shape[1], np.nan)
    pvalues = np.full(matrix.shape[1], np.nan)
    estimates[kept] = coefficients[kept]
    pvalues[kept] = 2.0 * norm.sf(np.abs(coefficients[kept] / standard_errors))
    return _BinomialFit(coefficients=estimates, pvalues=pvalues)


def _required_columns(frame: pd.DataFrame, name: str, required: Sequence[str]) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {', '.join(missing)}")


def _unique_id_rows(frame: pd.DataFrame, name: str) -> dict[object, pd.Series]:
    if frame["IDPERS"].isna().any():
        raise ValueError(f"{name} IDPERS must not contain missing values")
    if frame["IDPERS"].duplicated().any():
        raise ValueError(f"{name} IDPERS must be unique")
    return {row["IDPERS"]: row for _, row in frame.iterrows()}


def _age(value: object, label: str) -> float:
    """A finite numeric age; R's createdatadiscrete accepts fractional ages too."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{label} must be a numeric age")
    numeric = float(value)
    if not np.isfinite(numeric):
        raise ValueError(f"{label} must be a numeric age")
    return numeric


def _supplementary_value(value: object, label: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        return float(int(value))
    if not isinstance(value, Real) or not np.isfinite(float(value)):
        raise ValueError(f"supplementary variable {label} must be numeric and nonmissing")
    return float(int(float(value)))


def _person_period_frame(
    ids: list[object],
    data_rows: dict[object, pd.Series],
    variables: list[str],
    dependent: str,
    min_rows: dict[object, pd.Series],
    max_rows: dict[object, pd.Series],
    supplementary: list[str],
) -> pd.DataFrame:
    independent = [variable for variable in variables if variable != dependent]
    columns = ["IDPERS", "AGE", dependent, *independent, *supplementary]
    rows: list[list[object]] = []

    for position, identifier in enumerate(ids, start=1):
        if position % 1000 == 0:
            percent = round(position / len(ids) * 100, 2)
            print(f"{percent:g}% completed", file=sys.stderr)
        if identifier not in data_rows or identifier not in min_rows or identifier not in max_rows:
            continue
        source = data_rows[identifier]
        min_value = min_rows[identifier]["AGE"]
        max_value = max_rows[identifier]["AGE"]
        if pd.isna(min_value) or pd.isna(max_value):
            continue

        start = _age(min_value, f"agemin AGE for IDPERS {identifier}")
        end = _age(max_value, f"agemax AGE for IDPERS {identifier}")
        dependent_status = int(source[f"{dependent}ST"])
        if dependent_status == 1:
            end = _age(source[dependent], f"{dependent} for IDPERS {identifier}")
        if end < start:
            raise ValueError(
                f"observation end precedes start for IDPERS {identifier} and {dependent}"
            )

        # R: nl <- as.integer((agemaxi - agemini) + 1) rows, AGE agemini:agemaxi
        # (createdatadiscrete.R:50, 67), fractional ages included (TX-AUD-20261019-015).
        ages = [start + offset for offset in range(int(np.trunc(end - start + 1)))]
        identifier_value: object
        if isinstance(identifier, Real) and not isinstance(identifier, (bool, np.bool_)):
            identifier_value = float(identifier)
        else:
            identifier_value = identifier
        supplementary_values = [
            _supplementary_value(source[column], column) for column in supplementary
        ]
        for offset, age in enumerate(ages):
            indicators: list[float] = []
            for variable in independent:
                status = int(source[f"{variable}ST"])
                if status == 1:
                    event_age = _age(source[variable], f"{variable} for IDPERS {identifier}")
                    # R sets the rows from as.integer((age - agemini) + 1) on
                    # (createdatadiscrete.R:71-73); an index of 0 or less marks every row.
                    first_row = int(np.trunc(event_age - start + 1))
                    indicator = float(event_age <= end and offset + 1 >= first_row)
                else:
                    indicator = 0.0
                indicators.append(indicator)
            dependent_value = float(dependent_status == 1 and offset == len(ages) - 1)
            rows.append(
                [
                    identifier_value,
                    float(age),
                    dependent_value,
                    *indicators,
                    *supplementary_values,
                ]
            )

    frame = pd.DataFrame(rows, columns=columns)
    numeric_ids = all(
        isinstance(identifier, Real) and not isinstance(identifier, (bool, np.bool_))
        for identifier in ids
    )
    numeric_columns = columns if numeric_ids else columns[1:]
    for column in numeric_columns:
        frame[column] = frame[column].astype(float)
    return frame


def createdatadiscrete(
    ids: Sequence[object],
    data: pd.DataFrame,
    vars: Sequence[str],
    agemin: pd.DataFrame,
    agemax: pd.DataFrame,
    supvar: Sequence[str] | None = None,
) -> dict[str, object]:
    """Transform event times into one person-period frame per dependent event.

    The top-level and per-event insertion order follows TraMineRextras. Two
    unsafe upstream defects are corrected: one/multiple-variable inputs work,
    and IDs without complete observation windows are skipped instead of
    reusing another person's rows.

    Ages may be fractional, as in R: a person has ``trunc(end - start + 1)`` rows
    at ages ``start, start + 1, ...``, the dependent event sits on the last row,
    and an independent event at age ``a`` is set from row ``trunc(a - start + 1)``
    on (so an event at 21 with ``start = 19.5`` is set from the row of age 20.5).
    An event before ``start`` is refused (R stops on an index error).
    """
    if not all(isinstance(frame, pd.DataFrame) for frame in (data, agemin, agemax)):
        raise TypeError("data, agemin, and agemax must be pandas DataFrames")
    id_values = list(ids)
    if pd.Index(id_values).hasnans:
        raise ValueError("ids must not contain missing values")
    if pd.Index(id_values).has_duplicates:
        raise ValueError("ids must be unique")

    variables = [vars] if isinstance(vars, str) else list(vars)
    if not variables:
        raise ValueError("vars must contain at least one event variable")
    if any(not isinstance(variable, str) or not variable for variable in variables):
        raise ValueError("vars must contain nonempty strings")
    if len(set(variables)) != len(variables):
        raise ValueError("vars must be unique")
    supplementary = (
        [] if supvar is None else ([supvar] if isinstance(supvar, str) else list(supvar))
    )
    if any(not isinstance(variable, str) or not variable for variable in supplementary):
        raise ValueError("supvar must contain nonempty strings")
    if len(set(supplementary)) != len(supplementary):
        raise ValueError("supvar must be unique")
    reserved = {"IDPERS", "AGE", *variables, *(f"{item}ST" for item in variables)}
    overlap = [variable for variable in supplementary if variable in reserved]
    if overlap:
        raise ValueError(f"supvar overlaps reserved columns: {', '.join(overlap)}")

    _required_columns(
        data,
        "data",
        [
            "IDPERS",
            *variables,
            *(f"{variable}ST" for variable in variables),
            *supplementary,
        ],
    )
    _required_columns(agemin, "agemin", ["IDPERS", "AGE"])
    _required_columns(agemax, "agemax", ["IDPERS", "AGE"])
    data_rows = _unique_id_rows(data, "data")
    min_rows = _unique_id_rows(agemin, "agemin")
    max_rows = _unique_id_rows(agemax, "agemax")

    for variable in variables:
        status_column = f"{variable}ST"
        invalid = data[status_column].isna() | ~data[status_column].isin([0, 1])
        if invalid.any():
            raise ValueError(f"{status_column} must contain only 0 or 1")

    event_frames: dict[str, pd.DataFrame] = {}
    for dependent in variables:
        print(f"Creating data set for variable {dependent}", file=sys.stderr)
        event_frames[dependent] = _person_period_frame(
            id_values,
            data_rows,
            variables,
            dependent,
            min_rows,
            max_rows,
            supplementary,
        )
    return {"data": event_frames, "vars": [], "ids": id_values}


def _rule_events(label: str) -> list[str]:
    sequential = str(label).replace(",", ")-(")
    return [part.strip().strip("()") for part in sequential.split("-")]


def _r_number(value: float, digits: int) -> str:
    """R's ``paste(round(value, digits))`` (seqerulesdisc.R:89): fixed or scientific as R
    prints it, so 0.0002 reads ``2e-04`` (TX-AUD-20261019-043), rounded by R's ``round``
    (``r_round``: 6.165 gives 6.16, batch-11 review M-4); NaN reads ``NA``."""
    text = r_as_character(r_round(float(value), digits))
    return "NA" if text is None else text


def _empty_rules_frame() -> pd.DataFrame:
    result = pd.DataFrame()
    result.attrs["graph"] = {"events": [], "hr": [], "hrindep": [], "pvalues": []}
    return result


def _validate_rule_data(
    datadiscr: Mapping[str, object], tsef: pd.DataFrame
) -> Mapping[str, object]:
    if not isinstance(datadiscr, Mapping) or "data" not in datadiscr:
        raise ValueError("datadiscr must be a mapping with a 'data' mapping")
    event_frames = datadiscr["data"]
    if not isinstance(event_frames, Mapping):
        raise ValueError("datadiscr['data'] must be a mapping of event frames")
    if not isinstance(tsef, pd.DataFrame):
        raise TypeError("tsef must be a pandas DataFrame")
    if "IDPERS" not in tsef.columns:
        raise ValueError("tsef is missing required column: IDPERS")
    if tsef["IDPERS"].isna().any() or tsef["IDPERS"].duplicated().any():
        raise ValueError("tsef IDPERS must be unique and nonmissing")
    return event_frames


def _model_frame(
    frame: pd.DataFrame,
    tsef: pd.DataFrame,
    dependent: str,
    independent: str,
    supplementary: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    _required_columns(
        frame, f"datadiscr data for {dependent}", ["IDPERS", "AGE", dependent, independent]
    )
    model_frame = frame.copy()
    source = tsef.set_index("IDPERS")
    for variable in supplementary:
        model_frame[variable] = model_frame["IDPERS"].map(source[variable])
    model_columns = [independent, "AGE", *supplementary]
    numeric = model_frame[[dependent, *model_columns]].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError("discrete rule model columns must be numeric, finite, and nonmissing")
    response = numeric[dependent].to_numpy(dtype=float)
    design = np.column_stack(
        [np.ones(len(numeric), dtype=float), numeric[model_columns].to_numpy(dtype=float)]
    )
    return response, design


def seqerulesdisc(
    fsubseq: SubSeqElist,
    datadiscr: Mapping[str, object],
    tsef: pd.DataFrame,
    *,
    pvalue: float = 0.1,
    supvars: Sequence[str] | None = None,
    adjust: bool = True,
    topt: bool = False,
    link: Literal["cloglog", "logit"] = "cloglog",
    dep: str | None = None,
) -> pd.DataFrame:
    """Build discrete-time hazard rules from frequent event subsequences.

    This ports TraMineRextras ``seqerulesdisc`` while repairing its unusable
    supplementary-variable branch: controls are included as ordinary covariates
    in the same binomial GLM. The R function's final ``print("END")`` is display
    behavior and is intentionally not emitted by this programmatic API.
    """
    if not isinstance(fsubseq, SubSeqElist):
        raise TypeError("fsubseq must be a SubSeqElist")
    event_frames = _validate_rule_data(datadiscr, tsef)
    if isinstance(pvalue, (bool, np.bool_)) or not isinstance(pvalue, Real):
        raise ValueError("pvalue must be a number from 0 to 1")
    threshold = float(pvalue)
    if not np.isfinite(threshold) or threshold < 0.0 or threshold > 1.0:
        raise ValueError("pvalue must be a number from 0 to 1")
    if link not in ("cloglog", "logit"):
        raise ValueError("link must be 'cloglog' or 'logit'")
    if dep is not None and (not isinstance(dep, str) or not dep):
        raise ValueError("dep must be a nonempty event name or None")

    supplementary = (
        [] if supvars is None else ([supvars] if isinstance(supvars, str) else list(supvars))
    )
    if any(not isinstance(variable, str) or not variable for variable in supplementary):
        raise ValueError("supvars must contain nonempty strings")
    if len(set(supplementary)) != len(supplementary):
        raise ValueError("supvars must be unique")
    reserved_controls = {"IDPERS", "AGE"}
    if reserved_controls.intersection(supplementary):
        raise ValueError("supvars must not contain IDPERS or AGE")
    missing_controls = [variable for variable in supplementary if variable not in tsef.columns]
    if missing_controls:
        raise ValueError(f"supvars are missing from tsef: {', '.join(missing_controls)}")

    subsequences = list(fsubseq.subseq)
    if topt and subsequences:
        counts = seqeapplysub(fsubseq, method="presence", rules=True)
        subsequences = [
            label
            for position, label in enumerate(subsequences)
            if float(counts.iloc[:, position].sum()) == 1.0
        ]
    if not subsequences:
        return _empty_rules_frame()
    if adjust:
        threshold /= len(subsequences)
        print(
            f"Adjusted p-value threshold : {threshold:.15g}",
            file=sys.stderr,
        )

    hazards: list[str] = []
    graph_events: list[list[str]] = []
    graph_hr: list[list[float]] = []
    graph_hrindep: list[list[list[float]]] = []
    graph_pvalues: list[list[float]] = []

    for label in subsequences:
        events = _rule_events(label)
        if dep is not None and (not events or events[-1] != dep):
            continue
        if len(events) <= 1:
            continue

        rule_hr: list[float] = []
        rule_controls: list[list[float]] = []
        rule_pvalues: list[float] = []
        hazard = ""

        for edge in range(1, len(events)):
            dependent = events[edge]
            independent = events[edge - 1]
            if dependent not in event_frames or not isinstance(
                event_frames[dependent], pd.DataFrame
            ):
                raise ValueError(f"datadiscr is missing a DataFrame for event {dependent}")
            frame = event_frames[dependent].copy()

            if edge > 1:
                for earlier in range(edge - 1):
                    first = events[earlier]
                    second = events[earlier + 1]
                    _required_columns(tsef, "tsef", [first, second, f"{first}ST"])
                    ordered = (
                        pd.to_numeric(tsef[first], errors="coerce")
                        < pd.to_numeric(tsef[second], errors="coerce")
                    ) & (tsef[f"{first}ST"] == 1)
                    valid_ids = set(tsef.loc[ordered, "IDPERS"])
                    _required_columns(frame, f"datadiscr data for {dependent}", ["IDPERS", first])
                    frame = frame.loc[frame["IDPERS"].isin(valid_ids) & (frame[first] == 1)].copy()

            _required_columns(frame, f"datadiscr data for {dependent}", [dependent])
            outcome = pd.to_numeric(frame[dependent], errors="coerce")
            if outcome.isna().any():
                raise ValueError(f"datadiscr outcome {dependent} must be numeric and nonmissing")
            if float(outcome.sum()) > 1.0:
                response, design = _model_frame(frame, tsef, dependent, independent, supplementary)
                fit = _binomial_irls(response, design, link=link)
                if np.isnan(fit.coefficients[1]):
                    # R reports the next coefficient (AGE's) under the predictor's
                    # name (seqerulesdisc.R:81-87); without an estimate the p-value is
                    # missing and R's own filter leaves the rule out (:96-98).
                    warnings.warn(
                        f"seqerulesdisc: the predictor {independent} of rule {label} is aliased "
                        "(constant in its model data); the rule has no estimate and is left out",
                        UserWarning,
                        stacklevel=2,
                    )
                hazard_ratio = float(np.exp(fit.coefficients[1]))
                event_pvalue = float(fit.pvalues[1])
                rule_hr.append(hazard_ratio)
                rule_pvalues.append(event_pvalue)
                if supplementary:
                    rule_controls.append([float(value) for value in np.exp(fit.coefficients[3:])])
                hazard += (
                    f"{independent}(HR={_r_number(hazard_ratio, 2)},"
                    f"p={_r_number(event_pvalue, 4)}) =>"
                )
            else:
                hazard = "No occurences of "

        hazard = f"{hazard} {events[-1]}"
        missing_pvalue = any(np.isnan(value) for value in rule_pvalues)
        if missing_pvalue or any(value > threshold for value in rule_pvalues) or hazard in hazards:
            continue
        hazards.append(hazard)
        graph_events.append(events)
        if rule_hr:
            graph_hr.append(rule_hr)
            graph_pvalues.append(rule_pvalues)
            if supplementary:
                graph_hrindep.append(rule_controls)

    result = pd.DataFrame({"Hazards": hazards}) if hazards else pd.DataFrame()
    result.attrs["graph"] = {
        "events": graph_events,
        "hr": graph_hr,
        "hrindep": graph_hrindep,
        "pvalues": graph_pvalues,
    }
    return result
