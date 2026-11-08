"""Euclidean coordinates for longitudinal categorical timelines (EMLT)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from scipy import linalg

from .sequences import STSList
from .sequences import alphabet as sequence_alphabet


@dataclass
class EMLTResult:
    states: List[str]
    period: int
    situations: List[str]
    sit_freq: np.ndarray
    disjunctive: np.ndarray
    sit_transrate: np.ndarray
    sit_profil: np.ndarray
    distance_matrix: np.ndarray
    benz_covariance: np.ndarray
    sdev: np.ndarray
    loadings: np.ndarray
    coordinates: np.ndarray
    correlation: np.ndarray
    a: float
    b: float
    weighted: bool


def _situation_layout(n_states: int, n_periods: int) -> tuple[np.ndarray, np.ndarray]:
    """Return time-major time and state indexes for all situations."""
    times = np.repeat(np.arange(n_periods, dtype=int), n_states)
    states = np.tile(np.arange(n_states, dtype=int), n_periods)
    return times, states


def _situation_frequencies(
    data: np.ndarray,
    sit_time: np.ndarray,
    sit_state: np.ndarray,
) -> np.ndarray:
    frequencies = np.zeros(sit_time.size, dtype=float)
    for situation, (time, state) in enumerate(zip(sit_time, sit_state)):
        frequencies[situation] = np.count_nonzero(data[:, time] == state)
    return frequencies


def _disjunctive(
    data: np.ndarray,
    sit_time: np.ndarray,
    sit_state: np.ndarray,
) -> np.ndarray:
    encoded = np.zeros((data.shape[0], sit_time.size), dtype=float)
    for situation, (time, state) in enumerate(zip(sit_time, sit_state)):
        encoded[:, situation] = data[:, time] == state
    return encoded


def _transition_rates(
    disj: np.ndarray,
    weights: Optional[np.ndarray],
) -> np.ndarray:
    weighted = disj if weights is None else np.asarray(weights, dtype=float)[:, None] * disj
    n_situations = disj.shape[1]
    transitions = np.zeros((n_situations, n_situations), dtype=float)
    for source in range(n_situations):
        frequency = float(np.sum(weighted[:, source]))
        if frequency > 0:
            transitions[source, source:] = (weighted[:, source] @ weighted[:, source:]) / frequency
        else:
            transitions[source, source:] = np.nan
    return transitions


def _situation_profiles(
    transrate: np.ndarray,
    sit_time: np.ndarray,
    a: float,
    b: float,
) -> np.ndarray:
    profiles = np.zeros_like(transrate, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        for source in range(transrate.shape[0]):
            for target in range(source, transrate.shape[1]):
                transition = transrate[source, target]
                if np.isnan(transition):
                    profiles[source, target] = np.nan
                else:
                    delay = sit_time[target] - sit_time[source]
                    profiles[source, target] = transition / (a * delay + b)
            profiles[source] /= np.sum(profiles[source])
    return profiles


def _profile_distance(profiles: np.ndarray, active: np.ndarray) -> np.ndarray:
    selected = profiles[np.ix_(active, active)]
    column_sums = np.sum(selected, axis=0)
    compact = np.zeros((active.size, active.size), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        for left in range(active.size):
            for right in range(left + 1):
                delta = selected[left] - selected[right]
                value = float(np.dot(delta / column_sums, delta))
                compact[left, right] = value
                compact[right, left] = value
    full = np.full((profiles.shape[0], profiles.shape[0]), np.nan, dtype=float)
    full[np.ix_(active, active)] = compact
    return full


def _benzecri(dsub: np.ndarray) -> np.ndarray:
    row_means = np.mean(dsub, axis=1)
    overall = float(np.mean(row_means))
    return -0.5 * (dsub - row_means[None, :] - row_means[:, None] + overall)


def _correlations(benz: np.ndarray, active: np.ndarray, n_sit: int) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        compact = np.corrcoef(benz, rowvar=False)
    compact = np.atleast_2d(compact)
    result = np.full((n_sit, n_sit), np.nan, dtype=float)
    result[np.ix_(active, active)] = compact
    return result


def _princomp_cor(x: np.ndarray) -> dict:
    """Reproduce ``stats::princomp(cor=TRUE)`` and its ``cov.wt`` route."""
    matrix = np.asarray(x, dtype=float)
    n_observations = matrix.shape[0]
    weights = np.full(n_observations, 1.0 / n_observations)
    center = np.sum(weights[:, None] * matrix, axis=0)
    centered = matrix - center
    with np.errstate(divide="ignore", invalid="ignore"):
        weighted_centered = np.sqrt(weights)[:, None] * centered
        covariance = (weighted_centered.T @ weighted_centered) / (1.0 - np.sum(weights**2))
        covariance *= 1.0 - 1.0 / n_observations
        scale = np.sqrt(np.diag(covariance))
        correlation = covariance / np.outer(scale, scale)
    if not np.all(np.isfinite(correlation)):
        raise np.linalg.LinAlgError("Eigenvalues did not converge")
    eigenvalues, eigenvectors = linalg.eigh(
        correlation,
        lower=True,
        driver="evr",
        check_finite=False,
    )
    order = np.arange(eigenvalues.size - 1, -1, -1)
    eigenvalues = eigenvalues[order]
    loadings = eigenvectors[:, order]
    signs = np.where(loadings[0] < 0, -1.0, 1.0)
    loadings *= signs
    standardized = centered / scale
    scores = standardized @ loadings
    return {
        "sdev": np.sqrt(np.maximum(eigenvalues, 0.0)),
        "loadings": loadings,
        "scores": scores,
    }


def seqemlt(
    seqdata: STSList,
    a: float = 1.0,
    b: float = 1.0,
    weighted: bool = True,
) -> EMLTResult:
    """Compute TraMineRextras-compatible EMLT situation and case coordinates."""
    if not isinstance(seqdata, STSList):
        raise TypeError("seqdata must be an STSList state sequence object")
    states = list(sequence_alphabet(seqdata))
    state_codes = np.array([seqdata._code_for_state(state) for state in states], dtype=int)
    n_states = len(states)
    period = seqdata.seqlen
    data = np.asarray(seqdata.data, dtype=int)
    sit_time, sit_state = _situation_layout(n_states, period)
    situations = [f"{states[state]}{time + 1}" for time, state in zip(sit_time, sit_state)]
    frequencies = _situation_frequencies(data, sit_time, state_codes[sit_state])
    disjunctive = _disjunctive(data, sit_time, state_codes[sit_state])
    case_weights = None
    if weighted and seqdata.weights is not None:
        case_weights = np.asarray(seqdata.weights, dtype=float)
        if case_weights.shape != (seqdata.nseq,):
            raise ValueError("sequence weights must have one value per case")
    transitions = _transition_rates(disjunctive, case_weights)
    profiles = _situation_profiles(transitions, sit_time, float(a), float(b))
    active = np.flatnonzero(frequencies != 0)
    if active.size == 0:
        # No situation is observed (no sequence, or only missing cells): R's seqemlt
        # stops there (TX-REP-20261023-005).
        raise ValueError("subscript out of bounds: no state is observed in seqdata")
    distances = _profile_distance(profiles, active)
    benz = _benzecri(distances[np.ix_(active, active)])
    principal = _princomp_cor(benz)
    coordinates = disjunctive[:, active] @ principal["scores"]
    correlation = _correlations(benz, active, sit_time.size)
    return EMLTResult(
        states=states,
        period=period,
        situations=situations,
        sit_freq=frequencies,
        disjunctive=disjunctive,
        sit_transrate=transitions,
        sit_profil=profiles,
        distance_matrix=distances,
        benz_covariance=benz,
        sdev=principal["sdev"],
        loadings=principal["loadings"],
        coordinates=coordinates,
        correlation=correlation,
        a=a,
        b=b,
        weighted=bool(weighted),
    )
