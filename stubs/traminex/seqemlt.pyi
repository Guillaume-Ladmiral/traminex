from dataclasses import dataclass

import numpy as np

from .sequences import STSList as STSList

@dataclass
class EMLTResult:
    states: list[str]
    period: int
    situations: list[str]
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

def seqemlt(
    seqdata: STSList, a: float = 1.0, b: float = 1.0, weighted: bool = True
) -> EMLTResult: ...
