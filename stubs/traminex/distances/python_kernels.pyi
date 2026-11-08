import numpy as np

class CostDiagnosticMessage(UserWarning): ...

def ham_distance(a: np.ndarray, b: np.ndarray) -> float: ...
def om_distance(
    a: np.ndarray,
    b: np.ndarray,
    indel: float = 1.0,
    subm: float | np.ndarray = 2.0,
    *,
    indel_by_pos: np.ndarray | None = None,
    subm_by_pos: np.ndarray | None = None,
    indel_by_state: np.ndarray | None = None,
) -> float: ...
def lcs_distance(a: np.ndarray, b: np.ndarray) -> float: ...
def lcp_distance(a: np.ndarray, b: np.ndarray) -> float: ...
def rlcp_distance(a: np.ndarray, b: np.ndarray) -> float: ...
def dur_l1_distance(a: np.ndarray, b: np.ndarray, K: int) -> float: ...
def dur_cosine_distance(a: np.ndarray, b: np.ndarray, K: int) -> float: ...
def dur_chi2_distance(a: np.ndarray, b: np.ndarray, K: int, eps: float = 1e-12) -> float: ...
def jsd_transition_distance(a: np.ndarray, b: np.ndarray, K: int, eps: float = 0.0) -> float: ...
def dhd_distance(
    a: np.ndarray, b: np.ndarray, freq_by_pos: np.ndarray, *, variant: str = "mean"
) -> float: ...
def pairwise(seqmat: np.ndarray, method: str = "OM", **kwargs) -> np.ndarray: ...
def nmsmst_counts(
    a: np.ndarray,
    b: np.ndarray,
    dur_a: np.ndarray,
    dur_b: np.ndarray,
    ksize: int,
) -> np.ndarray: ...
def omslen_distance(
    a: np.ndarray,
    b: np.ndarray,
    *,
    indel: float | np.ndarray,
    subm: float | np.ndarray,
    link: str,
    dur_a: np.ndarray,
    dur_b: np.ndarray,
    h: float = 0.5,
) -> float: ...
def omspell_distance(
    a: np.ndarray,
    b: np.ndarray,
    *,
    indel: float | np.ndarray,
    subm: float | np.ndarray,
    timecost: float,
    tpow: float = 1.0,
) -> float: ...
def subseq_counts(
    a: np.ndarray, b: np.ndarray, ksize: int, *, prox: np.ndarray | None = None
) -> np.ndarray: ...
def svrspell_counts(
    a: np.ndarray,
    b: np.ndarray,
    dur_a: np.ndarray,
    dur_b: np.ndarray,
    ksize: int,
    softmatch: np.ndarray,
) -> np.ndarray: ...
def twed_distance(
    a: np.ndarray,
    b: np.ndarray,
    *,
    nu: float,
    lam: float,
    subm: float | np.ndarray = 2.0,
    indel: float | None = None,
) -> float: ...
