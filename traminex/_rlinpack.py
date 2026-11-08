"""R's LINPACK QR rank rule (private).

R's ``qr()``, ``lm`` and ``glm`` decide a design's rank with ``dqrdc2`` (R's ``src/appl/dqrdc2.f``, Ross Ihaka's
modification of LINPACK's ``dqrdc``): Householder QR with limited column pivoting, in which a column whose remaining
norm falls below ``tol`` times its own original norm is moved to the end, so the rank is the number of columns left
in front. ``qr()`` uses ``tol = 1e-7``; ``glm.fit`` calls it with ``min(1e-7, epsilon / 1000) = 1e-11``
(TX-REP-20261019-002).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.linalg.blas import dnrm2


def dqrdc2_pivot(matrix: np.ndarray, tol: float) -> tuple[int, list[int]]:
    """R's rank and column order for ``matrix``: ``(qr(x, tol)$rank, qr(x, tol)$pivot - 1)``.

    A line-for-line port of dqrdc2.f's decomposition loop; only the rank and the pivot are returned. The pivot lists
    0-based column positions: the first ``rank`` are the columns R keeps, in their order, and the rest those it
    moved to the end.
    """
    x = np.array(matrix, dtype=float, copy=True)
    if x.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    n, p = x.shape
    jpvt = list(range(p))
    qraux = np.zeros(p)
    work1 = np.zeros(p)
    work2 = np.zeros(p)
    if n > 0:
        for j in range(p):
            qraux[j] = dnrm2(x[:, j])
            work1[j] = qraux[j]
            work2[j] = qraux[j] if qraux[j] != 0.0 else 1.0
    k = p + 1
    for l in range(1, min(n, p) + 1):
        c = l - 1
        # Cycle the columns from l to p left-to-right until one with a non-negligible norm is in place.
        while not (l >= k or qraux[c] >= work2[c] * tol):
            x[:, c:] = np.roll(x[:, c:], -1, axis=1)
            jpvt[c:] = jpvt[c + 1 :] + [jpvt[c]]
            for values in (qraux, work1, work2):
                values[c:] = np.roll(values[c:], -1)
            k -= 1
        if l == n:
            continue
        # Householder transformation for column l.
        nrmxl = dnrm2(x[c:, c])
        if nrmxl == 0.0:
            continue
        if x[c, c] != 0.0:
            nrmxl = math.copysign(nrmxl, x[c, c])
        x[c:, c] *= 1.0 / nrmxl
        x[c, c] = 1.0 + x[c, c]
        # Apply it to the remaining columns, downdating their norms.
        for j in range(c + 1, p):
            t = -np.dot(x[c:, c], x[c:, j]) / x[c, c]
            x[c:, j] += t * x[c:, c]
            if qraux[j] != 0.0:
                tt = max(1.0 - (abs(x[c, j]) / qraux[j]) ** 2, 0.0)
                if abs(tt) >= 1e-6:
                    qraux[j] = qraux[j] * math.sqrt(tt)
                else:
                    qraux[j] = dnrm2(x[c + 1 :, j])
                    work1[j] = qraux[j]
        qraux[c] = x[c, c]
        x[c, c] = -nrmxl
    return min(k - 1, n), jpvt
