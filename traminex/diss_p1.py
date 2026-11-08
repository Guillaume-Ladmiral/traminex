"""Priority 1 dissimilarity functions: representative set selection."""

from __future__ import annotations

import math
from typing import Literal, Optional

import numpy as np
import pandas as pd

from ._deprecated import MISSING, resolve_deprecated_alias


def _exact_near_ties(score: np.ndarray, matrix: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Re-sum with correct rounding the scores that summation error could have split.

    ``score[i] = sum_j matrix[i, j] * weights[j]``. Floating-point summation can make two sums
    that are equal in exact arithmetic differ, and which one is larger then depends on the
    order of the additions (the BLAS build): on mvad R's OpenBLAS and NumPy ordered two such
    candidates differently (TX-AUD-20261019-017). A sum of ``n`` terms is within
    ``n * eps * sum|terms|`` of its exact value, so every score that close to another one is
    summed again with ``math.fsum``; an exact tie is then a tie, and the stable order keeps the
    lower row first, as R's ``order()`` does on a tie.
    """
    count = score.size
    if count < 2:
        return score
    magnitude = np.abs(matrix) @ np.abs(weights)
    error = matrix.shape[1] * np.finfo(float).eps * magnitude
    order = np.argsort(score, kind="stable")
    gaps = np.diff(score[order])
    close = gaps <= error[order][:-1] + error[order][1:]
    if not close.any():
        return score
    rows = set(order[:-1][close].tolist()) | set(order[1:][close].tolist())
    refined = score.copy()
    for row in rows:
        refined[row] = math.fsum((matrix[row] * weights).tolist())
    return refined


def _positive_weight_total(weights: np.ndarray) -> float:
    """The total of ``weights``, refused unless positive (2026-10-25 study, finding F-05).

    R stops too, on the NaN of a division by this total ("missing value where
    TRUE/FALSE needed", dissrep.R:104-105); this names the cause instead.
    """
    total = float(weights.sum())
    if not total > 0:
        raise ValueError(f"the weights must have a positive total, got {total!r}")
    return total


def dissrep(
    diss: np.ndarray,
    *,
    criterion: Literal["density", "freq", "dist", "random"] = "density",
    score: Optional[np.ndarray] = None,
    decreasing: bool = True,
    coverage: float = 0.25,
    nrep: Optional[int] = None,
    pradius: float = 0.10,
    dmax: Optional[float] = None,
    weights: Optional[np.ndarray] = None,
    random_state: int | np.random.Generator | None = None,
    trep=MISSING,
    tsim=MISSING,
) -> dict:
    """Extract a set of representative sequences from a dissimilarity matrix.

    Selects representative sequences based on various criteria (density, frequency,
    distance) with optional coverage or count constraints. Returns representatives
    and quality statistics.

    Parameters
    ----------
    diss : np.ndarray
        Dissimilarity matrix (square, shape (n, n)) or distance object
    criterion : {"density", "freq", "dist", "random"}, default "density"
        Selection criterion:
        - "density": neighborhood density (most neighbors within pradius)
        - "freq": frequency (number of exact duplicates)
        - "dist": minimum sum of distances
        - "random": random selection
    score : np.ndarray, optional
        Custom score array of length n for ranking candidates
    decreasing : bool, default True
        Sort order for scores (True for high-to-low, False for low-to-high)
    coverage : float, default 0.25
        Target coverage as proportion of total weight (used if nrep is None)
    nrep : int, optional
        Fixed number of representatives to select (overrides coverage)
    pradius : float, default 0.10
        Neighborhood radius as proportion of dmax (0 to 1)
    dmax : float, optional
        Maximum theoretical distance (default: max of diss)
    weights : np.ndarray, optional
        Sequence weights (default: uniform weights)
    random_state : int or numpy.random.Generator, optional
        Local random-number source for ``criterion="random"``. An integer gives
        reproducible independent calls and a generator is consumed in place.
        When omitted, the legacy NumPy global random stream is used.

    Returns
    -------
    dict
        Dictionary with keys:
        - representatives: array of selected representative indices (0-based)
        - statistics: DataFrame with quality measures per representative
          (na, na_pct, nb, nb_pct, SD, MD, DC, V, Q)
        - quality: overall quality measure (at most 1); NaN, as in R, when
          the total distance to the centre is 0 (e.g. all dissimilarities
          zero), which also makes the ``Q`` statistics NaN
        - distances: array of distances from each sequence to nearest representative
        - rep_group: array indicating which representative each sequence is assigned to
        - metadata: dict with criterion, dmax, pradius, weights_sum

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqdist, dissrep
    >>> data = np.array([
    ...     ['A', 'A', 'B'],
    ...     ['A', 'B', 'B'],
    ...     ['A', 'A', 'B'],
    ...     ['B', 'B', 'C'],
    ... ])
    >>> sts = seqdef(data)
    >>> D = seqdist(sts, method='OM', subm=2.0, indel=1.0)
    >>> result = dissrep(D, criterion='density', coverage=0.5)
    >>> result['representatives']
    array([0])
    """
    coverage = resolve_deprecated_alias(
        coverage,
        trep,
        canonical_name="coverage",
        alias_name="trep",
        function_name="dissrep",
        default=0.25,
    )
    pradius = resolve_deprecated_alias(
        pradius,
        tsim,
        canonical_name="pradius",
        alias_name="tsim",
        function_name="dissrep",
        default=0.10,
    )

    # Convert dist object to matrix if necessary
    if hasattr(diss, "shape") and len(diss.shape) == 1:
        # Condensed distance matrix - convert to square
        from scipy.spatial.distance import squareform

        diss = squareform(diss)
    elif not isinstance(diss, np.ndarray):
        diss = np.asarray(diss)

    if diss.ndim != 2 or diss.shape[0] != diss.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")

    nbobj = diss.shape[0]
    # R: diss0 <- all(diss == 0) (dissrep.R:16); see the nb rule below.
    diss0 = bool(np.all(diss == 0))

    # Initialize weights
    if weights is None:
        weights = np.ones(nbobj, dtype=float)
    else:
        weights = np.asarray(weights, dtype=float)
        if weights.shape[0] != nbobj:
            raise ValueError("weights must have length equal to number of objects")

    weights_sum = _positive_weight_total(weights)

    # Determine max distance
    if dmax is None:
        dmax = float(diss.max())

    # Validate and convert pradius to absolute value
    if pradius < 0 or pradius > 1:
        raise ValueError("pradius must be between 0 and 1")
    pradius_abs = dmax * pradius

    # Compute scores based on criterion
    if criterion == "density":
        # Neighborhood density: count neighbors within pradius
        neighbours = diss < pradius_abs
        score = neighbours @ weights
        decreasing = True

    elif criterion == "freq":
        # Frequency: count exact duplicates
        neighbours = diss == 0
        score = neighbours @ weights
        decreasing = True

    elif criterion == "dist":
        # Sum of distances
        score = diss @ weights
        decreasing = False

    elif criterion == "random":
        # Random order
        if random_state is None:
            score = np.random.rand(nbobj)
        else:
            rng = (
                random_state
                if isinstance(random_state, np.random.Generator)
                else np.random.default_rng(random_state)
            )
            score = np.asarray(rng.random(nbobj), dtype=float)
        decreasing = False

    elif score is not None:
        # Use provided score
        pass
    else:
        raise ValueError(f"Unknown criterion: {criterion}")

    if score is None:
        raise ValueError("score is None after criterion selection")

    if len(score) != nbobj:
        raise ValueError(f"score must have length {nbobj}, got {len(score)}")

    # Sort candidates by score, matching R's order(score, decreasing=) which
    # breaks ties by ascending original index (a stable sort), not numpy's
    # default unstable quicksort.
    score = np.asarray(score, dtype=float)
    if criterion in ("density", "freq"):
        score = _exact_near_ties(score, neighbours.astype(float), weights)
    elif criterion == "dist":
        score = _exact_near_ties(score, diss, weights)
    if decreasing:
        score_sort = np.argsort(-score, kind="stable")
    else:
        score_sort = np.argsort(score, kind="stable")

    # Reorder distance matrix
    rep_dist = diss[np.ix_(score_sort, score_sort)]

    # Select representatives
    idx = 0
    idxrep: list[int] = []

    if nrep is None:
        # Coverage-based selection
        pctrep = 0.0

        while pctrep < coverage and idx < nbobj:
            if idx == 0 or np.all(rep_dist[idx, idxrep] > pradius_abs):
                idxrep.append(idx)
                # Compute coverage
                tempm = rep_dist[:, idxrep]
                nbnear = ((tempm < pradius_abs).any(axis=1) * weights[score_sort]).sum()
                pctrep = nbnear / weights_sum
            idx += 1

        nbkeep = len(idxrep)

    else:
        # Fixed number of representatives
        repcount = 0

        while repcount < nrep and idx < nbobj:
            if idx == 0 or np.all(rep_dist[idx, idxrep] > pradius_abs):
                idxrep.append(idx)
                repcount += 1
            idx += 1

        nbkeep = len(idxrep)
        # Compute coverage for reporting
        if nbkeep > 0:
            tempm = rep_dist[:, idxrep]
            nbnear = ((tempm < pradius_abs).any(axis=1) * weights[score_sort]).sum()
            pctrep = nbnear / weights_sum
        else:
            pctrep = 0.0

    if nbkeep == 0:
        raise ValueError("No representatives selected")

    # Get distances to representatives (original ordering)
    dist_repseq = diss[:, score_sort[idxrep]]

    # Distance of each object to the global (weighted) centre, exactly as R
    # dissrep: dc.tot <- disscenter(diss, weights=weights).
    from .diss import disscenter

    dc_tot = np.asarray(disscenter(diss, weights=weights), dtype=float)

    if nbkeep == 1:
        # Single representative: every object is assigned to it.
        na_values = np.array([float(weights_sum)])
        sd_values = np.array([float((dist_repseq[:, 0] * weights).sum())])
        md_values = sd_values / weights_sum
        nb_values = np.array([float(((dist_repseq[:, 0] < pradius_abs) * weights).sum())])
        dc_values = np.array([float((dc_tot * weights).sum())])
        variance_values = dc_values / weights_sum
        minidx = np.zeros(nbobj, dtype=int)

    else:
        # Multiple representatives
        # Assign each sequence to nearest representative
        minidx = np.argmin(dist_repseq, axis=1)

        # Create masked distances (only nearest representative)
        dist_masked = dist_repseq.copy()
        for i in range(nbobj):
            mask = np.ones(nbkeep, dtype=bool)
            mask[minidx[i]] = False
            dist_masked[i, mask] = np.nan

        # Per-representative statistics
        na_values = np.zeros(nbkeep, dtype=float)
        sd_values = np.zeros(nbkeep, dtype=float)
        md_values = np.zeros(nbkeep, dtype=float)
        nb_values = np.zeros(nbkeep, dtype=float)
        dc_values = np.zeros(nbkeep, dtype=float)
        variance_values = np.zeros(nbkeep, dtype=float)

        for i in range(nbkeep):
            sel = minidx == i
            sel_idx = np.flatnonzero(sel)
            # R: na <- colSums((!is.na(dist.repseq))*weights) -> weighted count
            # of objects assigned (nearest) to representative i.
            na_values[i] = float(weights[sel].sum())
            sd_values[i] = float(np.nansum(dist_masked[:, i] * weights))
            # R: MD <- SD/na (dissrep.R:153), a plain division: a representative
            # whose assigned cases all weigh 0 has MD = 0/0 = NaN, not 0.
            with np.errstate(divide="ignore", invalid="ignore"):
                md_values[i] = np.float64(sd_values[i]) / np.float64(na_values[i])
            nb_values[i] = float(((dist_masked[:, i] < pradius_abs) * weights).sum())
            dc_values[i] = float((dc_tot[sel] * weights[sel]).sum())
            # R: V[i] <- mean(disscenter(diss[sel, sel], weights=weights[sel]))
            sub = diss[np.ix_(sel_idx, sel_idx)]
            variance_values[i] = float(np.mean(disscenter(sub, weights=weights[sel])))

    if diss0:
        # R: if (diss0) nb <- sum(weights, na.rm=TRUE) (dissrep.R:178): when all
        # objects are identical every one counts as similar to the
        # representative, although no distance is below a zero radius.
        nb_values = np.full(nbkeep, float(np.nansum(weights)))

    # Overall statistics
    na_total = float(na_values.sum())
    nb_total = float(nb_values.sum())
    sd_total = float(sd_values.sum())
    md_total = sd_total / weights_sum
    dc_total = float((dc_tot * weights).sum())
    variance_total = dc_total / weights_sum

    # Append totals
    na_values = np.append(na_values, na_total)
    nb_values = np.append(nb_values, nb_total)
    sd_values = np.append(sd_values, sd_total)
    md_values = np.append(md_values, md_total)
    dc_values = np.append(dc_values, dc_total)
    variance_values = np.append(variance_values, variance_total)

    # Compute percentages
    pcta = (na_values / weights_sum) * 100.0
    pctb = (nb_values / weights_sum) * 100.0

    # Quality measure Q = (DC - SD) / DC * 100 (dissrep.R:200). R divides
    # plainly, so a zero DC (e.g. all dissimilarities zero: nothing to explain)
    # gives NaN or -Inf, not 0; numpy's warning for that case is silenced.
    with np.errstate(divide="ignore", invalid="ignore"):
        quality_values = (dc_values - sd_values) / dc_values * 100.0

    # Build statistics DataFrame
    stats = pd.DataFrame(
        {
            "na": na_values,
            "na(%)": pcta,
            "nb": nb_values,
            "nb(%)": pctb,
            "SD": sd_values,
            "MD": md_values,
            "DC": dc_values,
            "V": variance_values,
            "Q": quality_values,
        },
        index=[f"r{i+1}" for i in range(nbkeep)] + ["Total"],
    )

    # Overall quality (dissrep.R:188), plain division as in R.
    with np.errstate(divide="ignore", invalid="ignore"):
        quality = float(np.float64(dc_total - sd_total) / np.float64(dc_total))

    # Get representative indices in original ordering
    representatives = score_sort[idxrep]

    # Get distances to nearest representative for each sequence
    dist_to_rep = np.min(dist_repseq, axis=1)

    # Build result dictionary
    result = {
        "representatives": representatives,
        "scores": score.copy(),
        "statistics": stats,
        "quality": quality,
        "distances": dist_to_rep,
        "rep_group": minidx,
        "metadata": {
            "criterion": criterion,
            "dmax": dmax,
            "pradius": pradius_abs,
            "weights_sum": weights_sum,
            "coverage": pctrep,
            "nrep": nbkeep,
        },
    }

    return result
