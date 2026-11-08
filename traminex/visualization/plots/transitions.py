from __future__ import annotations

from typing import cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ...sequences import STSList, _require_stslist
from ...stats_utils import seqtrate
from ..utils import _close_new_figures_on_error, _import_seaborn


@_close_new_figures_on_error
@_require_stslist
def transition_heatmap(
    x: STSList,
    *,
    normalize: bool = True,
    weighted: bool = True,
    figsize=(5, 4),
    seaborn: bool = False,
    ax=None,
):
    """Heatmap of the transition rates of :func:`seqtrate` (counts with ``normalize=False``).

    The cells are ``seqtrate(x, weighted=weighted)``, or its ``count=True``
    counts, over the states seqtrate reports (TX-AUD-20261011-028; the heatmap
    counted every code, the missing and void markers included, always weighted).
    """
    x = x._canonicalized()
    table = cast(pd.DataFrame, seqtrate(x, weighted=weighted, count=not normalize))
    M = table.to_numpy(dtype=float)
    labels = [str(state) for state in table.index]
    columns = [str(state) for state in table.columns]
    K = M.shape[0]
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    if seaborn:
        sns = _import_seaborn("transition_heatmap")
        sns.heatmap(M, cmap="magma", cbar=True, ax=ax, xticklabels=columns, yticklabels=labels)
        ax.set_xlabel("To state")
        ax.set_ylabel("From state")
        ax.set_title("Transition heatmap" + (" (row-normalized)" if normalize else ""))
        return fig, ax
    im = ax.imshow(M, interpolation="nearest", cmap="magma")
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(K))
    ax.set_yticklabels(labels)
    ax.set_xlabel("To state")
    ax.set_ylabel("From state")
    ax.set_title("Transition heatmap" + (" (row-normalized)" if normalize else ""))
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return fig, ax
