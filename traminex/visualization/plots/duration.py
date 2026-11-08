from __future__ import annotations

import inspect
from typing import Any, cast

import matplotlib.pyplot as plt
import numpy as np

from ...sequences import _real_state_codes, _require_stslist
from ..utils import _close_new_figures_on_error, _import_seaborn


def _boxplot(ax, values, *, tick_labels=None) -> None:
    """Draw a vertical boxplot without deprecated Matplotlib arguments."""
    parameters = inspect.signature(ax.boxplot).parameters
    kwargs: dict[str, Any] = {
        "patch_artist": True,
        ("orientation" if "orientation" in parameters else "vert"): (
            "vertical" if "orientation" in parameters else True
        ),
    }
    if tick_labels is not None:
        label_parameter = "tick_labels" if "tick_labels" in parameters else "labels"
        kwargs[label_parameter] = tick_labels
    cast(Any, ax.boxplot)(values, **kwargs)


@_close_new_figures_on_error
@_require_stslist
def duration_boxplot(
    x,
    labels: np.ndarray | None = None,
    *,
    by_state: bool = False,
    seaborn: bool = True,
    figsize=(8, 4),
):
    """Box/violin plots of durations per cluster or overall.

    ``labels`` must be one-dimensional with one value per sequence. Label values
    are preserved as supplied, including strings and NumPy scalar values.
    With ``by_state=True`` each state's distribution of per-sequence time
    spent in it (``seqistatd``) is drawn, split by cluster when ``labels`` is
    given: seaborn violins, or one matplotlib box per state (and cluster)
    when ``seaborn=False``. The missing and void markers are not states and
    are left out.
    """
    import pandas as pd

    from ...stats import seqistatd

    label_values = None
    if labels is not None:
        label_values = np.asarray(labels)
        if label_values.ndim != 1 or label_values.shape[0] != x.nseq:
            raise ValueError("labels must be one-dimensional with exactly one label per sequence")

    D = seqistatd(x, include_missing=False)
    # The marker columns are zero here and are not substantive states.
    states = [(k, x.alphabet[k]) for k in _real_state_codes(x)]
    df_rows = []
    for i in range(x.nseq):
        if by_state:
            for k, s in states:
                df_rows.append(
                    {
                        "seq": i,
                        "state": s,
                        "duration": int(D[i, k]),
                        "label": None if label_values is None else label_values[i],
                    }
                )
        else:
            df_rows.append(
                {
                    "seq": i,
                    "duration": int(D[i].sum()),
                    "label": None if label_values is None else label_values[i],
                }
            )
    df = pd.DataFrame(df_rows)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    if seaborn and by_state:
        sns = _import_seaborn("duration_boxplot")

        # seaborn 0.13 renamed violinplot's ``scale`` to ``density_norm`` and
        # deprecated ``scale``. seaborn 0.12, the declared floor, has no
        # ``density_norm`` and silently drops it through ``**kwargs``, leaving
        # the violins area-scaled (audit S7-11).
        width_norm = (
            {"density_norm": "width"}
            if "density_norm" in inspect.signature(sns.violinplot).parameters
            else {"scale": "width"}
        )
        if labels is not None:
            sns.violinplot(
                x="duration",
                y="state",
                hue="label",
                data=df,
                ax=ax,
                **width_norm,
            )
        else:
            sns.violinplot(x="duration", y="state", data=df, ax=ax, **width_norm)
        ax.set_title(
            "Durations"
            + (" by state" if by_state else "")
            + (" per cluster" if labels is not None else "")
        )
        return fig, ax
    # Matplotlib boxplots avoid seaborn's use of arguments deprecated by
    # Matplotlib 3.11 while remaining compatible with our Matplotlib 3.6 floor.
    if by_state:
        # One box per state -- per state and cluster with labels -- the groups
        # the seaborn path draws (2026-09-22 audit S6-04).
        cluster_labels = [None] if labels is None else sorted(df["label"].unique())
        groups = []
        tick_labels = []
        for _, state in states:
            for cluster in cluster_labels:
                selected = df["state"] == state
                if cluster is not None:
                    selected &= df["label"] == cluster
                groups.append(df[selected]["duration"].values)
                tick_labels.append(state if cluster is None else f"{state}\n{cluster}")
        _boxplot(ax, groups, tick_labels=tick_labels)
    elif labels is not None:
        cluster_labels = sorted(df["label"].unique())
        groups = [df[df["label"] == c]["duration"].values for c in cluster_labels]
        _boxplot(ax, groups, tick_labels=cluster_labels)
    else:
        _boxplot(ax, df["duration"].values)
    ax.set_title(
        "Durations"
        + (" by state" if by_state else "")
        + (" per cluster" if labels is not None else "")
    )
    return fig, ax
