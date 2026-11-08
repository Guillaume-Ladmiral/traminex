"""Reference renderer for TraMineRextras ``plot.SAMM``."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ...samm import SAMMResult, _subset_sequences, _validate_samm
from ..legend import _draw_state_legend
from ..utils import _close_new_figures_on_error
from .dispatch import seqplot

_SUPPORTED_TYPES = {"i", "I", "d", "Ht", "dH", "f", "mt", "ms", "r", "rf"}


def _origin_groups(x: SAMMResult):
    _validate_samm(x)
    if "transition" not in x.frame.columns:
        raise ValueError("SAMMResult frame is missing the transition column")
    transition = x.frame["transition"]
    if transition.isna().any():
        raise ValueError("SAMMResult transition column must not contain missing values")
    selected_source_rows = np.flatnonzero(transition.to_numpy(dtype=bool))
    if selected_source_rows.size == 0:
        raise ValueError("SAMMResult has no transition rows to plot")

    selected_sequences = _subset_sequences(x.sequences, selected_source_rows)
    origin_series = x.frame.iloc[selected_source_rows][x.state_columns[0]]
    categorical = pd.Categorical(origin_series)
    codes = np.asarray(categorical.codes, dtype=int)
    used = set(int(code) for code in codes if code >= 0)
    level_codes = [code for code in range(len(categorical.categories)) if code in used]
    groups = []
    for code in level_codes:
        positions = np.flatnonzero(codes == code)
        origin = categorical.categories[code]
        groups.append(
            (
                f"Transition out of {origin}",
                positions,
                _subset_sequences(selected_sequences, positions),
            )
        )
    return selected_source_rows, groups


def _resolve_axes(*, count: int, two_panel: bool, supplied, figsize):
    expected = count * (2 if two_panel else 1)
    if supplied is None:
        if two_panel:
            fig, axes = plt.subplots(
                count,
                2,
                figsize=figsize,
                squeeze=False,
                constrained_layout=True,
            )
        else:
            fig, axes = plt.subplots(
                count,
                1,
                figsize=figsize,
                squeeze=False,
                constrained_layout=True,
            )
        return fig, np.asarray(axes, dtype=object).reshape(-1)

    axes = np.asarray(supplied, dtype=object).reshape(-1)
    if axes.size != expected:
        raise ValueError(f"supplied axes must contain {expected} axes")
    if any(not hasattr(axis, "figure") for axis in axes):
        raise TypeError("supplied axes must contain matplotlib Axes objects")
    fig = axes[0].figure
    if any(axis.figure is not fig for axis in axes):
        raise ValueError("all supplied axes must belong to the same figure")
    return fig, axes


@_close_new_figures_on_error
def plot_SAMM(x: SAMMResult, type: str = "d", **kwargs: object):
    """Plot transition windows, grouped and labelled by their origin state."""
    if not isinstance(x, SAMMResult):
        raise TypeError("x must be a SAMMResult")
    if type not in _SUPPORTED_TYPES:
        raise ValueError(f"type must be one of {sorted(_SUPPORTED_TYPES)}")

    selected_source_rows, groups = _origin_groups(x)
    options = dict(kwargs)
    supplied_axes = options.pop("axes", None)
    supplied_ax = options.pop("ax", None)
    if supplied_axes is not None and supplied_ax is not None:
        raise ValueError("supply axes or ax, not both")
    if supplied_ax is not None:
        supplied_axes = [supplied_ax]
    figsize = options.pop("figsize", None)
    # R's plot.SAMM is seqplot(..., group = origin), which draws one state legend
    # for the figure, none for Ht (seqsamm.R:56-60, seqplot.R:210).
    with_legend = options.pop("with_legend", "auto")
    if with_legend not in (True, False, "auto"):
        raise ValueError("with_legend must be True, False or 'auto'")
    which_plot = options.get("which_plot", "both") if type == "rf" else None
    two_panel = type == "rf" and which_plot == "both"
    if figsize is None:
        figsize = (
            (10, max(3.0 * len(groups), 3.0))
            if two_panel
            else (
                8,
                max(3.0 * len(groups), 3.0),
            )
        )
    fig, axes = _resolve_axes(
        count=len(groups),
        two_panel=two_panel,
        supplied=supplied_axes,
        figsize=figsize,
    )

    distance = None
    if type in {"r", "rf"}:
        if "diss" not in options:
            raise ValueError(f"diss for transition rows is required when type={type!r}")
        distance = np.asarray(options.pop("diss"), dtype=float)
        expected = selected_source_rows.size
        if distance.shape != (expected, expected) or not np.all(np.isfinite(distance)):
            raise ValueError(
                "diss for transition rows must be a finite square matrix with "
                "one row per selected transition"
            )

    axis_offset = 0
    for label, positions, sequences in groups:
        group_options = dict(options)
        if distance is not None:
            group_options["diss"] = distance[np.ix_(positions, positions)]
        if type == "r" and sequences.nseq < 2:
            raise ValueError("type='r' requires at least two transition rows per origin group")
        if type == "rf" and group_options.get("k") is None and sequences.nseq < 10:
            raise ValueError("type='rf' requires an explicit k for origin groups smaller than 10")

        if two_panel:
            target_axes = (axes[axis_offset], axes[axis_offset + 1])
            axis_offset += 2
            _, rendered_axes = seqplot(
                sequences,
                type=type,
                ax=target_axes,
                with_legend=False,
                **group_options,
            )
            for axis in rendered_axes:
                axis.set_title(f"{label}: {axis.get_title()}")
        else:
            target_axis = axes[axis_offset]
            axis_offset += 1
            _, rendered_axis = seqplot(
                sequences,
                type=type,
                ax=target_axis,
                with_legend=False,
                **group_options,
            )
            rendered_axis.set_title(label)

    if with_legend and type != "Ht":
        _draw_state_legend(
            x.sequences,
            fig if supplied_axes is None else axes[0],
            options.get("cpal", options.get("cmap")),
        )
    returned_axes = axes.reshape(len(groups), 2) if two_panel else axes
    return fig, returned_axes
