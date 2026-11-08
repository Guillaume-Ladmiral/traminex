from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ...sequences import STSList, _alphabet_codes, _remap_issues, _require_stslist, alphabet
from ...stats import seqstatd
from ..ticks import _apply_time_axis_on_cells, apply_time_axis
from ..utils import _close_new_figures_on_error, _import_seaborn, determine_layout


def _chronogram_rows(x: STSList, include_missing: bool) -> tuple[list[str], np.ndarray]:
    """Row names and frequencies (states x time) of a chronogram: ``seqstatd``'s states.

    As R's ``seqstatd`` (seqstatd.R:3, 12): the real states, plus the missing
    state when ``include_missing`` (a zero row when the object has none); the
    void marker never (TX-AUD-20261019-052).
    """
    freq = np.asarray(seqstatd(x, include_missing=include_missing)["Frequencies"], dtype=float)
    rows = [
        freq[:, code] if code is not None else np.zeros(freq.shape[0])
        for code in _alphabet_codes(x, with_missing=include_missing)
    ]
    return alphabet(x, with_missing=include_missing), np.vstack(rows)


def _cluster_labels(x: STSList, labels) -> np.ndarray:
    """``labels`` as an array of one cluster label per sequence (TX-AUD-20261011-029)."""
    labels = np.asarray(labels)
    if labels.shape[:1] != (x.nseq,):
        raise ValueError("labels length must equal number of sequences")
    return labels


@_close_new_figures_on_error
@_require_stslist
def chronogram(
    x: STSList, *, include_missing: bool = False, figsize=(8, 3), seaborn: bool = False, ax=None
):
    """Chronogram: heatmap of state frequencies over time (states × time).

    The rows are ``seqstatd``'s states: the real states, plus the missing state
    when ``include_missing``.
    """
    x = x._canonicalized()
    names, freq = _chronogram_rows(x, include_missing)  # K × L
    K, L = freq.shape
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    if seaborn:
        sns = _import_seaborn("chronogram")
        sns.heatmap(freq, cmap="viridis", cbar=True, ax=ax)
        ax.set_ylabel("State")
        # seaborn draws time position j as the cell [j, j + 1]: centre the
        # labels on the cells (2026-09-22 audit S6-21).
        _apply_time_axis_on_cells(ax, x, L)
        ax.set_xlabel("Time")
        ax.set_yticks(np.arange(K) + 0.5)
        ax.set_yticklabels(names, rotation=0)
        ax.set_title("Chronogram (state frequency heatmap)")
        return fig, ax
    im = ax.imshow(freq, aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_yticks(np.arange(K))
    ax.set_yticklabels(names)
    apply_time_axis(ax, x, L)
    ax.set_xlabel("Time")
    ax.set_ylabel("State")
    ax.set_title("Chronogram (state frequency heatmap)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return fig, ax


@_close_new_figures_on_error
@_require_stslist
def cluster_chronograms(
    x: STSList,
    labels: np.ndarray,
    *,
    include_missing: bool = False,
    cols: int = 3,
    figsize=(8, 3),
    shared_scale: bool = True,
    layout: str = "column",
    nrows: int | None = None,
    ncols: int | None = None,
):
    """Grid of chronograms, one per cluster label.

    ``labels`` holds one label per sequence (else ValueError). If
    shared_scale=True, all subplots share a common color scale.
    """
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    labels = _cluster_labels(x, labels)
    uniq = np.unique(labels)
    nrows_eff, ncols_eff = determine_layout(
        len(uniq), layout=layout, nrows=nrows, ncols=ncols, default_cols=cols
    )
    if nrows_eff == 0 or ncols_eff == 0:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.axis("off")
        return fig, np.array([[ax]])
    fig, axes = plt.subplots(
        nrows_eff, ncols_eff, figsize=(ncols_eff * figsize[0] / 2.0, nrows_eff * figsize[1])
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([[axes]])
    axes = axes.reshape(nrows_eff, ncols_eff)
    # Precompute frequencies per cluster if shared scale
    freq_list: list[tuple[list[str], np.ndarray] | None] = []
    vmax = 0.0
    for lab in uniq:
        rows = np.where(labels == lab)[0]
        if rows.size == 0:
            freq_list.append(None)
            continue
        Xc = STSList(
            data=x.data[rows],
            alphabet=list(x.alphabet),
            labels=None if x.labels is None else list(x.labels),
            cpal=None if x.cpal is None else list(x.cpal),
            xtlab=None if x.xtlab is None else list(x.xtlab),
            xtstep=x.xtstep,
            tick_last=x.tick_last,
            nr=x.nr,
            void=x.void,
            nr_code=x.nr_code,
            void_code=x.void_code,
            weights=None if x.weights is None else x.weights[rows],
            start=x.start,
            ids=None if x.ids is None else [x.ids[int(i)] for i in rows],
            missing_color=x.missing_color,
            issues=_remap_issues(x.issues, rows),
        )
        names, freq = _chronogram_rows(Xc, include_missing)
        freq_list.append((names, freq))
        vmax = max(vmax, float(freq.max()))
    vmin = 0.0
    last_im = None
    for idx, lab in enumerate(uniq):
        ri, ci = divmod(idx, ncols_eff)
        ax = axes[ri, ci]
        item = freq_list[idx]
        if item is None:
            ax.axis("off")
            continue
        names, cluster_freq = item
        im = ax.imshow(
            cluster_freq,
            aspect="auto",
            interpolation="nearest",
            cmap="viridis",
            vmin=vmin if shared_scale else None,
            vmax=vmax if shared_scale else None,
        )
        last_im = im
        ax.set_yticks(np.arange(len(names)))
        ax.set_yticklabels(names)
        ax.set_title(f"Cluster {lab}")
        if ci == 0:
            ax.set_ylabel("State")
        if ri == nrows_eff - 1:
            apply_time_axis(ax, x, cluster_freq.shape[1])
            ax.set_xlabel("Time")
    # turn off any unused axes
    for j in range(len(uniq), nrows_eff * ncols_eff):
        ri, ci = divmod(j, ncols_eff)
        axes[ri, ci].axis("off")
    if shared_scale and last_im is not None:
        fig.colorbar(last_im, ax=axes.ravel().tolist(), fraction=0.046, pad=0.04)
    return fig, axes


@_close_new_figures_on_error
@_require_stslist
def cluster_chronograms_with_representatives(
    x: STSList,
    labels: np.ndarray,
    D: np.ndarray,
    *,
    include_missing: bool = False,
    cols: int = 3,
    figsize=(8, 3),
    shared_scale: bool = True,
    marker_size: int = 3,
    layout: str = "column",
    nrows: int | None = None,
    ncols: int | None = None,
):
    """Chronograms per cluster with overlaid representative (medoid) sequence markers.

    ``labels`` holds one label per sequence (else ValueError).
    """
    from ...cluster import medoid_index

    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    labels = _cluster_labels(x, labels)
    uniq = np.unique(labels)
    nrows_eff, ncols_eff = determine_layout(
        len(uniq), layout=layout, nrows=nrows, ncols=ncols, default_cols=cols
    )
    if nrows_eff == 0 or ncols_eff == 0:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.axis("off")
        return fig, np.array([[ax]])
    fig, axes = plt.subplots(
        nrows_eff, ncols_eff, figsize=(ncols_eff * figsize[0] / 2.0, nrows_eff * figsize[1])
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([[axes]])
    axes = axes.reshape(nrows_eff, ncols_eff)
    # Precompute frequencies and medoids
    freq_list: list[tuple[list[str], np.ndarray] | None] = []
    rep_idx = {}
    vmax = 0.0
    for lab in uniq:
        rows = np.where(labels == lab)[0]
        if rows.size == 0:
            freq_list.append(None)
            continue
        Xc = STSList(
            data=x.data[rows],
            alphabet=list(x.alphabet),
            labels=None if x.labels is None else list(x.labels),
            cpal=None if x.cpal is None else list(x.cpal),
            xtlab=None if x.xtlab is None else list(x.xtlab),
            xtstep=x.xtstep,
            tick_last=x.tick_last,
            nr=x.nr,
            void=x.void,
            nr_code=x.nr_code,
            void_code=x.void_code,
            weights=None if x.weights is None else x.weights[rows],
            start=x.start,
            ids=None if x.ids is None else [x.ids[int(i)] for i in rows],
            missing_color=x.missing_color,
            issues=_remap_issues(x.issues, rows),
        )
        names, freq = _chronogram_rows(Xc, include_missing)
        freq_list.append((names, freq))
        vmax = max(vmax, float(freq.max()))
        sub = D[np.ix_(rows, rows)]
        mlocal = medoid_index(sub)
        rep_idx[lab] = int(rows[mlocal])
    vmin = 0.0
    last_im = None
    for idx, lab in enumerate(uniq):
        ri, ci = divmod(idx, ncols_eff)
        ax = axes[ri, ci]
        item = freq_list[idx]
        if item is None:
            ax.axis("off")
            continue
        names, freq = item
        im = ax.imshow(
            freq,
            aspect="auto",
            interpolation="nearest",
            cmap="viridis",
            vmin=vmin if shared_scale else None,
            vmax=vmax if shared_scale else None,
        )
        last_im = im
        # Overlay representative markers
        rep = rep_idx.get(lab)
        if rep is not None:
            rep_row = x.data[rep]
            # plot points at (t, row of the state); a state without a row
            # (the void marker, or the missing one when excluded) is not marked
            row_of = {
                code: row
                for row, code in enumerate(_alphabet_codes(x, with_missing=include_missing))
                if code is not None
            }
            t = np.arange(x.seqlen)
            y = np.array([row_of.get(int(code), np.nan) for code in rep_row], dtype=float)
            ax.plot(
                t,
                y,
                marker="o",
                linestyle="None",
                markerfacecolor="white",
                markeredgecolor="black",
                markersize=marker_size,
            )
        ax.set_yticks(np.arange(len(names)))
        ax.set_yticklabels(names)
        ax.set_title(f"Cluster {lab}")
        if ci == 0:
            ax.set_ylabel("State")
        if ri == nrows_eff - 1:
            apply_time_axis(ax, x, freq.shape[1])
            ax.set_xlabel("Time")
    # turn off any unused axes
    for j in range(len(uniq), nrows_eff * ncols_eff):
        ri, ci = divmod(j, ncols_eff)
        axes[ri, ci].axis("off")
    if shared_scale and last_im is not None:
        fig.colorbar(last_im, ax=axes.ravel().tolist(), fraction=0.046, pad=0.04)
    return fig, axes
