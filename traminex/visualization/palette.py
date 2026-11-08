from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, ListedColormap, Normalize, is_color_like, to_rgba


def _state_code_norm(x) -> Normalize:
    """Fix an image's colour scale to the whole codebook of ``x``.

    ``imshow`` otherwise stretches the colormap over the codes that happen to be
    displayed, so a subset without the highest codes drew every state in a
    neighbour's colour (2026-09-22 audit S6-01). With ``[-0.5, K - 0.5]`` code
    ``c`` lands on palette entry ``c`` of a ``K``-colour listed colormap and on
    the centre of the ``c``-th band of any continuous one.
    """
    return Normalize(vmin=-0.5, vmax=len(x.alphabet) - 0.5)


def _state_color(cm, code: int):
    """One shared state-code -> color resolver for every renderer.

    Previously duplicated in five plot modules (cycle-8 audit F-33). The
    modulo keeps out-of-range codes drawable; renderers that accept declared
    external marker codes canonicalize their input first, so marker codes are
    in range and pick up the palette entries ``STSList._canonicalized``
    appends (``missing_color`` for nr, white for void).
    """
    if hasattr(cm, "colors"):
        return cm.colors[code % len(cm.colors)]
    return cm(code)


def ensure_cmap(x, cmap):
    """Resolve the state colours of ``x`` as a listed colormap, one entry per state.

    ``cmap`` may be ``None`` (the object's ``cpal``, else a default palette), a
    single colour, a sequence of colours (entry ``k`` colours the ``k``-th real
    state; the missing state takes the entry after the real states, else the
    object's ``missing_color``, and void is white, as in R; a list shorter than
    the real states, which R's Rd does not allow, is recycled over them so that
    each state keeps one colour), or a matplotlib colormap or its name, which is
    sampled at the centres of ``len(x.alphabet)`` equal bands -- the colours
    ``_state_code_norm`` gives the index plot. Every renderer therefore draws a
    state in the same colour; a named colormap used to hand the bar renderers
    its first, near-identical table entries (2026-09-22 audit S6-13).
    """
    state_count = max(len(x.alphabet), 1)
    if cmap is None and getattr(x, "cpal", None) is not None and len(x.cpal) >= len(x.alphabet):
        colors = [to_rgba(color) for color in x.cpal[: len(x.alphabet)]]
        return ListedColormap(colors)
    if cmap is None:
        if getattr(x, "cpal", None) is not None:
            warnings.warn(
                f"cpal has {len(x.cpal)} colors for {len(x.alphabet)} alphabet "
                "states; falling back to the tab20 palette",
                stacklevel=2,
            )
        # derive a default palette for the alphabet size
        base = plt.get_cmap("tab20")
        colors = [base(i % base.N) for i in range(len(x.alphabet))]
        return ListedColormap(colors)
    if isinstance(cmap, str):
        if is_color_like(cmap):
            return ListedColormap([to_rgba(cmap)])
        cmap = plt.get_cmap(cmap)
    if isinstance(cmap, Colormap):
        return ListedColormap([cmap((code + 0.5) / state_count) for code in range(state_count)])
    if not callable(cmap):
        try:
            colors = [to_rgba(color) for color in cmap]
        except TypeError:
            pass
        else:
            if colors:
                return ListedColormap(_user_state_colors(x, colors, state_count))
    return cmap


def _user_state_colors(x, colors, state_count: int) -> list:
    """Spread a user colour list over the real states; the markers keep their own colours.

    R appends ``missing.color`` after a user ``cpal`` and colours ``c(alphabet, nr)``
    with ``c(cpal, missing.color)`` in order (plot.stslist.R:83-90, seqplot.R:433-443):
    nr takes the entry after the real states when the list has one, else the
    object's missing colour, and void cells, which match no state, stay blank
    (TX-AUD-20261019-020). The list used to be recycled over every code, so
    missing cells took the first state's colour. R's Rd asks for one colour per
    state; given fewer, R recycles ``c(cpal, missing.color)`` cell by cell, so one
    state takes different colours (``AAA`` red, blue, darkgrey; batch-8 review
    M-3). Here a shorter list is recycled over the real states instead.
    """
    markers = {int(code) for code in (x.nr_code, x.void_code) if code is not None}
    real = [code for code in range(state_count) if code not in markers]
    by_code = {code: colors[rank % len(colors)] for rank, code in enumerate(real)}
    if x.nr_code is not None:
        by_code[int(x.nr_code)] = (
            colors[len(real)] if len(colors) > len(real) else to_rgba(x.missing_color)
        )
    if x.void_code is not None:
        by_code[int(x.void_code)] = to_rgba("white")
    return [by_code[code] for code in range(state_count)]
