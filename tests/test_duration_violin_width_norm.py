"""``duration_boxplot(by_state=True)`` must width-normalise its violins on every
supported seaborn (audit S7-11, 2026-09-22).

What was wrong: the plot passed ``density_norm="width"`` to
``seaborn.violinplot``. That keyword first exists in seaborn 0.13; the declared
floor is ``seaborn>=0.12``, whose ``violinplot`` accepts and silently ignores
unknown ``**kwargs``, so on 0.12.x the violins were area-scaled instead of
width-normalised (measured on seaborn 0.12.2: extents 0.287 / 0.8 / 0.268
against 0.8 / 0.8 / 0.8 on 0.13.2). seaborn 0.13 renamed ``scale`` to
``density_norm`` and deprecated ``scale``.

Source of truth: the two seaborn signatures (0.12: ``scale="area"`` default and
``**kwargs``; 0.13: ``density_norm``). The plot now passes whichever keyword
the installed ``violinplot`` declares.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("matplotlib")
sns = pytest.importorskip("seaborn")

import traminex as tx  # noqa: E402


def _sequences():
    rows = [["A"] * a + ["B"] * b + ["C"] * (10 - a - b) for a, b in ((1, 0), (5, 2), (8, 1))]
    rows += [["A"] * 3 + ["C"] * 7, ["B"] * 2 + ["A"] * 8]
    return tx.seqdef(np.array(rows, dtype=object))


def _recording_violinplot(calls):
    def violinplot_012(
        data=None,
        *,
        x=None,
        y=None,
        hue=None,
        order=None,
        hue_order=None,
        bw="scott",
        cut=2,
        scale="area",
        scale_hue=True,
        gridsize=100,
        width=0.8,
        inner="box",
        split=False,
        dodge=True,
        orient=None,
        linewidth=None,
        color=None,
        palette=None,
        saturation=0.75,
        ax=None,
        **kwargs,
    ):
        calls.append({"scale": scale, **kwargs})
        return ax

    return violinplot_012


@pytest.mark.parametrize("with_labels", [False, True])
def test_seaborn_012_signature_receives_scale_width(monkeypatch, with_labels):
    calls: list[dict] = []
    monkeypatch.setattr(sns, "violinplot", _recording_violinplot(calls))
    x = _sequences()
    labels = np.array([0, 0, 1, 1, 1]) if with_labels else None

    tx.duration_boxplot(x, labels, by_state=True)

    assert calls == [{"scale": "width"}]


def test_installed_seaborn_draws_width_normalised_violins():
    import warnings

    from matplotlib.collections import PolyCollection

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fig, ax = tx.duration_boxplot(_sequences(), by_state=True)
    # seaborn >= 0.13 deprecates ``scale``; the plot must use the spelling the
    # installed version declares.
    assert not [w for w in caught if "scale" in str(w.message)], [str(w.message) for w in caught]
    extents = [
        round(float(path.vertices[:, 1].max() - path.vertices[:, 1].min()), 6)
        for collection in ax.collections
        if isinstance(collection, PolyCollection)
        for path in collection.get_paths()
    ]

    assert len(extents) == 3
    assert len(set(extents)) == 1, extents
