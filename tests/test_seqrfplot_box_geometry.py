"""seqrfplot's boxes lie inside their bands, as R draws them (TX-AUD-20261019-018).

R draws the boxes at ``(cumsum(h) - h/2) / sum(h) * k`` with widths ``boxwex * h / max(h)``,
``boxwex = 0.8`` (dissrf.R:249; plot.seqrf.R:79-90; graphics::bxp), all on one 0..k axis, so
every box lies in its group's band. traminex draws the bands on a 0..sum(h) axis, so R's
geometry maps to widths ``0.8 * h / max(h) * sum(h) / k``. The boxes were 0.6 * k * h / sum(h)
tall on that axis and overlapped (0.6 for ten equal groups whose bands are 0.1 apart).
"""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

ROWS = [list(r) for r in ("AAAB", "AABB", "ABBB", "BBBB", "BBBA", "BBAA", "BAAA", "AAAA", "ABAB")]


def _boxes(k: int, **rf_kwargs):
    x = tx.seqdef(ROWS, alphabet=["A", "B"])
    diss = np.asarray(tx.seqdist(x, method="HAM"), dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig, (_, right) = tx.seqrfplot(x, diss, k=k, **rf_kwargs)
        rf = tx.seqrf(x, diss, k=k, **rf_kwargs)
    extents = [
        patch.get_path().get_extents()
        for patch in right.patches
        if isinstance(patch, mpatches.PathPatch)
    ]
    plt.close(fig)
    return np.asarray(rf["heights"], dtype=float), [(e.y0, e.y1) for e in extents]


@pytest.mark.parametrize("k, rf_kwargs", [(3, {}), (4, {"grp_meth": "first"})])
def test_every_box_lies_in_its_band_with_r_s_width(k: int, rf_kwargs: dict) -> None:
    heights, boxes = _boxes(k, **rf_kwargs)
    assert len(boxes) == k
    tops = np.cumsum(heights)
    bottoms = tops - heights
    widths = 0.8 * heights / heights.max() * heights.sum() / k
    for (y0, y1), bottom, top, width in zip(boxes, bottoms, tops, widths, strict=True):
        assert y1 - y0 == pytest.approx(width, abs=1e-12)
        assert bottom - 1e-12 <= y0 and y1 <= top + 1e-12
