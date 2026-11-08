"""plot_subseqelistchisq colours bars by R's residual bins (TX-AUD-20261010-045).

R's plot.subseqelistchisq colours each bar by the bin of its Pearson residual,
cut at +/- -qnorm(c(0.05, 0.01)) (Bonferroni-adjusted when used), with the five
Brewer "RdBu" colours, so a positive residual is blue, and titles each panel with
the group label (seqecmpgroup.R:105-142). traminex used a continuous red-blue
scale with positive red and titled "Group <g>". The horizontal layout stays a
recorded re-expression. Fixture: ``scripts/r_goldens/make_goldens_ecmpgroup_colours.R``
(colours and titles traced from R's barplot calls).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402

import traminex as tx  # noqa: E402

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(HERE / "ecmpgroup_colours.csv", dtype=str, keep_default_na=False)
TSE = pd.read_csv(HERE / "ecmpgroup_colours_tse.csv")
GROUP = ["g1"] * 20 + ["g2"] * 20


def _comparison(method: str):
    eseq = tx.seqecreate(id=TSE["id"], timestamp=TSE["time"], event=TSE["event"])
    subsequences = tx.seqefsub(eseq, pmin_support=0.05)
    return tx.seqecmpgroup(subsequences, group=GROUP, method=method)


@pytest.mark.parametrize("method", ["chisq", "bonferroni"])
def test_bar_colours_and_titles_are_r(method: str) -> None:
    fig, axes = tx.plot_subseqelistchisq(_comparison(method))
    r = GOLDEN[GOLDEN["method"] == method]
    for ax, (panel, rows) in zip(axes, r.groupby("panel")):
        rows = rows.sort_values("position", key=lambda c: c.astype(int))
        assert ax.get_title() == rows["title"].iloc[0]
        # The panels share the y axis; its labels show on the first panel.
        labels = [label.get_text() for label in axes[0].get_yticklabels()]
        assert labels == rows["subseq"].tolist()
        colours = [to_hex(patch.get_facecolor()).upper() for patch in ax.patches]
        assert colours == rows["colour"].str.upper().tolist()
    assert len(axes) == r["panel"].nunique()
    plt.close(fig)


def test_positive_residuals_are_blue() -> None:
    chisq = GOLDEN[GOLDEN["method"] == "chisq"]
    first = chisq[(chisq["title"] == "g1") & (chisq["subseq"] == "(a)")]["colour"].item()
    assert first.upper() == "#0571B0"  # g1 has more (a) than expected
