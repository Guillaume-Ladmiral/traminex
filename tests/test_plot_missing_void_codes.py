import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from traminex import seqdef, sequence_index_plot, state_distribution_plot


def test_sequence_index_plot_applies_overlays_for_missing_and_void_codes():
    arr = np.array(
        [
            ["A", "*", "B"],
            ["%", "A", "B"],
        ],
        dtype=object,
    )
    x = seqdef(arr, nr="*", void="%")
    fig, ax = sequence_index_plot(x)
    # Base image + missing overlay + void overlay
    assert len(ax.images) >= 3
    plt.close(fig)


def test_state_distribution_plot_uses_missing_and_void_alphas_by_code():
    arr = np.array(
        [
            ["A", "*", "B"],
            ["%", "A", "B"],
        ],
        dtype=object,
    )
    x = seqdef(arr, nr="*", void="%")
    fig, ax = state_distribution_plot(
        x,
        include_missing=True,
        mode="d",
        missing_alpha=0.11,
        void_alpha=0.22,
    )
    alphas = [p.get_alpha() for p in ax.patches if p.get_alpha() is not None]
    assert any(abs(a - 0.11) < 1e-12 for a in alphas)
    assert any(abs(a - 0.22) < 1e-12 for a in alphas)
    plt.close(fig)
