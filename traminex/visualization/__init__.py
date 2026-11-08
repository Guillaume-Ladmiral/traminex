from __future__ import annotations

from .legend import legend_figure, state_legend, state_legend_figure
from .plots.chronogram import (
    chronogram,
    cluster_chronograms,
    cluster_chronograms_with_representatives,
)
from .plots.distribution import state_distribution_plot
from .plots.duration import duration_boxplot
from .plots.frequency import seqfrequency_plot
from .plots.implicative import plot_seqimplic
from .plots.index import (
    overlay_representatives_on_index_plot,
    sequence_index_plot,
    sequence_index_plot_grouped,
)
from .plots.representatives import seqplot_rf
from .plots.samm import plot_SAMM
from .plots.summary import seqplot_tentrop, seqplot_tentrop_m
from .plots.survival import seqsplot
from .plots.transitions import transition_heatmap
from .utils import compute_even_ticks, determine_layout, figure_to_buffer, save_plot

# Convenience re-exports for importing this subpackage directly. This is NOT the
# package's public plotting contract: `traminex.plotting` is, with 37 names routed
# by `traminex/__init__.py::_LAZY_EXPORTS["plotting"]`, all 37 in `PUBLIC_API`
# (`traminex.plotting.__all__` additionally lists the non-public `save_plot`
# helper, giving it 38 entries). The list below holds 23 names -- it omits many
# of the 37 public plot names (seqplot, seqfplot, seqrfplot, the plot_* event
# renderers and the rest) and includes helpers that are deliberately not public
# (`compute_even_ticks`, `determine_layout`, `figure_to_buffer`, `legend_figure`,
# `save_plot`). Nothing user-visible depends on the membership, so it is left as
# it is; the earlier comment "Public visualization API: curated re-exports"
# overstated its role (cycle-6 §8.3; counts re-measured cycle-8 F-48).
__all__ = [
    # plots
    "sequence_index_plot",
    "sequence_index_plot_grouped",
    "overlay_representatives_on_index_plot",
    "state_distribution_plot",
    "duration_boxplot",
    "chronogram",
    "cluster_chronograms",
    "cluster_chronograms_with_representatives",
    "transition_heatmap",
    "seqfrequency_plot",
    "seqplot_tentrop",
    "seqplot_tentrop_m",
    "seqplot_rf",
    "plot_SAMM",
    "plot_seqimplic",
    "seqsplot",
    # legends/palette
    "state_legend",
    "state_legend_figure",
    "legend_figure",
    # utils
    "compute_even_ticks",
    "determine_layout",
    "figure_to_buffer",
    "save_plot",
]
