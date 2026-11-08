from io import BytesIO
from typing import Sequence

def compute_even_ticks(
    length: int, *, max_labels: int = 10, include_endpoints: bool = True
) -> list[int]: ...
def determine_layout(
    num_items: int,
    *,
    layout: str = "column",
    nrows: int | None = None,
    ncols: int | None = None,
    default_cols: int = 3,
) -> tuple[int, int]: ...
def legend_figure(
    labels: Sequence[str],
    colors: Sequence,
    *,
    ncol: int = 5,
    figsize: tuple[float, float] = (8, 1),
    fontsize: float | str = 10,
): ...
def figure_to_buffer(
    fig, *, fmt: str = "png", dpi: int = 200, bbox_inches: str = "tight"
) -> BytesIO: ...
def save_plot(
    fig_or_tuple, filepath: str, *, dpi: int = 300, bbox_inches: str | None = "tight", **kwargs
) -> None: ...
