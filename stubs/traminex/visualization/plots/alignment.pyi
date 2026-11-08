from ...align_p2 import seqalign as seqalign
from ..palette import ensure_cmap as ensure_cmap

OP_COLORS: dict[str, str]
GAP_COLOR: str

def seqalignplot(
    x,
    indices,
    *,
    indel=1.0,
    sm=None,
    include_missing: bool = False,
    cpal=None,
    figsize=(8, 3),
    ax=None,
): ...
