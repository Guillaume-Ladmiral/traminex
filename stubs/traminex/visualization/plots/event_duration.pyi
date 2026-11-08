from typing import Optional, Sequence

def seqedplot(
    eseq,
    *,
    type: str = "survival",
    group: Optional[Sequence] = None,
    ignore: Optional[Sequence[str]] = None,
    breaks: int = 20,
    ages: Optional[tuple[float, float]] = None,
    cpal=None,
    figsize=(8, 4),
    axes=None,
): ...
