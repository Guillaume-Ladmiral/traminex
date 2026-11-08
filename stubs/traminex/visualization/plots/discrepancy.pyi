from ...seqdiff import seqdiff as seqdiff

def seqdiffplot(
    x,
    group,
    *,
    method: str = "LCS",
    norm="auto",
    weighted: bool = True,
    stat: str = "Pseudo R2",
    figsize=(8, 4),
    ax=None,
    **kwargs,
): ...
