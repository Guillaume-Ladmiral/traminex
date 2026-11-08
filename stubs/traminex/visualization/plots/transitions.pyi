from ...sequences import STSList as STSList

def transition_heatmap(
    x: STSList,
    *,
    normalize: bool = True,
    weighted: bool = True,
    figsize=(5, 4),
    seaborn: bool = False,
    ax=None,
): ...
