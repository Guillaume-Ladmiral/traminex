from .sequences import STSList as STSList

def frequent_subsequences(
    x: STSList,
    *,
    max_len: int = 3,
    contiguous: bool = True,
    min_support: int = 2,
    support_type: str = "occurrences",
) -> list[tuple[tuple[str, ...], int, float]]: ...
