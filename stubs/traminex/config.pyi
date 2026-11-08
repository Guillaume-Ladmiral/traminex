import os
from typing import Any

def load_calibrated_params(
    paths: list[str] | str | os.PathLike[str] | None = None,
) -> dict[str, Any]: ...
