"""Import-order helpers for developer commands that must use the current checkout."""

from __future__ import annotations

import sys
from importlib.machinery import PathFinder


def prefer_checkout_imports() -> None:
    """Make ordinary ``sys.path`` lookup precede editable-install redirects."""
    redirect_finders = []
    for finder in sys.meta_path:
        finder_type = type(finder)
        identifier = f"{finder_type.__module__}.{finder_type.__qualname__}".lower()
        if "editable" in identifier or "redirect" in identifier:
            redirect_finders.append(finder)

    reordered_finders = []
    path_finder_seen = False
    for finder in sys.meta_path:
        if finder is PathFinder:
            reordered_finders.append(finder)
            reordered_finders.extend(redirect_finders)
            path_finder_seen = True
        elif finder not in redirect_finders:
            reordered_finders.append(finder)

    if not path_finder_seen:  # pragma: no cover - a broken interpreter configuration
        raise RuntimeError("Python PathFinder is unavailable")
    sys.meta_path[:] = reordered_finders
