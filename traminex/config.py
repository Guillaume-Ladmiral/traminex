from __future__ import annotations

import json
import math
import os
from json import JSONDecodeError
from typing import Any, Dict, Iterable, Mapping

_ENV_CALIBRATED_PATHS = "TRAMINEX_CALIBRATED_PARAMS"


def _normalize_paths(
    paths: str | os.PathLike[str] | Iterable[str | os.PathLike[str]] | None,
) -> list[str]:
    """Normalize optional path input into a deterministic path list."""
    if paths is None:
        raw = os.environ.get(_ENV_CALIBRATED_PATHS, "").strip()
        if not raw:
            return []
        return [p.strip() for p in raw.split(os.pathsep) if p.strip()]
    if isinstance(paths, (str, os.PathLike)):
        # A single path object is one path, not an iterable (audit S2-17).
        return [os.fspath(paths)]
    return [str(p).strip() for p in paths if str(p).strip()]


def _validate_calibrated_payload(payload: Dict[str, Any], fp: str) -> None:
    """Refuse a TRATE smoothing the cost builder would silently misuse.

    A negative, NaN or infinite value (or a JSON boolean) otherwise surfaced
    later as an unrelated cost error (audit S2-17).
    """
    trate = payload.get("trate")
    if not isinstance(trate, dict) or "smoothing" not in trate:
        return
    value = trate["smoothing"]
    try:
        smoothing = math.nan if isinstance(value, bool) else float(value)
    except (TypeError, ValueError):
        smoothing = math.nan
    if not math.isfinite(smoothing) or smoothing < 0:
        raise ValueError(
            f"trate.smoothing in calibrated params file {fp} must be a finite "
            f"non-negative number, got {value!r}"
        )


#: The keys each calibrated profile reads from its section of the file, as
#: ``scripts/calibrate_parity.py`` writes them.
_PROFILE_KEYS: Dict[str, tuple[str, ...]] = {
    "trate": ("smoothing", "mode", "normalize"),
    "dhd": ("variant", "smooth_window"),
}


def _load_first(paths: object) -> tuple[Dict[str, Any] | None, str | None, list[str]]:
    """The first calibration file that exists: ``(payload, file, paths tried)``."""
    tried: list[str] = []
    for p in _normalize_paths(paths):  # type: ignore[arg-type]
        fp = os.path.abspath(os.path.expanduser(p))
        tried.append(fp)
        if not os.path.exists(fp):
            continue
        try:
            with open(fp, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except FileNotFoundError:
            continue
        except JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in calibrated params file: {fp}") from exc
        except OSError as exc:
            raise ValueError(f"Unable to read calibrated params file: {fp}") from exc
        if not isinstance(payload, dict):
            raise ValueError(f"Calibrated params file must contain a JSON object: {fp}")
        _validate_calibrated_payload(payload, fp)
        return payload, fp, tried
    return None, None, tried


def _calibrated_profile_values(
    kind: str, explicit: Mapping[str, object], paths: object = None
) -> Dict[str, object]:
    """Resolve a calibrated ``trate`` or ``dhd`` profile (TX-AUD-20261010-032).

    A calibrated profile is an explicit request, so it fails loudly instead of
    falling back to the defaults: when a value it needs is not given in
    ``explicit``, the first existing file of ``paths`` (or of
    ``TRAMINEX_CALIBRATED_PARAMS``) must load and hold the profile's section,
    and that section may hold only the profile's keys. A value given in
    ``explicit`` overrides the file's, for both profiles. Keys found in neither
    are left out, so the caller's own defaults apply.
    """
    keys = _PROFILE_KEYS[kind]
    values: Dict[str, object] = {key: explicit[key] for key in keys if key in explicit}
    if len(values) == len(keys):
        return values
    payload, fp, tried = _load_first(paths)
    if payload is None:
        raise ValueError(
            f"calibrated {kind} profile requested but no calibration file loaded "
            f"(paths tried: {tried}); pass calibrated_params_path or set "
            f"{_ENV_CALIBRATED_PATHS}"
        )
    section = payload.get(kind)
    if not isinstance(section, dict):
        raise ValueError(f"calibration file {fp} has no '{kind}' section")
    unknown = sorted(set(section) - set(keys))
    if unknown:
        raise ValueError(
            f"unknown key(s) {unknown} in the '{kind}' section of {fp}; "
            f"expected only {list(keys)}"
        )
    for key in keys:
        if key not in values and key in section:
            values[key] = section[key]
    return values


def load_calibrated_params(
    paths: list[str] | str | os.PathLike[str] | None = None,
) -> Dict[str, Any]:
    """Load calibrated parameters JSON from explicit paths only.

    Deterministic behavior:
    - No implicit search in CWD, project root, or test folders.
    - ``paths`` is one path (a string or a path object) or a list of paths;
      the first file that exists is loaded.
    - If ``paths`` is None, uses env var ``TRAMINEX_CALIBRATED_PARAMS`` (pathsep-separated).
    - If no explicit paths are provided, returns {}.
    - Returns {} if no file loads successfully.
    - A ``trate.smoothing`` value that is not a finite non-negative number
      is refused with a ``ValueError`` naming the file.
    """
    payload, _, _ = _load_first(paths)
    return {} if payload is None else payload
