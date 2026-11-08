"""Reproducible distance-kernel benchmarks with correctness checks.

Agreement is evaluated before timing. Timing samples therefore contain only the
selected backend call and never the comparison or report serialization work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

from traminex.distances import python_kernels as kernels

# ``DUR`` (with ``metric="chi2"``) and ``JSD_TRANS`` are library Numba routes that
# the CLI previously could not address -- and they are precisely the two whose
# agreement is non-trivial, because their kernels reassociate logs and divisions
# (cycle-4b F-016). Without them the agreement-before-timing contract exercised
# only routes that are byte-identical by construction.
METHODS = ("DHD", "DUR", "HAM", "JSD_TRANS", "LCP", "LCS", "OM", "RLCP")
BACKENDS = ("auto", "python", "numba", "pybind11")
NUMBA_METHODS = frozenset({"DHD", "DUR", "JSD_TRANS", "LCP", "LCS", "OM", "RLCP"})
PYBIND11_METHODS = frozenset({"DHD", "HAM", "LCP", "LCS", "OM", "RLCP"})
# The options each timed kernel runs with; the report records them as
# config.kernel_options (TX-AUD-20261011-007). DHD: the legacy mean-frequency
# variant is the shared Python/Numba/pybind11 contract; the default TraMineR DHD
# has no compiled kernel contract. DUR: only the chi2 metric has a Numba kernel;
# the l1/cos metrics do not, so the default would silently benchmark an
# unaccelerated route.
KERNEL_OPTIONS: dict[str, dict[str, Any]] = {
    "DHD": {"variant": "mean"},
    "DUR": {"metric": "chi2"},
    "OM": {"indel": 1.0, "subm": 2.0},
}


class BackendUnavailableError(RuntimeError):
    """Raised when a requested benchmark backend cannot honor the method."""


@dataclass(frozen=True)
class BenchmarkConfig:
    method: str = "OM"
    backends: tuple[str, ...] = ("auto",)
    n_sequences: int = 200
    sequence_length: int = 30
    n_states: int = 5
    seed: int = 42
    warmups: int = 1
    repeats: int = 3
    check_agreement: bool = True

    def __post_init__(self) -> None:
        method = self.method.upper()
        backends = tuple(backend.lower() for backend in self.backends)
        if method not in METHODS:
            raise ValueError(f"method must be one of: {', '.join(METHODS)}")
        if not backends or any(backend not in BACKENDS for backend in backends):
            raise ValueError(f"backend must be one of: {', '.join(BACKENDS)}")
        if self.n_sequences < 1:
            raise ValueError("N must be at least 1")
        if self.sequence_length < 1:
            raise ValueError("L must be at least 1")
        if self.n_states < 2:
            raise ValueError("K must be at least 2")
        if self.warmups < 0:
            raise ValueError("warmups must be non-negative")
        if self.repeats < 1:
            raise ValueError("repeats must be at least 1")
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "backends", backends)


def load_pybind11() -> Any | None:
    """Return the optional compiled module without hiding non-import failures."""
    try:
        from traminex import _pb_kernels
    except ImportError:
        return None
    return _pb_kernels


def generate_input(config: BenchmarkConfig) -> np.ndarray:
    """Generate the deterministic categorical-state trajectory matrix."""
    rng = np.random.default_rng(config.seed)
    return rng.integers(
        0,
        config.n_states,
        size=(config.n_sequences, config.sequence_length),
        dtype=np.int32,
    )


def input_metadata(matrix: np.ndarray, config: BenchmarkConfig) -> dict[str, Any]:
    contiguous = np.ascontiguousarray(matrix)
    return {
        "shape": list(contiguous.shape),
        "n_sequences": config.n_sequences,
        "sequence_length": config.sequence_length,
        "n_states": config.n_states,
        "seed": config.seed,
        "dtype": str(contiguous.dtype),
        "sha256": hashlib.sha256(contiguous.tobytes()).hexdigest(),
    }


def resolve_backend(requested: str, method: str) -> str:
    """Resolve a requested backend or raise an actionable availability error."""
    requested = requested.lower()
    method = method.upper()
    if requested == "python":
        return "python"
    if requested == "numba":
        if method not in NUMBA_METHODS:
            raise BackendUnavailableError(f"Numba backend does not support method {method}")
        if not kernels._HAVE_NUMBA:
            if kernels._NUMBA_IMPORTED:
                # numba is installed; TRAMINEX_NUMBA=off disabled dispatch (S7-20).
                raise BackendUnavailableError(
                    "Numba backend is unavailable: numba is installed but dispatch "
                    "is disabled by TRAMINEX_NUMBA=off"
                )
            raise BackendUnavailableError("Numba backend is unavailable; install the 'perf' extra")
        return "numba"
    if requested == "pybind11":
        if method not in PYBIND11_METHODS:
            raise BackendUnavailableError(f"pybind11 backend does not support method {method}")
        if load_pybind11() is None:
            raise BackendUnavailableError(
                "pybind11 backend is unavailable; build/install the explicit compiled profile"
            )
        return "pybind11"
    if requested != "auto":
        raise BackendUnavailableError(f"unknown backend: {requested}")
    if method in PYBIND11_METHODS and load_pybind11() is not None:
        return "pybind11"
    if method in NUMBA_METHODS and kernels._HAVE_NUMBA:
        return "numba"
    return "python"


def _python_kernel(matrix: np.ndarray, method: str) -> np.ndarray:
    kwargs = KERNEL_OPTIONS.get(method, {})
    return np.asarray(kernels.pairwise(matrix, method=method, **kwargs), dtype=np.float64)


def _python_with_numba_mode(matrix: np.ndarray, method: str, *, enabled: bool) -> np.ndarray:
    previous = kernels._HAVE_NUMBA
    kernels._HAVE_NUMBA = enabled
    try:
        return _python_kernel(matrix, method)
    finally:
        kernels._HAVE_NUMBA = previous


def _pybind11_kernel(module: Any, matrix: np.ndarray, method: str) -> np.ndarray:
    if method == "OM":
        options = KERNEL_OPTIONS["OM"]
        result = module.pairwise_om(matrix, options["indel"], options["subm"])
    elif method == "HAM":
        result = module.pairwise_ham_unit(matrix)
    elif method == "LCS":
        result = module.pairwise_lcs(matrix)
    elif method == "LCP":
        result = module.pairwise_lcp(matrix)
    elif method == "RLCP":
        result = module.pairwise_rlcp(matrix)
    elif method == "DHD":
        frequency = kernels._freq_by_pos(matrix, smooth_window=None)
        result = module.pairwise_dhd(matrix, np.asarray(frequency, dtype=np.float64), 0)
    else:  # guarded by resolve_backend
        raise BackendUnavailableError(f"pybind11 backend does not support method {method}")
    return np.asarray(result, dtype=np.float64)


def execute_backend(matrix: np.ndarray, method: str, *, backend: str) -> np.ndarray:
    """Execute one resolved backend call without timing or agreement work."""
    if backend == "python":
        return _python_with_numba_mode(matrix, method, enabled=False)
    if backend == "numba":
        return _python_with_numba_mode(matrix, method, enabled=True)
    if backend == "pybind11":
        module = load_pybind11()
        if module is None:
            raise BackendUnavailableError("pybind11 backend became unavailable")
        return _pybind11_kernel(module, matrix, method)
    raise BackendUnavailableError(f"backend must be resolved before execution: {backend}")


def measure_call(call: Callable[[], np.ndarray], *, warmups: int, repeats: int) -> list[float]:
    """Time only ``call`` after the requested number of untimed warmups."""
    for _ in range(warmups):
        call()
    timings: list[float] = []
    for _ in range(repeats):
        started = time.perf_counter()
        call()
        timings.append(time.perf_counter() - started)
    return timings


def _summary(timings: Sequence[float]) -> dict[str, float | int]:
    return {
        "count": len(timings),
        "min": min(timings),
        "max": max(timings),
        "mean": statistics.fmean(timings),
        "median": statistics.median(timings),
        "pstdev": statistics.pstdev(timings),
    }


def _agreement(reference: np.ndarray, candidate: np.ndarray, backend: str) -> dict[str, Any]:
    difference = np.abs(reference - candidate)
    max_difference = float(difference.max()) if difference.size else 0.0
    if backend == "pybind11":
        comparison = "array_equal"
        agrees = bool(np.array_equal(reference, candidate))
    else:
        # Numba's JSD_TRANS route is not byte-identical to the reference kernel
        # (it reassociates logs and divisions; measured max |delta| ~1e-16), so
        # Numba routes are compared with a tolerance that leaves four orders of
        # magnitude of margin while still failing on a real divergence.
        comparison = "allclose(rtol=1e-10, atol=1e-12)"
        agrees = bool(np.allclose(reference, candidate, rtol=1e-10, atol=1e-12))
    return {
        "checked": True,
        "reference_backend": "python",
        "comparison": comparison,
        "agrees": agrees,
        "max_abs_difference": max_difference,
    }


def _backend_metadata(backend: str) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "numba_available": bool(kernels._HAVE_NUMBA),
        "omp_num_threads_env": os.environ.get("OMP_NUM_THREADS"),
        "numba_num_threads_env": os.environ.get("NUMBA_NUM_THREADS"),
    }
    if backend == "numba":
        try:
            import numba

            metadata["numba_version"] = numba.__version__
            metadata["numba_threads"] = int(numba.get_num_threads())
        except ImportError:  # resolve_backend already provides the user-facing error
            pass
    if backend == "pybind11":
        module = load_pybind11()
        metadata["openmp_enabled"] = bool(
            module is not None
            and hasattr(module, "accel_openmp_enabled")
            and module.accel_openmp_enabled()
        )
        if module is not None and hasattr(module, "accel_max_threads"):
            metadata["openmp_max_threads"] = int(module.accel_max_threads())
    else:
        metadata["openmp_enabled"] = False
    return metadata


def _cpu_model() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name"):
                return line.partition(":")[2].strip()
    return platform.processor() or "unknown"


def _package_identity() -> dict[str, Any]:
    """The imported traminex: its file, version and a digest of its Python sources.

    The digest is SHA-256 over each ``.py`` file under the package directory, in
    order of relative POSIX path, of ``path NUL bytes NUL``; it tells apart two
    trees that report the same version (TX-AUD-20261011-007).
    """
    import traminex

    root = Path(traminex.__file__).parent
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*.py"), key=lambda item: item.relative_to(root).as_posix()):
        digest.update(path.relative_to(root).as_posix().encode() + b"\0")
        digest.update(path.read_bytes() + b"\0")
    return {
        "name": "traminex",
        "version": traminex.__version__,
        "file": traminex.__file__,
        "source_sha256": digest.hexdigest(),
    }


def _environment_metadata(matrix: np.ndarray, config: BenchmarkConfig) -> dict[str, Any]:
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "package": _package_identity(),
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "cpu": {"model": _cpu_model(), "logical_count": os.cpu_count()},
        "input": input_metadata(matrix, config),
    }


def run_suite(config: BenchmarkConfig) -> dict[str, Any]:
    """Run all requested backends, refusing any agreement mismatch."""
    matrix = generate_input(config)
    results: list[dict[str, Any]] = []
    for requested in config.backends:
        resolved = resolve_backend(requested, config.method)

        # Correctness is deliberately outside measure_call.
        if config.check_agreement:
            reference = execute_backend(matrix, config.method, backend="python")
            candidate = execute_backend(matrix, config.method, backend=resolved)
            agreement = _agreement(reference, candidate, resolved)
            if not agreement["agrees"]:
                raise RuntimeError(
                    f"{resolved} disagrees with Python for {config.method}; "
                    f"max_abs_difference={agreement['max_abs_difference']}"
                )
        else:
            agreement = {"checked": False}

        timings = measure_call(
            lambda: execute_backend(matrix, config.method, backend=resolved),
            warmups=config.warmups,
            repeats=config.repeats,
        )
        results.append(
            {
                "requested_backend": requested,
                "resolved_backend": resolved,
                "agreement": agreement,
                "timings_seconds": timings,
                "summary_seconds": _summary(timings),
                "backend_metadata": _backend_metadata(resolved),
            }
        )

    return {
        "schema_version": 1,
        "config": {**asdict(config), "kernel_options": KERNEL_OPTIONS.get(config.method, {})},
        "metadata": _environment_metadata(matrix, config),
        "results": results,
    }


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _state_count(value: str) -> int:
    parsed = int(value)
    if parsed < 2:
        raise argparse.ArgumentTypeError("must be at least 2")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark traminex distance backends with untimed agreement checks"
    )
    parser.add_argument(
        "--method",
        choices=METHODS,
        default="OM",
        type=str.upper,
        help=(
            "distance to time (default: OM, indel 1 and substitution cost 2). DHD times the"
            " frequency variant variant='mean', not R's default DHD; DUR the chi2 metric."
            " The report records these as config.kernel_options."
        ),
    )
    parser.add_argument(
        "--backend",
        action="append",
        choices=BACKENDS,
        dest="backends",
        help="backend to run; repeat for comparisons (default: auto)",
    )
    parser.add_argument("--N", "--n-sequences", type=_positive_int, default=200)
    parser.add_argument("--L", "--sequence-length", type=_positive_int, default=30)
    parser.add_argument("--K", "--n-states", type=_state_count, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--warmups", type=_nonnegative_int, default=1)
    parser.add_argument("--repeats", type=_positive_int, default=3)
    parser.add_argument("--no-agreement", action="store_true")
    parser.add_argument("--json", metavar="PATH", help="write JSON report; use '-' for stdout")
    return parser


def _print_human(report: dict[str, Any]) -> None:
    config = report["config"]
    shape = report["metadata"]["input"]["shape"]
    print(
        f"method={config['method']} input={shape[0]}x{shape[1]} "
        f"K={report['metadata']['input']['n_states']} seed={config['seed']}"
    )
    for result in report["results"]:
        summary = result["summary_seconds"]
        agreement = result["agreement"]
        print(
            f"backend={result['requested_backend']} resolved={result['resolved_backend']} "
            f"mean={summary['mean']:.6f}s min={summary['min']:.6f}s "
            f"agreement={agreement.get('agrees', 'not-checked')}"
        )


def _write_json(report: dict[str, Any], destination: str) -> None:
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if destination == "-":
        print(rendered, end="")
        return
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        config = BenchmarkConfig(
            method=args.method,
            backends=tuple(args.backends or ("auto",)),
            n_sequences=args.N,
            sequence_length=args.L,
            n_states=args.K,
            seed=args.seed,
            warmups=args.warmups,
            repeats=args.repeats,
            check_agreement=not args.no_agreement,
        )
        report = run_suite(config)
    except (BackendUnavailableError, ValueError, RuntimeError) as error:
        parser.error(str(error))
    if args.json:
        _write_json(report, args.json)
    else:
        _print_human(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
