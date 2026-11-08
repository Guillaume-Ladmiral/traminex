from __future__ import annotations

import json

import numpy as np
import pytest

from benchmarks import run_benchmarks as runner


def _config(**overrides):
    values = {
        "method": "HAM",
        "backends": ("python",),
        "n_sequences": 4,
        "sequence_length": 5,
        "n_states": 3,
        "seed": 17,
        "warmups": 0,
        "repeats": 1,
        "check_agreement": True,
    }
    values.update(overrides)
    return runner.BenchmarkConfig(**values)


def test_input_generation_is_deterministic_and_records_identity():
    config = _config()
    first = runner.generate_input(config)
    second = runner.generate_input(config)

    assert np.array_equal(first, second)
    assert first.dtype == np.int32
    metadata = runner.input_metadata(first, config)
    assert metadata["shape"] == [4, 5]
    assert metadata["n_states"] == 3
    assert metadata["seed"] == 17
    assert len(metadata["sha256"]) == 64


def test_python_benchmark_records_agreement_outside_timing(monkeypatch):
    events: list[str] = []
    original_execute = runner.execute_backend
    original_measure = runner.measure_call

    def execute(*args, **kwargs):
        events.append(f"execute:{kwargs['backend']}")
        return original_execute(*args, **kwargs)

    def measure(call, *, warmups, repeats):
        events.append("measure:start")
        result = original_measure(call, warmups=warmups, repeats=repeats)
        events.append("measure:end")
        return result

    monkeypatch.setattr(runner, "execute_backend", execute)
    monkeypatch.setattr(runner, "measure_call", measure)

    report = runner.run_suite(_config())
    result = report["results"][0]

    assert result["agreement"]["agrees"] is True
    assert result["agreement"]["reference_backend"] == "python"
    assert result["summary_seconds"]["count"] == 1
    assert events[:2] == ["execute:python", "execute:python"]
    assert events[2] == "measure:start"
    assert events[-1] == "measure:end"


def test_auto_resolution_and_unavailable_backends_are_explicit(monkeypatch):
    monkeypatch.setattr(runner, "load_pybind11", lambda: None)
    monkeypatch.setattr(runner.kernels, "_HAVE_NUMBA", False)

    assert runner.resolve_backend("auto", "OM") == "python"
    with pytest.raises(runner.BackendUnavailableError, match="pybind11 backend is unavailable"):
        runner.resolve_backend("pybind11", "OM")
    with pytest.raises(runner.BackendUnavailableError, match="Numba backend is unavailable"):
        runner.resolve_backend("numba", "OM")
    with pytest.raises(runner.BackendUnavailableError, match="does not support method HAM"):
        runner.resolve_backend("numba", "HAM")


def _numba_skip_reason() -> str | None:
    """Why the Numba comparison cannot run here, or ``None`` if it can.

    ``_HAVE_NUMBA`` is ``_NUMBA_IMPORTED and _resolve_numba_mode() == "auto"``
    (``traminex/distances/python_kernels.py``), so it is false both when numba is
    absent *and* when the supported ``TRAMINEX_NUMBA=off`` is set. Reporting
    "Numba is not installed" for the second case tells a log reader the
    environment lacks numba when it is in fact a deliberate configuration --
    the reference-kernel suite reported exactly that while numba was installed
    (cycle-6 F-007).
    """
    if not runner.kernels._NUMBA_IMPORTED:
        return "Numba is not installed"
    if not runner.kernels._HAVE_NUMBA:
        return "Numba dispatch is disabled by TRAMINEX_NUMBA=off"
    return None


@pytest.mark.parametrize(
    ("imported", "dispatching", "expected"),
    [
        (False, False, "Numba is not installed"),
        (True, False, "Numba dispatch is disabled by TRAMINEX_NUMBA=off"),
        (True, True, None),
    ],
)
def test_the_numba_skip_reason_distinguishes_absent_from_disabled(
    monkeypatch, imported, dispatching, expected
):
    monkeypatch.setattr(runner.kernels, "_NUMBA_IMPORTED", imported, raising=True)
    monkeypatch.setattr(runner.kernels, "_HAVE_NUMBA", dispatching, raising=True)
    assert _numba_skip_reason() == expected


def test_numba_agrees_with_python_when_available():
    reason = _numba_skip_reason()
    if reason is not None:
        pytest.skip(reason)

    report = runner.run_suite(_config(method="LCS", backends=("numba",), warmups=1))
    result = report["results"][0]
    assert result["resolved_backend"] == "numba"
    assert result["agreement"]["agrees"] is True
    assert result["agreement"]["comparison"] == "allclose(rtol=1e-10, atol=1e-12)"


@pytest.mark.skipif(runner.load_pybind11() is None, reason="pybind11 module is not built")
@pytest.mark.parametrize("method", ["HAM", "OM", "LCS", "RLCP", "DHD"])
def test_pybind11_promised_methods_are_byte_identical(method):
    report = runner.run_suite(_config(method=method, backends=("pybind11",)))
    result = report["results"][0]
    assert result["resolved_backend"] == "pybind11"
    assert result["agreement"] == {
        "checked": True,
        "reference_backend": "python",
        "comparison": "array_equal",
        "agrees": True,
        "max_abs_difference": 0.0,
    }


def test_cli_writes_machine_readable_report(tmp_path):
    output = tmp_path / "nested" / "benchmark.json"
    exit_code = runner.main(
        [
            "--method",
            "HAM",
            "--backend",
            "python",
            "--N",
            "3",
            "--L",
            "4",
            "--K",
            "2",
            "--warmups",
            "0",
            "--repeats",
            "1",
            "--seed",
            "23",
            "--json",
            str(output),
        ]
    )

    assert exit_code == 0
    report = json.loads(output.read_text())
    assert report["schema_version"] == 1
    assert report["config"]["method"] == "HAM"
    assert report["metadata"]["package"]["version"] == "0.2.3"
    assert report["metadata"]["python"]["executable"]
    assert report["metadata"]["platform"]["system"] == "Linux"
    assert report["metadata"]["cpu"]["logical_count"]
    assert report["metadata"]["input"]["shape"] == [3, 4]
    assert report["results"][0]["timings_seconds"]


@pytest.mark.parametrize(
    "arguments",
    [
        ["--N", "0"],
        ["--L", "0"],
        ["--K", "1"],
        ["--warmups", "-1"],
        ["--repeats", "0"],
    ],
)
def test_cli_rejects_unsafe_or_meaningless_dimensions(arguments):
    with pytest.raises(SystemExit):
        runner.main(arguments)


@pytest.mark.parametrize(
    ("imported", "expected"),
    [
        (False, "install the 'perf' extra"),
        (True, "TRAMINEX_NUMBA=off"),
    ],
)
def test_the_numba_refusal_names_the_actual_cause(monkeypatch, imported, expected):
    """``--backend numba`` must not tell a user with numba to install it (S7-20).

    With ``TRAMINEX_NUMBA=off`` numba is importable but dispatch is disabled, and
    the benchmark still said "install the 'perf' extra" (audit S7-20,
    2026-09-22). ``_NUMBA_IMPORTED`` distinguishes the two causes.
    """
    monkeypatch.setattr(runner.kernels, "_NUMBA_IMPORTED", imported, raising=True)
    monkeypatch.setattr(runner.kernels, "_HAVE_NUMBA", False, raising=True)

    with pytest.raises(
        runner.BackendUnavailableError, match="Numba backend is unavailable"
    ) as caught:
        runner.resolve_backend("numba", "OM")

    assert expected in str(caught.value)


# TX-AUD-20261011-007: a report did not say which kernel options it timed (DHD times the
# frequency variant "mean", not R's default DHD) nor which traminex tree was imported.


def _cli_report(tmp_path, method):
    output = tmp_path / f"{method}.json"
    arguments = ["--method", method, "--backend", "python", "--N", "3", "--L", "4"]
    arguments += ["--K", "2", "--warmups", "0", "--repeats", "1", "--json", str(output)]
    assert runner.main(arguments) == 0
    return json.loads(output.read_text())


@pytest.mark.parametrize(
    ("method", "options"),
    [
        ("DHD", {"variant": "mean"}),
        ("DUR", {"metric": "chi2"}),
        ("OM", {"indel": 1.0, "subm": 2.0}),
        ("HAM", {}),
    ],
)
def test_the_report_records_the_kernel_options_it_timed(tmp_path, method, options):
    assert _cli_report(tmp_path, method)["config"]["kernel_options"] == options


def test_the_report_identifies_the_imported_traminex_source(tmp_path):
    import hashlib
    from pathlib import Path

    import traminex

    package = _cli_report(tmp_path, "HAM")["metadata"]["package"]

    root = Path(traminex.__file__).parent
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*.py"), key=lambda p: p.relative_to(root).as_posix()):
        digest.update(path.relative_to(root).as_posix().encode() + b"\0")
        digest.update(path.read_bytes() + b"\0")
    assert package["file"] == traminex.__file__
    assert package["version"] == traminex.__version__
    assert package["source_sha256"] == digest.hexdigest()


def test_the_help_names_the_kernels_timed():
    text = " ".join(runner.build_parser().format_help().split())

    assert "DHD times the frequency variant variant='mean', not R's default DHD" in text
    assert "DUR the chi2 metric" in text
