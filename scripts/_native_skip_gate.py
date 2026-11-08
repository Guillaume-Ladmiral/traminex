"""A pytest plugin for the native (pybind11) gate: a compiled-kernel skip fails the run.

docs/DEV_BUILD.md: "A Python fallback or skipped compiled differential module is a failed
pybind11 gate even when pytest exits zero." Nothing enforced it (TX-AUD-20261011-075,
TX-AUD-20261019-077): with the module missing, those tests skip and pytest exits 0.
``scripts/run_pybind11_tests.sh`` loads this plugin. Skips for other reasons (an optional
tool such as Graphviz or Rscript, a platform property) keep pytest's meaning.
"""

from __future__ import annotations

import pytest

#: The skip reasons in tests/ that mean the compiled kernels are unbuilt or not in use.
#: ``tests/test_native_gates.py`` requires this list to equal the reasons the tests use.
COMPILED_SKIP_REASONS = (
    "compiled _pb_kernels not built",
    "pybind11 module is not built",
    "active backend is not the compiled kernel",
)


def _skip_reason(report) -> str:
    longrepr = report.longrepr
    message = (
        str(longrepr[2]) if isinstance(longrepr, tuple) and len(longrepr) == 3 else str(longrepr)
    )
    return message.removeprefix("Skipped: ")


class NativeSkipGate:
    """Record compiled-kernel skips; turn a passing session into a failed one if any occurred."""

    def __init__(self) -> None:
        self.compiled_skips: list[tuple[str, str]] = []

    def _record(self, report) -> None:
        if report.skipped and not hasattr(report, "wasxfail"):
            reason = _skip_reason(report)
            if reason in COMPILED_SKIP_REASONS:
                self.compiled_skips.append((report.nodeid, reason))

    def pytest_runtest_logreport(self, report) -> None:
        self._record(report)

    def pytest_collectreport(self, report) -> None:
        self._record(report)

    def pytest_terminal_summary(self, terminalreporter) -> None:
        if self.compiled_skips:
            terminalreporter.section("native gate: compiled-kernel tests skipped")
            for nodeid, reason in self.compiled_skips:
                terminalreporter.line(f"{nodeid}: {reason}")

    @pytest.hookimpl(trylast=True)
    def pytest_sessionfinish(self, session, exitstatus) -> None:
        if self.compiled_skips and exitstatus == pytest.ExitCode.OK:
            session.exitstatus = pytest.ExitCode.TESTS_FAILED
