import os
import sys

import pytest

# Ensure project root is on sys.path so `import traminex` works when running tests directly
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def pytest_report_header(config):
    """Name the ICU library traminex loads (TX-AUD-20261019-053, owner decision of 2026-10-22).

    The automatic state order follows that library, as R's follows its own; a version that differs from R's is
    shown here, not refused. Never raises, so an unexpected library or lookup cannot stop the run.
    """
    try:
        import ctypes.util

        library = ctypes.util.find_library("icui18n")
    except Exception as exc:
        return f"ICU: lookup failed ({type(exc).__name__}); traminex needs ICU (docs/DEV_BUILD.md)"
    if library is None:
        return "ICU: not found (traminex needs it for R's label order; docs/DEV_BUILD.md)"
    try:
        from traminex._categorical_order import _icu_version

        version = _icu_version(library)
    except Exception:
        return f"ICU: {library} (version unreadable; traminex cannot use it; docs/DEV_BUILD.md)"
    return (
        f"ICU: {library} (major version {version}); the automatic state order follows it, and an "
        "explicit alphabet= keeps R's order under any ICU (docs/DEV_BUILD.md)"
    )


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Show the ICU line at the end when ``-q`` hid the header.

    ``make test-python`` and ``make test-pybind11`` pass ``-q``, which hides the header above, so their runs
    named no ICU library (2026-10-25 study, finding F-03). A run shows the line once, in either place.
    """
    if config.option.verbose < 0:
        terminalreporter.write_line(pytest_report_header(config))


@pytest.fixture(autouse=True)
def _close_matplotlib_figures():
    """Close all matplotlib figures after each test.

    The RP7 visualization tests open many figures; without this the suite trips
    matplotlib's "More than 20 figures have been opened" resource warning. Closing
    per test keeps memory bounded and the warning summary clean. No-op when
    matplotlib is not importable.
    """
    yield
    try:
        import matplotlib.pyplot as plt

        plt.close("all")
    except Exception:
        pass
